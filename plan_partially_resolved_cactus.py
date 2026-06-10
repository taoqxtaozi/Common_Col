#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plan_partially_resolved_cactus.py

A readable planner that:
1. reads a rooted / meaningfully rooted Newick tree,
2. conditionally keeps or reroots the input tree according to the reference rule,
3. auto-names unnamed internal nodes,
4. detects hierarchical unresolved regions,
5. creates nested folders + `*_aln*.txt`,
6. writes an `instruction.txt` that lists the Cactus commands in dependency order.

Important conventions
---------------------
- Rerooting rule:
    * if the input tree is fully binary, keep the input rooting unchanged;
    * else, if the reference is already a direct child of the input top root,
      keep the input rooting unchanged;
    * otherwise, reroot the tree so that the reference becomes the top-level outgroup.
- The root is treated as the final alignment level, but it is not handled as an ordinary
  non-root internal vertex when deciding biological sub-tasks.
- At the top root, if the planned tree has the form (reference, X) and X is an
  internal node, X itself is not materialized as a separate task. X.children enter
  the root-level alignment directly, and X.name can still be preserved as the
  ingroup-side top label in generated guide trees.
- For non-root internal nodes:
    * exactly two direct children -> fixed single-run task
    * more than two direct children -> consensus task
- Fixed tasks do not add the reference unless the reference is already one of their true
  direct descendants in the fixed tree.
- Consensus tasks always include the global reference as an outgroup for
  reference-coordinate consensus extraction.
- RunPipelineUseThis.sh is called with --ifRefNonOutgroup 1 for every consensus
  task. The pipeline therefore emits one nested HAL of the form
  (reference,(ingroup)Anc_name_ofIngroup)Anc_name; the ingroup ancestor remains
  the task-level FASTA used by higher levels.

Notes
-----
- Consensus tasks with exactly 3 sampled elements enumerate all 3 binary topologies internally.
- Consensus tasks with >3 sampled elements call an external guide-tree generator:
      generate_random_guidetrees_2models_2modes_finalver.py
- The consensus-extraction / ancestor-inference step is written as a real
  RunPipelineUseThis.sh command in instruction.txt.

Assumptions
-----------
- Taxon names and internal names do not contain commas, semicolons, parentheses.
- Branch lengths and support values are ignored when parsing trees.
- The path file has two whitespace-separated columns:
      TAXON   PATH
- If a path is relative, it is resolved relative to the path file directory.
"""

from __future__ import annotations

import argparse
import copy
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


# =============================================================================
# 0. Tree data structure and Newick utilities
# =============================================================================

@dataclass
class Node:
    name: Optional[str] = None
    children: List["Node"] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def clone(self) -> "Node":
        return copy.deepcopy(self)


_TOKEN_RE = re.compile(r"\s*([(),;]|[^(),;:\s]+)(?::[^(),;]+)?\s*")


def tokenize_newick(s: str) -> List[str]:
    tokens = _TOKEN_RE.findall(s.strip())
    if not tokens or tokens[-1] != ";":
        raise ValueError("Newick string must end with ';'")
    return tokens


def parse_newick(s: str) -> Node:
    tokens = tokenize_newick(s)
    idx = 0

    def parse_subtree() -> Node:
        nonlocal idx
        if tokens[idx] == "(":
            idx += 1
            children: List[Node] = []
            while True:
                children.append(parse_subtree())
                if tokens[idx] == ",":
                    idx += 1
                    continue
                if tokens[idx] == ")":
                    idx += 1
                    break
                raise ValueError(f"Unexpected token near {tokens[idx]!r}")
            name = None
            if idx < len(tokens) and tokens[idx] not in [",", ")", ";"]:
                name = tokens[idx]
                idx += 1
            return Node(name=name, children=children)

        if tokens[idx] in [",", ")", ";"]:
            raise ValueError(f"Unexpected token near {tokens[idx]!r}")
        leaf = Node(name=tokens[idx], children=[])
        idx += 1
        return leaf

    root = parse_subtree()
    if tokens[idx] != ";":
        raise ValueError("Unexpected trailing tokens in Newick")
    return root


def to_newick(node: Node) -> str:
    if node.is_leaf():
        return node.name or ""
    return "(" + ",".join(to_newick(c) for c in node.children) + ")" + (node.name or "")


def collapse_unary_nodes(node: Node) -> Node:
    """
    Recursively remove redundant unary internal nodes.

    Example:
        (((A,B)X)Y)  ->  ((A,B)X)
    """
    if node.is_leaf():
        return node.clone()

    new_children = [collapse_unary_nodes(c) for c in node.children]
    new_node = Node(name=node.name, children=new_children)

    while len(new_node.children) == 1:
        new_node = new_node.children[0]

    return new_node


def _ascii_tree_lines(node: Node, reference: str, prefix: str = "", is_last: bool = True) -> List[str]:
    label = node.name or "Unnamed"
    if node.is_leaf() and label == reference:
        label += " [reference]"

    if prefix == "":
        lines = [label]
    else:
        connector = "└── " if is_last else "├── "
        lines = [prefix + connector + label]

    if node.children:
        child_prefix = prefix + ("    " if is_last else "│   ")
        for i, child in enumerate(node.children):
            child_is_last = (i == len(node.children) - 1)
            lines.extend(_ascii_tree_lines(child, reference, child_prefix, child_is_last))
    return lines


def ascii_tree_view(root: Node, reference: str) -> str:
    return "\n".join(_ascii_tree_lines(root, reference, prefix="", is_last=True))


def find_leaf(root: Node, leaf_name: str) -> Node:
    stack = [root]
    found: List[Node] = []
    while stack:
        n = stack.pop()
        if n.is_leaf() and n.name == leaf_name:
            found.append(n)
        stack.extend(n.children)
    if not found:
        raise ValueError(f"Reference leaf {leaf_name!r} not found in tree")
    if len(found) > 1:
        raise ValueError(f"Leaf name {leaf_name!r} appears multiple times")
    return found[0]



def reroot_at_reference(root: Node, reference: str) -> Node:
    """
    Reroot a tree so that the reference leaf is the top-level outgroup.

    The returned root has the form:

        (reference, X)<old_root_name_or_auto_root>

    where X is the non-reference side.  If the reference was already attached
    directly to a multifurcating root, the non-reference children are wrapped
    into a new unnamed ingroup holder.  That holder is auto-named later.  This
    deliberately produces the top-level (reference, X) structure expected by
    the planner, while the root-level task will still use X.children directly
    and will not materialize X as a separate task.
    """
    root = root.clone()
    ref_leaf = find_leaf(root, reference)

    id_to_node: Dict[int, Node] = {}
    adj: Dict[int, List[int]] = {}

    stack = [root]
    while stack:
        n = stack.pop()
        nid = id(n)
        id_to_node[nid] = n
        adj.setdefault(nid, [])
        for c in n.children:
            cid = id(c)
            id_to_node[cid] = c
            adj.setdefault(cid, [])
            adj[nid].append(cid)
            adj[cid].append(nid)
            stack.append(c)

    ref_id = id(ref_leaf)
    ref_neighbors = adj[ref_id]
    if len(ref_neighbors) != 1:
        raise ValueError(f"Reference leaf {reference!r} should have exactly one neighbor")

    old_root_id = id(root)
    ingroup_neighbor = ref_neighbors[0]

    def rebuild(current_id: int, parent_id: int) -> Node:
        current = id_to_node[current_id]
        downstream = [x for x in adj[current_id] if x != parent_id]
        if not downstream:
            return Node(name=current.name, children=[])
        return Node(
            name=current.name,
            children=[rebuild(nb, current_id) for nb in downstream],
        )

    # If the reference is already a direct child of the original top root,
    # avoid duplicating the old root name on both the new outer root and the
    # ingroup holder.  The holder is left unnamed and auto-named later.
    if ingroup_neighbor == old_root_id:
        ingroup_children: List[Node] = []
        for c in root.children:
            if c.is_leaf() and c.name == reference:
                continue
            ingroup_children.append(rebuild(id(c), old_root_id))

        if len(ingroup_children) == 1:
            ingroup_side = ingroup_children[0]
        else:
            ingroup_side = Node(name=None, children=ingroup_children)

        return Node(
            name=root.name,
            children=[Node(name=reference, children=[]), ingroup_side],
        )

    return Node(
        name=root.name,
        children=[
            Node(name=reference, children=[]),
            rebuild(ingroup_neighbor, ref_id),
        ],
    )


def node_contains_leaf(node: Node, leaf_name: str) -> bool:
    stack = [node]
    while stack:
        n = stack.pop()
        if n.is_leaf() and n.name == leaf_name:
            return True
        stack.extend(n.children)
    return False


def validate_existing_node_names(root: Node) -> None:
    leaf_counts: Dict[str, int] = {}
    internal_counts: Dict[str, int] = {}

    def walk(n: Node) -> None:
        if n.name:
            if n.is_leaf():
                leaf_counts[n.name] = leaf_counts.get(n.name, 0) + 1
            else:
                internal_counts[n.name] = internal_counts.get(n.name, 0) + 1
        for c in n.children:
            walk(c)

    walk(root)

    duplicate_leaf_names = sorted(name for name, count in leaf_counts.items() if count > 1)
    duplicate_internal_names = sorted(name for name, count in internal_counts.items() if count > 1)
    leaf_internal_overlaps = sorted(set(leaf_counts) & set(internal_counts))

    if duplicate_leaf_names or duplicate_internal_names or leaf_internal_overlaps:
        messages = ["Input tree contains duplicate named nodes."]
        if duplicate_leaf_names:
            messages.append(
                "Duplicate leaf taxon names: " + ", ".join(repr(x) for x in duplicate_leaf_names)
            )
        if duplicate_internal_names:
            messages.append(
                "Duplicate internal node names: " + ", ".join(repr(x) for x in duplicate_internal_names)
            )
        if leaf_internal_overlaps:
            messages.append(
                "Names used by both leaf and internal node: "
                + ", ".join(repr(x) for x in leaf_internal_overlaps)
            )
        raise ValueError("\n".join(messages))


def auto_name_internal_nodes(root: Node) -> None:
    used_names: Set[str] = set()

    def collect_used_names(n: Node) -> None:
        if n.name:
            used_names.add(n.name)
        for c in n.children:
            collect_used_names(c)

    collect_used_names(root)

    counter = 1

    def next_internal_name() -> str:
        nonlocal counter
        while True:
            candidate = f"Anc_in{counter}"
            counter += 1
            if candidate not in used_names:
                used_names.add(candidate)
                return candidate

    def next_root_name() -> str:
        if "Root" not in used_names:
            used_names.add("Root")
            return "Root"
        suffix = 1
        while True:
            candidate = f"Root_{suffix}"
            suffix += 1
            if candidate not in used_names:
                used_names.add(candidate)
                return candidate

    def postorder(n: Node, is_root: bool = False) -> None:
        for c in n.children:
            postorder(c, False)
        if not n.is_leaf():
            if is_root:
                if not n.name:
                    n.name = next_root_name()
            else:
                if not n.name:
                    n.name = next_internal_name()

    postorder(root, True)


# =============================================================================
# 1. Input path file
# =============================================================================

def read_taxon_paths(path_file: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    base = path_file.resolve().parent
    with open(path_file, "r", encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"Invalid path file line {lineno}: {line!r}")
            taxon = parts[0]
            p = Path(parts[1])
            if not p.is_absolute():
                p = (base / p).resolve()
            mapping[taxon] = str(p)
    return mapping


# =============================================================================
# 2. Task-level tree helpers
# =============================================================================

def root_has_direct_reference(root: Node, reference: str) -> bool:
    return any(c.is_leaf() and c.name == reference for c in root.children)


def is_fully_binary_tree(root: Node) -> bool:
    """
    Return True only when every internal node has exactly two children.

    The check is applied after collapsing unary nodes in main(), so it treats
    the input as a conventional rooted binary tree whenever all remaining
    internal vertices are bifurcating.
    """
    if root.is_leaf():
        return True
    if len(root.children) != 2:
        return False
    return all(is_fully_binary_tree(c) for c in root.children)


def direct_alignment_children(node: Node, reference: str, is_root: bool) -> Tuple[List[Node], Optional[Node]]:
    """
    Return the direct children that define the current alignment decision.

    For ordinary non-root nodes, these are simply node.children.

    After the conditional rerooting step, at the top root
    there are two distinct reference cases:

    1. Explicit root split: (reference, X)Root, where X is a single internal
       holder.  In this case the reference is treated as an external top-level
       reference/outgroup, X is not materialized as a separate task, and
       X.children are used for the root-level decision.

    2. Root polytomy with reference as one direct child:
       (reference, A, B, ...)Root.  In this case the reference is a true direct
       alignment child and participates in the root-level consensus sampling.
       Therefore node.children, including the reference, are used directly.

    Example:
        (GALGA,CATAU,(CHLUN,(COLLI,MESUN)))Root

    with GALGA as reference is a 3-way root-level consensus task, because the
    direct alignment children are GALGA, CATAU, and the internal CHLUN/COLLI/MESUN
    side.
    """
    if is_root:
        ref_children = [c for c in node.children if c.is_leaf() and c.name == reference]
        non_ref = [c for c in node.children if not (c.is_leaf() and c.name == reference)]

        # Only the explicit binary root split treats the reference as external
        # and flattens the single ingroup holder.
        if len(ref_children) == 1 and len(non_ref) == 1 and not non_ref[0].is_leaf():
            return non_ref[0].children, non_ref[0]

        # In a root polytomy, the reference remains a true direct child.
        return node.children, None

    return node.children, None

def is_consensus_target(node: Node, reference: str, is_root: bool) -> bool:
    children, _ = direct_alignment_children(node, reference, is_root)
    return len(children) > 2


def relabel_tree_internal_nodes(
    root: Node,
    top_name: Optional[str],
    internal_prefix: str,
    protected_names: Optional[Set[str]] = None,
) -> Node:
    root = root.clone()
    protected_names = protected_names or set()
    counter = 1

    def next_label() -> str:
        nonlocal counter
        while True:
            label = f"{internal_prefix}_{counter}"
            counter += 1
            if label not in protected_names:
                return label

    def postorder(n: Node, is_top: bool = False) -> None:
        for c in n.children:
            postorder(c, False)
        if not n.is_leaf():
            if n.name in protected_names:
                return
            if is_top and top_name is not None:
                n.name = top_name
            else:
                n.name = next_label()

    postorder(root, True)
    return root

def relabel_external_outgroup_tree(
    full_root: Node,
    reference: str,
    outer_root_name: str,
    ingroup_top_name: Optional[str],
    internal_prefix: str,
) -> Node:
    """
    Relabel a generated tree of the form (reference, ingroup_subtree).
    The outer root receives outer_root_name. The non-reference child receives
    ingroup_top_name if provided; otherwise it is named by internal_prefix.
    """
    full_root = full_root.clone()
    if full_root.is_leaf():
        raise ValueError("External-outgroup guide tree unexpectedly has a leaf root.")

    ref_children = [c for c in full_root.children if c.is_leaf() and c.name == reference]
    non_ref_children = [c for c in full_root.children if not (c.is_leaf() and c.name == reference)]
    if len(ref_children) != 1 or len(non_ref_children) != 1:
        raise ValueError(
            "Expected generated external-outgroup tree to have one reference child "
            "and one non-reference ingroup child."
        )

    full_root.name = outer_root_name
    relabeled_ingroup = relabel_tree_internal_nodes(
        non_ref_children[0], top_name=ingroup_top_name, internal_prefix=internal_prefix
    )
    full_root.children = [ref_children[0].clone(), relabeled_ingroup]
    return full_root


Shape = Union[str, Tuple["Shape", "Shape"]]


# =============================================================================
# 3. Parameters and task records
# =============================================================================

@dataclass
class GuideParams:
    generator: Path
    num_trees: Optional[int] = None
    model: Optional[str] = None
    rf_threshold: Optional[str] = None
    t_threshold: Optional[str] = None
    max_tries: Optional[int] = None
    restarts: Optional[int] = None
    seed: Optional[int] = None
    python_bin: str = sys.executable

    def command_for(self, taxa: List[str], outgroup: Optional[str] = None) -> List[str]:
        cmd = [self.python_bin, str(self.generator), "--taxa", *taxa]
        if outgroup is not None:
            cmd += ["--outgroup", outgroup]
        if self.num_trees is not None:
            cmd += ["--num_trees", str(self.num_trees)]
        if self.model is not None:
            cmd += ["--model", self.model]
        if self.rf_threshold is not None:
            cmd += ["--rf-threshold", str(self.rf_threshold)]
        if self.t_threshold is not None:
            cmd += ["--t-threshold", str(self.t_threshold)]
        if self.max_tries is not None:
            cmd += ["--max-tries", str(self.max_tries)]
        if self.restarts is not None:
            cmd += ["--restarts", str(self.restarts)]
        if self.seed is not None:
            cmd += ["--seed", str(self.seed)]
        return cmd


@dataclass
class PipelineParams:
    pipeline: Path
    threads: int
    common_workers: int
    model_file: Path
    maf_to_concat_fasta: Path
    sep_length: Optional[int] = None
    chrom_length_threshold: Optional[int] = None
    global_num: Optional[int] = None
    separate_workers: Optional[int] = None


@dataclass
class Task:
    name: str
    kind: str
    node: Node
    is_root: bool
    folder: Path
    parent_task: Optional["Task"] = None
    dependencies: List["Task"] = field(default_factory=list)
    alignment_files: List[Path] = field(default_factory=list)
    output_fa: Optional[Path] = None
    descendant_consensus_names: Set[str] = field(default_factory=set)
    genome_txt: Optional[Path] = None

    def __hash__(self) -> int:
        return hash((self.name, str(self.folder)))


# =============================================================================
# 4. Planner
# =============================================================================

class Planner:
    def __init__(
        self,
        root: Node,
        reference: str,
        taxon_paths: Dict[str, str],
        outdir: Path,
        guide_params: GuideParams,
        pipeline_params: PipelineParams,
    ):
        self.root = root
        self.reference = reference
        self.taxon_paths = taxon_paths
        self.outdir = outdir.resolve()
        self.guide_params = guide_params
        self.pipeline_params = pipeline_params
        self.tasks_by_name: Dict[str, Task] = {}
        self.root_task: Optional[Task] = None

    # -------------------------------------------------------------------------
    # Basic paths
    # -------------------------------------------------------------------------

    def task_folder(self, node: Node, parent_task: Optional[Task]) -> Path:
        if parent_task is None:
            return self.outdir / node.name
        return parent_task.folder / node.name

    def task_output_fa(self, task: Task) -> Path:
        return task.folder / f"{task.name}.fa"

    def task_genome_txt(self, task: Task) -> Path:
        return task.folder / f"{task.name}.genome.txt"

    def task_primary_hal(self, task: Task) -> Path:
        if task.kind == "consensus":
            # RunPipelineUseThis.sh is called with --pre <task.folder>/consensus.
            return task.folder / "consensus.hal"
        return task.folder / f"{task.name}_aln1.hal"

    def task_primary_maf(self, task: Task) -> Path:
        if task.kind == "consensus":
            return task.folder / "consensus.maf"
        return task.folder / f"{task.name}_aln1.maf"

    def task_primary_fasta(self, task: Task) -> Path:
        if task.kind == "consensus":
            return task.folder / "consensus.fasta"
        return task.output_fa if task.output_fa is not None else (task.folder / f"{task.name}.fa")

    def task_primary_stem(self, task: Task) -> str:
        return self.task_primary_hal(task).stem

    def root_alltaxa_script(self) -> Path:
        return self.outdir / "finalization.sh"

    def bundled_hal_append_subtree(self) -> Path:
        return self.pipeline_params.pipeline.resolve().parent / "tool_used" / "bin" / "halAppendSubtree"

    def root_alltaxa_hal(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.root_task.folder / "FinalResultAlignment.hal"

    def root_alltaxa_maf(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.root_task.folder / "FinalResultAlignment.maf"

    def root_alltaxa_fasta(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.root_task.folder / "FinalResultAlignment.fasta"

    def final_output_dir(self) -> Path:
        return self.outdir / "final_output"

    def final_output_hal(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.final_output_dir() / self.root_alltaxa_hal().name

    def final_output_maf(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.final_output_dir() / self.root_alltaxa_maf().name

    def final_output_fasta(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.final_output_dir() / self.root_alltaxa_fasta().name

    def final_output_ancestor_fasta(self) -> Path:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        return self.final_output_dir() / self.root_task.output_fa.name

    def render_run_path(self, p: Path) -> str:
        p = p.resolve()
        try:
            rel = p.relative_to(self.outdir)
        except ValueError:
            return shlex.quote(str(p))
        if str(rel) == ".":
            return "${RUN_PATH}"
        return f"${{RUN_PATH}}/{rel.as_posix()}"

    # -------------------------------------------------------------------------
    # Task creation and dependency planning
    # -------------------------------------------------------------------------

    def create_task(self, node: Node, parent_task: Optional[Task], is_root: bool, kind: str) -> Task:
        if node.name in self.tasks_by_name:
            t = self.tasks_by_name[node.name]
            if parent_task is not None and t not in parent_task.dependencies:
                parent_task.dependencies.append(t)
            return t

        t = Task(
            name=node.name,
            kind=kind,
            node=node,
            is_root=is_root,
            folder=self.task_folder(node, parent_task),
            parent_task=parent_task,
        )
        t.output_fa = self.task_output_fa(t)
        if kind == "consensus":
            t.genome_txt = self.task_genome_txt(t)
        self.tasks_by_name[node.name] = t
        if parent_task is not None and t not in parent_task.dependencies:
            parent_task.dependencies.append(t)
        return t

    def task_direct_children(self, task: Task) -> Tuple[List[Node], Optional[Node]]:
        return direct_alignment_children(task.node, self.reference, task.is_root)

    def reference_is_direct_alignment_child(self, task: Task) -> bool:
        children, _ = self.task_direct_children(task)
        return any(c.is_leaf() and c.name == self.reference for c in children)

    def child_containing_reference(self, task: Task) -> Optional[Node]:
        children, _ = self.task_direct_children(task)
        for c in children:
            if node_contains_leaf(c, self.reference):
                return c
        return None

    def reference_role_for_consensus(self, task: Task) -> Tuple[str, Optional[Node]]:
        """
        Return one of:
            direct   : reference is a true direct alignment child of this task.
            nested   : reference is inside one direct child of this task.
            external : reference is outside this task and is added only as an external reference.
        """
        if self.reference_is_direct_alignment_child(task):
            return "direct", None
        child = self.child_containing_reference(task)
        if child is not None:
            return "nested", child
        return "external", None

    def consensus_reference_is_technical_only(self, task: Task) -> bool:
        return task.kind == "consensus" and not self.reference_is_direct_alignment_child(task)

    def ensure_task_for_materialized_node(self, node: Node, parent_task: Task, is_root: bool = False) -> Task:
        kind = "consensus" if is_consensus_target(node, self.reference, is_root) else "fixed"
        task = self.create_task(node, parent_task, is_root, kind)

        if kind == "consensus":
            children, _ = direct_alignment_children(node, self.reference, is_root)
            for u in children:
                if not u.is_leaf():
                    self.ensure_task_for_materialized_node(u, task, False)
        else:
            self.collect_consensus_descendants_for_fixed(node, task)
        return task

    def collect_consensus_descendants_for_fixed(self, node: Node, current_task: Task) -> None:
        def walk(n: Node) -> None:
            if n is not node and (not n.is_leaf()) and is_consensus_target(n, self.reference, False):
                current_task.descendant_consensus_names.add(n.name)
                self.ensure_task_for_materialized_node(n, current_task, False)
                return
            for c in n.children:
                walk(c)

        for c in node.children:
            walk(c)

    def plan(self) -> Task:
        root_kind = "consensus" if is_consensus_target(self.root, self.reference, True) else "fixed"
        self.root_task = self.create_task(self.root, None, True, root_kind)

        if root_kind == "consensus":
            children, _ = direct_alignment_children(self.root, self.reference, True)
            for u in children:
                if not u.is_leaf():
                    self.ensure_task_for_materialized_node(u, self.root_task, False)
        else:
            self.collect_consensus_descendants_for_fixed(self.root, self.root_task)

        return self.root_task

    # -------------------------------------------------------------------------
    # Guide-tree strings
    # -------------------------------------------------------------------------

    def unit_label(self, unit: Node) -> str:
        if not unit.name:
            raise ValueError("Encountered unnamed direct alignment child after auto-naming.")
        return unit.name

    def rec_fixed_subtree(self, n: Node, task: Task) -> str:
        materialized = set(task.descendant_consensus_names)
        if n.is_leaf():
            return n.name
        if n.name in materialized and n.name != task.name:
            return n.name
        inside = ",".join(self.rec_fixed_subtree(c, task) for c in n.children)
        return f"({inside}){n.name}"

    def fixed_tree_string(self, task: Task) -> str:
        # Root special case only for an explicit root split: (reference, X)Root,
        # where X is an internal holder. In that case X is not materialized as a
        # separate task, but the root-level fixed tree keeps reference outside
        # the holder-side subtree.
        if task.is_root and self.reference_is_direct_alignment_child(task):
            children, holder = self.task_direct_children(task)
            if holder is not None:
                if len(children) == 1:
                    return f"({self.reference},{self.rec_fixed_subtree(children[0], task)}){task.name};"
                if len(children) == 2:
                    left = self.rec_fixed_subtree(children[0], task)
                    right = self.rec_fixed_subtree(children[1], task)
                    ingroup_name = holder.name if holder.name else f"{task.name}_ingroup"
                    return f"({self.reference},({left},{right}){ingroup_name}){task.name};"

        return self.rec_fixed_subtree(task.node, task) + ";"


    def topology_shapes_for_three(self, labels: List[str]) -> List[Shape]:
        if len(labels) != 3:
            raise ValueError("topology_shapes_for_three requires exactly three labels")
        a, b, c = labels
        return [((a, b), c), ((a, c), b), ((b, c), a)]

    def build_named_tree_from_shape(self, shape: Shape, top_name: Optional[str], internal_prefix: str) -> str:
        counter = 1

        def rec(x: Shape, is_top: bool = False) -> str:
            nonlocal counter
            if isinstance(x, str):
                return x
            left = rec(x[0], False)
            right = rec(x[1], False)
            if is_top and top_name is not None:
                label = top_name
            else:
                label = f"{internal_prefix}_{counter}"
                counter += 1
            return f"({left},{right}){label}"

        return rec(shape, True)

    def parse_generator_newicks(self, cmd: List[str]) -> List[str]:
        try:
            res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.strip()
            stdout = e.stdout.strip()
            msg = ["Guide-tree generator failed.", "Command: " + " ".join(shlex.quote(x) for x in cmd)]
            if stdout:
                msg.append("STDOUT:\n" + stdout)
            if stderr:
                msg.append("STDERR:\n" + stderr)
            raise RuntimeError("\n".join(msg)) from e

        newick_lines: List[str] = []
        for raw in res.stdout.splitlines():
            s = raw.strip()
            if not s:
                continue
            if s.endswith(";"):
                newick_lines.append(s)

        if not newick_lines:
            raise RuntimeError(
                "Guide-tree generator returned no Newick trees. "
                f"Command was: {' '.join(shlex.quote(x) for x in cmd)}"
            )
        return newick_lines
    
    def generated_trees_no_outgroup(
        self,
        labels: List[str],
        top_name: str,
        internal_prefix: str,
        protected_names: Optional[Set[str]] = None,
    ) -> List[str]:
        cmd = self.guide_params.command_for(labels, outgroup=None)
        newick_lines = self.parse_generator_newicks(cmd)
        out: List[str] = []
        for s in newick_lines:
            parsed = parse_newick(s)
            relabeled = relabel_tree_internal_nodes(
                parsed,
                top_name=top_name,
                internal_prefix=internal_prefix,
                protected_names=protected_names,
            )
            out.append(to_newick(relabeled))
        return out    

    def generated_trees_with_external_outgroup(
        self,
        labels: List[str],
        outgroup: str,
        outer_root_name: str,
        ingroup_top_name: Optional[str],
        internal_prefix: str,
    ) -> List[str]:
        cmd = self.guide_params.command_for(labels, outgroup=outgroup)
        newick_lines = self.parse_generator_newicks(cmd)
        out: List[str] = []
        for s in newick_lines:
            parsed = parse_newick(s)
            relabeled = relabel_external_outgroup_tree(
                parsed,
                reference=outgroup,
                outer_root_name=outer_root_name,
                ingroup_top_name=ingroup_top_name,
                internal_prefix=internal_prefix,
            )
            out.append(to_newick(relabeled))
        return out

    def consensus_tree_strings(self, task: Task) -> List[str]:
        children, holder = self.task_direct_children(task)
        labels = [self.unit_label(u) for u in children]
        role, _ = self.reference_role_for_consensus(task)

        # Case 1: reference is a true direct child of the current polytomy.
        # It must participate in guide-tree sampling and must not be forced to
        # the external outgroup position.
        if role == "direct":
            top_name = task.name
            internal_prefix = task.name

            if len(labels) == 3:
                return [
                    self.build_named_tree_from_shape(
                        shape,
                        top_name=top_name,
                        internal_prefix=internal_prefix,
                    ) + ";"
                    for shape in self.topology_shapes_for_three(labels)
                ]

            return [s + ";" for s in self.generated_trees_no_outgroup(
                labels=labels,
                top_name=top_name,
                internal_prefix=internal_prefix,
                protected_names={self.reference},
            )]

        # Case 2: explicit root split or ordinary external-reference consensus.
        # The reference is a technical coordinate anchor / external outgroup.
        # It is added outside the sampled ingroup tree.
        if task.is_root and holder is not None and holder.name:
            ingroup_top_name = holder.name
            internal_prefix = holder.name
        elif task.is_root:
            ingroup_top_name = f"{task.name}_ingroup"
            internal_prefix = ingroup_top_name
        else:
            ingroup_top_name = task.name
            internal_prefix = task.name

        outer_root_name = task.name if task.is_root else f"{task.name}_{self.reference}"

        if len(labels) == 3:
            ingroups = [
                self.build_named_tree_from_shape(
                    shape,
                    top_name=ingroup_top_name,
                    internal_prefix=internal_prefix,
                )
                for shape in self.topology_shapes_for_three(labels)
            ]
            return [f"({self.reference},{ingroup}){outer_root_name};" for ingroup in ingroups]

        return [s + ";" for s in self.generated_trees_with_external_outgroup(
            labels=labels,
            outgroup=self.reference,
            outer_root_name=outer_root_name,
            ingroup_top_name=ingroup_top_name,
            internal_prefix=internal_prefix,
        )]

    # -------------------------------------------------------------------------
    # Alignment text files
    # -------------------------------------------------------------------------

    def leaf_paths_for_task(self, task: Task) -> Dict[str, str]:
        mapping: Dict[str, str] = {}

        def add_taxon_or_ancestor(label: str) -> None:
            if label in self.taxon_paths:
                mapping[label] = str(Path(self.taxon_paths[label]).resolve())
            elif label in self.tasks_by_name:
                dep_task = self.tasks_by_name[label]
                mapping[label] = str(dep_task.output_fa.resolve())
            else:
                raise ValueError(
                    f"Cannot find path for label {label!r}. It is neither a taxon in the path file "
                    "nor a planned ancestor task."
                )

        if task.kind == "consensus":
            # Consensus extraction always needs the reference as a coordinate anchor.
            add_taxon_or_ancestor(self.reference)
            children, _ = self.task_direct_children(task)
            for u in children:
                add_taxon_or_ancestor(self.unit_label(u))
        else:
            materialized = set(task.descendant_consensus_names)

            def rec(n: Node) -> None:
                if not n.is_leaf() and n.name in materialized and n.name != task.name:
                    add_taxon_or_ancestor(n.name)
                    return
                if n.is_leaf():
                    add_taxon_or_ancestor(n.name)
                    return
                for c in n.children:
                    rec(c)

            rec(task.node)

        return mapping

    def write_task_files(self, task: Task) -> None:
        task.folder.mkdir(parents=True, exist_ok=True)

        trees = self.consensus_tree_strings(task) if task.kind == "consensus" else [self.fixed_tree_string(task)]
        leaf_path_map = self.leaf_paths_for_task(task)
        task.alignment_files = []

        for i, tree_str in enumerate(trees, 1):
            aln_txt = task.folder / f"{task.name}_aln{i}.txt"
            task.alignment_files.append(aln_txt)
            with open(aln_txt, "w", encoding="utf-8") as fw:
                fw.write(tree_str.rstrip() + "\n")
                for label in sorted(leaf_path_map):
                    fw.write(f"{label} {leaf_path_map[label]}\n")

        if task.kind == "consensus" and task.genome_txt is not None:
            with open(task.genome_txt, "w", encoding="utf-8") as fw:
                for label in sorted(leaf_path_map):
                    fw.write(f"{label} {leaf_path_map[label]}\n")

    def write_all_task_files(self) -> None:
        for task in self.topologically_sorted_tasks():
            self.write_task_files(task)

    # -------------------------------------------------------------------------
    # Dependency sorting
    # -------------------------------------------------------------------------

    def task_dependency_graph(self) -> Dict[str, Set[str]]:
        graph: Dict[str, Set[str]] = {}
        for name, task in self.tasks_by_name.items():
            graph.setdefault(name, set())
            for dep in task.dependencies:
                graph[name].add(dep.name)
        return graph

    def topologically_sorted_task_levels(self) -> List[List[Task]]:
        graph = self.task_dependency_graph()
        indeg: Dict[str, int] = {n: len(deps) for n, deps in graph.items()}
        reverse: Dict[str, Set[str]] = {n: set() for n in graph}
        for n, deps in graph.items():
            for d in deps:
                reverse.setdefault(d, set()).add(n)

        levels: List[List[Task]] = []
        ready = sorted([n for n, deg in indeg.items() if deg == 0])
        seen: Set[str] = set()

        while ready:
            this_level_names = ready[:]
            ready = []
            level_tasks = [self.tasks_by_name[n] for n in this_level_names]
            levels.append(level_tasks)

            for n in this_level_names:
                seen.add(n)
                for m in sorted(reverse.get(n, [])):
                    indeg[m] -= 1
                    if indeg[m] == 0:
                        ready.append(m)

        if len(seen) != len(graph):
            raise RuntimeError("Task graph contains a cycle, which should not happen.")
        return levels

    def topologically_sorted_tasks(self) -> List[Task]:
        out: List[Task] = []
        for level in self.topologically_sorted_task_levels():
            out.extend(level)
        return out

    # -------------------------------------------------------------------------
    # Regrafting and final output
    # -------------------------------------------------------------------------

    def all_leaf_taxa(self) -> List[str]:
        taxa: List[str] = []

        def walk(n: Node) -> None:
            if n.is_leaf():
                if n.name not in taxa:
                    taxa.append(n.name)
                return
            for c in n.children:
                walk(c)

        walk(self.root)
        non_ref = sorted([x for x in taxa if x != self.reference])
        return [self.reference] + non_ref

    def regraft_order_within_task(self, task: Task) -> List[Task]:
        out: List[Task] = []
        seen: Set[str] = set()

        def append_task(child_task: Task) -> None:
            if child_task.name in seen:
                return
            seen.add(child_task.name)
            out.append(child_task)
            for deeper in self.regraft_order_within_task(child_task):
                if deeper.name not in seen:
                    seen.add(deeper.name)
                    out.append(deeper)

        if task.kind == "consensus":
            children, _ = self.task_direct_children(task)
            for u in children:
                if not u.is_leaf() and u.name in self.tasks_by_name:
                    append_task(self.tasks_by_name[u.name])
        else:
            materialized = set(task.descendant_consensus_names)

            def walk_fixed(n: Node) -> None:
                if not n.is_leaf() and n.name in materialized and n.name != task.name and n.name in self.tasks_by_name:
                    append_task(self.tasks_by_name[n.name])
                    return
                for c in n.children:
                    walk_fixed(c)

            for c in task.node.children:
                walk_fixed(c)

        return out

    def final_regraft_tasks(self) -> List[Task]:
        if self.root_task is None:
            return []
        return self.regraft_order_within_task(self.root_task)

    def regraft_cmd_for_script(self, destination_hal: Path, task: Task) -> str:
        return (
            f"\"$patched_halAppendSubtree\" {self.render_run_path(destination_hal)} "
            f"{self.render_run_path(self.task_primary_hal(task))} {shlex.quote(task.name)} {shlex.quote(task.name)} --merge"
        )

    def root_alltaxa_hal2maf_cmd(self) -> str:
        if self.root_task is None:
            raise RuntimeError("Root task has not been planned yet")
        jobstore = self.root_task.folder / f"jobstorehal2maf_{self.root_task.name}_allTaxa"
        return (
            f"cactus-hal2maf --dupeMode single --chunkSize 500000 "
            f"--refGenome {shlex.quote(self.reference)} --noAncestors "
            f"{self.render_run_path(jobstore)} {self.render_run_path(self.root_alltaxa_hal())} "
            f"{self.render_run_path(self.root_alltaxa_maf())}"
        )

    def root_alltaxa_fasta_cmd(self) -> str:
        taxa = self.all_leaf_taxa()
        taxa_args = " ".join(shlex.quote(t) for t in taxa)
        script = self.render_run_path(self.pipeline_params.maf_to_concat_fasta)
        return (
            f"taxa_args=({taxa_args})\n"
            f"python {script} \"${{taxa_args[@]}}\" < {self.render_run_path(self.root_alltaxa_maf())} > "
            f"{self.render_run_path(self.root_alltaxa_fasta())}"
        )

    def write_finalization_script(self) -> Optional[Path]:
        """
        Write the finalization script.

        This script is needed even when there are no descendant HALs to regraft,
        because a fully fixed/binary root alignment still needs to be exported
        from HAL to the final MAF and concatenated FASTA outputs.
        """
        regraft_tasks = self.final_regraft_tasks()
        script_path = self.root_alltaxa_script()

        if self.root_task is None:
            if script_path.exists():
                script_path.unlink()
            return None

        bundled_hal_append_subtree = self.bundled_hal_append_subtree().resolve()

        with open(script_path, "w", encoding="utf-8") as fw:
            fw.write("#!/bin/bash\n")
            fw.write("set -euo pipefail\n\n")
            fw.write("RUN_PATH=" + shlex.quote(str(self.outdir)) + "\n")

            if regraft_tasks:
                fw.write("patched_halAppendSubtree=" + shlex.quote(str(bundled_hal_append_subtree)) + "\n\n")
                fw.write('for exe in "$patched_halAppendSubtree"; do\n')
                fw.write('    if [[ ! -f "$exe" ]]; then\n')
                fw.write('        echo "Error: bundled executable not found: $exe"\n')
                fw.write("        exit 1\n")
                fw.write("    fi\n")
                fw.write('    chmod +x "$exe"\n')
                fw.write('    if [[ ! -x "$exe" ]]; then\n')
                fw.write('        echo "Error: bundled executable is still not executable: $exe"\n')
                fw.write("        exit 1\n")
                fw.write("    fi\n")
                fw.write("done\n\n")

            fw.write(f'root_alltaxa_hal="{self.render_run_path(self.root_alltaxa_hal())}"\n')
            fw.write('# Copy the root-level primary HAL to a stable final-working name.\n')
            fw.write(f'cp {self.render_run_path(self.task_primary_hal(self.root_task))} "$root_alltaxa_hal"\n')

            for task in regraft_tasks:
                fw.write(self.regraft_cmd_for_script(self.root_alltaxa_hal(), task) + "\n")

            fw.write(self.root_alltaxa_hal2maf_cmd() + "\n")
            fw.write(self.root_alltaxa_fasta_cmd() + "\n")
            fw.write(f'mkdir -p {self.render_run_path(self.final_output_dir())}\n')
            fw.write(f'mv "$root_alltaxa_hal" {self.render_run_path(self.final_output_hal())}\n')
            fw.write(f'mv {self.render_run_path(self.root_alltaxa_maf())} {self.render_run_path(self.final_output_maf())}\n')
            fw.write(f'mv {self.render_run_path(self.root_alltaxa_fasta())} {self.render_run_path(self.final_output_fasta())}\n')
            fw.write(f'mv {self.render_run_path(self.root_task.output_fa)} {self.render_run_path(self.final_output_ancestor_fasta())}\n')

        script_path.chmod(0o755)
        return script_path

    def write_final_output_summary(self, fw) -> None:
        fw.write("\n")
        fw.write("############################################################\n")
        fw.write("# Final output summary\n")
        fw.write(f"# RUN_PATH = {self.outdir}\n")

        if self.root_task is None:
            fw.write("# Root task was not created.\n")
            fw.write("############################################################\n")
            return

        final_hal = self.final_output_hal().resolve()
        final_maf = self.final_output_maf().resolve()
        final_fasta = self.final_output_fasta().resolve()
        final_ancestor_fasta = self.final_output_ancestor_fasta().resolve()

        fw.write(f"# Final result alignment in HAL   : {final_hal}\n")
        fw.write(f"# Final result alignment in MAF   : {final_maf}\n")
        fw.write(f"# Final result alignment in FASTA : {final_fasta}\n")
        fw.write(f"# Final inferred genome FASTA for {self.root_task.name}: {final_ancestor_fasta}\n")
        fw.write("############################################################\n")

    # -------------------------------------------------------------------------
    # Commands
    # -------------------------------------------------------------------------

    def processed_tree_newick(self) -> str:
        return to_newick(self.root) + ";"

    def processed_tree_ascii(self) -> str:
        return ascii_tree_view(self.root, self.reference)

    def cactus_cmd_for_fixed(self, task: Task) -> str:
        aln_txt = task.alignment_files[0]
        hal = task.folder / f"{task.name}_aln1.hal"
        jobstore = task.folder / f"jobstorecluster_{task.name}_aln1"
        cmd1 = (
            f"cactus --root {shlex.quote(task.name)} "
            f"{self.render_run_path(jobstore)} {self.render_run_path(aln_txt)} {self.render_run_path(hal)}"
        )
        cmd2 = (
            f"hal2fasta {self.render_run_path(hal)} {shlex.quote(task.name)} > "
            f"{self.render_run_path(task.output_fa)}"
        )
        return cmd1 + "\n" + cmd2

    def cactus_cmds_for_consensus_alignments(self, task: Task) -> List[str]:
        cmds: List[str] = []
        for i, aln_txt in enumerate(task.alignment_files, 1):
            hal = task.folder / f"{task.name}_aln{i}.hal"
            maf = task.folder / f"{task.name}_aln{i}.maf"
            jobstore = task.folder / f"jobstorecluster_{task.name}_aln{i}"
            jobstore_hal2maf = task.folder / f"jobstorehal2maf_{task.name}_aln{i}"
            cmd1 = f"cactus {self.render_run_path(jobstore)} {self.render_run_path(aln_txt)} {self.render_run_path(hal)}"
            cmd2 = (
                f"cactus-hal2maf --dupeMode single --chunkSize 500000 "
                f"--refGenome {shlex.quote(self.reference)} --noAncestors "
                f"{self.render_run_path(jobstore_hal2maf)} {self.render_run_path(hal)} {self.render_run_path(maf)}"
            )
            cmds.append(cmd1 + "\n" + cmd2)
        return cmds

    def consensus_pipeline_cmd(self, task: Task) -> str:
        maf_paths = [
            self.render_run_path(task.folder / f"{task.name}_aln{i}.maf")
            for i in range(1, len(task.alignment_files) + 1)
        ]
        pre = self.render_run_path(task.folder / "consensus")
        genome_txt = self.render_run_path(task.genome_txt) if task.genome_txt is not None else ""
        pipeline = self.render_run_path(self.pipeline_params.pipeline)
        model_file = self.render_run_path(self.pipeline_params.model_file)

        _, holder = self.task_direct_children(task)
        role, _ = self.reference_role_for_consensus(task)

        # --ifRefNonOutgroup convention used by RunPipelineUseThis.sh:
        #   0 : reference is a true direct child of the consensus task and is
        #       allowed to move among the sampled guide trees.
        #   1 : reference is an external/top-level outgroup used as the
        #       coordinate anchor; --Anc_name_ofIngroup names the ingroup
        #       ancestor used by higher levels.
        if role == "direct":
            anc_name_for_pipeline = task.name
            extra_pipeline_options: List[Tuple[str, str]] = [
                ("--ifRefNonOutgroup", "0"),
            ]
        else:
            if task.is_root:
                anc_name_for_pipeline = task.name
                ingroup_name = holder.name if holder is not None and holder.name else f"{task.name}_ingroup"
            else:
                anc_name_for_pipeline = f"{task.name}_{self.reference}"
                ingroup_name = task.name

            extra_pipeline_options = [
                ("--ifRefNonOutgroup", "1"),
                ("--Anc_name_ofIngroup", ingroup_name),
            ]

        option_lines = [
            ("--input", " ".join(maf_paths)),
            ("--reference", shlex.quote(self.reference)),
            ("--pre", pre),
            ("--threads", str(self.pipeline_params.threads)),
            ("--common_workers", str(self.pipeline_params.common_workers)),
            ("--genome_path", genome_txt),
            ("--Anc_name", shlex.quote(anc_name_for_pipeline)),
            ("--ModelFile", model_file),
            ("--ifDeleteImmdiFiles", "0"),
        ]
        option_lines.extend(extra_pipeline_options)

        if self.pipeline_params.sep_length is not None:
            option_lines.append(("--sep_length", str(self.pipeline_params.sep_length)))
        if self.pipeline_params.chrom_length_threshold is not None:
            option_lines.append(("--chrom_length_threshold", str(self.pipeline_params.chrom_length_threshold)))
        if self.pipeline_params.global_num is not None:
            option_lines.append(("--global_num", str(self.pipeline_params.global_num)))
        if self.pipeline_params.separate_workers is not None:
            option_lines.append(("--separate_workers", str(self.pipeline_params.separate_workers)))

        lines = [f"bash {pipeline} \\"]
        for i, (opt, val) in enumerate(option_lines):
            suffix = " \\" if i < len(option_lines) - 1 else ""
            lines.append(f"  {opt} {val}{suffix}")

        return "\n".join(lines)


    def write_instruction(self) -> Path:
        instruction = self.outdir / "instruction.txt"
        levels = self.topologically_sorted_task_levels()
        cmd_no = 1
        regraft_tasks = self.final_regraft_tasks()

        with open(instruction, "w", encoding="utf-8") as fw:
            fw.write("## Processed tree used for planning (conditionally rerooted; auto-named internal nodes if needed):\n")
            fw.write(self.processed_tree_newick() + "\n\n")
            fw.write("## Tree view:\n")
            fw.write(self.processed_tree_ascii() + "\n\n")
            fw.write("## Copy and run all lines below that do not start with ##\n")
            fw.write("## Command0\n")
            fw.write("RUN_PATH=" + shlex.quote(str(self.outdir)) + "\n\n")
            fw.write("## Run the following commands:\n")
            fw.write("## Commands in the same block can be run in parallel if they are marked so.\n")
            fw.write("## You may need to adjust Cactus resource parameters yourself, e.g. --defaultMemory --maxCores --maxMemory.\n\n")

            for level_idx, level in enumerate(levels, 1):
                fixed_tasks = [t for t in level if t.kind == "fixed"]
                consensus_tasks = [t for t in level if t.kind == "consensus"]

                n_parallel = len(fixed_tasks) + sum(len(t.alignment_files) for t in consensus_tasks)
                if n_parallel > 1:
                    fw.write(f"## Level {level_idx}: Command{cmd_no}~Command{cmd_no + n_parallel - 1} can be run in parallel\n")

                for task in fixed_tasks:
                    fw.write(f"## Command{cmd_no}  [{task.name}; fixed single-run]\n")
                    fw.write(self.cactus_cmd_for_fixed(task) + "\n\n")
                    cmd_no += 1

                consensus_align_ranges: List[Tuple[Task, int, int]] = []
                for task in consensus_tasks:
                    start_here = cmd_no
                    cmds = self.cactus_cmds_for_consensus_alignments(task)
                    for c in cmds:
                        fw.write(f"## Command{cmd_no}  [{task.name}; guide-tree-specific alignment]\n")
                        fw.write(c + "\n\n")
                        cmd_no += 1
                    end_here = cmd_no - 1
                    consensus_align_ranges.append((task, start_here, end_here))

                for task, start_here, end_here in consensus_align_ranges:
                    if start_here == end_here:
                        fw.write(f"## After Command{start_here} is done, run Command{cmd_no}\n")
                    else:
                        fw.write(f"## After Command{start_here}~Command{end_here} are done, run Command{cmd_no}\n")
                    fw.write(f"## Command{cmd_no}  [{task.name}; extract consensus and infer ancestor]\n")
                    fw.write(self.consensus_pipeline_cmd(task) + "\n\n")
                    cmd_no += 1

            if self.root_task is not None:
                if regraft_tasks:
                    fw.write(f"## After Command{cmd_no - 1} is done, run descendant-HAL regrafting and final alignment export:\n")
                else:
                    fw.write(f"## After Command{cmd_no - 1} is done, run final alignment export:\n")
                fw.write("bash ${RUN_PATH}/finalization.sh\n\n")

            fw.write("## The planning step is done. The final root-level workspace is under:\n")
            if self.root_task is not None:
                fw.write(f"## {self.render_run_path(self.root_task.folder)}\n")

            self.write_final_output_summary(fw)

        return instruction


# =============================================================================
# 5. CLI
# =============================================================================

def build_arg_parser() -> argparse.ArgumentParser:
    script_dir = Path(__file__).resolve().parent
    default_generator = script_dir / "generate_random_guidetrees_2models_2modes_finalver.py"
    default_pipeline = script_dir / "RunPipelineUseThis.sh"
    default_model_file = script_dir / "tryMLstartree.mod"
    default_maf_to_concat_fasta = script_dir / "maf_to_concat_fasta.py"

    ap = argparse.ArgumentParser(description="Plan hierarchical Cactus runs from a partially resolved tree.")
    ap.add_argument("--tree", required=True, help="Input Newick tree string. It is rerooted to the reference unless the reference is already a direct child of the input top root.")
    ap.add_argument("--reference", required=True, help="Reference / outgroup taxon name used for reference-coordinate consensus extraction")
    ap.add_argument("--paths", required=True, type=Path, help="Two-column file: taxon path")
    ap.add_argument("--outdir", default=".", type=Path, help="Output directory (default: current working directory)")

    ap.add_argument("--generator", default=str(default_generator), help="Path to generate_random_guidetrees_2models_2modes_finalver.py")
    ap.add_argument("--guide-num-trees", type=int, default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-model", choices=["uniform", "yule"], default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-rf-threshold", default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-t-threshold", default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-max-tries", type=int, default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-restarts", type=int, default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--guide-seed", type=int, default=None, help="Pass through to the guide-tree generator")
    ap.add_argument("--generator-python", default=sys.executable, help="Python executable used to run the guide-tree generator")

    ap.add_argument("--pipeline", default=str(default_pipeline), help="Path to RunPipelineUseThis.sh")
    ap.add_argument("--maf_to_concat_fasta", default=str(default_maf_to_concat_fasta), help="Path to maf_to_concat_fasta.py")
    ap.add_argument("--threads", required=True, type=int, help="Pass through to RunPipelineUseThis.sh --threads")
    ap.add_argument("--common_workers", required=True, type=int, help="Pass through to RunPipelineUseThis.sh --common_workers")
    ap.add_argument("--ModelFile", default=str(default_model_file), type=Path, help="Path to the base model file used by RunPipelineUseThis.sh; default: tryMLstartree.mod in the same directory as this planner")
    ap.add_argument("--sep_length", type=int, default=None, help="Pass through to RunPipelineUseThis.sh")
    ap.add_argument("--chrom_length_threshold", type=int, default=None, help="Pass through to RunPipelineUseThis.sh")
    ap.add_argument("--global_num", type=int, default=None, help="Pass through to RunPipelineUseThis.sh")
    ap.add_argument("--separate_workers", type=int, default=None, help="Pass through to RunPipelineUseThis.sh")
    return ap


def main() -> None:
    ap = build_arg_parser()
    args = ap.parse_args()

    taxon_paths = read_taxon_paths(args.paths)
    if args.reference not in taxon_paths:
        raise ValueError(f"Reference {args.reference!r} is not present in the path file")

    raw_root = parse_newick(args.tree)
    validate_existing_node_names(raw_root)
    find_leaf(raw_root, args.reference)

    # Conditional rerooting rule:
    #   - If the input tree is already fully binary, keep the input rooting.
    #     A fully resolved tree should be treated as a fixed-task hierarchy,
    #     not rewritten by moving the reference to the top.
    #   - Else, if the reference is already a direct child of the input top root,
    #     keep the input rooting.  This preserves cases such as
    #         (GALGA,CATAU,(CHLUN,(COLLI,MESUN)));
    #     as a true root-level 3-way consensus task with GALGA participating.
    #   - Otherwise, reroot to the reference, which is the original planner
    #     behavior for non-binary partially resolved trees whose top root does
    #     not already contain the reference as a direct child.
    input_root = collapse_unary_nodes(raw_root)
    if is_fully_binary_tree(input_root):
        root = input_root
    elif root_has_direct_reference(input_root, args.reference):
        root = input_root
    else:
        root = collapse_unary_nodes(reroot_at_reference(input_root, args.reference))

    auto_name_internal_nodes(root)

    generator = Path(args.generator).resolve()
    if not generator.is_file():
        raise FileNotFoundError(f"Guide-tree generator script not found: {generator}")

    pipeline = Path(args.pipeline).resolve()
    if not pipeline.is_file():
        raise FileNotFoundError(f"Consensus pipeline script not found: {pipeline}")

    model_file = Path(args.ModelFile).resolve()
    if not model_file.is_file():
        raise FileNotFoundError(f"Model file not found: {model_file}")

    maf_to_concat_fasta = Path(args.maf_to_concat_fasta).resolve()
    if not maf_to_concat_fasta.is_file():
        raise FileNotFoundError(f"maf_to_concat_fasta.py not found: {maf_to_concat_fasta}")

    guide_params = GuideParams(
        generator=generator,
        num_trees=args.guide_num_trees,
        model=args.guide_model,
        rf_threshold=args.guide_rf_threshold,
        t_threshold=args.guide_t_threshold,
        max_tries=args.guide_max_tries,
        restarts=args.guide_restarts,
        seed=args.guide_seed,
        python_bin=args.generator_python,
    )

    pipeline_params = PipelineParams(
        pipeline=pipeline,
        threads=args.threads,
        common_workers=args.common_workers,
        model_file=model_file,
        maf_to_concat_fasta=maf_to_concat_fasta,
        sep_length=args.sep_length,
        chrom_length_threshold=args.chrom_length_threshold,
        global_num=args.global_num,
        separate_workers=args.separate_workers,
    )

    planner = Planner(
        root=root,
        reference=args.reference,
        taxon_paths=taxon_paths,
        outdir=args.outdir,
        guide_params=guide_params,
        pipeline_params=pipeline_params,
    )
    planner.plan()
    planner.write_all_task_files()
    planner.write_finalization_script()
    instruction = planner.write_instruction()

    print("Done.")
    print(f"Instruction file: {instruction}")
    if planner.root_task is not None:
        print(f"Root task folder: {planner.root_task.folder}")


if __name__ == "__main__":
    main()
