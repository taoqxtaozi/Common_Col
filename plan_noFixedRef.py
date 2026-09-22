#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plan_noFixedRef.py

Planner for hierarchical Cactus alignment of a rooted, partially resolved tree.

Core conventions
----------------
1. The input rooting is preserved. There is no global reference and no
   reference-driven rerooting.
2. Internal nodes with exactly two direct children are fixed single-run tasks.
3. Internal nodes with more than two direct children are consensus tasks.
4. Guide-tree generation and Cactus alignment do not use a reference.
5. For each consensus task, the local reference is selected at run time from
   the genomes directly participating in that task:
       - original taxon FASTA files, and/or
       - ancestor FASTA files produced by lower-level tasks.
   The genome with the greatest total FASTA sequence length is selected by
   select_longest_reference.py.
6. The selected local reference is used only for:
       - cactus-hal2maf --refGenome
       - RunPipelineUseThis.sh --reference
   and every consensus task uses --ifRefNonOutgroup 0.
7. The root-level final HAL -> MAF export also requires a coordinate reference.
   That export-only reference is supplied explicitly with --final_reference.
   It is used only for the final full-alignment export and does not affect
   planning, guide trees, Cactus alignments, or consensus construction.

Notes
-----
- A 3-way polytomy enumerates all three binary guide trees internally.
- A polytomy with >3 direct children uses
  generate_random_guidetrees_2models_2modes_finalver.py without an outgroup.
- Task dependencies are planned statically. Local reference selection is
  deferred until the shell command that actually needs a reference is executed.

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


def _ascii_tree_lines(node: Node, prefix: str = "", is_last: bool = True) -> List[str]:
    label = node.name or "Unnamed"

    if prefix == "":
        lines = [label]
    else:
        connector = "└── " if is_last else "├── "
        lines = [prefix + connector + label]

    if node.children:
        child_prefix = prefix + ("    " if is_last else "│   ")
        for i, child in enumerate(node.children):
            child_is_last = (i == len(node.children) - 1)
            lines.extend(_ascii_tree_lines(child, child_prefix, child_is_last))
    return lines


def ascii_tree_view(root: Node) -> str:
    return "\n".join(_ascii_tree_lines(root, prefix="", is_last=True))


def leaf_names(root: Node) -> List[str]:
    names: List[str] = []

    def walk(n: Node) -> None:
        if n.is_leaf():
            if not n.name:
                raise ValueError("Encountered an unnamed leaf.")
            names.append(n.name)
            return
        for c in n.children:
            walk(c)

    walk(root)
    return names


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

def is_consensus_target(node: Node) -> bool:
    """Return True when the current node is a polytomy."""
    return (not node.is_leaf()) and len(node.children) > 2


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

    def command_for(self, taxa: List[str]) -> List[str]:
        # Local references do not participate in guide-tree generation.
        # Therefore the guide-tree generator is never given an outgroup.
        cmd = [self.python_bin, str(self.generator), "--taxa", *taxa]
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
    reference_selector: Path
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
        final_reference: Optional[str],
        taxon_paths: Dict[str, str],
        paths_file: Path,
        outdir: Path,
        guide_params: GuideParams,
        pipeline_params: PipelineParams,
    ):
        self.root = root
        self.final_reference = final_reference
        self.taxon_paths = taxon_paths
        self.paths_file = paths_file.resolve()
        self.outdir = outdir.resolve()
        self.guide_params = guide_params
        self.pipeline_params = pipeline_params
        self.final_reference_was_user_provided = final_reference is not None
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

    def reference_selector_expression(self, genome_txt: Path) -> str:
        """
        Return a shell command-substitution expression that resolves the
        longest-genome label at run time.

        Example:
            "$(python /path/select_longest_reference.py task.genome.txt)"
        """
        selector = self.render_run_path(self.pipeline_params.reference_selector)
        genome_file = self.render_run_path(genome_txt)
        return f'"$(python {selector} {genome_file})"'

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

    def task_direct_children(self, task: Task) -> List[Node]:
        return task.node.children

    def ensure_task_for_materialized_node(
        self,
        node: Node,
        parent_task: Task,
        is_root: bool = False,
    ) -> Task:
        kind = "consensus" if is_consensus_target(node) else "fixed"
        task = self.create_task(node, parent_task, is_root, kind)

        if kind == "consensus":
            for child in node.children:
                if not child.is_leaf():
                    self.ensure_task_for_materialized_node(child, task, False)
        else:
            self.collect_consensus_descendants_for_fixed(node, task)

        return task

    def collect_consensus_descendants_for_fixed(
        self,
        node: Node,
        current_task: Task,
    ) -> None:
        def walk(n: Node) -> None:
            if n is not node and (not n.is_leaf()) and is_consensus_target(n):
                current_task.descendant_consensus_names.add(n.name)
                self.ensure_task_for_materialized_node(n, current_task, False)
                return

            for c in n.children:
                walk(c)

        for c in node.children:
            walk(c)

    def plan(self) -> Task:
        root_kind = "consensus" if is_consensus_target(self.root) else "fixed"
        self.root_task = self.create_task(self.root, None, True, root_kind)

        if root_kind == "consensus":
            for child in self.root.children:
                if not child.is_leaf():
                    self.ensure_task_for_materialized_node(
                        child,
                        self.root_task,
                        False,
                    )
        else:
            self.collect_consensus_descendants_for_fixed(
                self.root,
                self.root_task,
            )

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
        cmd = self.guide_params.command_for(labels)
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

    def consensus_tree_strings(self, task: Task) -> List[str]:
        children = self.task_direct_children(task)
        labels = [self.unit_label(u) for u in children]

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

        return [
            s + ";"
            for s in self.generated_trees_no_outgroup(
                labels=labels,
                top_name=top_name,
                internal_prefix=internal_prefix,
            )
        ]

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
            # Only genomes directly participating in this consensus task are
            # listed. The local reference is selected from this exact set at run time.
            for u in self.task_direct_children(task):
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

    def resolve_final_reference_for_single_task(self) -> None:
        """
        Resolve the final export reference when no explicit reference is provided.

        In --noFixedRef mode, the final HAL-to-MAF export reference is selected
        automatically as the longest genome from the original input path file.
        This reference is only used for final export and does not affect guide
        trees or consensus construction.
        """
        if self.final_reference is not None:
            return

        selector = self.pipeline_params.reference_selector
        result = subprocess.run(
            [
                "python",
                str(selector),
                str(self.paths_file),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        self.final_reference = result.stdout.strip()

    def is_single_level_polytomy_workflow(self) -> bool:
        """
        Return True only for a root-level polytomy whose direct children are all
        leaves. This is the only case where the automatically selected reference
        is also the consensus extraction reference and therefore does not need an
        extra final-export note in instruction.txt.
        """
        return (
            is_consensus_target(self.root)
            and all(child.is_leaf() for child in self.root.children)
        )

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
        """
        Return extant taxa with --final_reference first, followed by the
        remaining taxa in deterministic lexical order.
        """
        taxa = sorted(leaf_names(self.root))
        others = [taxon for taxon in taxa if taxon != self.final_reference]
        return [self.final_reference] + others

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
            for u in self.task_direct_children(task):
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

        jobstore = (
            self.root_task.folder
            / f"jobstorehal2maf_{self.root_task.name}_allTaxa"
        )

        return (
            f"cactus-hal2maf --outType single --chunkSize 500000 "
            f"--refGenome {shlex.quote(self.final_reference)} --noAncestors "
            f"{self.render_run_path(jobstore)} "
            f"{self.render_run_path(self.root_alltaxa_hal())} "
            f"{self.render_run_path(self.root_alltaxa_maf())}"
        )

    def root_alltaxa_fasta_cmd(self) -> str:
        taxa_args = " ".join(shlex.quote(t) for t in self.all_leaf_taxa())
        script = self.render_run_path(self.pipeline_params.maf_to_concat_fasta)

        return (
            f"taxa_args=({taxa_args})\n"
            f"python {script} \"${{taxa_args[@]}}\" < "
            f"{self.render_run_path(self.root_alltaxa_maf())} > "
            f"{self.render_run_path(self.root_alltaxa_fasta())}"
        )

    def write_finalization_script(self) -> Optional[Path]:
        """
        Write the finalization script.

        For workflows without regrafting, RunPipelineUseThis.sh already
        produces the final consensus HAL/MAF/FASTA outputs. This script only
        copies them into final_output with standardized names.

        For workflows with regrafting, the final HAL is assembled first and
        then exported to MAF and FASTA.
        """
        regraft_tasks = self.final_regraft_tasks()
        script_path = self.root_alltaxa_script()

        if self.root_task is None:
            if script_path.exists():
                script_path.unlink()
            return None

        with open(script_path, "w", encoding="utf-8") as fw:
            fw.write("#!/bin/bash\n")
            fw.write("set -euo pipefail\n\n")
            fw.write("RUN_PATH=" + shlex.quote(str(self.outdir)) + "\n\n")
            fw.write(f'mkdir -p {self.render_run_path(self.final_output_dir())}\n\n')

            if regraft_tasks:
                fw.write("# Regrafting workflow\n\n")
                bundled_hal_append_subtree = self.bundled_hal_append_subtree().resolve()
                fw.write("patched_halAppendSubtree=" + shlex.quote(str(bundled_hal_append_subtree)) + "\n")
                fw.write('if [[ ! -f "$patched_halAppendSubtree" ]]; then\n')
                fw.write('    echo "Error: bundled executable not found"\n')
                fw.write("    exit 1\n")
                fw.write("fi\n")
                fw.write('chmod +x "$patched_halAppendSubtree"\n\n')

                fw.write(f'root_alltaxa_hal="{self.render_run_path(self.root_alltaxa_hal())}"\n')
                fw.write(f'cp {self.render_run_path(self.task_primary_hal(self.root_task))} "$root_alltaxa_hal"\n')

                for task in regraft_tasks:
                    fw.write(self.regraft_cmd_for_script(self.root_alltaxa_hal(), task) + "\n")

                fw.write("\n")
                fw.write(self.root_alltaxa_hal2maf_cmd() + "\n")
                fw.write(self.root_alltaxa_fasta_cmd() + "\n\n")

                fw.write(f'mv {self.render_run_path(self.root_alltaxa_hal())} {self.render_run_path(self.final_output_hal())}\n')
                fw.write(f'mv {self.render_run_path(self.root_alltaxa_maf())} {self.render_run_path(self.final_output_maf())}\n')
                fw.write(f'mv {self.render_run_path(self.root_alltaxa_fasta())} {self.render_run_path(self.final_output_fasta())}\n')

            else:
                fw.write("# No regrafting is required. Final output handling depends on the root task type.\n\n")

                if self.root_task.kind == "consensus":
                    fw.write("# Root consensus task: RunPipelineUseThis.sh already produced the final consensus HAL/MAF/FASTA.\n")
                    fw.write("# Only copy/rename outputs into final_output.\n\n")

                    fw.write(f'cp {self.render_run_path(self.task_primary_hal(self.root_task))} {self.render_run_path(self.final_output_hal())}\n')
                    fw.write(f'cp {self.render_run_path(self.task_primary_maf(self.root_task))} {self.render_run_path(self.final_output_maf())}\n')
                    fw.write(f'cp {self.render_run_path(self.task_primary_fasta(self.root_task))} {self.render_run_path(self.final_output_fasta())}\n')

                else:
                    fw.write("# Root fixed task: the root alignment is produced by cactus as <Root>_aln1.hal.\n")
                    fw.write("# Export HAL -> MAF using the final export reference, then convert MAF -> FASTA.\n\n")

                    fw.write(f'cp {self.render_run_path(self.task_primary_hal(self.root_task))} {self.render_run_path(self.root_alltaxa_hal())}\n')
                    fw.write(self.root_alltaxa_hal2maf_cmd() + "\n")
                    fw.write(self.root_alltaxa_fasta_cmd() + "\n\n")

                    fw.write(f'mv {self.render_run_path(self.root_alltaxa_hal())} {self.render_run_path(self.final_output_hal())}\n')
                    fw.write(f'mv {self.render_run_path(self.root_alltaxa_maf())} {self.render_run_path(self.final_output_maf())}\n')
                    fw.write(f'mv {self.render_run_path(self.root_alltaxa_fasta())} {self.render_run_path(self.final_output_fasta())}\n')

            fw.write(f'cp {self.render_run_path(self.root_task.output_fa)} {self.render_run_path(self.final_output_ancestor_fasta())}\n')

        script_path.chmod(0o755)
        return script_path

    def write_final_output_summary(self, fw) -> None:
        fw.write("\n")
        fw.write("############################################################\n")
        fw.write("# Final output summary\n")
        fw.write(f"# RUN_PATH = {self.outdir}\n")
        fw.write(
            f"# Final export reference (if HAL-to-MAF export was required): {self.final_reference}\n"
        )

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
        return ascii_tree_view(self.root)

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
        if task.genome_txt is None:
            raise RuntimeError(
                f"Consensus task {task.name!r} has no genome-path file."
            )

        reference_expr = self.reference_selector_expression(task.genome_txt)
        cmds: List[str] = []

        for i, aln_txt in enumerate(task.alignment_files, 1):
            hal = task.folder / f"{task.name}_aln{i}.hal"
            maf = task.folder / f"{task.name}_aln{i}.maf"
            jobstore = task.folder / f"jobstorecluster_{task.name}_aln{i}"
            jobstore_hal2maf = (
                task.folder / f"jobstorehal2maf_{task.name}_aln{i}"
            )

            cmd1 = (
                f"cactus {self.render_run_path(jobstore)} "
                f"{self.render_run_path(aln_txt)} "
                f"{self.render_run_path(hal)}"
            )

            cmd2 = (
                f"cactus-hal2maf --outType single --chunkSize 500000 "
                f"--refGenome {reference_expr} --noAncestors "
                f"{self.render_run_path(jobstore_hal2maf)} "
                f"{self.render_run_path(hal)} "
                f"{self.render_run_path(maf)}"
            )

            cmds.append(cmd1 + "\n" + cmd2)

        return cmds

    def consensus_pipeline_cmd(self, task: Task) -> str:
        if task.genome_txt is None:
            raise RuntimeError(
                f"Consensus task {task.name!r} has no genome-path file."
            )

        maf_paths = [
            self.render_run_path(
                task.folder / f"{task.name}_aln{i}.maf"
            )
            for i in range(1, len(task.alignment_files) + 1)
        ]

        pre = self.render_run_path(task.folder / "consensus")
        genome_txt = self.render_run_path(task.genome_txt)
        pipeline = self.render_run_path(self.pipeline_params.pipeline)
        model_file = self.render_run_path(self.pipeline_params.model_file)
        reference_expr = self.reference_selector_expression(task.genome_txt)

        option_lines = [
            ("--input", " ".join(maf_paths)),
            ("--reference", reference_expr),
            ("--pre", pre),
            ("--threads", str(self.pipeline_params.threads)),
            ("--common_workers", str(self.pipeline_params.common_workers)),
            ("--genome_path", genome_txt),
            ("--Anc_name", shlex.quote(task.name)),
            ("--ModelFile", model_file),
            ("--ifDeleteImmdiFiles", "0"),
            ("--ifRefNonOutgroup", "0"),
        ]

        if self.pipeline_params.sep_length is not None:
            option_lines.append(
                ("--sep_length", str(self.pipeline_params.sep_length))
            )
        if self.pipeline_params.chrom_length_threshold is not None:
            option_lines.append(
                (
                    "--chrom_length_threshold",
                    str(self.pipeline_params.chrom_length_threshold),
                )
            )
        if self.pipeline_params.global_num is not None:
            option_lines.append(
                ("--global_num", str(self.pipeline_params.global_num))
            )
        if self.pipeline_params.separate_workers is not None:
            option_lines.append(
                (
                    "--separate_workers",
                    str(self.pipeline_params.separate_workers),
                )
            )

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
            fw.write("## Input tree used for planning (rooting preserved; auto-named internal nodes if needed):\n")
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

            if (
                not self.final_reference_was_user_provided
                and not self.is_single_level_polytomy_workflow()
            ):
                fw.write(
                    "## The final HAL-to-MAF export reference was automatically selected as the longest genome from the input genome path file.\n"
                )
                fw.write(
                    "## If a different export reference is preferred, manually modify the --refGenome option in finalization.sh.\n\n"
                )

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

    default_generator = (
        script_dir
        / "generate_random_guidetrees_2models_2modes_finalver.py"
    )
    default_pipeline = script_dir / "RunPipelineUseThis.sh"
    default_reference_selector = script_dir / "select_longest_reference.py"
    default_model_file = script_dir / "tryMLstartree.mod"
    default_maf_to_concat_fasta = script_dir / "maf_to_concat_fasta.py"

    ap = argparse.ArgumentParser(
        description=(
            "Plan hierarchical Cactus runs from a rooted, partially "
            "resolved tree using per-consensus-task local references."
        )
    )

    ap.add_argument(
        "--tree",
        required=True,
        help=(
            "Input rooted Newick tree string. The input rooting is preserved; "
            "the planner does not reroot to any reference."
        ),
    )
    ap.add_argument(
        "--final_reference",
        required=False,
        default=None,
        help=(
            "Reference genome used only when exporting the final complete "
            "consensus HAL to MAF. Required for hierarchical workflows with "
            "multiple tasks. For a single-task workflow, the automatically "
            "selected local reference is reused as the final export reference."
        ),
    )
    ap.add_argument(
        "--paths",
        required=True,
        type=Path,
        help="Two-column file: taxon path",
    )
    ap.add_argument(
        "--outdir",
        default=".",
        type=Path,
        help="Output directory (default: current working directory)",
    )

    ap.add_argument(
        "--generator",
        default=str(default_generator),
        help="Path to generate_random_guidetrees_2models_2modes_finalver.py",
    )
    ap.add_argument(
        "--guide-num-trees",
        type=int,
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-model",
        choices=["uniform", "yule"],
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-rf-threshold",
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-t-threshold",
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-max-tries",
        type=int,
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-restarts",
        type=int,
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--guide-seed",
        type=int,
        default=None,
        help="Pass through to the guide-tree generator",
    )
    ap.add_argument(
        "--generator-python",
        default=sys.executable,
        help="Python executable used to run the guide-tree generator",
    )

    ap.add_argument(
        "--pipeline",
        default=str(default_pipeline),
        help="Path to RunPipelineUseThis.sh",
    )
    ap.add_argument(
        "--reference-selector",
        default=str(default_reference_selector),
        help=(
            "Path to select_longest_reference.py. Default: script with this "
            "name in the same directory as the planner."
        ),
    )
    ap.add_argument(
        "--maf_to_concat_fasta",
        default=str(default_maf_to_concat_fasta),
        help="Path to maf_to_concat_fasta.py",
    )
    ap.add_argument(
        "--threads",
        required=True,
        type=int,
        help="Pass through to RunPipelineUseThis.sh --threads",
    )
    ap.add_argument(
        "--common_workers",
        required=True,
        type=int,
        help="Pass through to RunPipelineUseThis.sh --common_workers",
    )
    ap.add_argument(
        "--ModelFile",
        default=str(default_model_file),
        type=Path,
        help=(
            "Base model file used by RunPipelineUseThis.sh; "
            "default: tryMLstartree.mod in the planner directory"
        ),
    )
    ap.add_argument(
        "--sep_length",
        type=int,
        default=None,
        help="Pass through to RunPipelineUseThis.sh",
    )
    ap.add_argument(
        "--chrom_length_threshold",
        type=int,
        default=None,
        help="Pass through to RunPipelineUseThis.sh",
    )
    ap.add_argument(
        "--global_num",
        type=int,
        default=None,
        help="Pass through to RunPipelineUseThis.sh",
    )
    ap.add_argument(
        "--separate_workers",
        type=int,
        default=None,
        help="Pass through to RunPipelineUseThis.sh",
    )

    return ap


def main() -> None:
    ap = build_arg_parser()
    args = ap.parse_args()

    taxon_paths = read_taxon_paths(args.paths)

    raw_root = parse_newick(args.tree)
    validate_existing_node_names(raw_root)

    # No global reference exists in the planning algorithm. Preserve the
    # supplied rooting and only collapse redundant unary nodes.
    root = collapse_unary_nodes(raw_root)
    auto_name_internal_nodes(root)

    missing_taxa = sorted(
        taxon
        for taxon in set(leaf_names(root))
        if taxon not in taxon_paths
    )
    if missing_taxa:
        raise ValueError(
            "The following leaf taxa are missing from --paths: "
            + ", ".join(repr(x) for x in missing_taxa)
        )

    leaf_taxa = set(leaf_names(root))
    if args.final_reference is not None:
        if args.final_reference not in leaf_taxa:
            raise ValueError(
                f"Final reference {args.final_reference!r} is not a leaf taxon "
                "in the input tree."
            )

        if args.final_reference not in taxon_paths:
            raise ValueError(
                f"Final reference {args.final_reference!r} is not present in --paths."
            )

    generator = Path(args.generator).resolve()
    if not generator.is_file():
        raise FileNotFoundError(
            f"Guide-tree generator script not found: {generator}"
        )

    pipeline = Path(args.pipeline).resolve()
    if not pipeline.is_file():
        raise FileNotFoundError(
            f"Consensus pipeline script not found: {pipeline}"
        )

    reference_selector = Path(args.reference_selector).resolve()
    if not reference_selector.is_file():
        raise FileNotFoundError(
            f"Reference selector script not found: {reference_selector}"
        )

    model_file = Path(args.ModelFile).resolve()
    if not model_file.is_file():
        raise FileNotFoundError(f"Model file not found: {model_file}")

    maf_to_concat_fasta = Path(args.maf_to_concat_fasta).resolve()
    if not maf_to_concat_fasta.is_file():
        raise FileNotFoundError(
            f"maf_to_concat_fasta.py not found: {maf_to_concat_fasta}"
        )

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
        reference_selector=reference_selector,
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
        final_reference=args.final_reference,
        taxon_paths=taxon_paths,
        paths_file=args.paths,
        outdir=args.outdir,
        guide_params=guide_params,
        pipeline_params=pipeline_params,
    )

    planner.plan()
    planner.resolve_final_reference_for_single_task()
    planner.write_all_task_files()
    planner.write_finalization_script()
    instruction = planner.write_instruction()

    print("Done.")
    print(f"Instruction file: {instruction}")
    if planner.root_task is not None:
        print(f"Root task folder: {planner.root_task.folder}")


if __name__ == "__main__":
    main()
