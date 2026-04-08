# -*- coding: utf-8 -*-
import argparse
import random
import re
import tempfile
from math import comb, ceil, floor, prod
from typing import List, Set, Dict, Tuple
from ete3 import Tree
import textwrap
# =============================================================================
# 0. Basic utilities
# =============================================================================

def number_of_rooted_binary_trees(n: int) -> int:
    if n < 2:
        return 1
    m = 2 * n - 3
    return prod(range(m, 0, -2))

def gen_placeholders(n: int) -> List[str]:
    return [f"T{i+1}" for i in range(n)]

def replace_placeholders_raw(nwk: str, mapping: Dict[str, str]) -> str:
    """Token-wise boundary replacement with a fallback scan to avoid chain replacements or missed replacements."""
    for tok in sorted(mapping.keys(), key=len, reverse=True):
        pat = re.compile(r'(?<=\(|,)\s*' + re.escape(tok) + r'\s*(?=[:),;])')
        nwk = pat.sub(mapping[tok], nwk)
    leftovers = re.findall(r'(?<=\(|,)\s*(T\d+|__OUTGROUP__)\s*(?=[:),;])', nwk)
    if leftovers:
        token_pat = re.compile(r'(?<=\(|,)\s*([^:),;]+?)\s*(?=[:),;])')
        def _repl(m: re.Match) -> str:
            tok = m.group(1).strip()
            return mapping.get(tok, tok)
        nwk = token_pat.sub(_repl, nwk)
    return nwk

try:
    _ = (1).bit_count
    def popcount(x: int) -> int:
        return x.bit_count()
except AttributeError:
    def popcount(x: int) -> int:
        return bin(x).count("1")

# =============================================================================
# 1. Random tree generation (Uniform / Yule)
# =============================================================================

def generate_uniform_tree(ingroup: List[str]) -> Tree:
    if not ingroup:
        return Tree()
    labs = random.sample(ingroup, len(ingroup))
    if len(labs) == 1:
        return Tree(name=labs[0])
    root = Tree()
    root.add_child(name=labs[0])
    root.add_child(name=labs[1])
    nodes = [root] + root.get_descendants()
    for lab in labs[2:]:
        chosen = random.choice(nodes)
        if chosen.is_root():
            new_root = Tree()
            new_root.add_child(chosen)
            new_root.add_child(name=lab)
            root = new_root
            nodes.append(new_root)
        else:
            parent = chosen.up
            new_internal = Tree()
            chosen.detach()
            new_internal.add_child(chosen)
            new_internal.add_child(name=lab)
            parent.add_child(new_internal)
            nodes.append(new_internal)
    return root

def generate_yule_tree(ingroup: List[str]) -> Tree:
    if not ingroup:
        return Tree()
    nodes = [Tree(name=n) for n in ingroup]
    random.shuffle(nodes)
    while len(nodes) > 1:
        a = nodes.pop(random.randrange(len(nodes)))
        b = nodes.pop(random.randrange(len(nodes)))
        p = Tree()
        p.add_child(a); p.add_child(b)
        nodes.append(p)
    return nodes[0] if nodes else Tree()

def generate_random_tree(ingroup: List[str], model: str = "uniform") -> Tree:
    m = (model or "uniform").lower()
    if m == "uniform": return generate_uniform_tree(ingroup)
    if m == "yule":    return generate_yule_tree(ingroup)
    raise ValueError(f"Unknown model type: {model}")

def attach_outgroup(rooted_ingroup: Tree, outgroup: str) -> Tree:
    if outgroup in rooted_ingroup.get_leaf_names():
        raise ValueError(f"Outgroup '{outgroup}' cannot also be present in the ingroup.")
    final_tree = Tree()
    final_tree.add_child(name=outgroup)
    final_tree.add_child(rooted_ingroup.copy(method="deepcopy"))
    return final_tree

# =============================================================================
# 2. RF utilities
# =============================================================================

def canonical_newick_rooted(t: Tree) -> str:
    def enc(n: Tree) -> str:
        if n.is_leaf(): return n.name
        parts = sorted([enc(ch) for ch in n.get_children()])
        return "(" + ",".join(parts) + ")"
    return enc(t)

def _leaf_index(leaf_names: List[str]) -> Dict[str, int]:
    return {name: i for i, name in enumerate(leaf_names)}

def clade_masks_bitset(t: Tree, leaf_index: Dict[str, int], n: int) -> Set[int]:
    masks: Set[int] = set()
    stack: List[Tuple[Tree, bool]] = [(t, False)]
    acc: Dict[int, int] = {}
    while stack:
        node, vis = stack.pop()
        if not vis:
            stack.append((node, True))
            for ch in node.get_children():
                stack.append((ch, False))
        else:
            if node.is_leaf():
                mask = 1 << leaf_index[node.name]
            else:
                mask = 0
                for ch in node.get_children():
                    mask |= acc[id(ch)]
                k = popcount(mask)
                if 1 < k < n:
                    masks.add(mask)
            acc[id(node)] = mask
    return masks

def allowed_shared_clades(n: int, tau: float) -> int:
    return floor((1.0 - tau) * (n - 2))

def shared_exceeds(candidate_masks: Set[int], kept_masks: Set[int], s_max: int) -> bool:
    shared = 0
    if len(candidate_masks) <= len(kept_masks):
        itA, itB = candidate_masks, kept_masks
    else:
        itA, itB = kept_masks, candidate_masks
    for m in itA:
        if m in itB:
            shared += 1
            if shared > s_max:
                return True
    return False

# =============================================================================
# 3. Built-in Triplet distance
#    Idea: for each unordered leaf pair {i,j}, take the LCA-subtree masks M1 and M2 in the two trees.
#          The number of k values for which both trees classify the pair as (i,j|k) is n - |M1 union M2|.
#          Consistent triplets are summed over all leaf pairs; raw = C(n,3) - consistent.
# =============================================================================

def preprocess_pair_lca_masks(t: Tree, leaf_index: Dict[str, int]) -> Tuple[Dict[Tuple[int,int], int], int]:
    """Return pair_mask[(i,j)] = the LCA(i,j) subtree mask (i<j), and the number of leaves n."""
    root = t
    parent: Dict[int, int] = {}
    depth: Dict[int, int] = {id(root): 0}
    node_mask: Dict[int, int] = {}
    leaf_node_id: Dict[int, int] = {}

    # Preorder: parent / depth
    stack = [root]
    while stack:
        node = stack.pop()
        for ch in node.get_children():
            parent[id(ch)] = id(node)
            depth[id(ch)] = depth[id(node)] + 1
            stack.append(ch)

    # Postorder: subtree leaf masks
    stack = [(root, False)]
    while stack:
        node, vis = stack.pop()
        if not vis:
            stack.append((node, True))
            for ch in node.get_children():
                stack.append((ch, False))
        else:
            if node.is_leaf():
                m = 1 << leaf_index[node.name]
                leaf_node_id[leaf_index[node.name]] = id(node)
            else:
                m = 0
                for ch in node.get_children():
                    m |= node_mask[id(ch)]
            node_mask[id(node)] = m

    n = len(leaf_index)

    def lca_mask(id_a: int, id_b: int) -> int:
        a, b = id_a, id_b
        da, db = depth[a], depth[b]
        while da > db:
            a = parent[a]; da -= 1
        while db > da:
            b = parent[b]; db -= 1
        while a != b:
            a = parent[a]; b = parent[b]
        return node_mask[a]

    pair_mask: Dict[Tuple[int,int], int] = {}
    idxs = list(range(n))
    for i in idxs:
        for j in range(i+1, n):
            m = lca_mask(leaf_node_id[i], leaf_node_id[j])
            pair_mask[(i,j)] = m
    return pair_mask, n

def triplet_distance_raw_from_pairmasks(pairA: Dict[Tuple[int,int], int],
                                        pairB: Dict[Tuple[int,int], int],
                                        n: int) -> int:
    agree = 0
    for (i,j), m1 in pairA.items():
        m2 = pairB[(i,j)]
        agree += (n - popcount(m1 | m2))
    total = comb(n, 3)
    return total - agree

def triplet_distance_raw_between_trees(t1: Tree, t2: Tree) -> Tuple[int, float]:
    """Given two trees (with identical leaf sets), return (raw, normalized)."""
    leaves1 = set(t1.get_leaf_names())
    leaves2 = set(t2.get_leaf_names())
    if leaves1 != leaves2:
        raise ValueError("The two trees have different leaf sets, so the triplet distance cannot be computed.")
    # Build the index in a stable order (lexicographic sorting, to ensure consistency across the two trees)
    ordered = sorted(leaves1)
    li = _leaf_index(ordered)
    pm1, n1 = preprocess_pair_lca_masks(t1, li)
    pm2, n2 = preprocess_pair_lca_masks(t2, li)
    assert n1 == n2
    raw = triplet_distance_raw_from_pairmasks(pm1, pm2, n1)
    denom = comb(n1, 3) if n1 >= 3 else 0
    norm = (raw / denom) if denom > 0 else 0.0
    return raw, norm

# =============================================================================
# 4. Single-attempt sampling logic (Worker)
# =============================================================================

def sample_far_trees_single_attempt(
    ingroup: List[str], num_trees: int,
    rf_mode: str, rf_tau: float,
    t_norm_threshold: float,
    max_tries: int, model_type: str, seed: int = None,
    progress_every: int = 10000
) -> List[Tree]:

    if seed is not None:
        random.seed(seed)

    n = len(ingroup)
    triplet_raw_min = ceil(t_norm_threshold * comb(n, 3))
    max_unique_topologies = number_of_rooted_binary_trees(n)

    kept_trees: List[Tree] = []
    visited_topologies: Set[str] = set()

    leaf_index = _leaf_index(ingroup)
    if rf_mode == "strict":
        union_clades: Set[int] = set()
    else:
        kept_clade_sets: List[Set[int]] = []
        s_max = allowed_shared_clades(n, rf_tau)

    # Cache pair->mask for the trees already selected
    kept_pair_masks: List[Dict[Tuple[int,int], int]] = []

    tries = 0
    early_stop = False
    while len(kept_trees) < num_trees and tries < max_tries:
        tries += 1
        if progress_every and tries % progress_every == 0:
            print(f"    - [Progress] {tries} attempts tried, still searching...", flush=True)

        if len(visited_topologies) >= max_unique_topologies:
            early_stop = True
            break

        cand = generate_random_tree(ingroup, model=model_type)
        topo_key = canonical_newick_rooted(cand)
        if topo_key in visited_topologies:
            continue
        visited_topologies.add(topo_key)

        # ---- RF filtering ----
        cand_masks = clade_masks_bitset(cand, leaf_index, n)
        rf_ok = True
        if rf_mode == "strict":
            if union_clades and not cand_masks.isdisjoint(union_clades):
                rf_ok = False
        else:
            for kept_masks in kept_clade_sets:
                if shared_exceeds(cand_masks, kept_masks, s_max):
                    rf_ok = False
                    break
        if not rf_ok:
            continue

        # ---- Triplet filtering (built-in algorithm) ----
        cand_pairmask, _ = preprocess_pair_lca_masks(cand, leaf_index)

        ok = True
        for pm_kept in kept_pair_masks:
            raw = triplet_distance_raw_from_pairmasks(cand_pairmask, pm_kept, n)
            if raw < triplet_raw_min:
                ok = False
                break
        if not ok:
            continue

        # Accepted: add it to the selected set and cache its pair->mask
        kept_trees.append(cand)
        kept_pair_masks.append(cand_pairmask)
        if rf_mode == "strict":
            union_clades |= cand_masks
        else:
            kept_clade_sets.append(cand_masks)

        print(f"    - (After {tries} attempts) found tree #{len(kept_trees)}...", flush=True)

    if len(kept_trees) < num_trees:
        if early_stop:
            print(f"    - [Info] This attempt has already exhausted all {max_unique_topologies} topologies; no more trees can be found.")
        else:
            print(f"    - [Info] This attempt has reached the maximum number of tries: {max_tries}.")

    return kept_trees

# =============================================================================
# 5. Restart-enabled logic (Manager)
# =============================================================================

def run_with_restarts(
    ingroup_internal: List[str], outgroup_internal: str, num_trees: int, restarts: int,
    rf_mode: str, rf_tau: float, t_norm_threshold: float, max_tries: int,
    model_type: str, seed: int = None
) -> List[Tree]:

    best: List[Tree] = []

    base_seed = seed if seed is not None else random.randrange(1 << 30)
    print(f"[Info] Using base seed: {base_seed} (it will increase by 1 for each restart)")

    for attempt in range(restarts + 1):
        print(f"\n--- Attempt {attempt + 1}/{restarts + 1} (model={model_type.upper()}, mode={rf_mode}, RF threshold={rf_tau}) ---")
        cur_seed = base_seed + attempt

        found = sample_far_trees_single_attempt(
            ingroup=ingroup_internal, num_trees=num_trees,
            rf_mode=rf_mode, rf_tau=rf_tau,
            t_norm_threshold=t_norm_threshold,
            max_tries=max_tries,
            model_type=model_type, seed=cur_seed
        )

        if len(found) > len(best):
            best = found
        if len(found) == num_trees:
            print(f"\n[Success] Found {num_trees} trees in attempt {attempt + 1}.")
            return best

    print(f"\n[Info] All {restarts + 1} attempts have finished, but {num_trees} trees could not be fully collected.")
    print(f"Returning the best result set found ({len(best)} trees).")
    return best

# =============================================================================
# 6. Command-line interface
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Efficiently sample a set of rooted binary trees with sufficiently large topological distances (random generation + RF filtering + Triplet filtering + deduplication + early stopping + multiple restarts).\n"
            "Supports Uniform / Yule random models; configurable RF threshold; configurable Triplet distance threshold.\n"
            "The default evolutionary model is yule, with strict filtering, RF=1 and triplet_distance=2/3.\n"
            "Under the default settings, when the number of ingroup taxa n is less than or equal to 10, at most n valid trees can usually be produced; when n > 10, --num_trees 10 is recommended."
        ),
        epilog = textwrap.dedent("""\
        Examples:
        --taxa A B C D E G H I --outgroup F --num_trees 8  --rf-threshold 1   --model uniform   # strict
        --taxa A B C D E G H I --outgroup F --num_trees 7 --rf-threshold 0.8 --model yule      # relax(0.8)
        --taxa "(A,B)" C D "(E,(G,(H,I)))" --outgroup F
        """),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--taxa", nargs="+", required=True, help="List of ingroup taxa (original labels; spaces/parentheses allowed). If a taxon is a fixed group, wrap it in quotes.")
    parser.add_argument("--outgroup", required=True, help="Outgroup label used for rooting (original label)")
    parser.add_argument("--num_trees", type=int, default=4, help="Number of trees to generate (default: 4)")

    parser.add_argument("--model", choices=['uniform', 'yule'], default='yule',
                        help='Model used to generate random trees ("uniform" or "yule", default: yule)')
    # <<< Replacement: remove --mode and add --rf-threshold >>>
    parser.add_argument("--rf-threshold", type=float, default=1.0,
                        help="Normalized RF distance threshold tau in (0,1]. =1 means strict (RFnorm=1.0), "
                             "0<tau<1 means relax (RFnorm>=tau). Default: 1.0.")

    parser.add_argument("--t-threshold", type=float, default=2/3,
                        help="Normalized Triplet distance threshold (default: 2/3)")
    parser.add_argument("--max-tries", type=int, default=200000,
                        help="Maximum number of tries per attempt (default: 200000)")
    parser.add_argument("--restarts", type=int, default=10,
                        help="Maximum number of random restarts when enough trees cannot be found (default: 10)")
    parser.add_argument("--seed", type=int, default=None,
                        help="Initial random seed for reproducibility (default: None)")

    args = parser.parse_args()

    # <<< New logic: determine rf_mode and rf_tau from --rf-threshold >>>
    if args.rf_threshold is None or args.rf_threshold >= 1.0:
        rf_mode = "strict"
        rf_tau  = 1.0
    else:
        if not (0.0 < args.rf_threshold < 1.0):
            raise ValueError("--rf-threshold must be in (0,1], for example 1 or 0.8")
        rf_mode = "relax"
        rf_tau  = float(args.rf_threshold)

    placeholders = gen_placeholders(len(args.taxa))
    OUT_INT = "__OUTGROUP__"
    int2ext = {placeholders[i]: ext for i, ext in enumerate(args.taxa)}
    mapping = {**int2ext, OUT_INT: args.outgroup}

    try:
        internal_ingroup_trees = run_with_restarts(
            ingroup_internal=placeholders,
            outgroup_internal=OUT_INT,
            num_trees=args.num_trees,
            restarts=args.restarts,
            rf_mode=rf_mode,
            rf_tau=rf_tau,
            t_norm_threshold=args.t_threshold,
            max_tries=args.max_tries,
            model_type=args.model,
            seed=args.seed
        )
    except (RuntimeError, ValueError) as e:
        print(f"\n[Error] Program execution failed: {e}")
        return

    print("\n--- Final generated guide trees ---")
    if not internal_ingroup_trees:
        print("Failed to generate any tree that satisfies the criteria.")
        return

    # Save the internal trees with outgroup attached so tree_1 vs tree_2 distance can be computed later
    full_internal_trees: List[Tree] = []

    print(f"Finally found {len(internal_ingroup_trees)} trees (target: {args.num_trees}):")
    for i, ingroup_tree in enumerate(internal_ingroup_trees, 1):
        try:
            full_tree_internal = attach_outgroup(ingroup_tree, OUT_INT)  # use the internal placeholder name
            full_internal_trees.append(full_tree_internal.copy(method="deepcopy"))

            nwk_internal = full_tree_internal.write(format=9)
            nwk_final = replace_placeholders_raw(nwk_internal, mapping)  # replace back to original names before output
            print(f"\n>tree_{i}\n{nwk_final}")
        except Exception as e:
            print(f"\n[Error] An error occurred while outputting tree_{i}: {e}")

    # ===== Extra: compute and print the triplet distance between tree_1 and tree_2 (using the built-in algorithm) =====
    if len(full_internal_trees) >= 2:
        try:
            t1 = internal_ingroup_trees[0]
            t2 = internal_ingroup_trees[1]
            raw, norm = triplet_distance_raw_between_trees(t1, t2)
            print(f"\n[Extra] Triplet distance for tree_1 vs tree_2: raw={raw}, normalized={norm:.12f}")
        except Exception as e:
            print(f"\n[Extra] Failed to compute the triplet distance for tree_1 vs tree_2: {e}")

if __name__ == "__main__":
    main()
