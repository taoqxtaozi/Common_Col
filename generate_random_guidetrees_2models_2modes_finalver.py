# -*- coding: utf-8 -*-
import argparse
import random
import re
from math import comb, ceil, floor, prod
from typing import List, Set, Dict, Tuple, Optional
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
    """
    Token-wise boundary replacement with a fallback scan to avoid chain replacements
    or missed replacements.

    This allows fixed groups such as "(A,B)" to be passed as one sampled element.
    For example, if T1 maps to "(M,TAEGU)", then:
        ((T1,T2),T3);
    becomes:
        (((M,TAEGU),D),E);
    """
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
        p.add_child(a)
        p.add_child(b)
        nodes.append(p)

    return nodes[0] if nodes else Tree()


def generate_random_tree(ingroup: List[str], model: str = "uniform") -> Tree:
    m = (model or "uniform").lower()

    if m == "uniform":
        return generate_uniform_tree(ingroup)

    if m == "yule":
        return generate_yule_tree(ingroup)

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
        if n.is_leaf():
            return n.name
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
#
# Idea:
# For each unordered leaf pair {i,j}, take the LCA-subtree masks M1 and M2
# in the two trees. The number of k values for which both trees classify the
# pair as (i,j|k) is n - |M1 union M2|.
# Consistent triplets are summed over all leaf pairs.
# raw = C(n,3) - consistent.
# =============================================================================

def preprocess_pair_lca_masks(t: Tree, leaf_index: Dict[str, int]) -> Tuple[Dict[Tuple[int, int], int], int]:
    """
    Return pair_mask[(i,j)] = the LCA(i,j) subtree mask (i<j),
    and the number of leaves n.
    """
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
            a = parent[a]
            da -= 1

        while db > da:
            b = parent[b]
            db -= 1

        while a != b:
            a = parent[a]
            b = parent[b]

        return node_mask[a]

    pair_mask: Dict[Tuple[int, int], int] = {}
    idxs = list(range(n))

    for i in idxs:
        for j in range(i + 1, n):
            m = lca_mask(leaf_node_id[i], leaf_node_id[j])
            pair_mask[(i, j)] = m

    return pair_mask, n


def triplet_distance_raw_from_pairmasks(
    pairA: Dict[Tuple[int, int], int],
    pairB: Dict[Tuple[int, int], int],
    n: int
) -> int:
    agree = 0

    for (i, j), m1 in pairA.items():
        m2 = pairB[(i, j)]
        agree += (n - popcount(m1 | m2))

    total = comb(n, 3)
    return total - agree


def triplet_distance_raw_between_trees(t1: Tree, t2: Tree) -> Tuple[int, float]:
    """
    Given two rooted trees with identical leaf sets, return (raw, normalized).
    """
    leaves1 = set(t1.get_leaf_names())
    leaves2 = set(t2.get_leaf_names())

    if leaves1 != leaves2:
        raise ValueError("The two trees have different leaf sets, so the triplet distance cannot be computed.")

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
# 4. Single-attempt sampling logic
# =============================================================================

def sample_far_trees_single_attempt(
    ingroup: List[str],
    num_trees: int,
    rf_mode: str,
    rf_tau: float,
    t_norm_threshold: float,
    max_tries: int,
    model_type: str,
    seed: Optional[int] = None,
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

    kept_pair_masks: List[Dict[Tuple[int, int], int]] = []

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

        # ---- Triplet filtering ----
        cand_pairmask, _ = preprocess_pair_lca_masks(cand, leaf_index)

        ok = True

        for pm_kept in kept_pair_masks:
            raw = triplet_distance_raw_from_pairmasks(cand_pairmask, pm_kept, n)
            if raw < triplet_raw_min:
                ok = False
                break

        if not ok:
            continue

        kept_trees.append(cand)
        kept_pair_masks.append(cand_pairmask)

        if rf_mode == "strict":
            union_clades |= cand_masks
        else:
            kept_clade_sets.append(cand_masks)

        print(f"    - (After {tries} attempts) found tree #{len(kept_trees)}...", flush=True)

    if len(kept_trees) < num_trees:
        if early_stop:
            print(
                f"    - [Info] This attempt has already exhausted all "
                f"{max_unique_topologies} topologies; no more trees can be found."
            )
        else:
            print(f"    - [Info] This attempt has reached the maximum number of tries: {max_tries}.")

    return kept_trees


# =============================================================================
# 5. Restart-enabled logic
# =============================================================================

def run_with_restarts(
    ingroup_internal: List[str],
    num_trees: int,
    restarts: int,
    rf_mode: str,
    rf_tau: float,
    t_norm_threshold: float,
    max_tries: int,
    model_type: str,
    seed: Optional[int] = None
) -> List[Tree]:

    best: List[Tree] = []

    base_seed = seed if seed is not None else random.randrange(1 << 30)
    print(f"[Info] Using base seed: {base_seed} (it will increase by 1 for each restart)")

    for attempt in range(restarts + 1):
        print(
            f"\n--- Attempt {attempt + 1}/{restarts + 1} "
            f"(model={model_type.upper()}, mode={rf_mode}, RF threshold={rf_tau}) ---"
        )

        cur_seed = base_seed + attempt

        found = sample_far_trees_single_attempt(
            ingroup=ingroup_internal,
            num_trees=num_trees,
            rf_mode=rf_mode,
            rf_tau=rf_tau,
            t_norm_threshold=t_norm_threshold,
            max_tries=max_tries,
            model_type=model_type,
            seed=cur_seed
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
            "Efficiently sample a set of rooted binary guide trees with sufficiently large "
            "topological distances. The script supports Uniform / Yule random models, RF "
            "filtering, Triplet-distance filtering, deduplication, early stopping, and "
            "multiple restarts.\n\n"
            "If --outgroup is provided, the script first samples rooted binary trees among "
            "--taxa and then attaches the outgroup outside the sampled tree.\n"
            "If --outgroup is omitted, the script directly samples rooted binary trees among "
            "all elements given by --taxa.\n\n"
            "Fixed groups can be passed as quoted Newick-like elements, e.g. \"(M,TAEGU)\"."
        ),
        epilog=textwrap.dedent("""\
        Examples:

        # 1. Old behavior: reference/outgroup is fixed outside the sampled ingroup tree
        python generate_random_guidetrees.py \\
            --taxa A B C D E G H I \\
            --outgroup F \\
            --num_trees 8 \\
            --rf-threshold 1 \\
            --model uniform

        # 2. No outgroup: all elements are sampled as ordinary elements
        python generate_random_guidetrees.py \\
            --taxa A B C TAEGU \\
            --num_trees 4

        # 3. Nested-reference case handled by the caller/planner:
        #    keep (M,TAEGU) as a fixed group and do not provide --outgroup
        python generate_random_guidetrees.py \\
            --taxa "(M,TAEGU)" D E \\
            --num_trees 3

        # 4. Fixed groups are allowed
        python generate_random_guidetrees.py \\
            --taxa "(A,B)" C D "(E,(G,(H,I)))" \\
            --outgroup F
        """),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--taxa",
        nargs="+",
        required=True,
        help=(
            "List of elements to be sampled. Each element can be a taxon, an inferred ancestor, "
            "or a fixed group such as \"(A,B)\". If a fixed group contains parentheses, wrap it in quotes."
        )
    )

    parser.add_argument(
        "--outgroup",
        default=None,
        help=(
            "Optional outgroup label. If provided, the outgroup is attached outside the sampled "
            "tree. If omitted, no external outgroup is attached."
        )
    )

    parser.add_argument(
        "--num_trees",
        type=int,
        default=4,
        help="Number of trees to generate (default: 4)"
    )

    parser.add_argument(
        "--model",
        choices=["uniform", "yule"],
        default="yule",
        help='Model used to generate random trees ("uniform" or "yule", default: yule)'
    )

    parser.add_argument(
        "--rf-threshold",
        type=float,
        default=1.0,
        help=(
            "Normalized RF distance threshold tau in (0,1]. "
            "=1 means strict (RFnorm=1.0); 0<tau<1 means relax (RFnorm>=tau). "
            "Default: 1.0."
        )
    )

    parser.add_argument(
        "--t-threshold",
        type=float,
        default=2 / 3,
        help="Normalized Triplet distance threshold (default: 2/3)"
    )

    parser.add_argument(
        "--max-tries",
        type=int,
        default=200000,
        help="Maximum number of tries per attempt (default: 200000)"
    )

    parser.add_argument(
        "--restarts",
        type=int,
        default=10,
        help="Maximum number of random restarts when enough trees cannot be found (default: 10)"
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Initial random seed for reproducibility (default: None)"
    )

    args = parser.parse_args()

    if len(args.taxa) == 0:
        raise ValueError("--taxa must contain at least one element.")

    if args.outgroup is not None and args.outgroup in args.taxa:
        raise ValueError(
            f"--outgroup {args.outgroup!r} is also present as an exact element in --taxa. "
            "If the reference should behave as an ordinary element, omit --outgroup."
        )

    if args.rf_threshold is None or args.rf_threshold >= 1.0:
        rf_mode = "strict"
        rf_tau = 1.0
    else:
        if not (0.0 < args.rf_threshold < 1.0):
            raise ValueError("--rf-threshold must be in (0,1], for example 1 or 0.8")
        rf_mode = "relax"
        rf_tau = float(args.rf_threshold)

    placeholders = gen_placeholders(len(args.taxa))
    int2ext = {placeholders[i]: ext for i, ext in enumerate(args.taxa)}

    OUT_INT = "__OUTGROUP__"

    if args.outgroup is not None:
        mapping = {**int2ext, OUT_INT: args.outgroup}
    else:
        mapping = dict(int2ext)

    try:
        internal_sampled_trees = run_with_restarts(
            ingroup_internal=placeholders,
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

    if not internal_sampled_trees:
        print("Failed to generate any tree that satisfies the criteria.")
        return

    if args.outgroup is not None:
        print("[Info] External outgroup mode: the sampled tree is generated among --taxa, then --outgroup is attached.")
    else:
        print("[Info] No-outgroup mode: the sampled tree is generated directly among all --taxa elements.")

    print(f"Finally found {len(internal_sampled_trees)} trees (target: {args.num_trees}):")

    output_internal_trees: List[Tree] = []

    for i, sampled_tree in enumerate(internal_sampled_trees, 1):
        try:
            if args.outgroup is not None:
                output_tree_internal = attach_outgroup(sampled_tree, OUT_INT)
            else:
                output_tree_internal = sampled_tree.copy(method="deepcopy")

            output_internal_trees.append(output_tree_internal.copy(method="deepcopy"))

            nwk_internal = output_tree_internal.write(format=9)
            nwk_final = replace_placeholders_raw(nwk_internal, mapping)

            print(f"\n>tree_{i}\n{nwk_final}")

        except Exception as e:
            print(f"\n[Error] An error occurred while outputting tree_{i}: {e}")
'''
    # Extra: compute and print the triplet distance between tree_1 and tree_2.
    # This is computed over the sampled elements only. If an external outgroup is attached,
    # the outgroup is not included in this distance, matching the filtering step.
    if len(internal_sampled_trees) >= 2:
        try:
            t1 = internal_sampled_trees[0]
            t2 = internal_sampled_trees[1]
            raw, norm = triplet_distance_raw_between_trees(t1, t2)

            print(
                f"\n[Extra] Triplet distance for tree_1 vs tree_2 "
                f"over sampled elements: raw={raw}, normalized={norm:.12f}"
            )

        except Exception as e:
            print(f"\n[Extra] Failed to compute the triplet distance for tree_1 vs tree_2: {e}")
'''	

if __name__ == "__main__":
    main()