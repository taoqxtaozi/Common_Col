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
# 0. 基础工具
# =============================================================================

def number_of_rooted_binary_trees(n: int) -> int:
    if n < 2:
        return 1
    m = 2 * n - 3
    return prod(range(m, 0, -2))

def gen_placeholders(n: int) -> List[str]:
    return [f"T{i+1}" for i in range(n)]

def replace_placeholders_raw(nwk: str, mapping: Dict[str, str]) -> str:
    """逐 token 的边界替换 + 兜底扫描，避免连环/漏替。"""
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
# 1. 随机树生成（Uniform / Yule）
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
    raise ValueError(f"未知的模型类型: {model}")

def attach_outgroup(rooted_ingroup: Tree, outgroup: str) -> Tree:
    if outgroup in rooted_ingroup.get_leaf_names():
        raise ValueError(f"外群 '{outgroup}' 不能同时存在于内群中。")
    final_tree = Tree()
    final_tree.add_child(name=outgroup)
    final_tree.add_child(rooted_ingroup.copy(method="deepcopy"))
    return final_tree

# =============================================================================
# 2. RF 工具
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
# 3. 内置 Triplet 距离
#    思路：对每个无序叶对 {i,j}，取两树 LCA 子树掩码 M1, M2。
#         这对在两树同时判为 (i,j|k) 的 k 的个数 = n - |M1 ∪ M2|。
#         一致 triplet = 对所有叶对求和；raw = C(n,3) - 一致。
# =============================================================================

def preprocess_pair_lca_masks(t: Tree, leaf_index: Dict[str, int]) -> Tuple[Dict[Tuple[int,int], int], int]:
    """返回 pair_mask[(i,j)] = LCA(i,j) 子树掩码（i<j），以及叶子数 n。"""
    root = t
    parent: Dict[int, int] = {}
    depth: Dict[int, int] = {id(root): 0}
    node_mask: Dict[int, int] = {}
    leaf_node_id: Dict[int, int] = {}

    # 前序：parent / depth
    stack = [root]
    while stack:
        node = stack.pop()
        for ch in node.get_children():
            parent[id(ch)] = id(node)
            depth[id(ch)] = depth[id(node)] + 1
            stack.append(ch)

    # 后序：子树叶掩码
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
    """给两棵树（叶集必须相同）返回 (raw, normalized)。"""
    leaves1 = set(t1.get_leaf_names())
    leaves2 = set(t2.get_leaf_names())
    if leaves1 != leaves2:
        raise ValueError("两棵树的叶集不同，无法计算 triplet 距离。")
    # 用稳定顺序构造索引（按字典序排序，确保两树一致）
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
# 4. 单轮采样逻辑 (Worker)
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

    # 缓存“已选”的 pair->mask
    kept_pair_masks: List[Dict[Tuple[int,int], int]] = []

    tries = 0
    early_stop = False
    while len(kept_trees) < num_trees and tries < max_tries:
        tries += 1
        if progress_every and tries % progress_every == 0:
            print(f"    - [进度] 已尝试 {tries} 次，正在寻找...", flush=True)

        if len(visited_topologies) >= max_unique_topologies:
            early_stop = True
            break

        cand = generate_random_tree(ingroup, model=model_type)
        topo_key = canonical_newick_rooted(cand)
        if topo_key in visited_topologies:
            continue
        visited_topologies.add(topo_key)

        # ---- RF 过滤 ----
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

        # ---- Triplet 过滤（内置算法）----
        cand_pairmask, _ = preprocess_pair_lca_masks(cand, leaf_index)

        ok = True
        for pm_kept in kept_pair_masks:
            raw = triplet_distance_raw_from_pairmasks(cand_pairmask, pm_kept, n)
            if raw < triplet_raw_min:
                ok = False
                break
        if not ok:
            continue

        # 通过：加入已选，并缓存其 pair->mask
        kept_trees.append(cand)
        kept_pair_masks.append(cand_pairmask)
        if rf_mode == "strict":
            union_clades |= cand_masks
        else:
            kept_clade_sets.append(cand_masks)

        print(f"    - (尝试 {tries} 次后) 找到第 {len(kept_trees)} 棵树...", flush=True)

    if len(kept_trees) < num_trees:
        if early_stop:
            print(f"    - [信息] 本轮尝试已遍历全部 {max_unique_topologies} 种拓扑，无法找到更多树。")
        else:
            print(f"    - [信息] 本轮尝试已达到最大次数 {max_tries}。")

    return kept_trees

# =============================================================================
# 5. 带重启功能 (Manager)
# =============================================================================

def run_with_restarts(
    ingroup_internal: List[str], outgroup_internal: str, num_trees: int, restarts: int,
    rf_mode: str, rf_tau: float, t_norm_threshold: float, max_tries: int,
    model_type: str, seed: int = None
) -> List[Tree]:

    best: List[Tree] = []

    base_seed = seed if seed is not None else random.randrange(1 << 30)
    print(f"[信息] 使用基础种子: {base_seed} (每次重启时会递增)")

    for attempt in range(restarts + 1):
        print(f"\n--- 第 {attempt + 1}/{restarts + 1} 轮尝试 (使用 {model_type.upper()} 模型, 模式={rf_mode}, RF阈值 {rf_tau}) ---")
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
            print(f"\n[成功] 在第 {attempt + 1} 轮尝试中找到了 {num_trees} 棵树。")
            return best

    print(f"\n[信息] 所有 {restarts + 1} 轮尝试都已完成，未能找齐 {num_trees} 棵树。")
    print(f"将返回找到最多的一组结果 ({len(best)} 棵)。")
    return best

# =============================================================================
# 6. 命令行接口
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "高效采样一组拓扑距离足够远的有根二叉树（随机+RF筛选+Triplet筛选+去重+早停+多次重启）。\n"
            "支持 Uniform / Yule 随机模型；RF 阈值可配；Triplet distance 阈值可配。\n"
            "默认进化模型为yule，进行严格筛选，RF=1, triplet_distance=2/3。\n"
            "默认设置下，对于内群taxa数量n小于等于10时，最多可产生n棵符合要求的树，n大于10时，建议--num_trees 10"
        ),
        epilog = textwrap.dedent("""\
        示例：
        --taxa A B C D E G H I --outgroup F --num_trees 8  --rf-threshold 1   --model uniform   # strict
        --taxa A B C D E G H I --outgroup F --num_trees 7 --rf-threshold 0.8 --model yule      # relax(0.8)
        --taxa "(A,B)" C D "(E,(G,(H,I)))" --outgroup F
        """),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("--taxa", nargs="+", required=True, help="内群的 taxa 列表（原始标签，可含空格/括号）,如果是fixed group, 需要用引号引起来")
    parser.add_argument("--outgroup", required=True, help="用于定根的外群标签（原始标签）")
    parser.add_argument("--num_trees", type=int, default=4, help="希望生成的树的数量 (默认: 4)")

    parser.add_argument("--model", choices=['uniform', 'yule'], default='yule',
                        help='用于生成随机树的模型 ("uniform" 或 "yule", 默认: yule)')
    # <<< 替换：去掉 --mode，新增 --rf-threshold >>>
    parser.add_argument("--rf-threshold", type=float, default=1.0,
                        help="归一化 RF 距离阈值 τ ∈ (0,1]。=1 表示 strict（RFnorm=1.0），"
                             "0<τ<1 表示 relax（RFnorm≥τ）。默认 1.0。")

    parser.add_argument("--t-threshold", type=float, default=2/3,
                        help="归一化 Triplet 距离阈值 (默认: 2/3)")
    parser.add_argument("--max-tries", type=int, default=200000,
                        help="每轮尝试的最大次数 (默认: 200000)")
    parser.add_argument("--restarts", type=int, default=10,
                        help="在未找齐树时，最大随机重启次数 (默认: 10)")
    parser.add_argument("--seed", type=int, default=None,
                        help="用于复现结果的初始随机种子 (默认: None)")

    args = parser.parse_args()

    # <<< 新逻辑：根据 --rf-threshold 决定 rf_mode 与 rf_tau >>>
    if args.rf_threshold is None or args.rf_threshold >= 1.0:
        rf_mode = "strict"
        rf_tau  = 1.0
    else:
        if not (0.0 < args.rf_threshold < 1.0):
            raise ValueError("--rf-threshold 必须在 (0,1]，例如 1 或 0.8")
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
        print(f"\n[错误] 程序执行失败: {e}")
        return

    print("\n--- 最终生成的指导树 ---")
    if not internal_ingroup_trees:
        print("未能生成任何满足条件的树。")
        return

    # 保存“带外群”的内部树以便后面算 tree_1 vs tree_2 距离
    full_internal_trees: List[Tree] = []

    print(f"最终找到 {len(internal_ingroup_trees)} 棵树 (目标是 {args.num_trees} 棵):")
    for i, ingroup_tree in enumerate(internal_ingroup_trees, 1):
        try:
            full_tree_internal = attach_outgroup(ingroup_tree, OUT_INT)  # 用内部占位名
            full_internal_trees.append(full_tree_internal.copy(method="deepcopy"))

            nwk_internal = full_tree_internal.write(format=9)
            nwk_final = replace_placeholders_raw(nwk_internal, mapping)  # 输出前替换回原名
            print(f"\n>tree_{i}\n{nwk_final}")
        except Exception as e:
            print(f"\n[错误] 在输出 tree_{i} 时发生: {e}")

    # ===== 额外：计算并输出 tree_1 vs tree_2 的 triplet 距离（用内置算法） =====
    if len(full_internal_trees) >= 2:
        try:
            t1 = internal_ingroup_trees[0]
            t2 = internal_ingroup_trees[1]
            raw, norm = triplet_distance_raw_between_trees(t1, t2)
            print(f"\n[附加] tree_1 vs tree_2 的 Triplet 距离：raw={raw}, normalized={norm:.12f}")
        except Exception as e:
            print(f"\n[附加] 计算 tree_1 vs tree_2 Triplet 距离失败：{e}")

if __name__ == "__main__":
    main()
