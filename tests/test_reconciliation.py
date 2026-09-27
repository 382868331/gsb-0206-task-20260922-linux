"""reconciliation 核心测试。

期望值来源：
1. 手工算例硬编码（speciation 树、复制树、单支丢失、多次复制）；
2. tests/reference.py 中不调用被测核心的独立参考实现，在一批
   确定性穷举的小树（叶标签允许重复，即一物种多基因）上交叉核对；
3. 往返/聚合自洽仅作补充断言。
"""

from __future__ import annotations

import math
import unittest
from itertools import product

from reference import reference_reconcile
from reference import reference_best_candidates

from reconciliation import (
    BUDGET_EXHAUSTED,
    DUPLICATION,
    OPTIMAL_CONFIRMED,
    SPECIATION,
    InvalidInputError,
    LeafMappingError,
    TreeStructureError,
    build_tree,
    reconcile,
    reconcile_candidates,
)

# ---------- 常用物种树 ----------
#  y:(A,B), x:(y,C), r:(x,D)
SP = {"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}
#  x:(A,B), r:(x,C)
SP3 = {"x": ("A", "B"), "r": ("x", "C")}


def rename_children(children, rename):
    return {
        rename.get(p, p): tuple(rename.get(c, c) for c in pair)
        for p, pair in children.items()
    }


class HandCheckedTests(unittest.TestCase):
    """手工可复验的小算例，期望值硬编码。"""

    def test_pure_speciation_tree(self):
        """基因树与物种树同构：全部 speciation，零丢失。"""
        g = {
            "gy": ("a", "b"),
            "gx": ("gy", "c"),
            "gr": ("gx", "d"),
        }
        sp = build_tree(SP, "r", name="species")
        gn = build_tree(g, "gr", name="gene")
        res = reconcile(sp, gn, {"a": "A", "b": "B", "c": "C", "d": "D"})

        self.assertEqual(res.mapping["a"], "A")
        self.assertEqual(res.mapping["gy"], "y")
        self.assertEqual(res.mapping["gx"], "x")
        self.assertEqual(res.mapping["gr"], "r")
        self.assertTrue(all(e.event_type == SPECIATION for e in res.events))
        self.assertEqual(res.speciations, 3)
        self.assertEqual(res.duplications, 0)
        self.assertEqual(res.total_losses, 0)
        self.assertTrue(all(e.loss_count == 0 for e in res.edge_losses))
        self.assertEqual(res.losses_by_lineage, {})

    def test_single_lineage_loss(self):
        """基因 (a,c)：路径 r->x->A 上只有 B 支系缺失，恰 1 次丢失。"""
        sp = build_tree(SP3, "r")
        gn = build_tree({"gr": ("a", "c")}, "gr")
        res = reconcile(sp, gn, {"a": "A", "c": "C"})

        self.assertEqual(res.duplications, 0)
        self.assertEqual(res.speciations, 1)
        self.assertEqual(res.total_losses, 1)
        by_edge = {(e.parent_gene, e.child_gene): e for e in res.edge_losses}
        e_ac = by_edge[("gr", "a")]
        self.assertEqual((e_ac.path_length, e_ac.loss_count), (2, 1))
        self.assertEqual([l.species_node for l in e_ac.lost_lineages], ["B"])
        self.assertEqual(e_ac.lost_lineages[0].descendant_leaves, frozenset({"B"}))
        self.assertEqual(by_edge[("gr", "c")].loss_count, 0)
        self.assertEqual(res.losses_by_lineage, {"B": 1})

    def test_one_duplication_with_losses(self):
        """根复制：gr 子支 g2 与 gr 同映到 x；复制边丢掉整个 C 支系。"""
        #  g1:(a1,b1)->y ; g2:(a2,c1)->x ; gr:(g1,g2)->x=dup
        g = {
            "g1": ("a1", "b1"),
            "g2": ("a2", "c1"),
            "gr": ("g1", "g2"),
        }
        lmap = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        sp = build_tree(SP, "r")
        gn = build_tree(g, "gr")
        res = reconcile(sp, gn, lmap)

        ev = {e.gene_node: e for e in res.events}
        self.assertEqual(ev["g1"].event_type, SPECIATION)
        self.assertEqual(ev["g2"].event_type, SPECIATION)
        self.assertEqual(ev["gr"].event_type, DUPLICATION)
        self.assertEqual((res.duplications, res.speciations), (1, 2))

        by_edge = {(e.parent_gene, e.child_gene): e for e in res.edge_losses}
        # 复制父 gr：x->y 路径长 1，不减 1 => 丢 1（兄弟 C 整支）
        self.assertEqual(by_edge[("gr", "g1")].loss_count, 1)
        self.assertEqual(
            [l.species_node for l in by_edge[("gr", "g1")].lost_lineages], ["C"]
        )
        self.assertEqual(by_edge[("gr", "g1")].lost_lineages[0].descendant_leaves,
                         frozenset({"C"}))
        # x->x 路径长 0，复制也只丢 0
        self.assertEqual(by_edge[("gr", "g2")].loss_count, 0)
        # g2 为 speciation：x->A 路径长 2，丢 1（跳过起点兄弟 C，记 B）
        e_a2 = by_edge[("g2", "a2")]
        self.assertEqual((e_a2.path_length, e_a2.loss_count), (2, 1))
        self.assertEqual([l.species_node for l in e_a2.lost_lineages], ["B"])
        # 其余边为 0
        for key in [("g1", "a1"), ("g1", "b1"), ("g2", "c1")]:
            self.assertEqual(by_edge[key].loss_count, 0, key)
        # C、B 各丢一次
        self.assertEqual(res.losses_by_lineage, {"B": 1, "C": 1})
        self.assertEqual(res.total_losses, 2)

    def test_multiple_duplications(self):
        """两次复制：种内复制 p:(a1,a2) 与上层复制 h:(p,q)。"""
        # p:(a1,a2)->A dup ; q:(a3,b1)->y spec ; h:(p,q)->y, q 映 y => dup
        g = {
            "p": ("a1", "a2"),
            "q": ("a3", "b1"),
            "h": ("p", "q"),
        }
        lmap = {"a1": "A", "a2": "A", "a3": "A", "b1": "B"}
        sp = build_tree({"y": ("A", "B"), "r": ("y", "C")}, "r")
        gn = build_tree(g, "h")
        res = reconcile(sp, gn, lmap)

        ev = {e.gene_node: e for e in res.events}
        self.assertEqual(ev["p"].event_type, DUPLICATION)
        self.assertEqual(ev["q"].event_type, SPECIATION)
        self.assertEqual(ev["h"].event_type, DUPLICATION)
        self.assertEqual((res.duplications, res.speciations), (2, 1))
        by_edge = {(e.parent_gene, e.child_gene): e for e in res.edge_losses}
        self.assertEqual(by_edge[("p", "a1")].loss_count, 0)
        self.assertEqual(by_edge[("p", "a2")].loss_count, 0)
        self.assertEqual(by_edge[("h", "p")].loss_count, 1)
        self.assertEqual(
            [l.species_node for l in by_edge[("h", "p")].lost_lineages], ["B"]
        )
        self.assertEqual(by_edge[("h", "q")].loss_count, 0)
        self.assertEqual(res.total_losses, 1)


class PermutationInvarianceTests(unittest.TestCase):
    def test_leaf_permutation_invariance(self):
        """基因叶重命名（双射置换）+ 内节点改名后，归因结果随改名一致。"""
        g = {
            "g1": ("a1", "b1"),
            "g2": ("a2", "c1"),
            "gr": ("g1", "g2"),
        }
        lmap = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        sp = build_tree(SP, "r")
        res1 = reconcile(sp, build_tree(g, "gr"), lmap)

        # 叶子置换 + 内节点改名（纯改名，拓扑不变；左/右顺序保持）。
        rename = {"a1": "z9", "b1": "z2", "a2": "z5", "c1": "z0",
                  "g1": "P", "g2": "Q", "gr": "R"}
        g2 = rename_children(g, rename)
        lmap2 = {rename[k]: v for k, v in lmap.items()}
        res2 = reconcile(sp, build_tree(g2, "R"), lmap2)

        self.assertEqual(
            (res1.duplications, res1.speciations, res1.total_losses),
            (res2.duplications, res2.speciations, res2.total_losses),
        )
        self.assertEqual(res1.losses_by_lineage, res2.losses_by_lineage)
        # 映射随基因节点改名一致
        for old, new in rename.items():
            self.assertEqual(res1.mapping[old], res2.mapping[new])
        # 每个基因内节点的事件类型随改名一致
        ev1 = {e.gene_node: e.event_type for e in res1.events}
        ev2 = {e.gene_node: e.event_type for e in res2.events}
        for old, new in rename.items():
            if old in ev1:
                self.assertEqual(ev1[old], ev2[new])

    def test_left_right_swap_invariance(self):
        """交换每个内节点的左右子（同时整体镜像）不改变计数与支系归因。"""
        g = {
            "g1": ("a1", "b1"),
            "g2": ("c1", "a2"),
            "gr": ("g2", "g1"),  # 交换左右
        }
        # g1 自身也镜像
        g["g1"] = ("b1", "a1")
        lmap = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        sp = build_tree(SP, "r")
        gn = build_tree(g, "gr")
        res = reconcile(sp, gn, lmap)
        self.assertEqual((res.duplications, res.speciations, res.total_losses), (1, 2, 2))
        self.assertEqual(res.losses_by_lineage, {"B": 1, "C": 1})


# ---------- 与独立参考的穷举交叉核对 ----------

def all_shapes(lo, hi):
    """位置区间 [lo,hi) 上全部满二叉形状（叶子为整数位置）。"""
    if hi - lo == 1:
        yield lo
        return
    for k in range(lo + 1, hi):
        for left in all_shapes(lo, k):
            for right in all_shapes(k, hi):
                yield (left, right)


def to_children(shape, prefix):
    """把嵌套 tuple 形状转成 children dict，返回 (children, root)。"""
    counter = [0]
    children = {}

    def build(node):
        if isinstance(node, int):
            return f"leaf{node}"
        left, right = node
        lc, rc = build(left), build(right)
        name = f"{prefix}n{counter[0]}"
        counter[0] += 1
        children[name] = (lc, rc)
        return name

    root = build(shape)
    return children, root


class CrossCheckReferenceTests(unittest.TestCase):
    """独立参考交叉核对：3/4 叶物种树全部形状 × 2..4 基因叶全部标记与形状。"""

    def _run_family(self, sp_leaves):
        k = len(sp_leaves)
        checked = 0
        for sp_shape in all_shapes(0, k):
            sp_ch, sp_root = to_children(sp_shape, "s")
            # 用物种叶位置重映射到真实标签
            sp_ch = {
                p: tuple(
                    sp_leaves[int(c[4:])] if c.startswith("leaf") else c
                    for c in pair
                )
                for p, pair in sp_ch.items()
            }
            sp_tree = build_tree(sp_ch, sp_root, name="species")
            for n in (2, 3, 4):
                for shape in all_shapes(0, n):
                    g_ch, g_root = to_children(shape, "g")
                    leaf_names = [f"leaf{i}" for i in range(n)]
                    for labels in product(sp_leaves, repeat=n):
                        lmap = dict(zip(leaf_names, labels))
                        g_tree = build_tree(g_ch, g_root, name="gene")
                        res = reconcile(sp_tree, g_tree, lmap)
                        rm, rev, reloss, resibs, rtotal = reference_reconcile(
                            sp_ch, sp_root, g_ch, g_root, lmap
                        )
                        # 映射一致
                        self.assertEqual(dict(res.mapping), rm)
                        # 事件一致
                        core_ev = {e.gene_node: e.event_type for e in res.events}
                        self.assertEqual(core_ev, rev)
                        # 逐边计数与缺失支系顺序一致
                        for edge in res.edge_losses:
                            key = (edge.parent_gene, edge.child_gene)
                            self.assertEqual(edge.loss_count, reloss[key], (key, labels))
                            self.assertEqual(
                                [l.species_node for l in edge.lost_lineages],
                                resibs[key],
                                (key, labels),
                            )
                        self.assertEqual(res.total_losses, rtotal)
                        # 聚合自洽（补充）
                        self.assertEqual(
                            res.total_losses,
                            sum(e.loss_count for e in res.edge_losses),
                        )
                        self.assertEqual(
                            res.total_losses,
                            sum(res.losses_by_lineage.values()),
                        )
                        checked += 1
        return checked

    def test_species_3_leaves(self):
        n = self._run_family(("A", "B", "C"))
        self.assertGreater(n, 100)

    def test_species_4_leaves(self):
        n = self._run_family(("A", "B", "C", "D"))
        self.assertGreater(n, 5000)


class BoundaryAndStructureTests(unittest.TestCase):
    def test_trivial_single_leaf_trees(self):
        sp = build_tree({}, "A", name="species")
        gn = build_tree({}, "a", name="gene")
        res = reconcile(sp, gn, {"a": "A"})
        self.assertEqual(res.mapping, {"a": "A"})
        self.assertEqual(res.events, ())
        self.assertEqual(res.edge_losses, ())
        self.assertEqual(res.total_losses, 0)

    def test_node_limit_accept_199_reject_201(self):
        # 满二叉树节点数必为奇数；100 叶 => 199 节点（<=200 接受）。
        def comb(k):
            ch = {}
            prev = "L1"
            for i in range(2, k + 1):
                node = f"u{i}"
                ch[node] = (prev, f"L{i}")
                prev = node
            return ch, prev

        ch100, root100 = comb(100)
        sp = build_tree(ch100, root100, name="species")
        gn = build_tree(ch100, root100, name="gene")
        lmap = {f"L{i}": f"L{i}" for i in range(1, 101)}
        res = reconcile(sp, gn, lmap)
        self.assertEqual(res.total_losses, 0)
        self.assertEqual(res.duplications, 0)

        ch101, root101 = comb(101)  # 201 节点
        with self.assertRaises(InvalidInputError):
            build_tree(ch101, root101)

    def test_cyclic_tree_rejected(self):
        # a -> b -> a 的回边
        with self.assertRaises(TreeStructureError) as cm:
            build_tree({"a": ("b", "c"), "b": ("d", "a")}, "a", name="gene")
        self.assertEqual(cm.exception.tree, "gene")
        self.assertIsNotNone(cm.exception.node)

    def test_cyclic_species_tree_rejected(self):
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("x", "A"), "x": ("r", "B")}, "r", name="species")

    def test_multiple_parents_rejected(self):
        with self.assertRaises(TreeStructureError) as cm:
            build_tree({"r": ("a", "b"), "x": ("a", "c")}, "r")
        self.assertEqual(cm.exception.node, "a")

    def test_unreachable_component_rejected(self):
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a", "b"), "x": ("c", "d")}, "r")

    def test_non_binary_rejected(self):
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a",)}, "r")
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a", "b", "c")}, "r")

    def test_equal_children_and_self_loop(self):
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a", "a")}, "r")
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("r", "a")}, "r")

    def test_ambiguous_and_unknown_root(self):
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a", "b"), "x": ("c", "d")})  # 两个无父内节点
        with self.assertRaises(TreeStructureError):
            build_tree({"r": ("a", "b")}, "zzz")

    def test_empty_children_requires_root(self):
        with self.assertRaises(InvalidInputError):
            build_tree({})

    def test_bad_label_types(self):
        with self.assertRaises(InvalidInputError):
            build_tree({1: ("a", "b")}, 1)
        with self.assertRaises(InvalidInputError):
            build_tree({"r": ("a", "")}, "r")
        with self.assertRaises(InvalidInputError):
            build_tree([("r", ("a", "b"))], "r")

    def test_numeric_input_rejects_bool_nan_inf(self):
        good = {"r": ("a", "b")}
        with self.assertRaises(InvalidInputError):
            build_tree(good, "r", max_nodes=True)   # bool 不是整数参数
        with self.assertRaises(InvalidInputError):
            build_tree(good, "r", max_nodes=False)
        for bad in (float("nan"), float("inf"), float("-inf"), 1.5, "200", None):
            with self.assertRaises(InvalidInputError):
                build_tree(good, "r", max_nodes=bad)
        for bad in (0, -1):
            with self.assertRaises(InvalidInputError):
                build_tree(good, "r", max_nodes=bad)
        with self.assertRaises(InvalidInputError):
            build_tree(good, "r", max_nodes=201)    # 超过题面 200 上限
        sp = build_tree(good, "r")
        with self.assertRaises(InvalidInputError):
            reconcile(sp, sp, {"a": "a", "b": "b"}, max_nodes=True)
        # math.isnan / isinf 覆盖的显式语义检查
        self.assertTrue(math.isnan(float("nan")))
        self.assertTrue(math.isinf(float("inf")))


class MappingErrorTests(unittest.TestCase):
    def setUp(self):
        self.sp = build_tree(SP3, "r")
        self.gn = build_tree({"gr": ("a", "b")}, "gr")

    def test_non_leaf_mapping_rejected(self):
        with self.assertRaises(LeafMappingError) as cm:
            reconcile(self.sp, self.gn, {"a": "x", "b": "B"})  # x 是物种内节点
        self.assertEqual(cm.exception.gene_leaf, "a")

    def test_missing_leaf(self):
        with self.assertRaises(LeafMappingError):
            reconcile(self.sp, self.gn, {"a": "A"})

    def test_unknown_species_target(self):
        with self.assertRaises(LeafMappingError) as cm:
            reconcile(self.sp, self.gn, {"a": "A", "b": "ZZZ"})
        self.assertEqual(cm.exception.gene_leaf, "b")

    def test_key_not_gene_leaf(self):
        # 内节点做键
        with self.assertRaises(LeafMappingError):
            reconcile(self.sp, self.gn, {"a": "A", "b": "B", "gr": "A"})
        # 根本不在基因树里
        with self.assertRaises(LeafMappingError):
            reconcile(self.sp, self.gn, {"a": "A", "b": "B", "zzz": "A"})

    def test_bad_mapping_types(self):
        sp = build_tree(SP3, "r")
        with self.assertRaises(InvalidInputError):
            reconcile(sp, sp, [("a", "A")])
        with self.assertRaises(InvalidInputError):
            reconcile("not a tree", sp, {"a": "A"})
        with self.assertRaises(InvalidInputError):
            reconcile(sp, object(), {"a": "A"})
        with self.assertRaises(LeafMappingError):
            reconcile(self.sp, self.gn, {"a": "A", "b": 3})
        with self.assertRaises(LeafMappingError):
            reconcile(self.sp, self.gn, {"a": "A", 7: "B"})


class EventPathReplayTests(unittest.TestCase):
    """补充：用结果中的映射与子节点关系回放，事件路径可独立复验。"""

    def test_replay_from_result(self):
        g = {
            "g1": ("a1", "b1"),
            "g2": ("a2", "c1"),
            "gr": ("g1", "g2"),
        }
        lmap = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        sp = build_tree(SP, "r")
        gn = build_tree(g, "gr")
        res = reconcile(sp, gn, lmap)

        # 不依赖核心的事件判定，仅用 mapping + 基因树结构复验事件类型
        for ev in res.events:
            c0, c1 = ev.child_genes
            expect_dup = (
                res.mapping[c0] == ev.mapped_species
                or res.mapping[c1] == ev.mapped_species
            )
            self.assertEqual(ev.event_type, DUPLICATION if expect_dup else SPECIATION)

        # 每条边：路径长度可由两物种节点深度关系复验（LCA 映射保证祖先关系）
        for edge in res.edge_losses:
            # 独立爬升求路径长
            x, steps = edge.child_species, 0
            while x != edge.parent_species:
                x = sp.parent[x]
                steps += 1
            self.assertEqual(steps, edge.path_length)
            expect_losses = (
                edge.path_length if edge.parent_gene in
                {e.gene_node for e in res.events if e.event_type == DUPLICATION}
                else edge.path_length - 1
            )
            self.assertEqual(edge.loss_count, expect_losses)
            self.assertEqual(edge.loss_count, len(edge.lost_lineages))


class CandidateSelectionTests(unittest.TestCase):
    """候选叶选择：手工算例、独立穷举交叉核对、并列顺序、额度、错误语义。"""

    def test_hand_checked_candidate_choice(self):
        """物种树 A,B|C；基因 (x,y)：x 候选 A/B，y=C 固定。

        x=A: (A,C) 为分化，路径经根，丢 B 支系 1 次 => 代价 0D+1L=1。
        x=B: 同理丢 A 支系 => 代价 1。两个解并列，按选定元组 (A,C)<(B,C) 排序。
        """
        sp = build_tree({"r": ("x", "C"), "x": ("A", "B")}, "r")
        gn = build_tree({"gr": ("gx", "gy")}, "gr")
        out = reconcile_candidates(
            sp, gn, {"gx": ["B", "A"], "gy": "C"}
        )
        self.assertEqual(out.status, OPTIMAL_CONFIRMED)
        self.assertEqual(out.best_cost, 1)
        self.assertEqual(out.evaluated, 2)
        self.assertEqual(out.total_combinations, 2)
        chosen = [tuple(s.leaf_mapping[g] for g in sorted(s.leaf_mapping))
                  for s in out.solutions]
        # 输入候选顺序 [B,A] 不影响结果；按物种标签稳定排列。
        self.assertEqual(chosen, [("A", "C"), ("B", "C")])
        for sol in out.solutions:
            self.assertEqual(sol.result.duplications, 0)
            self.assertEqual(sol.result.total_losses, 1)
            self.assertEqual((sol.duplication_cost, sol.loss_cost), (0, 1))

    def test_weights_passthrough_and_zero_cost_optimum(self):
        """权重由调用方透传；零复制零丢失映射在任何正权重下都最优。"""
        # 物种树 (((A,B),C),D) 同 SP；基因树 ((a,b),c)，a 固定 A，b 固定 B，
        # c 候选 C 或 D：c=C 时 (y,C) 分化无复制，边到 y 丢 0；c=D 丢更多。
        sp = build_tree(SP, "r")
        gn = build_tree({"gy": ("a", "b"), "gr": ("gy", "c")}, "gr")
        candidates = {"a": "A", "b": "B", "c": ["C", "D"]}
        out = reconcile_candidates(sp, gn, candidates)
        self.assertEqual(out.best_cost, 0)
        self.assertEqual([s.leaf_mapping["c"] for s in out.solutions], ["C"])
        # 加权后各解代价仍单调：C 解为 0，D 解必为正，结论不变且权重透传。
        out_w = reconcile_candidates(
            sp, gn, candidates, duplication_weight=3, loss_weight=5
        )
        self.assertEqual(out_w.best_cost, 0)
        self.assertEqual(out_w.duplication_weight, 3)
        self.assertEqual(out_w.loss_weight, 5)

    def test_duplication_vs_loss_tradeoff_weighted(self):
        """构造复制/丢失取舍：不同权重选出不同最优。"""
        # 物种树 (A,B) 根 r；基因树 gr:(a1,a2)，两叶均可映 A 或 B。
        sp = build_tree({"r": ("A", "B")}, "r")
        gn = build_tree({"gr": ("a1", "a2")}, "gr")
        candidates = {"a1": ["A", "B"], "a2": ["A", "B"]}
        # (A,B)/(B,A)：分化 0D0L，代价 0，恒最优；(A,A)/(B,B)：复制 1D0L。
        out = reconcile_candidates(sp, gn, candidates)
        self.assertEqual(
            [tuple(sorted(s.leaf_mapping.items())) for s in out.solutions],
            sorted([(("a1", "A"), ("a2", "B")),
                    (("a1", "B"), ("a2", "A"))]),
        )
        for sol in out.solutions:
            self.assertEqual(sol.result.duplications, 0)

    def test_all_fixed_matches_reconcile(self):
        """全部叶归属明确（str 候选）时，与直接 reconcile 的结果完全一致。"""
        g = {"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")}
        lmap = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        sp = build_tree(SP, "r")
        gn = build_tree(g, "gr")
        direct = reconcile(sp, gn, lmap)
        out = reconcile_candidates(sp, gn, lmap)
        self.assertEqual(out.status, OPTIMAL_CONFIRMED)
        self.assertEqual(out.total_combinations, 1)
        self.assertEqual(len(out.solutions), 1)
        sol = out.solutions[0]
        self.assertEqual(sol.leaf_mapping, lmap)
        self.assertEqual(sol.result.mapping, direct.mapping)
        self.assertEqual(sol.result.events, direct.events)
        self.assertEqual(sol.result.edge_losses, direct.edge_losses)
        self.assertEqual(sol.total_cost,
                         direct.duplications + direct.total_losses)

    def test_brute_force_reference_cross_check(self):
        """小树穷举：核心选优与独立参考（不调用核心）逐解核对代价与并列顺序。"""
        sp_children = {"x": ("A", "B"), "r": ("x", "C")}
        sp = build_tree(sp_children, "r")
        g_children = {"g1": ("a", "b"), "gr": ("g1", "c")}
        gn = build_tree(g_children, "gr")
        candidates = {
            "a": ["A", "B", "C"],
            "b": ["A", "C"],
            "c": ["B", "C"],
        }
        for dup_w, loss_w in ((1, 1), (1, 4), (3, 2)):
            out = reconcile_candidates(
                sp, gn, candidates,
                duplication_weight=dup_w, loss_weight=loss_w,
            )
            ref_cost, ref_tied = reference_best_candidates(
                sp_children, "r", g_children, "gr", candidates,
                dup_w, loss_w,
            )
            self.assertEqual(out.best_cost, ref_cost)
            self.assertEqual(len(out.solutions), len(ref_tied))
            got = [
                (tuple(s.leaf_mapping[leaf] for leaf in sorted(s.leaf_mapping)),
                 s.result.duplications, s.result.total_losses)
                for s in out.solutions
            ]
            self.assertEqual(got, ref_tied)
            self.assertEqual(out.status, OPTIMAL_CONFIRMED)
            self.assertEqual(out.evaluated,
                             3 * 2 * 2)

    def test_search_budget_exhausted(self):
        """额度不足：只评估预算份数，状态标记尚未确认最优。"""
        sp = build_tree({"r": ("A", "B")}, "r")
        gn = build_tree({"gr": ("a1", "a2")}, "gr")
        candidates = {"a1": ["A", "B"], "a2": ["A", "B"]}
        out = reconcile_candidates(sp, gn, candidates, search_budget=2)
        self.assertEqual(out.status, BUDGET_EXHAUSTED)
        self.assertFalse(out.optimal_confirmed)
        self.assertEqual(out.evaluated, 2)
        self.assertEqual(out.total_combinations, 4)
        # 额度 1：唯一评估的解就是当前返回解。
        out1 = reconcile_candidates(sp, gn, candidates, search_budget=1)
        self.assertEqual(out1.status, BUDGET_EXHAUSTED)
        self.assertEqual(len(out1.solutions), 1)
        # 不限额度时确认最优。
        full = reconcile_candidates(sp, gn, candidates, search_budget=None)
        self.assertEqual(full.status, OPTIMAL_CONFIRMED)
        self.assertTrue(full.optimal_confirmed)

    def test_budget_exactly_covers_space(self):
        """额度恰好等于组合数时视为已确认最优。"""
        sp = build_tree({"r": ("A", "B")}, "r")
        gn = build_tree({"gr": ("a1", "a2")}, "gr")
        out = reconcile_candidates(
            sp, gn, {"a1": ["A", "B"], "a2": ["A", "B"]}, search_budget=4
        )
        self.assertEqual(out.status, OPTIMAL_CONFIRMED)

    def test_candidate_validation_errors(self):
        """候选输入错误全部在搜索前拒绝。"""
        sp = build_tree({"r": ("A", "B")}, "r")
        gn = build_tree({"gr": ("a", "b")}, "gr")
        # 缺叶
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": "A"})
        # 多余/非基因叶键
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B", "gr": "A"})
        # 空候选
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": [], "b": "B"})
        # 候选重复
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": ["A", "A"], "b": "B"})
        # 非叶目标（映到物种内节点）
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": ["A", "r"], "b": "B"})
        # 未知物种
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": ["A", "Z"], "b": "B"})
        # 候选类型非法
        with self.assertRaises(LeafMappingError):
            reconcile_candidates(sp, gn, {"a": 3, "b": "B"})
        # 权重非法
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B"},
                                 duplication_weight=0)
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B"},
                                 loss_weight=True)
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B"},
                                 loss_weight=float("nan"))
        # 额度非法
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B"},
                                 search_budget=0)
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, {"a": "A", "b": "B"},
                                 search_budget=1.5)
        # 非 Mapping 入参
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, gn, [("a", "A"), ("b", "B")])
        # 错树类型
        with self.assertRaises(InvalidInputError):
            reconcile_candidates(sp, object(), {"a": "A", "b": "B"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
