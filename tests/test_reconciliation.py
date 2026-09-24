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

from reference import reference_check_transfer, reference_reconcile

from reconciliation import (
    DUPLICATION,
    SPECIATION,
    TRANSFER_ANCESTOR_RELATED,
    TRANSFER_NO_TIME_OVERLAP,
    TRANSFER_ROOT_BRANCH,
    TRANSFER_SAME_BRANCH,
    InvalidInputError,
    LeafMappingError,
    TreeStructureError,
    build_tree,
    check_transfer,
    reconcile,
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


class TransferHandCheckedTests(unittest.TestCase):
    """候选横向转移检查：手工可复验算例，期望值硬编码。

    物种树 SP：r:(x,D)，x:(y,C)，y:(A,B)。
    深度：r=0；x=1, D=1；y=2, C=2；A=3, B=3。
    分支时间区间 = [depth(父), depth(节点)]（闭区间）。
    """

    @classmethod
    def setUpClass(cls):
        cls.sp = build_tree(SP, "r", name="species")

    def test_feasible_full_overlap_between_siblings(self):
        # A、B 为姐妹分支，区间同为 [2,3]，完全重叠。
        chk = check_transfer(self.sp, "A", "B")
        self.assertTrue(chk.feasible)
        self.assertEqual(chk.interval, (2, 3))
        self.assertIsNone(chk.reason)
        self.assertIsNone(chk.reason_code)

    def test_feasible_partial_overlap(self):
        # 更深的树：p:(u,v)，u:(A,w)，w:(B,C)，v:(D,E)。
        # 深度：p=0；u=1, v=1；A=2, w=2, D=2, E=2；B=3, C=3。
        sp = build_tree(
            {"p": ("u", "v"), "u": ("A", "w"), "w": ("B", "C"), "v": ("D", "E")},
            "p",
        )
        # 供体 w 区间 [1,2]，受体 D 区间 [1,2] -> 完全重叠 (1,2)。
        chk = check_transfer(sp, "w", "D")
        self.assertTrue(chk.feasible)
        self.assertEqual(chk.interval, (1, 2))
        # 供体 B 区间 [2,3]，受体 A 区间 [1,2] -> 仅端点 t=2 相触。
        chk2 = check_transfer(sp, "B", "A")
        self.assertTrue(chk2.feasible)
        self.assertEqual(chk2.interval, (2, 2))

    def test_endpoint_touch_is_feasible(self):
        # C 区间 [1,2]，B 区间 [2,3]：仅在端点 t=2 相触，闭区间视为相交。
        chk = check_transfer(self.sp, "C", "B")
        self.assertTrue(chk.feasible)
        self.assertEqual(chk.interval, (2, 2))
        # 方向对调结论相同（区间相交对称）。
        chk_rev = check_transfer(self.sp, "B", "C")
        self.assertTrue(chk_rev.feasible)
        self.assertEqual(chk_rev.interval, (2, 2))
        # D 区间 [0,1]，C 区间 [1,2]：端点 t=1 相触。
        chk2 = check_transfer(self.sp, "D", "C")
        self.assertTrue(chk2.feasible)
        self.assertEqual(chk2.interval, (1, 1))

    def test_ancestor_relation_rejected_both_directions(self):
        # y 是 A 的祖先（供体为祖先）。
        chk = check_transfer(self.sp, "y", "A")
        self.assertFalse(chk.feasible)
        self.assertIsNone(chk.interval)
        self.assertEqual(chk.reason_code, TRANSFER_ANCESTOR_RELATED)
        self.assertIn("祖先", chk.reason)
        # 受体为祖先（方向对调同样拒绝）。
        chk2 = check_transfer(self.sp, "A", "y")
        self.assertFalse(chk2.feasible)
        self.assertEqual(chk2.reason_code, TRANSFER_ANCESTOR_RELATED)
        # 供体是受体的直接父分支（x -> y）也属于祖先关系。
        chk3 = check_transfer(self.sp, "x", "y")
        self.assertFalse(chk3.feasible)
        self.assertEqual(chk3.reason_code, TRANSFER_ANCESTOR_RELATED)

    def test_ancestor_check_precedes_time_overlap(self):
        # x 是 A 的祖先，且区间 [0,1] 与 [2,3] 本就不相交；
        # 固定检查顺序下应报祖先关系而非时间不相交。
        chk = check_transfer(self.sp, "x", "A")
        self.assertFalse(chk.feasible)
        self.assertEqual(chk.reason_code, TRANSFER_ANCESTOR_RELATED)

    def test_no_time_overlap_rejected(self):
        # D 区间 [0,1]，A 区间 [2,3]：互无祖先关系但时间不相交。
        chk = check_transfer(self.sp, "D", "A")
        self.assertFalse(chk.feasible)
        self.assertIsNone(chk.interval)
        self.assertEqual(chk.reason_code, TRANSFER_NO_TIME_OVERLAP)
        self.assertIn("不相交", chk.reason)

    def test_root_branch_rejected(self):
        for donor, recipient in [("r", "A"), ("A", "r"), ("r", "r")]:
            chk = check_transfer(self.sp, donor, recipient)
            self.assertFalse(chk.feasible, (donor, recipient))
            self.assertEqual(chk.reason_code, TRANSFER_ROOT_BRANCH)
            self.assertIsNone(chk.interval)

    def test_same_branch_rejected(self):
        chk = check_transfer(self.sp, "A", "A")
        self.assertFalse(chk.feasible)
        self.assertEqual(chk.reason_code, TRANSFER_SAME_BRANCH)
        # 同一内节点分支也一样。
        chk2 = check_transfer(self.sp, "y", "y")
        self.assertEqual(chk2.reason_code, TRANSFER_SAME_BRANCH)

    def test_result_is_pure_and_repeatable(self):
        a = check_transfer(self.sp, "C", "B")
        b = check_transfer(self.sp, "C", "B")
        self.assertEqual(a, b)
        # 对合法与非法各取一次，不影响后续复演结果。
        res1 = reconcile(
            self.sp,
            build_tree({"g1": ("a1", "b1"), "gr": ("g1", "c1")}, "gr"),
            {"a1": "A", "b1": "B", "c1": "C"},
        )
        check_transfer(self.sp, "D", "A")
        res2 = reconcile(
            self.sp,
            build_tree({"g1": ("a1", "b1"), "gr": ("g1", "c1")}, "gr"),
            {"a1": "A", "b1": "B", "c1": "C"},
        )
        self.assertEqual(res1, res2)


class TransferInputErrorTests(unittest.TestCase):
    """输入本身不合法（非候选不可行）时必须抛 InvalidInputError。"""

    @classmethod
    def setUpClass(cls):
        cls.sp = build_tree(SP, "r", name="species")

    def test_non_tree_rejected(self):
        with self.assertRaises(InvalidInputError):
            check_transfer("not a tree", "A", "B")
        with self.assertRaises(InvalidInputError):
            check_transfer(None, "A", "B")

    def test_unknown_node_rejected(self):
        with self.assertRaises(InvalidInputError):
            check_transfer(self.sp, "ZZZ", "A")
        with self.assertRaises(InvalidInputError):
            check_transfer(self.sp, "A", "ZZZ")

    def test_bad_label_types(self):
        for bad in (3, None, "", True, 1.5):
            with self.assertRaises(InvalidInputError):
                check_transfer(self.sp, bad, "A")
            with self.assertRaises(InvalidInputError):
                check_transfer(self.sp, "A", bad)


class TransferCrossCheckTests(unittest.TestCase):
    """与独立参考的穷举交叉核对：3/4 叶物种树全部形状 × 全部有序节点对。"""

    def _run_family(self, sp_leaves):
        k = len(sp_leaves)
        checked = 0
        for sp_shape in all_shapes(0, k):
            sp_ch, sp_root = to_children(sp_shape, "s")
            sp_ch = {
                p: tuple(
                    sp_leaves[int(c[4:])] if c.startswith("leaf") else c
                    for c in pair
                )
                for p, pair in sp_ch.items()
            }
            sp_tree = build_tree(sp_ch, sp_root, name="species")
            nodes = sorted(sp_tree.nodes)
            for donor in nodes:
                for recipient in nodes:
                    chk = check_transfer(sp_tree, donor, recipient)
                    ref = reference_check_transfer(sp_ch, sp_root, donor, recipient)
                    self.assertEqual(
                        (chk.feasible, chk.interval, chk.reason_code),
                        ref,
                        (sp_leaves, sp_shape, donor, recipient),
                    )
                    # 聚合自洽（补充）：可行 <=> 有区间且无原因码。
                    self.assertEqual(chk.feasible, chk.interval is not None)
                    self.assertEqual(chk.feasible, chk.reason_code is None)
                    self.assertEqual(chk.feasible, chk.reason is None)
                    checked += 1
        return checked

    def test_species_3_leaves_all_pairs(self):
        n = self._run_family(("A", "B", "C"))
        self.assertGreater(n, 40)  # 2 形状 × 5 节点² = 50

    def test_species_4_leaves_all_pairs(self):
        n = self._run_family(("A", "B", "C", "D"))
        self.assertGreater(n, 200)  # 5 形状 × 7 节点² = 245


if __name__ == "__main__":
    unittest.main(verbosity=2)
