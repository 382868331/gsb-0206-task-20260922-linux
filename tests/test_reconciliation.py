"""必要测试。

期望值全部手工推导（见各用例注释），不调用被测核心产生。
``_reference_reconcile`` 是用不同算法（祖先集合求 LCA、直接按题面公式计数）
写的独立参考，只做交叉核对（补充手段）；``Tree.lca`` / ``Tree.path`` /
``reconcile`` 均未被参考实现调用。
"""

import math
import unittest

from reconciliation import (
    EVENT_DUPLICATION,
    EVENT_SPECIATION,
    MAX_NODES,
    MappingError,
    Node,
    Reconciliation,
    Tree,
    TreeError,
    newick,
    node_from_spec,
    reconcile,
)


# ---------------------------------------------------------------------------
# 独立参考实现：直接在嵌套元组 spec 上工作，键为标签（标签全树唯一），
# 用祖先集合求 LCA——与被测实现的深度爬升法完全不同。
# ---------------------------------------------------------------------------
def _reference_reconcile(species_spec, gene_spec, leaf_map):
    def parse(spec):
        if isinstance(spec, str):
            return spec, None, None
        if len(spec) == 1:
            return spec[0], None, None
        label, left, right = spec
        return label, parse(left), parse(right)

    def index(root):
        parent, kids, leaves, post = {}, {}, [], []

        def walk(t, p):
            label, left, right = t
            parent[label] = p
            if left is None:
                kids[label] = ()
                leaves.append(label)
            else:
                kids[label] = (left[0], right[0])
                walk(left, label)
                walk(right, label)
            post.append(label)

        walk(root, None)
        return parent, kids, leaves, post

    sp, sk, sl, _ = index(parse(species_spec))
    gp, gk, gl, gpost = index(parse(gene_spec))

    def ancestors(x, parent):
        out = set()
        while x is not None:
            out.add(x)
            x = parent[x]
        return out

    def lca(a, b):
        aa = ancestors(a, sp)
        x = b
        while x not in aa:
            x = sp[x]
        return x

    def path(anc, desc):
        chain = [desc]
        while chain[-1] != anc:
            chain.append(sp[chain[-1]])
        chain.reverse()
        return chain

    gmap = {}
    for g in gpost:
        if g in gl:
            gmap[g] = leaf_map[g]
        else:
            gmap[g] = lca(gmap[gk[g][0]], gmap[gk[g][1]])

    events, counts, nd, ns, nl = [], {}, 0, 0, 0
    for g in gpost:
        if g in gl:
            continue
        m = gmap[g]
        is_dup = gmap[gk[g][0]] == m or gmap[gk[g][1]] == m
        nd += is_dup
        ns += not is_dup
        for c in gk[g]:
            p = path(m, gmap[c])
            k = len(p) - 1
            n = k if is_dup else k - 1  # 题面公式：D 不减 1，S 路径长度减 1
            counts[(g, c)] = n
            nl += n
        events.append((g, m, EVENT_DUPLICATION if is_dup else EVENT_SPECIATION))
    return gmap, events, counts, ns, nd, nl


class _SpecCase(unittest.TestCase):
    """对 spec 跑被测核心，并用独立参考交叉核对（补充断言）。"""

    def assertRefAgrees(self, species_spec, gene_spec, leaf_map):
        result = reconcile(
            Tree(node_from_spec(species_spec), name="species"),
            Tree(node_from_spec(gene_spec), name="gene"),
            leaf_map,
        )
        gmap, events, counts, ns, nd, nl = _reference_reconcile(
            species_spec, gene_spec, leaf_map
        )
        self.assertEqual(result.mapping, gmap)
        self.assertEqual(
            sorted((e.gene, e.species, e.kind) for e in result.events),
            sorted(events),
        )
        self.assertEqual(result.edge_loss_counts, counts)
        self.assertEqual(result.n_speciations, ns)
        self.assertEqual(result.n_duplications, nd)
        self.assertEqual(result.n_losses, nl)
        return result


# 物种树：R -> (X -> (Y -> (A,B), Z -> (C,D)), E)
SPECIES5 = ("R", ("X", ("Y", "A", "B"), ("Z", "C", "D")), "E")
# 三物种小树：R -> (Y -> (A,B), C)
SPECIES3 = ("R", ("Y", "A", "B"), "C")


class TestHandComputed(_SpecCase):
    def test_pure_speciation_no_loss(self):
        # ((a,b),(c,d))：g1=LCA(A,B)=Y，g2=LCA(C,D)=Z，g0=LCA(Y,Z)=X
        gene = ("g0", ("g1", "a", "b"), ("g2", "c", "d"))
        lmap = dict(a="A", b="B", c="C", d="D")
        r = self.assertRefAgrees(SPECIES5, gene, lmap)

        self.assertEqual(
            r.mapping, dict(a="A", b="B", c="C", d="D", g1="Y", g2="Z", g0="X")
        )
        self.assertEqual(r.n_speciations, 3)
        self.assertEqual(r.n_duplications, 0)
        self.assertEqual(r.n_losses, 0)
        self.assertTrue(all(e.kind == EVENT_SPECIATION for e in r.events))
        self.assertTrue(all(v == 0 for v in r.edge_loss_counts.values()))
        # 事件表的子映射与逐边丢失
        ev = {e.gene: e for e in r.events}
        self.assertEqual(
            ev["g0"].children, (("g1", "Y"), ("g2", "Z"))
        )
        self.assertEqual(
            ev["g0"].edge_losses,
            (("g1", ()), ("g2", ())),
        )

    def test_one_species_multiple_genes_duplication(self):
        # 一物种多基因：((a1,a2),b)，A 内两个基因
        gene = ("g0", ("g1", "a1", "a2"), "b")
        lmap = dict(a1="A", a2="A", b="B")
        r = self.assertRefAgrees(SPECIES3, gene, lmap)

        self.assertEqual(r.mapping["g1"], "A")  # LCA(A,A)=A
        self.assertEqual(r.mapping["g0"], "Y")  # LCA(A,B)=Y
        ev = {e.gene: e for e in r.events}
        self.assertEqual(ev["g1"].kind, EVENT_DUPLICATION)
        self.assertEqual(ev["g0"].kind, EVENT_SPECIATION)
        self.assertEqual(r.n_duplications, 1)
        self.assertEqual(r.n_speciations, 1)
        self.assertEqual(r.n_losses, 0)

    def test_single_lineage_loss(self):
        # (a1,c)：g0=LCA(A,C)=R；R->Y->A 路径上 A 的兄弟 B 是唯一丢失
        gene = ("g0", "a1", "c")
        lmap = dict(a1="A", c="C")
        r = self.assertRefAgrees(SPECIES3, gene, lmap)

        self.assertEqual(r.n_speciations, 1)
        self.assertEqual(r.n_duplications, 0)
        self.assertEqual(r.n_losses, 1)
        self.assertEqual(r.edge_loss_counts[("g0", "a1")], 1)
        self.assertEqual(r.edge_loss_counts[("g0", "c")], 0)
        self.assertEqual(r.edge_losses[("g0", "a1")], ("B",))
        self.assertEqual(r.edge_losses[("g0", "c")], ())

    def test_whole_sister_clade_lost(self):
        # (a,e) on 5 物种：R->X->Y->A，S 计数 3-1=2：
        # 在 X 丢兄弟整支 Z（C,D），在 Y 丢 B
        gene = ("g0", "a", "e")
        lmap = dict(a="A", e="E")
        r = self.assertRefAgrees(SPECIES5, gene, lmap)

        self.assertEqual(r.n_losses, 2)
        self.assertEqual(r.edge_losses[("g0", "a")], ("Z", "B"))
        self.assertEqual(r.edge_losses[("g0", "e")], ())

    def test_duplication_no_subtract_rule(self):
        # g1=(a1,b)->Y(S)；g0=(g1,a2)，LCA(Y,A)=Y，且 g1 映射==Y，故 g0 为 D。
        # 边 g0->a2：Y->A，路径长度 1；D 不减 1，故丢 1 个兄弟 B。
        # 若错按 S 公式（k-1）会得到 0——本用例专门钉住该差异。
        gene = ("g0", ("g1", "a1", "b"), "a2")
        lmap = dict(a1="A", b="B", a2="A")
        r = self.assertRefAgrees(SPECIES5, gene, lmap)

        self.assertEqual(r.n_duplications, 1)
        self.assertEqual(r.n_speciations, 1)
        self.assertEqual(r.n_losses, 1)
        self.assertEqual(r.edge_loss_counts[("g0", "a2")], 1)
        self.assertEqual(r.edge_losses[("g0", "a2")], ("B",))
        self.assertEqual(r.edge_loss_counts[("g0", "g1")], 0)

    def test_multiple_duplications(self):
        # (((a1,a2),(a3,a4)),e)：g1、g2、h 都映到 A，全部为 D，共 3 个 D；
        # g0=LCA(A,E)=R(S)，边 g0->h 走 R->X->Y->A，丢整支 Z 与 B。
        gene = ("g0", ("h", ("g1", "a1", "a2"), ("g2", "a3", "a4")), "e")
        lmap = dict(a1="A", a2="A", a3="A", a4="A", e="E")
        r = self.assertRefAgrees(SPECIES5, gene, lmap)

        self.assertEqual(r.n_duplications, 3)
        self.assertEqual(r.n_speciations, 1)
        self.assertEqual(r.n_losses, 2)
        self.assertEqual(r.edge_losses[("g0", "h")], ("Z", "B"))
        self.assertEqual(
            sorted(e.kind for e in r.events),
            [EVENT_DUPLICATION] * 3 + [EVENT_SPECIATION],
        )

    def test_combined_case_hand_values(self):
        # g0=(h,e)，h=(g1,g2)，g1=(a1,a2)，g2=(a3,c)
        # 手算：
        #   g1=LCA(A,A)=A (D)
        #   g2=LCA(A,C)=X (S)：边 g2->a3 丢 B，边 g2->c 丢 D
        #   h =LCA(A,X)=X，g2 映射==X -> D：边 h->g1 走 X->Y->A，
        #     D 不减 1 => 在 X 丢兄弟整支 Z、在 Y 丢 B，共 2
        #   g0=LCA(X,E)=R (S)，无丢失
        gene = (
            "g0",
            ("h", ("g1", "a1", "a2"), ("g2", "a3", "c")),
            "e",
        )
        lmap = dict(a1="A", a2="A", a3="A", c="C", e="E")
        r = self.assertRefAgrees(SPECIES5, gene, lmap)

        self.assertEqual(
            r.mapping,
            dict(
                a1="A", a2="A", a3="A", c="C", e="E",
                g1="A", g2="X", h="X", g0="R",
            ),
        )
        kinds = {e.gene: e.kind for e in r.events}
        self.assertEqual(
            kinds, dict(g1=EVENT_DUPLICATION, g2=EVENT_SPECIATION,
                        h=EVENT_DUPLICATION, g0=EVENT_SPECIATION)
        )
        self.assertEqual(r.edge_losses[("g2", "a3")], ("B",))
        self.assertEqual(r.edge_losses[("g2", "c")], ("D",))
        self.assertEqual(r.edge_losses[("h", "g1")], ("Z", "B"))
        self.assertEqual(r.n_duplications, 2)
        self.assertEqual(r.n_speciations, 2)
        self.assertEqual(r.n_losses, 4)
        # 丢失支系合计的具体标签
        all_lost = [x for tup in r.edge_losses.values() for x in tup]
        self.assertEqual(sorted(all_lost), ["B", "B", "D", "Z"])


class TestPermutationInvariance(_SpecCase):
    def test_child_swap_invariance(self):
        # 每个内部节点左右互换（叶子在树中的位置随之镜像），
        # 归因结果（映射/事件/逐边计数）必须不变。
        lmap = dict(a="A", b="B", c="C", d="D")
        g1 = ("g0", ("g1", "a", "b"), ("g2", "c", "d"))
        g2 = ("g0", ("g2", "d", "c"), ("g1", "b", "a"))
        r1 = self.assertRefAgrees(SPECIES5, g1, lmap)
        r2 = self.assertRefAgrees(SPECIES5, g2, lmap)
        self.assertEqual(r1.mapping, r2.mapping)
        self.assertEqual(r1.edge_losses, r2.edge_losses)
        self.assertEqual(r1.edge_loss_counts, r2.edge_loss_counts)
        self.assertEqual(
            sorted((e.gene, e.species, e.kind) for e in r1.events),
            sorted((e.gene, e.species, e.kind) for e in r2.events),
        )
        self.assertEqual(
            (r1.n_speciations, r1.n_duplications, r1.n_losses),
            (r2.n_speciations, r2.n_duplications, r2.n_losses),
        )

    def test_leaf_assignment_permutation_within_species(self):
        # 同一物种的多个基因互换名字（连同叶映射一起置换），事件统计不变。
        gene = ("g0", ("g1", "a1", "a2"), ("g2", "a3", "c"))
        base = dict(a1="A", a2="A", a3="A", c="C")
        perm = dict(a1="a3", a2="a1", a3="a2", c="c")  # 新名 -> 旧名
        # 构造置换后的基因树与映射
        renamed = ("g0", ("g1", "x1", "x2"), ("g2", "x3", "xc"))
        lmap2 = {"x1": base[perm["a1"]], "x2": base[perm["a2"]],
                 "x3": base[perm["a3"]], "xc": "C"}
        r1 = self.assertRefAgrees(SPECIES5, gene, base)
        r2 = self.assertRefAgrees(SPECIES5, renamed, lmap2)
        self.assertEqual(
            (r1.n_speciations, r1.n_duplications, r1.n_losses),
            (r2.n_speciations, r2.n_duplications, r2.n_losses),
        )


class TestRandomizedCrossCheck(_SpecCase):
    """确定性随机抽样：被测核心与独立参考必须在大量随机树上一致。

    参考实现不调用被测核心；这里还对每棵基因树做整体镜像，
    核对映射/计数/统计的左右置换不变性。
    """

    @staticmethod
    def _random_tree(rng, prefix, n_leaves):
        labels = [f"{prefix}{i}" for i in range(n_leaves)]
        counter = [0]

        def merge(nodes):
            while len(nodes) > 1:
                i = rng.randrange(len(nodes))
                a = nodes.pop(i)
                j = rng.randrange(len(nodes))
                b = nodes.pop(j)
                counter[0] += 1
                nodes.append((f"{prefix}n{counter[0]}", a, b))
            return nodes[0]

        return merge(labels)

    @staticmethod
    def _mirror(spec):
        if isinstance(spec, str):
            return spec
        label, left, right = spec
        return label, TestRandomizedCrossCheck._mirror(right), TestRandomizedCrossCheck._mirror(left)

    def test_200_random_pairs(self):
        import random

        rng = random.Random(20260923)
        for case in range(200):
            n_sp = rng.randrange(1, 8)
            n_gn = rng.randrange(1, 9)
            sp = self._random_tree(rng, "s", n_sp)
            gn = self._random_tree(rng, "g", n_gn)
            sp_leaves = [f"s{i}" for i in range(n_sp)]
            gl = [f"g{i}" for i in range(n_gn)]
            lmap = {g: rng.choice(sp_leaves) for g in gl}
            r1 = self.assertRefAgrees(sp, gn, lmap)
            r2 = self.assertRefAgrees(sp, self._mirror(gn), lmap)
            self.assertEqual(r1.mapping, r2.mapping)
            self.assertEqual(
                sorted(r1.edge_loss_counts.values()),
                sorted(r2.edge_loss_counts.values()),
            )
            self.assertEqual(
                (r1.n_speciations, r1.n_duplications, r1.n_losses),
                (r2.n_speciations, r2.n_duplications, r2.n_losses),
            )


class TestValidation(unittest.TestCase):
    def _tree(self, root, **kw):
        return Tree(root, name="t", **kw)

    def test_reject_non_leaf_mapping(self):
        # 基因叶映到物种内部节点必须拒绝
        sp = self._tree(node_from_spec(SPECIES3))
        gt = self._tree(node_from_spec(("g0", "a", "c")))
        with self.assertRaisesRegex(MappingError, "物种内部节点"):
            reconcile(sp, gt, dict(a="Y", c="C"))

    def test_reject_missing_and_extra_mapping(self):
        sp = self._tree(node_from_spec(SPECIES3))
        gt = self._tree(node_from_spec(("g0", "a", "c")))
        with self.assertRaisesRegex(MappingError, "缺少映射"):
            reconcile(sp, gt, dict(a="A"))
        with self.assertRaisesRegex(MappingError, "多余键"):
            reconcile(sp, gt, dict(a="A", c="C", z="A"))

    def test_reject_unknown_species_label(self):
        sp = self._tree(node_from_spec(SPECIES3))
        gt = self._tree(node_from_spec(("g0", "a", "c")))
        with self.assertRaisesRegex(MappingError, "不存在的物种标签"):
            reconcile(sp, gt, dict(a="A", c="ZZ"))

    def test_reject_cyclic_tree(self):
        # 手工制造环：n1 -> n2 -> n1
        leaf_l = Node("L")
        leaf_r = Node("R")
        n1 = Node("n1")
        n2 = Node("n2")
        n1.left = n2
        n1.right = leaf_r
        n2.left = n1
        n2.right = leaf_l
        with self.assertRaisesRegex(TreeError, "环"):
            self._tree(n1)

    def test_reject_shared_subtree(self):
        # 同一叶对象被两个父节点共享（多父/DAG）
        shared = Node("s")
        left = Node("u", shared, Node("x"))
        right = Node("v", shared, Node("y"))
        root = Node("r", left, right)
        with self.assertRaisesRegex(TreeError, "多个父节点"):
            self._tree(root)

    def test_reject_unary_node(self):
        root = Node("r", Node("only"), None)
        with self.assertRaisesRegex(TreeError, "不是满二叉"):
            self._tree(root)

    def test_reject_duplicate_labels(self):
        spec = ("r", ("u", "a", "b"), ("v", "a", "c"))
        with self.assertRaisesRegex(TreeError, "标签重复"):
            self._tree(node_from_spec(spec))

    def test_reject_same_object_children(self):
        leaf = Node("a")
        root = Node("r", leaf, leaf)
        with self.assertRaises(TreeError):
            self._tree(root)

    def test_reject_bad_node_type_and_label(self):
        with self.assertRaisesRegex(TreeError, "根必须是 Node"):
            self._tree("not a node")  # type: ignore[arg-type]
        root = Node("", Node("a"), Node("b"))
        with self.assertRaisesRegex(TreeError, "非空 str"):
            self._tree(root)

    def test_max_nodes_bool_and_floats_rejected(self):
        good = node_from_spec(SPECIES3)
        with self.assertRaisesRegex(TreeError, "bool"):
            Tree(good, max_nodes=True)  # type: ignore[arg-type]
        for bad in (1.0, float("nan"), float("inf"), -float("inf")):
            with self.subTest(bad=bad):
                with self.assertRaises(TreeError):
                    Tree(good, max_nodes=bad)  # type: ignore[arg-type]
        # 纯 NaN/Inf 不接受；整数边界也钉一下
        self.assertTrue(math.isnan(float("nan")))
        with self.assertRaisesRegex(TreeError, "200"):
            Tree(node_from_spec(SPECIES3), max_nodes=201)

    def test_node_count_boundary_199_ok_201_rejected(self):
        # 满二叉树节点数恒为奇数。199 <= 200 接受；201 超过上限拒绝。
        def balanced(n_leaves, prefix, counter):
            if n_leaves == 1:
                counter[0] += 1
                return f"{prefix}_leaf{counter[0]}"
            half = n_leaves // 2
            counter[0] += 1
            lab = f"{prefix}_int{counter[0]}"
            return (
                lab,
                balanced(half, prefix, counter),
                balanced(n_leaves - half, prefix, counter),
            )

        spec199 = balanced(100, "s", [0])
        self.assertEqual(Tree(node_from_spec(spec199)).leaves.__len__(), 100)
        spec201 = balanced(101, "s", [0])
        with self.assertRaisesRegex(TreeError, "超过上限 200"):
            Tree(node_from_spec(spec201))
        self.assertEqual(MAX_NODES, 200)

    def test_bad_leaf_map_types(self):
        sp = self._tree(node_from_spec(SPECIES3))
        gt = self._tree(node_from_spec(("g0", "a", "c")))
        with self.assertRaises(MappingError):
            reconcile(sp, gt, [("a", "A"), ("c", "C")])  # type: ignore[arg-type]
        with self.assertRaisesRegex(MappingError, "值必须是物种叶标签"):
            reconcile(sp, gt, dict(a=1, c="C"))  # type: ignore[dict-item]

    def test_wrong_tree_types(self):
        sp = self._tree(node_from_spec(SPECIES3))
        with self.assertRaises(TypeError):
            reconcile(sp, "x", {})  # type: ignore[arg-type]
        with self.assertRaises(TypeError):
            reconcile("x", sp, {})  # type: ignore[arg-type]


class TestMisc(unittest.TestCase):
    def test_newick_and_node_spec(self):
        t = Tree(node_from_spec(SPECIES3))
        self.assertEqual(newick(t), "((A,B)Y,C)R;")
        # 单叶树合法
        one = Tree(node_from_spec("A"))
        self.assertTrue(one.root.is_leaf)

    def test_result_is_frozen_and_label_keyed(self):
        sp = Tree(node_from_spec(SPECIES3))
        gt = Tree(node_from_spec(("g0", "a", "c")))
        r = reconcile(sp, gt, dict(a="A", c="C"))
        self.assertIsInstance(r, Reconciliation)
        with self.assertRaises(Exception):
            r.n_losses = 99  # type: ignore[misc]
        # 所有键/值均为标签（str），不暴露内部 Node 对象
        self.assertTrue(all(isinstance(k, str) for k in r.mapping))
        self.assertTrue(all(isinstance(v, str) for v in r.mapping.values()))

    def test_inputs_not_mutated(self):
        sp_root = node_from_spec(SPECIES3)
        gt_root = node_from_spec(("g0", "a", "c"))
        lmap = dict(a="A", c="C")
        lmap_before = dict(lmap)
        reconcile(Tree(sp_root), Tree(gt_root), lmap)
        self.assertEqual(lmap, lmap_before)
        self.assertEqual(newick(sp_root), "((A,B)Y,C)R;")


if __name__ == "__main__":
    unittest.main()
