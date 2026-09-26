"""不预设复制/丢失权重的 Pareto 前沿入口测试。

期望值来源：
1. 手工可复验小算例（含“多复制少丢失 / 少复制多丢失”双向非支配对）；
2. tests/reference.py 中穷举**全部合法映射**的独立参考（不调用被测核心）；
3. (1,1) 加权结果与旧固定权重入口 reconcile 的逐字段一致性。
"""

from __future__ import annotations

import math
import unittest
from itertools import combinations_with_replacement

from reference import reference_pareto_frontier

from reconciliation import (
    DUPLICATION,
    SPECIATION,
    BudgetExhaustedError,
    InvalidInputError,
    _gene_postorder,
    build_tree,
    reconcile,
    reconcile_pareto,
    reconcile_weighted,
)


SP3 = {"y": ("B", "C"), "r": ("A", "y")}


def binary_shapes(labels):
    labels = tuple(labels)
    if len(labels) == 1:
        yield (labels[0],)
        return
    for k in range(1, len(labels)):
        for left in binary_shapes(labels[:k]):
            for right in binary_shapes(labels[k:]):
                yield (left, right)


def shape_to_tree(shape, prefix, species_labels, *, keep_leaf_labels=False):
    children = {}
    counter = [0]
    leaf_names = {}

    def conv(node):
        if len(node) == 1:
            label = node[0]
            if keep_leaf_labels:
                return label
            name = f"{prefix}_leaf_{label}_{counter[0]}"
            counter[0] += 1
            leaf_names[name] = label
            return name
        name = f"{prefix}_n{counter[0]}"
        counter[0] += 1
        a = conv(node[0])
        b = conv(node[1])
        children[name] = (a, b)
        return name

    root = conv(shape)
    return children, root, leaf_names


def canonical_key(gene_tree, mapping):
    return tuple(mapping[g] for g in _gene_postorder(gene_tree))


class HandCheckedParetoTests(unittest.TestCase):
    def setUp(self):
        self.sp = build_tree(SP3, "r", name="species")
        # 基因树 ((b1,b2),c)：b1,b2 -> B；c -> C
        self.gn = build_tree(
            {"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr", name="gene"
        )
        self.lm = {"b1": "B", "b2": "B", "c": "C"}

    def test_two_way_tradeoff_pair(self):
        """前沿恰为 {(0,1),(1,0)}：少复制多丢失 vs 多复制少丢失。"""
        res = reconcile_pareto(self.sp, self.gn, self.lm)
        self.assertTrue(res.complete)
        self.assertEqual(res.counts(), ((0, 1), (1, 0)))

    def test_each_point_has_full_attribution(self):
        by_count = {(p.duplications, p.losses): p for p in res0_points(self)}
        zero_one = by_count[(0, 1)]
        one_zero = by_count[(1, 0)]
        self.assertEqual(zero_one.result.total_losses, 1)
        self.assertEqual(zero_one.result.duplications, 0)
        self.assertEqual(zero_one.mapping["gr"], "r")
        self.assertEqual(zero_one.mapping["g1"], "y")
        self.assertEqual(one_zero.result.duplications, 1)
        self.assertEqual(one_zero.result.total_losses, 0)
        self.assertEqual(one_zero.mapping["gr"], "y")
        self.assertEqual(one_zero.mapping["g1"], "B")
        for point in (zero_one, one_zero):
            self.assertEqual(
                point.result.duplications,
                sum(e.event_type == DUPLICATION for e in point.result.events),
            )
            self.assertEqual(
                point.result.total_losses,
                sum(e.loss_count for e in point.result.edge_losses),
            )

    def test_no_transfer_events(self):
        """本题物种归属固定：事件类型只有 speciation / duplication。"""
        res = reconcile_pareto(self.sp, self.gn, self.lm)
        for point in res.points:
            self.assertTrue(
                all(
                    e.event_type in (SPECIATION, DUPLICATION)
                    for e in point.result.events
                )
            )

    def test_equal_count_single_witness(self):
        """相同计数只保留规范序最小的一份见证，不出现重复取舍。"""
        res = reconcile_pareto(self.sp, self.gn, self.lm)
        counts = [canonical_key(self.gn, p.mapping) for p in res.points]
        self.assertEqual(len(counts), len(set(counts)))
        # 与独立穷举的规范见证逐键一致。
        front, witnesses = reference_pareto_frontier(
            SP3, "r", {"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr", self.lm
        )
        self.assertEqual(front, ((0, 1), (1, 0)))
        res_by = {(p.duplications, p.losses): canonical_key(self.gn, p.mapping) for p in res.points}
        for count, ref_key in witnesses.items():
            self.assertEqual(res_by[count], ref_key, count)

    def test_single_leaf_gene_tree(self):
        sp = build_tree({"r": ("A", "B")}, "r", name="species")
        gn = build_tree({}, "a", name="gene")
        res = reconcile_pareto(sp, gn, {"a": "A"})
        self.assertTrue(res.complete)
        self.assertEqual(res.counts(), ((0, 0),))
        weighted = reconcile_weighted(
            sp, gn, {"a": "A"}, dup_weight=3, loss_weight=2
        )
        self.assertEqual(weighted.objective, 0)
        self.assertEqual(weighted.tied_counts, ((0, 0),))


def res0_points(case):
    return reconcile_pareto(case.sp, case.gn, case.lm).points


class ExhaustiveSmallTreeTests(unittest.TestCase):
    """小树枚举全部合法映射：DP 前沿与独立穷举参考完全一致。"""

    def test_frontier_matches_bruteforce(self):
        cases = 0
        tradeoff_cases = 0
        for sp_labels in (("A", "B"), ("A", "B", "C")):
            for sp_shape in binary_shapes(sp_labels):
                sp_ch, sp_root, _ = shape_to_tree(
                    sp_shape, "s", sp_labels, keep_leaf_labels=True
                )
                sp = build_tree(sp_ch, sp_root, name="species")
                for ngl in range(2, 5):
                    for multis in combinations_with_replacement(sp_labels, ngl):
                        for g_shape in set(binary_shapes(multis)):
                            g_ch, g_root, leaf_names = shape_to_tree(
                                g_shape, "g", sp_labels
                            )
                            gn = build_tree(g_ch, g_root, name="gene")
                            got = reconcile_pareto(sp, gn, leaf_names)
                            self.assertTrue(got.complete)
                            ref_front, ref_witness = reference_pareto_frontier(
                                sp_ch, sp_root, g_ch, g_root, leaf_names
                            )
                            self.assertEqual(got.counts(), ref_front)
                            got_witness = {
                                (p.duplications, p.losses): canonical_key(gn, p.mapping)
                                for p in got.points
                            }
                            self.assertEqual(got_witness, ref_witness)
                            cases += 1
                            if any(a[0] > b[0] and a[1] < b[1] for a in ref_front for b in ref_front):
                                tradeoff_cases += 1
        self.assertGreater(cases, 50)
        # 必须真的核对到双向非支配（多复制少丢失 & 少复制多丢失）。
        self.assertGreater(tradeoff_cases, 0)

    def test_weighted_1_1_selects_frontier_optimum_and_legacy_is_legal(self):
        """(1,1) 从前沿选 D+L 最小；旧 LCA 复演的计数必为合法计数。

        字面事件语义下抬升映射可能合并复制，故 LCA 映射不必最优，
        但其 (D, L) 一定不被前沿点严格支配之外的非法状态污染——
        即 LCA 计数必须等于某个合法计数（此处用独立穷举参考核对）。
        """
        sp_shapes = list(binary_shapes(("A", "B", "C")))
        checked = 0
        lca_was_optimal = 0
        lca_was_suboptimal = 0
        for sp_shape in sp_shapes:
            sp_ch, sp_root, _ = shape_to_tree(
                sp_shape, "s", ("A", "B", "C"), keep_leaf_labels=True
            )
            sp = build_tree(sp_ch, sp_root, name="species")
            for multis in combinations_with_replacement(("A", "B", "C"), 3):
                for g_shape in set(binary_shapes(multis)):
                    g_ch, g_root, leaf_names = shape_to_tree(g_shape, "g", ("A", "B", "C"))
                    gn = build_tree(g_ch, g_root, name="gene")
                    weighted = reconcile_weighted(
                        sp, gn, leaf_names, dup_weight=1, loss_weight=1
                    )
                    legacy = reconcile(sp, gn, leaf_names)
                    # 加权最优计数必须在 tied_counts 中且目标值最小。
                    self.assertIn(
                        (weighted.best.duplications, weighted.best.total_losses),
                        weighted.tied_counts,
                    )
                    self.assertEqual(
                        weighted.objective,
                        min(
                            p.duplications + p.losses
                            for p in weighted.frontier.points
                        ),
                    )
                    # LCA 复演计数必须是合法计数：不被任何前沿点严格支配
                    # （若被严格支配说明 LCA 次优；若既非前沿点也不被支配则矛盾）。
                    legacy_count = (legacy.duplications, legacy.total_losses)
                    dominated = any(
                        p.duplications <= legacy_count[0]
                        and p.losses <= legacy_count[1]
                        and (p.duplications, p.losses) != legacy_count
                        for p in weighted.frontier.points
                    )
                    on_frontier = legacy_count in weighted.frontier.counts()
                    self.assertTrue(on_frontier or dominated)
                    if dominated:
                        lca_was_suboptimal += 1
                    else:
                        lca_was_optimal += 1
                    checked += 1
        self.assertGreater(checked, 10)
        # 穷举集合里两种情形都必须真实出现，避免测试空洞。
        self.assertGreater(lca_was_optimal, 0)
        self.assertGreater(lca_was_suboptimal, 0)

    def test_weights_select_different_frontier_points(self):
        sp = build_tree(SP3, "r", name="species")
        gn = build_tree({"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr")
        lm = {"b1": "B", "b2": "B", "c": "C"}
        cheap_loss = reconcile_weighted(sp, gn, lm, dup_weight=10, loss_weight=1)
        cheap_dup = reconcile_weighted(sp, gn, lm, dup_weight=1, loss_weight=10)
        self.assertEqual(cheap_loss.tied_counts, ((0, 1),))
        self.assertEqual(cheap_dup.tied_counts, ((1, 0),))
        self.assertEqual(cheap_loss.objective, 1)
        self.assertEqual(cheap_dup.objective, 1)


class BudgetTests(unittest.TestCase):
    def setUp(self):
        self.sp = build_tree(SP3, "r", name="species")
        self.gn = build_tree({"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr")
        self.lm = {"b1": "B", "b2": "B", "c": "C"}

    def test_exhaustion_marked_incomplete(self):
        res = reconcile_pareto(self.sp, self.gn, self.lm, budget=1)
        self.assertFalse(res.complete)
        self.assertEqual(res.budget, 1)
        self.assertGreaterEqual(res.combinations_explored, 2)
        # 明确：预算耗尽时前沿可能不全；点之间仍互不支配。
        for i, p in enumerate(res.points):
            for q in res.points[i + 1 :]:
                self.assertFalse(
                    (p.duplications <= q.duplications and p.losses <= q.losses)
                    or (q.duplications <= p.duplications and q.losses <= p.losses)
                )

    def test_weighted_raises_on_exhaustion_with_partial(self):
        with self.assertRaises(BudgetExhaustedError) as ctx:
            reconcile_weighted(
                self.sp,
                self.gn,
                self.lm,
                dup_weight=1,
                loss_weight=1,
                budget=1,
            )
        self.assertFalse(ctx.exception.partial.complete)
        self.assertIn("\u524d\u6cbf", str(ctx.exception))  # “前沿”

    def test_generous_budget_completes(self):
        res = reconcile_pareto(self.sp, self.gn, self.lm, budget=1_000_000)
        self.assertTrue(res.complete)


class ParameterValidationTests(unittest.TestCase):
    def setUp(self):
        self.sp = build_tree(SP3, "r", name="species")
        self.gn = build_tree({"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr")
        self.lm = {"b1": "B", "b2": "B", "c": "C"}

    def test_bad_budget(self):
        for bad in (0, -1, 1.0, True, False, math.nan, math.inf, "100"):
            with self.assertRaises(InvalidInputError):
                reconcile_pareto(self.sp, self.gn, self.lm, budget=bad)

    def test_bad_weights(self):
        for kwargs in (
            {"dup_weight": 0, "loss_weight": 1},
            {"dup_weight": 1, "loss_weight": -2},
            {"dup_weight": 1.5, "loss_weight": 1},
            {"dup_weight": True, "loss_weight": 1},
            {"dup_weight": 1, "loss_weight": math.nan},
        ):
            with self.assertRaises(InvalidInputError):
                reconcile_weighted(self.sp, self.gn, self.lm, **kwargs)

    def test_weights_required_as_keyword(self):
        with self.assertRaises(TypeError):
            reconcile_weighted(self.sp, self.gn, self.lm, 1, 1)

    def test_leaf_mapping_validation_shared(self):
        with self.assertRaises(InvalidInputError):
            reconcile_pareto("not a tree", self.gn, self.lm)


if __name__ == "__main__":
    unittest.main()
