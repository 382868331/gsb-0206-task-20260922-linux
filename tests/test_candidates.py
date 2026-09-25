"""候选物种叶集合联合归属的测试。

期望值来源：
1. 手工构造的小算例（含“各自看似最优、合起来更差”、零权重平局、
   一物种多基因）硬编码；
2. 独立暴力枚举（itertools.product 全组合 + 公开 reconcile 重算成本）
   作为参照，与被测 DP 结果逐字段交叉核对——参照只调用公开结果接口，
   不复用被测搜索代码；
3. 输入校验与预算耗尽（incomplete）边界。
"""

from __future__ import annotations

import unittest
from itertools import product

from reconciliation import (
    DUPLICATION,
    SPECIATION,
    CandidateMappingError,
    InvalidInputError,
    LeafMappingError,
    build_tree,
    reconcile,
    reconcile_with_candidates,
)

# 物种树：ab:(A,B), cd:(C,D), root:(ab,cd)
SP4 = {"ab": ("A", "B"), "cd": ("C", "D"), "root": ("ab", "cd")}
# 物种树：y:(A,B), x:(y,C), r:(x,D)
SP = {"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}


def make(children, root, name):
    return build_tree(children, root, name=name)


def brute_force(sp, gn, candidates, dup_w, loss_w):
    """独立参照：枚举全部候选组合，用公开 reconcile 计算成本，取
    (成本, 基因叶ID序的物种ID序列) 字典序最小。"""
    leaves = sorted(gn.leaves)
    best = None
    for combo in product(*(candidates[g] for g in leaves)):
        mapping = dict(zip(leaves, combo))
        rec = reconcile(sp, gn, mapping)
        cost = dup_w * rec.duplications + loss_w * rec.total_losses
        key = (cost, combo)
        if best is None or key < best[0]:
            best = (key, mapping, rec)
    return best


class JointChoiceTests(unittest.TestCase):
    """两个叶各自看似最优、合起来事件更多的经典陷阱。"""

    def test_joint_beats_greedy_per_leaf(self):
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        # 每叶候选 {A,B}：逐叶独立选“最近”无法判断，若都选 A（或都选 B），
        # 根映射到 A（或 B）=> duplication（成本含 dup_weight）；
        # 联合最优是 x->A, y->B（或镜像）：根映射 ab => speciation，零丢失。
        res = reconcile_with_candidates(
            sp, gn, {"x": {"A", "B"}, "y": {"A", "B"}}
        )
        self.assertTrue(res.complete)
        self.assertEqual(res.optimal_cost, 0)
        self.assertEqual(res.reconciliation.duplications, 0)
        self.assertEqual(res.reconciliation.total_losses, 0)
        self.assertEqual(res.reconciliation.events[0].event_type, SPECIATION)
        # 平局消解：序列 ("A","B") < ("B","A")
        self.assertEqual(res.assignments, {"x": "A", "y": "B"})

    def test_greedy_same_choice_is_worse(self):
        """验证陷阱真实存在：两叶同选 A 时确实产生 duplication。"""
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        rec = reconcile(sp, gn, {"x": "A", "y": "A"})
        self.assertEqual(rec.duplications, 1)
        self.assertEqual(rec.events[0].event_type, DUPLICATION)

    def test_zero_weights_tiebreak(self):
        """零权重下所有组合成本同为 0，按基因叶 ID 序的物种 ID 序列取最小。"""
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        res = reconcile_with_candidates(
            sp,
            gn,
            {"x": {"B", "A", "D"}, "y": {"C", "A"}},
            dup_weight=0,
            loss_weight=0,
        )
        self.assertTrue(res.complete)
        self.assertEqual(res.optimal_cost, 0)
        self.assertEqual(res.assignments, {"x": "A", "y": "A"})

    def test_one_species_multiple_genes(self):
        """一物种多基因：两个基因叶都可合法地落到同一物种叶。"""
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        res = reconcile_with_candidates(
            sp, gn, {"x": {"A"}, "y": {"A", "B"}}, dup_weight=1, loss_weight=1
        )
        # x 固定 A；y 选 A 则根 dup（成本 1），y 选 B 则 speciation 零成本。
        self.assertEqual(res.assignments, {"x": "A", "y": "B"})
        self.assertEqual(res.optimal_cost, 0)
        # 强制同种时结果也可验证。
        res2 = reconcile_with_candidates(sp, gn, {"x": {"A"}, "y": {"A"}})
        self.assertEqual(res2.assignments, {"x": "A", "y": "A"})
        self.assertEqual(res2.reconciliation.duplications, 1)
        self.assertEqual(res2.optimal_cost, 1)

    def test_weights_change_optimum(self):
        """同一候选集，权重不同则最优归属不同。"""
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        # x 固定 A；y 候选 {A, D}：
        #   y=A -> duplication（0 丢失），成本 = dup_weight；
        #   y=D -> speciation，两条边 d=2 各丢 1，成本 = 2 * loss_weight。
        res = reconcile_with_candidates(
            sp, gn, {"x": "A", "y": {"A", "D"}}, dup_weight=1, loss_weight=0
        )
        self.assertEqual(res.assignments["y"], "D")  # 0 < 1
        self.assertEqual(res.optimal_cost, 0)
        res2 = reconcile_with_candidates(
            sp, gn, {"x": "A", "y": {"A", "D"}}, dup_weight=3, loss_weight=1
        )
        self.assertEqual(res2.assignments["y"], "D")  # 2 < 3
        self.assertEqual(res2.optimal_cost, 2)
        res3 = reconcile_with_candidates(
            sp, gn, {"x": "A", "y": {"A", "D"}}, dup_weight=1, loss_weight=1
        )
        self.assertEqual(res3.assignments["y"], "A")  # 1 < 2
        self.assertEqual(res3.optimal_cost, 1)


class EquivalenceTests(unittest.TestCase):
    """单归属输入与旧接口结果完全一致。"""

    def test_single_assignment_matches_reconcile(self):
        sp = make(SP, "r", "species")
        gn = make(
            {"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")},
            "gr",
            "gene",
        )
        mapping = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}
        old = reconcile(sp, gn, mapping)
        new = reconcile_with_candidates(sp, gn, mapping)
        self.assertTrue(new.complete)
        self.assertEqual(new.assignments, mapping)
        self.assertEqual(new.reconciliation, old)
        self.assertEqual(
            new.optimal_cost, old.duplications + old.total_losses
        )

    def test_candidate_wrapping_matches_scalar(self):
        sp = make(SP, "r", "species")
        gn = make({"g": ("a", "b")}, "g", "gene")
        r1 = reconcile_with_candidates(sp, gn, {"a": "A", "b": "B"})
        r2 = reconcile_with_candidates(sp, gn, {"a": {"A"}, "b": {"B"}})
        self.assertEqual(r1.assignments, r2.assignments)
        self.assertEqual(r1.reconciliation, r2.reconciliation)

    def test_duplicate_candidates_deduped(self):
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        r1 = reconcile_with_candidates(sp, gn, {"x": ["A", "A", "B"], "y": ("B", "A")})
        r2 = reconcile_with_candidates(sp, gn, {"x": {"A", "B"}, "y": {"A", "B"}})
        self.assertEqual(r1.assignments, r2.assignments)
        self.assertEqual(r1.optimal_cost, r2.optimal_cost)


class CrossCheckTests(unittest.TestCase):
    """与独立暴力枚举在小树全组合上交叉核对最优值与平局映射。"""

    def test_bruteforce_crosscheck(self):
        sp = make(SP4, "root", "species")
        gene_trees = [
            ({"g": ("x", "y")}, "g"),
            ({"g1": ("x", "y"), "g": ("g1", "z")}, "g"),
        ]
        cand_pool = ["A", "B", "C", "D"]
        for children, root in gene_trees:
            gn = make(children, root, "gene")
            leaves = sorted(gn.leaves)
            # 确定性生成一批候选组合（含单子集与多子集）。
            cases = []
            for mask in range(1, 4 ** len(leaves)):
                cands = {}
                for i, g in enumerate(leaves):
                    pick = (mask // (3 ** i)) % 3
                    if pick == 0:
                        cands[g] = cand_pool[:2]
                    elif pick == 1:
                        cands[g] = cand_pool[1:3]
                    else:
                        cands[g] = cand_pool[2:]
                cases.append(cands)
            for cands in cases:
                for dup_w, loss_w in ((1, 1), (0, 1), (1, 0), (2, 3), (0, 0)):
                    res = reconcile_with_candidates(
                        sp, gn, cands, dup_weight=dup_w, loss_weight=loss_w
                    )
                    (exp_key, exp_map, _), = [brute_force(
                        sp, gn,
                        {g: sorted(set(cands[g])) for g in leaves},
                        dup_w, loss_w,
                    )]
                    self.assertTrue(res.complete)
                    self.assertEqual(
                        res.optimal_cost, exp_key[0],
                        msg=f"{cands} w=({dup_w},{loss_w})",
                    )
                    self.assertEqual(
                        res.assignments, exp_map,
                        msg=f"{cands} w=({dup_w},{loss_w})",
                    )
                    # 账本可验证：返回成本与最优值一致。
                    self.assertEqual(res.best_cost, res.optimal_cost)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.sp = make(SP4, "root", "species")
        self.gn = make({"g": ("x", "y")}, "g", "gene")

    def test_empty_candidate_set_rejected(self):
        with self.assertRaises(CandidateMappingError):
            reconcile_with_candidates(self.sp, self.gn, {"x": [], "y": "A"})

    def test_nonleaf_candidate_rejected(self):
        with self.assertRaises(CandidateMappingError) as cm:
            reconcile_with_candidates(self.sp, self.gn, {"x": {"ab"}, "y": "A"})
        self.assertEqual(cm.exception.gene_leaf, "x")

    def test_unknown_candidate_rejected(self):
        with self.assertRaises(CandidateMappingError):
            reconcile_with_candidates(self.sp, self.gn, {"x": {"ZZ"}, "y": "A"})

    def test_missing_leaf_rejected(self):
        with self.assertRaises(CandidateMappingError):
            reconcile_with_candidates(self.sp, self.gn, {"x": "A"})

    def test_unknown_gene_leaf_rejected(self):
        with self.assertRaises(CandidateMappingError):
            reconcile_with_candidates(
                self.sp, self.gn, {"x": "A", "y": "B", "zz": "A"}
            )

    def test_candidate_error_is_leaf_mapping_error(self):
        with self.assertRaises(LeafMappingError):
            reconcile_with_candidates(self.sp, self.gn, {"x": set(), "y": "A"})

    def test_bad_weight_types(self):
        for bad in (True, 1.5, float("nan"), float("inf"), -1, "1"):
            with self.assertRaises(InvalidInputError, msg=f"dup_weight={bad!r}"):
                reconcile_with_candidates(
                    self.sp, self.gn, {"x": "A", "y": "B"}, dup_weight=bad
                )
            with self.assertRaises(InvalidInputError, msg=f"loss_weight={bad!r}"):
                reconcile_with_candidates(
                    self.sp, self.gn, {"x": "A", "y": "B"}, loss_weight=bad
                )

    def test_bad_candidate_container(self):
        for bad in (123, None, b"AB"):
            with self.assertRaises(CandidateMappingError, msg=f"{bad!r}"):
                reconcile_with_candidates(self.sp, self.gn, {"x": bad, "y": "A"})

    def test_non_mapping_input(self):
        with self.assertRaises(InvalidInputError):
            reconcile_with_candidates(self.sp, self.gn, [("x", "A")])


class IncompleteTests(unittest.TestCase):
    """预算耗尽：明确 incomplete，保留最优已知映射与可验证账本。"""

    def test_budget_exhausted_returns_incomplete(self):
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        cands = {"x": {"A", "B", "C", "D"}, "y": {"A", "B", "C", "D"}}
        res = reconcile_with_candidates(sp, gn, cands, max_relaxations=1)
        self.assertFalse(res.complete)
        self.assertIsNone(res.optimal_cost)
        # 最优已知映射仍可完整复演验证。
        rec = reconcile(sp, gn, res.assignments)
        self.assertEqual(rec, res.reconciliation)
        self.assertEqual(
            res.best_cost, rec.duplications + rec.total_losses
        )

    def test_budget_result_is_deterministic(self):
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        cands = {"x": {"B", "A"}, "y": {"D", "C"}}
        r1 = reconcile_with_candidates(sp, gn, cands, max_relaxations=1)
        r2 = reconcile_with_candidates(sp, gn, cands, max_relaxations=1)
        self.assertEqual(r1.assignments, r2.assignments)

    def test_bad_budget_rejected(self):
        sp = make(SP4, "root", "species")
        gn = make({"g": ("x", "y")}, "g", "gene")
        for bad in (0, -3, True, 2.5):
            with self.assertRaises(InvalidInputError, msg=f"{bad!r}"):
                reconcile_with_candidates(
                    sp, gn, {"x": "A", "y": "B"}, max_relaxations=bad
                )


if __name__ == "__main__":
    unittest.main()
