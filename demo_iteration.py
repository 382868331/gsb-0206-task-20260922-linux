"""离线演示：不预设复制/丢失权重的 Pareto 前沿协调。

运行：python demo_iteration.py
全程真实计算，无网络、无账号/密钥、无 sleep、无预录结果；
在本机（Python 3.14.7）实测 < 1 秒（要求不超过 90 秒）。
"""

from __future__ import annotations

from reconciliation import (
    BudgetExhaustedError,
    build_tree,
    reconcile,
    reconcile_pareto,
    reconcile_weighted,
)


def print_attribution(point):
    print(f"    复制={point.duplications}  丢失={point.losses}")
    print("    映射：", dict(sorted(point.mapping.items())))
    for edge in point.result.edge_losses:
        lost = [lineage.species_node for lineage in edge.lost_lineages]
        print(
            f"      边 {edge.parent_gene}->{edge.child_gene}: "
            f"{edge.parent_species}->{edge.child_species}, "
            f"路径长={edge.path_length}, 丢失={edge.loss_count}, 缺失支系={lost}"
        )


def main():
    print("=" * 72)
    print("部分一：Pareto 前沿——多复制少丢失 vs 少复制多丢失")
    print("=" * 72)
    # 物种树:  r = (A, y), y = (B, C)
    sp = build_tree({"y": ("B", "C"), "r": ("A", "y")}, "r", name="species")
    # 基因树:  gr = (g1, c), g1 = (b1, b2)；B 有两个旁系拷贝
    gn = build_tree(
        {"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr", name="gene"
    )
    leaf_mapping = {"b1": "B", "b2": "B", "c": "C"}
    pareto = reconcile_pareto(sp, gn, leaf_mapping)
    print("不被支配的计数取舍：", pareto.counts())
    print(f"前沿完整：{pareto.complete}，DP 展开候选组合 {pareto.combinations_explored} 个")
    for point in pareto.points:
        print_attribution(point)

    print("=" * 72)
    print("部分二：调用方给正整数权重，从前沿选优（不重新搜索）")
    print("=" * 72)
    cheap_dup = reconcile_weighted(sp, gn, leaf_mapping, dup_weight=1, loss_weight=10)
    cheap_loss = reconcile_weighted(sp, gn, leaf_mapping, dup_weight=10, loss_weight=1)
    print(
        "权重 (复制=1, 丢失=10) -> 最优并列",
        cheap_dup.tied_counts,
        "目标值",
        cheap_dup.objective,
    )
    print(
        "权重 (复制=10, 丢失=1) -> 最优并列",
        cheap_loss.tied_counts,
        "目标值",
        cheap_loss.objective,
    )

    print("=" * 72)
    print("部分三：旧固定权重入口（LCA 复演）仍是合法见证")
    print("=" * 72)
    legacy = reconcile(sp, gn, leaf_mapping)
    print(f"LCA 复演：复制={legacy.duplications}，丢失={legacy.total_losses}")
    print("LCA 映射：", dict(sorted(legacy.mapping.items())))
    print("该计数在前沿上：", (legacy.duplications, legacy.total_losses) in pareto.counts())

    print("=" * 72)
    print("部分四：搜索预算耗尽（与“无解”严格区分，前沿明确标注可能不全）")
    print("=" * 72)
    partial = reconcile_pareto(sp, gn, leaf_mapping, budget=1)
    print(
        f"budget=1 -> complete={partial.complete}, "
        f"已展开 {partial.combinations_explored} 个候选, 已知点={partial.counts()}"
    )
    try:
        reconcile_weighted(
            sp, gn, leaf_mapping, dup_weight=1, loss_weight=1, budget=1
        )
    except BudgetExhaustedError as exc:
        print("加权入口抛出 BudgetExhaustedError（不假装已证明最优）：")
        print("   ", str(exc))

    print("=" * 72)
    print("部分五：小树穷举自检（枚举全部合法映射核对前沿）")
    print("=" * 72)
    # 直接复用测试中的独立穷举参考，在 3 叶物种树上实时核对一个算例。
    import os
    import sys

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "tests"))
    from reference import reference_pareto_frontier

    ref_front, _ = reference_pareto_frontier(
        {"y": ("B", "C"), "r": ("A", "y")},
        "r",
        {"g1": ("b1", "b2"), "gr": ("g1", "c")},
        "gr",
        leaf_mapping,
    )
    print("独立穷举参考前沿：", ref_front)
    assert ref_front == pareto.counts()
    print("一致：DP 前沿 == 全部合法映射穷举前沿。")


if __name__ == "__main__":
    main()
