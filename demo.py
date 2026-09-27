"""固定输入演示：真实计算一次含复制与丢失的复演，并实际触发一个拒绝边界。

运行：python demo.py
数据全部在本文件内确定性给出，无网络、无外部数据、无 sleep、无预录结果。
"""

from __future__ import annotations

from reconciliation import (
    BUDGET_EXHAUSTED,
    OPTIMAL_CONFIRMED,
    DUPLICATION,
    SPECIATION,
    LeafMappingError,
    TreeStructureError,
    build_tree,
    reconcile,
    reconcile_candidates,
)


def show_tree(title, tree):
    print(f"{title}（根 {tree.root}，{len(tree.internal)} 内节点 / {len(tree.leaves)} 叶）：")
    for node in sorted(tree.internal):
        a, b = tree.children[node]
        print(f"    {node} -> ({a}, {b})")


def main():
    print("=" * 72)
    print("部分一：正常复演（含 1 次复制、2 次丢失，事件路径可手工复验）")
    print("=" * 72)

    # 物种树:
    #        r
    #       / \\
    #      x   D
    #     / \\
    #    y   C
    #   / \\
    #  A   B
    species_children = {"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}
    # 基因树（同一物种 A 有两个旁系同源拷贝 a1/a2）:
    #       gr
    #      /  \\
    #    g1    g2
    #   / \\   / \\
    #  a1 b1 a2  c1
    gene_children = {
        "g1": ("a1", "b1"),
        "g2": ("a2", "c1"),
        "gr": ("g1", "g2"),
    }
    leaf_mapping = {"a1": "A", "b1": "B", "a2": "A", "c1": "C"}

    sp = build_tree(species_children, "r", name="species")
    gn = build_tree(gene_children, "gr", name="gene")
    show_tree("物种树", sp)
    show_tree("基因树", gn)
    print("基因叶 -> 物种叶：")
    for g, s in sorted(leaf_mapping.items()):
        print(f"    {g} -> {s}")

    res = reconcile(sp, gn, leaf_mapping)

    print("\n[1] LCA 映射（内节点 = 两个子映射的 LCA）：")
    for g, s in res.mapping.items():
        kind = "叶" if g in gn.leaves else "内节点"
        print(f"    {g}（{kind}） -> {s}")

    print("\n[2] 事件表（任一子映射 == 父映射 => duplication）：")
    for ev in res.events:
        c0, c1 = ev.child_genes
        m0, m1 = ev.child_maps
        flag = "  <== 复制" if ev.event_type == DUPLICATION else ""
        print(
            f"    {ev.gene_node}: LCA({m0},{m1}) = {ev.mapped_species}；"
            f"子 {c0}->{m0}, {c1}->{m1}；{ev.event_type}{flag}"
        )

    print("\n[3] 逐基因父子边的丢失计数：")
    print("    （d=物种树路径边数；speciation 丢 d-1，duplication 丢 d）")
    for edge in res.edge_losses:
        lost = (
            "缺失支系: " + ", ".join(
                f"{l.species_node}(叶 {sorted(l.descendant_leaves)})"
                for l in edge.lost_lineages
            )
            if edge.lost_lineages
            else "无缺失"
        )
        print(
            f"    {edge.parent_gene}->{edge.child_gene}: "
            f"{edge.parent_species}->{edge.child_species}, d={edge.path_length}, "
            f"丢失={edge.loss_count}；{lost}"
        )

    print("\n[4] 累计：")
    print(f"    speciation={res.speciations}, duplication={res.duplications}, "
          f"总丢失={res.total_losses}")
    print(f"    按缺失支系统计: {dict(res.losses_by_lineage)}")

    print(
        """
[5] 手工复验路径：
    g1=LCA(A,B)=y -> speciation；g2=LCA(A,C)=x -> speciation。
    gr=LCA(y,x)=x；子 g2 也映到 x（子映射==父映射）=> gr 为 duplication。
    边 gr->g1（复制，不减 1）：x->y 路径 d=1，丢 1 = C 整支缺失。
    边 gr->g2（复制）：x->x 路径 d=0，丢 0。
    边 g2->a2（分化）：x->y->A 路径 d=2，丢 d-1=1；
        起点 x 的兄弟 C 被另一子支 c1 占据（不计），y 处兄弟 B 缺失 => 丢 B。
    其余边 d=1（分化）或 d=0，均丢 0。合计 B、C 各 1，总丢失 2。
"""
    )

    print("=" * 72)
    print("部分二：实际触发的拒绝边界（非叶映射 / 循环树，均不静默纠正）")
    print("=" * 72)

    # 边界 1：基因叶映到物种“内节点” x —— 必须拒绝，且可定位到基因叶。
    bad_mapping = dict(leaf_mapping)
    bad_mapping["c1"] = "x"  # x 是物种内节点，不是物种叶
    try:
        reconcile(sp, gn, bad_mapping)
    except LeafMappingError as exc:
        print(f"    非叶映射被拒绝：{type(exc).__name__}: {exc}")
        print(f"    （定位字段 gene_leaf={exc.gene_leaf!r}，调用方可知错在哪个基因叶）")

    # 边界 2：循环树 —— b 的子是 a，a 的子是 b。
    try:
        build_tree({"a": ("b", "c"), "b": ("d", "a")}, "a", name="gene")
    except TreeStructureError as exc:
        print(f"    循环树被拒绝：{type(exc).__name__}: {exc}")
        print(f"    （定位字段 tree={exc.tree!r}, node={exc.node!r}）")

    print("=" * 72)
    print("部分三：候选叶选择（最小化复制+丢失总代价，并列解稳定返回）")
    print("=" * 72)

    # 物种树: r -> (x, C), x -> (A, B)；基因树 gr -> (u, v)。
    # u 归属未定：候选 {A, B}；v 固定为 C。
    sp3 = build_tree({"x": ("A", "B"), "r": ("x", "C")}, "r", name="species")
    gn3 = build_tree({"gr": ("u", "v")}, "gr", name="gene")
    candidates = {"u": ["B", "A"], "v": "C"}  # 故意打乱候选顺序
    print("    候选（输入顺序不影响结果排序）：", candidates)

    search = reconcile_candidates(
        sp3, gn3, candidates, duplication_weight=1, loss_weight=1
    )
    print(
        f"    候选空间 {search.total_combinations} 个完整映射，"
        f"已评估 {search.evaluated} 个，状态={search.status}（"
        + ("已确认全局最优" if search.status == OPTIMAL_CONFIRMED else "尚未确认最优")
        + ")"
    )
    print(f"    最低总代价 = {search.best_cost}"
          f"（复制权重 {search.duplication_weight}，丢失权重 {search.loss_weight}）")
    for sol in search.solutions:
        chosen = ", ".join(f"{g}->{s}" for g, s in sorted(sol.leaf_mapping.items()))
        print(
            f"      解 [{chosen}]：复制 {sol.result.duplications}、"
            f"丢失 {sol.result.total_losses}，"
            f"代价拆分 复制 {sol.duplication_cost} + 丢失 {sol.loss_cost} "
            f"= {sol.total_cost}"
        )
    print(
        """
    手工复验：u=A 时 gr=LCA(A,C)=r 为分化，边 gr->u（r->x->A, d=2）丢 1 = B 支系；
    u=B 时同理丢 1 = A 支系。两解代价均为 1，按选定物种元组 (A,C)<(B,C)
    稳定排列；事件判定与丢失计数完全沿用部分一的同一套 reconcile 规则。
"""
    )

    # 额度不足的实际触发：4 个完整映射只允许评估 2 个。
    sp2 = build_tree({"r": ("A", "B")}, "r", name="species")
    gn2 = build_tree({"gr": ("a1", "a2")}, "gr", name="gene")
    limited = reconcile_candidates(
        sp2, gn2, {"a1": ["A", "B"], "a2": ["A", "B"]}, search_budget=2
    )
    print(
        f"    额度演示：候选空间 {limited.total_combinations}、仅评估 "
        f"{limited.evaluated}，状态={limited.status}；"
        + ("已确认最优。" if limited.status == OPTIMAL_CONFIRMED
           else "调用方应据此知晓返回解尚未确认全局最优。")
    )

    print("\n演示结束：以上结果均由本库当场计算，未使用预录输出。")


if __name__ == "__main__":
    main()
