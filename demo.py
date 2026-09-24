"""固定输入演示：真实计算一次含复制与丢失的复演、候选横向转移合法性检查，
并实际触发拒绝边界。

运行：python demo.py
数据全部在本文件内确定性给出，无网络、无外部数据、无 sleep、无预录结果。
"""

from __future__ import annotations

from reconciliation import (
    DUPLICATION,
    SPECIATION,
    InvalidInputError,
    LeafMappingError,
    TreeStructureError,
    build_tree,
    check_candidate_transfer,
    reconcile,
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
    print("部分二：候选横向转移合法性检查（在上面的复演结果上，叶映射固定）")
    print("=" * 72)
    print(
        "    物种树分支以下端节点标识，分支 (父 p, 子 c) 的时间区间为 [time(p), time(c)]。"
        "\n    缺省时间：根 r=0，每条父子边 +1"
        "（x=1, y=2, C=2, D=1, A=B=3）。"
        "\n    合法条件：两条不同分支、互不为祖先、闭时间区间有公共时刻。"
    )

    def show_check(title, out):
        print(f"\n  {title}")
        print(f"    供体 {out.donor_branch} 区间 {out.donor_interval}；"
              f"受体 {out.recipient_branch} 区间 {out.recipient_interval}")
        if out.feasible:
            lo, hi = out.feasible_interval
            note = "（端点相交：零长度窗口，转移只能发生在该时刻）" if lo == hi else ""
            print(f"    => 可行，可转移时间区间 {out.feasible_interval} {note}")
        else:
            print(f"    => 拒绝，原因 {out.reason}")

    # [正常] 严格重叠窗口：自定义时间下 A 支 (7,8)、C 支 (6,9)。
    wide_times = {"r": 0, "x": 6, "y": 7, "A": 8, "B": 9, "C": 9, "D": 7}
    show_check(
        "[正常] 堂表分支 A -> C（自定义时间，区间严格重叠）",
        check_candidate_transfer(res, "A", "C", node_times=wide_times),
    )

    # [边界] 端点相交：缺省时间下 B 支 (2,3) 与 C 支 (1,2) 仅共享 t=2。
    show_check(
        "[边界] 堂表分支 B -> C（缺省时间，仅在端点 t=2 相交）",
        check_candidate_transfer(res, "B", "C"),
    )

    # [失败 1] 祖先/后裔分支：x 支 (0,1) 与 y 支 (1,2) 在 t=1 端点相接，
    #          但祖先判定优先，仍拒绝。
    show_check(
        "[失败] x -> y（祖先/后裔分支，拒绝优先于时间判定）",
        check_candidate_transfer(res, "x", "y"),
    )

    # [失败 2] 互不为祖先但时间不相交：B 支 (2,3) 与 D 支 (0,1)。
    show_check(
        "[失败] B -> D（互不为祖先，但时间区间不相交）",
        check_candidate_transfer(res, "B", "D"),
    )

    # [失败 3] 同一分支。
    show_check(
        "[失败] A -> A（同一分支，没有供体/受体之分）",
        check_candidate_transfer(res, "A", "A"),
    )

    # [实际错误] 根没有入边，不能作为分支 —— 输入非法，明确报错。
    try:
        check_candidate_transfer(res, "r", "A")
    except InvalidInputError as exc:
        print(f"\n  [输入错误] 根节点充当分支被拒绝：{type(exc).__name__}: {exc}")

    # [实际错误] 节点时间不满足沿父子边严格递增（y=2 不大于 x=2）。
    bad_times = {"r": 0, "x": 2, "y": 2, "A": 3, "B": 3, "C": 3, "D": 3}
    try:
        check_candidate_transfer(res, "A", "C", node_times=bad_times)
    except InvalidInputError as exc:
        print(f"  [输入错误] 非严格递增时间被拒绝：{type(exc).__name__}: {exc}")

    print(
        """
  手工复验：
    缺省时间下 B 支为 [t(y),t(B)]=[2,3]，C 支为 [t(x),t(C)]=[1,2]，
    交集为 [max(2,1), min(3,2)]=[2,2] —— 端点相交可行。
    B 支 [2,3] 与 D 支 [t(r),t(D)]=[0,1] 交集为空（2>1）=> 时间不相交。
    x、y 在同一条根-叶脉络上（y 是 x 的子节点）=> 祖先分支拒绝。
  以上检查只读复演结果，叶映射与事件结果均未改变。"""
    )

    print("\n" + "=" * 72)
    print("部分三：原有复演的拒绝边界（非叶映射 / 循环树，均不静默纠正）")
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

    print("\n演示结束：以上结果均由本库当场计算，未使用预录输出。")


if __name__ == "__main__":
    main()
