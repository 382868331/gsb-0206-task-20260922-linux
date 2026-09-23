#!/usr/bin/env python3
"""固定输入演示：LCA 基因树/物种树事件归因。

运行：python demo.py（约 8 秒内完成，无网络、无第三方依赖）。

前半部分打印一个同时含物种分化、复制与多支丢失的正常归因结果，
并逐条给出可手工复验的物种树路径；后半部分实际触发一个拒绝边界
（基因叶映射到物种内部节点 -> MappingError）。
"""

from reconciliation import (
    EVENT_DUPLICATION,
    EVENT_SPECIATION,
    MappingError,
    Tree,
    newick,
    node_from_spec,
    reconcile,
)

# 物种树（9 节点，<=200）：
#        R
#       / \
#       X  E
#      / \
#     Y   Z
#    / \ / \
#    A B C D
SPECIES = ("R", ("X", ("Y", "A", "B"), ("Z", "C", "D")), "E")

# 基因树（9 节点）：
#        g0
#       /  \
#       h   e
#      / \
#     g1  g2
#    / \  / \
#   a1 a2 a3 c
GENE = (
    "g0",
    ("h", ("g1", "a1", "a2"), ("g2", "a3", "c")),
    "e",
)
LEAF_MAP = {"a1": "A", "a2": "A", "a3": "A", "c": "C", "e": "E"}


def main() -> None:
    species_tree = Tree(node_from_spec(SPECIES), name="species")
    gene_tree = Tree(node_from_spec(GENE), name="gene")

    print("=" * 68)
    print("正常结果：LCA 事件归因")
    print("=" * 68)
    print("物种树:", newick(species_tree))
    print("基因树:", newick(gene_tree))
    print("基因叶 -> 物种叶:", LEAF_MAP)

    result = reconcile(species_tree, gene_tree, LEAF_MAP)

    print("\n[1] 全节点 LCA 映射（后序推导：叶子直映，内节点=子映射的 LCA）")
    for label, mapped in result.mapping.items():
        print(f"    {label:>3} -> {mapped}")

    print("\n[2] 事件表（任一子映射 == 父映射 即为 D，否则 S）")
    for ev in result.events:
        ch = ", ".join(f"{g}->{s}" for g, s in ev.children)
        print(f"    {ev.gene:>3} @ {ev.species}  [{ev.kind}]  子映射: {ch}")

    print("\n[3] 逐边丢失计数（S: 路径长度-1；D: 路径长度不减 1）")
    for (parent, child), lost in result.edge_losses.items():
        count = result.edge_loss_counts[(parent, child)]
        detail = ", ".join(lost) if lost else "-"
        print(f"    边 {parent:>2} -> {child:<3} 丢失 {count} 支: {detail}")

    print(
        "\n[4] 合计: "
        f"S={result.n_speciations}  D={result.n_duplications}  "
        f"丢失={result.n_losses}"
    )

    # ---- 可手工复验的事件路径：直接从物种树结构重算路径，不读结果表 ----
    print("\n[5] 可复验事件路径（物种树上的真实路径，对照 [2][3] 手工核对）")
    gnode = {n.label: n for n in gene_tree.postorder}
    for ev in result.events:
        g = gnode[ev.gene]
        for child_gene, child_species in ev.children:
            chain = species_tree.path(
                species_tree.label_to_node[ev.species],
                species_tree.label_to_node[child_species],
            )
            k = len(chain) - 1
            expect = k if ev.kind == EVENT_DUPLICATION else k - 1
            arrow = " -> ".join(n.label for n in chain)
            print(
                f"    {ev.gene}({ev.kind}) -> {child_gene}: {arrow} | "
                f"路径长度={k}, 本边丢失应={expect}"
            )

    print(
        "\n    手算要点: g1=LCA(A,A)=A(D); "
        "g2=LCA(A,C)=X(S), 去 A 路上丢 B、去 C 路上丢 D; "
        "h=LCA(A,X)=X 且 g2 已映到 X, 故 h 为 D, 去 A 不减 1 => 再丢 Z 整支和 B; "
        "g0=LCA(X,E)=R(S), 无丢失。合计 S=2 D=2 丢失=4。"
    )

    # ---- 实际触发的拒绝边界 ----
    print("\n" + "=" * 68)
    print("拒绝边界：基因叶映射到物种内部节点（非叶映射必须被拒绝）")
    print("=" * 68)
    bad_map = dict(LEAF_MAP)
    bad_map["c"] = "Z"  # Z 是物种内部节点 (C,D)Z，不是叶子
    print(f"leaf_map = {bad_map}")
    try:
        reconcile(species_tree, gene_tree, bad_map)
    except MappingError as exc:
        print(f"如预期抛出 MappingError: {exc}")
        print("（归因未执行，无任何部分结果产生——reconcile 是纯函数。）")
    else:
        raise AssertionError("非法的非叶映射未被拒绝")


if __name__ == "__main__":
    main()
