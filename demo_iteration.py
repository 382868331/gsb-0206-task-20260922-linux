"""离线演示计数前沿 DP、加权选择与预算耗尽标记。

运行：python demo_iteration.py
"""

from __future__ import annotations

from itertools import product

from reconciliation import (
    build_tree,
    reconcile,
    reconcile_pareto_frontier,
    reconcile_weighted,
)


def shapes(leaves):
    if leaves == 1:
        yield "leaf"
        return
    for split in range(1, leaves):
        for left in shapes(split):
            for right in shapes(leaves - split):
                yield left, right


def make_tree(shape, prefix):
    children = {}
    leaf_id = 0
    internal_id = 0

    def build(node_shape):
        nonlocal leaf_id, internal_id
        if node_shape == "leaf":
            name = f"{prefix}_leaf{leaf_id}"
            leaf_id += 1
            return name
        left = build(node_shape[0])
        right = build(node_shape[1])
        name = f"{prefix}_n{internal_id}"
        internal_id += 1
        children[name] = left, right
        return name

    root = build(shape)
    return children, root


def postorder(children, root):
    result = []
    stack = [(root, False)]
    while stack:
        node, expanded = stack.pop()
        if node not in children:
            continue
        if expanded:
            result.append(node)
            continue
        left, right = children[node]
        stack.extend(((node, True), (right, False), (left, False)))
    return result


def ancestors(tree, node):
    result = {node}
    while node != tree.root:
        node = tree.parent[node]
        result.add(node)
    return result


def descendants(tree):
    result = {leaf: {leaf} for leaf in tree.leaves}
    for node in postorder(tree.children, tree.root):
        left, right = tree.children[node]
        result[node] = {node, *result[left], *result[right]}
    return result


def enumerate_frontier_sizes(species_tree, species_children):
    """对 2..4 个基因叶的小树独立枚举祖先映射，并统计前沿大小。"""
    desc = descendants(species_tree)
    hist = {}
    checked = 0
    for gene_leaves in range(2, 5):
        for shape in shapes(gene_leaves):
            gene_children, gene_root = make_tree(shape, "g")
            gene_tree = build_tree(gene_children, gene_root, name="gene")
            internals = postorder(gene_children, gene_root)
            clade_leaves = {}
            for gene in internals:
                left, right = gene_children[gene]
                clade_leaves[gene] = (
                    (clade_leaves.get(left) or {left})
                    | (clade_leaves.get(right) or {right})
                )
            for labels in product(sorted(species_tree.leaves), repeat=gene_leaves):
                leaf_mapping = dict(zip(sorted(gene_tree.leaves), labels))
                allowed = {
                    gene: set.intersection(
                        *(ancestors(species_tree, leaf_mapping[leaf]) for leaf in clade_leaves[gene])
                    )
                    for gene in internals
                }
                counts = set()
                for values in product(*(allowed[gene] for gene in internals)):
                    mapping = dict(leaf_mapping)
                    mapping.update(zip(internals, values))
                    duplication_count = loss_count = 0
                    valid = True
                    for gene in internals:
                        mapped = mapping[gene]
                        left_gene, right_gene = gene_children[gene]
                        left_map = mapping[left_gene]
                        right_map = mapping[right_gene]
                        speciation = False
                        if mapped in species_children:
                            side_a, side_b = species_children[mapped]
                            speciation = (
                                left_map != mapped
                                and right_map != mapped
                                and (
                                    (left_map in desc[side_a] and right_map in desc[side_b])
                                    or (left_map in desc[side_b] and right_map in desc[side_a])
                                )
                            )
                        if not speciation and not (
                            left_map in desc[mapped] and right_map in desc[mapped]
                        ):
                            valid = False
                            break
                        duplication_count += 0 if speciation else 1
                        for child_map in (left_map, right_map):
                            edge_losses = species_tree.depth[child_map] - species_tree.depth[mapped]
                            loss_count += edge_losses - (1 if speciation else 0)
                    if valid:
                        counts.add((duplication_count, loss_count))
                frontier = {
                    count
                    for count in counts
                    if not any(
                        other[0] <= count[0]
                        and other[1] <= count[1]
                        and other != count
                        for other in counts
                    )
                }
                hist[len(frontier)] = hist.get(len(frontier), 0) + 1
                checked += 1
    return checked, hist


def main():
    species_children = {"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}
    species_tree = build_tree(species_children, "r", name="species")
    gene_children = {"g1": ("a", "b"), "g2": ("c", "d"), "root": ("g1", "g2")}
    gene_tree = build_tree(gene_children, "root", name="gene")
    leaf_mapping = {"a": "A", "b": "B", "c": "A", "d": "C"}

    baseline = reconcile(species_tree, gene_tree, leaf_mapping)
    frontier = reconcile_pareto_frontier(species_tree, gene_tree, leaf_mapping)
    weighted = reconcile_weighted(species_tree, gene_tree, leaf_mapping, 2, 1)
    exhausted = reconcile_pareto_frontier(species_tree, gene_tree, leaf_mapping, budget=1)

    print("固定 LCA 入口:", baseline.duplications, "次复制,", baseline.total_losses, "次丢失")
    print(
        "Pareto 前沿:",
        [(solution.duplications, solution.losses) for solution in frontier.solutions],
        "complete =",
        frontier.complete,
    )
    print(
        "权重 (2, 1) 最优代价:",
        weighted.optimal_cost,
        [(solution.duplications, solution.losses) for solution in weighted.solutions],
    )
    print(
        "budget=1:",
        "complete =",
        exhausted.complete,
        "remaining =",
        exhausted.budget_remaining,
        "solutions =",
        len(exhausted.solutions),
    )

    checked, histogram = enumerate_frontier_sizes(species_tree, species_children)
    print("独立小树枚举:", checked, "个固定叶归属输入，前沿大小分布 =", dict(sorted(histogram.items())))
    print("结论：固定叶归属且只允许复制/丢失时，LCA 同时最小化复制数和丢失数。")


if __name__ == "__main__":
    main()
