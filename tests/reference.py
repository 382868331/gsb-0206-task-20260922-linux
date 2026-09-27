"""独立参考实现（仅测试使用）。

注意：本参考不导入、不调用 reconciliation 核心的任何函数来产生期望值，
以“祖先集合 + 朴素路径爬升”的另一种写法独立实现 LCA 复演。
核心测试拿它的输出做交叉核对；手工算例另由硬编码期望值保证。
"""

from __future__ import annotations

from itertools import product


class RefTree:
    def __init__(self, children, root):
        self.children = dict(children)
        self.root = root
        self.parent = {}
        for p, (a, b) in self.children.items():
            self.parent[a] = p
            self.parent[b] = p
        self.nodes = {root}
        frontier = [root]
        while frontier:
            n = frontier.pop()
            if n in self.children:
                frontier.extend(self.children[n])
                self.nodes.update(self.children[n])
        self.leaves = frozenset(n for n in self.nodes if n not in self.children)

    def ancestors(self, n):
        """含自身的祖先集合（朴素实现，不做深度优化）。"""
        out = []
        while True:
            out.append(n)
            if n == self.root:
                return frozenset(out)
            n = self.parent[n]

    def lca(self, a, b):
        anc = self.ancestors(a)
        x = b
        while x not in anc:
            x = self.parent[x]
        return x

    def path(self, a, b):
        """a 为 b 祖先时返回 [a..b]。"""
        up = []
        x = b
        while x != a:
            up.append(x)
            x = self.parent[x]
        return [a] + up[::-1]


def reference_reconcile(sp_children, sp_root, g_children, g_root, leaf_map):
    """返回 (mapping, {gene: event}, {(p,c): loss_count}, {(p,c): [lost_siblings]}, total)。"""
    sp = RefTree(sp_children, sp_root)
    gn = RefTree(g_children, g_root)

    # 迭代后序
    order = []
    stack = [gn.root]
    seen = set()
    while stack:
        n = stack.pop()
        if n in gn.children:
            if n in seen:
                order.append(n)
            else:
                seen.add(n)
                a, b = gn.children[n]
                stack.append(n)
                stack.append(b)
                stack.append(a)
    mapping = dict(leaf_map)
    for n in order:
        a, b = gn.children[n]
        mapping[n] = sp.lca(mapping[a], mapping[b])

    events = {}
    for n in order:
        a, b = gn.children[n]
        events[n] = (
            "duplication"
            if mapping[a] == mapping[n] or mapping[b] == mapping[n]
            else "speciation"
        )

    edge_loss = {}
    edge_sibs = {}
    for n in order:
        dup = events[n] == "duplication"
        for c in gn.children[n]:
            path = sp.path(mapping[n], mapping[c])
            d = len(path) - 1
            sibs = []
            for i in range(d):
                ca, cb = sp.children[path[i]]
                nxt = path[i + 1]
                sib = cb if nxt == ca else ca
                if dup or i > 0:
                    sibs.append(sib)
            edge_loss[(n, c)] = d if dup else d - 1
            edge_sibs[(n, c)] = sibs

    total = sum(edge_loss.values())
    return mapping, events, edge_loss, edge_sibs, total


def reference_best_candidates(
    sp_children, sp_root, g_children, g_root, candidates, dup_weight=1, loss_weight=1
):
    """独立穷举候选完整映射并打分（只依赖本参考的 reference_reconcile）。

    candidates: {基因叶: (候选物种叶, ...)}。
    返回 (最优代价, [(选定元组（按叶标签字典序）, 复制数, 丢失数), ...])，
    并列解按选定元组字典序排列。
    """
    leaves = sorted(candidates)
    best = None
    tied = []
    for chosen in product(*(candidates[g] for g in leaves)):
        leaf_map = dict(zip(leaves, chosen))
        _, events, edge_loss, _, _ = reference_reconcile(
            sp_children, sp_root, g_children, g_root, leaf_map
        )
        dups = sum(1 for e in events.values() if e == "duplication")
        losses = sum(edge_loss.values())
        cost = dup_weight * dups + loss_weight * losses
        if best is None or cost < best:
            best = cost
            tied = [(tuple(chosen), dups, losses)]
        elif cost == best:
            tied.append((tuple(chosen), dups, losses))
    tied.sort(key=lambda item: item[0])
    return best, tied
