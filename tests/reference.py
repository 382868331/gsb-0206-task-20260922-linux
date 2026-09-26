"""独立参考实现（仅测试使用）。

注意：本参考不导入、不调用 reconciliation 核心的任何函数来产生期望值，
以“祖先集合 + 朴素路径爬升”的另一种写法独立实现 LCA 复演。
核心测试拿它的输出做交叉核对；手工算例另由硬编码期望值保证。
"""

from __future__ import annotations


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


def reference_pareto_frontier(sp_children, sp_root, g_children, g_root, leaf_map):
    """穷举全部合法映射求 (复制数, 丢失数) Pareto 前沿（独立参考，不调用被测核心）。

    合法映射：每个基因内节点映到其两个子映射的共同祖先（含自身）。
    事件与丢失口径与核心一致：任一子映射 == 父映射 => duplication
    （逐边路径长 d 计 d），否则 speciation（计 d-1）。
    返回按 (D, L) 升序的非支配计数元组，以及 {计数: 规范见证映射}：
    规范序 = 基因树后序（先左后右）节点上的物种标签序列字典序最小。
    """
    sp = RefTree(sp_children, sp_root)
    gn = RefTree(g_children, g_root)

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

    ancestors = {n: sp.ancestors(n) for n in sp.nodes}
    depth = {}
    frontier = [sp.root]
    depth[sp.root] = 0
    while frontier:
        n = frontier.pop()
        for c in sp.children.get(n, ()):
            depth[c] = depth[n] + 1
            frontier.append(c)

    candidates = {}
    for leaf in gn.leaves:
        candidates[leaf] = (leaf_map[leaf],)
    for n in order:
        a, b = gn.children[n]
        common = set()
        for ma in candidates[a]:
            for mb in candidates[b]:
                common.update(ancestors[ma] & ancestors[mb])
        candidates[n] = tuple(common)

    keys = order
    results = {}
    import itertools

    for vals in itertools.product(*[sorted(candidates[n]) for n in keys]):
        mapping = dict(zip(keys, vals))
        mapping.update(leaf_map)
        if not all(
            mapping[n] in ancestors[mapping[c]]
            for n in keys
            for c in gn.children[n]
        ):
            continue
        dups = losses = 0
        for n in keys:
            a, b = gn.children[n]
            s, ma, mb = mapping[n], mapping[a], mapping[b]
            dup = ma == s or mb == s
            da = depth[ma] - depth[s]
            db = depth[mb] - depth[s]
            dups += 1 if dup else 0
            losses += (da if dup else da - 1) + (db if dup else db - 1)
        count = (dups, losses)
        witness = tuple(mapping[n] for n in list(order) + sorted(gn.leaves - set(order)))
        # 规范见证键只含后序内节点序列；叶子固定，不含在比较键中。
        key = tuple(mapping[n] for n in order)
        if count not in results or key < results[count]:
            results[count] = key

    counts = set(results)
    front = tuple(
        sorted(
            p
            for p in counts
            if not any(q != p and q[0] <= p[0] and q[1] <= p[1] for q in counts)
        )
    )
    return front, {c: results[c] for c in front}
