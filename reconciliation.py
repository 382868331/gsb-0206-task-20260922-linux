"""基因树与物种树的事件归因（LCA 复演）核心库。

仅依赖 Python 3.14 标准库。公开接口见模块文档与 README：

- build_tree(children, root=None, *, name=..., max_nodes=200) -> Tree
- reconcile(species_tree, gene_tree, leaf_mapping, *, max_nodes=200) -> ReconciliationResult
- 数据类型: Tree / Event / LostLineage / EdgeLoss / ReconciliationResult
- 错误类型:
    ReconciliationError
      ├─ InvalidInputError     # 类型/取值非法（含 bool 充整数、NaN/Infinity、超节点上限）
      ├─ TreeStructureError    # 树结构非法（多父、循环、不可达、非二叉等），可定位到节点
      └─ LeafMappingError      # 基因叶→物种叶映射非法（含非叶映射），可定位到基因叶

契约摘要（以 README 为准）：
- 两棵树均为有根二叉树，节点标签为非空 str，单棵树节点数 <= max_nodes（默认 200）。
- 基因内节点映射 = 两个子映射在物种树中的 LCA；
  任一子映射等于父映射 => duplication，否则 speciation。
- 基因父子边 (g, c) 的丢失数：物种树路径长度 d（父映射到子映射的边数）；
  g 为 speciation 时丢失 d-1，g 为 duplication 时丢失 d（不减 1）。
- “具体缺失支系”：沿物种树路径每个分叉点上、不在路径上的兄弟支系；
  speciation 跳过起点处的兄弟（该支系被另一个基因子支占据）。

本模块为纯函数式、无状态：结果全部在校验通过后一次性构造，
非法输入在任何结果产生前抛错，不存在部分变更。
"""

from __future__ import annotations

import math
from collections import deque
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Mapping, Sequence

__all__ = [
    "Tree",
    "Event",
    "LostLineage",
    "EdgeLoss",
    "ReconciliationResult",
    "CandidateReconciliationResult",
    "build_tree",
    "reconcile",
    "reconcile_with_candidates",
    "ReconciliationError",
    "InvalidInputError",
    "TreeStructureError",
    "LeafMappingError",
    "CandidateMappingError",
    "SPECIATION",
    "DUPLICATION",
]

SPECIATION = "speciation"
DUPLICATION = "duplication"

# 题面给定的接口边界（单棵树节点数）。
DEFAULT_MAX_NODES = 200


class ReconciliationError(ValueError):
    """本库所有错误的基类（同时是 ValueError）。"""


class InvalidInputError(ReconciliationError):
    """参数类型或取值非法（bool 充整数、NaN/Infinity、超过节点上限等）。"""


class TreeStructureError(ReconciliationError):
    """树结构非法；:attr:`tree` 与 :attr:`node` 用于定位。"""

    def __init__(self, message: str, *, tree: str | None = None, node: str | None = None):
        self.tree = tree
        self.node = node
        where = []
        if tree is not None:
            where.append(f"树 {tree!r}")
        if node is not None:
            where.append(f"节点 {node!r}")
        prefix = f"[{' / '.join(where)}] " if where else ""
        super().__init__(prefix + message)


class LeafMappingError(ReconciliationError):
    """基因叶映射非法；:attr:`gene_leaf` 用于定位。"""

    def __init__(self, message: str, *, gene_leaf: str | None = None):
        self.gene_leaf = gene_leaf
        prefix = f"[基因叶 {gene_leaf!r}] " if gene_leaf is not None else ""
        super().__init__(prefix + message)


class CandidateMappingError(LeafMappingError):
    """候选物种叶集合非法（空集、含非叶/未知节点等）；仍是 LeafMappingError。"""


def _check_max_nodes(max_nodes: object) -> int:
    """整数参数校验：拒绝 bool、拒绝非整数、拒绝 NaN/Infinity（浮点一律拒绝）。"""
    if isinstance(max_nodes, bool) or not isinstance(max_nodes, int):
        if isinstance(max_nodes, float):
            if math.isnan(max_nodes):
                raise InvalidInputError("max_nodes 不能为 NaN")
            if math.isinf(max_nodes):
                raise InvalidInputError("max_nodes 不能为 Infinity/-Infinity")
        raise InvalidInputError(
            f"max_nodes 必须为 int 且不能是 bool，实际类型为 {type(max_nodes).__name__}"
        )
    if max_nodes < 1:
        raise InvalidInputError(f"max_nodes 必须为正整数，实际为 {max_nodes}")
    if max_nodes > DEFAULT_MAX_NODES:
        raise InvalidInputError(
            f"max_nodes 超过题面接口上限 {DEFAULT_MAX_NODES}，实际为 {max_nodes}"
        )
    return max_nodes


def _require_label(value: object, where: str) -> str:
    if not isinstance(value, str) or isinstance(value, bool):
        raise InvalidInputError(f"{where} 必须为非空 str，实际类型为 {type(value).__name__}")
    if value == "":
        raise InvalidInputError(f"{where} 不能为空字符串")
    return value


class Tree:
    """已校验的有根二叉树（不可变语义，请勿修改其属性）。

    Attributes:
        name: 树名（仅用于错误定位，如 "species" / "gene"）。
        root: 根节点标签。
        children: 内节点 -> 两个子节点标签（有序，按构造顺序）。
        parent: 子节点 -> 父节点（根不在表中）。
        depth: 节点 -> 根到该节点的边数。
        leaves: 全部叶节点标签。
        internal: 全部内节点标签。
        leaf_descendants: 节点 -> 其后裔叶节点集合。
    """

    __slots__ = (
        "name",
        "root",
        "children",
        "parent",
        "depth",
        "leaves",
        "internal",
        "leaf_descendants",
    )

    def __init__(
        self,
        *,
        name: str,
        root: str,
        children: Mapping[str, tuple[str, str]],
        parent: Mapping[str, str],
        depth: Mapping[str, int],
        leaves: frozenset[str],
        internal: frozenset[str],
        leaf_descendants: Mapping[str, frozenset[str]],
    ):
        self.name = name
        self.root = root
        self.children = dict(children)
        self.parent = dict(parent)
        self.depth = dict(depth)
        self.leaves = leaves
        self.internal = internal
        self.leaf_descendants = dict(leaf_descendants)

    @property
    def nodes(self) -> frozenset[str]:
        return self.leaves | self.internal

    def __repr__(self) -> str:
        return (
            f"Tree(name={self.name!r}, root={self.root!r}, "
            f"{len(self.leaves)} leaves, {len(self.internal)} internal)"
        )


def build_tree(
    children: Mapping[str, Sequence[str]],
    root: str | None = None,
    *,
    name: str = "tree",
    max_nodes: int = DEFAULT_MAX_NODES,
) -> Tree:
    """由 ``{内节点: (左子, 右子)}`` 构造并校验一棵有根二叉树。

    叶子 = 作为子节点出现、但自身不在 ``children`` 键中的节点。
    允许 ``children={}`` 且显式给定 ``root`` 的单叶（平凡）树。

    Args:
        children: 内节点到恰好两个**不同**子节点的映射。
        root: 根标签；省略时要求“恰有一个不作为任何子节点出现的键”，否则报错。
        name: 树名，用于错误信息定位。
        max_nodes: 节点数上限（正整数，<= 200，默认 200；拒绝 bool）。

    Raises:
        InvalidInputError: 类型/取值非法或超过节点上限。
        TreeStructureError: 多父、循环、不可达、孤键、非二叉等结构错误（可定位节点）。
    """
    limit = _check_max_nodes(max_nodes)
    _require_label(name, "树名 name")

    if not isinstance(children, Mapping):
        raise InvalidInputError(
            f"children 必须为 Mapping，实际类型为 {type(children).__name__}"
        )

    clean_children: dict[str, tuple[str, str]] = {}
    universe: set[str] = set()

    if not children:
        # 平凡单叶树：仅允许显式给定 root。
        if root is None:
            raise InvalidInputError(
                f"children 为空（单叶树）时必须显式给定 root 参数（树 {name!r}）"
            )
        only = _require_label(root, "root")
        return Tree(
            name=name,
            root=only,
            children={},
            parent={},
            depth={only: 0},
            leaves=frozenset({only}),
            internal=frozenset(),
            leaf_descendants={only: frozenset({only})},
        )

    for key, pair in children.items():
        node = _require_label(key, "children 的键")
        if not isinstance(pair, Sequence) or isinstance(pair, (str, bytes)):
            raise TreeStructureError(
                "子节点必须是长度为 2 的序列 (list/tuple)", tree=name, node=node
            )
        if len(pair) != 2:
            raise TreeStructureError(
                f"二叉树内节点必须恰有 2 个子节点，实际有 {len(pair)} 个",
                tree=name,
                node=node,
            )
        c0 = _require_label(pair[0], f"节点 {node!r} 的第 1 个子节点")
        c1 = _require_label(pair[1], f"节点 {node!r} 的第 2 个子节点")
        if c0 == c1:
            raise TreeStructureError(
                f"两个子节点不能相同（均为 {c0!r}）", tree=name, node=node
            )
        if node in (c0, c1):
            # 直接自环，三色遍历也能查，这里给出更明确的定位信息。
            raise TreeStructureError("检测到自环（节点是自己的子节点）", tree=name, node=node)
        clean_children[node] = (c0, c1)
        universe.add(node)
        universe.add(c0)
        universe.add(c1)

    internal = frozenset(clean_children)

    # 推断 / 校验根。
    if root is None:
        if not clean_children:
            raise InvalidInputError(
                f"children 为空（单叶树）时必须显式给定 root 参数（树 {name!r}）"
            )
        parentless = [k for k in clean_children if not any(k in pair for pair in clean_children.values())]
        if len(parentless) == 0:
            raise TreeStructureError(
                "无法推断根：所有内节点都作为某个子节点出现（可能存在环）", tree=name
            )
        if len(parentless) > 1:
            raise TreeStructureError(
                f"无法唯一推断根：存在多个无父内节点 {sorted(parentless)}，请显式给定 root",
                tree=name,
            )
        root_label = parentless[0]
    else:
        root_label = _require_label(root, "root")
        if root_label not in universe:
            raise TreeStructureError(
                f"root {root_label!r} 未出现在 children 的任何节点中", tree=name, node=root_label
            )

    # 唯一性父节点检查（一个节点最多一个父）。
    parent: dict[str, str] = {}
    for node, (c0, c1) in clean_children.items():
        for child in (c0, c1):
            if child in parent:
                raise TreeStructureError(
                    f"节点有多个父节点：{parent[child]!r} 与 {node!r}",
                    tree=name,
                    node=child,
                )
            parent[child] = node

    # 三色遍历：0=未访问, 1=在当前递归栈上, 2=完成。可定位环与不可达节点。
    state: dict[str, int] = {n: 0 for n in universe}
    depth: dict[str, int] = {}
    postorder: list[str] = []
    # (node, 是否已展开)
    stack: list[tuple[str, bool]] = [(root_label, False)]
    state[root_label] = 1
    depth[root_label] = 0
    while stack:
        node, expanded = stack.pop()
        if expanded:
            state[node] = 2
            postorder.append(node)
            continue
        pair = clean_children.get(node)
        if pair is None:
            # 叶子：完成并记入后序。
            state[node] = 2
            postorder.append(node)
            continue
        stack.append((node, True))
        for child in pair:
            s = state.get(child, 0)
            if s == 1:
                raise TreeStructureError(
                    f"检测到环：{node!r} -> {child!r} 的回边", tree=name, node=child
                )
            if s == 0:
                state[child] = 1
                depth[child] = depth[node] + 1
                stack.append((child, False))
            # s == 2 不可能出现（多父已被拦截），不特殊处理。

    unreachable = sorted(n for n, s in state.items() if s != 2)
    if unreachable:
        raise TreeStructureError(
            f"存在从根不可达的节点（游离子树/环分量）：{unreachable}",
            tree=name,
            node=unreachable[0],
        )

    leaves = frozenset(n for n in universe if n not in internal)

    node_count = len(universe)
    if node_count > limit:
        raise InvalidInputError(
            f"树 {name!r} 节点数为 {node_count}，超过上限 {limit}"
        )

    # 每个节点的后裔叶集合（用于标出具体缺失支系）。
    leaf_descendants: dict[str, frozenset[str]] = {}
    for node in postorder:
        pair = clean_children.get(node)
        if pair is None:
            leaf_descendants[node] = frozenset({node})
        else:
            leaf_descendants[node] = leaf_descendants[pair[0]] | leaf_descendants[pair[1]]

    return Tree(
        name=name,
        root=root_label,
        children=clean_children,
        parent=parent,
        depth=depth,
        leaves=leaves,
        internal=internal,
        leaf_descendants=leaf_descendants,
    )


@dataclass(frozen=True)
class Event:
    """一个基因内节点上的事件。

    Attributes:
        gene_node: 基因内节点标签。
        event_type: ``speciation`` 或 ``duplication``。
        mapped_species: 该基因节点映射到的物种树节点。
        child_genes: 两个基因子节点（按树构造顺序）。
        child_maps: 两个子节点各自映射到的物种树节点。
    """

    gene_node: str
    event_type: str
    mapped_species: str
    child_genes: tuple[str, str]
    child_maps: tuple[str, str]


@dataclass(frozen=True)
class LostLineage:
    """一条具体缺失支系：物种树节点 :attr:`species_node` 所代表的整个支系缺失。"""

    species_node: str
    descendant_leaves: frozenset[str]


@dataclass(frozen=True)
class EdgeLoss:
    """一条基因父子边上的丢失计数与具体缺失支系。

    Attributes:
        parent_gene / child_gene: 基因树边的两端。
        parent_species / child_species: 两端映射到的物种树节点。
        path_length: 物种树上父映射到子映射的路径边数 d。
        loss_count: 丢失数（speciation 父节点为 d-1，duplication 为 d）。
        lost_lineages: 具体缺失支系（顺序沿物种树路径自上而下）。
    """

    parent_gene: str
    child_gene: str
    parent_species: str
    child_species: str
    path_length: int
    loss_count: int
    lost_lineages: tuple[LostLineage, ...]


@dataclass(frozen=True)
class ReconciliationResult:
    """复演结果（不可变）。

    Attributes:
        mapping: 全部基因节点（叶+内节点）-> 物种树节点。
        events: 事件表（按基因树后序，子节点事件在前）。
        edge_losses: 逐基因父子边的丢失计数（后序父节点、先左后右）。
        losses_by_lineage: 各缺失物种支系 -> 被计为丢失的次数（按支系标签排序）。
        duplications / speciations / total_losses: 累计计数。
    """

    mapping: Mapping[str, str]
    events: tuple[Event, ...]
    edge_losses: tuple[EdgeLoss, ...]
    losses_by_lineage: Mapping[str, int]
    duplications: int
    speciations: int
    total_losses: int


def _lca(species: Tree, a: str, b: str) -> str:
    """物种树两节点的最近共同祖先（按深度抬升）。"""
    depth = species.depth
    parent = species.parent
    while depth[a] > depth[b]:
        a = parent[a]
    while depth[b] > depth[a]:
        b = parent[b]
    while a != b:
        a = parent[a]
        b = parent[b]
    return a


def reconcile(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_mapping: Mapping[str, str],
    *,
    max_nodes: int = DEFAULT_MAX_NODES,
) -> ReconciliationResult:
    """对基因树与物种树做 LCA 复演并归因事件与丢失。

    Args:
        species_tree: ``build_tree`` 构造的物种树。
        gene_tree: ``build_tree`` 构造的基因树。
        leaf_mapping: 每个基因叶 -> 一个**物种叶**；同一物种叶可对应多个基因叶。
        max_nodes: 保留的接口边界参数（树在构造时已受限，此处再次校验参数本身）。

    Returns:
        ReconciliationResult：映射、事件表、逐边计数与累计值。

    Raises:
        InvalidInputError: 参数类型非法（含 bool 充整数、NaN/Infinity）。
        TreeStructureError: 传入对象不是本库 Tree。
        LeafMappingError: 映射不完整、含未知键、目标不是物种叶（非叶映射）等。
    """
    _check_max_nodes(max_nodes)  # 显式拒绝 bool/NaN 等
    if not isinstance(species_tree, Tree):
        raise InvalidInputError(
            f"species_tree 必须为 build_tree 产出的 Tree，实际为 {type(species_tree).__name__}"
        )
    if not isinstance(gene_tree, Tree):
        raise InvalidInputError(
            f"gene_tree 必须为 build_tree 产出的 Tree，实际为 {type(gene_tree).__name__}"
        )
    if not isinstance(leaf_mapping, Mapping):
        raise InvalidInputError(
            f"leaf_mapping 必须为 Mapping，实际类型为 {type(leaf_mapping).__name__}"
        )

    # ---- 映射校验（全部在任何计算之前完成） ----
    gene_leaves = gene_tree.leaves
    raw_mapping: dict[str, str] = {}
    for key, value in leaf_mapping.items():
        # 键本身不是合法标签时无法定位到具体基因叶，gene_leaf=None。
        if not isinstance(key, str) or isinstance(key, bool) or key == "":
            raise LeafMappingError(
                f"leaf_mapping 的键必须是非空 str，实际为 {key!r}（{type(key).__name__}）"
            )
        g = key
        if not isinstance(value, str) or isinstance(value, bool) or value == "":
            raise LeafMappingError(
                f"基因叶 {g!r} 的映射目标必须是非空 str，实际为 {value!r}"
                f"（{type(value).__name__}）",
                gene_leaf=g,
            )
        s = value
        if g not in gene_leaves:
            if g in gene_tree.internal:
                raise LeafMappingError(
                    "映射键必须是基因叶，不能是基因内节点", gene_leaf=g
                )
            raise LeafMappingError("不是基因树中的节点", gene_leaf=g)
        if s not in species_tree.nodes:
            raise LeafMappingError(
                f"映射目标 {s!r} 不是物种树中的节点", gene_leaf=g
            )
        if s in species_tree.internal:
            # 题面硬性要求：拒绝非叶映射。
            raise LeafMappingError(
                f"映射目标 {s!r} 是物种内节点，基因叶只能映射到物种叶（拒绝非叶映射）",
                gene_leaf=g,
            )
        raw_mapping[g] = s

    missing = sorted(gene_leaves - raw_mapping.keys())
    if missing:
        raise LeafMappingError(f"存在未给出映射的基因叶：{missing}", gene_leaf=missing[0])

    # ---- 后序求 LCA 映射 ----
    g_children = gene_tree.children
    postorder: list[str] = []
    stack: list[tuple[str, bool]] = [(gene_tree.root, False)]
    while stack:
        node, expanded = stack.pop()
        pair = g_children.get(node)
        if pair is None:
            continue
        if expanded:
            postorder.append(node)
            continue
        stack.append((node, True))
        stack.append((pair[1], False))
        stack.append((pair[0], False))

    mapping: dict[str, str] = dict(raw_mapping)
    for node in postorder:
        c0, c1 = g_children[node]
        mapping[node] = _lca(species_tree, mapping[c0], mapping[c1])

    # ---- 事件归因 ----
    events: list[Event] = []
    event_type: dict[str, str] = {}
    for node in postorder:
        c0, c1 = g_children[node]
        m = mapping[node]
        m0, m1 = mapping[c0], mapping[c1]
        etype = DUPLICATION if (m0 == m or m1 == m) else SPECIATION
        event_type[node] = etype
        events.append(
            Event(
                gene_node=node,
                event_type=etype,
                mapped_species=m,
                child_genes=(c0, c1),
                child_maps=(m0, m1),
            )
        )

    # ---- 逐边丢失与具体缺失支系 ----
    s_children = species_tree.children
    s_parent = species_tree.parent
    edge_losses: list[EdgeLoss] = []
    lineage_counts: dict[str, int] = {}

    for node in postorder:
        c0, c1 = g_children[node]
        m_parent = mapping[node]
        parent_is_dup = event_type[node] == DUPLICATION
        for child in (c0, c1):
            m_child = mapping[child]
            # 收集 m_parent -> m_child 的物种树路径（自上而下）。
            upwards: deque[str] = deque()
            cur = m_child
            while cur != m_parent:
                upwards.appendleft(cur)
                cur = s_parent[cur]
            path = [m_parent, *upwards]  # m_parent 是 m_child 的祖先（LCA 性质保证）
            d = len(path) - 1

            lost: list[LostLineage] = []
            # 路径上每个分叉点 p_i -> p_{i+1}，其兄弟支系即潜在缺失支系。
            for i in range(d):
                p_i, p_next = path[i], path[i + 1]
                a, b = s_children[p_i]
                sibling = b if p_next == a else a
                # speciation：起点处的兄弟支系由另一个基因子支占据，不计丢失。
                if (not parent_is_dup) and i == 0:
                    continue
                lineage = LostLineage(
                    species_node=sibling,
                    descendant_leaves=species_tree.leaf_descendants[sibling],
                )
                lost.append(lineage)
                lineage_counts[sibling] = lineage_counts.get(sibling, 0) + 1

            loss_count = d if parent_is_dup else d - 1
            # 与枚举结果自检一致（纯整数算术，不允许偏差）。
            assert loss_count == len(lost), (node, child, loss_count, len(lost))
            edge_losses.append(
                EdgeLoss(
                    parent_gene=node,
                    child_gene=child,
                    parent_species=m_parent,
                    child_species=m_child,
                    path_length=d,
                    loss_count=loss_count,
                    lost_lineages=tuple(lost),
                )
            )

    total_losses = sum(e.loss_count for e in edge_losses)
    duplications = sum(1 for e in events if e.event_type == DUPLICATION)
    speciations = sum(1 for e in events if e.event_type == SPECIATION)

    return ReconciliationResult(
        mapping=dict(sorted(mapping.items())),
        events=tuple(events),
        edge_losses=tuple(edge_losses),
        losses_by_lineage=dict(sorted(lineage_counts.items())),
        duplications=duplications,
        speciations=speciations,
        total_losses=total_losses,
    )


# ---------------------------------------------------------------------------
# 候选物种叶集合：联合选择每个基因叶的归属，最小化加权 duplication+loss。
# ---------------------------------------------------------------------------

# DP 松弛工作量默认预算（sum over 内节点 of |maps(左子)| * |maps(右子)|）。
DEFAULT_MAX_RELAXATIONS = 4_000_000


@dataclass(frozen=True)
class CandidateReconciliationResult:
    """候选归属搜索结果（不可变）。

    Attributes:
        complete: True 表示已证明最优；False 表示搜索预算耗尽（incomplete），
            此时返回的是最优已知映射，optimal_cost 为 None。
        assignments: 每个基因叶最终选定的物种叶。
        optimal_cost: 已证明的最优目标值（complete 时非 None）。
        best_cost: 返回映射实际达到的目标值（可验证，complete 时 == optimal_cost）。
        reconciliation: 对 assignments 做标准 LCA 复演的完整事件账本。
        relaxation_work: 实际消耗的 DP 松弛工作量。
    """

    complete: bool
    assignments: Mapping[str, str]
    optimal_cost: int | None
    best_cost: int
    reconciliation: ReconciliationResult
    relaxation_work: int


class _BudgetExceeded(Exception):
    """内部控制流：DP 松弛工作量超过预算。"""

    def __init__(self, work: int):
        self.work = work
        super().__init__(f"relaxation work exceeded budget: {work}")


def _check_weight(value: object, name: str) -> int:
    """权重校验：非负 int，拒绝 bool 与一切非 int（含 NaN/Infinity）。"""
    if isinstance(value, bool) or not isinstance(value, int):
        if isinstance(value, float):
            if math.isnan(value):
                raise InvalidInputError(f"{name} 不能为 NaN")
            if math.isinf(value):
                raise InvalidInputError(f"{name} 不能为 Infinity/-Infinity")
        raise InvalidInputError(
            f"{name} 必须为非负 int 且不能是 bool，实际类型为 {type(value).__name__}"
        )
    if value < 0:
        raise InvalidInputError(f"{name} 必须为非负整数，实际为 {value}")
    return value


def _check_max_relaxations(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidInputError(
            f"max_relaxations 必须为正 int 且不能是 bool，实际类型为 {type(value).__name__}"
        )
    if value < 1:
        raise InvalidInputError(f"max_relaxations 必须为正整数，实际为 {value}")
    return value


def _normalize_candidates(
    species_tree: Tree, gene_tree: Tree, leaf_candidates: Mapping[str, object]
) -> dict[str, tuple[str, ...]]:
    """校验并规范化候选输入：每个基因叶 -> 去重后排序的候选物种叶元组。

    单条 str 视为单候选（与旧行为一致）；其余必须是可迭代的非空 str 集合。
    """
    if not isinstance(leaf_candidates, Mapping):
        raise InvalidInputError(
            f"leaf_candidates 必须为 Mapping，实际类型为 {type(leaf_candidates).__name__}"
        )
    gene_leaves = gene_tree.leaves
    normalized: dict[str, tuple[str, ...]] = {}
    for key, value in leaf_candidates.items():
        if not isinstance(key, str) or isinstance(key, bool) or key == "":
            raise CandidateMappingError(
                f"leaf_candidates 的键必须是非空 str，实际为 {key!r}（{type(key).__name__}）"
            )
        g = key
        if g not in gene_leaves:
            if g in gene_tree.internal:
                raise CandidateMappingError(
                    "候选键必须是基因叶，不能是基因内节点", gene_leaf=g
                )
            raise CandidateMappingError("不是基因树中的节点", gene_leaf=g)
        if isinstance(value, str):
            raw = [value]
        else:
            if not isinstance(value, Iterable) or isinstance(value, (bytes,)):
                raise CandidateMappingError(
                    f"候选必须是非空 str 或 str 的可迭代集合，实际为 {value!r}"
                    f"（{type(value).__name__}）",
                    gene_leaf=g,
                )
            raw = list(value)
        if not raw:
            raise CandidateMappingError("候选集合不能为空", gene_leaf=g)
        deduped = sorted(set(raw))
        for s in deduped:
            if not isinstance(s, str) or isinstance(s, bool) or s == "":
                raise CandidateMappingError(
                    f"候选目标必须是非空 str，实际为 {s!r}（{type(s).__name__}）",
                    gene_leaf=g,
                )
            if s not in species_tree.nodes:
                raise CandidateMappingError(
                    f"候选目标 {s!r} 不是物种树中的节点", gene_leaf=g
                )
            if s in species_tree.internal:
                raise CandidateMappingError(
                    f"候选目标 {s!r} 是物种内节点，候选必须全是物种叶",
                    gene_leaf=g,
                )
        normalized[g] = tuple(deduped)
    missing = sorted(gene_leaves - normalized.keys())
    if missing:
        raise CandidateMappingError(
            f"存在未给出候选的基因叶：{missing}", gene_leaf=missing[0]
        )
    return normalized


def _gene_internal_postorder(gene_tree: Tree) -> list[str]:
    postorder: list[str] = []
    stack: list[tuple[str, bool]] = [(gene_tree.root, False)]
    while stack:
        node, expanded = stack.pop()
        pair = gene_tree.children.get(node)
        if pair is None:
            continue
        if expanded:
            postorder.append(node)
            continue
        stack.append((node, True))
        stack.append((pair[1], False))
        stack.append((pair[0], False))
    return postorder


def _dp_optimal_cost(
    species_tree: Tree,
    gene_tree: Tree,
    postorder: list[str],
    candidates: Mapping[str, tuple[str, ...]],
    fixed: Mapping[str, str],
    dup_weight: int,
    loss_weight: int,
    budget: int,
) -> tuple[int, int]:
    """树 DP：给定部分叶已固定，求子树最小加权成本与松弛工作量。

    状态 = 基因节点映射到的物种节点（只保留可达 LCA，故每节点状态数
    以物种树节点数为界）；转移时把事件成本（本节点）与两条子边丢失
    一并计入。预算耗尽抛 _BudgetExceeded。
    """
    g_children = gene_tree.children
    s_depth = species_tree.depth
    work = 0
    # node -> {species_node: min_cost}
    table: dict[str, dict[str, int]] = {}
    for leaf in gene_tree.leaves:
        if leaf in fixed:
            table[leaf] = {fixed[leaf]: 0}
        else:
            table[leaf] = {c: 0 for c in candidates[leaf]}
    for node in postorder:
        c0, c1 = g_children[node]
        m0, m1 = table[c0], table[c1]
        work += len(m0) * len(m1)
        if work > budget:
            raise _BudgetExceeded(work)
        cell: dict[str, int] = {}
        for t0, cost0 in m0.items():
            for t1, cost1 in m1.items():
                s = _lca(species_tree, t0, t1)
                is_dup = t0 == s or t1 == s
                edge0 = s_depth[t0] - s_depth[s]
                edge1 = s_depth[t1] - s_depth[s]
                cost = cost0 + cost1
                if is_dup:
                    cost += dup_weight + loss_weight * (edge0 + edge1)
                else:
                    cost += loss_weight * (edge0 + edge1 - 2)
                prev = cell.get(s)
                if prev is None or cost < prev:
                    cell[s] = cost
        table[node] = cell
    return min(table[gene_tree.root].values()), work


def reconcile_with_candidates(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_candidates: Mapping[str, object],
    *,
    dup_weight: int = 1,
    loss_weight: int = 1,
    max_nodes: int = DEFAULT_MAX_NODES,
    max_relaxations: int = DEFAULT_MAX_RELAXATIONS,
) -> CandidateReconciliationResult:
    """为带候选物种叶集合的基因叶联合选择归属，最小化加权事件成本。

    目标：min  dup_weight * duplications + loss_weight * total_losses；
    平局按基因叶 ID（字典序）对应的物种 ID 序列取字典序最小者。
    选择完成后仍按原 LCA 规则复演（不能逐叶独立选最近物种——
    候选选择会改变祖先 LCA，必须联合优化）。

    Args:
        species_tree / gene_tree: ``build_tree`` 构造的树。
        leaf_candidates: 每个基因叶 -> 单个物种叶 str（与旧输入等价）
            或候选物种叶的非空可迭代集合（重复候选自动去重）。
        dup_weight / loss_weight: 非负 int 权重（拒绝 bool/浮点/负数）。
        max_nodes: 接口边界参数（同 reconcile）。
        max_relaxations: DP 松弛工作量预算；耗尽时不报错，而是返回
            ``complete=False`` 的结果（最优已知映射 + 可验证事件账本）。

    Returns:
        CandidateReconciliationResult。

    Raises:
        InvalidInputError: 参数类型/取值非法。
        CandidateMappingError: 候选集合非法（LeafMappingError 子类）。
    """
    _check_max_nodes(max_nodes)
    dup_w = _check_weight(dup_weight, "dup_weight")
    loss_w = _check_weight(loss_weight, "loss_weight")
    budget = _check_max_relaxations(max_relaxations)
    if not isinstance(species_tree, Tree):
        raise InvalidInputError(
            f"species_tree 必须为 build_tree 产出的 Tree，实际为 {type(species_tree).__name__}"
        )
    if not isinstance(gene_tree, Tree):
        raise InvalidInputError(
            f"gene_tree 必须为 build_tree 产出的 Tree，实际为 {type(gene_tree).__name__}"
        )
    candidates = _normalize_candidates(species_tree, gene_tree, leaf_candidates)
    postorder = _gene_internal_postorder(gene_tree)

    def fallback_result(work: int) -> CandidateReconciliationResult:
        # 预算耗尽：最优已知映射取每叶字典序最小候选，账本仍可完整验证。
        assignments = {g: cs[0] for g, cs in candidates.items()}
        rec = reconcile(species_tree, gene_tree, assignments, max_nodes=max_nodes)
        best = dup_w * rec.duplications + loss_w * rec.total_losses
        return CandidateReconciliationResult(
            complete=False,
            assignments=dict(sorted(assignments.items())),
            optimal_cost=None,
            best_cost=best,
            reconciliation=rec,
            relaxation_work=work,
        )

    try:
        best_cost, work = _dp_optimal_cost(
            species_tree, gene_tree, postorder, candidates, {}, dup_w, loss_w, budget
        )
    except _BudgetExceeded as exc:
        return fallback_result(exc.work)

    # 平局消解：按基因叶 ID 顺序逐叶固定为仍能达到最优值的最小候选，
    # 精确得到“物种 ID 序列字典序最小”的最优映射（逐叶可行性由 DP 重算保证）。
    fixed: dict[str, str] = {}
    for g in sorted(gene_tree.leaves):
        cs = candidates[g]
        if len(cs) == 1:
            continue
        for c in cs:
            fixed[g] = c
            cost, _ = _dp_optimal_cost(
                species_tree, gene_tree, postorder, candidates, fixed,
                dup_w, loss_w, budget,
            )
            if cost == best_cost:
                break
        else:  # pragma: no cover - 不变式保护，理论上不可达
            raise AssertionError(f"基因叶 {g!r} 无可行候选保持最优值")

    assignments = {g: fixed.get(g, candidates[g][0]) for g in gene_tree.leaves}
    for g, cs in candidates.items():
        if len(cs) == 1:
            assignments[g] = cs[0]
    rec = reconcile(species_tree, gene_tree, assignments, max_nodes=max_nodes)
    final_cost = dup_w * rec.duplications + loss_w * rec.total_losses
    assert final_cost == best_cost, (final_cost, best_cost)
    return CandidateReconciliationResult(
        complete=True,
        assignments=dict(sorted(assignments.items())),
        optimal_cost=best_cost,
        best_cost=final_cost,
        reconciliation=rec,
        relaxation_work=work,
    )
