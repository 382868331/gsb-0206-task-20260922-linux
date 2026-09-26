"""基因树与物种树的事件归因（LCA 复演）核心库。

仅依赖 Python 3.14 标准库。公开接口见模块文档与 README：

- build_tree(children, root=None, *, name=..., max_nodes=200) -> Tree
- reconcile(species_tree, gene_tree, leaf_mapping, *, max_nodes=200) -> ReconciliationResult
- reconcile_pareto_frontier(..., budget=None) -> FrontierResult
- reconcile_weighted(..., duplication_weight, loss_weight, budget=None) -> WeightedReconciliationResult
- 数据类型: Tree / Event / LostLineage / EdgeLoss / ReconciliationResult /
    FrontierSolution / FrontierResult / WeightedReconciliationResult
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
from dataclasses import dataclass
from typing import Mapping, Sequence

__all__ = [
    "Tree",
    "Event",
    "LostLineage",
    "EdgeLoss",
    "ReconciliationResult",
    "Count",
    "FrontierSolution",
    "FrontierResult",
    "WeightedReconciliationResult",
    "build_tree",
    "reconcile",
    "reconcile_pareto_frontier",
    "reconcile_weighted",
    "ReconciliationError",
    "InvalidInputError",
    "TreeStructureError",
    "LeafMappingError",
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


Count = tuple[int, int]


@dataclass(frozen=True)
class FrontierSolution:
    """一个非支配计数及其一份规范最小事件归因。"""

    duplications: int
    losses: int
    reconciliation: ReconciliationResult

    @property
    def count(self) -> Count:
        return self.duplications, self.losses


@dataclass(frozen=True)
class FrontierResult:
    """复制/丢失计数 Pareto 前沿。

    ``complete=False`` 表示搜索预算在 DP 合并阶段耗尽；此时 ``solutions``
    只是已探索状态中的非支配解，不能宣称是全局完整前沿。
    """

    solutions: tuple[FrontierSolution, ...]
    complete: bool
    budget_remaining: int | None


@dataclass(frozen=True)
class WeightedReconciliationResult:
    """正整数权重下的最优解；不同计数若代价相同会全部保留。"""

    duplication_weight: int
    loss_weight: int
    optimal_cost: int
    solutions: tuple[FrontierSolution, ...]


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


def _validate_leaf_mapping(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_mapping: Mapping[str, str],
) -> dict[str, str]:
    if not isinstance(leaf_mapping, Mapping):
        raise InvalidInputError(
            f"leaf_mapping 必须为 Mapping，实际类型为 {type(leaf_mapping).__name__}"
        )

    gene_leaves = gene_tree.leaves
    raw_mapping: dict[str, str] = {}
    for key, value in leaf_mapping.items():
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
            raise LeafMappingError(
                f"映射目标 {s!r} 是物种内节点，基因叶只能映射到物种叶（拒绝非叶映射）",
                gene_leaf=g,
            )
        raw_mapping[g] = s

    missing = sorted(gene_leaves - raw_mapping.keys())
    if missing:
        raise LeafMappingError(f"存在未给出映射的基因叶：{missing}", gene_leaf=missing[0])
    return raw_mapping


def _gene_postorder(gene_tree: Tree) -> list[str]:
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


def _gene_preorder(gene_tree: Tree) -> list[str]:
    order: list[str] = []
    stack = [gene_tree.root]
    while stack:
        node = stack.pop()
        children = gene_tree.children.get(node)
        if children is None:
            continue
        order.append(node)
        stack.append(children[1])
        stack.append(children[0])
    return order


def _species_preorder(species_tree: Tree) -> list[str]:
    order: list[str] = []
    stack = [species_tree.root]
    while stack:
        node = stack.pop()
        order.append(node)
        children = species_tree.children.get(node)
        if children is not None:
            stack.append(children[1])
            stack.append(children[0])
    return order


def _species_postorder(species_tree: Tree) -> list[str]:
    order: list[str] = []
    stack: list[tuple[str, bool]] = [(species_tree.root, False)]
    while stack:
        node, expanded = stack.pop()
        children = species_tree.children.get(node)
        if children is None or expanded:
            order.append(node)
            continue
        stack.append((node, True))
        stack.append((children[1], False))
        stack.append((children[0], False))
    return order


def _lca_mapping(species_tree: Tree, gene_tree: Tree, raw_mapping: Mapping[str, str]) -> dict[str, str]:
    mapping = dict(raw_mapping)
    for node in _gene_postorder(gene_tree):
        c0, c1 = gene_tree.children[node]
        mapping[node] = _lca(species_tree, mapping[c0], mapping[c1])
    return mapping


def _build_reconciliation(
    species_tree: Tree,
    gene_tree: Tree,
    mapping: Mapping[str, str],
) -> ReconciliationResult:
    g_children = gene_tree.children
    postorder = _gene_postorder(gene_tree)
    events: list[Event] = []
    event_type: dict[str, str] = {}
    for node in postorder:
        c0, c1 = g_children[node]
        mapped = mapping[node]
        m0, m1 = mapping[c0], mapping[c1]
        etype = DUPLICATION if (m0 == mapped or m1 == mapped) else SPECIATION
        event_type[node] = etype
        events.append(
            Event(
                gene_node=node,
                event_type=etype,
                mapped_species=mapped,
                child_genes=(c0, c1),
                child_maps=(m0, m1),
            )
        )

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
            upwards: deque[str] = deque()
            cur = m_child
            while cur != m_parent:
                upwards.appendleft(cur)
                cur = s_parent[cur]
            path = [m_parent, *upwards]
            d = len(path) - 1
            lost: list[LostLineage] = []
            for i in range(d):
                p_i, p_next = path[i], path[i + 1]
                a, b = s_children[p_i]
                sibling = b if p_next == a else a
                if not parent_is_dup and i == 0:
                    continue
                lineage = LostLineage(
                    species_node=sibling,
                    descendant_leaves=species_tree.leaf_descendants[sibling],
                )
                lost.append(lineage)
                lineage_counts[sibling] = lineage_counts.get(sibling, 0) + 1
            loss_count = d if parent_is_dup else d - 1
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

    return ReconciliationResult(
        mapping=dict(sorted(mapping.items())),
        events=tuple(events),
        edge_losses=tuple(edge_losses),
        losses_by_lineage=dict(sorted(lineage_counts.items())),
        duplications=sum(1 for event in events if event.event_type == DUPLICATION),
        speciations=sum(1 for event in events if event.event_type == SPECIATION),
        total_losses=sum(edge.loss_count for edge in edge_losses),
    )


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
    raw_mapping = _validate_leaf_mapping(species_tree, gene_tree, leaf_mapping)

    mapping = _lca_mapping(species_tree, gene_tree, raw_mapping)
    return _build_reconciliation(species_tree, gene_tree, mapping)


def _validate_reconciliation_inputs(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_mapping: Mapping[str, str],
    max_nodes: int,
) -> dict[str, str]:
    _check_max_nodes(max_nodes)
    if not isinstance(species_tree, Tree):
        raise InvalidInputError(
            f"species_tree 必须为 build_tree 产出的 Tree，实际为 {type(species_tree).__name__}"
        )
    if not isinstance(gene_tree, Tree):
        raise InvalidInputError(
            f"gene_tree 必须为 build_tree 产出的 Tree，实际为 {type(gene_tree).__name__}"
        )
    return _validate_leaf_mapping(species_tree, gene_tree, leaf_mapping)


def _check_budget(budget: object) -> int | None:
    if budget is None:
        return None
    if isinstance(budget, bool) or not isinstance(budget, int):
        raise InvalidInputError(
            f"budget 必须为正整数 int 或 None，实际类型为 {type(budget).__name__}"
        )
    if budget < 1:
        raise InvalidInputError(f"budget 必须为正整数，实际为 {budget}")
    return budget


def _prune_local_states(
    states: dict[Count, tuple[str, ...]],
    species_order: Mapping[str, int],
) -> dict[Count, tuple[str, ...]]:
    kept: dict[Count, tuple[str, ...]] = {}
    for count, witness in states.items():
        duplicates, losses = count
        if any(
            other_duplications <= duplicates
            and other_losses <= losses
            and (other_duplications, other_losses) != count
            for other_duplications, other_losses in states
        ):
            continue
        old = kept.get(count)
        if old is None or tuple(species_order[node] for node in witness) < tuple(
            species_order[node] for node in old
        ):
            kept[count] = witness
    return kept


def reconcile_pareto_frontier(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_mapping: Mapping[str, str],
    *,
    max_nodes: int = DEFAULT_MAX_NODES,
    budget: int | None = None,
) -> FrontierResult:
    """枚举所有合法内部映射的 Pareto 计数前沿。

    树 DP 在每个基因节点和映射物种节点下保存非支配
    ``(duplications, losses)``；不通过反复调用固定权重搜索前沿。
    """
    raw_mapping = _validate_reconciliation_inputs(
        species_tree, gene_tree, leaf_mapping, max_nodes
    )
    remaining_budget = _check_budget(budget)
    canonical_species_order = {
        node: i for i, node in enumerate(_species_postorder(species_tree))
    }
    postorder = _gene_postorder(gene_tree)
    all_descendants: dict[str, set[str]] = {}
    for node in _species_preorder(species_tree).__reversed__():
        all_descendants[node] = {node}
        children = species_tree.children.get(node)
        if children is not None:
            all_descendants[node].update(all_descendants[children[0]])
            all_descendants[node].update(all_descendants[children[1]])

    states_by_gene: dict[str, dict[str, dict[Count, tuple[str, ...]]]] = {}
    complete = True
    for gene in gene_tree.leaves:
        states_by_gene[gene] = {raw_mapping[gene]: {(0, 0): ()}}

    def consume(cost_units: int = 1) -> bool:
        nonlocal remaining_budget
        if remaining_budget is None:
            return True
        remaining_budget -= cost_units
        return remaining_budget >= 0

    for gene in postorder:
        left_gene, right_gene = gene_tree.children[gene]
        left_states = states_by_gene[left_gene]
        right_states = states_by_gene[right_gene]
        local_by_species: dict[str, dict[Count, tuple[str, ...]]] = {}

        for mapped in _species_preorder(species_tree):
            candidates: dict[Count, tuple[str, ...]] = {}
            mapped_descendants = all_descendants[mapped]
            is_speciation_parent = mapped in species_tree.children
            if is_speciation_parent:
                child_a, child_b = species_tree.children[mapped]
                side_a = all_descendants[child_a]
                side_b = all_descendants[child_b]

            for left_map, left_counts in left_states.items():
                if left_map not in mapped_descendants:
                    continue
                for right_map, right_counts in right_states.items():
                    if right_map not in mapped_descendants:
                        continue
                    speciation = (
                        is_speciation_parent
                        and left_map != mapped
                        and right_map != mapped
                        and (
                            (left_map in side_a and right_map in side_b)
                            or (left_map in side_b and right_map in side_a)
                        )
                    )
                    for (left_count, left_witness) in left_counts.items():
                        for (right_count, right_witness) in right_counts.items():
                            if not consume():
                                complete = False
                                break
                            left_edge_losses = species_tree.depth[left_map] - species_tree.depth[mapped]
                            right_edge_losses = species_tree.depth[right_map] - species_tree.depth[mapped]
                            if speciation:
                                left_edge_losses -= 1
                                right_edge_losses -= 1
                            count = (
                                left_count[0] + right_count[0] + (0 if speciation else 1),
                                left_count[1]
                                + right_count[1]
                                + left_edge_losses
                                + right_edge_losses,
                            )
                            witness = (*left_witness, *right_witness, mapped)
                            old = candidates.get(count)
                            if old is None or tuple(canonical_species_order[node] for node in witness) < tuple(
                                canonical_species_order[node] for node in old
                            ):
                                candidates[count] = witness
                        if not complete:
                            break
                    if not complete:
                        break
                if not complete:
                    break

            if not complete:
                break
            pruned = _prune_local_states(candidates, canonical_species_order)
            if pruned:
                local_by_species[mapped] = pruned

        if not complete:
            return FrontierResult(
                solutions=(),
                complete=False,
                budget_remaining=remaining_budget,
            )
        states_by_gene[gene] = local_by_species

    root_candidates: dict[Count, tuple[str, ...]] = {}
    for mapped_counts in states_by_gene[gene_tree.root].values():
        root_candidates.update(mapped_counts)
    root_frontier = _prune_local_states(root_candidates, canonical_species_order)

    solutions: list[FrontierSolution] = []
    for count, witness in sorted(root_frontier.items()):
        mapping = dict(raw_mapping)
        mapping.update(zip(postorder, witness))
        solutions.append(
            FrontierSolution(
                duplications=count[0],
                losses=count[1],
                reconciliation=_build_reconciliation(species_tree, gene_tree, mapping),
            )
        )
    return FrontierResult(
        solutions=tuple(solutions),
        complete=True,
        budget_remaining=remaining_budget,
    )


def _check_positive_weight(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise InvalidInputError(
            f"{name} 必须为正整数 int，实际类型为 {type(value).__name__}"
        )
    if value < 1:
        raise InvalidInputError(f"{name} 必须为正整数，实际为 {value}")
    return value


def reconcile_weighted(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_mapping: Mapping[str, str],
    duplication_weight: int,
    loss_weight: int,
    *,
    max_nodes: int = DEFAULT_MAX_NODES,
    budget: int | None = None,
) -> WeightedReconciliationResult:
    """从 Pareto 前沿选择正整数权重下的全部并列最优方案。"""
    duplication_weight = _check_positive_weight(duplication_weight, "duplication_weight")
    loss_weight = _check_positive_weight(loss_weight, "loss_weight")
    frontier = reconcile_pareto_frontier(
        species_tree,
        gene_tree,
        leaf_mapping,
        max_nodes=max_nodes,
        budget=budget,
    )
    if not frontier.complete:
        raise InvalidInputError(
            "搜索预算已耗尽，无法保证加权最优；请增大 budget 或使用 budget=None"
        )
    optimal_cost = min(
        solution.duplications * duplication_weight + solution.losses * loss_weight
        for solution in frontier.solutions
    )
    optimal_solutions = tuple(
        solution
        for solution in frontier.solutions
        if solution.duplications * duplication_weight
        + solution.losses * loss_weight
        == optimal_cost
    )
    return WeightedReconciliationResult(
        duplication_weight=duplication_weight,
        loss_weight=loss_weight,
        optimal_cost=optimal_cost,
        solutions=optimal_solutions,
    )
