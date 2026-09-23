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
from dataclasses import dataclass
from typing import Mapping, Sequence

__all__ = [
    "Tree",
    "Event",
    "LostLineage",
    "EdgeLoss",
    "ReconciliationResult",
    "build_tree",
    "reconcile",
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
