"""基因树 / 物种树的 LCA 事件归因（reconciliation）。

本模块只实现题面契约规定的数学子集：

* 输入为有根、满二叉（内部节点恰有两个子节点）的物种树与基因树；
* 基因叶子映射到物种叶子，允许一个物种对应多个基因；
* 基因内部节点映射为其子节点映射在物种树上的最近公共祖先（LCA）；
* 任一子节点的映射等于父节点映射，该内部节点记为 duplication，否则记为 speciation；
* 每条基因树父子边 (g, c) 的丢失数为物种树上 map(g)→map(c) 的路径长度减 1；
  当 g 是 duplication 时不减 1（即按完整路径长度计）。
  丢失的具体支系取路径上每个“未被选择的兄弟子树”，以其根节点标签标识。

模块为纯标准库、无浮点运算、无可变状态：``reconcile`` 不修改任何输入，
要么返回完整结果，要么抛异常，不存在部分变更。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

#: 接口边界：每棵树允许的最大节点数。
MAX_NODES = 200

#: 物种分化事件。
EVENT_SPECIATION = "S"
#: 复制事件。
EVENT_DUPLICATION = "D"


class ReconciliationError(ValueError):
    """本模块所有错误的基类（消息均可定位到具体标签/参数）。"""


class TreeError(ReconciliationError):
    """树结构非法：非二叉、重名、成环/多父、超出节点上限等。"""


class MappingError(ReconciliationError):
    """基因叶→物种叶映射非法：漏配、多配、目标不存在或不是物种叶。"""


@dataclass(eq=False)
class Node:
    """有根二叉树节点。

    ``left`` 与 ``right`` 同时为 None 表示叶子，同时为 Node 表示内部节点，
    只给一个子节点会在 :class:`Tree` 校验时被拒绝。节点按对象身份区分，
    因此同一个 Node 不能在树中出现两次（多父/成环均会被拒绝）。
    """

    label: str
    left: Optional["Node"] = None
    right: Optional["Node"] = None

    @property
    def children(self) -> tuple["Node", ...]:
        """叶子返回空元组，内部节点返回 ``(left, right)``。"""
        if self.left is None and self.right is None:
            return ()
        return (self.left, self.right)  # type: ignore[return-value]

    @property
    def is_leaf(self) -> bool:
        return self.left is None and self.right is None


def node_from_spec(spec) -> Node:
    """由确定性的嵌套元组构造 :class:`Node`。

    * 叶子：``"a"`` 或 ``("a",)``；
    * 内部节点：``("label", 左子树spec, 右子树spec)``。

    解析阶段不做完整校验（成环由元组结构天然不可能），完整校验在 :class:`Tree`。
    """
    if isinstance(spec, str):
        return Node(spec)
    if isinstance(spec, tuple):
        if len(spec) == 1 and isinstance(spec[0], str):
            return Node(spec[0])
        if len(spec) == 3 and isinstance(spec[0], str):
            return Node(spec[0], node_from_spec(spec[1]), node_from_spec(spec[2]))
    raise TreeError(
        "非法的树 spec：叶子须为 str，内部节点须为 (label, left, right)；"
        f"实际得到 {spec!r}"
    )


def newick(tree: "Tree | Node") -> str:
    """以最简 Newick 风格（无枝长）渲染树，仅用于打印/调试。"""

    root = tree.root if isinstance(tree, Tree) else tree

    def walk(n: Node) -> str:
        if n.is_leaf:
            return n.label
        return f"({walk(n.left)},{walk(n.right)}){n.label}"  # type: ignore[arg-type]

    return walk(root) + ";"


class Tree:
    """经过完整校验的有根满二叉树。

    校验内容：根为 Node；标签为非空 str 且全树唯一；每个节点要么是叶子要么
    恰有两个 Node 子节点；无多父、无环；节点数不超过 ``max_nodes``（默认
    :data:`MAX_NODES`）。校验失败抛 :class:`TreeError`，对象不会半成品对外暴露
    （构造失败即不存在）。
    """

    def __init__(self, root: Node, *, max_nodes: int = MAX_NODES, name: str = "tree"):
        if isinstance(max_nodes, bool) or not isinstance(max_nodes, int):
            # 本模块不接受任何浮点数值；bool 是 int 的子类，必须显式拒绝。
            raise TreeError(
                f"max_nodes 必须是 int 且不能是 bool，实际得到 {type(max_nodes).__name__}"
            )
        if not 1 <= max_nodes <= MAX_NODES:
            raise TreeError(
                f"max_nodes 必须落在 [1, {MAX_NODES}]，实际得到 {max_nodes}"
            )
        if not isinstance(name, str) or not name:
            raise TreeError("name 必须是非空 str")
        if not isinstance(root, Node):
            raise TreeError(f"{name} 的根必须是 Node，实际得到 {type(root).__name__}")

        self.max_nodes = max_nodes
        self.name = name
        self.root = root
        self.parent: dict[Node, Node] = {}
        self.depth: dict[Node, int] = {}
        self.label_to_node: dict[str, Node] = {}
        self.postorder: list[Node] = []
        self.leaves: list[Node] = []

        self._validate()

    def _validate(self) -> None:
        white, gray, black = 0, 1, 2
        color: dict[Node, int] = {self.root: white}

        def visit(n: Node, parent: Optional[Node], depth: int) -> None:
            state = color.get(n, white)
            if state == gray:
                raise TreeError(
                    f"{self.name} 存在环：节点 {n.label!r} 再次出现在自己的祖先链上"
                )
            if state == black:
                # 同一节点被第二次访问：DAG（多父），作为树同样拒绝。
                prev = self.parent.get(n)
                prev_label = prev.label if prev is not None else "<root>"
                par_label = parent.label if parent is not None else "<root>"
                raise TreeError(
                    f"{self.name} 节点 {n.label!r} 有多个父节点 "
                    f"({prev_label!r} 与 {par_label!r})，树不得共享子树"
                )
            color[n] = gray

            if not isinstance(n.label, str) or not n.label:
                raise TreeError(
                    f"{self.name} 节点标签必须是非空 str，实际得到 {n.label!r}"
                )
            if n.label in self.label_to_node:
                raise TreeError(f"{self.name} 节点标签重复：{n.label!r}")
            self.label_to_node[n.label] = n
            self.depth[n] = depth
            if parent is not None:
                self.parent[n] = parent

            kids = n.left, n.right
            if kids == (None, None):
                self.leaves.append(n)
            elif None in kids:
                raise TreeError(
                    f"{self.name} 节点 {n.label!r} 不是满二叉：恰有一个子节点"
                )
            else:
                left, right = kids
                if not isinstance(left, Node) or not isinstance(right, Node):
                    raise TreeError(
                        f"{self.name} 节点 {n.label!r} 的子节点必须是 Node"
                    )
                if left is right:
                    raise TreeError(
                        f"{self.name} 节点 {n.label!r} 的两个子节点是同一个对象"
                    )
                visit(left, n, depth + 1)
                visit(right, n, depth + 1)

            color[n] = black
            self.postorder.append(n)

        visit(self.root, None, 0)
        if len(self.postorder) > self.max_nodes:
            raise TreeError(
                f"{self.name} 节点数 {len(self.postorder)} 超过上限 {self.max_nodes}"
            )

    def lca(self, a: Node, b: Node) -> Node:
        """物种树上两节点的最近公共祖先。"""
        while self.depth[a] > self.depth[b]:
            a = self.parent[a]
        while self.depth[b] > self.depth[a]:
            b = self.parent[b]
        while a is not b:
            a = self.parent[a]
            b = self.parent[b]
        return a

    def path(self, ancestor: Node, descendant: Node) -> list[Node]:
        """物种树上 ancestor 到 descendant 的自上而下节点序列（含两端）。"""
        if self.depth[ancestor] > self.depth[descendant]:
            raise TreeError("path() 的祖先比后代更深，不可能在有根树上成立")
        chain = [descendant]
        cur = descendant
        while cur is not ancestor:
            cur = self.parent.get(cur)  # type: ignore[assignment]
            if cur is None:
                raise TreeError("path() 的两个节点不在同一条祖先链上")
            chain.append(cur)
        chain.reverse()
        return chain


@dataclass(frozen=True)
class Event:
    """一个基因内部节点上的归因事件。

    ``children`` 为 ``(基因子节点标签, 该子节点映射到的物种标签)``，保持基因树
    左、右子顺序；``edge_losses`` 为每条父子边上标出的具体缺失支系（物种兄弟
    子树的根标签），自上而下排列。
    """

    gene: str
    species: str
    kind: str  # EVENT_SPECIATION / EVENT_DUPLICATION
    children: tuple[tuple[str, str], tuple[str, str]]
    edge_losses: tuple[tuple[str, tuple[str, ...]], tuple[str, tuple[str, ...]]]


@dataclass(frozen=True)
class Reconciliation:
    """归因结果（全部以标签为键，可直接序列化/打印）。"""

    mapping: dict[str, str]
    events: tuple[Event, ...]
    edge_losses: dict[tuple[str, str], tuple[str, ...]]
    edge_loss_counts: dict[tuple[str, str], int]
    n_speciations: int
    n_duplications: int
    n_losses: int


def reconcile(
    species_tree: Tree,
    gene_tree: Tree,
    leaf_map: dict[str, str],
) -> Reconciliation:
    """对一对物种树/基因树执行 LCA 归因。

    ``leaf_map`` 将每个基因叶子标签映射到一个**物种叶子**标签：必须恰好覆盖
    全部基因叶（不多不少），目标必须存在且为物种叶，否则抛
    :class:`MappingError`。树结构问题抛 :class:`TreeError`。
    """
    if not isinstance(species_tree, Tree):
        raise TypeError("species_tree 必须是 Tree 实例")
    if not isinstance(gene_tree, Tree):
        raise TypeError("gene_tree 必须是 Tree 实例")
    if not isinstance(leaf_map, dict):
        raise MappingError(
            f"leaf_map 必须是 dict，实际得到 {type(leaf_map).__name__}"
        )

    for key, value in leaf_map.items():
        if not isinstance(key, str):
            raise MappingError(f"leaf_map 的键必须是基因叶标签 str，实际得到 {key!r}")
        if not isinstance(value, str):
            raise MappingError(
                f"leaf_map[{key!r}] 的值必须是物种叶标签 str，实际得到 {value!r}"
            )

    species_leaves = {n.label for n in species_tree.leaves}
    gene_leaves = {n.label for n in gene_tree.leaves}

    missing = sorted(gene_leaves - leaf_map.keys())
    if missing:
        raise MappingError(f"基因叶缺少映射：{missing}")
    extra = sorted(leaf_map.keys() - gene_leaves)
    if extra:
        raise MappingError(f"leaf_map 含有不是基因叶的多余键：{extra}")

    for gene_label in sorted(gene_leaves):
        target = leaf_map[gene_label]
        target_node = species_tree.label_to_node.get(target)
        if target_node is None:
            raise MappingError(
                f"基因叶 {gene_label!r} 映射到不存在的物种标签 {target!r}"
            )
        if not target_node.is_leaf:
            raise MappingError(
                f"基因叶 {gene_label!r} 映射到物种内部节点 {target!r}："
                "只允许映射到物种叶子"
            )

    # 1) LCA 映射（后序，叶子先行）。
    gmap: dict[Node, Node] = {}
    for gn in gene_tree.postorder:
        if gn.is_leaf:
            gmap[gn] = species_tree.label_to_node[leaf_map[gn.label]]
        else:
            gmap[gn] = species_tree.lca(gmap[gn.left], gmap[gn.right])  # type: ignore[index]

    # 2) 事件判定与逐边丢失。
    events: list[Event] = []
    edge_losses: dict[tuple[str, str], tuple[str, ...]] = {}
    edge_loss_counts: dict[tuple[str, str], int] = {}
    n_dup = n_spec = n_loss = 0

    for gn in gene_tree.postorder:
        if gn.is_leaf:
            continue
        mapped = gmap[gn]
        child_events = []
        child_losses = []
        is_dup = any(gmap[kid] is mapped for kid in (gn.left, gn.right))
        kind = EVENT_DUPLICATION if is_dup else EVENT_SPECIATION
        if is_dup:
            n_dup += 1
        else:
            n_spec += 1

        for kid in (gn.left, gn.right):
            path = species_tree.path(mapped, gmap[kid])  # type: ignore[arg-type]
            k = len(path) - 1
            # S：路径长度 k 减 1，跳过 v0 处的兄弟（另一基因子占据该分支）；
            # D：不减 1，v0 处未被选择的兄弟也是丢失。
            start = 0 if is_dup else 1
            lost = []
            for i in range(start, k):
                vi, nxt = path[i], path[i + 1]
                sibling = vi.right if vi.left is nxt else vi.left
                lost.append(sibling.label)  # type: ignore[union-attr]
            lost_t = tuple(lost)
            edge_key = (gn.label, kid.label)
            edge_losses[edge_key] = lost_t
            edge_loss_counts[edge_key] = len(lost_t)
            n_loss += len(lost_t)
            child_events.append((kid.label, gmap[kid].label))  # type: ignore[union-attr]
            child_losses.append((kid.label, lost_t))

        events.append(
            Event(
                gene=gn.label,
                species=mapped.label,
                kind=kind,
                children=(child_events[0], child_events[1]),
                edge_losses=(child_losses[0], child_losses[1]),
            )
        )

    return Reconciliation(
        mapping={gn.label: gmap[gn].label for gn in gene_tree.postorder},
        events=tuple(events),
        edge_losses=edge_losses,
        edge_loss_counts=edge_loss_counts,
        n_speciations=n_spec,
        n_duplications=n_dup,
        n_losses=n_loss,
    )
