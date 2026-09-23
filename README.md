# 基因树 / 物种树事件归因（LCA reconciliation）

离线、纯标准库的可复用 Python 库：给定有根二叉**物种树**与**基因树**，
以及基因叶→物种叶的映射，计算每个基因节点的 LCA 映射、事件类型
（物种分化 S / 复制 D）与每条基因树父子边上的丢失数及具体缺失支系。

## 运行环境

- Python **3.14.7**（3.10+ 亦可运行；开发与验收在 3.14.7）
- Linux 原生运行；**仅用标准库，无需安装任何第三方依赖**，无网络、账号、
  密钥、数据库、中间件、云 API 或 Docker。

## 准确命令

```bash
# 必要测试（含手工期望值与独立参考交叉核对）
python -m unittest discover -s tests -v

# 固定输入演示（正常结果 + 一个实际触发的拒绝边界，秒级完成，远低于 8 秒）
python demo.py
```

## 快速上手

```python
from reconciliation import Tree, node_from_spec, reconcile

# 叶子写字符串；内部节点写 (标签, 左子树, 右子树)
species = Tree(node_from_spec(("R", ("Y", "A", "B"), "C")))
gene    = Tree(node_from_spec(("g0", "a1", "c")))
result  = reconcile(species, gene, {"a1": "A", "c": "C"})

result.mapping            # {'a1': 'A', 'c': 'C', 'g0': 'R'}
result.n_losses           # 1
result.edge_losses        # {('g0', 'a1'): ('B',), ('g0', 'c'): ()}
```

## 公开接口与数据结构

均定义于 `reconciliation.py`，以**标签（str）**为公开键，不暴露内部节点对象。

- `Node(label, left=None, right=None)`：树节点；两个子节点同为 `None` 为叶，
  同为 `Node` 为内部节点。节点按对象身份区分，禁止共享子树。
- `node_from_spec(spec)`：由嵌套元组确定性构造。叶 `"A"` / `("A",)`；
  内部节点 `("R", left_spec, right_spec)`。
- `Tree(root, *, max_nodes=200, name="tree")`：构造时完成全部校验
  （见“输入上限与错误语义”），并提供 `lca` / `path` / `leaves` /
  `postorder` / `label_to_node` / `depth` / `parent`。
- `newick(tree)`：最简 Newick 文本（无枝长），仅用于打印。
- `reconcile(species_tree, gene_tree, leaf_map) -> Reconciliation`：纯函数，
  不修改任何输入，要么返回完整结果要么抛异常，无部分结果。
- `Reconciliation`（frozen dataclass）：
  - `mapping: dict[str, str]`：基因节点标签 → 物种节点标签（全节点 LCA 映射）；
  - `events: tuple[Event, ...]`：每个基因内部节点一条事件；
  - `edge_losses: dict[(父, 子), tuple[str, ...]]`：每条基因树父子边上标出的
    **具体缺失支系**，以丢失发生处物种兄弟子树的根标签标识，自上而下排列；
  - `edge_loss_counts`：同一批边的丢失计数（= 标签元组长度）；
  - `n_speciations / n_duplications / n_losses`：合计。
- `Event` 字段：`gene, species, kind, children, edge_losses`；
  `kind` 为 `"S"`（`EVENT_SPECIATION`）或 `"D"`（`EVENT_DUPLICATION`）；
  `children` / `edge_losses` 保持基因树左、右子顺序。

## 算法契约（严格按题面）

1. 基因叶子直接映到 `leaf_map` 指定的物种叶子；
2. 基因内部节点 `g` 映到两个子映射在物种树上的最近公共祖先（LCA）；
3. **任一**子节点映射等于 `g` 的映射 ⇒ `g` 为 duplication（D），否则为
   speciation（S）；
4. 对每条基因树父子边 `(g, c)`，设物种树上 `map(g) → map(c)` 的路径长度为
   `k`（边数）：
   - `g` 为 S：丢失数 = `k - 1`；
   - `g` 为 D：丢失数 = `k`（不减 1）；
   - 具体缺失支系：沿路径在每个处分叉点取“未被选择的兄弟子树”，
     以该兄弟子树根标签标识；S 时跳过路径起点处的兄弟（另一基因子占据该
     物种分支，不是丢失），D 时不跳过（被复制的另一拷贝映射在原处，该兄弟
     在本边上确为丢失）。

树遍历为迭代式后序，LCA 用深度爬升，时间复杂度 O((Nₛ+N₉)·树高)，
空间 O(Nₛ+N₉)。

## 输入上限

- 每棵树节点数 **≤ 200**（`Tree(..., max_nodes=...)` 可收紧，不可放宽）；
- 满二叉：内部节点恰有两个子节点，单叶树（1 节点）允许；
- 标签必须为非空 `str` 且全树唯一；
- `leaf_map` 必须恰好覆盖全部基因叶子（不多不少）。

## 错误语义（均可定位，不静默纠正）

- `TreeError`：根不是 `Node`、标签空/重名、单子（非满二叉）、子节点类型错、
  同一子对象被两个父节点引用（DAG）、成环（沿祖先链回到灰色节点）、
  节点数超上限、`max_nodes` 非 `int`/是 `bool`/越界；
- `MappingError`：`leaf_map` 不是 dict、键值不是 str、基因叶漏配、
  存在非基因叶的多余键、目标物种标签不存在、**目标是物种内部节点
  （拒绝非叶映射）**；
- 树参数类型错误抛 `TypeError`。
- 本库无浮点输入与浮点结果，故不存在 NaN/Infinity 问题；唯一数值参数
  `max_nodes` 拒绝 `float`（含 NaN/Infinity）与 `bool`（`bool` 虽是
  `int` 子类仍显式拒绝）。
- `reconcile` 为纯函数，不存在批次/有状态部分变更；本算法是确定性 LCA
  映射，不存在搜索过程，因此没有“预算耗尽”这一结果类别——输入合法即得
  唯一解，不合法即上述异常。

## 浮点与误差口径

无任何浮点运算：映射、事件、计数全部基于整数树拓扑。**结果精确，无误差、
无容差**。

## 测试与验收要点

- `tests/test_reconciliation.py` 的期望值全部**手工推导**（用例注释给出
  LCA 与每条边的手算过程），覆盖：纯物种分化零丢失、一物种多基因复制、
  单支丢失、整支姐妹群丢失、D 边“不减 1”与 S 边“减 1”的差异用例、
  多次复制、组合情形、左右子互换与叶名置换不变性、输入不可变、
  结果冻结，以及全部拒绝边界（含手工构造的环与多父 DAG、199 节点接受 /
  201 节点拒绝）。
- 文件内另有**独立参考实现** `_reference_reconcile`：直接在嵌套元组上用
  “祖先集合相交”求 LCA（不同于被测的深度爬升法），并按题面公式独立计数，
  仅作交叉核对的**补充**手段；它不 import、不调用被测核心的任何函数来
  产生期望值。

## 设计取舍

- 丢失支系用“分叉点处未选择的兄弟子树**根标签**”标识：这是路径上唯一确定
  且标签稳定的表示；一次标注代表该兄弟整支在本边上缺失（其规模可在物种树
  上展开，模块不重复罗列其全部叶子）。
- 仅实现题面契约的二叉 LCA 归因：不做多叉树、不做 HGT/ILS、不做代价最优
  调解搜索、不做 Newick 解析器（`newick` 只单向输出用于演示）。
- `Tree` 构造即校验并缓存 `parent/depth/后序序列`，使后续 LCA、路径与
  归因均为简单查询；节点按对象身份而非标签区分（标签仅用于对外报告）。

## 已知限制

- 不是行业标准 Newick/NHX 工具，不解析枝长、自举值、基因名分隔符；
- 不支持非二叉（多叉）树与无根树；
- 不区分“同一复制事件在两条子边上的丢失归属”之外的历史解释，输出严格对应
  题面定义的 LCA（最简约）归因；
- 节点上限 200 是接口边界而非性能上限。

## 文件

- `reconciliation.py`：核心库；
- `tests/test_reconciliation.py`：必要测试（手工期望值 + 独立参考）；
- `demo.py`：固定输入演示，真实计算并打印可复验事件路径，以及一个实际
  触发的 `MappingError` 拒绝边界。
