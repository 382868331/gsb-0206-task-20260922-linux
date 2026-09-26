# 基因树与物种树的事件归因（LCA 复演）

离线、零第三方依赖的可复用 Python 库：输入有根二叉**物种树**与**基因树**及
基因叶→物种叶映射，输出 LCA 映射、事件表（speciation / duplication）、
逐基因父子边的丢失计数与**具体缺失支系**。

## 运行环境

- Python **3.14.7**（开发与验证版本；仅用标准库，`from __future__ import annotations`，
  3.10+ 亦可导入，但只承诺 3.14.7）。
- Linux 原生运行，装好 Python 后**无需安装任何其他依赖**（`requirements.txt` 为空占位）。
- 完全离线：无网络、账号、密钥、数据库、中间件、云 API 或 Docker。

## 快速命令

```bash
# 推荐：用项目目录内隔离 venv；不访问网络、不安装第三方包
python3.14 -m venv .venv
. .venv/bin/activate
python --version  # 应为 Python 3.14.7
```

```bash
# 必要测试（30 个用例，含基线参考与合法映射穷举交叉核对）
python -m unittest discover -s tests -v

# 固定输入演示：一个正常复演结果 + 两个实际触发的拒绝边界
python demo.py

# 增量演示：计数前沿 DP、正整数权重选择、预算耗尽与小树独立枚举
python demo_iteration.py
```

两个演示均真实计算、无 `sleep`、无预录结果、无账号/密钥/网络依赖，实测均远低于 90 秒。

## 公开接口

```python
from reconciliation import build_tree, reconcile

sp  = build_tree({"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}, "r",
                 name="species")
gn  = build_tree({"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")},
                 "gr", name="gene")
res = reconcile(sp, gn, {"a1": "A", "b1": "B", "a2": "A", "c1": "C"})
```

不预先固定复制/丢失权重时使用：

```python
from reconciliation import reconcile_pareto_frontier, reconcile_weighted

frontier = reconcile_pareto_frontier(sp, gn, leaf_mapping)
for item in frontier.solutions:
    item.duplications, item.losses
    item.reconciliation.mapping
    item.reconciliation.events

best = reconcile_weighted(sp, gn, leaf_mapping, duplication_weight=2, loss_weight=1)
```

### `build_tree(children, root=None, *, name="tree", max_nodes=200) -> Tree`

- `children`: `{内节点: (左子, 右子)}`，映射类型（dict）；两个子节点必须是
  互不相同、且不等于父节点的非空 `str`。
- 叶子 = 作为子节点出现、但自身不在键中的节点。
- `root`: 可省略；省略时必须能唯一推断（恰有一个不作为任何子节点出现的键），
  否则抛 `TreeStructureError`。`children={}` 时必须显式给 `root`（单叶平凡树）。
- 返回不可变语义的 `Tree`（`.root/.children/.parent/.depth/.leaves/.internal/
  .leaf_descendants/.nodes`）。

### `reconcile(species_tree, gene_tree, leaf_mapping, *, max_nodes=200) -> ReconciliationResult`

- `leaf_mapping`: `{基因叶: 物种叶}`，**每个基因叶恰好一条**；同一物种叶可对应
  多个基因叶（一物种多基因）。目标必须是物种**叶**，非叶映射一律拒绝。
- `ReconciliationResult` 字段：
  - `mapping: dict[str,str]` — 全部基因节点（叶+内节点）→ 物种节点；
  - `events: tuple[Event, ...]` — 事件表（基因树后序），`Event` 含
    `gene_node / event_type / mapped_species / child_genes / child_maps`；
  - `edge_losses: tuple[EdgeLoss, ...]` — 逐边计数，`EdgeLoss` 含
    基因边两端、物种映射两端、`path_length`、`loss_count`、
    `lost_lineages`（`LostLineage(species_node, descendant_leaves)`，自上而下有序）；
  - `losses_by_lineage: dict[str,int]` — 每个缺失物种支系被计丢的次数；
  - `duplications / speciations / total_losses` — 累计计数。
- `reconcile_pareto_frontier(..., budget=None) -> FrontierResult`：树 DP 返回全部
  非支配 `(duplications, losses)`；每个计数只保留一份按规范映射顺序最小的见证。
  - `FrontierResult.solutions: tuple[FrontierSolution, ...]`，按计数排序；
  - `FrontierSolution.duplications/losses/reconciliation`；
  - `complete=False` 且 `budget_remaining<0` 表示预算耗尽，前沿可能不全，不能当作“证明无解”。
- `reconcile_weighted(..., duplication_weight, loss_weight, budget=None)`：
  两个权重必须是正整数 `int`（拒绝 `bool`、0、负数、浮点），从已算前沿选出代价并列最小的方案；
  预算耗尽时拒绝给可能错误的加权最优。
- 常量：`SPECIATION="speciation"`、`DUPLICATION="duplication"`。

## 算法规约（与题面逐条对应）

1. **映射**：基因叶映射由输入给定；每个基因内节点映射到两个子映射在物种树中的
   LCA（深度抬升法求祖先）。
2. **事件**：任一子映射等于父映射 ⇒ `duplication`，否则 `speciation`。
3. **逐边丢失数**：对基因父子边 `(g,c)`，令 `d` 为物种树上 `map(g)→map(c)` 的
   路径边数（LCA 性质保证前者是后者祖先）：
   - `g` 为 speciation：丢失 `d-1`；
   - `g` 为 duplication：丢失 `d`（不减 1）。
4. **具体缺失支系**：沿物种树路径每个分叉点 `p_i→p_{i+1}`，不在路径上的兄弟
   子树即一条缺失支系（附带其后裔叶集合，便于核对“缺了谁”）；
   speciation 时跳过**起点**处的兄弟——该支系由同一基因节点的另一个子支占据，
   不是丢失。duplication 不跳过（两个子支都从同一物种节点长出，兄弟支系无人占据）。
   枚举到的支系数与计数公式逐边 `assert` 自检一致。
5. 复杂度：校验与复演均为线性级（LCA 按深度爬升，节点 ≤ 200）。

## 不预设权重的计数前沿

前沿 DP 自底向上处理基因树；对每个“基因节点 → 物种节点”保存一组
`(复制数, 丢失数)`。合并两个子树时枚举该基因节点上的一次合法 speciation 或
duplication，把子边丢失数加总，然后严格删除被其他状态逐项不大于、且至少一项更小的状态。
它不是循环替换多组权重；权重只在完整前沿得到后用于线性选择，因此并列方案会同时保留。

相同计数的等价路径只保留一份见证。规范映射顺序为：先按基因内节点后序
（左子树先于右子树，子先于父），再按物种树左右优先的后序排列；LCA 位于祖先前，
因而会稳定选择最深的合法 LCA 映射。

`budget` 按实际执行的子状态合并计数。预算耗尽时立即停止，并返回
`FrontierResult(solutions=(), complete=False, budget_remaining=-1)`；
调用方必须把它解释为“可能未探索完”，而不是空前沿或证明无解。

### 适用范围中的一个重要数学结论

本项目的输入叶物种归属固定，事件只包含 speciation/duplication/loss，
不增加水平转移。把任何基因内节点从其 LCA 抬到祖先，无法减少复制数；
其下行边还会增加或保持丢失数。故 LCA 映射同时最小化复制数和丢失数：
在本严格模型内，不会出现“一个方案多复制少丢失、另一个少复制多丢失”的非支配对。
小树测试独立枚举 3–4 叶物种树、2–4 个基因叶的全部形状和全部固定叶归属，
核对 DP 前沿并验证该结论。若需要这类真正的多目标权衡，模型必须放宽为转移、
不固定叶归属、多叉/不定根等额外事件集合；本库不把这些模型伪装成当前能力。

## 输入上限

- 单棵树节点数 **≤ 200**（题面接口边界），`max_nodes` 只能为 `1≤n≤200` 的 `int`。
- 满二叉树节点数恒为奇数：199 节点（100 叶）接受、201 节点拒绝，测试已覆盖。
- 标签为非空 `str`；左右子必须不同；不支持多边、非二叉、森林。

## 错误语义（可定位，不静默纠正）

所有错误继承 `ReconciliationError`（也是 `ValueError`）：

| 异常 | 触发 | 定位字段 |
|---|---|---|
| `InvalidInputError` | 类型/取值非法：非映射入参、标签非非空 str、`max_nodes` 为 bool/浮点/NaN/Infinity/越界等 | — |
| `TreeStructureError` | 自环、回边（循环）、多父、多根/根不可推断、根未知、不可达分量、非二叉、左右子相同 | `.tree`、`.node` |
| `LeafMappingError` | 映射缺失/多余键、键非基因叶、目标非物种节点、**非叶映射（映到物种内节点）**、键值类型非法 | `.gene_leaf`（键本身非法时为 `None`） |

- 数值输入：唯一数值参数 `max_nodes` 拒绝 `bool`（`True/False` 不当 1/0）、
  拒绝 `float`（含 `NaN`、`+/-Infinity`，信息中区分 NaN 与 Infinity）、
  拒绝 ≤0 与 >200。
- **原子性/无状态**：库不持有任何可变状态；`reconcile` 先完成全部校验，
  再一次性构造结果，非法输入时不会产生任何部分结果或外部变更。
- 固定 LCA 入口无回溯/搜索过程；增量前沿搜索接受显式正整数 `budget`，
  耗尽时以 `FrontierResult.complete=False` 标明“可能不全”，与非法输入和“证明无解”区分。

## 浮点与误差口径

本库是离散算法：**全部输入、中间值与输出都是整数或字符串，没有任何浮点运算**，
因此不存在浮点误差与容差——计数（路径长、丢失数、累计）是精确整数，
实现中以整数恒等断言自洽（逐边枚举支系数 == 公式计数；
`total_losses == Σ逐边 == Σ按支系`）。不会用宽容差掩盖边界错误。

## 设计取舍

- 树以 `{父: (子, 子)}` 邻接字典输入，标签即节点身份；构造时一次性算出
  父表、深度、后裔叶集，复演阶段只做查表与爬升。
- 循环检测用三色迭代 DFS（白/灰/黑），遇灰节点即报告回边并给出节点定位；
  不使用递归，避免深树触发递归上限（虽然 200 节点内不会发生）。
- 事件表与逐边结果按基因树后序排列（子先于父、先左后右），输出稳定、可 diff。
- “缺失支系”按物种树兄弟子树枚举，并给出后裔叶，使每条丢失都可人工追溯；
  同一支系在不同边上可被计多次（独立丢失事件），`losses_by_lineage` 给累计次数。
- 只实现本题契约（有根二叉、LCA 复演、D/S 两类事件），不宣称兼容
  完整行业标准（不做 HGT、转移-丢失联合模型、非二叉/无定根树和解）。

## 已知限制

- 仅支持有根二叉满树与 LCA（复制-丢失）复演；不支持水平基因转移、
  基因树不定根、多叉节点、分支长度、概率模型或转移-复制-丢失联合搜索。
- 前沿搜索适用于复制/丢失两类事件；固定叶归属且无转移时，标准 LCA 解是唯一计数最优点。
- 叶映射必须全给定，不做“最适合映射”推断。
- 单文件模块布局（`reconciliation.py`），未打包为可安装发行版；
  使用时把该文件放入项目即可导入。
- 左右子顺序只影响输出排列，不影响计数（镜像不变性有测试覆盖）。

## 文件清单

```
reconciliation.py        # 核心库（唯一交付模块，纯标准库）
demo.py                  # 固定输入演示
demo_iteration.py        # 前沿 DP、加权选择、预算与小树枚举演示
tests/test_reconciliation.py  # 30 个必要测试（手工硬编码 + 穷举 + 边界）
tests/reference.py       # 独立参考实现（不导入/不调用被测核心求期望值）
README.md, requirements.txt
```

测试期望值来源：手工算例硬编码（纯分化树、单支丢失、单次/多次复制）；
以及 `tests/reference.py` 中以“祖先集合 + 朴素爬升”另写的独立参考，
在 3 叶与 4 叶物种树的**全部二叉形状 × 2..4 基因叶全部形状 × 全部叶标记组合**
（允许一物种多基因）上逐字段交叉核对；新增测试还独立枚举全部合法内部节点祖先映射，
核对计数前沿、相同计数见证去重、权重选择和预算耗尽语义。往返/聚合自洽断言仅作补充。
