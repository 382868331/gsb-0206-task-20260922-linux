# 基因树与物种树的事件归因（LCA 复演 + 候选横向转移检查）

离线、零第三方依赖的可复用 Python 库：输入有根二叉**物种树**与**基因树**及
基因叶→物种叶映射，输出 LCA 映射、事件表（speciation / duplication）、
逐基因父子边的丢失计数与**具体缺失支系**；并在既有复演结果上检查
**候选横向转移（HGT）**的供体/受体分支是否合法，返回可行转移时间区间或
结构化拒绝原因（叶映射固定，不修改复演结果）。

## 运行环境

- Python **3.14.7**（开发与验证版本；仅用标准库，`from __future__ import annotations`，
  3.10+ 亦可导入，但只承诺 3.14.7）。
- Linux 原生运行，装好 Python 后**无需安装任何其他依赖**（`requirements.txt` 为空占位）。
- 完全离线：无网络、账号、密钥、数据库、中间件、云 API 或 Docker。

## 快速命令

```bash
# 全部测试（44 个用例，含与独立参考的穷举交叉核对）
python -m unittest discover -s tests -v

# 固定输入演示：正常复演 + 候选转移（正常/端点边界/三类拒绝/两类输入错误）+ 原有拒绝边界
python demo.py
```

演示真实计算、无 `sleep`、无预录结果，实测 < 0.1 秒（要求 90 秒内）。

## 公开接口

```python
from reconciliation import build_tree, reconcile, check_candidate_transfer

sp  = build_tree({"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}, "r",
                 name="species")
gn  = build_tree({"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")},
                 "gr", name="gene")
res = reconcile(sp, gn, {"a1": "A", "b1": "B", "a2": "A", "c1": "C"})

# 在复演结果上检查候选横向转移（叶映射沿用 res 中固定的映射）
chk = check_candidate_transfer(res, donor="B", recipient="C")
chk.feasible, chk.feasible_interval, chk.reason
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
  - `duplications / speciations / total_losses` — 累计计数；
  - `species_tree / leaf_mapping` — 输入物种树与**固定**叶映射的不可变快照，
    供候选转移检查等后续步骤复用（既有位置字段顺序不变，旧调用不受影响）。
- 常量：`SPECIATION="speciation"`、`DUPLICATION="duplication"`。

### `check_candidate_transfer(result, donor, recipient, *, node_times=None) -> TransferCheckResult`

在既有 `ReconciliationResult` 上检查一条候选横向转移
（`donor` 供体分支 → `recipient` 受体分支）是否合法；**叶映射固定**为
`result.leaf_mapping`，函数纯读取、不修改任何状态。

- **分支标识**：物种树的分支（两节点之间的边）以其**下端（子）节点**标识；
  根节点没有入边，不能作为分支（传根抛 `InvalidInputError`）。
- `node_times: {节点: int} | None`：节点时间，根最老、沿每条父子边
  **严格递增**；必须覆盖全部节点且不能含树外节点。省略时取缺省时间
  “根 0、每条父子边 +1”（即节点深度）。时间只接受 `int`（拒绝 `bool`/`float`）。
- 分支 `(父 p, 子 c)` 的时间区间为闭区间 **`[time(p), time(c)]`**。
- 合法当且仅当同时满足：
  1. `donor != recipient`（否则原因 `same_branch`）；
  2. 两条分支**互不为祖先**——两个下端节点不在同一条根-叶脉络上
     （否则原因 `ancestor_branches`，**该判定优先于时间判定**）；
  3. 两个闭时间区间**有公共时刻**，即
     `max(下界) <= min(上界)`（否则原因 `disjoint_intervals`）。
- `TransferCheckResult` 字段：
  - `feasible: bool`、`reason: str | None`；
  - `donor_branch / recipient_branch`；
  - `donor_interval / recipient_interval: tuple[int,int]`；
  - `feasible_interval: tuple[int,int] | None` — 合法时为交集
    `(max(下界), min(上界))`；**端点相交给零长度窗口 `(t, t)`**，
    表示转移只能发生在该时刻。拒绝原因为结构化字符串常量：
    `REASON_SAME_BRANCH="same_branch"`、
    `REASON_ANCESTOR_BRANCHES="ancestor_branches"`、
    `REASON_DISJOINT_INTERVALS="disjoint_intervals"`。
- 供体/受体互换不影响合法性与窗口（`feasible_interval` 对称）。

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
6. **候选横向转移检查**（在复演结果上，叶映射固定）：
   - 分支以下端节点标识，`(父 p, 子 c)` 的时间区间为闭区间 `[time(p), time(c)]`；
   - 祖先关系：从一个下端节点沿 `parent` 上溯能到达另一个，即二者在同一条
     根-叶脉络上（含端点相接的父子分支，如 `x` 支与 `y` 支）；
   - 判定顺序：同分支 → 祖先分支 → 区间不相交，全部通过才合法；
   - 可行窗口：`[max(两下界), min(两上界)]`，端点相交给零长度窗口；
   - 分支区间与祖先判定均为对树表的常数/深度级查表，单次检查线性级以内。

## 输入上限

- 单棵树节点数 **≤ 200**（题面接口边界），`max_nodes` 只能为 `1≤n≤200` 的 `int`。
- 满二叉树节点数恒为奇数：199 节点（100 叶）接受、201 节点拒绝，测试已覆盖。
- 标签为非空 `str`；左右子必须不同；不支持多边、非二叉、森林。

## 错误语义（可定位，不静默纠正）

所有错误继承 `ReconciliationError`（也是 `ValueError`）：

| 异常 | 触发 | 定位字段 |
|---|---|---|
| `InvalidInputError` | 类型/取值非法：非映射入参、标签非非空 str、`max_nodes` 为 bool/浮点/NaN/Infinity/越界、候选分支不是物种树非根节点、`node_times` 缺键/含树外键/非 int/非严格递增等 | — |
| `TreeStructureError` | 自环、回边（循环）、多父、多根/根不可推断、根未知、不可达分量、非二叉、左右子相同 | `.tree`、`.node` |
| `LeafMappingError` | 映射缺失/多余键、键非基因叶、目标非物种节点、**非叶映射（映到物种内节点）**、键值类型非法 | `.gene_leaf`（键本身非法时为 `None`） |

- 数值输入：唯一浮点风险参数 `max_nodes` 拒绝 `bool`（`True/False` 不当 1/0）、
  拒绝 `float`（含 `NaN`、`+/-Infinity`，信息中区分 NaN 与 Infinity）、
  拒绝 ≤0 与 >200；`node_times` 同样拒绝 `bool` 与非 `int`。
- **原子性/无状态**：库不持有任何可变状态；`reconcile` 先完成全部校验，
  再一次性构造结果，非法输入时不会产生任何部分结果或外部变更。
  `check_candidate_transfer` 只读入参与复演结果（`ReconciliationResult` 为
  不可变 dataclass），**永远不改叶映射或事件结果**，拒绝只体现在返回值
  （`feasible=False` + `reason`），非法入参才抛错；不存在部分更新。
- 本任务无回溯/搜索过程，不存在“搜索预算耗尽”这一结果；输入非法即确定性拒绝，
  与“证明无解”不存在混淆面，故不设置预算类返回值。

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
- 只实现本题契约（有根二叉、LCA 复演、D/S 两类事件，以及在既有复演结果上的
  候选横向转移合法性检查），不宣称兼容完整行业标准（不做转移-丢失联合代价
  优化、不枚举多候选叶映射、不做非二叉/无定根树和解）。

## 已知限制

- 复演仍只支持有根二叉满树与 LCA（复制-丢失）归因，事件表只有 D/S 两类；
  横向转移部分只做**单条候选转移的合法性判定**（分支时间区间相交 + 互不为
  祖先，叶映射固定），不重算事件、不枚举可行转移集合、不做全局事件代价
  （D/T/L）联合优化与最优解搜索——这些不在本次迭代范围内。
- 节点时间为离散整数、沿父子边严格递增；不接受分支长度浮点时间或不确定时间带。
- 叶映射必须全给定，不做“最适合映射”推断。
- 单文件模块布局（`reconciliation.py`），未打包为可安装发行版；
  使用时把该文件放入项目即可导入。
- 左右子顺序只影响输出排列，不影响计数（镜像不变性有测试覆盖）。

## 文件清单

```
reconciliation.py        # 核心库（唯一交付模块，纯标准库）
demo.py                  # 固定输入演示（复演 + 候选转移检查）
tests/test_reconciliation.py  # 44 个测试（手工硬编码 + 性质 + 边界 + 穷举交叉核对）
tests/reference.py       # 独立参考实现（不导入/不调用被测核心求期望值）
README.md, requirements.txt
```

测试期望值来源：手工算例硬编码（纯分化树、单支丢失、单次/多次复制；
候选转移的端点相交、姐妹/堂表分支重叠窗口、同分支/祖先/不相交三类拒绝）；
以及 `tests/reference.py` 中以“祖先集合 + 朴素爬升/朴素区间相交”另写的
独立参考，在 3 叶与 4 叶物种树的**全部二叉形状 × 2..4 基因叶全部形状 ×
全部叶标记组合**（允许一物种多基因）上逐字段交叉核对复演结果，
并在 **3/4 叶物种树全部分支有序对 × 13 组确定性合法时间**上交叉核对
候选转移判定。往返/聚合自洽断言仅作补充。
