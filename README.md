# 基因树与物种树的事件归因（LCA 复演）

离线、零第三方依赖的可复用 Python 库：输入有根二叉**物种树**与**基因树**及
基因叶→物种叶映射，输出 LCA 映射、事件表（speciation / duplication）、
逐基因父子边的丢失计数与**具体缺失支系**。

扩展能力：部分或全部基因叶可携带**候选物种叶集合**，库联合选择每叶归属，
最小化 `dup_weight * duplications + loss_weight * total_losses`
（非负整数权重），平局按基因叶 ID 序的物种 ID 序列取字典序最小；
选定后仍按原 LCA 规则复演。单归属输入与旧接口结果完全一致。

## 运行环境

- Python **3.14.7**（开发与验证版本；仅用标准库，`from __future__ import annotations`，
  3.10+ 亦可导入，但只承诺 3.14.7）。
- Linux 原生运行，装好 Python 后**无需安装任何其他依赖**（`requirements.txt` 为空占位）。
- 完全离线：无网络、账号、密钥、数据库、中间件、云 API 或 Docker。

## 快速命令

```bash
# 必要测试（47 个用例，含与独立参考/暴力枚举的穷举交叉核对）
python -m unittest discover -s tests -v

# 固定输入演示：正常复演 + 候选联合归属 + 实际触发的拒绝/预算边界
python demo.py
```

演示真实计算、无 `sleep`、无预录结果，实测 < 0.1 秒（要求约 8 秒内）。

## 公开接口

```python
from reconciliation import build_tree, reconcile

sp  = build_tree({"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}, "r",
                 name="species")
gn  = build_tree({"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")},
                 "gr", name="gene")
res = reconcile(sp, gn, {"a1": "A", "b1": "B", "a2": "A", "c1": "C"})
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
- 常量：`SPECIATION="speciation"`、`DUPLICATION="duplication"`。

### `reconcile_with_candidates(species_tree, gene_tree, leaf_candidates, *, dup_weight=1, loss_weight=1, max_nodes=200, max_relaxations=4_000_000) -> CandidateReconciliationResult`

- `leaf_candidates`: `{基因叶: 物种叶 str 或候选 str 的非空可迭代集合}`。
  单条 str 等价于单候选（与旧 `leaf_mapping` 输入一致，结果不变）；
  集合中重复候选自动去重；候选必须**全是现有物种叶且非空**，否则抛
  `CandidateMappingError`（`LeafMappingError` 子类，可定位基因叶）。
- `dup_weight / loss_weight`：非负 `int` 权重（拒绝 bool/浮点/负数，
  含 NaN/Infinity），目标为
  `min dup_weight * duplications + loss_weight * total_losses`。
- **平局消解**：达到同一最优成本的所有归属中，按基因叶 ID（字典序）排列的
  物种 ID 序列取字典序最小者（精确实现：逐叶固定为仍可达最优的最小候选）。
- `max_relaxations`：树 DP 松弛工作量预算（正 int）。耗尽时**不报错**，
  返回 `complete=False` 的 incomplete 结果。
- `CandidateReconciliationResult` 字段：
  - `complete: bool` — True 表示已证明最优；
  - `assignments: dict[str,str]` — 每个基因叶选定的物种叶；
  - `optimal_cost: int | None` — 已证明最优值（incomplete 时为 None）；
  - `best_cost: int` — 返回映射实际达到的成本（可验证）；
  - `reconciliation: ReconciliationResult` — 对 `assignments` 的标准
    LCA 复演完整事件账本；
  - `relaxation_work: int` — 实际消耗的松弛工作量。

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

## 候选归属搜索算法规约

1. **为什么不能逐叶独立选最近物种**：候选选择会改变祖先节点的 LCA，
   进而改变事件归因与丢失计数，必须联合优化。
2. **树 DP**：按基因树后序，状态 = 基因节点映射到的物种节点（只保留
   子树候选可达的 LCA，故每节点状态数以物种树节点数 ≤200 为界）。
   转移时把本节点事件成本（dup_weight，仅 duplication 计入）与两条
   子边的加权丢失（speciation 丢 d-1、duplication 丢 d，d 为物种树
   路径边数）一并累加；根节点最小值即最优成本。该 DP 是精确多项式
   算法，不做逐叶贪心。
3. **平局消解**：求出最优成本 C* 后，按基因叶 ID 字典序逐叶尝试候选
   （升序），固定为第一个仍使 DP 最优值等于 C* 的候选——精确得到
   “物种 ID 序列字典序最小”的最优归属，而非近似。
4. **预算与 incomplete**：DP 松弛工作量（逐内节点累计
   `|maps(左子)| * |maps(右子)|`）超过 `max_relaxations` 时停止搜索，
   返回 `complete=False`、`optimal_cost=None`，并保留最优已知映射
   （每叶取字典序最小候选）及其可验证的完整事件账本；这与“证明最优”
   严格区分，也与输入非法的报错路径区分。
5. 零权重合法：`dup_weight=loss_weight=0` 时所有归属成本同为 0，
   结果完全由平局消解规则决定（有测试覆盖）。

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
| `CandidateMappingError` | 候选集合为空、含非物种叶/未知节点、键非法、基因叶缺候选等（`LeafMappingError` 子类） | `.gene_leaf` |

- 数值输入：唯一数值参数 `max_nodes` 拒绝 `bool`（`True/False` 不当 1/0）、
  拒绝 `float`（含 `NaN`、`+/-Infinity`，信息中区分 NaN 与 Infinity）、
  拒绝 ≤0 与 >200。
- **原子性/无状态**：库不持有任何可变状态；`reconcile` 先完成全部校验，
  再一次性构造结果，非法输入时不会产生任何部分结果或外部变更。
- 单归属 `reconcile` 无搜索过程，输入非法即确定性拒绝。
  `reconcile_with_candidates` 的搜索预算耗尽与“证明最优”严格分开：
  前者返回 `complete=False` 的最优已知映射，后者返回 `complete=True`
  且 `optimal_cost` 非 None；两者都不是错误，非法输入仍然确定性抛错。

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
  基因树不定根、多叉节点、分支长度、概率模型。
- 候选归属只做本题契约内的联合最优选择（dup/loss 加权和 + 字典序平局），
  不做基因转移、时间估计，也不接外部生物数据库。
- 候选搜索的平局消解需逐叶重跑 DP，最坏情形为 `O(叶数 × 候选数)` 次
  DP；节点 ≤ 200 的题面规模内实测无感，但超大规模候选集应调低
  `max_relaxations` 以换取 incomplete 快速返回。
- 单文件模块布局（`reconciliation.py`），未打包为可安装发行版；
  使用时把该文件放入项目即可导入。
- 左右子顺序只影响输出排列，不影响计数（镜像不变性有测试覆盖）。

## 文件清单

```
reconciliation.py        # 核心库（唯一交付模块，纯标准库）
demo.py                  # 固定输入演示
tests/test_reconciliation.py  # 26 个必要测试（手工硬编码 + 性质 + 边界）
tests/test_candidates.py      # 21 个候选归属测试（含暴力枚举交叉核对）
tests/reference.py       # 独立参考实现（不导入/不调用被测核心求期望值）
README.md, requirements.txt
```

测试期望值来源：手工算例硬编码（纯分化树、单支丢失、单次/多次复制）；
以及 `tests/reference.py` 中以“祖先集合 + 朴素爬升”另写的独立参考，
在 3 叶与 4 叶物种树的**全部二叉形状 × 2..4 基因叶全部形状 × 全部叶标记组合**
（允许一物种多基因）上逐字段交叉核对。往返/聚合自洽断言仅作补充。
候选归属部分另以 `itertools.product` 暴力枚举全部候选组合、用公开
`reconcile` 重算成本作独立参照，在多种候选形状与 5 组权重（含零权重）
下交叉核对最优值与平局映射。
