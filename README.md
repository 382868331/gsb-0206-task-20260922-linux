# 基因树与物种树的事件归因（LCA 复演 + 复制/丢失 Pareto 前沿）

离线、零第三方依赖的可复用 Python 库：输入有根二叉**物种树**与**基因树**及
基因叶→物种叶映射，输出事件表（speciation / duplication）、逐基因父子边的
丢失计数与**具体缺失支系**。除旧的 LCA 固定权重入口外，本次新增：

- **不预先固定复制/丢失权重**的协调入口，一次树 DP 返回全部不被支配的
  `(复制数, 丢失数)` 取舍及各自一份完整事件归因；
- 调用方随后给**正整数权重**时，直接从前沿选加权最优（含并列计数），
  不重新搜索。

## 运行环境（隔离 venv 准备）

- Python **3.14.7**（开发与验证版本；仅用标准库，`from __future__ import annotations`，
  3.10+ 亦可导入，但只承诺 3.14.7）。
- Linux 原生运行；建议在仓库目录建隔离虚拟环境（标准库即可，venv 内无需装包）：

```bash
python3.14 -m venv .venver
. .venver/bin/activate
python -m ensurepip --upgrade || true   # 本库不需要 pip；仅为隔离环境一致性可选执行
python -c "import sys; print(sys.version)"   # 应显示 3.14.7
```

  不建 venv 也可以：装好 Python 后**无需安装任何其他依赖**
  （`requirements.txt` 为空占位）；venv 只是推荐的隔离方式。
- 完全离线：无网络、账号、密钥、数据库、中间件、云 API 或 Docker。

## 快速命令

```bash
# 全部必要测试（基线 26 个 + 新增 15 个，含与独立参考的穷举交叉核对）
python -m unittest discover -s tests -v

# 旧固定输入演示：一个正常 LCA 复演结果 + 两个实际触发的拒绝边界
python demo.py

# 本轮迭代演示：Pareto 前沿、正整数权重选优、预算耗尽边界、穷举自检
# 实测 < 0.1 秒（要求不超过 90 秒）；无 sleep、无预录结果
python demo_iteration.py
```

两个演示均为真实计算、无 `sleep`、无预录结果。

## 公开接口

```python
from reconciliation import build_tree, reconcile

sp  = build_tree({"y": ("A", "B"), "x": ("y", "C"), "r": ("x", "D")}, "r",
                 name="species")
gn  = build_tree({"g1": ("a1", "b1"), "g2": ("a2", "c1"), "gr": ("g1", "g2")},
                 "gr", name="gene")
res = reconcile(sp, gn, {"a1": "A", "b1": "B", "a2": "A", "c1": "C"})
```

前沿/权重入口的最小调用例：

```python
from reconciliation import build_tree, reconcile_pareto, reconcile_weighted

sp = build_tree({"y": ("B", "C"), "r": ("A", "y")}, "r")
gn = build_tree({"g1": ("b1", "b2"), "gr": ("g1", "c")}, "gr")
lm = {"b1": "B", "b2": "B", "c": "C"}

front = reconcile_pareto(sp, gn, lm)
front.counts()          # ((0, 1), (1, 0))：少复制多丢失 / 多复制少丢失
front.complete          # True（预算足够时前沿完整）

pick = reconcile_weighted(sp, gn, lm, dup_weight=10, loss_weight=1)
pick.tied_counts        # ((0, 1),)
pick.best.total_losses  # 1（ReconciliationResult，字段与旧入口完全一致）
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

### `reconcile_pareto(species_tree, gene_tree, leaf_mapping, *, budget=100_000_000) -> ParetoResult`

不预先固定复制/丢失权重；树 DP 一次算出全部**不被支配**的 `(复制数 D, 丢失数 L)`：
点 `(d,l)` 被支配当且仅当存在 `(d',l')` 满足 `d'≤d`、`l'≤l` 且至少一处严格小。

- `ParetoResult.points`：`FrontierPoint` 元组，按 `(D, L)` 升序；每个点带
  - `duplications / losses`：整数计数；
  - `mapping`：一份完整合法映射（基因节点→物种节点）；
  - `result`：由该映射构造的 `ReconciliationResult`（与旧入口同构的完整归因）。
- `ParetoResult.complete`：前沿是否完整。`False` 表示**搜索预算耗尽**，
  `points` 只是彼此不支配的已知子集、**可能不全**；这与“证明无解”严格区分
  （整棵基因树映到物种树根恒为合法映射，本题不存在无解）。
  `budget` 与 `combinations_explored` 给出预算口径与实际展开的候选组合数。
- 相同 `(D, L)` 只保留**一份**见证：规范序 = 基因树后序（先左后右）节点上的
  物种标签序列字典序最小者；多条等价路径不构成多个取舍。

### `reconcile_weighted(..., *, dup_weight, loss_weight, budget=...) -> WeightedReconciliation`

调用方给**正整数**权重（拒绝 `bool`/浮点/≤0），从前沿选 `dup_weight*D +
loss_weight*L` 最小的方案，**不重新搜索**：

- `objective`：最优目标值；`tied_counts`：达到该值的全部不同 `(D, L)`；
- `best`：规范序（先按 `(D,L)`、再按映射序）选定的那一份 `ReconciliationResult`；
- `frontier`：决策所用前沿。前沿不完整时抛 `BudgetExhaustedError`，
  异常的 `.partial` 携带不完整 `ParetoResult`，绝不把“没算完”伪装成最优。

与旧入口的关系（重要）：旧 `reconcile` 固定走 **LCA 映射**。在字面事件语义下，
把映射抬到 LCA 之上有时能把两次复制合并为一次（如三个同属一个物种的拷贝），
因此 LCA 映射**不一定**是加权最优点；旧入口行为保持不变，其 `(D,L)`
始终是一个合法计数。新入口负责在整个前沿上按权重选优。

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
5. 复杂度：校验与 LCA 复演均为线性级（LCA 按深度爬升，节点 ≤ 200）。

### 前沿树 DP（本次新增，不是“换权重循环碰运气”）

1. **合法映射**：叶映射固定（本题物种归属固定，不增加转移/HGT 事件）；
   每个基因内节点可映到其两个子映射的任意**共同祖先**（含 LCA 本身与
   LCA 之上的祖先）。事件仍按上面的字面规则判定（任一子映射 == 父映射即
   duplication，否则 speciation），逐边丢失仍为 `d`（dup）/`d-1`（spec）。
2. **状态**：对每个 `(基因节点 g, 物种节点 s)` 维护一张
   `(复制数, 丢失数)` 局部前沿：`g` 映到 `s`、子树各自取最优计数组合。
   叶子只在其归属物种叶处有 `(0,0)`。
3. **转移**：枚举两个子前沿的候选组合，父事件代价只依赖三条映射的物种树
   深度差；每个局部表生成后立即做**严格支配剪枝**——被 `(≤,≤)` 且一处严格
   的计数盖住的状态删除。同计数只留规范见证键（子树后序物种标签序列）最小者。
4. **根前沿合并**：根在各物种节点的局部前沿合并后再支配剪枝，即全局前沿。
5. **预算**：预算以“展开的候选子状态组合数”计；一旦超限，当前结果标记
   `complete=False`（可能少算点，尤其可能漏掉多复制少丢失的方案），
   加权入口据此拒绝宣称最优。
6. 测试中的 `tests/reference.py:reference_pareto_frontier` **枚举全部合法映射**
   独立求前沿与规范见证（不调用被测核心），在 2/3 叶物种树的全部形状 ×
   2..4 基因叶全部形状 × 全部叶标记组合上与 DP 逐计数、逐见证一致，
   并显式断言穷举集合中存在“多复制少丢失 / 少复制多丢失”的双向非支配对。

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
| `BudgetExhaustedError` | 加权入口在前沿枚举完整前耗尽 `budget`，无法证明最优（**不是无解**） | `.partial`（不完整 `ParetoResult`） |

- 数值输入：唯一数值参数 `max_nodes` 拒绝 `bool`（`True/False` 不当 1/0）、
  拒绝 `float`（含 `NaN`、`+/-Infinity`，信息中区分 NaN 与 Infinity）、
  拒绝 ≤0 与 >200。
- **原子性/无状态**：库不持有任何可变状态；`reconcile` 先完成全部校验，
  再一次性构造结果，非法输入时不会产生任何部分结果或外部变更。
- 旧 LCA 入口无搜索过程，不存在预算问题；**新增前沿入口是树 DP 搜索**，
  `ParetoResult.complete=False` 与 `BudgetExhaustedError` 把“预算耗尽、前沿
  可能不全”和“证明无解”严格分开（本题合法映射恒存在，无解面不存在）。

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

- 仅支持有根二叉满树与复制-丢失（D/L）复演；**不支持水平基因转移（HGT）**、
  基因树不定根、多叉节点、分支长度、概率/成本最优搜索。物种归属由叶映射
  固定输入，前沿搜索不引入任何转移事件。
- 前沿的“合法映射”允许把基因节点映到 LCA 之上的祖先；事件判定仍严格沿用
  基线字面规则（子映射==父映射即复制）。这意味着抬升映射在某些算例上可合并
  复制（例如同一物种三个拷贝），旧 LCA 入口不利用这种抬升、行为保持不变。
- 权重只接受正整数；加权选择是对已算前沿的线性目标比较，不做连续/概率权重。
- 叶映射必须全给定，不做“最适合映射”推断。
- 单文件模块布局（`reconciliation.py`），未打包为可安装发行版；
  使用时把该文件放入项目即可导入。
- 左右子顺序只影响输出排列，不影响计数（镜像不变性有测试覆盖）。

## 文件清单

```
reconciliation.py            # 核心库（唯一交付模块，纯标准库）
demo.py                      # 旧固定输入演示（LCA 复演 + 拒绝边界）
demo_iteration.py            # 本轮迭代演示（Pareto 前沿/权重/预算/穷举自检）
tests/test_reconciliation.py # 基线 26 个测试（全部保留为回归）
tests/test_pareto.py         # 新增 15 个测试（手工 + 全合法映射穷举 + 预算 + 校验）
tests/reference.py           # 独立参考（LCA 复演 + 全合法映射前沿穷举，均不调被测核心）
README.md, requirements.txt
```

测试期望值来源：手工算例硬编码（纯分化树、单支丢失、单次/多次复制、
双向非支配前沿）；以及 `tests/reference.py` 中以“祖先集合 + 朴素爬升”
另写的两个独立参考（LCA 复演、全合法映射前沿穷举），在 2/3/4 叶物种树的
**全部二叉形状 × 2..4 基因叶全部形状 × 全部叶标记组合**（允许一物种多基因）
上逐字段/逐计数/逐见证交叉核对。往返/聚合自洽断言仅作补充。
