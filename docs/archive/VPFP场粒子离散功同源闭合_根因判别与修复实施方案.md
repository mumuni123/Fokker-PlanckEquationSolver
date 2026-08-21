# VPFP 场粒子离散功同源闭合：根因判别与修复实施方案

> 执行对象：自动编码智能体；按“较弱模型也能机械执行”的标准书写  
> 性质：代码实施、验收和停止规则，不是讨论稿。  
> 依据：`VPFP公共能量残差_离散功配对审计与修复实施方案.md` 及 Gate A--H。
> 最后更新：2026-08-16。Gate I 已完成，唯一根因为 `C`（时间中心）；当前只允许执行 `TASK_ID=JC`。

## 0. 智能体执行入口（最高优先级）

> **本节优先级高于本文其他所有章节。** 本文是总设计和任务索引，不是一条可一次性执行的提示词。
> 每次智能体任务只能执行一个明确的 `TASK_ID`。即使用户只说“执行本文档”，智能体也只能从
> §0.3 状态表中选择第一个满足前置条件且尚未完成的任务，完成后立即停止并报告，禁止自动继续下一项。

### 0.1 一次任务的固定输入

执行者开始修改前必须在回复或工作日志中原样给出：

```text
TASK_ID=<I0|I1|I2|I3|I4|I5A|I5B|I6|JA|JB|JC|JD|JE|JF|K1|K2|K3|K4>
BASELINE=<本次开始时的 commit；无 Git 时写 unavailable，不得因此判失败>
PRECONDITION_EVIDENCE=<前置 Gate 的 result 文件及 PASS 字段>
ALLOWED_FILES=<本任务允许修改的文件>
FORBIDDEN_FILES=<本任务不得修改的文件>
ACCEPTANCE_COMMANDS=<本任务结束前必须执行的命令>
STOP_CONDITION=<PASS 后停止，或首次失败项>
```

若用户没有给出 `TASK_ID`：

1. 读取 §0.3 状态表和现有 `.result`/日志；
2. 选择**第一个证据完整、前置条件满足、状态不是 PASS** 的任务；
3. 先声明所选 `TASK_ID`，不得同时选择多个任务；
4. 若无法判断状态，输出 `BLOCKED_BY_MISSING_EVIDENCE`，列出缺失文件，不得猜测。

### 0.2 单任务事务协议

每个 `TASK_ID` 必须严格执行以下七步，不能合并：

1. **要求摘录**：仅摘录该任务的公式、接口、允许文件、验收字段和禁止项。
2. **当前实现审计**：对每项要求给出 `SATISFIED / MISSING / INCORRECT / UNVERIFIED`，并引用文件和行号。
3. **修改计划**：只包含 `MISSING/INCORRECT` 项；禁止顺手重构、性能优化或修改下一个 Gate。
4. **实施**：使用生产已有量和生产算子；测试不得复写被测公式制造自洽。
5. **静态检查**：编译目标、`git diff --check`、Python `py_compile`；无 Git 环境只跳过 Git 项。
6. **运行验收**：执行该任务规定的全部正测试和负测试；保存 `.result`，不得只引用终端文本。
7. **停止报告**：输出 §0.8 模板，然后停止。PASS 也不得自动进入下一任务。

任何失败只允许修复**首次失败的本任务实现**。禁止通过以下方式获得 PASS：

- 放宽容差、删除断言或把 FAIL 改成 warning；
- 以累计正负抵消替代逐步闭合；
- 用测试重算量替代生产实际量；
- 令诊断参与推进、接受条件或状态修正；
- 因暂时无法运行 MPI/长跑而伪造 PASS。

### 0.3 唯一状态表

状态只能取：`NOT_STARTED`、`IMPLEMENTED_UNVERIFIED`、`PASS`、`FAIL`、`BLOCKED`。
`代码已写`只能标记为 `IMPLEMENTED_UNVERIFIED`。只有规定结果文件存在且全部验收字段通过，才能标记 `PASS`。

| TASK_ID | 内容 | PASS 的唯一证据 | 当前状态 | 权威证据/限制 |
|---|---|---|---|---|
| I0 | 可空只读接口和 audit on/off 等价 | bitwise 状态/RNG/ledger/checkpoint | **PASS** | Gate I 只读链已通过 |
| I1 | bulk x 实际 swept number | x 正/负/混合/边界/zero-dt 全 PASS | **PASS** | `x_transport_unit.result: status=PASS` |
| I2 | Beam/Tail 轨迹连续性及边界源 | 所有 accepted row `continuity_pass=1` | **PASS** | coarse/fine breakdown 无 bad row |
| I3 | bulk/Tail/Beam 逐 cell 实际功 | `local_work_ledger_pass=1` | **PASS** | `force_work.result: status=PASS` |
| I4 | 四恒等式、实际 `E_pair`、accepted 诊断 | 单 rank 配对与重构 PASS | **PASS** | coarse/fine 均 `reconstruction_pass=1` |
| I5A | 六类 injected-fault 负测试 | 所有故障均被检出 | **PASS** | `injected_faults.result: status=PASS` |
| I5B | 2/5-rank MPI owner 与确定性 | 两档 MPI `.result` 均 PASS | **PASS** | `pairing_mpi_n2/n5.result: status=PASS` |
| I6 | 115 fs 粗/细档唯一根因 | `PASS_ROOT_CAUSE_UNIQUE` | **PASS** | `root_cause=C`, candidate count `1` |
| JC | 唯一根因 C：时间中心修复 | 仅实施 §6.3，局部回归通过 | **NOT_STARTED** | **当前唯一可执行任务** |
| JA/JB/JD/JE/JF | 未命中分支 | 不适用 | **BLOCKED** | 不得修改这些分支的生产离散 |
| K1 | 结构回归 | JC PASS 后 Gate I 与原结构测试全 PASS | **BLOCKED** | 等待 JC |
| K2 | 宏观物理回归 | K1 PASS | **BLOCKED** | 等待 K1 |
| K3 | restart 回归 | K2 PASS | **BLOCKED** | 等待 K2 |
| K4 | 性能回归 | K3 PASS | **BLOCKED** | 等待 K3 |

普通执行者不得无证据把状态写成 PASS。上表 I0--I6 的 PASS 来自已保存的
`output/vpfp_pairing_gate_i/*.result` 和 coarse/fine accepted-step 诊断，不是代码存在即 PASS。

### 0.4 Gate I 任务卡与文件白名单

#### TASK I0：只读接口

- 允许文件：审计结构声明、子算子可空参数、I0 等价测试、`CMakeLists.txt`。
- 禁止文件/行为：生产公式、时间层、边界条件、MPI 通信协议、容差。
- 完成定义：`audit=nullptr` 与 `audit=enabled` 的接受态、RNG、ledger 和 checkpoint bitwise 相同。
- PASS 后动作：报告 I0，停止；不得开始 I1。

#### TASK I1：bulk x swept number

- 允许文件：`conservative_ppm_remap.*`、`vlasov_split_step.*`、
  `field_particle_power_audit.*` 中的 x 容器、`vpfp_x_transport_flux_audit_test.cpp`。
- 必须读取生产 `work_.swept_mass` 的最终接受通量；不得在测试中重新积分 PPM。
- 必须覆盖正/负速度、混合速度、reservoir inflow、左右 outflow、zero-dt 和 audit-null。
- PASS 后动作：报告 I1，停止；不得开始 PIC 接线。

#### TASK I2：PIC 轨迹连续性

- 允许文件：`beam_pic.*`、`background_tail_pic.*`、轨迹审计容器和专用测试。
- 只能快照现有 `ChargeConservingTrajectory1D/FaceCurrentAccumulator`；禁止新增 $qnv$ 电流。
- Beam/Tail 必须分别覆盖域内轨迹、MPI face、左右物理边界和源项符号。
- PASS 后动作：报告 I2，停止。

#### TASK I3：逐 cell 受力功

- 允许文件：`conservative_ppm_remap.*` 的最终 $u_\parallel$ 通量记账、
  `beam_pic.*`/`background_tail_pic.*` kick 记账、专用测试。
- cell work 必须从最终实际通量/动能差累加，所有 limiter/转换决定之后再记录。
- 三物种 cell 求和分别等于修改前已存在的全局 ledger；正负场均测试。
- PASS 后动作：报告 I3，停止。

#### TASK I4：公共配对计算器

- 允许文件：`field_particle_power_audit.*`、`open_electrostatic_solver.*` 的只读 helper、
  `vpfp_integrator.*` 的采集接线、`vpfp_diagnostics.*`、单 rank 配对测试。
- `E_pair` 必须来自 Poisson 离散势--电荷恒等式，禁止复制 `Ex^{n+1}` 或自动选择残差最小候选。
- finalize 必须发生在接受态提交、能量 ledger 和源项完成之后；trial 不得写 accepted 文件。
- 四项重构必须逐步闭合，且输出首个失败 rank/index、物种和区域。
- PASS 后动作：报告 I4，停止。

#### TASK I5A：负测试

- 只允许修改测试、测试夹具和因负测试暴露出的 I0--I4 明确缺陷。
- 六类故障必须逐个具有独立 `detected` 和 `first_bad_index`；只检查总 FAIL 不算通过。
- 测试夹具必须遵守当前生产数组布局；禁止依赖旧 ghost/moment 布局。
- PASS 后动作：报告 6/6，停止。

#### TASK I5B：MPI 测试

- 只允许修改 MPI 测试、face owner/guard 暴露及已定位的 MPI 接线。
- 必须真实使用 2 和 5 rank 分区；每 rank 复制同一全域状态不算 MPI 测试。
- 共享 face 只计一次，物理左右端独立；测试输出引用值与分布式值。
- 本地 fake/stub MPI 只能用于语法检查，不能标记 I5B PASS。
- PASS 后动作：报告两个结果文件，停止。

#### TASK I6：固定 checkpoint 根因选择

- 禁止修改任何生产公式；只允许修复诊断缺列、分析器格式错误或实验配置错误。
- coarse/fine 必须来自同一 checkpoint、同一物理窗口、同一物理配置。
- 输出 multiple/unexplained 时必须停止；不得选择“看起来最大”的分支。
- 只有 `PASS_ROOT_CAUSE_UNIQUE` 才允许用户另开任务执行对应 J 分支。

### 0.5 Gate J/K 防跨阶段规则

1. J 不是一个任务集合；I6 只允许激活 `JA/JB/JC/JD/JE/JF` 中的**一个**。
2. 未激活的 J 分支全部视为禁止文件范围，不能“顺便一起修”。
3. J 实施完成后只运行该分支的局部回归；通过后停止，再单独执行 K1。
4. K1、K2、K3、K4 是四次独立任务，不能在同一修改中混合结构、物理、restart 和性能优化。
5. K2/K4 失败不得回头放宽 Gate I 容差；必须回到产生回归的具体 J 代码。

### 0.6 测试与生产代码隔离规则

测试允许：构造输入状态、调用生产算子、读取生产审计出口、计算误差范数和判断门。

测试禁止：

1. 复制 PPM swept-mass、Poisson、PIC deposition 或功电流公式作为“参考实现”；
2. 使用与生产相同的错误公式同时生成输入和期望值；
3. 只验证有限性、非零或累计抵消，却声称恒等式通过；
4. 用固定的 `=0`、空数组或单 rank 复制模拟 MPI 参考值；
5. 把 bitwise hash 不同自动视为物理失败，或把物理范数接近自动视为 bitwise 相同。

每个测试结果至少包含：

```text
schema_version
task_id
case_count_expected
case_count_executed
positive_pass_count
negative_detected_count
first_failed_case
first_bad_rank
first_bad_index
status
```

缺少任一必需 case 时，状态必须是 `INCOMPLETE`，不能写 `PASS`。

### 0.7 集群与本地证据边界

- 本地无真实 MPI 环境：只能报告 `single_rank_local=PASS`、`mpi=UNVERIFIED`。
- 没有 115 fs checkpoint：只能报告 I0--I5B 状态，I6 为 `BLOCKED`。
- 集群没有 Git：Git hash 写 `unavailable`，不影响物理门；编译目标和结果 schema 仍必须验收。
- 长跑未完成：不得将短跑趋势外推为 0--120 fs 已修复。
- 编译成功不是 Gate PASS，代码审阅也不是运行证据。

### 0.8 每次任务的强制报告模板

```text
TASK_ID:
STATUS: PASS | FAIL | IMPLEMENTED_UNVERIFIED | BLOCKED

Requirements audited:
- requirement -> SATISFIED/MISSING/INCORRECT/UNVERIFIED -> file:line evidence

Files changed:
- path -> reason

Commands executed:
- exact command -> exit code

Result evidence:
- result path
- required field=value

First failure, if any:
- case/rank/index
- absolute residual, scale, tolerance

Forbidden changes check:
- production formula unchanged: yes/no
- threshold unchanged: yes/no
- next Gate touched: yes/no

NEXT_TASK_ALLOWED:
- exactly one task id, or NONE/BLOCKED
```

报告完成后必须结束当前任务。`NEXT_TASK_ALLOWED` 只是供用户决定下一次任务，不是继续执行授权。

### 0.9 用户发给智能体的推荐提示词

不要使用“执行整份文档”。每次使用以下模板，只替换 `TASK_ID`：

```text
阅读 docs/VPFP场粒子离散功同源闭合_根因判别与修复实施方案.md。
本次只执行 TASK_ID=<填写一个任务，例如 I2>。

必须优先遵守文档第0节。先审计当前实现并逐项给出
SATISFIED/MISSING/INCORRECT/UNVERIFIED；只修复本任务缺口。
不得执行下一TASK，不得修改未授权生产公式，不得放宽阈值。

完成后运行本TASK规定的全部验收命令，读取.result字段，按§0.8报告并停止。
集群测试无法在当前环境运行时标记IMPLEMENTED_UNVERIFIED，不得声明PASS。
```

当任务是 I6 时，再增加：

```text
本次禁止修改生产离散。若根因不唯一，输出INCONCLUSIVE或BLOCKED并停止，
不得自行选择Gate J分支。
```

当前 I6 已结束。下一个智能体任务必须使用：

```text
阅读 docs/VPFP场粒子离散功同源闭合_根因判别与修复实施方案.md。
本次只执行 TASK_ID=JC。
PRECONDITION_EVIDENCE=output/vpfp_pairing_gate_i/root_cause_selection.result:
status=PASS_ROOT_CAUSE_UNIQUE,root_cause=C,root_cause_candidate_count=1。
只实施文档 §6.3，禁止执行 JA/JB/JD/JE/JF，禁止自动进入 Gate K。
```

## 1. 结论与目标

目标是修复：

$$
R_{fp}=\Delta U_E+\Delta K_{bulk,force}+\Delta K_{tail,kick}
+\Delta K_{beam,kick}-W_{electrode}.
$$

Gate H 表明 $R_{fp}$ 解释约 99% 的完整残差，定位了上一级根因：
**场--粒子功配对不闭合**。Gate I 随后使用生产 accepted state、同一个 115 fs checkpoint
和相同物理窗口完成粗/细档判别，已将下一级根因唯一定位为 **C：时间中心不一致**。

历史上在 Gate I 之前不能唯一证明“Poisson 时间层错误”，原因是：

1. 其他阶段闭合后，`field_particle_pair` 接近完整余额具有分解恒等式性质。
2. 真正的 `midpoint_poisson + force/kick + final_poisson` stage 组合只解释 72.2%/76.2%，
   未达到原文 80% 单根因门。
3. `conversion_after_force` 仍占约 24%--28%。
4. $x$ remap 动能闭合不代表其 swept-mass 电荷电流与速度受力功电流同源。
5. Gate F 只证明 Poisson 空间恒等式，不证明输运电流、场能和受力功三方兼容。

完成 Gate I 后的权威根因表述：

> 生产离散使用的 midpoint 场与电荷输运所对应的离散梯度场没有形成同一个
> 收敛的场--粒子时间中心状态。连续性、局部受力功、Poisson 空间恒等式、PIC 伴随、
> Bulk/Tail 转换和边界主导已从本轮主根因中排除。修复必须重构统一时间中心的事务 trial，
> 不得修改 Poisson 空间算子或加入事后能量补丁。

最新权威结果：

```text
output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps.result: status=PASS
output/vpfp_pairing_gate_i/checkpoint115_dt025_20steps.result: status=PASS
output/vpfp_pairing_gate_i/root_cause_selection.result:
  status=PASS_ROOT_CAUSE_UNIQUE
  root_cause=C
  root_cause_candidate_count=1
  same_initial_physical_state=1
  same_physical_window=1
```

### 1.1 本任务的明确交付物

以下是整个项目的总交付 backlog，不是一次智能体任务的交付要求。每次只能按 §0 执行一个
`TASK_ID`，不得因为本清单列出了后续工作而跨 Gate 实施：

1. [已完成] 生产实际量的只读审计数据结构和采集接口；
2. [已完成] bulk、Tail、Beam 三种表示统一量纲的 face swept charge 和 cell force work；
3. [已完成] 四个离散恒等式的公共计算器；
4. [已完成] accepted-step-only 的全局与 profile 诊断；
5. [已完成] C++ 正/负/MPI 测试和 Python 分析器；
6. [已完成] 115 fs checkpoint 的粗细档根因选择：`C`；
7. [当前任务] 仅针对 C 时间中心分支的生产修复；
8. 修复前后结构、宏观物理、restart 和性能 A/B；
9. 一份按 §11 模板填写的最终报告。

Gate I 已结束。现在只能开始交付物 7 的 `JC`；其他 J 分支仍然禁止。

## 2. 已排除项

- 已知源账无遗漏或重复：`known_source_minus_accounted=0`。
- Poisson 空间恒等式通过，残差约 $10^{-38}\ \mathrm{J/m^2}$。
- bulk 最终 $u_\parallel$ 通量、Tail kick、Beam kick 的局部动能恒等式通过。
- 碰撞、$x$ remap、H10 的阶段账本残差为舍入级。
- 115--117 fs 的残差功率对 $dt$ 减半近零阶；全局继续减步长不是修复。

### 2.1 证据边界

| 已有结果 | 能证明 | 不能证明 |
|---|---|---|
| Poisson work identity PASS | 固定两个电荷状态时，空间 Poisson/边界功恒等式正确 | 生产空间输运是否提供了与之兼容的电荷通量 |
| bulk u-face identity PASS | 最终 u-face 通量与 bulk 动能差一致 | 该功电流是否等于 x-remap 电荷电流 |
| Tail/Beam kick identity PASS | kick 前后粒子动能差记账正确 | gather 是否与轨迹 deposition 构成伴随 |
| x stage balance 舍入级 | 空间输运不直接改变粒子动能 | x swept charge 是否正确驱动场能变化 |
| source ownership PASS | 已列源项无漏记/重复 | 状态更新内部是否有结构缺口 |
| field-particle 解释约 99% | 缺口位于公共场粒子余额 | 缺口必然是时间中心 |
| Gate I I0--I5B PASS | 连续性、局部功、重构、MPI owner 可用于根因判别 | 时间中心修复已完成 |
| Gate I I6 `root_cause=C` | 当前只允许修复统一时间中心 | JC 能保证 0--120 fs 宏观和性能验收 |

后续任何结论必须落在“能证明”列内。禁止把单个 PASS 扩大解释为整个耦合系统已闭合。

## 3. 硬约束

1. 不修改 Poisson 空间差分、Dirichlet 边界、场能定义和 Gauss 门。
2. 不做全局能量补丁、电场/动量缩放、投影或零模修正。
3. 不用 $qnv$ 或重算矩替代生产实际 swept mass、PIC 轨迹和最终受力通量。
4. 不修改碰撞、flux-interface、H10/controller，除非新审计唯一定位到对应阶段。
5. 诊断只读；开关前后接受态、RNG、ledger、checkpoint 必须 bitwise 相同。
6. MPI 共享 face 只有一个 owner；左入/出、右入/出独立记账，不周期回卷。
7. Gate I 已完成且唯一命中 C；只允许修改 JC 授权的时间中心代码。
8. 大数组仅在 `diagnostic-level=2` 分配；level 0/1 不增加逐步 I/O、全扫描或 collective。
9. 使用 `--restart-dir`；不存在 `--restart` 和 `--overwrite-output`。

### 3.1 固定执行顺序

历史全流程为：

1. **I0 接口接入**：只增加可空审计接口；默认路径 bitwise 回归。
2. **I1 bulk x-flux**：先使 bulk swept-mass 更新恒等式通过。
3. **I2 PIC trajectory**：复用 Beam/Tail 现有轨迹接口，完成物种连续性。
4. **I3 local force work**：完成逐 cell 功与现有全局功的一致性。
5. **I4 pairing calculator**：实现四个总恒等式和区域分解。
6. **I5A unit/negative tests**：六类注入故障逐项被检测，然后停止。
7. **I5B MPI tests**：真实 2/5-rank owner 测试，然后停止。
8. **I6 checkpoint 10/20 步**：唯一根因选择，然后停止。
9. **JC 时间中心生产修复**：当前新任务只实施 C 分支。
10. **K1--K4**：结构、宏观、restart、性能分别作为四个新任务执行。

当前 1--8 已完成，从第 9 步开始。不得因为 JC 未实现而重跑 1--8。

禁止把 I0--I4 和 Jx 合成一次修改。每阶段失败时只修该阶段首次失败项。

### 3.2 文件职责矩阵

| 文件 | 允许修改 | 禁止内容 |
|---|---|---|
| `conservative_ppm_remap.*` | 暴露最终 x swept mass；输出逐 cell u 功 | 新建另一套 PPM；修改生产通量 |
| `vlasov_split_step.*` | 透传可空审计指针 | 改 Strang 顺序 |
| `beam_pic.*` | 只读导出既有轨迹电流和逐 cell kick work | 新建 qv 沉积；改变注入/出流 |
| `background_tail_pic.*` | 复用 `ChargeConservingTrajectory1D` 和 `FaceCurrentAccumulator`；输出逐 cell work | 第二套轨迹算法；改变粒子权重 |
| `field_particle_power_audit.*` | 统一量纲、恒等式、区域归并 | 推进状态；修正生产数组 |
| `vpfp_integrator.*` | Gate I 采集接线；Gate J 命中后实施单分支 | Gate I 阶段改物理 |
| `vpfp_diagnostics.*` | accepted-step-only 写文件 | 写 trial；参与接受条件 |
| `open_electrostatic_solver.*` | 只读 helper；仅 A/C 分支按证据扩展 | 修改 `solve()` 空间算子 |
| `main_vpfp.cpp` | 配置诊断和 Gate J A/B CLI | 硬编码输出目录或测试参数 |
| `CMakeLists.txt` | 新增目标并链接生产实现 | 复制生产公式到测试 |

### 3.3 数据生命周期和所有权

1. `VpfpIntegrator` 持有每步审计 workspace，避免 remap/PIC 子对象互相持有指针。
2. `advance()` 开始时仅在 level 2 清零；试探失败时丢弃，接受后随 `VpfpStepResult` 输出标量。
3. profile 大数组不放入 checkpoint，不进入 physical config hash，也不影响 restart。
4. 第一、第二 x 半步、第一、第二 PIC drift 的积分量分别保存；最终计算器才组合。
5. rank-local face 数组不 gather 成全域数组。标量用一次 packed reduction；profile 各 rank 独立写，
   rank 0 写 manifest。
6. 共享 face 的权威值在子算子完成通信后读取，不允许审计层再次平均。

## 4. Gate I：生产实际量的局部同源审计

### 4.1 新增文件和结构

新增 `src/field_particle_power_audit.h/.cpp`：

```cpp
struct XFaceTransportAudit {
    bool enabled;
    double substep_dt;
    std::vector<double> bulk_number_swept_face;
    std::vector<double> tail_number_swept_face;
    std::vector<double> beam_number_swept_face;
};
struct CellForceWorkAudit {
    bool enabled;
    std::vector<double> bulk_delta_ke_cell;
    std::vector<double> tail_delta_ke_cell;
    std::vector<double> beam_delta_ke_cell;
};
struct FieldParticlePowerAuditResult {
    bool valid;
    int first_bad_rank, first_bad_index;
    double continuity_linf;
    double poisson_transport_residual;
    double force_work_residual;
    double current_pair_residual;
    double conversion_transfer_residual;
    double boundary_residual;
    double reconstructed_full_residual;
    double reconstruction_mismatch;
};
```

数组只保存本 rank 的 `nx_local+1` face 或 `nx_local` cell，不复制完整状态。

同时定义稳定枚举，结果文件只写整数值，Python 再映射名称：

```cpp
enum PairingRootCauseMask {
    PAIRING_CAUSE_NONE       = 0,
    PAIRING_CAUSE_TRANSPORT  = 1 << 0,
    PAIRING_CAUSE_WORK       = 1 << 1,
    PAIRING_CAUSE_TIME       = 1 << 2,
    PAIRING_CAUSE_PIC        = 1 << 3,
    PAIRING_CAUSE_CONVERSION = 1 << 4,
    PAIRING_CAUSE_BOUNDARY   = 1 << 5
};
```

结构默认构造必须将 bool 设为 false、索引设为 -1、累计量设为 0；不得用未初始化的 `inf` 作为
“未执行”标记。未执行状态通过独立 `valid/evaluated` 字段表示。

量纲契约：

| 名称 | 位置 | 单位 | 是否除以 dt |
|---|---|---|---|
| `*_number_swept_face` | face | $\mathrm{m^{-2}}$ | 否 |
| `current_face_x` | face | $\mathrm{A/m^2}$ | 已是完整步平均 |
| `*_delta_ke_cell` | cell | $\mathrm{J/m^2}$ | 否 |
| `E_pair` | face | $\mathrm{V/m}$ | 不适用 |
| pairing residual | 全局/区域 | $\mathrm{J/m^2}$ | 否 |
| residual power | 全局/区域 | $\mathrm{J/(m^2s)}$ | 是 |

任何函数入口和输出文件都必须沿用该表，禁止同名字段在不同测试中改变语义。

### 4.2 bulk 空间实际通量

修改 `conservative_ppm_remap.h/.cpp`、`vlasov_split_step.h/.cpp`：

1. 给 `advect_x()` 增加可空 `XFaceTransportAudit* audit=NULL`。
2. 在现有速度循环内累加**最终用于更新 output.f** 的 `work_.swept_mass[f]`。
3. 权重与 `Species::density` 完全一致；若 `f` 已是 cell-integrated mass，不得重复乘速度体积。
4. 两个 $x$ 半步分别保存，不能提前相加。
5. MPI face 以交换后的权威 left flux 为准；内部 face owner 固定为右 rank。
6. `audit==NULL` 时不分配、不清零、不归约，生产状态 bitwise 不变。

接口必须最终写成：

```cpp
RemapDiagnostics ConservativePpmRemap::advect_x(
    const Species& input,
    Species& output,
    double dt,
    double time,
    const OpenBackgroundBoundary& boundary,
    int mpi_rank,
    int mpi_size,
    XFaceTransportAudit* audit = NULL);
```

`VlasovSplitStep::first_x_half()` 和 `second_x_half()` 在参数末尾增加同一可空指针并原样透传。
不允许保留“测试专用 advect_x 复制版本”。调用方必须为第一、第二半步提供两个不同对象。

`bulk_number_swept_face[f]` 的规范语义是：在该子步内穿过 face $f$、按正 $x$ 方向为正的电子
**数量面密度**，单位 $\mathrm{m^{-2}}$，尚未乘电子电荷、尚未除子步时间。转换为电流时统一使用：

$$
J_{bulk,f}=q_e Q_{bulk,f}/\Delta t_{full}.
$$

两个半步的 $Q$ 先相加，再除完整物理步长；不得分别除半步长后再直接相加，否则会多一倍。

### 4.3 PIC 实际轨迹通量

修改 `beam_pic.h/.cpp`、`background_tail_pic.h/.cpp`：

1. Beam 复用生产 `current_face_x` 或其 shape-difference 来源，不用末态 $qv$ 重沉积。
2. Tail 两段 drift 输出同语义 swept number；跨 rank/边界轨迹截断在真实 crossing。
3. 注入/出流写 boundary source，不混入内部 face divergence。
4. 审计不得改变粒子、RNG、迁移缓冲和 ledger。

禁止新增第二套 PIC 沉积器。必须复用：

- Beam：`begin_current_interval()`、`predict_to_midpoint()`、`finish_from_midpoint()`、
  `finalize_charge_conserving_current()` 及最终 `current_face_x`；
- Tail：`ChargeConservingTrajectory1D::deposit_segment()`、`FaceCurrentAccumulator`、
  `drift_half()`、`finalize_trajectory_current()` 及最终 `current_face_x`。

为了区分两个 drift 半步，只允许在现有累加器上做只读阶段快照：

```cpp
struct PicTrajectoryStageAudit {
    std::vector<double> after_first_drift_current_face;
    std::vector<double> after_second_drift_current_face;
    double left_boundary_number;
    double right_boundary_number;
};
```

第二半步贡献定义为 final accumulator 减 first snapshot，不得重新追踪粒子。Beam 若内部以
shape-density delta 重构最终电流，也必须快照该生产 delta/累加器，而不是用粒子末态重算。

PIC `current_face_x` 已是 $\mathrm{A/m^2}$ 的完整步平均电流；不得再乘 $q_e$ 或再除 $dt$。
bulk raw swept number 与 PIC current 只能在 `FieldParticlePowerAudit` 中完成统一量纲。

### 4.4 按 cell 的实际受力功

1. bulk：在 `advect_u_parallel()` 现有离散动能差位置，由最终 `upar_swept_` 累加
   `bulk_delta_ke_cell[ix]`，不得重构候选通量。
2. Tail/Beam：kick 前后按实际 shape 权重累计动能差到 cell。
3. cell 求和必须分别等于现有 `bulk_upar_face_work`、`tail_kick_work`、`beam_kick_work`。

接口要求：

```cpp
RemapDiagnostics ConservativePpmRemap::advect_u_parallel(
    /* existing arguments */,
    std::vector<double>* local_delta_ke_by_x = NULL);

void BackgroundTailPIC::kick(
    /* existing arguments */,
    double* local_kinetic_work = NULL,
    std::vector<double>* local_delta_ke_by_x = NULL);
```

Beam 不得再 push 一次。应在 `finish_from_midpoint()` 的现有 kick 循环中可选累加
`local_delta_ke_by_x`，或将它保存到 `BeamPIC` 的只读 per-step audit workspace。无论采用哪种
接口，`last_field_work()` 的现有语义不得改变。

cell 归属规则：

1. bulk 按本地物理 x cell；
2. PIC 按 kick 瞬间 $x^{n+1/2}$ 的 CIC 两个 cell 权重分配；
3. 落在物理域外的 shape 份额进入 boundary work audit，不得重归一化回域内；
4. MPI guard 份额交换到 owner 后再求和；
5. 三个数组单位均为 $\mathrm{J/m^2}$ 的本步动能变化，不是功率。

### 4.5 必须闭合的四个恒等式

先用原始积分量确认符号和量纲，再构造电流：

$$
\rho_s^{n+1}-\rho_s^n+D_xQ_{s,swept}=S_{s,boundary},
$$

$$
\Delta U_E-W_{electrode}
=-\Delta t\langle E_{pair},J_{charge}\rangle_h+R_{P\leftrightarrow J},
$$

$$
\Delta K_{s,force}
=\Delta t\langle E_{pair},J_{s,work}\rangle_h+R_{K_s},
$$

$$
R_{fp}=R_{P\leftrightarrow J}
+\Delta t\langle E_{pair},J_{work}-J_{charge}\rangle_h
+\sum_sR_{K_s}+R_{conversion}+R_{boundary}.
$$

`E_pair` 不预设为普通平均。必须从已通过的
`OpenPoissonWorkIdentity::potential_charge_work`、生产 $G/D/G^*$ 和连续性推导。审计端点
平均 face 场、生产 midpoint 场和离散势差配对场；face/cell 映射必须为离散伴随。

实现时必须保留“原始量”和“物理量”两层，避免单位错误：

1. 原始层：cell charge change、face swept charge、cell kinetic-energy change；
2. 物理层：除以完整 $dt$ 后的 $J_{charge}$、$J_{work}$ 和 residual power；
3. 所有离散内积显式乘本地 $\Delta x$，共享 face 使用已定义的 face quadrature 权重；
4. 开放端 face 不得套用内部 face 的半权重，权重必须从已通过的 Poisson identity 推导；
5. 电子电荷使用代码内唯一常量，禁止局部写 `-1.602e-19`；
6. 每个恒等式同时输出 signed residual、absolute contribution sum、roundoff tolerance 和 pass。

候选 `E_pair` 的选择过程必须可复核：

- `endpoint_average_face`：$0.5(E_f^n+E_f^{n+1})$；
- `production_midpoint_face`：当前 `midpoint_fields_.Ex_face`；
- `potential_discrete_gradient`：由 Poisson potential-charge work 和 $G/G^*$ 推导。

分析器不得选择“残差最小者”作为最终定义。只有数学推导和 manufactured test 同时通过的候选
才有效；其余候选只作对照。

### 4.6 输出

level 2 且命中输出步时写：

- `field_particle_power_pairing.dat`：逐步全局量；
- `field_particle_power_pairing_profile_<step>.dat`：指定步 profile。

profile 列：

```text
global_index x_m E_pair
J_charge_bulk J_charge_tail J_charge_beam
J_work_bulk J_work_tail J_work_beam
poisson_transport_density force_work_density pairing_residual_density region_id
```

固定汇总左 5%、核心 90%、右 5%，不得按结果动态划区。

文件格式必须固定：

`field_particle_power_pairing.dat` 首行以 `#` 写 schema/version，随后列顺序固定：

```text
step time_s dt_s accepted
continuity_bulk continuity_tail continuity_beam
poisson_transport_residual
force_work_bulk force_work_tail force_work_beam
current_pair_residual conversion_residual boundary_residual
full_residual reconstructed_residual reconstruction_mismatch
roundoff_tolerance root_cause_mask
```

要求：

1. 只写 `accepted=1` 的状态；失败 trial 另写 `trial_field_particle_power_pairing.dat`，默认关闭。
2. 每个物理步只允许一行；restart 后步号不得重复。
3. 所有数值使用 `std::setprecision(17)`。
4. rank 0 写标量；profile 每 rank 独立写并由 manifest 给出 `rank_count`、全局 index 范围和 hash。
5. level 0/1 不创建这些文件；level 2 但 interval 未命中时只保留接受所需标量，不生成 profile。
6. 诊断写失败不得改变已接受物理状态，但必须令作业以明确 I/O failure code 结束，不能静默丢失。

### 4.6.1 与生产积分器的接线位置

修改 `vpfp_integrator.h`：

1. 给 `VpfpStepResult` 增加一个标量 `FieldParticlePowerAuditResult pairing_audit` 和
   `bool pairing_audit_enabled`；
2. 给 `VpfpIntegrator` 增加持久 workspace，分别保存 bulk x1/x2、Tail/Beam drift1/drift2、
   三物种 cell work；
3. 增加 `set_field_particle_power_audit_enabled(bool)`，仅由 `main_vpfp.cpp` 根据 level 2 设置。

修改 `vpfp_integrator.cpp`，采集点固定为：

1. 调用 `first_x_half()` 时写 bulk x1；
2. Tail `drift_half(dt/2)` 和 Beam `predict_to_midpoint()` 后保存 first-drift snapshot；
3. midpoint Poisson 后只保存候选场，不计算最终结论；
4. `u_full()`、Tail kick、Beam `finish_from_midpoint()` 的原循环中采集逐 cell work；
5. force 后 conversion 前后分别保存 bulk/Tail $N,P_x,K$；
6. `second_x_half()` 写 bulk x2；第二 PIC drift 完成后保存 final trajectory accumulator；
7. final Poisson 完成且所有 ledger globalize 后调用一次 `finalize_field_particle_power_audit()`；
8. 只有 `result.accepted=true` 后才交给 `VpfpDiagnostics::write_accepted_step()`。

Beam 和 no-Beam 的两条生产路径必须调用同一个 finalize helper。禁止复制两份公式。

修改 `vpfp_diagnostics.h/.cpp`：

1. 新增 `write_field_particle_power_pairing_accepted_step(...)`；
2. 从 `write_accepted_step()` 内部调用，不能由 main 另走旁路；
3. level/interval 判断沿用现有 `VpfpDiagnostics` 成员；
4. profile 写出只接收 const workspace，不允许诊断层修改或重新计算物理状态。

修改 `main_vpfp.cpp`：

1. level 2 时启用审计，level 0/1 关闭；
2. Gate I 阶段不新增用户 CLI，避免改变 physical config hash；
3. 审计配置只能进入 diagnostic config hash，不进入 physical config hash。

`finalize_field_particle_power_audit()` 的实现顺序固定为：

```cpp
// Pseudocode: do not change the order without updating tests.
verify_workspace_shapes();
combine_bulk_x_half_swept_numbers();
snapshot_existing_beam_and_tail_face_currents();
convert_bulk_swept_number_to_full_step_current();
build_species_charge_change_from_accepted_n_and_candidate_np1();
evaluate_species_continuity_with_boundary_sources();

sum_local_force_work_by_species();
compare_local_work_sum_with_existing_global_ledgers();

evaluate_poisson_transport_pair_for_each_E_pair_candidate();
evaluate_force_work_pair_for_each_species_and_candidate();
evaluate_conversion_and_boundary_terms();
reconstruct_full_field_particle_residual();

pack_scalar_sums_maxima_and_first_bad_location();
MPI_Allreduce_once_for_sums_and_maxima();
resolve_first_bad_rank_and_index_deterministically();
classify_regions_and_root_cause_mask();
```

MPI 注意：

1. SUM、MAX 和首个失败位置若不能放在同一 datatype，可最多使用两个 collective；不得每字段归约。
2. 首个失败位置按最小 rank、再最小 global index 决定，不能依赖到达顺序。
3. rank 0 不得把已经 global 的 bulk work 再次参与 SUM；现有 Gate C 的 rank0-only 语义保持。
4. `all_finite` 在进入数值归约前检查；非有限值不得参与普通 SUM 后再查。
5. profile 不做全局 gather，避免 $N_x$ 大时占用 rank 0 内存。

### 4.7 新测试

新增文件：

- `tests/vpfp_x_transport_flux_audit_test.cpp`
- `tests/vpfp_field_particle_pairing_test.cpp`
- `tests/vpfp_field_particle_pairing_mpi_test.cpp`
- `tools/analyze_vpfp_field_particle_pairing.py`

测试必须链接生产实现，禁止复制生产公式。下面各项均为必须实现内容。

#### 4.7.1 共用测试夹具

新增 `tests/vpfp_field_particle_pairing_test_support.h`，只允许放置状态构造、结果比较和输出辅助，
不得在该文件重新实现 swept-mass、Poisson、PIC deposition 或功电流公式。

夹具至少提供：

```cpp
struct PairingTestState {
    SpatialGrid grid;
    CylindricalVelocityGrid velocity_grid;
    Species bulk_n;
    Species bulk_np1;
    BeamPIC beam_n;
    BeamPIC beam_np1;
    BackgroundTailPIC tail_n;
    BackgroundTailPIC tail_np1;
    EMFields field_n;
    EMFields field_np1;
    std::vector<double> ion_density;
};

PairingTestState make_smooth_bulk_case(...);
PairingTestState make_open_boundary_case(...);
PairingTestState make_beam_single_particle_case(...);
PairingTestState make_tail_single_particle_case(...);
bool bitwise_equal_physical_state(...);
double machine_scaled_tolerance(double scale, double unit_floor);
```

统一构造规则：

1. 小网格只用于单元测试，例如 $N_x=16$ 或 32；速度网格必须调用生产
   `CylindricalVelocityGrid` 初始化，不得自行构造等距替代网格。
2. bulk 平滑分布使用正值 Maxwellian 乘低幅空间扰动：
   $f=f_M[1+10^{-3}\cos(2\pi x/L)]$，保证没有 limiter、尾部转换和边界异常介入。
3. ion density 必须由初始电子密度构造，使初态满足离散 Poisson；不能直接假设连续解析中性。
4. 随机 PIC 测试使用固定 RNG state；能用单粒子确定性轨迹时不得使用统计噪声。
5. 每个测试先保存 bulk 数组、PIC 粒子字节、RNG、ledger 和场，再开启审计执行同一生产调用；
   关闭审计的对照路径也执行一次。最终逐字段或逐字节比较。
6. 结果文件固定为 `key=value`，最后一行必须为 `status=PASS|FAIL`；失败时写
   `first_failure_case`、`first_bad_rank`、`first_bad_index`、绝对残差、尺度和容差。

#### 4.7.2 `vpfp_x_transport_flux_audit_test.cpp`

该测试只验证 `ConservativePpmRemap::advect_x()` 暴露的通量就是实际更新通量，不测试
Poisson 或速度受力。

必须直接调用生产：

```cpp
ConservativePpmRemap remap;
remap.init(grid, velocity_grid);
RemapDiagnostics d = remap.advect_x(
    input, output, dt, time, boundary, rank, size, &audit);
```

若为避免破坏现有 API 而把审计参数放入重载，测试仍必须调用最终进入生产的同一个内部实现。

必须实现以下 case：

1. **uniform-positive-vx**：只占据一个 $u_\parallel>0$ 速度单元，周期内部不含物理边界。
   检查每个 cell：
   `output-input == swept_left-swept_right`。
2. **uniform-negative-vx**：只占据 $u_\parallel<0$ 单元，确认符号反转且同一恒等式成立。
3. **mixed-velocity**：同时填充正负速度，确认速度积分后的 face swept number 等于逐速度
   `work_.swept_mass` 的生产累加结果；测试只能读取生产调试出口，不能另写 PPM 积分。
4. **left-reservoir-inflow**：设置左侧 incoming 分布和正速度，验证左边界 source、内部 divergence
   和总粒子数变化闭合。
5. **left-outflow/right-outflow**：分别构造负速度左出和正速度右出，验证出流符号和边界 owner。
6. **zero-dt/zero-vx**：swept number 必须精确为 0，输出与输入 bitwise 相同。
7. **audit-null-equivalence**：`audit=NULL` 与启用审计所得 `output.f`、moments、诊断 ledger
   bitwise 相同。

逐 case 计算：

$$
R_i=M_i^{out}-M_i^{in}-Q_{i-1/2}+Q_{i+1/2},
$$

并输出 `mass_update_linf`、`global_mass_balance`、`boundary_flux_balance`、
`audit_null_bitwise_equal`。所有非零残差使用 §5 的机器缩放门。

#### 4.7.3 `vpfp_field_particle_pairing_test.cpp`

该测试验证单 rank 下四个恒等式和物种分解。必须调用生产
`VlasovSplitStep`、`OpenElectrostaticSolver`、`BeamPIC`、`BackgroundTailPIC` 和新增
`FieldParticlePowerAudit`；不得用测试代码手动推进状态。

必须实现以下 case：

1. **bulk-continuity-only**：$E=0$，只执行两个 bulk x 半步；连续性必须通过，受力功精确为 0。
2. **bulk-force-plus/minus**：固定平滑密度，分别施加 $+E$、$-E$，只执行生产
   `advect_u_parallel()`。按 cell 功求和与现有全局 bulk work 一致，反场后符号正确。
3. **bulk-poisson-force**：执行无碰撞、无转换、无 PIC 的完整 x-u-x + 两次 Poisson 链。
   输出三种候选 `E_pair` 的 $R_{P\leftrightarrow J}$、$R_K$ 和 $R_{fp}$，不得预先指定赢家。
4. **beam-single-particle**：单个 Beam 宏粒子不跨边界，执行 drift-kick-drift；轨迹连续性、
   `current_face_x`、kick 动能差和 gather/deposit 配对均通过。
5. **beam-boundary-crossing**：粒子在一步内从左右物理边界出流各一次；边界 source 与轨迹
   divergence 必须解释电荷变化，越界后不再产生内部电流。
6. **tail-single-particle**：与 Beam 类似，但调用 `BackgroundTailPIC` 生产 drift/kick。
7. **combined-no-transfer**：bulk、Tail、Beam 同时开启，碰撞、转换和 H10 关闭。完整重构误差
   必须为舍入级，各物种分量之和等于总量。
8. **conversion-accounting**：只启用一次已通过的 flux-interface 转换，验证转换残差是独立项且
   不与场功重复计数；转换前后 $N,P_x,K$ 的现有事务门必须通过。
9. **nonzero-dirichlet**：使用非零 $\phi_L,\phi_R$，验证
   `boundary_energy_work()` 的电极功进入且只进入一次。
10. **read-only-equivalence**：对上述 combined case 比较审计关闭/开启后的 bulk、Tail、Beam、
    fields、RNG 和 ledger，全部 bitwise 相同。

每个 case 必须输出：

```text
continuity_residual_abs continuity_scale continuity_tolerance continuity_pass
poisson_transport_residual_abs poisson_transport_scale poisson_transport_pass
force_work_residual_abs force_work_scale force_work_pass
current_pair_residual_abs current_pair_scale
conversion_transfer_residual_abs conversion_transfer_scale
boundary_residual_abs boundary_scale
full_residual_abs reconstructed_residual_abs reconstruction_mismatch
reconstruction_tolerance reconstruction_pass
```

#### 4.7.4 `vpfp_field_particle_pairing_mpi_test.cpp`

该测试必须使用与生产相同的空间分区，分别由 `yhrun -n 2` 和 `-n 5` 执行。禁止在测试中
把全局状态复制到每个 rank 后绕过 halo/face 协议。

必须实现：

1. 平滑 bulk 分布跨越所有 rank，正负速度均非零。
2. 至少一条 Beam 和一条 Tail 轨迹跨 MPI 内部边界，但不跨物理边界。
3. 另一个 case 同时包含左 reservoir inflow 和右物理 outflow。
4. 每个内部共享 face 输出两个 rank 的原始候选值、最终 owner 值和 owner rank。
5. 统计 `shared_face_duplicate_count`、`shared_face_mismatch_count` 和
   `shared_face_owner_error_count`，三者必须为 0。
6. 将 MPI 结果与同一全局初态的单 rank 参考比较。全局连续性、各物种总 swept number、功和
   $R_{fp}$ 必须在确定性归约容差内一致。
7. 所有局部标量先写入一个 packed `MPI_Allreduce`；测试不得因每项单独归约形成不同执行序。

结果至少包含：

```text
mpi_size
global_state_definition_hash
shared_face_count
shared_face_duplicate_count
shared_face_mismatch_count
shared_face_owner_error_count
single_rank_reference_available
continuity_mpi_vs_reference_abs
work_mpi_vs_reference_abs
full_residual_mpi_vs_reference_abs
status
```

#### 4.7.5 负测试必须能识别错误

`vpfp_field_particle_pairing_test.cpp` 增加 `--case injected-faults`，错误只能注入审计输入副本，
不得修改生产状态或生产算子：

1. 将 `J_charge` face 数组循环平移一格；
2. 翻转一个内部 face 的符号；
3. 将一个 MPI shared face 重复计数；
4. 从 boundary source 中删除一次右出流；
5. 将 `J_work_bulk` 错位到相邻 cell；
6. 将电极功重复加入一次。

每个错误必须使对应门 FAIL，并准确报告首个 bad index。若任一注入错误仍得到 PASS，则
`negative_test_effective=0`，整个 Gate I 失败。负测试不要求执行生产长跑。

#### 4.7.6 Python 分析器实现

`tools/analyze_vpfp_field_particle_pairing.py` 必须支持两种模式：

1. `--run <目录>`：分析单个短跑；
2. `--coarse <目录> --fine <目录> --source-checkpoint <源checkpoint目录>`：相同 checkpoint/物理窗口的
   联合根因选择。`--source-checkpoint` 指向真正的 `$CHECKPOINT_115`，分析器从它的 `manifest.txt`
   读取 `step`、`time`、`physical_config_hash`、`mpi_size`、`nx_global` 作为实验身份；禁止从
   coarse/fine 输出目录的 manifest 反推身份。时间步缩放正是本次对照变量，不得用 coarse/fine 的
   不同 `dt` 判定物理配置不同。

分析器必须：

1. 读取 accepted-step 文件，拒绝 trial/rejected rows；
2. 检查步号连续、接受步数、split/failure/rollback；
3. 用有符号累计、绝对累计和单位物理时间功率三种尺度报告每个残差；
4. 分别报告 bulk/Tail/Beam、左/核心/右区域；
5. 验证四项重构逐步成立，而不是只检查累计抵消；
6. 输出 §5 的 A--F 解释比例和 `root_cause_candidate_count`；
7. 只有一个候选满足阈值时输出 `PASS_ROOT_CAUSE_UNIQUE`；
8. 多个候选输出 `INCONCLUSIVE_MULTIPLE_CAUSES`，无候选输出
   `INCONCLUSIVE_UNEXPLAINED`；
9. 缺文件、NaN/Inf、不同 checkpoint、不同初态或不同物理窗口必须输出
   `INVALID_COMPARISON`，不得降级为物理 FAIL；
10. 退出码固定：PASS=0，物理 FAIL/INCONCLUSIVE=2，输入或格式错误=3。

#### 4.7.7 CMake 接入

`CMakeLists.txt` 必须：

1. 三个 C++ 目标链接对应生产 `.cpp`，并使用与 `fp_solver` 相同的速度网格编译定义；
2. 单 rank 测试可加入 CTest；MPI 2/5 rank 测试只给目标，不在登录节点自动运行；
3. 不把 `main_vpfp.cpp` 链入测试；
4. 不复制 `VPFP_SOURCES` 后再删除不确定文件。优先复用现有 `VPFP_CORE_SOURCES`；
5. Python 分析器和测试脚本先执行 `python3 -m py_compile`。

Gate I 的通过定义是：正测试全部 PASS、注入错误全部被拒绝、MPI 2/5 rank 均通过、审计路径
bitwise 只读，并且分析器能从 115 fs 粗细档唯一选择一个根因分支。

## 5. Gate I 验收和根因选择

舍入门：

$$
T_{round}=\max(10^{-12}S_{unit},512\epsilon_{double}S),
$$

$S$ 为恒等式绝对项之和，禁止除以近零量。

Poisson identity crosscheck 不得复用上式的 reconstruction scale。它是
$\Delta U_E-W_{electrode}-W_{\phi\Delta\rho}$ 中大量交换项相减的独立交叉检查，使用：

$$
T_P=\max\left(4T_{round},4096\epsilon_{double}S_P\right),
$$

其中 $S_P$ 是 Poisson identity 各绝对项的独立尺度。这是对审计量量纲/尺度的纠正，
不是放宽生产能量门；它不参与状态推进或接受条件。

必须输出并通过：

```text
audit_read_only_bitwise_equal=1
rng_bitwise_equal=1
ledger_bitwise_equal=1
shared_face_duplicate_count=0
continuity_identity_pass=1
local_work_sum_matches_existing_ledger=1
full_residual_reconstruction_pass=1
all_finite=1
```

唯一根因门：

- **A 输运/Poisson**：连续性通过，$|R_{P\leftrightarrow J}|$ 解释至少 50%。
- **B 功电流/电荷电流不同源**：Poisson 配对和局部 kick 通过，
  $\Delta t\langle E,J_{work}-J_{charge}\rangle$ 解释至少 50%。
- **C 时间中心**：同一空间电流下离散梯度场通过、生产 midpoint 场失败，差值解释至少 50%。
- **D PIC 伴随**：bulk-only 通过，Beam-only 或 Tail-only 失败。
- **E 转换**：独立转换残差解释至少 20%，且非重复项。
- **F 边界**：两端 5% 解释至少 80%，核心低于 10%。

无唯一分支时输出 `INCONCLUSIVE`，禁止混改多个模块。

### 5.1 根因选择的确定算法

Python 分析器必须按以下顺序执行，不允许智能体凭观察手工选分支：

1. 先验证输入有效、相同 checkpoint、相同初态和相同物理时间窗口；
2. 再验证 bitwise、连续性、局部功求和和完整重构结构门；
3. 结构门失败时输出 `INVALID_AUDIT_STRUCTURE`，不得进入 A--F；
4. 对每步计算各候选的 signed sum 和 absolute sum；
5. 独立交换尺度定义为
   $S_{exchange}=\sum(|A|+|B|+|D|+|E|+|R_{boundary}|)$。候选 C 是两个场配对候选之差，
   不是恒等式中的独立能量项，因此不进入分母；
6. 解释比例使用 $F_c=\sum|R_c|/S_{exchange}$，并同时报告 signed sum、absolute sum
   和单位物理时间功率。因此 $F_C>1$ 在数学上允许，不是分析器错误；
7. 候选必须在 coarse/fine 两档都超过门槛，且主导区域/物种一致；
8. 若 A 与 C 同时命中，先检查同一个 `J_charge` 下不同 `E_pair` 的差值：
   - 空间 Poisson 配对对所有场候选都失败，选 A；
   - 离散梯度候选通过而 production midpoint 失败，选 C；
9. 若 B 与 D 同时命中，先做 species-only：
   - bulk-only 失败，保留 B；
   - 仅 PIC 物种失败，选 D；
10. E 只有在 conversion residual 不包含于 A/B/C 的重构项时才能命中；
11. F 是位置限定符。若 F 与 A--E 同时出现，输出例如 `A_BOUNDARY`，不是第二独立根因；
12. 两个无法按上述规则消歧的候选均超过门槛时，输出
    `INCONCLUSIVE_MULTIPLE_CAUSES`；
13. 根因结果必须列出最坏 10 步，不能只列最大一步；
14. `order_A/B/C/full` 必须报告但不单独决定分支。本轮 `order_C=-0.534`
    说明旧时间离散不呈现正常的缩步收敛，正是 JC 需要修复的风险；它不否定
    coarse/fine 两档均唯一命中 C 的结构判定。

### 5.2 Gate I 阶段性通过标准

| 子门 | 通过条件 | 失败后只允许修改 |
|---|---|---|
| I0 | audit NULL/on bitwise 相同 | 审计接口和 workspace |
| I1 | bulk 每 cell 更新恒等式通过 | x swept 暴露、owner、权重 |
| I2 | Beam/Tail 连续性通过 | 轨迹快照和边界 source |
| I3 | cell work 求和等于旧全局功 | 局部功累加 |
| I4 | 四项重构逐步通过 | pairing calculator、单位、伴随 |
| I5A | 六类注入故障全部被检测 | 测试或已定位的单 rank 接口 |
| I5B | 真实 2/5-rank owner/边界测试通过 | MPI face owner、guard、测试分区 |
| I6 | 唯一根因 | 不修改代码；补充判别实验 |

本轮实际验收结果：

| 项目 | coarse (`dt-scale=0.5`) | fine (`dt-scale=0.25`) | 判断 |
|---|---:|---:|---|
| accepted steps | 10 | 20 | 相同物理窗口 |
| run status | PASS | PASS | 结构有效 |
| $F_A$ | 0.499481 | 0.501486 | 未在两档稳定超门 |
| $F_B$ | 0.500519 | 0.498514 | 未在两档稳定超门 |
| $F_C$ | 1.789266 | 2.226579 | 两档显著命中 |
| $F_D$ | $7.63\times10^{-15}$ | $8.78\times10^{-15}$ | 排除 PIC 功账主导 |
| $F_E$ | $3.02\times10^{-17}$ | $2.20\times10^{-17}$ | 排除 conversion 主导 |
| boundary fraction $F_F$ | 0.086683 | 0.086568 | 排除边界主导 |

`same_initial_physical_state=1`、`same_physical_window=1`，且结果为：

```text
status=PASS_ROOT_CAUSE_UNIQUE
root_cause=C
root_cause_candidate_count=1
```

因此 I0--I6 全部 PASS。这个结论授权进入 `JC`，但不等于 JC 已实现，也不等于
0--120 fs 能量问题已修复。

## 6. Gate J：按唯一分支修复

**当前分支锁：`JC`**。§6.1、§6.2 和 §6.4 仅作历史设计索引，本轮不得实施。
任何对 Poisson 空间算子、电荷连续性电流、PIC gather/deposit、conversion 或边界的修改都视为跨分支。

### 6.1 A：输运/Poisson

以两个 bulk x 半步和 PIC 轨迹为唯一 `J_charge`；在
`OpenElectrostaticSolver` 新增只读配对 helper，不改 `solve()`。只修命中的 owner、boundary
source、量纲或 $G/G^*$。禁止给电流加常数或投影。

具体实施顺序：

1. 先修 raw swept charge 的 cell divergence，使连续性逐 cell 机器闭合；
2. 再由连续性和势差推导 `E_pair` 内积，禁止先调场；
3. 如果仅 MPI face 失败，只改 face owner/交换；如果仅物理边界失败，只改 boundary source；
4. `OpenElectrostaticSolver::solve()` 保持 bitwise，不得把审计电流传入 Poisson 求解器；
5. 修复后由用户分别创建 I1、I4、I5A、I5B、I6 回归任务；当前 J 任务不得连续代跑。
   这些独立回归均通过且 A 不再主导，才允许创建 Gate K 任务。

允许修改文件仅为命中错误的 remap/PIC owner、`field_particle_power_audit.*` 和只读 Poisson
helper。若需要改 `solve()`，立即停止并重新证明 Gate F。

### 6.2 B：功电流/电荷电流

连续性权威为生产 `J_charge`，功权威为生产动能差。修改受力 gather，使同一 face `E_pair`
进入 bulk remap、Tail kick、Beam kick。bulk 使用经测试的 $G^*$；PIC gather 是 shape deposit
的伴随。禁止事后把 `J_work` 替换为 `J_charge`。

具体要求：

1. 在测试中显式构造 gather 矩阵 $S$ 和 deposit 的转置作用，验证
   $\langle E,J\rangle_f=\langle S E,j_p\rangle_p$；
2. bulk 使用生产网格权重验证 $G/G^*$，不能复用旧周期求解器算子；
3. PIC 单粒子测试分别覆盖粒子位于 cell center、face、MPI face 和物理边界附近；
4. 修改 gather 后，密度和轨迹 deposition 不得变化；
5. 受力后的粒子动量、bulk 分布变化应来自新 gather 自然产生，不允许二次动量修正；
6. B 修复必须同时让正场/负场测试通过，防止只修一个符号。

### 6.3 C：时间中心

#### 6.3.0 当前任务定义

```text
TASK_ID=JC
PRECONDITION=output/vpfp_pairing_gate_i/root_cause_selection.result
REQUIRED_PRECONDITION_FIELDS=status=PASS_ROOT_CAUSE_UNIQUE,root_cause=C,
                             root_cause_candidate_count=1
ALLOWED_PRODUCTION_FILES=src/vpfp_integrator.h,src/vpfp_integrator.cpp,
                         src/main_vpfp.cpp,src/vpfp_checkpoint.h,
                         src/vpfp_checkpoint.cpp,src/vpfp_diagnostics.h,
                         src/vpfp_diagnostics.cpp
ALLOWED_TEST_FILES=tests/vpfp_field_particle_time_center_test.cpp,
                   tests/checkpoint_restart_equivalence_test.cpp,CMakeLists.txt
FORBIDDEN_PHYSICS=Poisson spatial operator,x remap flux,Beam/Tail deposition,
                  collision kernels,conversion/return,boundary conditions
STOP_CONDITION=JC local tests PASS or first JC failure; do not run Gate K
```

只在 C 命中时修改 `vpfp_integrator.h/.cpp`，新增事务 trial：

```cpp
bool evaluate_field_particle_trial(
    const Species& collision_input,
    const BeamPIC& beam_n,
    const BackgroundTailPIC& tail_n,
    const EMFields& field_n,
    const EMFields& field_pair_guess,
    FieldParticleTrial& trial,
    FieldParticlePowerAuditResult* audit);
```

要求：

1. collision half1、注入 schedule 和随机结果在接受步内冻结，迭代不消耗新 RNG。
2. 第一 x 半步只算一次。
3. 用 guess 完成 force/kick、转换、第二 drift 和候选末态。
4. 候选末态解生产 Poisson，构造离散梯度 `field_pair_new`。
5. face 场残差和局部功残差同时收敛才接受。
6. trial 失败完整回滚所有状态和 ledger。
7. 初版仅有界阻尼 Picard，不同时加入 Aitken/Anderson。
8. 禁止软接受未收敛状态。

初始结构门：`max_iters=12`、`field_rel_tol=1e-8`、功残差相对交换尺度 `<=1e-8`。

必须新增：

- `enum class FieldParticleCouplingMode { LEGACY, DISCRETE_GRADIENT };`
- `set_field_particle_coupling_mode(...)`；
- CLI `--field-particle-coupling legacy|discrete-gradient`；
- checkpoint manifest 中记录该物理配置；restart 模式不一致必须拒绝，除非专门的测试型 override。

`FieldParticleTrial` 必须拥有候选 bulk/Tail/Beam/final field 和 trial ledger；不得引用接受态的可写
数组。迭代流程必须保证：

1. 每次 trial 从同一个冻结输入恢复；
2. MPI 所有 rank 使用相同 iteration count 和 accept decision；
3. 任一 rank 非有限时 packed consensus 后全体回滚；
4. 收敛量至少包括 `field_rel_l2`、`field_rel_linf`、`pairing_abs`、
   `pairing_relative_to_exchange`；
5. 达到迭代上限返回明确 failure code，不走当前软接受路径；
6. accepted result 记录迭代次数和耗时，供 Gate K 性能门使用。

禁止把 collision 或 stochastic Tail collision 放进 trial 反复执行。若当前 split 不能冻结相关状态，
先新增事务快照/恢复单元测试，不得依赖“同 seed 应该相同”。

#### 6.3.1 实施顺序（不得跳步）

1. 新增 `FieldParticleTrial`，使其完整拥有 bulk/Tail/Beam/field/ledger 候选态；
2. 将 collision half1、Beam injection schedule、Tail collision RNG 固化为 accepted-step 不变输入；
3. 将“给定 `field_pair_guess` 评估候选末态”抽出为无副作用的
   `evaluate_field_particle_trial()`；
4. 由候选末态的生产 Poisson 解构造离散梯度 `field_pair_new`；
5. 使用有界阻尼 Picard 迭代 `field_pair_guess -> field_pair_new`，不加 Aitken/Anderson；
6. 同时检查场残差与局部 pairing 残差，两者通过才 commit trial；
7. 失败时验证状态、RNG、ledger 和 checkpoint hash 未变，返回专用 failure code；
8. 最后增加 CLI/manifest/checkpoint 往返，但不改默认 `legacy` 行为，直到 Gate K 通过。

#### 6.3.2 JC 局部验收标准

- 零场候选一次收敛，状态不变；
- 正/负场均满足 `field_rel_l2 <= 1e-8` 和 `pairing_relative_to_exchange <= 1e-8`；
- 第二次以相同冻结输入评估必须 bitwise 相同；
- 注入故障必须全 rank consensus rollback，不得保留 partial ledger/RNG 消耗；
- legacy mode 必须与 JC 前基线 bitwise 一致；
- coupling mode 的 CLI、manifest、physical config hash、checkpoint read/write 往返通过；
- 只运行 JC 单元/事务测试，PASS 后停止，由用户另行授权 K1。

### 6.4 D/E/F

- D：统一 PIC 密度、轨迹电流和 gather 的 shape 伴随；相对论 kick 使用离散梯度速度。
- E：直接比较转换前后 bulk+Tail 的 $N,P_x,K$，一个 owner 记一次接口能量；不得改动量配平。
- F：分离 reservoir、物理出流、MPI guard 和内部 face；电极功仅由
  `boundary_energy_work()` 给出；核心不做全局修正。

细化约束：

- D 修改后必须通过单粒子 bitwise 可重复测试和 2/5 rank crossing；不得依赖大量粒子统计平均。
- E 的转换前后比较必须使用同一 cell-volume 矩定义；不能混用 cell-center 近似矩。
- E 若误差只处于既有转换容差内，不得为了降低全局能量残差收紧或重写转换。
- F 只允许改边界两侧固定区域；任何作用到核心 face 的修正均视为全局补丁并拒绝。

### 6.5 Gate J 代码提交边界

一个 Gate J 提交只能包含：

1. 一个根因分支的生产修改；
2. 对应的单元/MPI 测试；
3. A/B CLI 和 manifest 字段；
4. 必要诊断更新；
5. 不得顺带做性能重构、格式化全文件或更换数据布局。

提交前后应输出 `git diff --stat` 和按文件的修改理由。若修改超过 8 个生产源文件，必须暂停并
解释为何根因不是局部的；未经用户确认不得继续扩大范围。

## 7. Gate K 验收

重跑 Gate I、Poisson identity、force-work、x-remap 单/MPI、Beam/Tail pusher、checkpoint restart。
从同一 115 fs checkpoint 运行旧模式 10 步、修复模式 10 步、修复模式 dt/2 20 步。

通过标准：

1. $R_{fp}$ 残差功率至少下降 10 倍；
2. 完整重构误差为机器缩放级；
3. 核心残差至少下降 10 倍，不能搬到边界；
4. 完整能量残差功率不高于旧值 20%；
5. 粗细两档不呈显著负阶；
6. 无 split/failure/rollback；
7. Tail/Beam 事务门不退化；
8. 10 步宏观量变化超过 0.5% 必须解释，超过 2% FAIL；
9. level 0 的 100 步 wall time增加不超过 10%；超过 25% FAIL。

### 7.1 K1：结构回归

执行顺序：

1. Gate I 所有正测试；
2. injected-faults 负测试；
3. Poisson、force-work、x-remap、PIC pusher 原有单元测试；
4. 2/5 rank MPI；
5. 115 fs 修复模式 10/20 步。

K1 必须证明：

- 四项恒等式逐步闭合，不依赖累计正负抵消；
- 原 Gate A--H 已通过项没有退化；
- 修复后的主导分支解释比例低于 10%，且没有新分支超过 20%；
- `root_cause_candidate_count=0`。若仍为 1，说明修复无效；若大于 1，说明引入新缺陷。

### 7.2 K2：宏观物理回归

从同一 checkpoint、同一物理窗口比较 legacy/fixed。分析器必须使用相同网格和相同输出时间，
不能按数组行号直接比较错位时间。

至少输出：

$$
\frac{\|E_x^{fixed}-E_x^{legacy}\|_2}{\|E_x^{legacy}\|_2},\quad
\frac{\|n_e^{fixed}-n_e^{legacy}\|_2}{\|n_e^{legacy}\|_2},
$$

以及 field energy、bulk/Tail/Beam kinetic energy、combined spectrum、总粒子数、左右边界累计
通量的相对差。分母接近零时改用绝对尺度门并标记 `relative_valid=0`。

门槛：

- 10 步差异 $\le0.5\%$：PASS；
- 0.5%--2%：`REVIEW_REQUIRED`，必须证明差异与被修复的局部功配对一致；
- $>2\%$：FAIL。

该门不是要求永远贴合 legacy，而是防止结构修复在极短时间内压制真实波形。

### 7.3 K3：checkpoint/restart

执行连续 20 步与“10 步 checkpoint + restart 10 步”。必须比较：

- bulk 完整 `f`；
- Tail/Beam 全粒子状态；
- RNG/persistent state；
- fields；
- coupling mode 和其参数；
- energy、conversion、boundary ledger；
- pairing audit 的累计标量。

物理状态必须 bitwise 相同；若 MPI reduction 导致诊断文本末位不同，可允许机器缩放差异，但
checkpoint 二进制物理状态不能只用宏观量验收。

### 7.4 K4：性能

使用同一 checkpoint、同一 100 步和 `diagnostic-level=0`：

1. legacy 运行一次；
2. fixed 运行一次；
3. 记录 wall seconds、每阶段 timing、MaxRSS、Tail 粒子数和 MPI collective seconds；
4. 比较 wall/accepted-step，而不是作业总排队时间；
5. level 2 另跑 10 步，只用于评估审计开销。

通过判定：

- fixed level 0 增幅 $\le10\%$：PASS；
- 10%--25%：结构可接受但禁止长跑，先优化已定位热点；
- $>25\%$：FAIL；
- level 2/level 0 比值 $>3$：详细审计实现不合格，但不否定物理修复。

### 7.5 Gate K 失败回退表

| 失败项 | 回到 | 不得做 |
|---|---|---|
| 结构重构不闭合 | Gate I4/J 当前分支 | 放宽能量门 |
| 原测试退化 | 产生退化的具体文件 | 同时换另一根因分支 |
| 宏观变化 >2% | Gate J 离散实现 | 用 legacy 波形作强制拟合 |
| restart 不等价 | checkpoint/config/RNG 接线 | 忽略 hash 或只比宏观量 |
| 性能 >25% | workspace/扫描/MPI 优化 | 降低物理精度或放宽收敛 |

## 8. 编译与单元测试命令

> **命令调度规则：本节不是一个可整体复制执行的脚本。** 构建命令可以一次构建全部 Gate I 目标，
> 但运行命令必须按当前 `TASK_ID` 选择。执行 I1 时不得因为下一行存在 I4/I5B 命令而继续运行或修复它们。
> 每个任务写入独立结果文件；禁止以后一个任务 PASS 掩盖前一个任务 FAIL。

### 8.1 公共构建（不代表任何 Gate 通过）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j4 --target fp_solver \
  vpfp_x_transport_flux_audit_test vpfp_field_particle_pairing_test \
  vpfp_field_particle_pairing_mpi_test vpfp_poisson_work_identity_test \
  vpfp_force_work_audit_test conservative_x_remap_test \
  conservative_x_remap_mpi_test checkpoint_restart_equivalence_test
```

构建失败时，当前任务状态为 `FAIL_BUILD`。只修复当前 `TASK_ID` 引入或暴露的首个编译错误；
不得以“顺便整理 CMake”为由修改无关目标。

### 8.2 TASK I1 命令：仅 bulk x 通量

```bash
mkdir -p ./output/vpfp_pairing_gate_i
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/vpfp_x_transport_flux_audit_test \
  --case all --result ./output/vpfp_pairing_gate_i/x_transport_unit.result || exit 11
```

验收时必须读取 `x_transport_unit.result`，确认预期 case 全部执行、`zero_dt_bitwise=1`、
`audit_null_bitwise_equal=1` 和 `status=PASS`。随后停止。

### 8.3 TASK I2/I3/I4 命令：单 rank 配对正测试

I2、I3、I4 可以使用同一测试二进制，但每次任务只能读取属于该任务的 case 和字段。
若测试程序当前只有 `--case all`，必须在进入 I2 前补充独立 case 选择；禁止用一个总 `status=PASS`
代替物种连续性、逐 cell 功和四恒等式的分别验收。

```bash
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/vpfp_field_particle_pairing_test \
  --case all --result ./output/vpfp_pairing_gate_i/pairing_unit.result || exit 12
```

I4 还必须运行已有 Poisson 和 force-work 独立门：

```bash
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/vpfp_poisson_work_identity_test \
  --case all --result ./output/vpfp_pairing_gate_i/poisson_identity.result || exit 14
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/vpfp_force_work_audit_test \
  --case all --result ./output/vpfp_pairing_gate_i/force_work.result || exit 15
```

### 8.4 TASK I5A 命令：六类注入故障

```bash
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/vpfp_field_particle_pairing_test \
  --case injected-faults \
  --result ./output/vpfp_pairing_gate_i/injected_faults.result || exit 16
```

必须逐项检查 `fault_1_detected` 至 `fault_6_detected`，以及
`injected_detected_count=6`。进程退出码为 0 只表示测试程序完成了负测试判定，不能替代字段检查。

### 8.5 TASK I5B 命令：真实 2/5-rank MPI

```bash
for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores ./build/vpfp_field_particle_pairing_mpi_test \
    --case all --result "./output/vpfp_pairing_gate_i/pairing_mpi_n$NP.result" || exit 13
done
```

两个文件都必须包含真实 `mpi_size`、共享面引用值、分布式值、物理边界值和 `status=PASS`。
若测试内部把所有 rank 初始化为 `(rank=0,size=1)`，无论输出为何均判 `INVALID_TEST`。

### 8.6 静态检查（每个修改任务均执行）

```bash
python3 -m py_compile tools/analyze_vpfp_field_particle_pairing.py

ctest --test-dir build --output-on-failure \
  -R 'open_electrostatic_solver_test|vpfp_poisson_work_identity_test|vpfp_force_work_audit_test|conservative_x_remap_test|beam_leapfrog_test|background_tail_pusher_test'
```

注意：

1. CTest 是回归补充，不替代当前任务的专用 `.result`，也不替代 2/5 rank 显式 MPI 测试。
2. 若测试程序没有文档中的 `--result` 或独立 `--case` 参数，必须先补齐；不得依赖终端文本验收。
3. 运行了后续任务命令不代表获得执行后续任务的授权；发现后续任务失败时只记录，不得跨任务修复。

### 8.7 TASK JC 实现后的局部命令

JC 必须新增 CMake 目标 `vpfp_field_particle_time_center_test`。该测试必须直接调用
`evaluate_field_particle_trial()` 生产实现，禁止复写简化时间离散。

编译：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target fp_solver \
  vpfp_field_particle_time_center_test checkpoint_restart_equivalence_test
```

局部测试：

```bash
mkdir -p ./output/vpfp_pairing_gate_jc

for CASE in zero-field signed-field deterministic-replay rollback; do
  yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_time_center_test \
    --case "$CASE" \
    --result "./output/vpfp_pairing_gate_jc/${CASE}.result" || exit 51
done

mkdir -p ./output/vpfp_pairing_gate_jc/checkpoint_tmp
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/checkpoint_restart_equivalence_test \
  --case all \
  --workdir ./output/vpfp_pairing_gate_jc/checkpoint_tmp \
  --result ./output/vpfp_pairing_gate_jc/checkpoint_restart.result || exit 52
```

字段验收：

```text
zero-field.result: status=PASS,iterations=1,state_unchanged=1
signed-field.result: status=PASS,positive_field_pass=1,negative_field_pass=1,
                     field_rel_l2<=1e-8,pairing_relative<=1e-8
deterministic-replay.result: status=PASS,state_bitwise_equal=1,rng_bitwise_equal=1,
                             ledger_bitwise_equal=1
rollback.result: status=PASS,all_rank_consensus=1,state_bitwise_equal=1,
                 rng_bitwise_equal=1,ledger_bitwise_equal=1
checkpoint_restart.result: status=PASS,coupling_mode_roundtrip=1,
                           mismatch_rejected=1
```

上述字段中任何一项缺失都是 `INVALID_TEST`，不能仅凭进程退出码判定 PASS。
JC 局部测试全部通过后立即停止；不在同一任务中执行 §9A 的 K1 生产短跑。

## 9. TASK I6：115 fs Gate I 短跑命令

> **状态：已完成（历史复现命令）。** I0--I5B、粗档、细档和根因选择均已通过。
> 不要为了执行 JC 重跑本节；只有 JC 改动了 Gate I 审计 schema/生产采集语义时，才在 K1 重跑。

先设置实际 checkpoint：

```bash
export OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close
CHECKPOINT_115=/absolute/path/to/checkpoint_target115fs
CP_STEP=5071
```

若 manifest 中 step 不是 5071，必须以 manifest 为准。粗档：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" --restart-source-dt-scale 0.5 \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 --tail-population-control-interval 0 \
  --dt-scale 0.5 --stop-after-steps "$((CP_STEP+10))" --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir ./output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps || exit 21
```

细档必须完整执行以下命令，不要手工编辑粗档输出目录：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" --restart-source-dt-scale 0.5 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 --tail-population-control-interval 0 \
  --dt-scale 0.25 --stop-after-steps "$((CP_STEP+20))" --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir ./output/vpfp_pairing_gate_i/checkpoint115_dt025_20steps || exit 22
```

分析：

```bash
python3 tools/analyze_vpfp_field_particle_pairing.py \
  --run ./output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps \
  --expected-accepted-steps 10 --require-no-split \
  --result ./output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps.result || exit 23

python3 tools/analyze_vpfp_field_particle_pairing.py \
  --run ./output/vpfp_pairing_gate_i/checkpoint115_dt025_20steps \
  --expected-accepted-steps 20 --require-no-split \
  --result ./output/vpfp_pairing_gate_i/checkpoint115_dt025_20steps.result || exit 24

python3 tools/analyze_vpfp_field_particle_pairing.py \
  --coarse ./output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps \
  --fine ./output/vpfp_pairing_gate_i/checkpoint115_dt025_20steps \
  --source-checkpoint "$CHECKPOINT_115" \
  --coarse-dt-scale 0.5 --fine-dt-scale 0.25 \
  --result ./output/vpfp_pairing_gate_i/root_cause_selection.result || exit 25
```

已获得 `status=PASS_ROOT_CAUSE_UNIQUE`，且 `root_cause=C`。本节的两次生产短跑不需要再执行。

`root_cause_selection.result` 至少必须包含：

```text
same_initial_physical_state
same_physical_window
source_checkpoint_path
source_checkpoint_step
source_checkpoint_time
source_physical_config_hash
abs_A_coarse/fine
abs_B_coarse/fine
abs_C_coarse/fine
abs_D_coarse/fine
abs_E_coarse/fine
signed_A_coarse/fine
signed_B_coarse/fine
signed_C_coarse/fine
fraction_A_coarse/fine
fraction_B_coarse/fine
fraction_C_coarse/fine
fraction_D_coarse/fine
fraction_E_coarse/fine
fraction_F_coarse/fine
exchange_scale_coarse/fine
order_A/order_B/order_C/order_full_residual
coarse_worst_steps
fine_worst_steps
root_cause_candidate_count
root_cause
status
```

coarse/fine 各自的结构通过证据保存在
`checkpoint115_dt050_10steps.result` 和 `checkpoint115_dt025_20steps.result`，不得伪造未输出的汇总字段。

## 9A. 当前 JC 与后续 Gate K 命令

> I6 已结束且唯一选中 `JC`。当前先实施 §6.3 和 JC 局部测试；下面的生产短跑属于
> **K1，不得在 JC 实施任务内自动执行**。

JC 修复完成后必须新增显式 A/B CLI；K1 使用以下命令：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" --restart-source-dt-scale 0.5 \
  --field-particle-coupling discrete-gradient \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 --tail-population-control-interval 0 \
  --dt-scale 0.5 --stop-after-steps "$((CP_STEP+10))" --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir ./output/vpfp_pairing_gate_k/fixed_dt050_10steps || exit 31

python3 tools/analyze_vpfp_field_particle_pairing.py \
  --run ./output/vpfp_pairing_gate_k/fixed_dt050_10steps \
  --baseline ./output/vpfp_pairing_gate_i/checkpoint115_dt050_10steps \
  --expected-accepted-steps 10 --require-no-split \
  --result ./output/vpfp_pairing_gate_k/fixed_dt050_10steps.result || exit 32
```

本轮不存在替换为 A/B/D/E/F 的权限。CLI 名称固定为 `discrete-gradient`。

### 9A.1 Gate K 的命令生成规则

Gate K 不得临时手写另一套物理参数。以 §9 粗档完整命令为模板，只允许下表中的替换：

| 测试 | 允许增加/替换 | 输出目录 |
|---|---|---|
| K1 fixed dt | 增加命中分支 CLI | `vpfp_pairing_gate_k/fixed_dt050_10steps` |
| K1 fixed dt/2 | 再增加 dt override，20步 | `fixed_dt025_20steps` |
| K2 legacy | coupling=`legacy`，10步，level2 | `macro_legacy_10steps` |
| K2 fixed | coupling=修复模式，10步，level2 | `macro_fixed_10steps` |
| K4 legacy | coupling=`legacy`，100步，level0 | `perf_legacy_100steps` |
| K4 fixed | coupling=修复模式，100步，level0 | `perf_fixed_100steps` |

除表中项目外，field boundary、Beam、collision、Tail、return、MPI/OpenMP、checkpoint 和初始 step
必须完全相同。每个输出目录运行前必须不存在；禁止覆盖旧结果后进行 A/B。

修复模式 dt/2 命令必须在 §9 细档完整命令上增加：

```bash
--field-particle-coupling discrete-gradient \
--output-dir ./output/vpfp_pairing_gate_k/fixed_dt025_20steps
```

实现该 CLI 时必须同步更新 usage、parse validation、manifest、physical config hash、
checkpoint write/read 和 restart mismatch 文本。任何缺失都使 JC 本地验收失败，不得留到 K3 才补。

### 9A.2 checkpoint/restart 单元命令

```bash
mkdir -p ./output/vpfp_pairing_gate_k/checkpoint_restart_tmp
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/checkpoint_restart_equivalence_test \
  --case all \
  --workdir ./output/vpfp_pairing_gate_k/checkpoint_restart_tmp \
  --result ./output/vpfp_pairing_gate_k/checkpoint_restart.result || exit 41
```

Gate J 新增配置后，必须扩展该测试检查 coupling mode 和参数往返；不能只运行旧测试而不增加断言。

### 9A.3 性能结果要求

若系统没有 `/usr/bin/time`，不得因此跳过性能门。由 `fp_solver` 自身或作业脚本使用 shell
`date +%s.%N` 记录 wall time，并从诊断读取：

```text
accepted_steps
wall_seconds
wall_seconds_per_accepted_step
max_rss_kib
mpi_collective_seconds
pairing_audit_seconds
field_particle_iterations_total
field_particle_iterations_max
```

K4 level 0 中 `pairing_audit_seconds` 应为 0。若不是 0，说明只读详细审计错误进入生产路径。

## 10. 停止规则

任一条件命中即停止：

1. Gate I 不能在机器缩放容差内重构完整残差；
2. 根因不唯一；
3. 审计改变状态/RNG/ledger；
4. 修复只把核心残差搬到边界或其他物种；
5. 修复依赖全局补丁；
6. 10 步宏观变化超过 2%；
7. level 0 性能下降超过 25%；
8. checkpoint/restart 不等价；
9. 无分支证据却要同时修改多个已通过模块。

输出 `BLOCKED_BY_EVIDENCE`、未闭合恒等式、首个 rank/index、绝对尺度和容差；不得改阈值重跑。

### 10.1 统一失败分类

新增审计不得复用现有物理推进 failure code。Python 和 C++ 结果统一使用：

| code | 名称 | 含义 | 是否允许 Gate J |
|---:|---|---|---|
| 0 | PASS | 当前子门通过 | 按顺序继续 |
| 101 | INVALID_INPUT | 文件、schema、checkpoint 或窗口无效 | 否，修测试输入 |
| 102 | AUDIT_MUTATED_STATE | 审计改变状态/RNG/ledger | 否，修只读路径 |
| 103 | CONTINUITY_MISMATCH | 实际输运不能重构密度变化 | 否，只修 I1/I2 |
| 104 | LOCAL_WORK_MISMATCH | 局部功和旧全局功不一致 | 否，只修 I3 |
| 105 | RESIDUAL_RECONSTRUCTION_MISMATCH | 四项不能重构完整残差 | 否，只修 I4 |
| 106 | MPI_FACE_OWNERSHIP_MISMATCH | 共享 face 重复/不一致 | 否，只修 MPI owner |
| 107 | INCONCLUSIVE | 根因不唯一或解释不足 | 否，补判别测试 |
| 108 | MACRO_REGRESSION | 修复压制宏观响应 | 否，回 Gate J |
| 109 | RESTART_REGRESSION | 连续/restart 不等价 | 否，修 checkpoint |
| 110 | PERFORMANCE_REGRESSION | level 0 超过性能门 | 否，做等价优化 |

生产 `fp_solver` 在 Gate I 的只读审计不因 103--107 拒绝物理步；这些状态写入诊断，由测试分析器
决定 Gate。只有非有限、数组越界、MPI 协议错误和 I/O 失败仍可终止生产短跑。Gate J 修复模式
可把结构不收敛作为 trial failure，但必须使用新的明确 code，不能报告 `unknown`。

### 10.2 禁止智能体自行推断的事项

若文档没有明确授权，智能体不得：

1. 自动选择残差最小的 `E_pair`；
2. 把 115 fs 结果外推为 0--120 fs 已修复；
3. 更改 `tail-convert-energy-mev=6.0`、return 参数或碰撞模型；
4. 调整网格、$dt$、MPI rank 数以帮助测试通过；
5. 将边界残差从分析中删除；
6. 把低于某个相对比例的绝对大残差写成 0；
7. 使用旧周期 Vlasov--Ampere 模块的 $G/G^*$ 代替当前开放 VPFP 算子；
8. 以 EPOCH 波形作为单元测试真值；
9. 修改验收阈值后继续同一 Gate；
10. 删除 legacy A/B 模式，直到 Gate K 全部通过。

### 10.3 实施完成前自检清单

```text
[x] I0 audit off/on bitwise
[x] I1 bulk raw swept continuity
[x] I2 Beam/Tail trajectory continuity
[x] I3 per-cell work sums
[x] I4 full residual reconstruction
[x] I5A positive/injected-fault tests
[x] I5B real 2/5-rank MPI tests
[x] I6 unique root cause = C
[ ] JC implementation and local tests
[ ] K1 structure
[ ] K2 macro
[ ] K3 restart
[ ] K4 performance
[ ] final report includes absolute scales and tolerances
```

任一项未勾选时，禁止报告“整个问题已修复”。

## 11. 智能体报告模板

```text
Gate:
修改文件:
明确未修改模块:
编译和测试命令:
结果文件:
bitwise验收:
连续性残差/容差:
Poisson-transport残差及解释比例:
work-charge残差及解释比例:
时间中心候选:
conversion独立残差:
边界/核心比例:
完整重构误差/容差:
唯一根因分支:
最小修复:
修复前后残差功率:
宏观量变化:
性能变化:
是否进入下一Gate:
剩余风险:
```

禁止只报告 PASS；比例必须给出有符号值、绝对值、尺度和零分母有效标志。
