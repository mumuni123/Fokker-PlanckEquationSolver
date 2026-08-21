# VPFP 联合相空间时间中心能量闭合重构实施方案

## 0. 任务身份与唯一目标

本方案用于替代当前 collisionless background 的 Strang 分裂：

$$
T_x(\Delta t/2)\circ T_u(E,\Delta t)\circ T_x(\Delta t/2).
$$

已确认该分裂的 PPM x swept current 与 u-face kinetic-energy work 不属于同一个离散能量系统。
本方案的唯一目标是在一个新、显式选择的模式中构造：

$$
\frac{f^{n+1}-f^n}{\Delta t}
+D_xF_x^{n+1/2}
+D_uF_u^{n+1/2}=0,
$$

并使同一个中点 phase-space flux 同时定义 density、charge current、u-force 动能功和 Poisson 场功。

当前模式必须保持默认：

```text
analytic-cell-center x velocity
existing Strang PPM x/u split
existing production output and checkpoint behavior
```

新模式名称固定为：

```text
--background-phase-space-mode joint-midpoint-energy
```

该模式在所有 J 阶段完成前严禁成为默认生产模式。

## 1. 已确认事实与非目标

### 1.1 已确认事实

1. Poisson、Gauss、端点 $dx/2$、$G/G^*$、MPI owner、Beam/Tail 边界功、conversion/return ledger 和逐组分连续性均已通过对应审计。
2. 单速度、低 limiter 覆盖状态已有 5%--7% 的功配对误差；PPM/FCT 不是首要根因。
3. 用能量共轭速度表替换 x transport speed 后，所有 table/continuity/MPI/checkpoint 审计通过，但单速度功误差没有下降。
4. 结论：根因是 **两个独立 split flux 没有共同的离散乘积法则**，不是单一速度表、边界或诊断项。

### 1.2 非目标

本方案不允许：

- 修改 `OpenElectrostaticSolver` 的 Poisson 算子、Dirichlet phi 边界、端点权重或物理电极功；
- 以缩放 $E$、$J$、$a_u$、$K$、$\rho$ 的方式补偿能量；
- 使用 Poisson 投影、净电荷扣除、平均场扣除或全局能量重分配；
- 在 J0/J1 修改 Beam、Tail、collision、FCT、PPM 或旧 Strang 生产路径；
- 因新模式未收敛而放宽旧物理门或接受 soft-unconverged 状态；
- 直接删除旧模式、旧 checkpoint reader 或既有测试。

## 2. 新离散系统

### 2.1 网格与未知量

背景 Bulk 分布仍存为 cell-integrated mass：

```text
M[i,j,k] = f_bar[i,j,k] * dx * du[j] * ring[k]
```

新模式只对物理 x cell 组装残差；ghost cell 只按既有 MPI/开放边界规则填充，绝不当作独立未知量。

每个时间步的未知量为：

```text
M_np1[physical x cells, all bulk-resolved u cells]
phi_np1[physical x cells plus Dirichlet endpoints]
```

使用：

$$
M^{n+1/2}=\frac12(M^n+M^{n+1}),
\qquad
\Phi^{n+1/2}=\mathcal G_P(\Phi^n,\Phi^{n+1}),
$$

其中 $\mathcal G_P$ 是当前已经通过 Gate F 的 production discrete-gradient pairing field，
不得改写为最终场、显式场或任意场平均。

### 2.2 x flux

对每个 velocity cell $\alpha=(j,k)$，定义单个中点 x flux：

$$
F^x_{i+1/2,\alpha}
=v^H_{\alpha}\,\widehat M^{n+1/2}_{i+1/2,\alpha}.
$$

初版 J1 必须采用一个明确、反对称、无 limiter 的中心 numerical trace：

$$
\widehat M_{i+1/2,\alpha}^{n+1/2}
=\frac12\left(M_{i,\alpha}^{n+1/2}+M_{i+1,\alpha}^{n+1/2}\right).
$$

禁止在 J1 中调用现有 PPM/FCT；其目的不是先追求正性，而是先验证联合离散能量恒等式。

$v^H_\alpha$ 必须由 J0 的离散 Hamiltonian/kinetic-energy compatibility test 决定。它不是简单复用现有
`vx`，也不是 P3-V.2 被否定的局部替代表。J0 必须给出该速度与 u flux kinetic-energy difference 的
离散弱形式证明；证明缺失时停止，禁止猜测速度投影。

### 2.3 u flux

对内部 u face，定义：

$$
F^u_{i,j+1/2,k}
=a^{n+1/2}_i\,widehat M^{n+1/2}_{i,j+1/2,k},
$$

$$
a^{n+1/2}_i=
\frac{qE^{n+1/2}_i}{m_ec}.
$$

其中 $E_i^{n+1/2}$ 必须由已验收的 face-to-cell gather $G$ 给出。J1 同样采用中心 trace：

$$
\widehat M_{i,j+1/2,k}^{n+1/2}
=\frac12\left(M_{i,j,k}^{n+1/2}+M_{i,j+1,k}^{n+1/2}\right).
$$

u 两端使用零 inflow；任何 outflow 先只进入既有速度边界 ledger。J1 不启用 Tail conversion sink。

### 2.4 同源 charge current 与功电流

charge current 必须由**同一** $F^x$ 定义：

$$
J_{i+1/2}^{n+1/2}
=\frac{q}{\Delta t}
\sum_\alpha F^x_{i+1/2,\alpha}.
$$

注意 $F^x$ 的单位必须在 J0 统一：若它已是每个 full-step swept mass，保留 $1/\Delta t$；
若它是 flux rate，禁止再次除以 $\Delta t$。每个字段的量纲必须写进类型注释和 result 文件 header。

u 功使用同一个 discrete kinetic energy difference：

$$
\Delta K_u
=\sum_{i,j,k}
K_{j,k}\left(M_{i,j,k}^{n+1}-M_{i,j,k}^{n}\right).
$$

J0/J1 的必须证明目标是：

$$
\Delta K_u
=\Delta t\sum_i dx\,E_i^{n+1/2}(GJ^{n+1/2})_i
+B_u^K.
$$

### 2.5 场方程与场能

每个 Newton residual evaluation 用当前候选 $M^{n+1}$ 计算 $\rho^{n+1}$，并调用**现有**
`OpenElectrostaticSolver` 得到 $\Phi^{n+1},E^{n+1}$。不在残差中更新 accepted `EMFields`。

已验收的离散 Poisson identity 必须保持：

$$
\Delta U_E-W_{\rm electrode}
=-\Delta t\sum_f w_f E_f^{n+1/2}J_f^{n+1/2}.
$$

将其与 2.4 合并，J1 的背景无碰撞、无速度边界 outflow 制造态必须满足：

$$
\Delta K_{\rm bulk}+\Delta U_E-W_{\rm electrode}=O(\epsilon).
$$

## 3. 文件与模式边界

### 3.1 新增文件

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
tests/joint_phase_space_midpoint_unit_test.cpp
tests/joint_phase_space_midpoint_mpi_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp
tools/analyze_joint_phase_space_energy.py
```

### 3.2 允许修改的现有文件

```text
src/vpfp_integrator.h / src/vpfp_integrator.cpp
src/main_vpfp.cpp
src/vpfp_checkpoint.h / src/vpfp_checkpoint.cpp
src/vpfp_diagnostics.h / src/vpfp_diagnostics.cpp
CMakeLists.txt
```

J0/J1 不允许修改：

```text
src/open_electrostatic_solver.*
src/conservative_ppm_remap.*
src/beam_pic.*
src/background_tail_pic.*
src/tail_bulk_return.*
src/cylindrical_fp_collision.*
```

### 3.3 新模式 dispatch

在 `VpfpIntegrator` 增加：

```cpp
enum class BackgroundPhaseSpaceMode {
    STRANG_PPM = 0,
    JOINT_MIDPOINT_ENERGY = 1
};
```

要求：

1. 默认永远是 `STRANG_PPM`；默认 hash、checkpoint、结果不能变化。
2. CLI 只接受：

   ```text
   --background-phase-space-mode strang-ppm|joint-midpoint-energy
   ```

3. mode 必须写入 run manifest、checkpoint manifest、physical configuration hash 和 accepted-step diagnostics。
4. restart mode 不同必须拒绝；不增加 production override。
5. 新模式 J1 只允许 `--collision-model none --background-tail-mode none --beam-enabled 0`。
   违反时 CLI 在启动前失败，不能静默退回旧路径。

## 4. 非线性求解器

### 4.1 不可使用的求解方式

禁止复用当前“先完整 x，再完整 u”的 Picard/Strang 试探态作为 joint mode residual，因为那会重新引入根因。
禁止对每个 x cell 单独缩放 $a_u$ 或做 u-flux 能量补丁。

### 4.2 J1 初版求解器

实现 distributed Newton-Krylov 外层，残差由 `JointPhaseSpaceMidpointOperator::evaluate_residual()` 唯一计算。

残差向量必须包含：

```text
R_M[i,j,k]  = joint midpoint FV residual
R_phi[i]    = current production Poisson residual / constrained solve residual
```

允许的初版线性策略：matrix-free GMRES + block-diagonal local phase-space preconditioner。
预条件器只用于线性求解，不能改变 residual、flux 或最终解。

Newton 每次 trial：

1. 深拷贝 candidate buffers；
2. 填充 ghost；
3. 计算中点 flux 和 Poisson candidate；
4. 计算 residual；
5. 仅当 global finite、residual 减小且 line search 接受时更新 candidate；
6. 失败时丢弃 candidate，不修改 accepted state/RNG/ledger/time。

J1 不允许 soft acceptance。所有 MPI rank 每轮必须进入同一 collective 序列；global finite、line-search
accept 和最终收敛均必须使用 MPI consensus。

## 5. 实施任务顺序

### J0：离散恒等式与小网格单元测试

只新增 `joint_phase_space_midpoint.*` 的纯算子辅助函数和 `joint_phase_space_midpoint_unit_test`。

必须完成：

1. 给定小的非均匀 u grid、任意 $M^{mid}$、任意离散 $E$，验证 x/u flux 的离散 summation-by-parts。
2. 验证 $G/G^*$、端点 $dx/2$、u kinetic-energy difference、cell volume 和单位。
3. 覆盖正/负场、单速度、两速度、u 边界零 inflow。
4. result 必须输出 `mass_residual`、`kinetic_work_residual`、`poisson_work_residual`、`combined_energy_residual`。

验收：所有恒等式到稳定求和误差。J0 不调用完整 `fp_solver`，但必须调用 production grid/Poisson helper；
不得复写这些公式。

#### J0 首轮 FAIL 的审计修复（2026-08-20，待重跑）

远端首轮 `output/joint_midpoint_j0.result` 的 mass、cell volume、Hamiltonian velocity、u boundary、
Poisson work 和 kinetic work 均已接近对应制造态的舍入尺度；`status=FAIL` 来自 J0 测试的两个审计错误，
不是联合 flux 公式已经被证伪：

1. $G/G^*$ residual 错把两个同号的 production work 相加：

   $$
   R_{G/G^*}^{old}=W_{pair}+W_{potential}.
   $$

   正确 residual 是：

   $$
   R_{G/G^*}=W_{pair}-W_{potential}.
   $$

2. kinetic identity 用 double 的不同累加顺序比较两个强抵消量，却以净 work 作为误差尺度。
   现改为从实际 `u_flux_rate` 与实际 $K_R-K_L$ 重算 long-double kinetic work，且用
   $\sum |\Delta t(K_R-K_L)F_u|$ 作为稳定求和尺度。

已修改：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
tests/joint_phase_space_midpoint_unit_test.cpp
```

该修复不改变 J0 flux、Poisson helper、任何生产推进或默认模式。必须重跑 J0 后才能决定 J1；
旧 `joint_midpoint_j0.result` 不得再用于判定联合离散本身失败。

#### J0 正式集群验收结果（2026-08-20）

集群最新结果：

```text
output/joint_midpoint_j0.result
```

四类 J0 场景均通过：

```text
case_positive_field_pass=1
case_negative_field_pass=1
case_single_velocity_pass=1
case_two_velocity_pass=1
status=PASS
```

验收摘要：

```text
mass_residual=998.421875
mass_roundoff_bound=2.174829026e7
kinetic_work_residual_max=1.728678234e-25
poisson_work_residual_max=2.407412431e-35
combined_energy_residual_max=1.728678234e-25
g_gstar_residual_max=4.786612997e-35
cell_volume_residual=0
hamiltonian_velocity_residual=3.155443621e-30
u_boundary_flux=0
```

`mass_residual` 的绝对值约为 `998`，但小于逐单元浮点舍入界，属于大质量状态的
可解释舍入误差，不是守恒缺陷。正/负场、单速度、两速度均覆盖，production
`OpenElectrostaticSolver` 的 Poisson/G-G* helper、端点 `dx/2` 权重和 u 端点零通量
均已参与测试。

因此 J0 正式状态为：

```text
J0=PASS
J0_audit=PASS
J0_physical_operator_failure=0
next_allowed_task=J1
default_Strang_PPM_mode_unchanged=1
```

该结果只允许进入 J1 background-only 联合中点 smoke；不得据此提前接入 Beam、Tail、
collision、FCT/PPM 或将 joint mode 设为默认模式。

### J1：background-only 联合中点模式

实现新 mode 的最小闭环，仅 background、无 Beam/Tail/collision、均匀或平滑小扰动初态。

验收：

```text
no soft acceptance
all finite
Newton residual decreases monotonically after accepted line-search steps
mass conservation = roundoff
Gauss/Poisson = existing tolerance
combined energy residual = roundoff-scaled tolerance
```

不得在 J1 使用 FCT/PPM。若中心通量出现负值，报告 `min_mass` 并停止；不要补 limiter。
每个 accepted 或 rejected J1 step 都写入 `joint_midpoint_iteration.dat`，记录
GMRES 维数、`R_M`、`R_phi`、line-search、最小 trial mass 和失败原因；失败态不得写入
`vpfp_step_diagnostics.dat`。

#### J1 首轮 FAIL 的求解器修复（2026-08-21，待重跑）

远端首轮结果：

```text
joint_midpoint_j1_unit.result:
  failure_code=73
  failure_stage=joint_midpoint_line_search
  finite=1
  iterations=1
  residual_linf=2.3754154384855132e-08
  poisson_residual_linf=0

joint_midpoint_j1_smoke:
  step=0
  failure_code=73
  failure_stage=joint_midpoint_line_search
```

这不是 Gauss、质量、能量或中心 flux 已被证明失败。根因是 `advance_joint_midpoint()` 的 directional
finite-difference probe 步长写反：旧实现使用

$$
h_{old}=\sqrt{\epsilon}\max(1,\|d\|_\infty),
$$

使 probe 位移 $h_{old}d$ 随 $\|d\|^2$ 增大。在大 cell mass 状态下，Jacobian probe 远离当前
candidate，导致所有 line-search trial 被拒绝。正确尺度是：

$$
h=\sqrt{\epsilon}
\frac{\max(1,\|M\|_\infty)}{\max(1,\|d\|_\infty)},
$$

使 $\|hd\|_\infty$ 为 candidate 规模的 $O(\sqrt{\epsilon})$ 相对扰动。

已完成本地修复：

1. 使用上述 MPI-global state/direction norm 构造 `fd`；
2. 初始 finite candidate 后立即回填 `joint_midpoint_residual_linf`、
   `joint_midpoint_poisson_residual_linf` 和 `gauss_ok`，避免 line-search failure 被错误报告为 `gauss_ok=0`。
3. 首轮修复后仍有 `failure_code=73`，其原因是 line search 对 exact-zero Maxwellian tail cell 使用
   fraction-to-boundary 上界；只要该 cell 的 Newton direction 为负，$\lambda_{max}=0$，所有 trial 都等于
   原 candidate。现改为：Newton/probe/line-search 阶段允许 signed finite candidate 用于求代数 residual；
   只有收敛后才用统一全局负值门验收。若

   $$
   \min(M^{n+1}) < -4096\epsilon\max(1,\|M^{n+1}\|_\infty),
   $$

   返回 `failure_code=76`、`joint_midpoint_negative_solution`，并记录首个负 x/u 格点。
   这不是接受负分布，也不是 FCT/limiter；它只把“Newton trial 可以带符号”和“最终物理解必须非负”分开。

已修改：

```text
src/vpfp_integrator.cpp
```

必须先重编译、重跑 J1 unit 和 10-step smoke。旧 J1 result 不得用于判定联合中点离散失败。

解释规则：

```text
failure_code=73: 仍是 Newton/line-search 机制错误；输出每次 lambda、trial residual、trial Poisson residual。
failure_code=76: 联合中心 trace 已收敛但得到显著负分布；这是中心 flux 非保正的真实算法问题，停止 J1，
                 不得加 FCT/放宽负值门/接受该状态。
failure_code=75: 联合代数解存在，但能量恒等式失败；定位 joint flux 公式，不得调 solver tolerance。
status=PASS:    才允许进入 J2。
```

#### J1 第二轮 FAIL：当前实现不是完整 JFNK（2026-08-21）

第二轮远端结果：

```text
unit:
  failure_code=73
  iterations=2
  residual_linf: 2.3754e-08 -> 2.1787e-08
  finite=1, gauss_ok=1, poisson_residual_linf=0

smoke:
  failure_code=74
  failure_stage=joint_midpoint_not_converged
  step=0, no accepted step
```

这证明首轮的 probe/zero-tail 锁死已缓解：smoke 已能继续进行 residual iteration，而不再立即在 line search
失败。但它也暴露当前代码没有实现本方案 §4.2 要求的 matrix-free GMRES：

```text
当前：preconditioned residual direction + 单个 directional Jv + 标量 beta
要求：preconditioned matrix-free Newton-Krylov / restarted GMRES
```

单方向 secant 只能沿一个近似方向最小化残差，不能求解 x/u/Poisson 耦合 Jacobian；unit 的第二轮 line-search
失败和 smoke 做满迭代后 74 都是该不完整线性求解器的表现，尚不是 joint flux 或能量恒等式失败。

`joint_midpoint_j1_smoke/vpfp_step_diagnostics.dat` 只有表头是**预期行为**：该文件仅写 accepted step，
而 smoke 在 step 0 失败，没有可写入的 accepted 物理状态。失败诊断应单独写入 failure/iteration 通道，
不能用伪 step 行填充该文件。

下一步必须替换 J1 的单方向 secant 为真正的分布式 restarted GMRES，并新增 accepted/trial 分离的
`joint_midpoint_iteration.dat`：每轮写 `iteration`、`residual_linf`、`poisson_linf`、`fd`、
`gmres_dimension`、`line_search_lambda`、`trial_min_mass`、`accepted` 和 failure reason。
在 GMRES 完成前，不得通过增大 `max_iterations`、放宽 tolerance 或 soft accept 强行推进。

#### J1 GMRES 符号修复与 smoke 输出语义（2026-08-21，待重跑）

当前 restarted GMRES 已接入，但首次实现将左预条件残差直接作为 Krylov RHS：

$$
P^{-1}J\delta=P^{-1}R.
$$

这与 Newton 方程相反。正确系统必须是：

$$
J\delta=-R,
\qquad
P^{-1}J\delta=-P^{-1}R.
$$

旧代码因此构造 $+J^{-1}R$ 的近似上升方向，unit 在 iteration 1 的 `failure_code=73` 是预期结果。
已在 `src/vpfp_integrator.cpp` 中在 GMRES normalization 前显式反转
`preconditioned_residual` 符号；不改变 Jv、preconditioner、flux、Poisson、容差或接受规则。

本轮 smoke 的 `vpfp_step_diagnostics.dat` 只有表头仍是正确现象：没有 accepted step 就不能写物理诊断。
但 smoke 目录同时缺少 `vpfp_failure.dat` 与 `joint_midpoint_iteration.dat`，无法视为已记录的数值失败。
当前实现的 `VpfpDiagnostics::write_failure()` 会在正常 rejected step 写 iteration 文件；重跑后：

```text
若有 rejection：必须同时存在 vpfp_failure.dat 和 joint_midpoint_iteration.dat。
若两者仍缺失：优先检查 scheduler/job-level termination，而不是继续修改数值算子。
```

重跑 J1 unit 与 smoke 前禁止继续修改 J1；新结果应先验证 GMRES RHS 符号修复是否使 residual 下降。

#### J1 第三轮修复：local slab 索引与 signed line-search trial（2026-08-21，待重跑）

对最新 unit/smoke 的只读审计发现两处仍未按 J1 设计接通的实现错误：

1. `advance_joint_midpoint()` 已改为每 rank 持有 local mass slab，但在 Poisson charge deposition 和
   accepted u-work 统计中仍用：

   ```cpp
   base = (grid_.ix_start + ix) * nq;
   ```

   索引 local `candidate`/`accepted_bundle`。除 rank 0 外这越过 local vector 范围，导致 smoke 在
   `evaluate_local_residual()` 路径 SIGSEGV，主循环来不及执行 `write_failure()`，因此没有
   `vpfp_failure.dat` 或 `joint_midpoint_iteration.dat`。现已改为：

   ```cpp
   base = ix * nq;
   ```

2. J1 的设计规定 signed Newton probe/line-search candidate 可用于 residual evaluation，最终才执行
   code-76 正性门。但 line-search 调用仍传入 `allow_negative_probe=false`，与 Jacobian probe 的
   `true` 不一致；含极小零尾部的 candidate 会被 residual evaluator提前拒绝，形成 unit 的
   `failure_code=73`。现已让 line-search trial 也使用 signed residual domain。

3. 最新 smoke 的项目根目录 `7074311.err` 证明没有 failure/iteration 文件是 SIGSEGV，而非诊断 writer
   未调用。`addr2line` 将栈定位到 `JointPhaseSpaceMidpointOperator::evaluate_local_residual()`。
   该函数计算 u-face kinetic-energy difference 时，错误地用包含 x 维的
   `cell_index(ix,j,k)` 索引只含 velocity slot 的 `vg.kinetic_energy[j,k]`。当 $ix>0$ 时越界，
   MPI 在 `advance()` 返回前终止，故不可能生成 `vpfp_failure.dat` 或 `joint_midpoint_iteration.dat`。
   现已分离：

   ```text
   distribution index = (ix * Nv + j) * Nmu + k
   kinetic index      = j * Nmu + k
   ```

   并修正该 u-face `delta_k` 读取。

两个修复只影响 joint J1 mode；不改 production Strang/PPM/PIC/Poisson 和默认 mode。

重跑验收顺序：

```text
1. J1 unit: 不得在 iteration 1 因 code 73 退出；必须输出真实 residual 进展。
2. J1 smoke: 不得 SIGSEGV；若 step 0 被拒绝，必须同时生成 vpfp_failure.dat 和
   joint_midpoint_iteration.dat。
3. 若出现 code 76，判定中心 trace 非保正，停止 J1；若出现 code 75，定位 flux 能量式；
   仅 status=PASS 才进入 J2。
```

为避免 unit 因四维 GMRES normal equation 近相关而重复 code 73，J1 另增加严格的 solver safeguard：
若 GMRES direction 的 12 档 line search 全部失败，再沿已带正确负号的 $-P^{-1}R$ 方向执行同一下降门。
该回退在 `joint_midpoint_iteration.dat` 中以 `gmres_dimension=-1` 标记；它不改变 residual、flux、Poisson、
容差或 accepted-state 规则。若该回退也失败，code 73 才表示两种下降方向均无可接受候选。

#### J1 第四轮修复：MPI local slab 与相对 Gauss 门（2026-08-21，待重跑）

最新 smoke 已不再 SIGSEGV，并正常写出 `vpfp_failure.dat` 与 `joint_midpoint_iteration.dat`。其 20 次迭代
的相空间 residual 从 $4.0512\times10^{-8}$ 缓慢降至 $3.7856\times10^{-8}$；但记录的绝对 Gauss residual
约 $1.4\times10^9$。这不是 Poisson 方程未满足，而是 J1 错把具有单位的

$$
\left|\partial_xE-\rho/\epsilon_0\right|
$$

与无量纲固定阈值 $10^{-8}$ 比较。等离子体密度下两项本身很大，双精度相减的绝对残差可以很大而相对残差
仍为舍入尺度。现已改为按 production Gauss gate 的局部物理尺度：

$$
r_G=
\frac{\|\partial_xE-\rho/\epsilon_0\|_\infty}
{\max(10^{-30},\|\partial_xE\|_\infty,
\|\rho/\epsilon_0\|_\infty)}.
$$

J1 Newton 只用 $r_G$ 与 `poisson_tolerance` 比较；诊断需明确这是无量纲相对量。

同时，J1 改为 local mass slab 后，两个位置仍错误使用 `grid_.ix_start + ix` 索引 local vector：

```text
candidate rho deposition
accepted u-work summation
```

已统一改为 `ix * nq`。此前该错误使非 rank-0 越界，SIGSEGV 发生在 `write_failure()` 前，故没有 failure/iteration 文件。

为避免 unit result 只报告最终 `73`，`joint_phase_space_midpoint_energy_test` 现写出每条
`JointPhaseSpaceIterationRecord`，包括 GMRES dimension、residual、relative Gauss residual、line-search alpha、
trial min mass、accepted 和 failure code。

重跑后验收：

```text
smoke 不得 SIGSEGV；若拒绝，failure 与 iteration 文件必须存在。
unit 的 iteration_log_count 必须 > 0。
relative Gauss residual 必须 <= 1e-8，不能再以 1e9 绝对值阻断。
若仍 code 73，读取 unit iteration_* 字段决定 GMRES 与 fallback 哪一个未下降。
若仍 code 74，确认相空间 residual 的收敛率后实施真正 restarted-GMRES 改进；不得放宽门或 soft accept。
```

### J2：开放背景边界和源项

仅在 J1 PASS 后，接入 production reservoir/absorbing ghost fill 与电极功。先无 Beam/Tail。

必须验证：

$$
\Delta K+\Delta U_E-W_{\rm electrode}=O(\epsilon)
$$

并分别覆盖 left reservoir、right absorbing、双端开放和 MPI shared-face。

### J3：Beam/Tail 和碰撞的逐项接入

按以下顺序，每次仅接入一个物种/算子：

```text
J3a: Beam charge-conserving trajectory current
J3b: Tail PIC kick/current/boundary work
J3c: Bulk-to-Tail conversion
J3d: H10 return
J3e: collision half steps
```

每一项接入前后都要通过 J1/J2 能量闭合、对应物种 continuity 和 checkpoint/restart 回归。
任何一项破坏闭合，立即停止在该子项，不得继续下一项。

## 6. 正性与性能原则

### 6.1 正性

J1 的中心时间格式优先证明能量结构；它不是最终生产正性方案。

若 J1 平滑小扰动已出现显著负值，说明基础 flux/solve 有错误，停止。
若仅强非线性生产态需要正性控制，后续另开任务实现**守恒且能量误差显式记账**的 limiter；
不得将旧 FCT 直接搬入 joint mode。

### 6.2 性能

- J0/J1 小网格完成前，禁止优化与并行重构；
- 所有 joint operator work arrays 在 init 时预分配；
- residual 评估只进行必要的 MPI reduction；
- diagnostics 仅在 accepted state 或测试模式输出；
- 不允许为改善性能降低 Newton 收敛门或启用 soft acceptance。

## 7. 编译与运行命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx

# J0
cmake --build build -j4 --target joint_phase_space_midpoint_unit_test
yhrun -N 1 -n 1 --cpu-bind=cores ./build/joint_phase_space_midpoint_unit_test \
  --case all --result ./output/joint_midpoint_j0.result

# J1 small-grid background-only smoke
cmake --build build -j4 --target fp_solver joint_phase_space_midpoint_energy_test
yhrun -N 1 -n 1 --cpu-bind=cores ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result ./output/joint_midpoint_j1_unit.result

yhrun -N 1 -n 1 --cpu-bind=cores ./build/fp_solver \
  --background-phase-space-mode joint-midpoint-energy \
  --background-x-boundary periodic \
  --collision-model none --background-tail-mode none --beam-enabled 0 \
  --stop-after-steps 10 --diagnostic-level 2 \
  --output-dir ./output/joint_midpoint_j1_smoke
```

J2/J3 命令只能在前一阶段 PASS 后补充；不得预先提交。

## 8. 验收与停止条件

J0/J1 任何一项失败时：

```text
STOP_KEEP_STRANG_PPM_DEFAULT
do not touch Beam/Tail/collision/FCT/production checkpoint path
report first residual term and first nonfinite location
```

只有 J3 全部通过，并且新模式在 K1 coarse/fine 同时满足：

```text
energy_residual_reduction=1
current_pair/full=roundoff scale
no soft acceptance
all existing continuity/Gauss/direct-face/restart tests PASS
```

才允许比较 K2 宏观物理响应。任何“账变好但波形恶化”的情况都必须保留两种模式结果，
不得将 joint mode 设为默认。
