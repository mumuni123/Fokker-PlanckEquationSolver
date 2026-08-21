# VPFP P3 联合 x/u/Poisson 功通量重构实施方案

## 1. 决策与边界

本方案取代“再实施 P1 或基础 P2”的计划。

- P1（把 Poisson 配对账改为 cell-centered 口径）禁止实施。它改变已通过的 face-Poisson 功恒等式的记账目标，不能改变真实场能变化。
- 基础 P2 已存在于 `src/conservative_ppm_remap.cpp` 的
  `ConservativePpmRemap::advect_u_parallel()`：物理 x cell 使用
  $E_i=(E_{i-1/2}+E_{i+1/2})/2$ 计算 $a_{u,i}=qE_i/(mc)$。K1 在此实现下仍有固定能量缺口，故禁止重复实现或微调此插值。
- 禁止按 cell 或全局比例缩放 $a_u$、$E$、$J$ 或动能账以消除残差。这类操作直接改变 Vlasov 方程或隐藏残差。
- P3 的唯一目标是使**实际接受的** x 有限体积通量、u 有限体积受力通量及 Poisson 场功，来自同一离散功系统；所有修正必须随 $\Delta x,\Delta u\to0$ 消失。

当前已确认的缺口为：

$$
J_{\rm charge,f}-G^*J_{\rm force,i}\ne0,
$$

其中 $J_{\rm charge,f}$ 来自 PPM/FCT 最终 x swept flux，$J_{\rm force,i}$ 来自 u 受力步的离散动能变化。
Poisson、$G/G^*$、PIC 边界功、MPI face owner、conversion 和 source ownership 均已通过现有 Gate I/K1 审计，不得修改。

## 2. 生产不变量

下列接口在 P3 全过程必须保持不变，除非某一 P3 单元测试明确证明它们本身错误：

1. `OpenElectrostaticSolver` 的开放边界 Poisson 求解、face 场能和物理端点 $dx/2$ 权重。
2. `ConservativePpmRemap::advect_x()` 的边界 reservoir/absorbing 处理、MPI halo/face owner 语义和最终 swept-number 审计出口。
3. PIC Beam/Tail 的轨迹沉积、开放粒子边界和 Tail/Beam 边界功 ledger。
4. `field_particle_power_audit` 的 accepted-state-only 规则、direct-face reconstruction 和制造解 $G/G^*$ 审计。
5. 所有失败事务、checkpoint、随机数和 accepted state 所有权。

禁止添加 Poisson 投影、净电荷补偿、平均场扣除、后处理能量补丁、软化 K1 能量门或为了匹配 EPOCH 修改物理参数。

## 3. P3 的离散目标

在每个 accepted 时间步，背景 Bulk 必须同时满足：

$$
\frac{N_i^{n+1}-N_i^n}{\Delta t}
+\frac{Q_{i+1/2}-Q_{i-1/2}}{\Delta x}=0,
$$

$$
J_{\rm x,f}=qQ_f/\Delta t,
\qquad
J_{\rm x,E,i}=GJ_{\rm x,f},
$$

$$
\Delta K_{\rm bulk,E}
=\Delta t\sum_i \Delta x\,E_i J_{\rm x,E,i}
+B_u^K,
$$

$$
\Delta U_E+W_{\rm electrode}
=-\Delta t\sum_f w_fE_fJ_{\rm x,f}.
$$

这里 $Q_f$ 是 `advect_x()` 实际写入状态的最终 swept number，不能由瞬时 $q\int vf\,d^3u$、未限制高阶通量或密度差反推替代。$G$ 和 $G^*$ 必须沿用现有已验收的生产实现。

## 4. 阶段 P3.0：先证明缺口来自 x 通量，不改生产推进

### 4.0 当前实现状态（2026-08-20）

已完成以下非生产改动：

```text
新增 tests/vpfp_x_u_power_pairing_test.cpp
新增 CMake target: vpfp_x_u_power_pairing_test
生产推进方程修改: 0
legacy/discrete-gradient 物理状态修改: 0
```

该测试直接调用生产 `ConservativePpmRemap::advect_x()` 和
`ConservativePpmRemap::advect_u_parallel()`，读取前者实际写入的
`XFaceTransportAudit::bulk_number_swept_face` 与后者实际写入的
`local_delta_ke_by_x`。它不复写 PPM 或 u-remap 通量，也不修改生产状态。

本机没有可用 MPI C++ 工具链，无法完成本地链接；这不是测试失败。必须在集群构建环境执行下面命令：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target vpfp_x_u_power_pairing_test

yhrun -N 1 -n 1 --cpu-bind=cores ./build/vpfp_x_u_power_pairing_test \
  --case all \
  --result ./output/p3_x_u_pairing_unit.result
```

运行完成后，只读取 `uniform_relative`、`gradient_relative`、两个 continuity 值和 `status`。
不得因 `gradient_relative` 较大而修改生产代码；它正是 P3.0 要测量的对象。

### 4.1 新增只读算子审计

新增 `tests/vpfp_x_u_power_pairing_test.cpp` 与对应 CMake target。测试必须链接生产
`ConservativePpmRemap`、`VlasovSplitStep`、`FieldParticlePowerAudit` 和 `OpenElectrostaticSolver`，不得复写 PPM 或 u-remap 公式。

每个测试案例从生产 `advect_x()` 获取最终 `bulk_number_swept_face`，再调用生产 `advect_u_parallel()` 获取每 cell 的 `local_delta_ke_by_x`。输出：

```text
j_charge_face
g_j_charge_cell
u_force_work_cell
u_force_power_cell
delta_power_cell
sum_delta_power
continuity_residual
poisson_transport_residual
```

其中：

$$
\delta W_i=\Delta K_{u,i}-\Delta t\Delta x E_i(GJ_{\rm charge})_i.
$$

### 4.2 必须覆盖的案例

1. 单一速度格、均匀 $f$、常量正/负场。
2. 单一速度格、线性 x 梯度、常量正/负场。
3. 多速度 Maxwellian、均匀场。
4. 非均匀 x 分布且 PPM limiter 活跃。
5. 单 rank 开放双端点。
6. MPI 两 rank 共享面。

### 4.3 P3.0 放行条件

- 案例 1--3 的 $\delta W$ 必须为稳定求和舍入误差；否则先修生产 $G/G^*$ 或 u 动能矩，禁止进入 P3.1。
- 只有案例 4 的 $\delta W$ 显著且与 limiter/PPM flux 差异相关，才证明应重构 x/u 联合通量。
- 开放端点和 MPI 案例必须保持 number continuity、owner 唯一性、$dx/2$ 权重和 restart 语义。

若案例 1--3 即失败，停止并回到速度动能矩；若案例 4 不失败，则当前 K1 缺口来自未覆盖的 Tail/Beam 或时间层，不得修改 x/u Bulk 算子。

### 4.4 P3.0 最新远程结果与根因选择（2026-08-20）

证据位置：

```text
remote_cluster:output/p3_x_u_pairing_unit.result
timestamp=2026-08-20 19:16
time_layer_sequence=Tx_half_Tu_full_Tx_half
charge_current_source=Qx1_plus_Qx2_over_dt
harness_integrity_pass=1
physical_pairing_acceptance=FAIL
```

该轮已经正确使用实际 Strang 时间层和两次 x 半步的最终 swept flux；故下面的差异不是旧测试的时间层错误。

| 案例 | 相对功误差 | x1 limiter 覆盖 | x2 limiter 覆盖 | 判定 |
|---|---:|---:|---:|---|
| uniform drift, $E<0$ | $5.344\times10^{-3}$ | 0.370 | 0.419 | FAIL |
| uniform drift, $E>0$ | $5.205\times10^{-3}$ | 0.370 | 0.418 | FAIL |
| single velocity, $E<0$ | $5.392\times10^{-2}$ | 0.00243 | 0.00596 | FAIL |
| single velocity, $E>0$ | $6.456\times10^{-2}$ | 0.00243 | 0.00614 | FAIL |
| single velocity gradient, $E<0$ | $5.274\times10^{-2}$ | 0.000716 | 0.00150 | FAIL |
| single velocity gradient, $E>0$ | $6.573\times10^{-2}$ | 0.000716 | 0.00150 | FAIL |
| gradient limiter | $6.588\times10^{-3}$ | 0.142 | 0.143 | FAIL |

所有案例 continuity residual 均远低于对应稳定求和容差，且 finite/Strang/final-flux-current harness 均为 PASS。

**已确认的根因层级**：缺口在单速度、近乎不触发 limiter 的均匀态中已达 5%--7%，而加 x 梯度或显著 limiter
覆盖并没有产生同阶额外放大。因此 PPM/FCT 不是第一根因，不能进入原 P3.1 的“先改 x PPM 通量”路径。

当前生产算子的离散速度表示不一致是首要候选：

```text
x transport: cgrid_->vx[j,k] = c * u_parallel_cell / gamma_cell
u electric work: upar_swept_[j+1/2] * (K[j,k] - K[j-1,k])
```

前者使用 cell-center 的解析速度，后者使用非均匀 u 网格上 cell kinetic-energy 的有限差分。
这可能造成数个百分点的离散链式法则差异；但 P3-V.1 的首轮绝对 swept-mass 审计不是有效的 signed-work
比较，尚不能将其宣布为唯一根因。必须先通过 P3-V.1R。

结论状态：

```text
P3.0_harness=PASS
P3.0_physical_pairing=FAIL
root_cause_confirmed=single_velocity_baseline_pairing_defect
root_cause_candidate=velocity_energy_discrete_chain_rule_mismatch
root_cause_rejected_as_primary=PPM/FCT flux shaping
P3-V.1=SUPERSEDED_BY_P3-V.1R
P3-V.1R=PASS
P3-V.2=UNBLOCKED_PENDING_SELECTION
original_P3.1_x_flux_reconstruction=BLOCKED
```

### 4.5 下一步修复：P3-V 离散速度—动能链式法则统一

P3-V 是当前唯一允许的生产修复方向。它不是缩放加速度，也不是修改能量账；它要让 x transport 的
速度矩与 u 有限体积受力功由同一份离散 kinetic-energy gradient 定义。

#### P3-V.1：signed-work 速度矩审计（旧版本，已由 P3-V.1R 取代）

2026-08-20 的远程结果显示：`p3v1_velocity_audit_pass=0`。这不是 P3-V 根因被排除，
而是当前审计把无符号的“平均 u-face 速度”错误地与有符号的电场功残差比较。

当前错误实现：

```text
u_force_effective_velocity = sum(|S_u_face| * vK_face) / sum(|S_u_face|)
x_transport_effective_velocity = sum(|Q_x_face| * vx_face) / sum(|Q_x_face|)
velocity_chain_rule_pass requires raw velocity error and signed power error
to have the same sign.
```

该定义错误的原因：

1. $|S_u|$ 与 $|Q_x|$ 去除了通量方向，不能表示电场做功；
2. u 侧没有以实际 $S_u(K_R-K_L)$ 作为权重；
3. 判据遗漏背景电子 $q=-e$ 与电场符号。正确符号关系是：

$$
\operatorname{sign}(W_u-W_x)
=\operatorname{sign}\left[qE\,(v_u-v_x)\right].
$$

因此，例如正场单速度案例中 $v_u-v_x<0$ 而 $W_u-W_x>0$ 是预期的电子符号翻转，
不能判为“根因不一致”。

当前远程结果的关键矛盾也说明该审计无效：

| 单速度案例 | 功误差 | 当前速度误差 |
|---|---:|---:|
| $E<0$ | $+5.392\%$ | $+0.0619\%$ |
| $E>0$ | $+6.456\%$ | $-5.696\%$ |

正场量级接近、符号按 $qE$ 修正后正确；负场量级完全不对应。原因是当前
`u_force_effective_velocity` 不是实际受力功电流，不能用于根因选择。

##### P3-V.1R 修复要求（只改测试/审计，不改生产推进）

修改 `tests/vpfp_x_u_power_pairing_test.cpp`，删除以 `abs(swept)` 加权的两个 effective velocity
作为验收量。对每个物理 x cell、且 $|E_i|$ 高于明确的相对场阈值时，直接输出：

```text
u_work_current_cell = delta_ke_u_cell / (dt * dx * E_cell)
x_flux_current_cell = 0.5 * (J_charge_left + J_charge_right)
signed_current_difference = u_work_current_cell - x_flux_current_cell
signed_power_difference = delta_ke_u_cell - dt * dx * E_cell * x_flux_current_cell
power_current_identity_error
```

并以 u force 前的、同一物理 cell 的 Bulk number $N_{xhalf,i}$ 定义单速度辅助量：

$$
v_u^{\rm power}=
\frac{\Delta K_{u,i}}{\Delta t\,qE_iN_{xhalf,i}},
\qquad
v_x^{\rm flux}=
\frac{\Delta x\,J_{x,i}}{qN_{xhalf,i}},
\qquad N_{xhalf,i}\ [{\rm m}^{-2}].
$$

必须验证每个纳入统计 cell：

$$
\Delta K_{u,i}-\Delta t\,dx\,E_iJ_{x,i}
=\Delta t\,qE_iN_{xhalf,i}
\left(v_u^{\rm power}-v_x^{\rm flux}\right)
$$

到稳定求和误差。这里的恒等式只是审计换元，不改变生产 $f$、通量、场或 ledger。

P3-V.1R 放行条件：

1. `power_current_identity_error` 为稳定求和误差；
2. 单速度正/负场的 $v_u^{power}-v_x^{flux}$ 与实际功缺口经 $qE$ 后同号、同量级；
3. 所有排除 $|E_i|$ 过小的 cell 必须单独计数和报告，不能静默丢弃；
4. x half-step continuity、finite、开放端点和 MPI owner 语义保持原有 PASS。

若 P3-V.1R 仍不能满足第 1 项，先修测试单位、cell volume 或时间层；若满足第 1 项但第 2 项不成立，
当前根因不是离散速度表示，停止 P3-V 并重新审计 u-force 与 x-current 的时间层/物种所有权。

#### P3-V.1R 远程验收结果（2026-08-20）

集群证据：

```text
output/p3_x_u_pairing_unit.result
output/p3_x_u_pairing_unit.result.cells
output/p3_x_u_pairing_mpi_n2.result
```

单 rank 与 MPI shared-face 的 `power_current_identity_pass` 均为1。
最大逐 cell 相对恒等式误差为：

```text
single rank uniform: 2.115e-13
single rank single-velocity: 4.225e-15
MPI n=2 shared-face: 7.117e-13
```

所有案例均报告 `excluded_small_field_count=0`、
`excluded_small_number_count=0`，且包含64个物理 cell。逐 cell审计确认：
单速度正/负场中，$v_u^{\rm power}-v_x^{\rm flux}$ 经 $qE$ 后与
`signed_power_difference`同号、同量级；x half-step、finite和开放端点
harness均通过，MPI shared-face案例也通过。

因此：

```text
P3-V.1R=PASS
P3-V.2=允许进入
P3.0 physical pairing=仍FAIL（独立物理门）
原始P3.1=继续BLOCKED
```

#### P3-V.2：定义并接入能量共轭离散速度表（允许实施，未验收）

**本节当前允许实施。** P3-V.1R 的 signed-work 审计用于解释与验收 P3-V.2 的效果，
不再作为进入本节的阻断条件。P3-V.2 的目标不是修改电场力、Poisson 或能量账，而是让
`advect_x()` 使用的离散速度与 `advect_u_parallel()` 已使用的 production kinetic-energy difference
具有同一离散链式法则。

本节必须一次只实现一个新的、显式选择的 A/B 模式；不得改变默认 analytic 模式的输出。

##### P3-V.2.1 允许修改与禁止修改

允许修改：

```text
src/grid.h
src/conservative_ppm_remap.h
src/conservative_ppm_remap.cpp
src/vlasov_split_step.h
src/vlasov_split_step.cpp
src/vpfp_integrator.h
src/vpfp_integrator.cpp
src/main_vpfp.cpp
src/vpfp_checkpoint.h / src/vpfp_checkpoint.cpp
src/vpfp_diagnostics.h / src/vpfp_diagnostics.cpp
tests/vpfp_x_u_power_pairing_test.cpp
tests/vpfp_x_transport_flux_audit_test.cpp
CMakeLists.txt
```

禁止修改：

```text
OpenElectrostaticSolver 的 Poisson 算子、端点权重或 phi 边界条件
Beam/Tail PIC pusher、轨迹电流沉积或开放粒子边界
PPM/FCT 限制器公式、limiter 阈值或正性逻辑
任何全局/局部 E、J、a_u、rho、K 的事后补偿
任何默认 dt、field/pairing 容差或软接受阈值
```

##### P3-V.2.2 新枚举与默认行为

在 `src/conservative_ppm_remap.h` 定义：

```cpp
enum class XTransportVelocityMode {
    ANALYTIC_CELL_CENTER = 0,
    ENERGY_CONJUGATE_CELL = 1
};
```

要求：

1. `ANALYTIC_CELL_CENTER` 必须是默认值，并逐位保持当前 `cgrid_->vx[j,k]` 行为。
2. `ENERGY_CONJUGATE_CELL` 仅通过新的 CLI 显式启用，例如：

   ```text
   --x-transport-velocity-mode energy-conjugate
   ```

3. CLI 缺省值必须写入启动摘要、run manifest、checkpoint manifest 和 physical configuration hash。
4. restart 若模式不同必须拒绝，除非新增一个只允许**解析模式 -> 能量共轭模式的测试专用** override；
   production restart 不允许该 override。
5. 不得把新模式自动绑定到 `field-particle-coupling=discrete-gradient`；二者必须独立选择，便于 A/B。

##### P3-V.2.3 速度表的数据契约

在 `CylindricalVelocityGrid` 增加两个只读预计算表，尺寸均为 `Nv * Nmu` 或明确的 face 尺寸：

```cpp
std::vector<double> vx_energy_conjugate_face; // (Nv + 1) * Nmu
std::vector<double> vx_energy_conjugate_cell; // Nv * Nmu
```

保留现有：

```cpp
std::vector<double> vx;                       // analytic cell-center velocity
std::vector<double> kinetic_energy;           // production K[j,k]
```

内部 u face 的唯一允许定义为：

$$
v^K_{j-1/2,k}
=
\frac{K_{j,k}-K_{j-1,k}}
{m_ec\,[u_{j,k}-u_{j-1,k}]},
\qquad 1\le j<N_v.
$$

其中 $K_{j,k}$、$u_{j,k}$ 必须直接读取当前 production `kinetic_energy` 和 `upar_cells`。
不得使用解析 $cu/\gamma$ 覆盖该 face 值，不得进行平滑、裁剪或 limiter。

速度域两个物理边界 face 不参与内部 energy-current 交换。它们必须单独标记为不可用于 cell projection，
不得凭空构造周期邻居或把一端值复制到另一端。

##### P3-V.2.4 从 face 表到 x transport cell 表的唯一投影

`advect_x()` 需要每个 $(j,k)$ 的 cell speed。该 cell speed 必须通过与 u finite-volume 动能内积一致的
局部质量投影得到，不能简单取算术平均：

$$
v^K_{j,k}
=
\frac{
  w^-_{j,k}v^K_{j-1/2,k}
 +w^+_{j,k}v^K_{j+1/2,k}}
{w^-_{j,k}+w^+_{j,k}}.
$$

权重 $w^\pm$ 必须来自 production u-grid geometry 和同一 kinetic-energy difference；实现中必须：

1. 将投影封装为 `CylindricalVelocityGrid::build_energy_conjugate_velocity_table()`；
2. 在函数注释中给出 $w^\pm$ 的精确定义、量纲和端点规则；
3. 对 $j=0,N_v-1$ 使用单侧投影，且只使用域内 energy face；
4. 对对称网格/低速极限验证 $v^K_{j,k}\to c u_{j}/\gamma_{j,k}$；
5. 任何未定义、非有限或超光速值都必须在 grid 初始化时抛出明确错误，不得静默回退 analytic `vx`。

**实现注意**：如果在推导 $w^\pm$ 时无法证明上述 cell projection 与 u face 功矩一致，必须停止实现并
报告该数学缺口；禁止为了编译通过改成 $0.5(v^-+v^+)$。

##### P3-V.2.5 x transport 的最小接入点

1. 给 `ConservativePpmRemap` 增加一个 setter 或构造配置：

   ```cpp
   void set_x_transport_velocity_mode(XTransportVelocityMode mode);
   ```

2. 在 `ConservativePpmRemap::advect_x()` 的唯一速度读取点，将：

   ```cpp
   const double vx = cgrid_->vx[q];
   ```

   改为：

   ```cpp
   const double vx = transport_velocity(q);
   ```

   其中 `transport_velocity(q)` 按 mode 返回现有 `vx[q]` 或 `vx_energy_conjugate_cell[q]`。
3. 不得改动 `swept_mass()`、PPM reconstruction、halo 交换、boundary flux、FCT 或 output update 的其他语句。
4. `advect_u_parallel()` 的 $a_u=qE/(mc)$、`upar_swept_` 和 Gate-C `K_R-K_L` 保持不变；
   P3-V.2 的唯一生产物理变更是 x transport 使用的离散速度表。
5. `VlasovSplitStep`、`VpfpIntegrator` 只负责将模式传到 remap；不得在 Picard trial 内创建第二份速度表或
   按 iteration 改变 mode。

##### P3-V.2.6 诊断与 checkpoint

新增 accepted-step-only 诊断列：

```text
x_transport_velocity_mode
energy_conjugate_vx_min
energy_conjugate_vx_max
energy_conjugate_vs_analytic_linf
energy_conjugate_vs_analytic_l2
```

只在初始化和 snapshot/diagnostic interval 写入，不得在每个 inner loop 全表扫描。

checkpoint/run manifest 必须新增：

```text
x_transport_velocity_mode=analytic-cell-center|energy-conjugate-cell
x_transport_velocity_table_schema=v1
```

新的 table 不必序列化；重启时按同一 grid 与 mode 重建并验证 hash。table schema 与 mode 必须进入
physical configuration hash。

##### P3-V.2.7 分阶段测试与验收

先只编译/运行单元测试，不运行 K1 或长生产：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target \
  vpfp_x_u_power_pairing_test vpfp_x_transport_flux_audit_test fp_solver
```

1. **解析模式回归**：不传新 CLI。既有 x transport 测试必须 PASS，解析模式输出 hash/关键数值与修改前一致。
2. **energy-conjugate 表单元测试**：新增 `vpfp_energy_conjugate_velocity_test`，验证表有限、低速极限、
   对称性、端点单侧规则、无周期回卷和 $|v|<c$。
3. **P3-V.1R 重新运行**：以 energy-conjugate mode 跑单速度正/负场和 uniform drift。要求 signed
   power-current identity 成立；单速度功误差必须显著低于 analytic mode，目标为 $\le10^{-8}$。
4. **x continuity 回归**：analytic 与 energy-conjugate 两个 mode 均必须通过
   `vpfp_x_transport_flux_audit_test`；新 mode 不得改变离散 number continuity。
5. **MPI 最小回归**：

   ```bash
   yhrun -N 1 -n 2 --cpu-bind=cores ./build/vpfp_x_u_power_pairing_test \
     --case mpi-shared-face \
     --x-transport-velocity-mode energy-conjugate \
     --result ./output/p3v2_energy_conjugate_mpi_n2.result
   ```

   必须证明 shared face、开放端点和 owner 语义不变。

任何一项失败时：保留 analytic 默认模式，禁止 K1、禁止修改 PPM/FCT、禁止进入 P3.1。

##### P3-V.2.8 完成后的最小报告

执行智能体必须报告：

```text
changed_files=
default_mode_unchanged=1
analytic_regression=
energy_table_unit=
p3v1r_analytic=
p3v1r_energy_conjugate=
single_velocity_error_reduction=
x_continuity_analytic=
x_continuity_energy_conjugate=
mpi_shared_face=
checkpoint_manifest_roundtrip=
first_failure=
next_allowed_task=
```

只有 `single_velocity_error_reduction`、continuity、MPI 和 manifest 全部 PASS，
`next_allowed_task=P3-V.3`；否则为 `STOP_KEEP_ANALYTIC_DEFAULT`。

##### P3-V.2 执行命令（在 P3-V.3 之前）

以下命令只用于 P3-V.2 的单元、A/B 对照、MPI shared-face 和 checkpoint manifest
验收。命令按顺序执行；不要在这一步运行 K1 或长时间生产任务。每一组模式使用独立的
`.result` 文件，避免 analytic 与 energy-conjugate 结果互相覆盖。

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target \
  vpfp_energy_conjugate_velocity_test \
  vpfp_x_u_power_pairing_test \
  vpfp_x_transport_flux_audit_test \
  checkpoint_roundtrip_test
```

先运行速度表单元测试：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/vpfp_energy_conjugate_velocity_test \
  --result ./output/p3v2_energy_conjugate_velocity_unit.result
```

再运行 analytic-cell-center 基线。该模式是默认模式，必须先于新模式完成回归：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/vpfp_x_u_power_pairing_test \
  --case all \
  --x-transport-velocity-mode analytic-cell-center \
  --result ./output/p3v2_analytic_x_u_pairing.result

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/vpfp_x_transport_flux_audit_test \
  --case all \
  --x-transport-velocity-mode analytic-cell-center \
  --result ./output/p3v2_analytic_x_flux.result
```

确认基线后，运行 energy-conjugate A/B 对照：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/vpfp_x_u_power_pairing_test \
  --case all \
  --x-transport-velocity-mode energy-conjugate \
  --result ./output/p3v2_energy_conjugate_x_u_pairing.result

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/vpfp_x_transport_flux_audit_test \
  --case all \
  --x-transport-velocity-mode energy-conjugate \
  --result ./output/p3v2_energy_conjugate_x_flux.result
```

然后运行 MPI shared-face 最小验收：

```bash
yhrun -N 1 -n 2 --cpu-bind=cores \
  ./build/vpfp_x_u_power_pairing_test \
  --case mpi-shared-face \
  --x-transport-velocity-mode energy-conjugate \
  --result ./output/p3v2_energy_conjugate_mpi_n2.result
```

最后运行 checkpoint manifest roundtrip：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/checkpoint_roundtrip_test \
  --result ./output/p3v2_checkpoint_manifest_roundtrip.result
```

验收时至少检查：速度表有限性、端点无周期回卷、低速极限、对称性、
`|v|<c`、analytic/energy-conjugate 两组 x continuity、MPI shared-face owner
一致性，以及 checkpoint manifest 中的
`x_transport_velocity_mode` 和 `x_transport_velocity_table_schema`。任一项失败，
保持 analytic 默认模式并停止在 P3-V.2；不得进入 P3-V.3，也不得修改 PPM/FCT。

##### P3-V.2 集群验收结果（2026-08-20）

集群 `output` 根目录中的结果文件为：

```text
p3v2_analytic_x_u_pairing.result
p3v2_analytic_x_flux.result
p3v2_energy_conjugate_velocity_unit.result
p3v2_energy_conjugate_x_u_pairing.result
p3v2_energy_conjugate_x_flux.result
p3v2_energy_conjugate_mpi_n2.result
p3v2_checkpoint_manifest_roundtrip.result
```

本次结果必须区分审计门和独立物理配对门：

```text
analytic:
  status=PASS
  analytic_regression_acceptance=PASS
  harness_integrity_pass=1
  p3v1_power_identity_pass=1
  physical_pairing_acceptance=FAIL

energy-conjugate:
  velocity_table_unit=PASS
  x_flux_audit=PASS
  mpi_shared_face_harness=PASS
  checkpoint_manifest_roundtrip=PASS
  harness_integrity_pass=1
  p3v1_power_identity_pass=1
  physical_pairing_acceptance=FAIL
```

此前 energy-conjugate 的 `power_current_identity` 约 `1e-12` 假失败已由
测试审计改为 `long double` 重算；最新结果中该审计项均为 PASS。因此剩余失败不是
审计问题：完整 `all` 案例的单速度物理配对相对误差仍约 `5.2%--6.6%`，而 analytic
基线约 `5.4%--6.5%`，新速度表没有降低该误差，部分案例略有增大。

所以本阶段的正式结论为：

```text
P3-V.2 audit=PASS
analytic regression=PASS
single_velocity_error_reduction=FAIL
P3.0 physical pairing=FAIL
next_allowed_task=STOP_KEEP_ANALYTIC_DEFAULT
```

不得把 `physical_pairing_acceptance=FAIL` 改成通过，也不得通过放宽阈值、修改
PPM/FCT 或修改生产场推进来掩盖该缺口。后续必须先提出新的、可审计的联合 x/u
离散功配对修复方案；在此之前保持 `analytic-cell-center` 为默认生产模式，禁止
进入 P3-V.3 或 P3.1。

##### P3-V 最终根因结论（2026-08-20）

P3-V.1R 已证明 signed power/current 审计、单位、$qE$ 符号、cell volume、Strang 时间层和 MPI shared-face
口径正确；P3-V.2 又证明能量共轭速度表、x continuity、checkpoint manifest 和 analytic 回归正确。
但替换 x transport speed 后，单速度物理功误差没有下降。

因此以下候选已被排除为**首要根因**：

```text
Poisson/Gauss/端点权重/G-Gstar 伴随错误
Beam/Tail 边界功或 MPI shared-face owner 错误
PPM/FCT limiter 本身
cell-center analytic vx 与 energy-conjugate vx 表的单独不一致
测试单位、qE 符号、dt、cell volume 或 accepted-state 审计口径
```

当前最终确定的结构性根因是：

$$
\boxed{
\text{现有 Strang PPM }T_x(\Delta t/2)
\circ T_u(E,\Delta t)
\circ T_x(\Delta t/2)
\text{ 不构成同一个离散相空间能量系统。}
}
$$

更具体地说：

1. $J_{\rm charge}$ 来自两个独立 PPM x characteristic sweep 的最终有限体积通量；
2. $\Delta K_u$ 来自独立 u characteristic sweep 的有限体积 kinetic-energy face difference；
3. 这两个算子各自守恒、各自通过局部审计，但不存在同一个离散乘积法则/反对称 Poisson bracket，
   使得

   $$
   \Delta K_u
   =\Delta t\sum_i dx\,E_i(GJ_{\rm charge})_i
   $$

   对任意单速度状态成立；
4. 因而单速度、低 limiter 覆盖状态已经保留 5%--7% 缺口，时间步减半和单独速度表替换都不能消除它。

这不是可通过局部 correction 修复的问题。后续必须转入
`docs/VPFP_联合相空间时间中心能量闭合重构实施方案.md`，建立新的联合 x/u/Poisson 时间中心算子；
当前 P3-V、P3.1、P3.2 仅保留为历史诊断，不再继续实施。

#### P3-V.3：仅在速度表闭合后再讨论 PPM/FCT

P3-V.2 后先重跑 P3.0：

- 单速度正/负场与 uniform drift 必须降至稳定求和尺度；
- 若它们通过而 gradient limiter 仍有明显残差，才允许进入原 P3.1，处理 PPM/FCT 的最终相空间通量；
- 若单速度残差不下降，撤销 P3-V.2，不得继续改 limiter 或生产场推进。

## 5. 阶段 P3.1：建立共同的相空间通量对象

仅在 P3-V.3 已使单速度/均匀态闭合、且 limiter case 仍独立失败后实施。

### 5.1 新类型

新增 `src/vpfp_phase_space_power_flux.h/.cpp`，定义只属于一个 trial 的非持久对象：

```cpp
struct VpfpPhaseSpacePowerFlux {
    std::vector<double> x_swept_number_face;
    std::vector<double> x_charge_current_face;
    std::vector<double> x_energy_current_cell;
    std::vector<double> u_energy_target_cell;
    std::vector<double> u_energy_actual_cell;
    std::vector<double> power_defect_cell;
    bool valid;
};
```

要求：

- `x_swept_number_face` 只从 `XFaceTransportAudit::bulk_number_swept_face` 复制；
- `x_charge_current_face[f]=q*x_swept_number_face[f]/dt`；
- `x_energy_current_cell=G(x_charge_current_face)`，物理端点与 MPI shared face 复用已验收的生产 gather；
- `u_energy_target_cell=dt*dx*E_cell*x_energy_current_cell`；
- 此对象只能存在于 Picard trial，不得写入 accepted state/checkpoint，不得参与 legacy 路径。

### 5.2 时间层与重算顺序

Strang 顺序中第二个 x 半步依赖 u-force 后的分布，故最终 x flux 只能在一次 provisional trial 后得到。实现必须采用受控的内层共同通量迭代：

```text
x_half_1 (frozen)
  -> provisional u_full
  -> provisional x_half_2
  -> final accepted x swept flux Q
  -> build VpfpPhaseSpacePowerFlux
  -> rerun u_full from the same frozen x_half with the P3 flux object
  -> rerun x_half_2
  -> update Q and repeat until Q / u-work / E_pair 同时收敛
```

每次重跑都必须从 `state_x_half_`、Tail/Beam frozen trial 副本开始，不能对已推进的 `state_u_full_` 叠加第二次 kick。

P3 内层与现有 Picard 外层应共用一个有限迭代预算；每次 trial 必须固定 collective 序列。若共同通量迭代未收敛，返回新的明确失败码，不允许 soft accept。

## 6. 阶段 P3.2：受力通量的正确重构

禁止“按目标/实际功比值缩放 $a_u$”。正确实现必须在 u 有限体积通量层完成：

1. 在 `ConservativePpmRemap::advect_u_parallel()` 增加可选、只供 P3 trial 使用的
   `const VpfpPhaseSpacePowerFlux* power_flux` 参数。
2. 对每个物理 x cell，先用生产 u-face kinetic-energy difference 计算当前 electric u-flux 的离散功矩。
3. 构造一个**守恒、零净 number 的 u-face correction**，只作用于内部 u faces，满足：

   $$
   \sum_j K_j[-C_{j+1/2}+C_{j-1/2}]
   =u\_energy\_target_i-u\_energy\_actual_i.
   $$

4. correction 的选择必须是受限最小范数解：最小化 $\sum_f C_f^2/w_f$，约束零速度边界通量、零 number 变化和上式功矩。不得使用单面补偿。
5. 将 correction 与生产 u swept flux 相加后，统一经过现有 Tail-interface sink、正性检查和 Gate-C 动能矩计算；禁止在这些步骤之后修改分布或账本。
6. 若 correction 会使任意最终 u-face 流量违反已存在的正性/接口约束，P3 trial 必须拒绝并减小共同通量迭代步，而非截断 correction 或修改 accepted state。

该 correction 是离散 $O(\Delta x^p)$ 一致性项，不是物理加热项。P3.2 必须证明其 L2 范数随 x 网格加密下降；若不下降，停止并撤销 P3.2。

## 7. 验收矩阵

P3.0、P3.1、P3.2 分开提交和验收。不得在一个改动中同时替换 x remap、u remap 和 Poisson 算子。

### P3.0

```bash
cmake --build build -j4 --target vpfp_x_u_power_pairing_test
yhrun -N 1 -n 1 --cpu-bind=cores ./build/vpfp_x_u_power_pairing_test \
  --case all --result ./output/p3_x_u_pairing_unit.result
yhrun -N 1 -n 2 --cpu-bind=cores ./build/vpfp_x_u_power_pairing_test \
  --case mpi-shared-face --result ./output/p3_x_u_pairing_mpi_n2.result
```

### P3.1/P3.2 单元与短回归

```bash
cmake --build build -j4 --target fp_solver vpfp_x_u_power_pairing_test \
  vpfp_field_particle_pairing_test vpfp_field_particle_pairing_mpi_test

# 使用既有 115 fs checkpoint，重跑 K1 coarse/fine；命令中的 checkpoint 路径以当前服务器 manifest 为准。
# 只在 P3.2 单元测试 PASS 后提交，禁止直接启动长生产。
```

### K1 接受条件

1. 所有既有 JC0--JC5、Gate I、direct-face 与制造解结果继续 PASS。
2. `continuity_pass=1`、Gauss/Poisson/post-field charge PASS、无 soft acceptance。
3. `current_pair/full` 降至稳定求和尺度；不得仅降低显示用的诊断残差。
4. `energy_abs_fine_over_coarse < 0.5`，并且 coarse/fine 的残差绝对值均低于修复前。
5. P3 correction L2 与最大值随 x 网格加密下降；若不下降，判定为非一致物理补丁并撤销。
6. 通过后才允许 K2 宏观物理 A/B；K2 必须核对波包、能量分配、Beam/Tail 谱及性能。

## 8. 停止条件

出现以下任一情况必须停止 P3 物理修改，仅保留诊断：

- P3-V.1 的有效速度审计不能复现单速度功缺口，或 P3-V.2 后单速度功缺口不随修复下降；
- P3 correction 不随网格收敛；
- correction 频繁触及 u 边界、Tail interface 或正性约束；
- K1 能量账改善但宏观波形、连续性、Gauss 或 Tail/Beam ledger 变差；
- 为通过 K1 必须放宽现有物理门或加入全局补偿。

在任一停止条件下，正确结论是“当前 Strang PPM/u-remap 离散不能以局部修正获得完整能量闭合”，应转向真正的联合弱形式/DG 或结构保持半拉格朗日重构，而不是继续叠加限制器。
