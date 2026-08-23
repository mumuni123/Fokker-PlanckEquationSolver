# VPFP F10 情形 A：连续性—Poisson 功配对严格修复实施方案

## 0. 用途与最高纪律

本文档处理情形 B 修复完成后，F10/J1 剩余的唯一失败：

```text
Poisson scalar identity PASS
Poisson-current pairing FAIL
```

本文档是交给 GPT-5.6 Luna、DeepSeek-v4-Flash 等代码代理的强制执行规格。代理不得自行选择算法、改变阶段顺序或扩大修改范围。

### 0.1 不可变物理约束

1. 严格保持 `docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md` 定义的物理模型。
2. 生产场保持非周期 Gauss/Poisson、`DIRICHLET_PHI` 和默认零端点电势。
3. 生产背景 bulk 保持 reservoir/open；Beam 保持开放注入/开放流出；Tail 保持同一背景电子物种的 PIC 高能表示。
4. J1 的 periodic x 仅是隔离测试拓扑，不得进入生产开放边界。
5. 不得改变 Poisson stencil、边界条件、场能定义、Beam、Tail、conversion、return、collision 或 checkpoint 语义。
6. 不得回退或继续修改已经通过的情形 B helper `build_periodic_x_adjoint_cell_field()`。

### 0.2 禁止的“修复”

- 能量投影；
- 电流或电场缩放；
- 分布函数裁剪；
- 全局补偿、零模扣除或强制平均场；
- 放宽 `failure_code=75` 的能量门；
- 把 FAIL 改成 PASS；
- 未经 A1--A-S0--A-S1--A-FS--A3R 证明就收紧全部容差；
- 测试复写生产公式后宣称生产实现通过。

### 0.3 执行约束

1. 所有文件修改在本地工作区完成。
2. 集群编译与测试由用户执行。
3. 没有用户提供的集群结果时，阶段报告必须写“待集群执行”，不得推测 PASS/FAIL。
4. 每次只执行一个阶段；当前阶段未 PASS 立即停止。
5. 所有报告写入 `docs/情形A执行情况.md`。
6. 不得删除、覆盖或清空旧 `output`；每次使用新的时间戳目录。

### 0.3.1 阶段门与原始求解器 `status` 必须分离

本文包含多种“定位阶段”和“局部修复阶段”。这些阶段的任务是验证某个结构恒等式或消除某一层误差，而不是提前解决最终总能量 code 75。因此必须区分：

```text
raw_solver_status：测试程序对完整J1最终接受条件的结果
stage_gate：当前文档阶段自己的结构验收结果
```

固定规则：

1. A-S0、A-S1、A-FS-R1 的前置检查必须读取对应执行报告中的 `阶段门`，不得直接用原始 `.result` 的 `status` 替代。
2. `smooth-perturbed-background` 在 A4 前保留 `status=FAIL/failure_code=75` 是预期现象，因为剩余 `W_transport` 尚未修复。
3. 若某阶段文档明确允许 code 75，同时其结构门全部通过，则该阶段必须记为 PASS。
4. 只有 A5/A6/A7/A8 等最终接受回归才要求相应正测试 `status=PASS/failure_code=0`。
5. 禁止形成以下逻辑死锁：

```text
必须先消除code75才允许进入A4，
但code75恰好需要A4才能消除。
```

6. A3R 的前置条件是：

```text
A-S0 stage_gate=PASS
A-S1 stage_gate=PASS
A-FS-R1 stage_gate=PASS
```

不是要求 A-S0/A-S1 的原始 F10 `.result status=PASS`。

### 0.4 2026-08-22 当前执行状态与修订结论

```text
A0 = PASS
A1 = PASS
A2 = PASS
A3 = 命中旧分类 A-S
```

A3 的依据是：

```text
charge_projection_mismatch_linf = 3.2623855132196505e-11 C/m2
旧 tau_C = 1.8189894035458565e-12 C/m2
u_boundary_charge_linf = 0
```

该旧 `tau_C` 只包含电荷变化和 current divergence 的小差量，没有包含构造 `rho=qe*(ni-ne)` 时参与近中性相减的父电荷尺度。因此旧 A-S 报告保留为历史执行记录，但其“真实 source”解释被本修订撤销。当前唯一允许的下一阶段是 A-S0。

修订后的顺序为：

```text
A0 -> A1 -> A2 -> A3
                 |
                 +-- A-N ------------------------> A4 -> A5N
                 |
                 +-- provisional A-S -> A-S0 -> A-S1 -> A-FS -> A3R
                                                            |
                                                            +-- assembly-only -> A5A
                                                            +-- A-N ----------> A4 -> A5N
                                                            +-- physical/mixed -> STOP

A5A or A5N -> A6 -> A7 -> A8
```

### 0.5 2026-08-22 A-FS 第二轮结果后的容差决策

最新 A-FS 结果确认：

```text
normal case PASS
zero-change PASS
stable_accumulation_used=1
J0/B4/A2/A-S1 无回归
large-baseline 正负 case 因 double rho_delta 输入量化失败
nonzero endpoint case 暴露当前生产模型之外的边界功缺口
F10 两个 case 尚未完成 A-FS 回归
```

用户允许在不掩盖累计漂移的前提下适当放宽舍入门。本文据此采用以下唯一调整：

1. 稳定求和后的 Poisson scalar identity 系数由 `8192` 调为 `16384`；只放宽 2 倍。
2. A4 的 pairing solver 最终保持 `1e-9`。A6-R2 已证明修复多步 signed-state 契约后，`dt`、`dt/2` 和10步累计都能在该门下通过；不再需要放宽到 `2e-9`。
3. 总能量接受门继续保持 `1e-8`，不得放宽。
4. 增加 10 步累计有符号残差门，快速检查单步残差是否立即同号积累；100 步降级为 A8 后的可选非阻断验证。
5. large-baseline case 必须修制造解输入，不能靠扩大容差通过。
6. 非零端点电势 case 标记为当前零端点生产模型之外的已知限制，不得伪装成 PASS，也不再阻断零端点 A-FS。

---

## 1. 已确认事实

### 1.1 情形 B 已经闭合

B8 给出：

```text
W_u = 17.779409055027578 J/m2
W_F = 17.779409055027571 J/m2
W_J = 17.779409055027124 J/m2
W_u-W_F = 7.1e-15 J/m2
R_FJ/scale = 2.52e-14
R_uJ = 4.5474735088646412e-13 J/m2
```

因此：

$$
W_u=W_F=W_J
$$

已达到舍入误差。不得重新怀疑 u-face 几何、Hamiltonian velocity、midpoint mass、periodic seam 加权转置或情形 B 的 current/force 映射。

### 1.2 当前剩余量

```text
field_energy_change = -17.779411022900604 J/m2
electrode_work = 0
poisson_potential_charge_work = -17.779411022868356 J/m2
W_J = 17.779409055027124 J/m2
R_PJ = -1.9678734801686915e-06 J/m2
energy_residual = -1.9678730254213406e-06 J/m2
relative_energy_residual = 1.107e-07
failure_code = 75
```

Poisson scalar identity定义为：

$$
R_P=\Delta U_E-W_{\mathrm{electrode}}-W_{\rho\phi}
$$

其中：

$$
W_{\rho\phi}=\sum_i\bar\phi_i\left(\rho_i^{n+1}-\rho_i^n\right)\Delta x
$$

当前未闭合的是：

$$
R_{PJ}=\Delta U_E-W_{\mathrm{electrode}}+\Delta t\left\langle E^{\mathrm{pair}},J^{\mathrm{charge}}\right\rangle_f
$$

所以只允许检查最终候选态的电荷变化、连续性电流及二者的时间层/权重/所有权。

---

## 2. 情形 A 的精确分解

### 2.1 相空间残差

必须直接使用 `evaluate_local_residual()` 返回的最终 `accepted_residual`：

$$
r^M_{i,j,k}=M^{n+1}_{i,j,k}-M^n_{i,j,k}
+\Delta t\left(F^x_{i+1/2,j,k}-F^x_{i-1/2,j,k}\right)
+\Delta t\left(F^u_{i,j+1/2,k}-F^u_{i,j-1/2,k}\right)
$$

不得在诊断里重建另一套通量。

### 2.2 两种电荷残差

电子电荷为：

$$
q_e=-e
$$

相空间残差的电荷矩：

$$
r^Q_i=q_e\sum_{j,k}r^M_{i,j,k}
$$

由最终 Poisson 电荷和生产 face current 独立构造：

$$
r^C_i=\left(\rho_i^{n+1}-\rho_i^n\right)\Delta x
+\Delta t\left(J_{i+1/2}-J_{i-1/2}\right)
$$

速度边界项：

$$
r^{u,\partial}_i=q_e\Delta t\sum_k
\left(F^u_{i,N_u+1/2,k}-F^u_{i,1/2,k}\right)
$$

正确关系为：

$$
r^C_i=r^Q_i-r^{u,\partial}_i
$$

J1 的速度边界通量应为零，因此预期：

$$
r^C_i=r^Q_i
$$

### 2.3 potential-weighted defect

必须使用 `evaluate_work_identity()` 相同的 cell-average potential：

$$
\bar\phi_i^s=\phi_i^s+\frac{\Delta x}{12}
\left(E_{i+1/2}^s-E_{i-1/2}^s\right),\qquad s\in\{n,n+1\}
$$

$$
\bar\phi_i=\frac{\bar\phi_i^n+\bar\phi_i^{n+1}}{2}
$$

定义：

$$
W_C=\sum_i\bar\phi_i r^C_i
$$

严格配对必须满足：

$$
R_{PJ}=R_P+W_C
$$

预测误差：

$$
R_A^{\mathrm{pred}}=R_{PJ}-R_P-W_C
$$

---

## 3. 根因选择树

当前首选假设是：J1 在归一化相空间残差约 `1.148e-13` 时停止，但残差经电荷矩和大电势加权后仍形成约 `1.968e-6 J/m2` 的功缺口。A1--A3 通过前不得把该假设当成结论。

### 3.1 A-N：非线性残差投影

仅当下列全部成立：

```text
poisson_scalar_identity_pass=1
charge_residual_projection_pass=1
poisson_current_prediction_pass=1
u_boundary_flux_pass=1
```

并且：

$$
R_{PJ}\approx R_P+W_C
$$

才允许进入 A4。

### 3.2 A-T：时间层不一致

scalar identity 与 charge projection 通过、预测恒等式失败，且只有某个 old/candidate/final 组合闭合时，输出 `A-T_TIME_LAYER_MISMATCH` 并停止。

### 3.3 A-O：face ownership/endpoint

单 rank 通过而多 rank 失败，或误差集中在 rank 接口/物理端点时，输出 `A-O_FACE_OWNERSHIP_OR_ENDPOINT`。只允许检查 endpoint half-weight、shared-face 单次 ownership 和 MPI 重复计数。

### 3.4 A-S：待拆分的 source/u 边界或电荷装配误差

若下式不在舍入误差内：

$$
r^C_i-r^Q_i+r^{u,\partial}_i
$$

输出 `A-S_SOURCE_OR_U_BOUNDARY_OR_ASSEMBLY`，但该分类只是进入 A-S0 的中间门，不能直接解释为真实物理 source。

J1 中 Beam/Tail/collision 均关闭，不应存在外部 source。生产 J1 当前分别构造：

```text
rho_n   = qe * (ni - ne_n)
rho_np1 = qe * (ni - ne_np1)
```

当离子与电子密度均约为背景密度、净扰动远小于两者时，两次近中性大数相减后再计算 `rho_np1-rho_n` 会损失有效数字。相比之下，`rQ` 直接对 `candidate-m_old` 和相空间 residual 求矩，数值路径更稳定。因此 A-S 必须再拆成：

```text
A-S-PHYSICAL_SOURCE_OR_BOUNDARY
A-S-CHARGE_ASSEMBLY_ROUNDOFF
A-S-MIXED
```

拆分只能由 A-S0 完成；A3 不得根据旧 `tau_C` 直接判定存在真实 source。

### 3.5 A-P：Poisson scalar identity

若 `OpenPoissonWorkIdentity.residual` 自身失败，输出 `A-P_POISSON_SCALAR_IDENTITY`。只检查 fixture、phi reconstruction 和只读 identity；禁止直接修改 Poisson stencil。

### 3.6 不唯一

多个分支同时命中必须输出 `INCONCLUSIVE`。

---

## 4. A0：冻结基线

不改源码。记录 commit、dirty 状态、相关文件 hash 和 B8 两个 result 的完整路径。B8 数值必须与第 1 节一致才允许 A1。

---

## 5. A1：只读残差分解诊断

### 5.1 白名单

```text
src/vpfp_integrator.h
src/vpfp_integrator.cpp
tests/joint_phase_space_midpoint_energy_test.cpp
docs/情形A执行情况.md
```

禁止修改 Poisson、joint operator、flux、Newton、line search、energy gate、acceptance、dt 和容差。

### 5.2 新增字段

```cpp
double joint_midpoint_poisson_scalar_identity_residual;
double joint_midpoint_continuity_charge_linf;
double joint_midpoint_continuity_charge_l1;
double joint_midpoint_residual_charge_linf;
double joint_midpoint_charge_projection_mismatch_linf;
double joint_midpoint_u_boundary_charge_linf;
double joint_midpoint_potential_weighted_continuity_defect;
double joint_midpoint_poisson_current_predicted_residual;
double joint_midpoint_poisson_current_prediction_error;
double joint_midpoint_continuity_roundoff_bound;
double joint_midpoint_prediction_roundoff_bound;
int joint_midpoint_continuity_first_bad_global_ix;
```

全部字段必须初始化，且只用于诊断。

### 5.3 数据来源

在最终 `poisson_work`、`accepted_bundle`、`accepted_residual` 已生成后、code 75 判断前计算。只使用：

```text
candidate
m_old
accepted_residual
accepted_bundle.x_flux_rate
accepted_bundle.u_flux_rate
accepted_bundle.charge_current_face
fields
candidate_fields
poisson_work
```

### 5.4 实现细则

每个 x cell 的 `accepted_residual` 电荷矩使用 `long double` 累加：

```cpp
long double sum_residual = 0.0L;
for (int q = 0; q < nq; ++q) {
    sum_residual += static_cast<long double>(
        accepted_residual[ix * nq + q]);
}
const long double rQ = -Const::qe * sum_residual;
```

独立连续性残差：

```cpp
const long double delta_rho_dx =
    static_cast<long double>(candidate_fields.rho[ng + ix] -
                             fields.rho[ng + ix]) * grid_.dx;
const long double current_div_dt =
    static_cast<long double>(dt) *
    (accepted_bundle.charge_current_face[ix + 1] -
     accepted_bundle.charge_current_face[ix]);
const long double rC = delta_rho_dx + current_div_dt;
```

u 边界项必须由 `u_flux_rate` 上下端实际值计算，不得假定为零。paired potential 必须逐字使用第 2.3 节公式。

本地 `sum_abs_rC`、`sum_phi_rC`、`sum_abs_phi_rC` 使用 `long double`。所有 rank 进入同样顺序的 collective；A1 首先只跑单 rank。

### 5.5 舍入界

$$
\tau_C=8192\epsilon_{\mathrm{mach}}\max\left(
1,\max_i|\Delta\rho_i\Delta x|,
\max_i|\Delta t\Delta J_i|,\max_i|r^Q_i|\right)
$$

$$
\tau_A=8192\epsilon_{\mathrm{mach}}\max\left(
1,|R_{PJ}|,|R_P|,|W_C|,
\sum_i|\bar\phi_i r^C_i|\right)
$$

### 5.6 A1 门

```text
全部A1字段存在且finite
scalar residual与poisson_work.residual一致
first_bad_global_ix合法或为-1
旧B8字段仍存在
failure_code=75时仍保持FAIL
```

A1 只验收诊断，不要求 J1 PASS。

---

## 6. A2：生产对象直调恒等式测试

### 6.1 白名单与 case

只允许修改现有 J0/J1 测试；无法复用 target 时才最小修改 CMake。新增：

1. `a2-zero-residual`；
2. `a2-injected-residual`，使用确定性、非对称、总电荷可控的小扰动，使 `W_C` 显著非零。

必须直接调用：

```text
JointPhaseSpaceMidpointOperator::evaluate_local_residual()
OpenElectrostaticSolver::solve()
OpenElectrostaticSolver::evaluate_work_identity()
OpenElectrostaticSolver::build_potential_pairing_field()
```

### 6.2 A2 门

```text
finite=1
poisson_scalar_identity_pass=1
charge_residual_projection_pass=1
poisson_current_prediction_pass=1
```

$$
|r^C-r^Q+r^{u,\partial}|_{\infty}\le\tau_C
$$

$$
|R_A^{\mathrm{pred}}|\le\tau_A
$$

注入 case 还要求：

$$
|W_C|>100\tau_A
$$

防止平凡零值掩盖错误。

---

## 7. A3：唯一根因分类

A3 不改代码，只读取 A1/A2，输出且只输出一个：

```text
A-N_NONLINEAR_RESIDUAL_PROJECTION
A-T_TIME_LAYER_MISMATCH
A-O_FACE_OWNERSHIP_OR_ENDPOINT
A-S_SOURCE_OR_U_BOUNDARY
A-P_POISSON_SCALAR_IDENTITY
INCONCLUSIVE
```

如果 A3 输出 `A-S_SOURCE_OR_U_BOUNDARY_OR_ASSEMBLY`，下一步只能依次进入 A-S0、A-S1、A-FS、A3R；不得提前进入 A4。

只有 A3 直接输出 A-N，或 A3R 在完成 A-S0/A-S1/A-FS 后输出 A-N，才允许进入 A4。

---

## 7A. A-S0：拆分真实输运残差与近中性电荷装配误差

### 7A.1 目标

把旧诊断中的：

$$
r^C_i-r^Q_i+r^{u,\partial}_i
$$

拆成两个独立对象：

1. 质量差与生产 x-current 的真实输运残差；
2. 从两份近中性 `rho` 状态相减产生的电荷装配误差。

本阶段只增加只读诊断，不改变候选态、Poisson、flux、Newton 或接受逻辑。

### 7A.2 白名单

```text
src/vpfp_integrator.h
src/vpfp_integrator.cpp
tests/joint_phase_space_midpoint_energy_test.cpp
docs/情形A执行情况.md
```

禁止修改：

```text
src/maxwell.cpp
src/open_electrostatic_solver.cpp
src/joint_phase_space_midpoint.cpp
任何rho装配公式
任何solver或阈值
```

### 7A.3 三个基本量

必须直接用 `m_old` 与最终 `candidate` 逐速度单元做差后再用 `long double` 累加：

$$
\Delta Q_i^{M}
=q_e\sum_{j,k}\left(M^{n+1}_{i,j,k}-M^n_{i,j,k}\right)
$$

禁止先分别求两个大总和再相减。实现形式必须是：

```cpp
long double delta_number = 0.0L;
long double sum_abs_old = 0.0L;
long double sum_abs_new = 0.0L;
for (int q = 0; q < nq; ++q) {
    const long double old_m = m_old[ix * nq + q];
    const long double new_m = candidate[ix * nq + q];
    delta_number += new_m - old_m;
    sum_abs_old += std::fabs(old_m);
    sum_abs_new += std::fabs(new_m);
}
const long double delta_q_mass = -Const::qe * delta_number;
```

由两个 Poisson 电荷态得到：

$$
\Delta Q_i^{\rho}
=\left(\rho_i^{n+1}-\rho_i^n\right)\Delta x
$$

定义装配误差：

$$
r_i^{\mathrm{assembly}}
=\Delta Q_i^{\rho}-\Delta Q_i^M
$$

定义不经过两次近中性 `rho` 相减的输运残差：

$$
r_i^{\mathrm{transport}}
=\Delta Q_i^M
+\Delta t\left(J_{i+1/2}-J_{i-1/2}\right)
$$

正确分解为：

$$
r_i^C=r_i^{\mathrm{assembly}}+r_i^{\mathrm{transport}}
$$

同时应有：

$$
r_i^{\mathrm{transport}}=r_i^Q-r_i^{u,\partial}
$$

### 7A.4 父操作数与舍入界

旧 `tau_C` 只按小差量定标，不能判断近中性装配误差。新增每个 cell 的父量：

$$
S_i^{\mathrm{parent}}
=|q_e n_i\Delta x|
+|q_e|\sum_{j,k}|M^n_{i,j,k}|
+|q_e|\sum_{j,k}|M^{n+1}_{i,j,k}|
$$

定义：

$$
\gamma_m=\frac{m\epsilon_{\mathrm{mach}}}{1-m\epsilon_{\mathrm{mach}}},
\qquad m=N_uN_{\perp}+8
$$

若分母非正或非有限，诊断直接失败。装配舍入界为：

$$
\tau_i^{\mathrm{assembly}}
=32\gamma_m\max\left(1,S_i^{\mathrm{parent}}\right)
$$

全局门取：

$$
\tau_{\mathrm{assembly}}
=\max_i\tau_i^{\mathrm{assembly}}
$$

该门只用于识别浮点装配误差，不得成为物理连续性接受阈值。

输运残差仍使用小差量尺度，但其舍入界必须包含逐速度单元差分的绝对和：

$$
S_i^{\mathrm{transport}}
=|q_e|\sum_{j,k}|M^{n+1}_{i,j,k}-M^n_{i,j,k}|
+|\Delta tJ_{i+1/2}|+|\Delta tJ_{i-1/2}|
$$

$$
\tau_{\mathrm{transport}}
=8192\epsilon_{\mathrm{mach}}
\max\left(1,\max_i S_i^{\mathrm{transport}}\right)
$$

### 7A.5 新增字段

```cpp
double joint_midpoint_density_assembly_mismatch_linf;
double joint_midpoint_density_assembly_mismatch_l1;
double joint_midpoint_density_assembly_roundoff_bound;
double joint_midpoint_mass_transport_charge_linf;
double joint_midpoint_mass_transport_roundoff_bound;
double joint_midpoint_transport_projection_mismatch_linf;
double joint_midpoint_parent_charge_scale_max;
double joint_midpoint_mass_delta_charge_linf;
int joint_midpoint_density_assembly_first_bad_global_ix;
int joint_midpoint_mass_transport_first_bad_global_ix;
```

还必须输出 potential-weighted 分解：

$$
W_{\mathrm{assembly}}=\sum_i\bar\phi_i r_i^{\mathrm{assembly}}
$$

$$
W_{\mathrm{transport}}=\sum_i\bar\phi_i r_i^{\mathrm{transport}}
$$

新增字段：

```cpp
double joint_midpoint_potential_weighted_assembly_defect;
double joint_midpoint_potential_weighted_transport_defect;
double joint_midpoint_weighted_defect_reconstruction_error;
```

并验证：

$$
W_C=W_{\mathrm{assembly}}+W_{\mathrm{transport}}
$$

### 7A.6 A-S0 分类门

A-S0 是诊断分类阶段。只要以下分类条件形成唯一结论，其 `stage_gate=PASS`；原始 F10 保留 code 75 不影响 A-S0 阶段门。

按以下固定规则分类：

#### A-S-CHARGE_ASSEMBLY_ROUNDOFF

```text
abs(assembly_mismatch_linf) <= assembly_roundoff_bound
mass_transport_charge_linf <= mass_transport_roundoff_bound
transport_projection_mismatch_linf <= mass_transport_roundoff_bound
u_boundary_charge_linf <= mass_transport_roundoff_bound
weighted_defect_reconstruction_error <= tau_A
```

并且：

```text
abs(W_assembly) >= 0.9 * abs(W_C)
```

#### A-S-PHYSICAL_SOURCE_OR_BOUNDARY

```text
assembly_mismatch_linf <= assembly_roundoff_bound
mass_transport_charge_linf > mass_transport_roundoff_bound
```

或 u 边界项真实超门。

#### A-S-MIXED

装配项和输运项都宏观可见，或分解重构自身不闭合。

不得仅凭 `first_bad_global_ix=0` 判定边界问题；必须比较完整空间分布和上述两个独立门。

---

## 7B. A-S1：使用数学等价的稳定增量电荷装配

本阶段仅当 A-S0 分类为 `A-S-CHARGE_ASSEMBLY_ROUNDOFF` 时执行。

### 7B.1 修复目标

J1 candidate 的离子密度固定，Beam/Tail/source 关闭，因此以下两式数学等价：

$$
\rho_i^{n+1}=q_e\left(n_i-n_{e,i}^{n+1}\right)
$$

$$
\rho_i^{n+1}=\rho_i^n
+\frac{q_e}{\Delta x}
\sum_{j,k}\left(M^{n+1}_{i,j,k}-M^n_{i,j,k}\right)
$$

第二式避免分别计算两个近中性大数再相减。该修改是稳定的代数等价实现，不是 source、能量补丁或状态修正。

### 7B.2 白名单

```text
src/vpfp_integrator.cpp
src/vpfp_integrator.h（仅新增必要诊断字段时）
tests/joint_phase_space_midpoint_energy_test.cpp
tests/joint_phase_space_midpoint_unit_test.cpp
docs/情形A执行情况.md
```

禁止修改 `EMFields::set_charge_density()` 的生产通用路径。本阶段只修改 J1 candidate evaluate，因为 J2 的 open/reservoir、Beam、Tail/source 需要单独推导，不能套用 background-only 增量式。

### 7B.3 candidate rho 的唯一实现

替换 J1 candidate evaluate 中：

```cpp
number = sum(state);
eval_fields.rho = qe * (ni - number / dx);
```

为：

```cpp
long double delta_number = 0.0L;
for (int q = 0; q < nq; ++q) {
    const size_t id = static_cast<size_t>(ix) * nq + q;
    delta_number += static_cast<long double>(state[id]) -
                    static_cast<long double>(m_old[id]);
}
const long double delta_rho =
    -static_cast<long double>(Const::qe) * delta_number /
    static_cast<long double>(grid_.dx);
eval_fields.rho[grid_.nghost + ix] =
    fields.rho[grid_.nghost + ix] + static_cast<double>(delta_rho);
```

符号必须与电子电荷一致。禁止从 current divergence 构造或覆盖 `rho`；rho 仍来自 candidate mass，只是采用稳定增量形式。

### 7B.4 同一 candidate 的一致性检查

只读计算 absolute form：

$$
\rho_i^{\mathrm{absolute}}=q_e\left(n_i-n_{e,i}^{n+1}\right)
$$

输出：

```text
candidate_rho_incremental
candidate_rho_absolute
candidate_rho_form_difference
candidate_rho_form_roundoff_bound
```

该差异只作诊断，不能用于修正增量结果。

### 7B.5 A-S1 门

重跑 A1/B8 case 后要求：

```text
density_assembly_mismatch_linf <= mass_transport_roundoff_bound
charge_projection_mismatch_linf <= mass_transport_roundoff_bound
weighted_defect_reconstruction_error <= tau_A
poisson scalar identity字段finite且预测分解闭合；其严格PASS由后续A-FS-R1负责
R_uJ保持roundoff
B4/J0全部无回归
```

A-S1 是电荷装配局部修复阶段。装配门通过但最终 F10 仍因 `W_transport` 返回 code 75 时，A-S1 必须记为 `stage_gate=PASS`。不得要求 A-S1 在 A4 前消除最终 code 75。

然后根据新的 `R_PJ` 决定：

- A-S1 完成后必须先执行 A-FS；禁止直接进入 A3R。
- A-FS PASS 后，若 F10 已满足现有 `1e-8` energy gate，进入 A3R，分类为 assembly-only resolved；不得为了得到更漂亮的账强制执行 A4。
- A-FS PASS 后，若 `R_PJ` 仍超能量门且剩余量由 `W_transport` 解释，进入 A3R 的 A-N 分支。
- 若增量/绝对形式差超过理论舍入界，A-S1 FAIL，禁止继续。

---

## 7C. A-FS：Poisson scalar identity 稳定求和

本阶段在 A-S1 完成后执行。当前严格 Gate F 数据为：

```text
poisson_identity_scale = 17.779361198365223 J/m2
poisson_identity_roundoff_bound = 3.23404696216407e-11 J/m2
abs(poisson_scalar_identity_residual) = 3.539213366821059e-11 J/m2
residual_to_bound_ratio = 1.0943605359560971
poisson_scalar_identity_pass = 0
```

同时：

```text
sum_abs_energy_before = 139160.18982782625 J/m2
sum_abs_energy_after = 139142.41046662789 J/m2
sum_abs_potential_charge = 26.669088288436217 J/m2
```

当前 `field_energy_change` 由两个约 `1.39e5 J/m2` 的独立 `double` 总和相减得到约 `17.78 J/m2`，存在明显抵消放大。A-FS 只修正 `evaluate_work_identity()` 的数值求和，不修改 Poisson 方程或场状态。

### 7C.1 唯一允许的选择与修订后 Gate F

执行稳定求和，并仅采用用户批准的 2 倍舍入裕度：

$$
\tau_P=16384\epsilon_{\mathrm{mach}}S_P
$$

其中：

$$
S_P=\max\left(
\mathrm{DBL\_MIN},
|\Delta U_E|,
|W_{\mathrm{electrode}}|,
|W_{\rho\phi}|
\right)
$$

不得把 `sum_abs_energy_before/after` 加入 Gate F 尺度；这些量只用于解释舍入来源。不得把 `16384` 再扩大。历史结果仍同时输出旧 `8192` 门的 ratio，便于比较，但阶段 PASS 使用本节新的 `16384` 门。

### 7C.2 白名单

```text
src/open_electrostatic_solver.h
src/open_electrostatic_solver.cpp
tests/vpfp_poisson_work_identity_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp（只允许输出新增诊断）
docs/情形A执行情况.md
```

只有新增测试 target 无法复用现有 `vpfp_poisson_work_identity_test` 时，才允许最小修改 `CMakeLists.txt`。优先扩展现有 target。

禁止修改：

```text
OpenElectrostaticSolver::solve()
reconstruct_phi()
build_potential_pairing_field()
Poisson stencil和边界
EMFields中的任何场值
J1 flux、Newton、G*、energy gate、dt和容差
```

### 7C.3 稳定累加器

在 `open_electrostatic_solver.cpp` 的匿名 namespace 增加只供本文件使用的 `LongDoubleNeumaierSum`：

```cpp
struct LongDoubleNeumaierSum {
    long double sum;
    long double correction;

    LongDoubleNeumaierSum() : sum(0.0L), correction(0.0L) {}

    void add(long double value) {
        const long double trial = sum + value;
        if (std::fabs(sum) >= std::fabs(value))
            correction += (sum - trial) + value;
        else
            correction += (value - trial) + sum;
        sum = trial;
    }

    long double value() const { return sum + correction; }
};
```

不得创建全局状态，不得改变 OpenMP 或 MPI 顺序。

### 7C.4 每 cell 场能与直接差分

对每个 cell，先把所有 face 值转换为 `long double`：

```cpp
const long double old_l = before.Ex_face[ix];
const long double old_r = before.Ex_face[ix + 1];
const long double new_l = after.Ex_face[ix];
const long double new_r = after.Ex_face[ix + 1];
```

场能 cell 项保持原离散定义：

$$
U_{E,i}
=\frac{\varepsilon_0\Delta x}{6}
\left(E_L^2+E_LE_R+E_R^2\right)
$$

`energy_before_sum` 和 `energy_after_sum` 分别用 `LongDoubleNeumaierSum` 累加，仅用于输出总场能。

关键的 `field_energy_change` 不得再由两个全局大总量相减得到。必须逐 cell 直接累加差量，并使用因式分解降低局部抵消：

$$
\Delta U_{E,i}
=\frac{\varepsilon_0\Delta x}{6}
\left[
(E_L^{n+1}-E_L^n)(E_L^{n+1}+E_L^n)
+(E_L^{n+1}-E_L^n)E_R^{n+1}
+E_L^n(E_R^{n+1}-E_R^n)
+(E_R^{n+1}-E_R^n)(E_R^{n+1}+E_R^n)
\right]
$$

用独立的 `field_delta_sum` 累加该项。

### 7C.5 potential-charge work

保留原 cell-average potential 公式，但所有中间运算使用 `long double`，并用独立的 `potential_charge_sum` 稳定累加：

$$
W_{\rho\phi}
=\sum_i
\frac{\bar\phi_i^n+\bar\phi_i^{n+1}}{2}
\Delta\rho_i\Delta x
$$

不得改变 `rho_delta`、phi reconstruction 或 pairing field。

### 7C.6 MPI 归约

本地至少归约以下五项：

```text
field_energy_before
field_energy_after
field_energy_change_direct
potential_charge_work
各项absolute-term-sum诊断
```

数值 identity 量必须使用 `MPI_LONG_DOUBLE` 和 `MPI_SUM`。所有 rank 必须执行相同数量、相同顺序的 collective。归约后再转换为 `double` 写入 `OpenPoissonWorkIdentity`。

不得先转成 `double` 再归约。

### 7C.7 结果字段语义

在 `OpenPoissonWorkIdentity` 中增加并在构造函数初始化：

```cpp
double field_energy_change_direct;
double field_energy_change_from_totals;
double field_energy_change_reconstruction_error;
double term_abs_sum_energy_before;
double term_abs_sum_energy_after;
double term_abs_sum_potential_charge;
bool stable_accumulation_used;
```

默认值全部为零，`stable_accumulation_used` 默认 `false`。只有完成全部 stable local sum 和 `MPI_LONG_DOUBLE` 归约后才置 `true`；任何输入尺寸或 finite 检查失败时保持 `false`。

```cpp
result.field_energy_before = static_cast<double>(global_before);
result.field_energy_after = static_cast<double>(global_after);
result.field_energy_change = static_cast<double>(global_delta_direct);
result.potential_charge_work = static_cast<double>(global_potential_work);
result.residual = result.field_energy_change -
                  result.electrode_work -
                  result.potential_charge_work;
```

字段输出语义：

```text
field_energy_change_direct
field_energy_change_from_totals
field_energy_change_reconstruction_error
term_abs_sum_energy_before
term_abs_sum_energy_after
term_abs_sum_potential_charge
stable_accumulation_used=1
poisson_identity_roundoff_bound_8192
poisson_identity_roundoff_bound_16384
poisson_identity_residual_to_bound_ratio_8192
poisson_identity_residual_to_bound_ratio_16384
poisson_scalar_identity_pass（按16384门）
```

兼容说明：已有单值 `poisson_identity_roundoff_bound` 和 `poisson_identity_residual_to_bound_ratio` 如继续保留，必须明确映射到新 `16384` 生产门；禁止同名字段在不同测试中使用不同系数。

其中：

$$
R_{\mathrm{reconstruct}}
=\Delta U_E^{\mathrm{direct}}
-\left(U_E^{n+1}-U_E^n\right)
$$

该 reconstruction error 只诊断抵消，不得进入 Gate F 或修正场能。

### 7C.8 制造解测试

扩展 `vpfp_poisson_work_identity_test`，必须至少覆盖：

1. 现有普通制造解；
2. `large-baseline-small-delta-positive`；
3. `large-baseline-small-delta-negative`；
4. 零变化场；
5. 非零 `DIRICHLET_PHI` endpoint 的边界功符号回归。

大基线 case 必须让：

```text
field_energy_before >= 1e4 * abs(field_energy_change)
field_energy_after  >= 1e4 * abs(field_energy_change)
```

同时必须满足输入可解析性：

```text
abs(expected_potential_charge_work) >= 128 * rho_delta_input_quantization_bound
abs(expected_field_energy_change) >= 128 * field_increment_quantization_bound
```

制造解的 `rho_delta` 不得由两个大 `double rho` 数组直接相减得到。必须从制造的场增量通过同一个离散 Gauss 差分构造：

$$
\Delta\rho_i
=\varepsilon_0
\frac{
\left(E_{i+1/2}^{n+1}-E_{i+1/2}^n\right)
-\left(E_{i-1/2}^{n+1}-E_{i-1/2}^n\right)
}{\Delta x}
$$

全部中间量使用 `long double`，最后一次转换为 API 所需的 `double`。这属于制造解输入修复，不得修改生产 solver。

并使用生产 `evaluate_work_identity()`，不得在测试中复制稳定求和公式。

测试调用前后必须验证 `before`、`after` 的 `Ex_face`、`Ex`、`phi`、`rho` 逐位不变，证明该函数只读。

非零 endpoint case 必须继续运行并输出：

```text
endpoint_nonzero_identity_supported=0
endpoint_nonzero_known_limitation=1
endpoint_nonzero_residual
endpoint_nonzero_residual_to_bound_ratio
```

当前不可变生产配置为 `phi_left=phi_right=0`，所以该 case 不再进入零端点生产门的总体 `status`。但不得把 `endpoint_boundary_work_pass` 强制写成 1；它必须保持真实 FAIL，并作为未来支持非零端点时的独立阻断项。

### 7C.9 A-FS 阶段门

零端点生产范围内以下全部满足才 PASS：

```text
所有制造解finite=1
stable_accumulation_used=1
poisson_scalar_identity_pass=1
poisson_identity_residual_to_bound_ratio_16384<=1
large-baseline正负case均PASS
零变化case严格为0或舍入量级
identity调用前后物理状态逐位一致
原Gate F普通case无回归
endpoint_nonzero_known_limitation=1且真实失败值完整输出
```

然后重跑 A-S1 的 J0/B4 与两个 F10 case，并要求：

```text
J0/B4全部PASS
smooth-background PASS
smooth-perturbed-background 的 poisson_scalar_identity_pass=1
R_uJ保持roundoff
assembly/transport/weighted decomposition保持闭合
```

稳定求和并采用 `16384` 门后，若零端点 case 仍有 `ratio>1`，A-FS FAIL。不得再扩大系数，不得修改 Poisson；必须输出首个失败 case 和全部 absolute-term-sum。

### 7C.10 A-FS-R1 当前唯一下一步

第二轮已经证明稳定求和代码被调用，但测试夹具和测试范围仍不符合当前生产模型。A-FS-R1 不是再次重写 Poisson identity，也不是继续调参；它只完成三件事：修复不可解析的制造解输入、隔离非零端点已知限制、补齐零端点 F10 回归。

#### 7C.10.1 A-FS-R1 白名单

默认只允许修改：

```text
tests/vpfp_poisson_work_identity_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp
docs/情形A执行情况.md
```

若现有结果结构还没有两套 Gate F 字段，才允许最小修改：

```text
src/open_electrostatic_solver.h
src/open_electrostatic_solver.cpp
src/vpfp_integrator.h
src/vpfp_integrator.cpp
```

这些源码修改只能增加或转发诊断字段，不得再次改变 7C.3--7C.7 已实现的稳定求和公式、Poisson 求解、场状态或 residual 定义。

禁止修改：

```text
OpenElectrostaticSolver::solve()
OpenElectrostaticSolver::reconstruct_phi()
OpenElectrostaticSolver::build_potential_pairing_field()
boundary_energy_work()
Poisson stencil和边界条件
J1 flux、G*、Newton、line search、energy gate、dt
任何Beam/Tail/collision/reservoir代码
```

#### 7C.10.2 修复 large-baseline fixture，而非修改生产 identity

两个 large-baseline case 的 `before` 与 `after` face field 必须由确定性解析模式生成。推荐形式：

$$
E_f^n=A\cos\left(2\pi x_f/L\right)
$$

$$
\delta E_f=\eta A
\left[
0.7\cos\left(2\pi x_f/L\right)
+0.3\sin\left(2\pi x_f/L\right)
\right]
$$

$$
E_f^{n+1}=E_f^n\pm\delta E_f
$$

其中正负 case 分别使用正号和负号。两个模式在完整离散域上的 trapezoid 均值必须达到舍入误差内的零，以满足零端点电势约束；测试必须显式输出：

```text
before_face_integral
after_face_integral
before_face_integral_tolerance
after_face_integral_tolerance
zero_endpoint_compatibility_pass
```

如果解析模式在当前离散网格上不能满足零积分，测试必须在构造阶段减去按同一 trapezoid face 权重计算的常数均值；该操作只用于制造满足零端点边界的测试场，不得进入生产代码。

#### 7C.10.3 `rho_delta` 的唯一构造方式

不得执行：

```cpp
rho_delta = after.rho - before.rho;
```

必须先在 `long double` 中形成 face 增量：

```cpp
const long double delta_e_left =
    static_cast<long double>(after.Ex_face[ix]) -
    static_cast<long double>(before.Ex_face[ix]);
const long double delta_e_right =
    static_cast<long double>(after.Ex_face[ix + 1]) -
    static_cast<long double>(before.Ex_face[ix + 1]);
const long double delta_rho =
    static_cast<long double>(Const::eps0) *
    (delta_e_right - delta_e_left) /
    static_cast<long double>(grid.dx);
rho_delta[ng + ix] = static_cast<double>(delta_rho);
```

`before.rho` 与 `after.rho` 仍可为状态一致性诊断而分别由 face Gauss 差分构造，但它们不得用于传给 `evaluate_work_identity()` 的 `rho_delta`。

#### 7C.10.4 输入量化界

对每个转换为 `double` 的 `rho_delta`，用 `std::nextafter()` 计算一半 ULP：

```cpp
const double x = rho_delta[ng + ix];
const double up = std::nextafter(
    x, std::numeric_limits<double>::infinity());
const double down = std::nextafter(
    x, -std::numeric_limits<double>::infinity());
const double half_ulp = 0.5 * std::max(up - x, x - down);
```

potential-work 输入量化界为：

$$
B_{\rho,\mathrm{quant}}
=\sum_i|\bar\phi_i|\frac{\mathrm{ULP}(\Delta\rho_i)}{2}\Delta x
$$

场增量量化界按相同原则对 `delta_e_left/right` 构造，并传播到逐 cell 场能差；输出：

```text
rho_delta_input_quantization_bound
field_increment_quantization_bound
expected_potential_charge_work
expected_field_energy_change
rho_delta_resolvable_pass
field_delta_resolvable_pass
```

必须满足：

$$
|W_{\rho\phi}^{\mathrm{expected}}|
\ge128B_{\rho,\mathrm{quant}}
$$

$$
|\Delta U_E^{\mathrm{expected}}|
\ge128B_{E,\mathrm{quant}}
$$

#### 7C.10.5 确定性选择扰动幅度

禁止运行时无限循环调参。测试只允许按顺序尝试以下固定列表：

```text
eta_candidates = {1e-4, 3e-5, 1e-5, 3e-6, 1e-6}
```

选择第一个同时满足以下条件的候选：

```text
field_energy_before >= 1e4*abs(field_energy_change)
field_energy_after  >= 1e4*abs(field_energy_change)
rho_delta_resolvable_pass=1
field_delta_resolvable_pass=1
zero_endpoint_compatibility_pass=1
```

如果没有候选满足，fixture 构造 FAIL，输出每个 eta 的四个比例并停止；不得扩大 Gate F 容差或降低 `128` 和 `1e4`。

#### 7C.10.6 两套 Gate F 字段

每个 zero-endpoint case 必须输出：

```text
poisson_identity_scale
poisson_identity_roundoff_bound_8192
poisson_identity_roundoff_bound_16384
poisson_identity_residual_to_bound_ratio_8192
poisson_identity_residual_to_bound_ratio_16384
poisson_scalar_identity_pass_8192
poisson_scalar_identity_pass_16384
poisson_scalar_identity_pass
```

兼容字段 `poisson_scalar_identity_pass` 必须等于 `_16384`。测试总体零端点生产门只使用 `_16384`，但旧 `_8192` 必须真实输出，不得覆盖。

#### 7C.10.7 非零 endpoint case 的隔离语义

该 case 必须继续运行，真实计算和输出：

```text
endpoint_nonzero_identity_pass
endpoint_nonzero_residual
endpoint_nonzero_roundoff_bound
endpoint_nonzero_residual_to_bound_ratio
endpoint_nonzero_known_limitation
endpoint_nonzero_in_production_scope
```

固定语义：

```text
endpoint_nonzero_in_production_scope=0
endpoint_nonzero_known_limitation=!endpoint_nonzero_identity_pass
```

禁止把 `endpoint_nonzero_identity_pass` 强制为 1。零端点总体状态计算必须显式排除该布尔量，但要求上述字段完整、finite。若未来生产配置允许非零端点电势，必须另立边界功方案，不能沿用本例外。

#### 7C.10.8 总体状态聚合

新增明确的聚合字段：

```text
zero_endpoint_case_normal_pass
zero_endpoint_case_zero_change_pass
zero_endpoint_case_large_positive_pass
zero_endpoint_case_large_negative_pass
zero_endpoint_production_gate_pass
nonzero_endpoint_diagnostic_complete
status
```

定义：

```cpp
zero_endpoint_production_gate_pass =
    normal_pass && zero_change_pass &&
    large_positive_pass && large_negative_pass &&
    state_unchanged_pass;

status = zero_endpoint_production_gate_pass &&
         nonzero_endpoint_diagnostic_complete;
```

其中 `nonzero_endpoint_diagnostic_complete` 只要求诊断完整、finite 和真实布尔值，不要求其 identity PASS。

#### 7C.10.9 F10 回归的特殊验收语义

必须补跑两个 F10 case。`smooth-background` 必须正常 PASS。

对 `smooth-perturbed-background`，A-FS-R1 只验收 Poisson scalar identity 与旧离散链未回归。由于剩余 `W_transport` 尚未通过 A4 修复，该 case 仍可能：

```text
status=FAIL
failure_code=75
failure_stage=joint_midpoint_energy_residual
```

这不自动导致 A-FS-R1 FAIL。A-FS-R1 对该 case 的门是：

```text
finite=1
gauss_ok=1
converged=1
poisson_scalar_identity_pass_16384=1
poisson_identity_residual_to_bound_ratio_16384<=1
R_uJ保持roundoff
assembly mismatch保持roundoff
weighted decomposition闭合到tau_A
W_transport仍完整解释剩余R_PJ
```

禁止为了让 F10 总体 `status=PASS` 提前实现 A4 或放宽总能量门。

#### 7C.10.10 A-FS-R1 最终阶段门

全部满足才 PASS：

1. zero-endpoint 四个制造解均 PASS；
2. large-baseline 两 case 同时满足 `1e4` 大基线比和 `128` 输入可解析性门；
3. 两套 Gate F 字段完整，生产 alias 使用 `16384`；
4. nonzero endpoint 真实失败值完整输出并明确 out-of-scope；
5. identity 调用前后状态逐位不变；
6. J0/B4/A2/A-S1 单元门无回归；
7. F10 smooth-background PASS；
8. F10 smooth-perturbed 满足 7C.10.9 的结构门，即使仍为 code 75；
9. 所有 `.result/.out/.err` 完整保存。

A-FS-R1 PASS 后只允许进入 A3R，不得自行进入 A4。

这里的 `A-FS-R1 PASS` 指 7C.10.10 的结构门 PASS；不要求 `smooth-perturbed-background` 的原始总体 `status=PASS`。

---

## 7D. A3R：修复装配与求和后重新分类

A3R 不改代码，读取 A-S0/A-S1/A-FS-R1 的阶段报告和底层结果。前置判断必须使用三个阶段报告中的 `stage_gate`，不得使用原始 F10 `status` 代替。A-FS-R1 结构门未 PASS 时禁止执行 A3R。输出一个且只输出一个：

```text
A-ASSEMBLY_ONLY_RESOLVED
A-N_NONLINEAR_RESIDUAL_PROJECTION
A-S-PHYSICAL_SOURCE_OR_BOUNDARY
A-S-MIXED
INCONCLUSIVE
```

机械规则：

1. A-S1 后 F10 PASS，`R_PJ` 满足现有 energy gate，scalar identity、`R_uJ`、J0/B4 均通过：输出 `A-ASSEMBLY_ONLY_RESOLVED`，跳过 A4，进入不含 pairing-negative hook 的 A5A 回归。
2. 装配误差已降至新 roundoff bound，但 `W_transport` 仍使 `R_PJ` 超门，且预测恒等式闭合：输出 A-N，进入 A4。
3. 真实输运/source/u 边界项超门：输出 A-S-PHYSICAL，停止并另写 source 方案。
4. 不满足唯一规则：输出 `INCONCLUSIVE`。

### 7D.1 当前结果的强制分类

当前已记录证据为：

```text
A-S0 stage_gate=PASS
  unique_classification=A-S-CHARGE_ASSEMBLY_ROUNDOFF
A-S1 stage_gate=PASS
  assembly mismatch已降至roundoff
A-FS-R1 stage_gate=PASS
  zero-endpoint scalar identity已通过
J0/B4/A2=PASS
R_uJ=roundoff
W_assembly可忽略
W_transport完整解释剩余R_PJ
mass_transport_charge_linf > mass_transport_roundoff_bound
smooth-perturbed raw status=FAIL/code75
```

最后一行是 A-N 尚未修复的结果，不是前置门失败。根据本节机械规则，唯一分类必须为：

```text
A-N_NONLINEAR_RESIDUAL_PROJECTION
```

下一步允许进入 A4。任何执行代理若因 A-S0/A-S1 原始 code 75 拒绝 A3R，属于阶段门语义错误，必须纠正报告，不得修改代码。

---

## 8. A4：把 Poisson-current metric 纳入收敛门

仅当 A3 直接为 A-N，或 A3R 为 A-N 时执行。若 A3R 为 `A-ASSEMBLY_ONLY_RESOLVED`，明确跳过 A4，禁止为追求更小能量数字而增加 solver 约束。

### 8.1 白名单与原则

只修改 `vpfp_integrator.h/.cpp`、J1 energy test 和执行报告。不改变联合中点方程，只改变何时宣布求解充分。禁止修改 Poisson、flux、G*、现有 `1e-8` energy gate、`max_iterations=20` 和物理模块。

### 8.2 candidate metric

candidate evaluate 返回：

```text
candidate_poisson_current_residual
candidate_poisson_current_scale
candidate_poisson_current_relative
candidate_poisson_scalar_residual
candidate_weighted_continuity_defect
```

必须来自同一 `state`、`eval_fields`、`bundle.charge_current_face` 和时间层 n 的 `fields`。

### 8.3 三重收敛门

```cpp
const double pairing_tolerance = 1.0e-9;
const bool phase_converged = accepted_norm <= residual_tolerance;
const bool poisson_converged =
    accepted_poisson_residual <= poisson_tolerance;
const bool pairing_converged =
    accepted_poisson_current_relative <= pairing_tolerance;
```

三个布尔量全真才收敛。`1e-9` 为 A6-R2 验证后的最终值，不得再次放宽。

### 8.4 line search

phase 未通过时保留现有逻辑。phase 已通过而 pairing 未通过时，仅允许：

```cpp
trial_norm <= residual_tolerance &&
trial_poisson <= poisson_tolerance &&
trial_poisson_current_relative < accepted_poisson_current_relative
```

两个 line-search 路径必须一致。不得接受违反 phase 或 Poisson 门的状态。

### 8.5 日志与停滞

`JointPhaseSpaceIterationRecord` 增加：

```cpp
double poisson_current_relative;
double poisson_current_residual;
double weighted_continuity_defect;
int phase_converged;
int poisson_converged;
int pairing_converged;
```

20 次后仍失败则不接受状态、不提高上限、不放宽阈值，返回：

```text
failure_stage=joint_midpoint_poisson_current_not_converged
```

并输出最近三次 pairing relative。

---

## 9. A5：正负测试

本阶段分两条互斥路径：

### 9.1 A5A：assembly-only 路径

仅当 A3R 为 `A-ASSEMBLY_ONLY_RESOLVED` 时执行。覆盖：

```text
smooth-background
smooth-perturbed-background
a2-zero-residual
a2-injected-residual
全部J0和B4
```

不得新增 `pairing-gate-negative`，因为 A4 未执行，也没有新的 pairing solver gate。

### 9.2 A5N：A4 solver-gate 路径

仅当 A4 已执行并通过构建时，覆盖 `smooth-background`、`smooth-perturbed-background` 和测试专用 `pairing-gate-negative`。测试 hook 不得进入生产 CLI。

正测试必须：

```text
accepted=1
finite=1
gauss_ok=1
converged=1
failure_code=0
phase_converged=1
poisson_converged=1
pairing_converged=1
```

$$
\frac{|R_{PJ}|}{\max(1,|\Delta U_E|,|W_J|)}\le10^{-9}
$$

$$
\frac{|R_E|}{\max(1,|W_u|,|\Delta U_E|)}\le10^{-8}
$$

仅 A5N 的负测试必须：

```text
accepted=0
failure_stage=joint_midpoint_poisson_current_not_converged
```

---

## 10. A6：步长缩放

J1 test 增加 `--dt-scale <positive finite>`，默认 1.0，只控制测试 dt。运行 1.0 和 0.5，两档均须独立通过 A5；仅 dt/2 通过则 A6 FAIL。

### 10.1 10 步累计残差门

J1 test 同时增加：

```text
--steps <positive integer>
```

默认值为 1，只用于测试，不进入生产 CLI。完成单步 `dt`/`dt/2` 后，在小网格 `smooth-perturbed-background` 上运行 10 个连续接受步，并输出：

```text
accepted_step_count
cumulative_signed_energy_residual
cumulative_absolute_energy_residual
cumulative_exchange_scale
cumulative_relative_energy_residual
max_step_relative_energy_residual
initial_domain_energy
cumulative_drift_budget
cumulative_drift_pass
```

定义：

$$
S_{\mathrm{exchange}}
=\sum_{n=0}^{N-1}
\max\left(
1,
|W_u^n|,
|\Delta U_E^n|,
|W_{\mathrm{electrode}}^n|
\right)
$$

$$
R_{\mathrm{cum}}=\sum_{n=0}^{N-1}R_E^n
$$

累计预算为：

$$
B_{\mathrm{cum}}
=\max\left(
10^{-8}S_{\mathrm{exchange}},
10^{-12}\max(1,|U_{\mathrm{domain}}^0|)
\right)
$$

验收必须同时满足：

```text
accepted_step_count=10
每一步finite=1、gauss_ok=1、failure_code=0
max_step_relative_energy_residual<=1e-8
abs(cumulative_signed_energy_residual)<=cumulative_drift_budget
cumulative_relative_energy_residual<=1e-8
```

该累计门用于保证适度放宽的舍入/solver 阈值不会产生宏观同号漂移。不得通过每步清零累计量或只输出最后一步残差通过。

100 步累计测试不再是 A6/A7/A8 的阻断门。只有 A8 全部通过且用户需要更长数值稳定性证据时才可选运行；未运行不得据此判定当前阶段 FAIL。

### 10.2 多步接受态的 signed residual 契约

A6 首次10步运行在第二步暴露：第一步最终 code-76 门允许舍入级 signed mass，
但下一步初始 residual 仍以 `allow_negative_probe=false` 拒绝任何负值，导致
`failure_code=71/joint_midpoint_initial_residual`。这是状态延续契约矛盾，
不是累计能量失败。

修复要求：`advance_joint_midpoint()` 对初始 candidate 调用 `evaluate()` 时必须
使用 signed residual domain，与 Jacobian probe 和 line-search trial 一致：

```cpp
evaluate(...,
         finite,
         true,
         accepted_e_local)
```

该修改不等于关闭正性门。最终候选仍必须经过原有 code-76：

```text
candidate_min_mass >= -negative_tolerance
```

禁止修改 `negative_tolerance`、禁止裁剪 accepted state、禁止将 code-76 改为软接受。

10步回归必须证明：上一接受态可以成为下一步初值；若最终 code-76 真实失败，
仍按原逻辑停止。

---

## 11. A7：MPI ownership

相同全局网格、初态、dt、case 运行 1/2/5 rank，比较：

```text
R_P
W_C
R_PJ
prediction_error
continuity_charge_linf
charge_projection_mismatch_linf
iterations
accepted
failure_code
```

三档均满足 A5 且：

$$
|R_A^{\mathrm{pred}}|\le\tau_A
$$

### 11.1 不得混用的两类 MPI 比较尺度

结构恒等式必须使用舍入尺度比较：

```text
prediction_error <= tau_A
charge_projection_mismatch_linf <= continuity_roundoff_bound
weighted_defect_reconstruction_error <= tau_A
Poisson scalar identity满足Gate F
```

非线性终止量必须使用 solver tolerance 比较，不能要求跨 rank 达到
`tau_A`。不同 MPI 分区改变归约顺序和有限差分 Jacobian 的末位，允许 Newton
在不同迭代次数、但同一个收敛球内停止。

定义每档：

$$
r_p^{(m)}
=\frac{|R_{PJ}^{(m)}|}{S_{PJ}^{(m)}}
$$

其中 `m` 为 rank 数，`S_PJ` 必须使用测试输出的
`candidate_poisson_current_scale`。A7 使用当前二进制实际输出的
`pairing_tolerance`；若旧结果来自更严格门，也允许按该旧门验收，但必须记录
二进制门值。

### 11.2 修订后的 A7 门

每个 rank 档必须：

```text
accepted=1
failure_code=0
phase_converged=1
poisson_converged=1
pairing_converged=1
candidate_poisson_current_relative<=pairing_tolerance
```

跨 rank 要求：

$$
\max_m r_p^{(m)}-\min_m r_p^{(m)}
\le\mathrm{pairing\_tolerance}
$$

$$
\max_m R_{PJ}^{(m)}-\min_m R_{PJ}^{(m)}
\le
\mathrm{pairing\_tolerance}
\max_m S_{PJ}^{(m)}+\tau_A
$$

迭代次数可以不同，不作为 FAIL 条件；必须全部在最大迭代数内正常收敛。

### 11.3 ownership 真失败的定义

只有出现以下任一情况才归因于 MPI ownership：

- 某 rank 档的 prediction identity 超过 `tau_A`；
- continuity/charge projection 超过各自舍入界；
- shared face 重复计数导致残差近似按 rank 数系统增长；
- 某档在相同 solver tolerance 下不收敛或产生不同 failure code；
- normalized residual spread 超过本节 solver-aware 门。

禁止仅因 `R_PJ` 绝对值不逐位一致或迭代次数不同判定 ownership FAIL。

---

## 12. A8：最终门

重跑全部 J0、B4 seam adjoint、F10 两个 case、A2 两个 case、A-S0/A-S1 分解、A-FS 全部制造解、A6 两档和 A7 三档。只有 A4 实际执行时才要求 A5N 负测试。全部通过后才可写：

```text
F10 Scenario A resolved:
Poisson scalar identity, charge continuity,
Poisson-current pairing and total J1 energy identity
close on the same converged midpoint candidate.
```

A8 PASS 后才允许更新根因总结并进入 J1 MPI 最终门。J2 必须为生产 open/reservoir 边界重新推导伴随映射和边界功，不得直接复用 J1 periodic 规则。

---

## 13. 文件矩阵

| 阶段 | 文件 | 内容 |
|---|---|---|
| A0 | `docs/情形A执行情况.md` | 基线 |
| A1 | `vpfp_integrator.h/.cpp` | 只读分解诊断 |
| A1 | J1 energy test | 输出字段 |
| A2 | 现有 J0/J1 test | 两个生产直调 case |
| A3 | 执行报告 | 唯一分类 |
| A-S0 | `vpfp_integrator.h/.cpp`、J1 test | 装配/输运残差拆分 |
| A-S1 | `vpfp_integrator.cpp`、J0/J1 test | J1稳定增量电荷装配 |
| A-FS | `open_electrostatic_solver.h/.cpp`、Poisson identity test | 稳定直接差分与长双精度求和 |
| A3R | 执行报告 | 修复后唯一分类 |
| A4 | `vpfp_integrator.h/.cpp` | candidate metric、收敛门、受限 line search |
| A5A | J0/J1 test | assembly-only 回归，无负hook |
| A5N | J1 energy test | A4正负门 |
| A6 | J1 energy test | `--dt-scale` |
| A7 | 无生产修改 | MPI 回归 |
| A8 | 无生产修改 | 全回归 |

需要修改表外文件时必须停止。

---

## 14. 集群命令

全部由用户在集群根目录执行：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
```

### 14.1 A0

```bash
git rev-parse HEAD
git status --short
md5sum src/vpfp_integrator.cpp src/vpfp_integrator.h \
  src/open_electrostatic_solver.cpp src/open_electrostatic_solver.h \
  src/joint_phase_space_midpoint.cpp \
  tests/joint_phase_space_midpoint_energy_test.cpp
```

### 14.2 A1

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a1_${RUN_ID}"
mkdir -p "$OUT"
yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/a1_decomposition.result" \
  > "$OUT/a1_decomposition.out" \
  2> "$OUT/a1_decomposition.err"
sed -n '1,520p' "$OUT/a1_decomposition.result"
sed -n '1,240p' "$OUT/a1_decomposition.err"
```

A1 可能因旧 code 75 非零退出；按第 5 节诊断门验收。

### 14.3 A2

```bash
cmake --build build --target joint_phase_space_midpoint_unit_test -j4
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a2_${RUN_ID}"
mkdir -p "$OUT"
yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all --result "$OUT/a2_identity.result" \
  > "$OUT/a2_identity.out" 2> "$OUT/a2_identity.err" || exit 1
sed -n '1,620p' "$OUT/a2_identity.result"
```

### 14.4 A3

A3 无运行命令。读取 A1/A2 结果并机械分类。当前已命中中间分类 `A-S_SOURCE_OR_U_BOUNDARY_OR_ASSEMBLY`，下一步执行 A-S0。

### 14.4A A-S0

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_as0_${RUN_ID}"
mkdir -p "$OUT"
yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/as0_charge_assembly.result" \
  > "$OUT/as0_charge_assembly.out" \
  2> "$OUT/as0_charge_assembly.err"
sed -n '1,680p' "$OUT/as0_charge_assembly.result"
sed -n '1,260p' "$OUT/as0_charge_assembly.err"
```

A-S0 仍可能因旧 code 75 非零退出。按第 7A.6 节分类，不得只看进程退出码。

### 14.4B A-S1

只有 A-S0 为 `A-S-CHARGE_ASSEMBLY_ROUNDOFF` 才运行：

```bash
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_as1_${RUN_ID}"
mkdir -p "$OUT"
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all --result "$OUT/as1_j0_all.result" \
  > "$OUT/as1_j0_all.out" 2> "$OUT/as1_j0_all.err" || exit 1
for CASE_NAME in smooth-background smooth-perturbed-background; do
  yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE_NAME" --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err"
done
```

### 14.4C A-FS

```bash
cmake --build build --target \
  vpfp_poisson_work_identity_test \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_afs_${RUN_ID}"
mkdir -p "$OUT"

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/vpfp_poisson_work_identity_test \
  --case all --result "$OUT/afs_poisson_identity.result" \
  > "$OUT/afs_poisson_identity.out" \
  2> "$OUT/afs_poisson_identity.err" || exit 1

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all --result "$OUT/afs_j0_all.result" \
  > "$OUT/afs_j0_all.out" 2> "$OUT/afs_j0_all.err" || exit 1

for CASE_NAME in smooth-background smooth-perturbed-background; do
  yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE_NAME" --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err"
done

sed -n '1,620p' "$OUT/afs_poisson_identity.result"
sed -n '1,620p' "$OUT/afs_j0_all.result"
grep -H -E '^(status=|failure_code=|joint_midpoint_poisson_scalar_identity_residual=|poisson_identity_scale=|poisson_identity_roundoff_bound_8192=|poisson_identity_roundoff_bound_16384=|poisson_identity_residual_to_bound_ratio_8192=|poisson_identity_residual_to_bound_ratio_16384=|poisson_scalar_identity_pass=|joint_midpoint_current_pair_residual=|joint_midpoint_weighted_defect_reconstruction_error=|endpoint_nonzero_)' \
  "$OUT"/smooth-*.result
```

需要回传该目录中的全部 `.result/.out/.err`。任何制造解或旧回归失败都禁止 A3R。

### 14.4D A3R

A3R 无命令。读取 A-S0/A-S1/A-FS `.result`，按第 7D 节输出唯一分类。

### 14.5 A4

```bash
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
```

### 14.6 A5

若 A3R 为 assembly-only，执行 A5A：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a5a_${RUN_ID}"
mkdir -p "$OUT"
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all --result "$OUT/j0_a2_all.result" \
  > "$OUT/j0_a2_all.out" 2> "$OUT/j0_a2_all.err" || exit 1
for CASE_NAME in smooth-background smooth-perturbed-background; do
  yhrun -N 1 -n 1 --cpu-bind=cores \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE_NAME" --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err" || exit 1
done
```

若 A4 已执行，执行 A5N：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a5n_${RUN_ID}"
mkdir -p "$OUT"
for CASE_NAME in smooth-background smooth-perturbed-background pairing-gate-negative; do
  yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE_NAME" --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err"
done
```

负测试预期非零退出；PASS 表示正确拒绝。

### 14.7 A6

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a6_${RUN_ID}"
mkdir -p "$OUT"
for DT_SCALE in 1.0 0.5; do
  TAG=$(printf '%s' "$DT_SCALE" | tr '.' 'p')
  yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case smooth-perturbed-background --dt-scale "$DT_SCALE" \
    --result "$OUT/dt_${TAG}.result" \
    > "$OUT/dt_${TAG}.out" 2> "$OUT/dt_${TAG}.err" || exit 1
done

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --dt-scale 1.0 --steps 10 \
  --result "$OUT/cumulative_10steps.result" \
  > "$OUT/cumulative_10steps.out" \
  2> "$OUT/cumulative_10steps.err" || exit 1

sed -n '1,620p' "$OUT/cumulative_10steps.result"
```

### 14.8 A7

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a7_${RUN_ID}"
mkdir -p "$OUT"
for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case smooth-perturbed-background --dt-scale 1.0 \
    --result "$OUT/n${NP}.result" \
    > "$OUT/n${NP}.out" 2> "$OUT/n${NP}.err" || exit 1
done
```

### 14.9 A8

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_a_a8_${RUN_ID}"
mkdir -p "$OUT"
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all --result "$OUT/j0_all.result" \
  > "$OUT/j0_all.out" 2> "$OUT/j0_all.err" || exit 1
for CASE_NAME in smooth-background smooth-perturbed-background; do
  yhrun -N 1 -n 1 --cpu-bind=cores \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE_NAME" --dt-scale 1.0 \
    --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err" || exit 1
done
```

---

## 15. 报告格式

在 `docs/情形A执行情况.md` 追加：

````markdown
# 阶段：A? 执行报告

## 当前阶段
`A?`

## 前置阶段
- `<PASS证据>`

## 实际修改文件
- `<绝对路径>`

## 修改内容
1. `<函数、位置、行为>`

## 明确未修改
- `<物理与算法对象>`

## 编译命令
```bash
...
```

## 测试命令
```bash
...
```

## 数值结果
```text
...
```

## 阶段门
PASS、FAIL 或待集群执行

- PASS：说明允许进入哪个阶段，但不要自行进入。
- FAIL：列首个失败字段、数值、源码位置和允许检查范围；停止。
- 待集群执行：列用户需运行的命令及需回传文件；不得推测。
````

---

## 16. 最终复核

提交前确认未：修改生产物理或边界、修改 Poisson stencil、修改情形 B `G*`、从密度差覆盖生产电流、用残差修正状态、放宽任何门、提高迭代上限隐藏停滞、A-S0/A-S1/A-FS/A3R 未完成就进入 A4、assembly-only 已解决仍强行执行 A4、A8 未通过就进入 J2/J3 或长跑。

违反任一项，当前阶段直接 FAIL。
