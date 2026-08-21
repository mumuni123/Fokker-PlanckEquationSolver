# VPFP 联合相空间时间中心能量闭合重构
# GitHub 源码实证强制修复实施方案 v3

> **适用仓库**：`mumuni123/GPT-fokker-planck-equation-solver`  
> **源码基线**：`main` @ `ad09bf2cb52f92b21a3a3cb71c75b6d798b3accd`  
> **基线提交说明**：`重构后未修复能量问题`  
> **文档性质**：执行规约，不是讨论稿。  
> **目标执行模型**：GPT-5.5 / GPT-5.6 Luna 或同等级代码执行模型。  
> **执行原则**：执行模型不得重新设计方案；只按本文顺序修改、编译、运行和验收。  
> **默认生产模式**：始终保持 `STRANG_PPM`。  
> **本文件取代**：此前的 `强制修复实施方案.md` 与 `源码适配强制修复实施方案_v2.md` 中与 J0/J1 核心算子有关的内容。

---

# 0. 最重要的结论变化

重新直接阅读当前 GitHub 源码后，必须修正此前对 J1 code-75 的解释。

此前根据运行结果认为：

```text
当前两个 center trace 虽共享 M_mid，
但没有共同 bracket，因此 center-trace 路线已经被 code-75 证伪。
```

这个结论现在**不能直接成立**。

原因是当前 production joint 路径实际并没有实现文档所描述的那个“正确的联合时间中心 center-trace 算子”。

当前源码中至少存在四个先于“中心 trace 是否具有结构”的确定性实现错误：

```text
A1. evaluate_local_residual() 使用 M_candidate，而不是 0.5*(M_old+M_candidate)
A2. u_flux_rate 对 cell-integrated M 重复乘 dx * uperp_ring_areas
A3. u force 使用 candidate 最终 Ex cell-center，而不是 production exact-dual E_pair
A4. J1 mass flux 硬编码周期回卷，而 Poisson/field work 使用 open/Dirichlet 拓扑
```

另外：

```text
A5. 旧 J0 PASS 没有验证 x-charge-current 与 u-kinetic-work 的交叉恒等式
```

所以当前正确的根因层级必须改成：

```text
第一层：已确认的源码实现契约错误 A1--A5
第二层：修复 A1--A5 后仍需验证的非均匀-u 离散伴随问题
第三层：若第二层仍失败，再讨论更一般的 phase-space bracket
```

本 v3 不允许直接跳到 v2 的 arbitrary common-corner formula。

---

# 1. 当前源码的权威数据契约

以下内容已经由当前 GitHub 源码确认，执行模型不得重新猜测。

## 1.1 Background `Species::f` 的真实含义

对背景电子：

```text
Species::cylindrical_mass_representation = true
```

`Species::f` 存的是 cylindrical phase-space **cell-integrated mass / number**：

$$
M_{i,j,k}
=
f_{3,i,j,k}
\Delta x
\Delta u_{\parallel,j}
A_{\perp,k}.
$$

其中：

$$
A_{\perp,k}
=
\pi
\left(
u_{\perp,k+1/2}^2
-
u_{\perp,k-1/2}^2
\right).
$$

源码对应：

```cpp
f[idx3(...)] =
    f3 * sgrid->dx * cgrid.cell_phase_volume(iv, imu);
```

以及：

```cpp
cell_phase_volume(j,k)
=
upar_widths[j] * uperp_ring_areas[k];
```

所以：

```text
M 的单位 = m^-2
M 已经包含 dx
M 已经包含 upar_width
M 已经包含 uperp_ring_area
```

这条是整个修复的最高优先级量纲契约。

---

## 1.2 Background number 与 kinetic energy

production：

```cpp
total_particle_number = sum(M)
```

production：

```cpp
total_kinetic_energy = sum(M * kinetic_energy[j,k])
```

因此不得在 joint energy ledger 中再次乘：

```text
dx
du
ring
```

到 `M*K` 上。

---

## 1.3 velocity-only table 布局

production：

```text
kinetic_energy[j,k]
vx[j,k]
```

索引：

```cpp
qv = j * Nmu + k;
```

distribution local slab：

```cpp
qd = (ix_local * Nv + j) * Nmu + k;
```

永久禁止：

```cpp
kinetic_energy[cell_index(ix,j,k)]
```

---

## 1.4 production u remap 的 mass-flux 语义

`ConservativePpmRemap::advect_u_parallel()` 已经给出权威范式：

```cpp
m = input.f[upar_index(ix,j,k)];
upar_average_[j] = m / upar_widths[j];
```

然后 finite-volume update：

```cpp
m_new =
    input_mass
    - upar_swept_[j+1]
    + upar_swept_[j];
```

这里没有：

```text
* dx
* uperp_ring_areas[k]
```

因为这两个量已经包含在 `M` 中。

因此 joint u-flux rate 的基本单位必须与：

```text
M / s
```

一致。

---

# 2. 当前 joint 源码中必须先修复的确定性错误

---

## 2.1 A1：local residual 实际不是 midpoint

当前：

```text
build_periodic_center_flux()
```

的 global/J0 helper 会先构造：

$$
M^H
=
\frac12(M^n+M^{n+1}).
$$

但 production joint 路径实际调用：

```text
evaluate_local_residual()
```

而该函数直接从：

```cpp
m_candidate_local
```

构造 x flux 与 u flux。

也就是说当前实际 J1 是：

$$
F(M^{n+1}),
$$

不是：

$$
F\!\left(
\frac{M^n+M^{n+1}}2
\right).
$$

### 强制修复

`evaluate_local_residual()` 内第一步必须构造：

```cpp
m_mid_local[q] =
    0.5 * (m_old_local[q] + m_candidate_local[q]);
```

所有：

```text
x trace
u trace
charge current
force-current diagnostic
```

只能读取 `m_mid_local`。

禁止任何 flux 继续读取 `m_candidate_local`。

---

## 2.2 A2：当前 u_flux_rate 重复乘 phase-space measure

当前 local joint 代码：

```cpp
ftrace =
    0.5 *
    (m_candidate[left]  / upar_widths[j-1]
   + m_candidate[right] / upar_widths[j]);

u_flux_rate =
    a * ftrace * sg.dx * uperp_ring_areas[k];
```

这是错误的。

因为：

```text
M 已含 dx * du * ring
M/du 已经是对 x 与 u_perp 积分后的 u_parallel line mass
```

正确的 center u-face rate 是：

$$
Q^u_{i,j+1/2,k}
=
a_i
\frac12
\left(
\frac{M^H_{i,j,k}}{\Delta u_j}
+
\frac{M^H_{i,j+1,k}}{\Delta u_{j+1}}
\right).
$$

即代码必须改为：

```cpp
u_flux_rate =
    a * 0.5 *
    (m_mid[left]  / upar_widths[j]
   + m_mid[right] / upar_widths[j+1]);
```

**不得再乘：**

```cpp
sg.dx
uperp_ring_areas[k]
```

---

## 2.3 A3：当前 force field 不是 production exact dual field

当前 `advance_joint_midpoint()` 的 candidate evaluation：

```text
candidate M
 -> candidate rho
 -> OpenElectrostaticSolver::solve(candidate_fields)
 -> e_local[ix] = candidate_fields.Ex[ng+ix]
 -> evaluate_local_residual(... e_local ...)
```

也就是说当前 u force 使用的是：

```text
candidate final Ex cell-center
```

但 production 已经有一个完全不同、专门为离散 Poisson 功恒等式构造的 exact dual：

```cpp
OpenElectrostaticSolver::build_potential_pairing_field(
    before,
    after,
    pairing_face,
    rank,
    size);
```

这个 pairing field 不是简单：

```text
0.5*(E^n + E^{n+1})
```

也不是简单：

```text
-(phi[i+1]-phi[i])/dx
```

production 实际使用：

$$
\bar\phi_i
=
\phi_i
+
\frac{\Delta x}{12}
(E_{i+1/2}-E_{i-1/2}),
$$

然后时间中心：

$$
\bar\phi_i^H
=
\frac12
(\bar\phi_i^n+\bar\phi_i^{n+1}).
$$

再构造 exact face dual。

内部 face：

$$
E^{pair}_{i+1/2}
=
\frac{
\bar\phi_i^H-\bar\phi_{i+1}^H
}{\Delta x}.
$$

物理端点不是这个公式，而是：

$$
E^{pair}_{1/2}
=
-\frac{2\bar\phi_0^H}{\Delta x},
$$

$$
E^{pair}_{N+1/2}
=
+\frac{2\bar\phi_{N-1}^H}{\Delta x}.
$$

MPI shared face 使用相邻 rank 的 paired cell-average potential。

### 强制修复

candidate residual evaluation 必须在 candidate Poisson solve 后调用：

```cpp
field_solver_.build_potential_pairing_field(
    fields_n,
    candidate_fields,
    pairing_face,
    mpi_rank,
    mpi_size);
```

然后 force 用 production gather：

$$
E^{pair}_{i}
=
\frac12
\left(
E^{pair}_{i-1/2}
+
E^{pair}_{i+1/2}
\right).
$$

禁止再把：

```cpp
candidate_fields.Ex[ng+ix]
```

作为 joint u-force 的权威场。

---

## 2.4 A4：当前 J1 混用 periodic mass topology 与 open field topology

当前 `joint_phase_space_midpoint_energy_test`：

```cpp
background boundary = PERIODIC
```

同时：

```cpp
field boundary = DIRICHLET_PHI
phi_left = 0
phi_right = 0
```

而 `OpenElectrostaticSolver` 是 open-domain solver。

更严重的是：

```cpp
evaluate_local_residual()
```

自己硬编码 MPI periodic wrap：

```text
rank 0 的 left neighbour = last rank
last rank 的 right neighbour = rank 0
```

所以即使上层 boundary config 改了，joint residual 仍然周期回卷。

### 结论

当前 code-75 的 J1 根本不是一个拓扑一致的：

```text
periodic Vlasov + periodic field
```

系统，也不是：

```text
open Vlasov + open field
```

系统。

它是：

```text
periodic mass flux
+
open/Dirichlet Poisson
```

的混合。

### 强制修复

J1 core test 从现在开始**禁止 periodic x**。

在 production open boundary 的完整能量推导完成前，J1 core 使用一个明确的：

```text
ADJOINT_ONE_SIDED_TEST
```

x-boundary closure。

该 closure 只用于 J0/J1 数学验证，不允许 production CLI 使用。

MPI 内部 shared face：

```text
正常非周期 neighbour exchange
```

global physical endpoints：

```text
不 wrap
```

具体 flux 见第 8 节。

---

## 2.5 A5：旧 J0 PASS 不是 joint-work PASS

当前 `joint_phase_space_midpoint_unit_test.cpp` 做了：

```text
1. kinetic state change vs u-flux kinetic work
2. Poisson field-energy identity
3. G/G* identity using x charge current
```

但 Poisson 的 `rho_delta` 在测试中是**直接从 x face current 重新构造的**：

```cpp
charge_delta =
    -dt * (J_right-J_left)/dx;
```

因此旧 J0 实际证明的是：

```text
u flux 自己的 kinetic telescoping 正确
Poisson 自己的 work identity 正确
x current 自己的 continuity/G-G* 正确
```

它没有证明：

$$
\boxed{
\Delta K_u
=
\Delta t
\langle E_{pair},J_x\rangle
}
$$

所以：

```text
J0=PASS
```

必须重新命名为：

```text
J0_component_self_consistency=PASS
```

而不是：

```text
joint Hamiltonian operator PASS
```

---

# 3. 修订后的根因判定

当前允许写进项目结论的内容只有：

## 3.1 已确定

```text
A1 midpoint implementation mismatch
A2 u-flux phase-space measure duplication
A3 wrong force-field time/dual object
A4 particle/field topology mismatch
A5 J0 missing cross-work gate
```

这五项由源码直接确定。

## 3.2 尚需重新验证

修复 A1--A5 后，当前非均匀 u 网格上的：

```text
build_hamiltonian_velocity()
```

是否与 center u-face mass trace互为离散伴随。

源码表明它当前用：

$$
v_j^{old}
=
\frac{
K_{j+1}-K_{j-1}
}{
m_ec(u_{j+1}-u_{j-1})
}.
$$

而 center u-flux 的 kinetic work 自然诱导的伴随速度是另一种公式。

所以结构性问题仍然高度可疑，但必须用新的 cross-work test 正式判定。

---

# 4. 本 v3 的最终 J1 core scheme

新 scheme 固定命名：

```text
ADJOINT_CENTER_MIDPOINT_V1
```

内部 manifest/schema 名称：

```text
joint_phase_space_scheme=adjoint-center-midpoint-v1
joint_phase_space_scheme_schema=3
```

在 J2 open-boundary 完成前不写 production checkpoint；这里只用于测试结果和启动摘要。

---

# 5. 新 scheme 的离散状态

## 5.1 时间中心 mass

唯一允许：

$$
M^H_{i,j,k}
=
\frac12
\left(
M^n_{i,j,k}
+
M^{n+1}_{i,j,k}
\right).
$$

---

## 5.2 production kinetic table

唯一允许：

$$
K_{j,k}
=
\texttt{cgrid.kinetic\_energy[j*Nmu+k]}.
$$

定义内部 u face：

$$
\Delta K_{j+1/2,k}
=
K_{j+1,k}-K_{j,k}.
$$

---

# 6. center u-trace 的严格离散伴随速度

这是 v3 相对 v2 最重要的数学更新。

不是继续使用当前：

```text
center-distance Hamiltonian velocity
```

也不是继续使用 P3-V.2 的 cell projection。

---

## 6.1 u-flux

内部 u face：

$$
Q^u_{i,j+1/2,k}
=
a_i^H
\frac12
\left(
\frac{M^H_{i,j,k}}{\Delta u_j}
+
\frac{M^H_{i,j+1,k}}{\Delta u_{j+1}}
\right),
$$

其中：

$$
a_i^H
=
\frac{q_e E^{pair}_i}{m_e c}.
$$

velocity endpoints：

$$
Q^u_{i,-1/2,k}=0,
\qquad
Q^u_{i,N_v+1/2,k}=0
$$

仅用于 J1 core。

---

## 6.2 从该 u-trace 代数推导 velocity

将 kinetic work 按 cell mass 重排。

对每个 velocity cell 定义：

$$
\Delta K^-_{j,k}
=
\begin{cases}
K_{j,k}-K_{j-1,k}, & j>0,\\
0,&j=0,
\end{cases}
$$

$$
\Delta K^+_{j,k}
=
\begin{cases}
K_{j+1,k}-K_{j,k},&j<N_v-1,\\
0,&j=N_v-1.
\end{cases}
$$

定义：

$$
\boxed{
v^{\dagger}_{j,k}
=
\frac{
\Delta K^-_{j,k}
+
\Delta K^+_{j,k}
}{
2m_ec\,\Delta u_j
}
}
$$

其中：

```text
Delta u_j = cgrid.upar_widths[j]
```

不是：

```text
upar_cells[j+1]-upar_cells[j-1]
```

---

## 6.3 为什么这个速度是唯一允许的 J1 core 速度

由 u flux：

$$
\dot K_{u,i}
=
\sum_{j+1/2,k}
\Delta K_{j+1/2,k}
Q^u_{i,j+1/2,k}.
$$

代入 center trace 并按 cell `j` 重排：

$$
\dot K_{u,i}
=
\frac{q_eE_i^{pair}}{m_ec}
\sum_{j,k}
\frac{M^H_{i,j,k}}{2\Delta u_j}
\left(
\Delta K^-_{j,k}
+
\Delta K^+_{j,k}
\right).
$$

所以：

$$
\boxed{
\dot K_{u,i}
=
E_i^{pair}\Delta x
J^{force}_i
}
$$

其中：

$$
\boxed{
J^{force}_i
=
\frac{q_e}{\Delta x}
\sum_{j,k}
v^\dagger_{j,k}M^H_{i,j,k}
}
$$

这是代数恒等式。

不是拟合，不是 correction，不是 energy scaling。

---

# 7. x flux 必须使用同一个 v-dagger

内部 x face：

$$
Q^x_{i+1/2,j,k}
=
v^\dagger_{j,k}
\frac{
M^H_{i,j,k}
+
M^H_{i+1,j,k}
}{
2\Delta x}.
$$

charge current：

$$
J^x_{i+1/2}
=
q_e
\sum_{j,k}
Q^x_{i+1/2,j,k}.
$$

因此内部 face 自动满足：

$$
\boxed{
J^x_{i+1/2}
=
\frac12
\left(
J^{force}_i
+
J^{force}_{i+1}
\right)
}
$$

这正是 production `G*` 在内部 face 的定义。

这就是本方案要修复的“共同离散乘积法则”。

---

# 8. J1 core 的 x physical endpoint closure

production `FieldParticlePowerAudit` 的 `G*` 在 global physical endpoint 定义为：

$$
(G^*J)_{1/2}=J_0,
$$

$$
(G^*J)_{N+1/2}=J_{N-1}.
$$

所以 J1 core test 必须使用一侧 trace：

左端：

$$
Q^x_{1/2,j,k}
=
v^\dagger_{j,k}
\frac{M^H_{0,j,k}}{\Delta x}.
$$

右端：

$$
Q^x_{N+1/2,j,k}
=
v^\dagger_{j,k}
\frac{M^H_{N-1,j,k}}{\Delta x}.
$$

于是：

$$
J^x_{1/2}=J^{force}_0,
$$

$$
J^x_{N+1/2}=J^{force}_{N-1}.
$$

从而所有 face：

$$
\boxed{
J_x=G^*J_{force}
}
$$

到舍入误差。

### 注意

这个 endpoint closure：

```text
ADJOINT_ONE_SIDED_TEST
```

只用于 J1 core 数学验证。

它不是最终 reservoir / absorbing production boundary。

production open boundary 推迟到 J2，并且 J2 在本文件末尾被明确 BLOCK。

---

# 9. production pairing field 与 G/G* 链

每个 candidate residual：

```text
fields_n = accepted old Poisson state
fields_np1 = candidate Poisson state
```

调用：

```cpp
build_potential_pairing_field(
    fields_n,
    fields_np1,
    e_pair_face,
    rank,
    size);
```

然后：

$$
E^{pair}_i
=
\frac12
(E^{pair}_{i-1/2}+E^{pair}_{i+1/2}).
$$

production 已有：

$$
\Delta U_E-W_{electrode}
=
W_{\phi\Delta\rho},
$$

以及：

$$
W_{\phi\Delta\rho}
=
-\Delta t
\langle E_{pair},J_{charge}\rangle_f.
$$

新 joint scheme 又强制：

$$
J_{charge}=G^*J_{force},
$$

所以：

$$
\Delta t
\langle E_{pair},J_{charge}\rangle_f
=
\Delta t
\sum_i
\Delta x E_i^{pair}J_i^{force}.
$$

由第 6 节：

$$
\Delta K_u
=
\Delta t
\sum_i
\Delta x E_i^{pair}J_i^{force}.
$$

最终：

$$
\boxed{
\Delta K_u
+
\Delta U_E
-
W_{electrode}
=
O(\epsilon_{sum})
}
$$

必须在 operator level 成立。

---

# 10. 新旧 velocity 的 A/B 测试

不得直接删除现有：

```cpp
build_hamiltonian_velocity()
```

先保留为：

```text
LEGACY_CENTER_DISTANCE_DIAGNOSTIC
```

新增：

```cpp
build_u_trace_adjoint_velocity(...)
```

对应：

```text
ADJOINT_U_TRACE
```

J0 必须同时跑：

```text
legacy velocity
adjoint velocity
```

预期：

```text
legacy:
  非均匀 u grid cross-work 可 FAIL

adjoint:
  cross-work 必须 PASS 到 stable-sum roundoff
```

只有 A/B 实际证明后，才允许把 J1 core 固定到：

```text
ADJOINT_U_TRACE
```

---

# 11. 需要修改的真实源码文件

第一轮只允许修改：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp

src/vpfp_integrator.h
src/vpfp_integrator.cpp

tests/joint_phase_space_midpoint_unit_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp

CMakeLists.txt
```

第一轮禁止修改：

```text
src/grid.h
src/species.*
src/conservative_ppm_remap.*
src/open_electrostatic_solver.*
src/field_particle_power_audit.*
src/open_boundary.*
src/beam_pic.*
src/background_tail_pic.*
src/cylindrical_fp_collision.*
src/vpfp_checkpoint.*
```

原因：

```text
这些 production 组件当前已经提供所需权威接口；
第一轮修复不能通过修改它们让 joint test 通过。
```

---

# 12. `joint_phase_space_midpoint.h` 的强制修改

---

## 12.1 新增 enum

```cpp
enum class JointVelocityMomentMode {
    LEGACY_CENTER_DISTANCE_DIAGNOSTIC = 0,
    ADJOINT_U_TRACE = 1
};

enum class JointXBoundaryClosure {
    ADJOINT_ONE_SIDED_TEST = 0
};
```

不新增 production reservoir mode。

---

## 12.2 扩展 `JointPhaseSpaceFluxBundle`

必须至少新增：

```cpp
std::vector<double> midpoint_mass;

std::vector<double> pairing_field_face;
std::vector<double> pairing_field_cell;

std::vector<double> force_current_cell;
std::vector<double> gstar_force_current_face;
std::vector<double> current_dual_face_residual;

double delta_ke_x_flux;
double delta_ke_u_flux;

double pairing_current_work;
double u_force_work;

double current_pair_residual;
double current_pair_abs_scale;
```

保留现有：

```text
x_flux_rate
u_flux_rate
charge_current_face
mass_delta_x
mass_delta_u
mass_delta_total
```

---

## 12.3 新 helper

```cpp
static std::vector<double> build_u_trace_adjoint_velocity(
    const CylindricalVelocityGrid& vg,
    double particle_mass);
```

---

## 12.4 新 G* helper

新增纯函数：

```cpp
static bool build_gstar_force_current(
    const SpatialGrid& sg,
    const std::vector<double>& force_current_cell,
    int mpi_rank,
    int mpi_size,
    std::vector<double>& gstar_face);
```

其实现必须逐字遵循 production `FieldParticlePowerAudit` 的当前定义：

```text
global left endpoint  -> first cell
global right endpoint -> last cell
local interior face   -> average adjacent local cells
MPI shared face       -> average neighbor boundary cells
```

禁止重新设计权重。

---

# 13. `evaluate_local_residual()` 新接口

旧接口：

```cpp
evaluate_local_residual(
    ...,
    e_cell_local,
    ...,
    allow_negative_probe);
```

改成：

```cpp
evaluate_local_residual(
    ...,
    pairing_face_local,
    particle_charge,
    particle_mass,
    JointVelocityMomentMode velocity_mode,
    JointXBoundaryClosure x_boundary_closure,
    ...,
    bool allow_negative_state);
```

注意：

```text
operator 只读 pairing_face
operator 不自己调用 Poisson
```

Poisson 仍由 integrator transaction 管理。

---

# 14. `evaluate_local_residual()` 固定执行顺序

严格按以下顺序：

```text
1. validate dimensions
2. validate all finite
3. build m_mid = 0.5*(m_old+m_candidate)
4. build velocity table
5. exchange m_mid internal MPI boundary slices with NONPERIODIC neighbors
6. gather pairing_face -> pairing_cell
7. build x_flux_rate from m_mid
8. build u_flux_rate from m_mid
9. build charge_current_face
10. build force_current_cell
11. build G*force_current_face
12. build current dual residual
13. build mass_delta_x
14. build mass_delta_u
15. build mass_delta_total
16. residual = candidate-old-mass_delta_total
17. compute residual norms
18. optional negative-state domain check
```

禁止：

```text
x flux 读 candidate
u flux 读 candidate
force current 用另一份 state
```

---

# 15. MPI mass halo：取消周期回卷

旧：

```text
rank 0 left <- last rank
last rank right <- rank 0
```

删除。

新：

```cpp
left =
    mpi_rank > 0
    ? mpi_rank - 1
    : MPI_PROC_NULL;

right =
    mpi_rank + 1 < mpi_size
    ? mpi_rank + 1
    : MPI_PROC_NULL;
```

只交换：

```text
internal MPI shared faces
```

global endpoint 使用：

```text
ADJOINT_ONE_SIDED_TEST
```

不读 ghost。

---

# 16. x flux 的实现

对 internal local face：

```cpp
trace =
    0.5 * (m_mid_left + m_mid_right);

x_flux_rate =
    v_adjoint * trace / dx;
```

MPI shared face：

```text
通过 exchange 得到远端 m_mid slice
两 rank 算出的 face flux 必须一致
```

global left endpoint：

```cpp
x_flux_rate =
    v_adjoint * m_mid_first_cell / dx;
```

global right endpoint：

```cpp
x_flux_rate =
    v_adjoint * m_mid_last_cell / dx;
```

---

# 17. u flux 的实现

```cpp
const double e_pair_cell =
    0.5 * (pairing_face[ix] + pairing_face[ix+1]);

const double a =
    particle_charge * e_pair_cell /
    (particle_mass * Const::c);

const double trace =
    0.5 *
    (m_mid[left]  / vg.upar_widths[jf-1]
   + m_mid[right] / vg.upar_widths[jf]);

u_flux_rate =
    a * trace;
```

明确禁止：

```cpp
* sg.dx
* vg.uperp_ring_areas[k]
```

---

# 18. force current 的实现

```cpp
J_force[ix] =
    particle_charge / sg.dx
    * sum_{j,k}(v_adjoint[j,k] * m_mid[ix,j,k]);
```

必须使用与 x flux **同一个**：

```text
v_adjoint
m_mid
```

---

# 19. charge current 的实现

```cpp
J_charge[face] =
    particle_charge
    * sum_{j,k}(x_flux_rate[face,j,k]);
```

禁止：

```text
从 rho_delta 反推 current
从 q∫vf 重新用 analytic vx 计算 current
从旧 PPM audit 复制 current
```

---

# 20. 当前对偶误差

逐 face：

```cpp
dual_face_residual[f] =
    J_charge[f] - Gstar_J_force[f];
```

全局：

```text
current_pair_linf
current_pair_l1 / abs scale
```

必须输出。

---

# 21. u kinetic work

直接从真正更新使用的 `u_flux_rate`：

$$
W_u
=
\Delta t
\sum_{i,j+1/2,k}
\Delta K_{j+1/2,k}
Q^u_{i,j+1/2,k}.
$$

使用 `long double` 累加。

---

# 22. pairing current work

必须使用 production face quadrature：

```text
global physical endpoint : dx/2
global interior face      : dx
MPI left-rank duplicate right face : 0
```

与 `FieldParticlePowerAudit` 完全一致。

$$
W_{EJ}
=
\Delta t
\sum_f
w_f
E^{pair}_f
J^{charge}_f.
$$

---

# 23. direct cross-work gate

定义：

$$
R_{cross}
=
W_u-W_{EJ}.
$$

这是新的核心 J0/J1 Gate。

不得再用：

```text
kinetic_work_residual + poisson_work_residual
```

替代。

---

# 24. stable-sum bound

定义：

$$
S_{cross}
=
\sum |u\text{-face kinetic work term}|
+
\sum |face E J term|.
$$

$$
T_{cross}
=
8192\epsilon_{double}S_{cross}.
$$

若：

```text
S_cross == 0
```

要求：

```text
R_cross == 0
```

不使用：

```text
max(1, S)
```

避免带单位残差与无量纲 `1` 混合。

---

# 25. 新 J0 测试结构

继续使用现有 target：

```text
joint_phase_space_midpoint_unit_test
```

不急着新增大量 target。

将其扩展为以下 case。

---

## J0-S1：source-unit-contract

验证：

```text
M [m^-2]
x_flux_rate [m^-2 s^-1]
u_flux_rate [m^-2 s^-1]
charge_current_face [A m^-2]
delta_ke [J]
kinetic work [J m^-2]
```

制造一个单一 `(ix,j,k)` mass，并显式检查：

```text
M / du
```

进入 u flux 后没有额外 `dx*ring`。

---

## J0-S2：local-midpoint-contract

单 rank。

构造：

```text
m_old != m_candidate
```

检查：

```text
local operator flux
```

必须等于显式传入：

```text
0.5*(old+candidate)
```

的 helper flux。

这个测试在旧源码上必须 FAIL。

修复后 PASS。

---

## J0-S3：pairing-field-contract

用：

```cpp
OpenElectrostaticSolver::build_potential_pairing_field()
```

得到 face field。

验证：

```text
operator 实际使用的 pairing_field_face
```

与 helper bitwise/roundoff一致。

禁止 operator 使用：

```text
candidate_fields.Ex
```

替代。

---

## J0-S4：legacy-cross-work

使用：

```text
LEGACY_CENTER_DISTANCE_DIAGNOSTIC
```

输出：

```text
legacy_current_pair_residual
legacy_current_pair_relative
```

不要求 PASS。

它只作为 A/B 基线。

---

## J0-S5：adjoint-cross-work

使用：

```text
ADJOINT_U_TRACE
```

要求：

```text
J_charge_face - G*J_force_face = stable-sum roundoff
W_u - W_EJ = stable-sum roundoff
```

必须 PASS。

---

## J0-S6：nonuniform-u

必须使用 production：

```text
CylindricalVelocityGrid::init_grid(...)
```

非均匀 `upar_widths`。

不能用手工 uniform grid 逃避问题。

---

## J0-S7：positive/negative field

分别测试：

```text
E_pair > 0
E_pair < 0
```

电子：

```text
charge = -Const::qe
```

---

## J0-S8：single/two/multivelocity

覆盖：

```text
single active velocity
two adjacent velocities
smooth Maxwellian-like M
```

---

# 26. J0 必须新增的输出字段

```text
source_mass_contract_pass
local_midpoint_contract_pass
u_flux_measure_contract_pass
pairing_field_contract_pass

legacy_current_pair_residual
legacy_current_pair_relative

adjoint_current_pair_residual
adjoint_current_pair_roundoff_bound
adjoint_current_pair_pass

adjoint_dual_face_linf
adjoint_dual_face_roundoff_bound
adjoint_dual_face_pass

u_work
pairing_current_work
cross_work_residual
cross_work_abs_scale
cross_work_roundoff_bound
cross_work_pass

status
```

---

# 27. J0 放行条件

只有：

```text
source_mass_contract_pass=1
local_midpoint_contract_pass=1
u_flux_measure_contract_pass=1
pairing_field_contract_pass=1
adjoint_current_pair_pass=1
adjoint_dual_face_pass=1
cross_work_pass=1
```

才能进入 J1。

`legacy` 可以 FAIL。

---

# 28. J1 integrator 的 field transaction

当前 `advance_joint_midpoint()` 已经正确采用：

```text
candidate_fields scratch
```

而不是直接改 accepted `fields`。

保留。

但 candidate evaluate 流程改成：

```text
candidate M
-> candidate rho
-> field_solver.solve(candidate_fields)
-> build_potential_pairing_field(fields_n,candidate_fields)
-> evaluate_local_residual(... pairing_face ...)
```

禁止：

```text
candidate Ex cell -> u force
```

---

# 29. J1 unknown

保持：

```text
M_candidate local physical slab
```

不新增 `phi` 为 Krylov unknown。

当前源码实际上已经是这种 Poisson-eliminated架构。

因此 v2 中“删除 phi unknown”的任务取消：

```text
当前源码无需执行该改动。
```

只需删除无用/误导的：

```text
phi_residual vector 作为 unknown 的语义
```

Gauss residual 保留为 diagnostic gate。

---

# 30. J1 Poisson/Gauss gate

保留当前相对定义：

$$
r_G
=
\frac{
\|D_xE-\rho/\epsilon_0\|_\infty
}{
\max(
\|D_xE\|_\infty,
\|\rho/\epsilon_0\|_\infty,
\text{safe floor}
)
}.
$$

不得恢复绝对：

```text
1e-8
```

比较。

---

# 31. J1 energy decomposition

在 `VpfpStepResult` 新增：

```cpp
double joint_midpoint_delta_ke_x_flux;
double joint_midpoint_delta_ke_u_flux;

double joint_midpoint_force_current_work;
double joint_midpoint_charge_current_work;

double joint_midpoint_current_pair_residual;

double joint_midpoint_poisson_potential_work;
double joint_midpoint_field_energy_change;
double joint_midpoint_electrode_work;
double joint_midpoint_poisson_identity_residual;

double joint_midpoint_state_flux_kinetic_residual;
double joint_midpoint_combined_energy_residual;
```

旧：

```cpp
joint_midpoint_energy_residual
```

保留兼容，但赋值为：

```text
combined_energy_residual
```

---

# 32. `delta_K_x_flux`

$$
\Delta K_x^{flux}
=
\sum_{i,j,k}
K_{j,k}
\Delta M^x_{i,j,k}.
$$

---

# 33. `delta_K_u_flux`

$$
\Delta K_u^{flux}
=
\sum_{i,j,k}
K_{j,k}
\Delta M^u_{i,j,k}.
$$

它应等于：

$$
W_u.
$$

---

# 34. state/flux kinetic gate

$$
\Delta K_{state}
=
K(M^{n+1})-K(M^n).
$$

定义：

$$
R_{state-flux}
=
\Delta K_{state}
-
\Delta K_x^{flux}
-
\Delta K_u^{flux}.
$$

nonlinear residual 收敛后必须到 stable-sum scale。

---

# 35. current-pair gate

$$
R_{current}
=
\Delta K_u^{flux}
-
W_{EJ}.
$$

必须 PASS。

这是当前 code-75 应该首先细分出的 residual。

---

# 36. Poisson gate

production：

```cpp
OpenPoissonWorkIdentity poisson_work =
    evaluate_work_identity(...);
```

直接保存：

```text
field_energy_change
electrode_work
potential_charge_work
residual
scale
```

要求：

```text
poisson_work.residual <= production Gate-F tolerance
```

禁止重新实现 Poisson energy identity。

---

# 37. current/Poisson bridge

production 已保证：

$$
W_{\phi\Delta\rho}
=
-\Delta t\langle E_{pair},J_{charge}\rangle
$$

当 continuity/current 来自同一 x flux。

新增：

$$
R_{transport}
=
W_{\phi\Delta\rho}
+
W_{EJ}.
$$

必须单独输出。

---

# 38. field-particle combined gate

定义：

$$
R_{FP}
=
\Delta K_u^{flux}
+
\Delta U_E
-
W_{electrode}.
$$

这才是 collisionless force/field exchange residual。

---

# 39. full kinetic + boundary transport gate

由于 J1 core 的 one-sided x endpoint 可能携带 kinetic energy：

$$
\Delta K_{state}
=
\Delta K_x^{flux}
+
\Delta K_u^{flux}.
$$

因此完整域能量必须报告：

$$
R_{domain}
=
\Delta K_{state}
+
\Delta U_E
-
W_{electrode}
-
\Delta K_x^{flux}.
$$

J1 core 要求：

```text
R_domain = roundoff
```

`delta_K_x_flux` 在这里就是 test boundary 的显式 kinetic transport source。

---

# 40. failure code 重构

保留现有：

```text
70 configuration
71 initial residual
72 linear/Jv/preconditioner
73 line search
74 not converged
75 legacy total energy
76 negative solution
```

新增：

```text
77 joint_midpoint_source_unit_contract
78 joint_midpoint_cross_work_identity
79 joint_midpoint_current_dual_identity
80 joint_midpoint_transport_poisson_identity
81 joint_midpoint_state_flux_identity
82 joint_midpoint_pairing_field_build
83 joint_midpoint_topology_mismatch
84 joint_midpoint_no_feasible_descent
```

---

# 41. periodic/open topology防误用

在 `advance_joint_midpoint()` 开头新增：

```text
若 background boundary = PERIODIC
且 field solver 使用当前 open/Dirichlet topology：
直接 failure_code=83
```

实际上当前项目只有：

```text
OpenElectrostaticSolver
```

没有 periodic electrostatic solver。

因此 J1 修复完成前：

```text
joint-midpoint-energy + background-x-boundary periodic
```

一律拒绝。

---

# 42. `joint_phase_space_midpoint_energy_test.cpp` 必须修改

删除：

```cpp
left_type = PERIODIC;
right_type = PERIODIC;
```

该 test 不再通过 production background boundary 驱动 x endpoint。

它必须显式启用：

```text
ADJOINT_ONE_SIDED_TEST
```

这是 test-only closure。

如果实现上需要 setter：

```cpp
integrator.set_joint_core_test_boundary_closure(
    JointXBoundaryClosure::ADJOINT_ONE_SIDED_TEST);
```

该 setter：

```text
不得被 main_vpfp.cpp 调用
```

---

# 43. main_vpfp.cpp 的约束

第一轮修复不允许 fp_solver 真正跑新的 joint production boundary。

若用户 CLI 指定：

```text
--background-phase-space-mode joint-midpoint-energy
```

但没有完成 J2 production boundary：

应在启动时明确失败：

```text
joint-midpoint-energy core scheme is not production-open-boundary complete
```

不要静默退回 Strang。

---

# 44. 第一轮不修改 checkpoint

当前 checkpoint 已经存：

```text
background_phase_space_mode
x_transport_velocity_mode
x_transport_velocity_table_schema
```

但在 J2 完成前：

```text
joint mode 不允许 production checkpoint
```

所以第一轮：

```text
禁止修改 vpfp_checkpoint.*
```

等 J2 production boundary 通过后再增加 scheme schema。

---

# 45. nonlinear solver：结构 PASS 前冻结

在 J0 新 cross-work PASS 前：

```text
禁止修改 GMRES
禁止修改 max_iterations
禁止修改 residual tolerance
禁止修改 line search
```

先修方程。

---

# 46. JFNK 当前额外问题

源码当前 GMRES probe 使用：

```cpp
fd =
    sqrt_eps * max(1.0, global_state_norm);
```

这仍存在：

```text
M 有量纲，但与无量纲 1 比较
```

不过 Krylov basis 已归一化，所以 direction denominator 不再是首要问题。

该项推迟到 operator PASS 后处理。

---

# 47. line-search 当前问题

当前 line-search trial 明确调用：

```text
allow_negative_probe=true
```

且只要 residual 下降就：

```text
candidate.swap(trial)
```

所以显著负 trial 会成为下一轮 Newton base state。

这与目前 smoke 的：

```text
min_mass = -1e-3
...
-1e8
```

一致。

J0/J1 operator PASS 后必须修。

---

# 48. feasible Newton gate

step-start：

$$
M_{scale}
=
\max |M^n|.
$$

$$
\tau_{neg}
=
4096\epsilon M_{scale}.
$$

不得：

```text
max(1, ||M||)
```

---

# 49. fraction-to-boundary

对：

```text
d[q] < 0
```

计算：

$$
\lambda_q
=
\frac{M_q+\tau_{neg}}{-d_q}.
$$

$$
\lambda_{pos}
=
0.99\min\lambda_q.
$$

初始：

```text
lambda = min(1, lambda_pos)
```

---

# 50. Armijo

固定：

```text
c1 = 1e-4
max_backtrack = 20
```

trial 必须同时：

```text
finite
min_mass >= -tau_neg
r_trial <= (1-c1*lambda)*r0
Gauss diagnostic not degraded beyond allowed gate
```

否则：

```text
lambda *= 0.5
```

失败：

```text
code=84
```

---

# 51. signed Jv 与 accepted iterate 分离

允许：

```text
JFNK finite-difference probe signed
```

禁止：

```text
accepted nonlinear iterate signed beyond tau_neg
```

代码中必须用不同语义变量，不能继续一个：

```text
allow_negative_probe
```

同时服务两者。

---

# 52. test Maxwellian exact-zero tail

当前 `Species::initialize_maxwellian()` 在 double 下可能产生 exact-zero tail。

在 feasible Newton unit 中使用 test-only floor：

$$
M_{floor}
=
10^{-14}M_{max}.
$$

然后保持 cell/global number重新归一化。

必须输出：

```text
floor_number_change
floor_energy_change
```

要求相对：

```text
< 1e-10
```

该 floor 不进入 production initialization。

---

# 53. 新 J1 运行顺序

只有 J0 全 PASS 后：

```text
J1-A  operator one candidate identity
J1-B  nonlinear single rank
J1-C  nonlinear MPI n=2
```

此阶段不运行：

```text
fp_solver production smoke
```

因为 production open boundary尚未完成。

---

# 54. J1-A：operator identity

不要求 nonlinear solve。

给定：

```text
M_old
M_candidate
fields_old
candidate Poisson fields
```

验证：

```text
M_mid contract
u measure
pairing field
J_charge = G*J_force
u_work = E_pair-J work
Poisson transport identity
```

全部 stable-sum PASS。

---

# 55. J1-B：single-rank nonlinear

要求：

```text
converged=1
finite=1
gauss_ok=1

current_pair_pass=1
transport_poisson_pass=1
state_flux_pass=1
field_particle_energy_pass=1
domain_energy_with_x_boundary_ledger_pass=1

minimum_accepted_iterate >= -tau_neg
```

---

# 56. J1-C：MPI n=2

必须验证：

```text
nonperiodic rank neighbour exchange
shared m_mid trace equal
shared x_flux equal
shared J_charge equal
shared E_pair equal
G* neighbour-cell gather equal
single owner face quadrature
current_pair global roundoff
```

---

# 57. MPI shared face owner

严格沿用 production `FieldParticlePowerAudit`：

```text
shared face 由右侧 rank 在 global face-work sum 中拥有
左 rank 的 local right shared face weight = 0
```

但：

```text
两 rank 的物理 x_flux/current 数值必须相同
```

owner 规则只作用于 ledger sum，不作用于物理 flux continuity。

---

# 58. 新 iteration record

扩展：

```cpp
JointPhaseSpaceIterationRecord
```

新增：

```text
fd_h
lambda_pos
tau_neg
trial_feasible
armijo_pass

current_pair_residual
transport_poisson_residual
state_flux_residual
field_particle_energy_residual
```

---

# 59. J1 result 必须输出

```text
source_commit

scheme=adjoint-center-midpoint-v1

midpoint_contract_pass
u_flux_measure_contract_pass
pairing_field_contract_pass
topology_contract_pass

delta_K_state
delta_K_x_flux
delta_K_u_flux

force_current_work
charge_current_work

current_pair_residual
current_pair_roundoff_bound

potential_charge_work
transport_poisson_residual
transport_poisson_roundoff_bound

field_energy_change
electrode_work
poisson_identity_residual

state_flux_kinetic_residual
field_particle_energy_residual
domain_energy_with_x_boundary_ledger_residual

minimum_accepted_iterate
tau_neg

converged
residual_linf
relative_gauss_residual
failure_code
status
```

---

# 60. 旧 `joint_midpoint_energy_residual` 的兼容语义

旧字段不得在 accepted 后被另一种“domain energy change”覆盖成不同含义。

固定：

```text
joint_midpoint_energy_residual
=
field_particle_energy_residual
```

若还需要：

```text
domain_energy_change
```

使用 ledger 现有字段，不复用同一个 joint 字段。

---

# 61. J0/J1 编译

继续使用当前项目真实 target：

```bash
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx

cmake --build build -j4 --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test
```

---

# 62. 第一轮执行命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

git rev-parse HEAD
git status --short

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target joint_phase_space_midpoint_unit_test

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result ./output/joint_midpoint_j0_v3.result
```

---

# 63. 第一轮只允许做到这里

第一次交给执行模型：

```text
只执行源码修改 + J0
```

禁止：

```text
修改 Newton
修改 GMRES
修改 line search
运行 joint nonlinear J1
运行 fp_solver
进入 Beam/Tail/collision
```

---

# 64. 第一轮 PASS 条件

```text
git_commit == expected or executor reports exact newer commit

source_mass_contract_pass=1
local_midpoint_contract_pass=1
u_flux_measure_contract_pass=1
pairing_field_contract_pass=1

adjoint_dual_face_pass=1
adjoint_current_pair_pass=1
cross_work_pass=1

status=PASS
```

---

# 65. 第一轮 FAIL 后唯一动作

若：

```text
source/unit/midpoint/measure/pairing
```

任何一项 FAIL：

```text
STOP
next_allowed_task=FIX_CURRENT_J0_ONLY
```

若：

```text
legacy cross-work FAIL
adjoint cross-work PASS
```

这是预期根因确认：

```text
next_allowed_task=J1_OPERATOR_INTEGRATION
```

若：

```text
adjoint cross-work 仍 FAIL
```

停止：

```text
next_allowed_task=STOP_REDERIVE_U_TRACE_ADJOINT
```

此时不得让执行模型自己发明 common-corner formula。

---

# 66. 第二轮：J1 operator integration

只有第一轮 PASS 才执行：

```bash
cmake --build build -j4 --target joint_phase_space_midpoint_energy_test

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result ./output/joint_midpoint_j1_v3.result
```

但第二轮先保持现有 nonlinear solver。

只看新的分解 residual。

---

# 67. 第二轮解释规则

如果：

```text
current_pair FAIL
```

停止在 operator。

如果：

```text
current_pair PASS
transport_poisson FAIL
```

检查：

```text
x continuity
pairing_field
face ownership
```

不改 u flux。

如果：

```text
current_pair PASS
transport_poisson PASS
state_flux FAIL
```

检查：

```text
mass residual
x/u flux divergence
state/flux kinetic bookkeeping
```

如果三者 PASS 但 nonlinear 74：

```text
才进入 feasible JFNK 修复。
```

---

# 68. 第三轮：feasible Newton

只在 operator identities PASS 后修改：

```text
fd scale
fraction-to-boundary
Armijo
accepted iterate positivity
```

---

# 69. 不再允许的伪修复

永久禁止：

```text
1. u_flux 继续乘 dx*ring
2. local residual 继续用 candidate 而非 midpoint
3. joint force 继续用 final candidate Ex
4. periodic mass + Dirichlet Poisson
5. 用旧 J0 PASS 宣称 cross-work PASS
6. 修改 Poisson 让 kinetic residual 变小
7. scale E
8. scale J
9. scale a_u
10. scale K
11. cell-wise energy patch
12. global energy patch
13. 放宽 1e-8 energy gate
14. 增大 max_iterations 代替结构修复
15. soft accept
16. 接受 -1e8 Newton iterate
17. converged 后 clip M
18. 把 PPM/FCT 搬进 J1 core
19. 再试 P3-V.2 velocity table作为生产修复
20. 从 rho delta 反推第二套 current
21. 修改 endpoint dx/2
22. 修改 OpenPoissonWorkIdentity
23. 修改 build_potential_pairing_field 以配合 joint
24. J1 core 未 PASS 就跑 production open boundary
```

---

# 70. 当前 production open boundary：明确 BLOCK

源码确认：

```text
OpenBackgroundBoundary
```

production reservoir/absorbing boundary 是基于：

```text
incoming/outgoing characteristic
```

处理的。

而 J1 core 的：

```text
ADJOINT_ONE_SIDED_TEST
```

不是 production boundary。

因此：

```text
J2 production open boundary
```

在本 v3 中明确：

```text
BLOCKED
```

直到 J1 core PASS。

原因不是缺代码，而是 boundary energy-compatible bracket 需要单独把：

```text
reservoir incoming current
absorbing outflow current
kinetic boundary flux
Poisson endpoint current
electrode work
```

放进同一个开放系统能量推导。

不允许 GPT-5.5/5.6 Luna在没有下一份边界推导文档时自行实现。

---

# 71. Beam/Tail/collision 同样 BLOCK

在 J2 production open boundary PASS 前：

```text
Beam BLOCKED
Tail BLOCKED
conversion BLOCKED
return BLOCKED
collision BLOCKED
```

---

# 72. default production path 回归

每轮代码变更必须确认：

```text
background_phase_space_mode default = STRANG_PPM
```

并运行已有至少：

```text
vpfp_x_transport_flux_audit_test
vpfp_x_u_power_pairing_test
vpfp_poisson_work_identity_test
```

这些旧测试不能因为 joint 修复改变结果。

---

# 73. checkpoint 第一阶段不改

第一轮/第二轮：

```text
checkpoint files unchanged
checkpoint schema unchanged
```

只有 joint production open boundary完成后，才增加：

```text
joint_phase_space_scheme
scheme_schema
```

---

# 74. 执行报告模板

每一轮执行模型只按以下模板报告：

```text
source_commit=
working_tree_before=

stage=
status=

changed_files=

build_status=

source_mass_contract=
local_midpoint_contract=
u_flux_measure_contract=
pairing_field_contract=
topology_contract=

legacy_current_pair_residual=
legacy_current_pair_relative=

adjoint_current_pair_residual=
adjoint_current_pair_roundoff_bound=

dual_face_linf=
dual_face_roundoff_bound=

u_force_work=
charge_current_work=
cross_work_residual=
cross_work_roundoff_bound=

transport_poisson_residual=
transport_poisson_roundoff_bound=

state_flux_kinetic_residual=
state_flux_kinetic_roundoff_bound=

field_particle_energy_residual=
field_particle_energy_roundoff_bound=

nonlinear_residual=
relative_gauss_residual=

minimum_accepted_iterate=
tau_neg=

default_strang_regression=

first_failed_gate=
first_failed_rank=
first_failed_ix=
first_failed_j=
first_failed_k=
first_failed_value=

next_allowed_task=
```

---

# 75. 执行模型不得输出推测结论

禁止：

```text
probably
likely
maybe
seems
应该
可能
大概
看起来像
```

只能：

```text
PASS
FAIL
BLOCKED
```

以及实际数值。

---

# 76. 本次修复真正的第一目标

不是：

```text
让 J1 smoke 跑完
```

不是：

```text
让 combined energy residual 数值变小
```

第一目标是同时证明：

$$
\boxed{
J_{charge}=G^*J_{force}
}
$$

以及：

$$
\boxed{
\Delta K_u
=
\Delta t
\langle E_{pair},J_{charge}\rangle_f
}
$$

到 stable-sum roundoff。

---

# 77. 本次修复真正的第二目标

只有第一目标 PASS 后才要求：

$$
\boxed{
\Delta K_u
+
\Delta U_E
-
W_{electrode}
=
O(\epsilon_{sum})
}
$$

---

# 78. 对此前 v2 的明确撤销项

以下 v2 内容作废：

```text
1. 不再直接实现 arbitrary bilinear common-corner mobility
2. 不再用 ΔPhi 代替 production E_pair
3. 不再假设 E_pair=-simple cell potential difference/dx at endpoints
4. 不再把 periodic J1 视为合法 topology
5. 不再在 M 已含 ring/dx 的前提下再次构造 g=M/(dx*du) 后随意乘 ring
6. 不再把 code-75 直接视作“center trace 已被证明结构失败”
```

v3 优先使用：

```text
当前 u center trace 的精确离散伴随
+
production G/G*
+
production exact-dual E_pair
```

这是对当前源码改动最小、证据最强的路径。

---

# 79. 为什么不直接采用旧 P3-V.2 velocity

旧 P3-V.2 / 当前 `vx_energy_conjugate_cell` 使用：

```text
center-to-center distance
```

构造 cell velocity。

而本 v3 的 `v_dagger` 是从当前 u-face center trace **反向做代数转置**得到：

$$
v^\dagger_j
=
\frac{
\Delta K^-_j+\Delta K^+_j
}{
2m_ec\,\Delta u_j
}.
$$

非均匀网格上：

```text
upar_width[j]
```

与：

```text
u[j+1]-u[j-1]
```

不是同一个量。

所以 v3 不是重复 P3-V.2。

---

# 80. 交给 GPT-5.5 / GPT-5.6 Luna 的最终执行指令

```text
你不是来重新分析算法的。

仓库基线：
mumuni123/GPT-fokker-planck-equation-solver
main
ad09bf2cb52f92b21a3a3cb71c75b6d798b3accd

首先核对实际 checkout commit。

第一轮只修改：
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
tests/joint_phase_space_midpoint_unit_test.cpp
必要的 CMake test wiring

第一轮必须修：
A1 midpoint
A2 u-flux measure
并新增 adjoint velocity + direct cross-work audit。

第一轮不得修改：
Newton
GMRES
line search
Poisson
PPM/FCT
Beam/Tail/collision
checkpoint

第一轮只运行 J0。

只有：
J_charge = G*J_force
以及
DeltaK_u = dt<E_pair,J_charge>
达到 stable-sum roundoff，
才进入第二轮 integrator field-pairing 接线。

任何 Gate FAIL：
立即停止在该 Gate。
不要自行更换方案。
```

---

# 81. 当前项目状态在执行前应记录为

```text
source_root_cause_A1_midpoint_mismatch=CONFIRMED
source_root_cause_A2_u_measure_duplication=CONFIRMED
source_root_cause_A3_wrong_force_field_object=CONFIRMED
source_root_cause_A4_topology_mismatch=CONFIRMED
source_root_cause_A5_J0_missing_cross_gate=CONFIRMED

structural_nonuniform_u_adjoint_mismatch=TO_BE_RETESTED
joint_core_energy_closed=NOT_YET
production_open_boundary_joint=BLOCKED
beam_tail_collision_joint=BLOCKED
default_strang_ppm=UNCHANGED
```

---

# 82. 最终停止条件

本 v3 文档的执行到以下状态后停止：

```text
J0 source contracts PASS
J0 adjoint cross-work PASS
J1 operator energy chain PASS
J1 nonlinear single-rank PASS
J1 MPI PASS
```

然后输出：

```text
next_required_document=
VPFP_joint_open_boundary_energy_compatible_derivation.md
```

不要自行进入 production open-boundary 实现。

