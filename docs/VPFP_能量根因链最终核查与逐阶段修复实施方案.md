# VPFP 能量根因链最终核查与逐阶段修复实施方案

## 文档编排说明（本次整理新增）

本次整理**不删除、不改写、不替换**本文既有任何技术结论、公式、代码片段、阶段门、禁止项或历史记录。原有内容仍完整保留在下方；本节以及文末“执行导航与命令附录”仅新增：

1. 原有章节之间的依赖关系；
2. 阶段执行时应当阅读的原文位置；
3. 当前已经存在的测试目标的可执行命令；
4. 未来 J2/J3 接入测试必须提供的命令契约。

因此，发生冲突时，原有 F0--F11、J1--J3、禁止项和阶段门的正文优先级高于本次新增的导航文字。命令中的输出目录均使用新的时间戳目录，**不得**使用删除旧输出或覆盖旧结果的方式复跑。

## 执行导航（本次整理新增）

本文原有内容可按三个层次阅读和执行：

| 层次 | 原有章节 | 目的 | 前置条件 | 进入下一层的必要条件 |
|---|---|---|---|---|
| 物理约束与根因证据 | §1--§9、§21、§24--§28 | 固定不可改变的 VPFP-PIC 物理模型，说明 J1 之前必须排除的确定性实现错误 | 无 | 已理解并接受“不修改 Poisson/Beam/Tail/碰撞物理”的约束 |
| 联合 bulk 核心闭合 | F0--F10、§11、§12--§13 | 先清除 J1 fixture、时间层、几何量纲、转置速度与 pairing field 错误，再证明 J0/J1 的三重离散功关系 | 上一层 | 单 rank 与 1/2/5 rank J1 均通过；若能量代数通过但 Newton 不通过，才进入 F11 |
| 生产物理逐项接回 | F11、J2、J3a--J3e、§20 | 在 bulk 核心闭合后，依次接入开放 background、Beam、Tail、转换、return、碰撞 | J1 MPI 门通过 | 每一项独立能量账通过后才允许接入下一项 |

**实际执行顺序**必须仍严格遵从原文 §23；本导航不授权跳过任何阶段。尤其 F1、F3、F4、F5、F6 是当前已从源码直接确认的修复顺序，F2/F7/F9 是对应的测试与诊断支撑，F11 不是 code 75 的提前修复手段。

## 原始审计结论（以下内容保持原文，未改写）

我已经把仓库当前 `main` 分支的源码与 `docs` 下除 `archive` 外的文档交叉核对了一遍，重点检查了 `joint_phase_space_midpoint.*`、`vpfp_integrator.*`、`species.*`、`grid.h`、`open_electrostatic_solver.*`、`field_particle_power_audit.*`、J0/J1 测试和 CMake 测试入口。

结论比现有 `J1联合中点当前失败分析_2026-08-21.md` 还要更具体：**原先对旧 Strang/PPM 路线的“场—粒子离散功不同源”根因判断仍然成立，但目前 J1 的** **`1.1132815e14 J/m²`** **不能再直接视为“联合中点格式本身仍然不保能”的干净证据，因为当前 J1 实现和测试中至少还存在 5 个必须先修掉的确定性缺陷。**

其中最严重的三个是：

1. **J1 测试初始化本身有错误**：初始化 Maxwellian 后没有调用 `compute_moments()`，就把仍为全零的 `electrons.number_density` 复制给 `ion_density`。初始 Poisson 因而看到零电荷，而第一次 Newton candidate 又用真实电子质量、零离子背景重算 Poisson，相当于一步内突然出现整个背景电子的净负电荷。当前 `1.113e14 J/m²` 因此受到测试夹具错误污染。
2. **`evaluate_local_residual()`** **根本没有使用真正的时间中点质量**。全局测试辅助函数 `evaluate_residual()` 会构造 `0.5*(M_old+M_candidate)`，但生产 J1 用的 `evaluate_local_residual()` 直接用 `M_candidate` 做 MPI ghost、x trace 和 u trace。这已经违反了联合中点方案自己的核心定义。
3. **u 方向通量存在明确的重复几何因子错误**。`Species::f` 在 cylindrical 模式保存的已经是完整 cell-integrated mass，初始化时已经乘过 `dx × Δu × ring`；但 J0/J1 又在 `M/Δu` 后额外乘了一次 `dx × ring`。因此当前 u 通量不仅不是正确量纲，而且会随 `dx` 和 `u_perp` ring 被错误二次缩放。

另外还有两个结构性缺陷：

4. J1 Newton 中实际传给 u-force 的是 candidate 的最终 `Ex`，不是已经由 Poisson 功恒等式验证过的离散梯度 pairing field。
5. 当前 `build_hamiltonian_velocity()` 在非均匀 `u_parallel` 网格上并不是 u-face 中心 trace 的严格离散转置，因此即便时间中点、量纲和场层全部修正，x 电流和 u 动能功仍未必严格相等。当前 J0 的能量测试又恰好只用“同一个 u 通量”验证自身求和分部关系，没有真正检验 x-current 与 u-work 的交叉配对，所以 J0 PASS 不能排除这个问题。

下面这份方案我按你要求写成了**可以直接交给 GPT-5.6 Luna / DeepSeek-v4-Flash 执行的机械式实施规格**。核心原则是：**执行模型不得自行推理、不得换算法、不得“顺手优化”、不得根据测试结果自行发明修复；每一步只执行明确写出的内容。**

---

# VPFP 能量根因链最终核查与逐阶段修复实施方案

## 0. 文档用途与最高执行纪律

本方案用于修复当前仓库 `Fokker-PlanckEquationSolver` 的离散能量闭合问题。

执行本方案的自动编码模型必须把本文视为**强制执行规格**，不得把它理解为“建议”“参考思路”或“可以自行改进的方案”。

### 0.1 绝对禁止自行思考

执行模型必须遵守以下规则：

1. **不得自行改变本文给出的算法。**
2. **不得自行选择另一种数值格式。**
3. **不得为了让测试通过而放宽阈值。**
4. **不得为了让能量残差变小而加入任何能量补丁。**
5. **不得自行修改 Poisson 方程。**
6. **不得自行修改 Beam 物理模型。**
7. **不得自行修改 Tail 物理模型。**
8. **不得自行修改碰撞物理模型。**
9. **不得自行改变边界条件。**
10. **不得删除失败测试。**
11. **不得把失败测试改成始终 PASS。**
12. **不得把 NaN、Inf、负质量或不守恒状态通过裁剪隐藏。**
13. **不得因为“认为某修改更好”而跳过本文步骤。**
14. **不得同时执行两个尚未通过阶段门的阶段。**
15. 当前阶段失败时，只允许检查和修改当前阶段明确允许修改的文件。
16. 当前阶段没有 PASS 前，禁止进入下一阶段。
17. 除本文明确要求外，不得重构无关代码。
18. 不得更换现有生产 Poisson 离散。
19. 不得恢复周期电场。
20. 不得恢复旧 Ampere 场推进。
21. 不得进行全局零模扣除。
22. 不得强制平均电场为零。
23. 不得改变电子束的开放注入和开放流出。
24. 不得把 background bulk 和 background tail 当作两个独立物种。

如果执行模型遇到本文没有说明的问题：

> **停止当前阶段，报告实际源码、错误信息、测试数值，不得猜测修复。**

---

# 1. 不可改变的最终物理模型

本节优先级高于后面所有数值修改。

仓库最终物理模型必须继续严格遵循：

`docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md`

该文档明确规定：bulk 与 tail 是同一背景电子物种的两种数值表示；转换只是内部表示变换；场必须保持非周期 Gauss/Poisson；Beam 保持开放注入和开放流出；禁止全局能量补丁、人工电流缩放、强制平均电场为零以及无记账裁剪。

## 1.1 背景电子

物理背景电子分布始终是

$$
f_e=f_{\mathrm{bulk}}+f_{\mathrm{tail}}.
$$

其中：

- bulk：Eulerian 1D2V；
- tail：PIC 3V；
- 二者属于**同一个背景电子物种**；
- bulk-to-tail 和 tail-to-bulk 是内部表示变化。

必须保持

$$
N_e=N_{\mathrm{bulk}}+N_{\mathrm{tail}},
$$

$$
P_{e,x}=P_{\mathrm{bulk},x}+P_{\mathrm{tail},x},
$$

$$
K_e=K_{\mathrm{bulk}}+K_{\mathrm{tail}}.
$$

Poisson 中必须使用

$$
n_e=n_{\mathrm{bulk}}+n_{\mathrm{tail}}.
$$

该单位和统一矩定义已经在主重构方案中明确规定。

---

## 1.2 Beam

Beam 必须继续作为独立外部注入电子物种。

Beam：

- 从左边界注入；
- 注入时间保持原方案；
- 使用原有 D-K-D；
- 从任意物理边界离开后删除；
- 不允许周期回卷；
- 不得为了能量闭合修改 Beam 注入统计、粒子权重、随机数序列或边界条件。

---

## 1.3 场

生产场必须继续使用：

- 非周期静电 Gauss/Poisson；
- 生产边界为 `DIRICHLET_PHI`；
- 当前生产值保持

```text
phi_left  = 0
phi_right = 0
```

不得：

- 改成周期 Poisson；
- 改用 Ampere 推进；
- 扣掉平均电荷；
- 扣掉平均电场；
- 修改 Poisson stencil 来迁就粒子算法。

现有 `OpenElectrostaticSolver` 的 Poisson 功恒等式和 pairing-field 构造已经存在，并且正是后续联合离散应该复用的对象。

---

## 1.4 空间边界

生产背景 bulk：

- reservoir/open boundary；
- 出流自然离开；
- 入流由 reservoir 分布给出；
- reservoir 交换进入数目、动量和能量账本。

Tail：

- 开放流出；
- 删除前必须完成域内轨迹、电流和出流账本。

Beam：

- 开放注入；
- 开放流出。

当前 J0/J1 中出现的 periodic x 只能作为**隔离数值测试拓扑**存在。

**绝对不得把 periodic x 变成生产物理边界。**

生产开放边界只能在 J2 阶段接入。主物理文档已经明确生产空间边界和背景 reservoir 语义。

---

# 2. 当前能量问题的完整根因链

## 2.1 第一层：旧 Strang 路线的历史根因仍然成立

旧算法把

```text
x 输运
u_parallel 场加速
Poisson
```

作为相互独立的离散算子构造。

因此：

- x 输运产生一个连续性电流；
- u 加速产生另一个动能功；
- Poisson 又有自己的离散功配对；
- 三者没有来自同一个离散 Hamiltonian/phase-space flux。

所以旧路线并不存在必然成立的离散恒等式

$$
\Delta K_{\mathrm{bulk}}+\Delta U_E-W_{\mathrm{electrode}}=0.
$$

这一历史根因判断仍然正确。

因此：

> **不得返回修补旧 Strang/PPM/FCT 的路线。**

联合相空间离散的方向是正确的。

---

# 3. 当前 J1 结果首先存在测试夹具错误

当前 `joint_phase_space_midpoint_energy_test.cpp` 的顺序是：

```cpp
electrons.initialize_maxwellian(0.0);

std::vector<double> ion_density = electrons.number_density;
```

但 `initialize_maxwellian()` 不会自动更新 `number_density`。

`Species::init()` 创建的 `number_density` 初始是零数组，而真正的 density 只有调用

```cpp
electrons.compute_moments();
```

之后才会由 `Species::f` 计算。

因此当前 J1 测试实际上执行了：

```text
真实 f_e != 0

但

electrons.number_density == 0
ion_density == 0
```

随后初始场又调用：

```cpp
fields.set_charge_density(...)
```

而该函数只是读取 `electrons.number_density`，不会重新计算 moments。

于是初始 Poisson 实际看到

$$
\rho^n=0.
$$

但是进入 Newton candidate 后，`advance_joint_midpoint()` 又直接由真实 `M` 求电子密度：

$$
n_e^{\mathrm{candidate}}=\frac{1}{\Delta x}\sum_{j,k}M_{i,j,k}.
$$

而固定离子仍然为零。

于是 candidate 突然变为

$$
\rho^{\mathrm{candidate}}=-e n_e.
$$

这不是所设计的平滑背景 VP 测试。

因此当前文档记录的

```text
energy_residual = 1.1132815186046597e14
failure_code    = 75
```

**不能再作为联合中点格式结构缺陷的定量证据。**

这里必须修正文档解释：

> 旧 Strang 不闭合的根因结论仍成立；但是当前 J1 的 `1.113e14 J/m²` 数值，在测试初始化修复前属于“受错误初态污染的结果”。

---

# 4. 当前 J1 实现没有真正使用时间中点质量

联合格式要求

$$
M^{n+1/2}=\frac{M^n+M^{n+1}}{2}.
$$

仓库中的全局测试函数 `evaluate_residual()` 确实这样做。

但是生产 J1 使用的是 `evaluate_local_residual()`。

当前 `evaluate_local_residual()`：

- MPI ghost 来自 `m_candidate_local`；
- x trace 来自 `m_candidate_local`；
- u trace来自 `m_candidate_local`。

也就是说当前实际上在用

$$
M^{n+1}
$$

构造通量，而不是

$$
M^{n+1/2}.
$$

源码已经能直接确认这一点。

因此目前所谓 `joint-midpoint-energy` 在生产 J1 residual 中实际上并非真正的 joint midpoint。

---

# 5. 当前 u 通量存在重复几何权重错误

这是本轮源码审计中新确认的确定性错误。

## 5.1 `Species::f` 的真实定义

背景 cylindrical 模式中存储的是 cell-integrated mass。

即

$$
M_{i,j,k}=\bar f_{i,j,k}\Delta x\Delta u_j A_k,
$$

其中

$$
A_k=\pi\left(u_{\perp,k+1/2}^2-u_{\perp,k-1/2}^2\right).
$$

源码初始化本身就是：

```cpp
f3 * sgrid->dx * cgrid.cell_phase_volume(iv, imu)
```

而

```cpp
cell_phase_volume = Δu_parallel * uperp_ring_area
```

因此 `M` 已经包含：

```text
dx
Δu_parallel
uperp ring area
```

这与头文件关于 cell-integrated mass 的契约一致。

---

## 5.2 正确的 u-face trace

定义

$$
\bar M_{i,j,k}=\frac{M^n_{i,j,k}+M^{n+1}_{i,j,k}}{2}.
$$

u-face 的中心 trace 应为

$$
T^u_{i,j+1/2,k}=\frac12\left(\frac{\bar M_{i,j,k}}{\Delta u_j}+\frac{\bar M_{i,j+1,k}}{\Delta u_{j+1}}\right).
$$

注意：

$$
\frac{M}{\Delta u}
$$

已经具有

```text
dx × uperp ring
```

因子。

因此不能再乘一次。

---

## 5.3 正确 u 通量

定义

$$
a_i^{n+1/2}=\frac{q_e E_i^{n+1/2}}{m_ec}.
$$

则正确质量通量率为

$$
F^u_{i,j+1/2,k}=a_i^{n+1/2}T^u_{i,j+1/2,k}.
$$

即代码必须是：

```cpp
u_flux_rate = a * ftrace;
```

而不是当前的：

```cpp
u_flux_rate =
    a * ftrace * sg.dx * vg.uperp_ring_areas[k];
```

当前代码在 `build_periodic_center_flux()` 和 `evaluate_local_residual()` 中都存在这一重复乘法。

---

# 6. 当前 u-work 与 x-current 仍然没有严格共轭

这是修完前面三个错误后仍必须修的结构问题。

对于内部 u-face，定义

$$
\Delta K_{j+1/2,k}=K_{j+1,k}-K_{j,k}.
$$

u 通量导致的离散动能变化为

$$
\Delta K_u=\Delta t\sum_{i,k}\sum_{j+1/2}\Delta K_{j+1/2,k}F^u_{i,j+1/2,k}.
$$

代入中心 trace：

$$
\Delta K_u=\Delta t\sum_i E_i^{\mathrm{pair}}q_e\sum_{j,k}v^H_{j,k}\bar M_{i,j,k}.
$$

为了使这个式子**逐代数恒等成立**，`vH` 不能再凭经验选择，而必须是 u trace 的离散转置。

## 6.1 必须使用的 J1 Hamiltonian velocity

初始化所有

```cpp
vH[j,k] = 0
```

随后遍历所有内部 u-face。

对于每一个

```text
jface = 1 ... Nv-1
```

定义：

$$
\Delta K_{j-1/2,k}=K_{j,k}-K_{j-1,k}.
$$

然后执行：

$$
v^H_{j-1,k}\mathrel{+}=\frac{\Delta K_{j-1/2,k}}{2m_ec\Delta u_{j-1}},
$$

$$
v^H_{j,k}\mathrel{+}=\frac{\Delta K_{j-1/2,k}}{2m_ec\Delta u_j}.
$$

对应 C++ 算法必须按下面结构实现：

```cpp
std::vector<double> velocity(nv * nmu, 0.0);

for (int jf = 1; jf < nv; ++jf) {
    const int jl = jf - 1;
    const int jr = jf;

    for (int k = 0; k < nmu; ++k) {
        const double delta_k =
            K[jr, k] - K[jl, k];

        velocity[jl, k] +=
            delta_k /
            (2.0 * me * c * du[jl]);

        velocity[jr, k] +=
            delta_k /
            (2.0 * me * c * du[jr]);
    }
}
```

### 禁止

不得在本阶段使用：

```cpp
vg.vx
```

不得使用：

```cpp
vg.vx_energy_conjugate_cell
```

不得使用解析式：

```cpp
c * u_parallel / gamma
```

不得继续使用当前：

```text
(K[j+1]-K[j-1]) /
(me*c*(u[j+1]-u[j-1]))
```

原因不是这些速度“物理上错误”，而是它们**不是当前 u-face trace 的严格离散转置**。

尤其当前网格是非均匀 `u_parallel` 网格，二者一般不相等。`grid.h` 中现有 energy-conjugate 表也有自己的非均匀投影定义，不能直接拿来替代 J1 所需要的转置。

---

# 7. x 通量必须和上述速度完全同源

定义

$$
F^x_{i+1/2,j,k}=v^H_{j,k}\frac{\bar M_{i,j,k}+\bar M_{i+1,j,k}}{2\Delta x}.
$$

电荷电流必须直接由该**同一个生产 x 通量**得到：

$$
J^{\mathrm{charge}}_{i+1/2}=q_e\sum_{j,k}F^x_{i+1/2,j,k}.
$$

禁止重新根据 distribution moments 构造另一个 J。

禁止调用 `Species::current_x` 来替代这个 current。

禁止用解析速度重新积分 current。

---

# 8. force current 与 charge current 的离散共轭关系

定义 cell force current：

$$
J_i^{\mathrm{force}}=\frac{q_e}{\Delta x}\sum_{j,k}v^H_{j,k}\bar M_{i,j,k}.
$$

内部空间 face 应满足

$$
J^{\mathrm{charge}}_{i+1/2}=\frac12\left(J_i^{\mathrm{force}}+J_{i+1}^{\mathrm{force}}\right).
$$

这不是额外物理假设。

它应该由上述 x 中心 flux 自动得到。

因此必须增加独立测试验证这个关系，而不能“如果不相等就平均修正”。

**绝对禁止修改** **`J_charge`** **来强迫它等于** **`J_force`****。**

正确做法是：

> 两者不相等则说明 flux 构造错误，测试失败，停止。

---

# 9. J1 必须使用 Poisson 的离散 pairing field

当前 Newton candidate 在 Poisson 解完后直接使用：

```cpp
eval_fields.Ex
```

作为 u force field。

这一行为必须修改。

现有 `OpenElectrostaticSolver` 已经提供：

```cpp
build_potential_pairing_field(
    before,
    after,
    pairing_face,
    mpi_rank,
    mpi_size
)
```

它由与 Poisson 能量恒等式相同的 cell-average potential 构造，并且已经正确处理：

- 非周期 Dirichlet 电势；
- endpoint half-weight；
- MPI shared face；
- 电极功。

源码中已经明确说明其目标恒等式。

因此 J1 必须复用它。

---

# 10. 完整修复步骤

以下阶段必须严格顺序执行。

---

# 阶段 F0：冻结物理模型和修改范围

## F0.1 本阶段不修改代码

只检查工作树。

允许读取任何文件。

禁止修改任何文件。

## F0.2 记录 Git 状态

执行：

```bash
git status --short
git rev-parse HEAD
```

记录：

```text
baseline_commit=
working_tree_clean=
```

如果工作树已有修改：

- 不得删除；
- 不得 reset；
- 不得覆盖；
- 报告已有修改。

---

## F0.3 初始允许修改文件白名单

F1 至 F8 只允许修改：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
src/vpfp_integrator.h
src/vpfp_integrator.cpp
tests/joint_phase_space_midpoint_unit_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp
CMakeLists.txt
```

其中 `CMakeLists.txt` 只有增加测试 target 时才允许修改。

本轮禁止修改：

```text
src/open_electrostatic_solver.*
src/conservative_ppm_remap.*
src/beam_pic.*
src/background_tail_pic.*
src/bulk_tail_converter.*
src/tail_bulk_return.*
src/cylindrical_fp_collision.*
src/background_tail_collision*
src/hybrid_collision_step.*
src/open_boundary.*
```

特别是：

> **不得修改** **`OpenElectrostaticSolver`****。**

---

# 阶段 F1：先修 J1 测试初始化

这是所有算法判断之前的强制阶段。

## F1.1 修改现有 `smooth-background`

在：

```text
tests/joint_phase_space_midpoint_energy_test.cpp
```

找到：

```cpp
electrons.initialize_maxwellian(0.0);
```

紧接着添加：

```cpp
electrons.compute_moments();
```

然后才允许执行：

```cpp
std::vector<double> ion_density = electrons.number_density;
```

顺序必须是：

```cpp
electrons.initialize_maxwellian(0.0);
electrons.compute_moments();

std::vector<double> ion_density = electrons.number_density;
```

不得反过来。

---

## F1.2 增加初始化防回归检查

在第一次 Poisson 之前计算：

```text
electron_density_max
electron_density_min
ion_density_max
ion_density_min
initial_rho_linf
initial_E_linf
initial_gauss_linf
```

`smooth-background` 的目的就是离散中性 Maxwellian。

因此至少要求：

```text
electron_density_max > 0
ion_density_max > 0
```

并且离子密度与离散电子密度相同至浮点精度。

---

## F1.3 此时只运行 J1 一次

不要修改任何联合算子。

重新编译：

```bash
make -j4 joint_phase_space_midpoint_energy_test
```

运行现有 J1。

记录：

```text
accepted
converged
residual_linf
poisson_residual_linf
energy_residual
failure_code
initial_rho_linf
initial_E_linf
```

### 本阶段非常重要的判断规则

如果仅修这个初始化以后，原来的

```text
1.1132815186046597e14
```

消失：

> 将原数字标记为“历史错误测试结果”。

不得继续在文档中把它称为已经确认的纯联合离散能量误差。

即使此时 J1 PASS，也不能结束修复，因为 `smooth-background` 是接近零场的平凡情况。

---

# 阶段 F2：增加非平凡 J1 场测试

不能让能量算法只由零场 Maxwellian 验证。

## F2.1 新增测试 case

增加：

```text
smooth-perturbed-background
```

电子密度设置为：

$$
n_e(x_i)=n_0\left[1+10^{-4}\cos\left(2\pi\frac{i_g+1/2}{N_x}\right)\right].
$$

其中：

```text
ig = grid.ix_start + ix
```

离子保持：

$$
n_i=n_0.
$$

实现必须：

```cpp
std::vector<double> electron_profile(grid.nx_local);

for (...) {
    electron_profile[ix] =
        Param::dens *
        (1.0 + 1.0e-4 *
         std::cos(
             2.0 * Const::pi *
             (static_cast<double>(ig) + 0.5) /
             static_cast<double>(grid.nx_global)));
}

electrons.initialize_maxwellian_profile(
    electron_profile, 0.0);

electrons.compute_moments();

std::vector<double> ion_density(
    grid.nx_local, Param::dens);
```

然后正常构造 Poisson。

---

## F2.2 此测试不得替换原测试

保留：

```text
smooth-background
```

再新增：

```text
smooth-perturbed-background
```

以后：

- equilibrium 测试检查零场平衡；
- perturbed 测试检查真正的 field-particle energy exchange。

---

# 阶段 F3：修正 u 通量量纲

## F3.1 修改 `build_periodic_center_flux()`

找到：

```cpp
bundle.u_flux_rate[...] =
    a * f_trace * sg.dx * vg.uperp_ring_areas[k];
```

修改为：

```cpp
bundle.u_flux_rate[...] =
    a * f_trace;
```

除此以外本步不得修改 flux 算法。

---

## F3.2 修改 `evaluate_local_residual()`

找到同样形式：

```cpp
bundle.u_flux_rate[...] =
    a * ftrace * sg.dx * vg.uperp_ring_areas[k];
```

修改为：

```cpp
bundle.u_flux_rate[...] =
    a * ftrace;
```

---

## F3.3 暂时不要修其他东西

本阶段不允许：

- 修改 Hamiltonian velocity；
- 修改 E pairing；
- 修改 Newton；
- 修改 Poisson；
- 修改 positivity。

---

## F3.4 新增 u-flux 单位测试

测试直接构造一个常数 phase density：

$$
M_j=f_0\Delta x\Delta u_j A_k.
$$

则

$$
\frac{M_j}{\Delta u_j}=f_0\Delta x A_k.
$$

对于左右状态相同：

$$
F^u=a f_0\Delta x A_k.
$$

测试必须比较生产 `u_flux_rate` 与这个解析结果。

至少测试：

```text
dx_1 != dx_2
ring_k1 != ring_k2
```

目的：

- 如果以后误加第二个 `dx`，测试失败；
- 如果以后误加第二个 ring，测试失败。

---

## F3.5 F3 阶段门

必须 PASS：

```text
joint_phase_space_midpoint_unit_test
```

且新 geometry test 相对误差达到浮点舍入量级。

失败时禁止进入 F4。

---

# 阶段 F4：修复真正的时间中点质量

## F4.1 修改 `evaluate_local_residual()`

函数开始验证完输入后，构造：

```cpp
std::vector<double> m_mid_local(count, 0.0);

for (size_t i = 0; i < count; ++i) {
    m_mid_local[i] =
        0.5 *
        (m_old_local[i] +
         m_candidate_local[i]);
}
```

必须检查 finite。

---

## F4.2 MPI ghost 必须来自 `m_mid_local`

当前所有从：

```cpp
m_candidate_local
```

交换空间 ghost 的代码改成：

```cpp
m_mid_local
```

单 rank 情况同样如此。

---

## F4.3 x trace 只能读取 midpoint

将 x trace 中所有：

```cpp
m_candidate_local[...]
```

改成：

```cpp
m_mid_local[...]
```

以及 midpoint ghost。

---

## F4.4 u trace 只能读取 midpoint

必须变为：

```cpp
const double ftrace =
    0.5 *
    (
        m_mid_local[left] /
            vg.upar_widths[jf - 1]
        +
        m_mid_local[right] /
            vg.upar_widths[jf]
    );
```

禁止读取 candidate 本身。

---

## F4.5 candidate 只能继续用于两处

`m_candidate_local` 此后只允许用于：

1. residual 左侧：

```cpp
candidate - old - delta
```

2. 非负性检查。

不得再直接用于 phase-space flux。

---

## F4.6 新增 local/global 一致性测试

单 rank 情况下，对完全相同：

```text
M_old
M_candidate
E
dt
```

分别调用：

```cpp
evaluate_residual()
evaluate_local_residual()
```

要求：

```text
x_flux_rate
u_flux_rate
mass_delta_x
mass_delta_u
mass_delta_total
residual
```

全部在浮点 roundoff 内一致。

如果不同：

> F4 FAIL。

不得进入下一阶段。

---

# 阶段 F5：让 u-work 定义自己的严格共轭速度

这是本轮核心离散修复。

## F5.1 重写 `build_hamiltonian_velocity()`

严格使用第 6 节给出的 face-to-cell transpose 算法。

不得保留当前 centered derivative 算法作为生产分支。

函数结果仍然可以叫：

```cpp
build_hamiltonian_velocity
```

这样可以减少接口修改。

---

## F5.2 必须在函数注释中写清楚

添加注释，意思必须包括：

```text
This velocity is not an analytic cell velocity.
It is the exact algebraic transpose of the J1 centered
u_parallel face trace under the cell-integrated-mass
representation.

Do not replace it with vg.vx,
vg.vx_energy_conjugate_cell,
or a centered derivative on the nonuniform grid.
```

---

## F5.3 增加严格 u-energy adjoint 测试

构造任意正的 deterministic `M_mid`。

不要使用 Maxwellian 对称性来隐藏错误。

例如让每个 cell 质量依赖：

```text
ix
j
k
```

但保持全部正值。

给定非零 `E_cell`。

生产代码计算：

$$
W_u=\Delta t\sum_{i,j+1/2,k}\Delta K_{j+1/2,k}F^u_{i,j+1/2,k}.
$$

独立按 `vH` 计算：

$$
W_J=\Delta t\sum_i E_i q_e\sum_{j,k}v^H_{j,k}M_{i,j,k}.
$$

要求：

$$
W_u-W_J
$$

只剩浮点求和误差。

测试不能调用两个完全相同的 helper 计算左右两边。

必须：

- 左边从实际 `u_flux_rate` 得到；
- 右边从实际 `vH` 和 `M_mid` 独立求和。

---

## F5.4 必须检查非均匀网格

测试必须使用生产 `CylindricalVelocityGrid` 的非均匀 `u_parallel` 网格。

不得只用 uniform grid。

---

## F5.5 额外不变量

要求 `vH`：

- 全 finite；
- 对称网格上关于 `u_parallel` 反对称；
- 不得出现明显超过光速的值。

但不得为了满足光速限制而 clip。

若超过：

> FAIL 并报告具体 `j,k,vH`。

---

# 阶段 F6：J1 改用 Poisson pairing field

## F6.1 保持 candidate Poisson 求解不变

继续：

1. 从 candidate `M` 求电子密度；
2. 与固定离子组成 `rho`;
3. 调用现有 `field_solver_.solve()`。

不得修改这部分 Poisson 算法。

---

## F6.2 Poisson 后构造 pairing face

在 `advance_joint_midpoint()` 的 candidate `evaluate` lambda 内，在 Poisson solve 后调用：

```cpp
std::vector<double> pairing_face;

const bool pairing_ok =
    field_solver_.build_potential_pairing_field(
        fields,
        eval_fields,
        pairing_face,
        mpi_rank,
        mpi_size);
```

这里：

```text
fields
```

必须是接受态：

$$
E^n,\phi^n.
$$

`eval_fields` 必须是当前 Newton candidate：

$$
E^{n+1},\phi^{n+1}.
$$

---

## F6.3 pairing 构造失败

如果：

```cpp
pairing_ok == false
```

则 candidate evaluation 失败。

不得 fallback 到：

```cpp
eval_fields.Ex
```

不得 fallback 到：

```cpp
0.5*(E_n+E_np1)
```

---

## F6.4 构造 cell force field

使用：

```cpp
e_pair_cell[ix] =
    0.5 *
    (
        pairing_face[ix]
        +
        pairing_face[ix + 1]
    );
```

然后：

```cpp
evaluate_local_residual(
    ...,
    e_pair_cell,
    ...
)
```

---

## F6.5 禁止再使用 final Ex 推动 J1 bulk

当前：

```cpp
eval_fields.Ex[ng + ix]
```

不再允许作为 J1 的 `u_parallel` force field。

它仍可以：

- 保存在 candidate field；
- 输出诊断；
- 用于比较。

但不允许参与 J1 phase-space residual。

---

# 阶段 F7：新增真正的 x/u/Poisson 三重能量诊断

当前 J0 audit 的最大问题是它只验证了 u flux 自己的求和分部，没有证明 x-current 和 u-work 是同一个离散功对象。

必须增加下面的诊断。

## F7.1 新增结果字段

至少增加：

```text
joint_midpoint_delta_k_x
joint_midpoint_delta_k_u
joint_midpoint_u_face_work
joint_midpoint_force_current_work
joint_midpoint_charge_current_work
joint_midpoint_charge_current_work_interior
joint_midpoint_charge_current_work_endpoint
joint_midpoint_poisson_potential_charge_work
joint_midpoint_poisson_transport_residual
joint_midpoint_current_pair_residual
joint_midpoint_force_charge_residual
joint_midpoint_energy_residual
joint_midpoint_domain_energy_change
```

不得继续让一个字段同时代表两个不同量。

---

## F7.2 charge-current work

使用生产 `accepted_bundle.charge_current_face`。

使用生产 `pairing_face`。

face quadrature 必须与 `FieldParticlePowerAudit` 已验证规则一致：

- 全局物理 endpoint：权重为半个 cell；
- 内部 face：一个完整 `dx`；
- MPI shared face：全局只计一次。

已有 Gate I 就是这样做的，不能发明另一套权重。

定义：

$$
W_J=\Delta t\left\langle E^{\mathrm{pair}},J^{\mathrm{charge}}\right\rangle_f.
$$

---

## F7.3 u-face kinetic work

定义：

$$
W_u=\Delta t\sum_{i,j+1/2,k}\Delta K_{j+1/2,k}F^u_{i,j+1/2,k}.
$$

这个量必须直接使用最终 production `u_flux_rate`。

禁止重新构造 u flux。

---

## F7.4 current-pair residual

定义：

$$
R_{uJ}=W_u-W_J.
$$

这就是最重要的新诊断。

同时必须输出：

$$
W_F=Delta tsum_{i,j,k}E_i^{\mathrm{pair}}q_ev^H_{j,k}M^{n+1/2}_{i,j,k},
$$

以及按当前 face ownership 拆分的：

```text
W_J_interior
W_J_endpoint
R_FJ = W_F - W_J
```

其中 `W_J_interior + W_J_endpoint = W_J`。这些字段只用于定位内部 x/u 配对错误与周期接缝/物理 endpoint 拓扑不一致；不得据此缩放、平均或重构 `J_charge`。

如果 F3、F4、F5 正确，background-only J1 的内部相空间离散应该使这个残差达到 roundoff。

---

## F7.5 Poisson transport residual

现有 Poisson 满足：

$$
\Delta U_E-W_{\mathrm{electrode}}+W_J\approx0.
$$

因此定义：

$$
R_{PJ}=\Delta U_E-W_{\mathrm{electrode}}+W_J.
$$

必须和：

```cpp
OpenElectrostaticSolver::evaluate_work_identity()
```

的结果交叉检查。

如果 `R_PJ` 大，而 Poisson scalar identity 自身正常：

> 检查 continuity current / pairing field 的使用。

**不得修改 Poisson。**

---

## F7.6 总能量残差

定义：

$$
R_E=R_{uJ}+R_{PJ}.
$$

因此：

$$
R_E=W_u+\Delta U_E-W_{\mathrm{electrode}}.
$$

对 background-only、periodic manufactured x、零碰撞、零 Beam、零 Tail 的 J1：

$$
R_E
\rightarrow
O(\epsilon_{\mathrm{machine}}).
$$

---

# 阶段 F8：修复当前 energy 字段被覆盖的问题

当前 `advance_joint_midpoint()` 在能量 gate 前把：

```cpp
joint_midpoint_energy_residual
```

设为联合离散残差。

但是 commit 后又执行：

```cpp
result.joint_midpoint_energy_residual =
    result.ledger.domain_energy_change;
```

这把该字段的含义改掉了。

必须删除这种覆盖。

保持：

```text
joint_midpoint_energy_residual = R_E
```

另存：

```text
joint_midpoint_domain_energy_change
```

用于：

$$
(K^{n+1}+U_E^{n+1})-(K^n+U_E^n).
$$

两个量不得混用。

---

# 阶段 F9：重新定义 J0 验收

当前 J0 PASS 不够。

J0 必须拆成至少以下六项。

## J0-A：cell mass 单位

验证：

$$
M=f\Delta x\Delta u A_\perp.
$$

---

## J0-B：u-flux geometry

验证：

$$
F_u=a\frac12\left(\frac{M_L}{\Delta u_L}+\frac{M_R}{\Delta u_R}\right).
$$

特别防止重复：

```text
dx
ring
```

---

## J0-C：midpoint 一致性

验证：

```text
global evaluate_residual
```

和

```text
single-rank evaluate_local_residual
```

完全使用同一：

$$
M^{n+1/2}.
$$

---

## J0-D：u-work / force-current adjoint

验证：

$$
W_u=W_{\mathrm{force}}.
$$

---

## J0-E：x-current / force-current adjoint

内部面验证：

$$
J^{\mathrm{charge}}_{i+1/2}=\frac{J_i^{\mathrm{force}}+J_{i+1}^{\mathrm{force}}}{2}.
$$

---

## J0-F：Poisson pairing

使用现有：

```cpp
build_potential_pairing_field()
```

验证：

$$
\Delta U_E-W_{\mathrm{electrode}}+W_J=O(\epsilon).
$$

最后才验证：

$$
R_E=R_{uJ}+R_{PJ}=O(\epsilon).
$$

### 只有六项全部 PASS 才允许称 J0 PASS。

---

# 阶段 F10：重新执行 J1

首先只执行单 rank。

推荐构建：

```bash
cmake .. \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON
```

然后：

```bash
make -j4 \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test \
  vpfp_poisson_work_identity_test
```

执行：

```bash
ctest -R \
"joint_phase_space_midpoint_unit_test|joint_phase_space_midpoint_energy_test|vpfp_poisson_work_identity_test" \
--output-on-failure
```

---

## F10.1 每个 J1 结果必须输出

至少输出：

```text
case
accepted
finite
gauss_ok
converged
iterations

residual_linf
poisson_residual_linf

delta_k_x
delta_k_u
u_face_work
charge_current_work

field_energy_change
electrode_work
poisson_potential_charge_work

current_pair_residual
poisson_transport_residual
energy_residual

domain_energy_change

min_mass
max_mass

failure_code
failure_stage
```

---

# 11. 根据残差自动决定“停止在哪里”

执行模型不得自行判断。

严格执行下面的决策树。

## 情形 A

如果：

$$
|R_{PJ}|
$$

大，而：

```text
OpenPoissonWorkIdentity.residual
```

小：

立即停止。

报告：

```text
Poisson scalar identity PASS
Poisson-current pairing FAIL
```

只检查：

- `charge_current_face`;
- x continuity；
- pairing field 使用；
- MPI face ownership。

不得修改 Poisson。

---

## 情形 B

如果：

$$
|R_{uJ}|
$$

大：

立即停止。

只检查：

- u trace；
- `vH`;
- x flux；
- midpoint mass；
- geometry factors。

不得修改 Newton。

不得修改 Poisson。

---

## 情形 C

如果：

$$
|R_{uJ}|
$$

小，

$$
|R_{PJ}|
$$

小，

但：

$$
|R_E|
$$

大：

这是符号、重复计数或诊断 time layer 错误。

检查：

```text
sign
dt
face weight
electrode work sign
double counting
```

不得改变物理算子。

---

## 情形 D

如果所有 energy identity 都 PASS，但 Newton 出现：

```text
code 73
code 74
code 76
```

这时才允许进入 F11。

---

# 阶段 F11：仅在能量代数已经 PASS 后修 Newton positivity

当前 line search 会允许带显著负质量的 Newton trial 被接受，只要 residual 下降。源码确实对 Jacobian probe 和实际 line-search candidate 都传入 `allow_negative_probe=true`，最终才用 code 76 检查。

这不是当前 code 75 的第一修复对象。

只有 F10 情形 D 才执行本阶段。

---

## F11.1 Jacobian probe

有限差分 Jacobian probe 继续允许 signed state。

即：

```text
Jacobian probe:
allow_negative_probe = true
```

不要修改。

因为这是代数导数采样，不是可接受物理状态。

---

## F11.2 实际 line-search trial

真实 line-search candidate 必须做 positivity-domain gate。

在调用完整 `evaluate(trial,...)` 之前检查。

定义全局质量尺度：

$$
M_{\mathrm{scale}}=\max |M|.
$$

允许的 roundoff 负值：

$$
M_{\mathrm{tol}}=4096\epsilon_{\mathrm{mach}}\max\left(1,M_{\mathrm{scale}}\right).
$$

如果某个 trial 有：

$$
M_i<-M_{\mathrm{tol}},
$$

则：

```text
该 lambda 直接拒绝
lambda *= 0.5
```

不要调用完整 candidate Poisson。

---

## F11.3 严禁 clip

不得写：

```cpp
trial[i] = std::max(0.0, trial[i]);
```

不得把负质量改成零。

因为这会破坏已经证明的 phase-space residual。

---

## F11.4 暂时不要改 GMRES 算法

当前 code 75 已经发生在 nonlinear residual 收敛之后，因此 GMRES 不是该能量失败的第一根因。

只有满足：

```text
R_uJ PASS
R_PJ PASS
R_E PASS
```

同时仍然发生 line-search failure，才允许单独评估 GMRES。

在此之前：

- 不增加 iteration；
- 不增加 restart dimension；
- 不放松 tolerance；
- 不换 nonlinear solver。

---

# 12. J1 最终阶段门

J1 必须同时通过：

## 12.1 平衡测试

```text
smooth-background
```

要求：

- finite；
- Gauss PASS；
- nonlinear PASS；
- 无显著负质量；
- energy PASS。

---

## 12.2 非平凡场测试

```text
smooth-perturbed-background
```

要求：

- 初始 `E != 0`；
- 有真实 field-particle energy exchange；
- `R_uJ` 达到浮点误差量级；
- `R_PJ` 达到 Poisson 已验证误差量级；
- 总能量 gate 至少满足现有生产要求：

$$
\frac{|R_E|}
{\max(
|W_u|,
|\Delta U_E|,
|W_{\mathrm{electrode}}|,
\epsilon_{\mathrm{floor}}
)}
\le10^{-8}.
$$

但单元测试应尽量达到明显优于此阈值的水平。

不得通过刚好把结果压在 `1e-8` 以下就停止分析。

---

# 13. J1 MPI 门

单 rank PASS 后才能执行。

必须执行：

```text
1 rank
2 ranks
5 ranks
```

比较：

```text
R_uJ
R_PJ
R_E
delta_K_u
W_J
field_energy_change
electrode_work
```

MPI 改变只能产生舍入顺序差。

不得随 rank 数产生系统性改变。

---

## 13.1 MPI face ownership

face 功积分必须严格复用当前 Gate I 的所有权原则：

- global physical endpoint：半权重；
- internal face：完整权重；
- shared face：只由一个 rank 计入。

现有 `FieldParticlePowerAudit` 已经实现了这一规则。

不得发明第二套 MPI face 规则。

---

# 14. J2：接入真正开放 background 边界

**只有 J1 全部 PASS 才能执行。**

这一步是保证最终模型没有被 J1 的 manufactured periodic test 偷换。

生产最终仍必须回到：

```text
开放 background
非周期 Poisson
```

---

## J2.1 不再周期交换物理 endpoint

内部 MPI face：

- 继续正常交换 midpoint state。

全局左、右物理边界：

- 使用现有 `OpenBackgroundBoundary`；
- 不允许 wrap 到另一端。

---

## J2.2 reservoir

左/右 reservoir 必须继续遵循主物理方案：

- incoming characteristic：reservoir；
- outgoing characteristic：域内状态；
- 不允许每步重置边界 cell；
- 不允许体源化 reservoir。

---

## J2.3 boundary kinetic flux 必须来自同一个 x flux

对于每一个物理 x-face，使用\*\*实际 production `F_x`\*\*计算能量流：

$$
F_K^x=\sum_{j,k}K_{j,k}F^x_{j,k}.
$$

不得由边界 cell 前后能量差反推。

这样：

- number boundary ledger；
- charge boundary current；
- kinetic boundary ledger；

都来自同一个 x flux。

---

## J2.4 J2 能量验收

完整开放 background 能量账必须包括：

- bulk kinetic change；
- field energy change；
- electrode work；
- reservoir inflow kinetic energy；
- background outflow kinetic energy。

这里不再要求“每个物理 endpoint 的 `J_charge` 必须等于内部 `G*J_force`”。

物理边界可能产生独立 boundary contribution。

正确要求是：

1. 核心区域的 x/u current pairing 闭合；
2. endpoint 差异被真实 boundary ledger 解释；
3. 完整开放系统能量账闭合。

不得为了消灭 endpoint residual 人工修改边界 current。

---

# 15. J3a：重新接入 Beam

只有 J2 PASS 后执行。

本阶段：

```text
Beam ON
Tail OFF
conversion OFF
collision OFF
```

不得修改 Beam 现有开放注入/流出物理。

Beam 已有 charge-conserving trajectory current 和 kick work 诊断，旧 Gate I 也已经把 PIC force work/current pairing 分开审计。

本阶段只做：

> 把 Beam 接入已经通过 J1/J2 的同一场—粒子 energy ledger。

检查：

- Beam injection energy；
- Beam outflow energy；
- Beam kick work；
- Beam trajectory current；
- field energy；
- electrode work。

失败时只处理 Beam coupling。

不得修改 bulk joint operator。

---

# 16. J3b：接入 background Tail PIC

只有 J3a PASS 后执行。

本阶段：

```text
Beam ON
Tail ON
conversion 暂时 OFF
collision OFF
```

Tail：

- 必须进入 Poisson density；
- 使用已有 PIC trajectory current；
- 保持开放流出；
- 流出前完成域内轨迹；
- Tail kick work 进入 energy ledger。

如果失败：

只检查 Tail gather/deposit/work pairing。

不得修改 J1 bulk。

---

# 17. J3c：接入 bulk-to-tail conversion

只有 J3b PASS。

转换是内部表示变化。

必须保持：

$$
\Delta N_{\mathrm{bulk}}^{\mathrm{conv}}
+
\Delta N_{\mathrm{tail}}^{\mathrm{conv}}
=0,
$$

$$
\Delta P_{\mathrm{bulk}}^{\mathrm{conv}}
+
\Delta P_{\mathrm{tail}}^{\mathrm{conv}}
=0,
$$

$$
\Delta K_{\mathrm{bulk}}^{\mathrm{conv}}
+
\Delta K_{\mathrm{tail}}^{\mathrm{conv}}
=0.
$$

主物理规格已经明确规定转换不能成为外部能量源。

因此：

```text
conversion_energy_source_external = 0
```

任何 conversion energy residual 都只能作为表示转换错误处理。

不得把它加到右端“补偿”。

---

# 18. J3d：接入 H10 Tail-to-Bulk return

只有 J3c PASS。

H10 仍然只是内部表示转换。

已有历史 A/B 已经证明：

- return ON/OFF 两档都存在近似同量级的公共能量余额；
- 两者只差约 0.6%；
- 因此公共能量缺陷不能归因于 H10；
- 禁止把公共 energy residual 塞进 H10 projection compensation。

因此不得修改已经通过的：

```text
H10 N/P/K transaction
local representation optimization
```

除非新的隔离测试明确证明它们回归。

---

# 19. J3e：最后接入碰撞

碰撞必须最后接入。

原因：

碰撞不是当前已经确认的无碰撞场—粒子结构根因。

接入顺序：

```text
collision half step
joint collisionless field-particle step
collision half step
```

碰撞自身：

- bulk-bulk；
- tail-tail；
- tail-bulk；
- bulk reaction；

继续使用当前主物理方案的 pair registry 和 reaction ledger。

不得为了总能量好看把碰撞交换写成全局补丁。

---

# 20. 最终生产能量账

完整模型最终必须验证：

$$
\Delta\left(K_{\mathrm{bulk}}+K_{\mathrm{tail}}+K_{\mathrm{beam}}+U_E\right)
$$

等于所有真实外部能量交换之和。

外部项包括：

- background reservoir inflow；
- background outflow；
- Tail outflow；
- Beam injection；
- Beam outflow；
- electrode work；
- 碰撞模型中明确存在的外部 reservoir exchange。

而：

```text
bulk -> tail
tail -> bulk
```

绝对不是外部能量源。

---

# 21. 明确禁止的“修复方法”

执行模型任何阶段都不得实施以下措施。

## 21.1 禁止能量投影

不得计算 residual 后执行：

```text
调整 E
调整 f
调整 current
```

使总能量强行相等。

---

## 21.2 禁止电流缩放

不得：

```cpp
J *= correction_factor;
```

---

## 21.3 禁止电场缩放

不得：

```cpp
E *= correction_factor;
```

---

## 21.4 禁止修改 Poisson 迎合粒子

不得因为：

```text
R_PJ != 0
```

就修改 `OpenElectrostaticSolver`。

先证明：

- charge current；
- continuity；
- pairing field；

使用正确。

---

## 21.5 禁止裁剪 distribution 修能量

不得：

```cpp
M = std::max(M, 0.0);
```

然后继续当作守恒解。

---

## 21.6 禁止重新引入 FCT 修 J1 能量

J1 的目标首先是证明离散 Hamiltonian/Poisson 能量结构。

FCT 会改变最终 flux，因此在基础能量恒等式证明之前不得重新加入。

---

## 21.7 禁止依靠减小 dt 宣称修复

减小 dt 可以降低误差，但不能替代：

$$
R_{uJ}=0
$$

和

$$
R_{PJ}=0
$$

的离散结构。

---

# 22. 每个阶段必须采用的工作模板

交给较弱模型时要求它在每一个阶段严格输出：

````markdown
## 当前阶段

阶段编号：

## 本阶段允许修改文件

- ...

## 本阶段实际修改文件

- ...

## 修改内容

1. ...
2. ...

## 编译命令

```bash
...
````

## 测试命令

```bash
...
```

## 数值结果

```text
...
```

## 阶段门

PASS / FAIL

## 若 FAIL

停止。
不得进入下一阶段。
列出实际失败数值和对应源码位置。

````

不得一次完成 F1-F11。

---

# 23. 对执行模型的强制工作顺序

最终完整顺序固定为：

```text
F0  冻结基线
 |
F1  修 J1 初始化测试
 |
F2  加非平凡场测试
 |
F3  修 u-flux 重复几何因子
 |
F4  修真正的 M midpoint
 |
F5  构造 u-trace 的精确转置 vH
 |
F6  使用 Poisson discrete-gradient pairing field
 |
F7  加 x/u/Poisson 三重功诊断
 |
F8  修 energy diagnostic 字段语义
 |
F9  强化 J0
 |
F10 重跑 J1
 |
 +---- energy identity FAIL ----> 只修对应离散身份
 |
 +---- energy PASS + Newton FAIL -> F11 positivity
 |
J1 1/2/5 rank PASS
 |
J2 开放 background
 |
J3a Beam
 |
J3b Tail
 |
J3c bulk->tail
 |
J3d tail->bulk
 |
J3e collisions
 |
短完整生产
 |
正式生产
````

**禁止跨阶段。**

---

# 24. 最关键的三个代数阶段门

较弱模型执行时，不需要“理解整个 Hamiltonian 理论”。

只要机械验证下面三个式子。

## 第一关：u flux 和 force current

必须有

$$
R_{uF}=\Delta t\sum \Delta K F_u-\Delta t\Delta x\sum_i E_i^{\mathrm{pair}}J_i^{\mathrm{force}}\approx0.
$$

失败只检查：

```text
u trace
u geometry
vH
```

---

## 第二关：force current 和 charge current

必须有

$$
R_{FJ}=\Delta t\Delta x\sum_i E_i^{\mathrm{pair}}J_i^{\mathrm{force}}-\Delta t\left\langle E^{\mathrm{pair}},J^{\mathrm{charge}}\right\rangle_f\approx0.
$$

失败只检查：

```text
x flux
face averaging
face ownership
G/G*
```

---

## 第三关：charge current 和 Poisson

必须有

$$
R_{PJ}=\Delta U_E-W_{\mathrm{electrode}}+\Delta t\left\langle E^{\mathrm{pair}},J^{\mathrm{charge}}\right\rangle_f\approx0.
$$

失败时：

> **禁止修改 Poisson。**

先检查是否真正使用了 Poisson 已提供的 pairing field。

---

## 三者合并

当三个离散关系都成立时：

$$
R_E=\Delta K_u+\Delta U_E-W_{\mathrm{electrode}}\approx0.
$$

这才是 J1 应当证明的东西。

---

# 25. 为什么这套方案比当前“继续调 J1 Newton”更优先

当前 J1 已经出现过：

```text
nonlinear residual converged
Gauss converged
energy code 75
```

所以单纯：

- 增加 Newton 次数；
- 增大 GMRES 维数；
- 放宽 residual；
- 放宽 energy gate；

都没有触及场—粒子离散身份。

更重要的是，现在源码又能确认：

```text
错误 J1 初态
+
错误 u-flux 几何权重
+
生产 local residual 没用 midpoint
+
force field 没用 Poisson pairing field
+
vH 不是 u-trace 精确转置
```

所以**现在不应该再把主要精力投入 nonlinear solver**。

正确顺序一定是：

> 先把离散代数闭合，再处理 nonlinear admissibility。

---

# 26. 对现有文档的修订要求

代码修复完成后才修改文档。

## 26.1 更新

```text
docs/J1联合中点当前失败分析_2026-08-21.md
```

必须保留历史数据，但明确标记：

```text
1.1132815186046597e14 J/m^2
```

来自旧 J1 fixture，其中离子 density 在 `compute_moments()` 前读取，因此该数值不能作为纯联合离散缺陷的定量证据。

不得直接删除历史记录。

---

## 26.2 更新

```text
docs/VPFP能量问题根因链_已确认结论.md
```

旧 Strang 根因继续保留。

新增当前 J1 实现层根因链：

```text
J1 fixture initialization
    ↓
cell-integrated mass geometry
    ↓
true temporal midpoint
    ↓
u-trace transpose velocity
    ↓
x charge-current / u-force-current adjoint
    ↓
Poisson discrete-gradient pairing field
    ↓
joint energy closure
```

---

## 26.3 不修改主物理规格

除非用户另行明确要求，否则本轮不要修改：

```text
docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
```

它继续作为最高物理模型规格。

---

# 27. 如果通过 SSH 在服务器执行

仓库自己的 SSH 协议明确要求：

- 默认工作目录为
  `/HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver`
- 集群上的修改操作必须先取得批准；
- 必须先观察再修改；
- 必须遵守调度器；
- 不得越出指定项目路径。

因此自动模型若需要远程操作：

```text
没有获得修改批准：
只允许读取、git status、编译已有内容等不改变仓库的操作。

获得批准以后：
才允许修改代码。
```

不得把本文当作自动绕过该审批规则的授权。

---

# 28. 最终结论

目前最合理的根因链不是简单的：

```text
Strang split
→ J1 center trace
→ 仍不保能
→ 需要发明更复杂 Hamiltonian 格式
```

而应修正为：

```text
旧生产 Strang/PPM
存在已确认的 x-current / u-work 不同源
        ↓
联合相空间方向正确
        ↓
当前 J1 测试初始化存在错误
        ↓
当前 local J1 residual 没真正使用 M midpoint
        ↓
当前 u flux 对 cell-integrated M 重复乘 dx × ring
        ↓
当前 u force 使用 final candidate Ex，
而不是 Poisson-compatible pairing field
        ↓
当前 vH 在非均匀 u 网格上不是 u trace 的精确转置
        ↓
当前 J0 audit 又没有直接测试
u-work ↔ x-charge-current 的交叉恒等式
        ↓
所以现有 J1 code75 结果尚不能用于判断
“联合中点思想本身失败”
        ↓
先完成 F1-F9
        ↓
重新得到干净 J1 数据
        ↓
只有此时仍有 R_uJ，才讨论更高阶离散括号问题
```

这里最重要的一点是：**不需要改变《开放Beam\_非周期场\_VPFP\_PIC完整重构方案.md》定义的任何物理模型。** 修复对象应严格限制在“同一个物理 VPFP-PIC 模型的离散场—粒子功共轭关系”。主方案要求的开放 Beam、非周期 Poisson、Eulerian bulk + PIC tail、内部 bulk/tail 转换、开放边界和碰撞结构全部保持不变。

而且从当前源码看，**下一步最先应该做的不是改 Newton，也不是继续尝试新的 Hamiltonian flux，而是严格按 F1 → F3 → F4 → F5 → F6 的顺序，把已经能从源码确定的实现错误全部清掉。** 在这些问题清掉之前，继续根据现在的 `1.113e14 J/m²` 推导新的数值格式，很容易在错误测试和错误量纲之上继续叠加复杂度。

---

# 29. 执行导航与命令附录（本次整理新增）

本节只补充执行顺序、命令和输出位置；不替代原有阶段的实现细节、数学定义或验收要求。所有命令默认在集群项目根目录执行：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
```

本附录遵循原文 §27：远程修改代码前仍须遵守 `docs/ssh远程连接协议.md` 的授权与先观察后修改规则。

## 29.1 统一运行前置条件

### 29.1.1 模块与构建目录

下面命令不删除 `build/`，避免破坏已有构建产物或无关测试。首次配置、CMake 变更后、或新增加测试 target 后，执行：

```bash
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON
```

每次只构建当前阶段需要的 target。除非 CMake 配置已经失效，不得使用会删除 `build/` 的重建脚本作为本方案的默认流程。

### 29.1.2 唯一输出目录

每次测试先创建唯一目录，不覆盖旧结果：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_${RUN_ID}"
mkdir -p "$OUT"
```

若调度环境不允许 `date`，手工指定不重复的 `RUN_ID`。所有 `.result`、日志与 MPI 子目录必须写入该 `$OUT`。

### 29.1.3 结果读取规则

每次命令结束后至少读取：

```bash
find "$OUT" -maxdepth 2 -type f -name '*.result' -print -exec sed -n '1,220p' {} \;
```

非零退出、`status=FAIL`、缺失预期 `.result`、或产生 NaN/Inf 都是当前阶段 FAIL。不得仅因程序退出码为零就判定 PASS。

## 29.2 F0：基线冻结命令

对应原文 F0。此阶段只记录，不修改源文件：

```bash
git status --short
git rev-parse HEAD
git diff --check
```

将三段输出原样写入当前阶段报告中的 `baseline_commit`、`working_tree_clean` 和 `diff_check`。工作树非空不等于失败，但必须停止并记录已有修改；不得 reset、checkout 或覆盖。

## 29.3 F1：fixture 初始化修复后的最小回归

对应原文 F1。此命令仅在完成 F1.1/F1.2 的源码修改后执行。当前 target 名称为 `joint_phase_space_midpoint_energy_test`。

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result "$OUT/f1_smooth_background.result" \
  > "$OUT/f1_smooth_background.out" \
  2> "$OUT/f1_smooth_background.err"
```

验收与 F1.3 完全一致：记录 `accepted`、`converged`、`residual_linf`、`poisson_residual_linf`、`energy_residual`、`failure_code`，并额外确认 F1.2 新增的初态字段已经实际输出。此步骤不能替代 F2 的非平凡场测试。

## 29.4 F2：非平凡场测试命令

对应原文 F2。**在 F2 源码尚未增加 `smooth-perturbed-background` CLI case 前，本命令必然返回 usage/exit 2，不得把这种结果记为数值 FAIL。** 完成 F2.1/F2.2 后：

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/f2_smooth_perturbed_background.result" \
  > "$OUT/f2_smooth_perturbed_background.out" \
  2> "$OUT/f2_smooth_perturbed_background.err"
```

除原文 F2 与 §12.2 指标外，结果必须显式出现非零初始场或等价的 `initial_E_linf`，否则该测试仍然只是在重复平衡态。

## 29.5 F3：u 通量量纲修复与单项测试命令

对应原文 F3。F3.4 明确要求新增 geometry test，但 F3.5 明确规定阶段门为既有 `joint_phase_space_midpoint_unit_test`。因此正确实现方式是：**只在** `tests/joint_phase_space_midpoint_unit_test.cpp` **内增加该 geometry case，并继续由既有 target 编译、运行和报告。** 不新增 executable，不新增 CMake target，不另写一套 flux。

新增 case 必须直接调用生产 `JointPhaseSpaceMidpointOperator::build_periodic_center_flux()`，以两个不同 `dx` 和两个不同 `uperp_ring_area` 检查 F3.4 的解析式；该 case 必须被 `--case all` 覆盖并在 `.result` 中有独立 PASS 字段。

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f3_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/f3_joint_midpoint_unit.result" \
  > "$OUT/f3_joint_midpoint_unit.out" \
  2> "$OUT/f3_joint_midpoint_unit.err"
```

F3 PASS 的唯一命令门是 `f3_joint_midpoint_unit.result` 中 `status=PASS`，且新增 geometry case 的相对误差达到浮点舍入量级。F3 PASS 后才能开始 F4。

## 29.6 F4：真实 midpoint mass 一致性测试命令

对应原文 F4。F4.6 的 local/global 一致性测试同样增加到既有 `tests/joint_phase_space_midpoint_unit_test.cpp`，由 `joint_phase_space_midpoint_unit_test --case all` 执行；不得新增独立 target。

新增 case 对同一 `M_old`、`M_candidate`、`E`、`dt` 直接调用 production `evaluate_residual()` 和 `evaluate_local_residual()`，并在 result 中独立输出 `x_flux_rate`、`u_flux_rate`、`mass_delta_x`、`mass_delta_u`、`mass_delta_total`、`residual` 的最大差值。不得复写替代 residual。

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f4_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/f4_joint_midpoint_unit.result" \
  > "$OUT/f4_joint_midpoint_unit.out" \
  2> "$OUT/f4_joint_midpoint_unit.err"
```

F4 PASS 的唯一命令门是 `f4_joint_midpoint_unit.result` 中 `status=PASS`，且 F4.6 的六类量全部在浮点 roundoff 内一致。失败时只允许修改原文 F4.1--F4.5 指定的 local residual 路径。

## 29.7 F5：u-work 严格转置速度测试命令

对应原文 F5。F5.3 的严格 u-energy adjoint 测试也增加到 `tests/joint_phase_space_midpoint_unit_test.cpp`，由既有 `joint_phase_space_midpoint_unit_test` 执行；不得新增 executable。

新增 case 必须用 production 非均匀 `CylindricalVelocityGrid`，构造任意正且依赖 `ix/j/k` 的 `M_mid` 和非零 `E_cell`。左式直接由生产 `u_flux_rate` 求和，右式直接由生产 `build_hamiltonian_velocity()` 的 `vH` 与 `M_mid` 求和；不得调用同一个 helper 计算两边。

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f5_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/f5_joint_midpoint_unit.result" \
  > "$OUT/f5_joint_midpoint_unit.out" \
  2> "$OUT/f5_joint_midpoint_unit.err"
```

F5 PASS 的唯一命令门是 `f5_joint_midpoint_unit.result` 中 `status=PASS`；其中新增 adjoint case 必须报告 `W_u-W_J`、`vH` finite/antisymmetry/speed 检查，并满足原文 F5.3--F5.5。F5 PASS 后才能开始 F6。

## 29.8 F6：Poisson pairing field 接入测试命令

对应原文 F6。F6 没有要求新建独立 Poisson 测试，也没有要求运行 `vpfp_poisson_work_identity_test`。本阶段只在既有 `joint_phase_space_midpoint_unit_test` 中增加/更新一个 J1 pairing-field 使用 case：直接检查 J1 传入 `evaluate_local_residual()` 的 cell force field 来自 `build_potential_pairing_field()` 的相邻 face 平均，而不是 final `Ex`。

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f6_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/f6_joint_midpoint_unit.result" \
  > "$OUT/f6_joint_midpoint_unit.out" \
  2> "$OUT/f6_joint_midpoint_unit.err"

```

F6 PASS 要求 `f6_joint_midpoint_unit.result` 为 `status=PASS`，并有独立 pairing-field case 证明 J1 复用 pairing field。F6 失败时禁止修改 `OpenElectrostaticSolver`；是否验证 Poisson scalar identity 属于后续 F7/F9 的诊断与 J0-F 汇总，不是 F6 的额外阶段门。

## 29.9 F7：三重功诊断字段验收命令

对应原文 F7。F7 的产物是写入 `VpfpStepResult` 和 J1 `.result` 的诊断字段，不要求也不授权新增独立 executable。修改字段并扩展现有 `joint_phase_space_midpoint_energy_test.cpp` 的输出后，直接运行既有 J1 非平凡场 case：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f7_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_energy_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/f7_smooth_perturbed_background.result" \
  > "$OUT/f7_smooth_perturbed_background.out" \
  2> "$OUT/f7_smooth_perturbed_background.err"
```

`code 75` 若仍发生，命令可非零退出，但 `.result` 仍必须写出 F7.1 的全部字段以及 `W_u`、`W_J`、`R_uJ`、`R_PJ`、`R_E`。F7 的验收重点是字段完整、同一 accepted/candidate 时间层可追溯，并能用结果复算 F7.4--F7.6；它不是在 F3--F6 未完成前要求能量 gate 必须 PASS 的阶段。

## 29.10 F8：能量字段语义回归命令

对应原文 F8。F8 不新增 target；扩展 `joint_phase_space_midpoint_energy_test.cpp` 的结果输出后，使用既有可接受的 `smooth-background` case 验证 accepted-step commit 之后字段未被覆盖：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f8_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_energy_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result "$OUT/f8_smooth_background.result" \
  > "$OUT/f8_smooth_background.out" \
  2> "$OUT/f8_smooth_background.err"
```

F8 PASS 要求该 result 同时输出 `joint_midpoint_energy_residual` 和 `joint_midpoint_domain_energy_change` 的独立字段，并能确认前者保持 F7.6 的 `R_E` 含义、后者只表示 domain energy change。不得因平衡态两者数值都为零，就省略字段或把同一内存字段复用。

## 29.11 F9：J0 六项汇总验收命令

对应原文 F9。J0-A 至 J0-F 都应作为 `joint_phase_space_midpoint_unit_test` 的独立 case 汇总；J0-F 直接检查该测试中实际调用的 `build_potential_pairing_field()`。F9 不要求也不授权将 `vpfp_poisson_work_identity_test` 作为额外阶段门，因此只运行既有 J0 target：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f9_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/f9_j0_all.result" \
  > "$OUT/f9_j0_all.out" \
  2> "$OUT/f9_j0_all.err"

```

F9 的 PASS 含义严格受原文 J0-A 至 J0-F 约束。只有 F3、F4、F5、F6、F7、F8 的各自阶段工作已完成，且 `f9_j0_all.result` 为 PASS、J0-A--J0-F 均有独立输出，才允许进入 F10；可选 `ctest` 汇总不能替代该显式结果文件。

## 29.12 F10：单 rank J1 回归命令

对应原文 F10、§11 和 §12。只在 F3--F9/J0 已通过后运行两个 J1 case：

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4

for CASE in smooth-background smooth-perturbed-background; do
  yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_phase_space_midpoint_energy_test \
    --case "$CASE" \
    --result "$OUT/j1_${CASE}_n1.result" \
    > "$OUT/j1_${CASE}_n1.out" \
    2> "$OUT/j1_${CASE}_n1.err" || exit 1
done
```

对每个 case 按原文 F10.1 完整输出字段验收，并按 §11 的 A--D 决策树处理。`R_uJ` 或 `R_PJ` 失败时只允许修改该决策树列出的对应离散对象；禁止借由 F11 放宽或隐藏 code 75。

## 29.13 F11：仅 energy PASS 后的 positivity 回归命令

对应原文 F11。F11 不新增独立物理测试；它必须重跑 29.6 的两个 J1 case，并额外检查每个迭代记录：

```text
line_search_alpha
trial_min_mass
accepted
failure_code
```

运行命令与 29.12 完全相同，但输出目录必须新建，例如：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/joint_energy_f11_${RUN_ID}"
mkdir -p "$OUT"
# 随后逐字执行 29.6 的命令。
```

只有原文 §11 情形 D 才允许使用这一节。不得把 F11 当作 F3--F10 能量代数失败的替代方案。

## 29.14 J1 MPI 门命令

对应原文 §13。前置条件是 29.6 的两个单 rank case PASS。

**实现前置条件：** 当前 J1 测试源码若仍将 `nx_global` 固定为 `4`，则不能以 `5 ranks` 正常分解。执行 5-rank 命令前，必须按 F2/F10 的测试实现将该 manufactured 网格显式扩展至 `nx_global >= 5`，且不改变测试的离散公式；否则 5-rank 的失败只能说明网格分解非法，不能说明 J1 MPI 配对失败。

完成该前置条件后，对每个 case 和 rank 数运行：

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4

for CASE in smooth-background smooth-perturbed-background; do
  for NP in 1 2 5; do
    yhrun -n "$NP" --cpu-bind=cores stdbuf -oL -eL \
      ./build/joint_phase_space_midpoint_energy_test \
      --case "$CASE" \
      --result "$OUT/j1_${CASE}_n${NP}.result" \
      > "$OUT/j1_${CASE}_n${NP}.out" \
      2> "$OUT/j1_${CASE}_n${NP}.err" || exit 1
  done
done
```

按原文 §13 比较 `R_uJ`、`R_PJ`、`R_E`、`delta_K_u`、`W_J`、`field_energy_change` 和 `electrode_work`。只允许可解释的 MPI 归约舍入差；系统性 rank 依赖为 FAIL。

## 29.15 J2/J3 的测试命令契约

原文 J2/J3a--J3e 定义了严格的接入顺序，但当前文档没有给出相应测试 target 名称。为避免未来执行者擅自复用 J1 periodic manufactured test，接入每一阶段时必须**先在 CMake 中新增一个直接调用生产 `VpfpIntegrator` 的回归 target**，再运行下列命令模板。

新增 target 的命名和最小覆盖必须为：

| 原有阶段 | 必须新增的 target | 必须禁用/启用的物理项 |
|---|---|---|
| J2 | `joint_open_background_energy_test` | Beam OFF，Tail OFF，conversion OFF，collision OFF；open/reservoir background + nonperiodic Poisson |
| J3a | `joint_beam_energy_test` | Beam ON，Tail OFF，conversion OFF，collision OFF |
| J3b | `joint_tail_energy_test` | Beam ON，Tail ON，conversion OFF，collision OFF |
| J3c | `joint_bulk_tail_conversion_energy_test` | conversion ON，return OFF，collision OFF |
| J3d | `joint_tail_return_energy_test` | return ON，collision OFF |
| J3e | `joint_collision_energy_test` | collisions ON，Strang collision half-step + joint collisionless step |

每个 target 的 result 至少包含原文对应阶段要求的所有能量 ledger、`accepted`、`finite`、`gauss_ok`、`failure_code`、MPI rank 数和本阶段明确要求的 source/transaction ledger。测试必须用新目录，命令结构固定为：

```bash
# 以 J2 为例；后续 J3a--J3e 仅替换 TARGET、RESULT 与 NP。
TARGET=joint_open_background_energy_test
RESULT="$OUT/j2_open_background_n1.result"

cmake --build build --target "$TARGET" -j4
yhrun -n 1 --cpu-bind=cores stdbuf -oL -eL \
  "./build/$TARGET" --case all --result "$RESULT" \
  > "$OUT/j2_open_background_n1.out" \
  2> "$OUT/j2_open_background_n1.err"
```

每个 J2/J3 target 的单 rank PASS 后，再执行 `-n 2` 与 `-n 5`，并按 J1 MPI 门同样比较数目、动量、动能和电极功 ledger。没有对应 target、没有真实 production operator 调用、或只复写公式的测试都不能作为该阶段 PASS 证据。

## 29.16 最终短生产与正式生产的命令约束

原文 §20 的完整模型只有在 J3e PASS 后才能运行。此时生产命令必须继续使用主物理规格的开放 Beam、非周期 `DIRICHLET_PHI`、reservoir/open background 和完整能量 ledger。不得借由 `joint-midpoint-energy` 的隔离 periodic J1 参数替换生产边界。

正式生产前先由执行者为“短完整生产”增加独立的生产配置文件或 CLI 参数集合，并将其原样保存至 `$OUT/production_command.txt`。该配置必须明确写出：

```text
field boundary
phi_left / phi_right
background boundary
Beam enabled and injection interval
Tail/conversion/return modes
collision model and integrator
dt
MPI ranks and OMP_NUM_THREADS
```

没有完成 J2--J3e 的逐项通过前，禁止给出或执行 120 fs 正式生产命令。

## 29.17 每阶段最小报告索引

为方便按原文 §22 输出，每阶段报告应引用下表，而不是重新组织或改写原有技术结论：

| 当前执行内容 | 必须引用的原文段落 | 最少结果文件 |
|---|---|---|
| F0 | F0.1--F0.3 | `git_status.txt`、`git_head.txt`、`git_diff_check.txt` |
| F1 | F1.1--F1.3 | `f1_smooth_background.result` |
| F2 | F2.1--F2.2 | `f2_smooth_perturbed_background.result` |
| F3--F9 | F3--F9、§24 | `j0_after_f3_to_f9.result`、`poisson_work_identity.result`、新增直接生产算子测试结果 |
| F10 | F10、§11、§12 | 两个 `j1_*_n1.result` |
| F11 | F11 | 新目录下两个 `j1_*_n1.result` 与迭代日志 |
| J1 MPI | §13 | 两个 case 的 `n1/n2/n5` 结果 |
| J2--J3e | §14--§20 | 与 29.9 target 一一对应的 n1/n2/n5 结果 |

此索引只降低阅读和执行成本，不改变原文 §22 的报告模板、§23 的强制顺序以及所有禁止项。

## 30. F10 实测根因更新（2026-08-22）

本节只记录 F10 的已获得证据，不改变原有 F0--F11 的算法约束或第 11 节决策树。

### 30.1 结果

单 rank `smooth-background` 已 PASS。单 rank `smooth-perturbed-background` 则：

```text
accepted=0
finite=1
gauss_ok=1
converged=1
failure_code=75
failure_stage=joint_midpoint_energy_residual
```

对应功分解：

$$
W_u=28.891863089349869,
$$

$$
W_F=28.891863089349886,
$$

$$
W_J=22.224220061661196,
$$

$$
R_{uJ}=6.6676430276886727\ \mathrm{J/m^2},
\qquad
R_{PJ}=-5.2950189122213942\times10^{-6}\ \mathrm{J/m^2}.
$$

其中 $W_F$ 是使用同一 accepted midpoint mass、同一 pairing cell field 和生产 `vH` 的独立 cell-force 功。故：

$$
W_u-W_F\simeq-1.7\times10^{-14}\ \mathrm{J/m^2},
$$

而：

$$
W_F-W_J=6.6676430276886904\ \mathrm{J/m^2}.
$$

### 30.2 强制分流结论

该结果严格落入第 11 节**情形 B**，不是情形 A、C 或 D：

- `R_uJ` 显著；
- `R_PJ` 相比主残差小六个数量级；
- u-face 功至 cell-force 功已经闭合；
- 能量 identity 未通过，故 F11 的 positivity 修复尚不允许进入。

因此下一次代码检查只能审计 `x_flux_rate`、`charge_current_face`、face averaging、midpoint trace、首末 cell/periodic seam 映射和 MPI face ownership。不得修改 Newton、Poisson、能量容差、Beam、Tail、碰撞或物理边界。

### 30.3 当前停止条件

在明确 $W_F\rightarrow W_J$ 的第一处离散不一致前：

```text
禁止 F11；
禁止 J1 MPI 门；
禁止 J2/J3；
禁止能量投影、电流缩放、field scaling 或全局补偿。
```
