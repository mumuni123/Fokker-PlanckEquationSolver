# VPFP F10 情形 B：`W_F-W_J` 根因定位与严格修复实施方案

> **方案性质**：本方案是对现有 `docs/VPFP_能量根因链最终核查与逐阶段修复实施方案.md` 第 11 节“情形 B”的补充执行规格。
> **适用状态**：F1–F9 已执行，F10 单 rank 的 `smooth-background` PASS，而 `smooth-perturbed-background` 以 `failure_code=75` 失败，并已经确认属于情形 B。
> **执行对象**：GPT-5.6 Luna、DeepSeek-v4-Flash 或能力相近的自动编码模型。
> **最高原则**：执行模型不得自行换算法、不得自行“优化”、不得根据经验猜修复。只能按照本文的阶段顺序机械执行。若实际源码与本文前提不一致，立即停止并报告，不得自行设计替代方案。

## 文档编排说明（本次整理新增）

本次整理**不删除、不改写、不替换**本文既有的根因结论、数学公式、代码片段、阶段门、禁止项或历史数值。原有第 0--55 节全部保留在下方；本节和文末第 56 节只增加：

1. 将原有零散章节归并为便于执行的五个大节；
2. 明确各大节的输入、产出、停止条件和后继关系；
3. 补齐每个实际测试阶段的构建、运行和结果读取命令。

发生冲突时，原有 B0--B8、阶段门和禁止项正文优先级高于本次导航与命令附录。每次运行必须使用新的输出目录；不得删除、覆盖或清空已有 `output` 结果。

## 五大执行节导航（本次整理新增）

| 大节 | 原有章节 | 执行目的 | 唯一产出 | 不通过时的停止位置 |
|---|---|---|---|---|
| A. 事实与边界 | §0--§7、§47--§51、§55 | 固定 F10 情形 B 的已知事实、物理边界与禁止项 | 只读基线与 seam 预测公式 | 不修改算法；回报源码/数据矛盾 |
| B. B1 根因确认 | B0 与 §9--§15 | 验证 $R_{FJ}$ 是否由 periodic seam / nonperiodic pairing gather 的解析余额完全解释 | `B1 PASS` 或 `B1 FAIL` | B1 FAIL 时禁止 B2--B8 |
| C. 加权伴随修复 | B2 与 §16--§42 | 只将 u-force cell field 改为现有 periodic x-current map 的 weighted transpose $G^*E$ | 唯一 helper、Newton 接线与最终诊断接线 | 仅修改 B2--B6 白名单；禁止改 current/Poisson |
| D. 单元与 J1 回归 | B4、B7、B8 与 §28--§46 | 先证实 seam adjoint 恒等式，再验证旧 J0 门和两个 F10 J1 case | J0-E2、全部 J0、F10 结果 | 依据 §46 重新进入 A/B/C/D 决策树 |
| E. 交付与文档 | §52--§54 | 记录阶段报告；仅在 B8 PASS 后更新根因文档 | `docs/执行情况.md` 与后续根因记录 | B8 前不得提前写修复成功 |

### 强制顺序

```text
B0
 ↓
B1 诊断与 seam 预测验收
 ├─ FAIL：停止；不得修改 flux / pairing field
 └─ PASS
     ↓
   B2--B3 helper 设计与实现
     ↓
   B4 isolated J0-E2
     ↓
   B5 Newton candidate 接线
     ↓
   B6 final diagnostic 接线
     ↓
   B7 全部 J0 回归
     ↓
   B8 F10 两个单 rank case
     ↓
   §46 重新分流
```

**B2--B6 是一个实现组，不是可任意拆开验收的独立物理阶段；B4、B7、B8 才是该组对应的三道测试门。**

---

# 0. 结论先行

当前 F10 的主问题已经不是原方案 F3、F4、F5 所修复的 u 方向通量、midpoint mass 或 Hamiltonian velocity。

当前最新结果为：

```text
W_u = 28.891863089349869
W_F = 28.891863089349886
W_J = 22.224220061661196

R_uJ = 6.6676430276886727 J/m²
R_FJ = 6.6676430276886904 J/m²
R_PJ = -5.2950189122213942e-6 J/m²
```

并且：

$$
W_u-W_F
\simeq
-1.7\times10^{-14}\ {\rm J/m^2}.
$$

因此 u-face work 与 cell-force work 已经闭合，当前首个确定的不闭合对象是：

$$
W_F-W_J.
$$

仓库的根因链第 9 节也已经把问题收窄到 `x_flux_rate`、`charge_current_face`、face averaging、periodic seam、first/last cell 映射以及 face ownership，而明确禁止修改 Newton 和 Poisson。 F10 执行报告给出了完全一致的判断。

本次源码核查进一步说明：

> **当前最值得首先验证的、而且从现有离散公式可以直接推导出的 B 情形根因，是：J1 的 x 通量采用 periodic seam，但 u-force 使用的** **`e_pair_cell`** **却仍采用适用于普通非周期相邻面的简单算术平均。二者不是同一个 face-to-cell 离散算子的严格转置。**

当前 `evaluate_local_residual()` 在 x 方向明确把最后一个 cell 与第一个 cell 周期连接，并且 `iface=0` 和 `iface=nx` 都表示周期接缝。

与此同时，当前 `advance_joint_midpoint()` 在构造 Poisson pairing face 后仍执行：

```cpp
e_pair_cell[ix] =
    0.5 *
    (pairing_face[ix] +
     pairing_face[ix + 1]);
```

然后把这个 field 用于 J1 u-force。

但是 `OpenElectrostaticSolver::build_potential_pairing_field()` 明确是一个**非周期 Poisson 的 face dual**，两个全局物理 endpoint 各自使用半 cell 权重，而且左、右端 pairing face 是两个独立的 face 值。

因此：

**periodic x current operator 与当前 naive face-to-cell field gather 并不是互为离散转置。**

这正好可以产生现在观察到的：

$$
W_F\neq W_J.
$$

---

# 1. 不得改变的物理模型

本方案不得以“修 J1”为理由修改最终生产物理模型。

最高物理规格仍然是：

```text
docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
```

该文档明确规定：

1. `bulk` 和 `tail` 是同一个背景电子物种的两种数值表示；
2. bulk-to-tail 和 tail-to-bulk 只是内部表示转换；
3. 场始终由**非周期 Gauss/Poisson**约束；
4. Beam 保持开放注入和开放流出；
5. 禁止周期电场、全局能量补丁、电流缩放、强制平均电场为零以及无记账裁剪。

因此本轮绝对禁止修改：

```text
src/open_electrostatic_solver.cpp
src/open_electrostatic_solver.h

src/beam_pic.*
src/background_tail_pic.*
src/bulk_tail_converter.*
src/tail_bulk_return.*
src/cylindrical_fp_collision.*
src/background_tail_collision*
src/hybrid_collision_step.*
src/open_boundary.*
```

尤其：

> **不得把生产 Poisson 改成周期 Poisson。**

> **不得为了让 J1 PASS 而把最终背景边界改成 periodic。**

当前 J1 的 periodic x 只允许作为隔离的 manufactured-test 拓扑存在。

真正生产模型在后续 J2 必须回到开放 background + 非周期 Poisson。

---

# 2. 本轮不得重新怀疑已经通过的部分

执行模型必须接受以下结论，不得重新从 F1 开始“优化”。

## 2.1 u-flux geometry 已通过

当前 `Species::f` 是 cell-integrated mass，F3 已经修复重复 `dx × ring`。

生产 u 通量为：

$$
F^u_{i,j+1/2,k}
=
a_i
\frac12
\left(
\frac{\bar M_{i,j,k}}{\Delta u_j}
+
\frac{\bar M_{i,j+1,k}}{\Delta u_{j+1}}
\right).
$$

不得重新加入：

```cpp
sg.dx
vg.uperp_ring_areas[k]
```

---

## 2.2 midpoint mass 已通过

必须继续使用：

$$
\bar M
=
\frac{M^n+M^{n+1}}{2}.
$$

不得改回 candidate-only。

---

## 2.3 `vH` 已通过直接 u-energy adjoint

生产 `build_hamiltonian_velocity()` 已按 u-face trace 的严格转置构造。

当前已经有：

$$
W_u-W_F
\simeq
10^{-14}\ {\rm J/m^2}.
$$

因此本轮：

```text
不得修改 build_hamiltonian_velocity()
不得换成 vg.vx
不得换成 vg.vx_energy_conjugate_cell
不得换成 c*u/gamma
```

---

## 2.4 Poisson scalar identity 不是当前根因

当前：

$$
|R_{PJ}|
\ll
|R_{FJ}|.
$$

而 Poisson pairing helper 自身已经有明确的有限体积 dual 定义。

所以：

```text
不得修改 Poisson stencil
不得修改 phi 边界
不得缩放 pairing field
不得重新定义 field energy
```

---

# 3. B 情形真正需要检查的离散链

定义 cell force current：

$$
J_i^F
=
\frac{q_e}{\Delta x}
\sum_{j,k}
v^H_{j,k}
\bar M_{i,j,k}.
$$

代码中电子电荷取：

$$
q_e=-e.
$$

当前内部 x-face 的中心通量自动给出：

$$
J^C_{i+1/2}
=
\frac12
\left(
J_i^F+J_{i+1}^F
\right).
$$

F9 的 J0-E 测试确实验证了这个关系。

但需要特别注意：

当前 J0-E 只遍历：

```cpp
for (int iface = 1; iface < nx; ++iface)
```

也就是说它只检查**内部 face**，没有检查 periodic seam。

所以：

> **J0-E PASS 并不能证明 first/last cell 与 periodic seam 的全局 adjoint 关系。**

这正是当前测试体系留下的盲区。

---

# 4. 当前 periodic seam 的实际离散

设全局共有 N 个 cell。

Poisson pairing face 记为：

$$
\mathcal E_f,
\qquad
f=0,1,\ldots,N.
$$

其中：

```text
f = 0    左物理 endpoint
f = N    右物理 endpoint
```

当前 J1 periodic x flux 在接缝处满足：

$$
J^C_0
=
J^C_N
=
\frac12
\left(
J^F_{N-1}+J^F_0
\right).
$$

这是因为当前代码把：

```text
iface = 0
iface = nx
```

都映射到：

```text
left  = last cell
right = first cell
```

作为同一个周期 face。

---

# 5. Poisson face 功使用的权重

当前 Poisson pairing 的 face inner product 使用：

$$
\langle \mathcal E,J\rangle_f
=
\Delta x
\left[
\frac12\mathcal E_0J_0
+
\sum_{f=1}^{N-1}\mathcal E_fJ_f
+
\frac12\mathcal E_NJ_N
\right].
$$

两个物理 endpoint 是半权重。

这是 `build_potential_pairing_field()` 的设计前提。源码注释已经明确指出 endpoint face weight 是：

```text
dx/2
```

并以此使：

$$
-\Delta t
\langle
E_{\rm pair},
J
\rangle_f
$$

与 Poisson potential-charge work 相同。

因此这部分**不能修改**。

---

# 6. 当前 `e_pair_cell` 为什么不再是正确的 adjoint

当前代码使用：

$$
E^{\rm naive}_{i}
=
\frac12
\left(
\mathcal E_i+\mathcal E_{i+1}
\right).
$$

对于普通内部 cell，这没有问题。

但是对于周期 seam：

第一 cell 当前用：

$$
E^{\rm naive}_0
=
\frac12
\left(
\mathcal E_0+\mathcal E_1
\right),
$$

最后一 cell 当前用：

$$
E^{\rm naive}_{N-1}
=
\frac12
\left(
\mathcal E_{N-1}+\mathcal E_N
\right).
$$

这实际上对应的是一种**非周期、两个 endpoint 分离的 face-to-cell gather**。

但当前 x current 却在这两个 endpoint 上使用**同一个 periodic seam current**：

$$
J^C_0=J^C_N.
$$

因此：

```text
x current operator G
```

和

```text
当前 face -> cell force field
```

不是：

$$
G
\quad\text{与}\quad
G^*.
$$

---

# 7. 当前 B 残差可以直接推导出的 seam 预测公式

在保持当前 naive cell field 时，可以严格推导：

$$
W_F^{\rm naive}
-
W_J
=
\frac{
\Delta t\Delta x
}{4}
\left(
\mathcal E_0-\mathcal E_N
\right)
\left(
J_0^F-J_{N-1}^F
\right).
$$

定义：

$$
R_{\rm seam}^{\rm pred}
=
\frac{
\Delta t\Delta x
}{4}
\left(
\mathcal E_0-\mathcal E_N
\right)
\left(
J_0^F-J_{N-1}^F
\right).
$$

如果当前分析正确，则旧 F10 必须满足：

$$
R_{FJ}
=
W_F^{\rm naive}-W_J
\simeq
R_{\rm seam}^{\rm pred}.
$$

当前结果已经知道：

$$
R_{FJ}
=
6.6676430276886904\ {\rm J/m^2}.
$$

但是当前 `.result` 没有同时保存：

```text
pairing_face_left
pairing_face_right
force_current_first_cell
force_current_last_cell
```

所以**现在不能跳过数值验证直接修改算法**。

必须先执行 B0。

---

# 8. 阶段 B0：冻结修改范围

## B0.1 允许修改文件

本轮 B0–B7 只允许修改：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp

src/vpfp_integrator.h
src/vpfp_integrator.cpp

tests/joint_phase_space_midpoint_unit_test.cpp
tests/joint_phase_space_midpoint_energy_test.cpp

docs/执行情况.md
```

只有所有 B 阶段 PASS 后才允许更新：

```text
docs/VPFP能量问题根因链_已确认结论.md
docs/VPFP_能量根因链最终核查与逐阶段修复实施方案.md
```

---

## B0.2 禁止修改

```text
OpenElectrostaticSolver
Newton tolerance
GMRES
line search
energy threshold
dt
Beam
Tail
collision
conversion
return
open boundary
```

---

## B0.3 Git 基线

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

不得：

```text
git reset --hard
git clean
删除未提交修改
覆盖未知修改
```

如果通过 SSH 操作集群，还必须遵守仓库 `ssh远程连接协议.md`：远程写操作需要批准，工作目录只能位于指定项目路径。

---

# 9. 阶段 B1：只增加 seam 根因诊断，不改变算法

**这一阶段绝对不能改变任何 flux、field 或 residual。**

目的只有一个：

验证：

$$
R_{FJ}
\simeq
R_{\rm seam}^{\rm pred}.
$$

---

## B1.1 增加诊断字段

在 `VpfpStepResult` 增加：

```cpp
double joint_midpoint_pairing_face_left;
double joint_midpoint_pairing_face_right;

double joint_midpoint_force_current_first_cell;
double joint_midpoint_force_current_last_cell;

double joint_midpoint_naive_force_current_work;

double joint_midpoint_seam_predicted_residual;
double joint_midpoint_seam_prediction_error;
```

初始化全部为 `0.0`。

不得复用已有其他字段。

---

# 10. B1.2 独立构造 cell force current

在最终 accepted candidate 的能量诊断处，已经有：

```cpp
hamiltonian_velocity
m_old
candidate
```

构造：

```cpp
std::vector<double> force_current_cell(
    static_cast<size_t>(grid_.nx_local), 0.0);
```

对每个 local cell：

```cpp
long double sum = 0.0L;

for (int q = 0; q < nq; ++q) {
    const double midpoint_mass =
        0.5 *
        (m_old[base + q] +
         candidate[base + q]);

    sum +=
        static_cast<long double>(
            hamiltonian_velocity[q]) *
        static_cast<long double>(
            midpoint_mass);
}

force_current_cell[ix] =
    (-Const::qe) *
    static_cast<double>(sum) /
    grid_.dx;
```

必须除：

```cpp
grid_.dx
```

因为定义是 current density。

---

# 11. B1.3 提取四个全局 endpoint 值

需要得到：

```text
E_pair_left
E_pair_right
J_force_first
J_force_last
```

本地数组：

```cpp
double endpoint_local[4] = {
    0.0, 0.0, 0.0, 0.0
};
```

rank 0：

```cpp
endpoint_local[0] = final_pairing_face.front();
endpoint_local[2] = force_current_cell.front();
```

最后一个 rank：

```cpp
endpoint_local[1] = final_pairing_face.back();
endpoint_local[3] = force_current_cell.back();
```

然后：

```cpp
double endpoint_global[4] = {};

MPI_Allreduce(
    endpoint_local,
    endpoint_global,
    4,
    MPI_DOUBLE,
    MPI_SUM,
    MPI_COMM_WORLD);
```

映射：

```text
endpoint_global[0] = E_left
endpoint_global[1] = E_right
endpoint_global[2] = J_first
endpoint_global[3] = J_last
```

---

# 12. B1.4 计算 seam 预测残差

严格按照：

```cpp
const double seam_pred =
    0.25 *
    dt *
    grid_.dx *
    (E_left - E_right) *
    (J_first - J_last);
```

即：

$$
R_{\rm seam}^{\rm pred}
=
\frac{
\Delta t\Delta x
}{4}
(\mathcal E_0-\mathcal E_N)
(J_0^F-J_{N-1}^F).
$$

保存：

```cpp
result.joint_midpoint_seam_predicted_residual =
    seam_pred;
```

然后：

```cpp
result.joint_midpoint_seam_prediction_error =
    result.joint_midpoint_force_charge_residual -
    seam_pred;
```

---

# 13. B1.5 测试输出必须增加

`joint_phase_space_midpoint_energy_test.cpp` 输出：

```text
joint_midpoint_pairing_face_left=
joint_midpoint_pairing_face_right=

joint_midpoint_force_current_first_cell=
joint_midpoint_force_current_last_cell=

joint_midpoint_force_charge_residual=
joint_midpoint_seam_predicted_residual=
joint_midpoint_seam_prediction_error=
```

---

# 14. B1.6 运行旧算法一次

只运行：

```bash
./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result output/b1_seam_diagnostic.result
```

不要期待测试 PASS。

仍允许：

```text
failure_code=75
```

---

# 15. B1 阶段门

计算：

$$
S=
\max
\left(
1,
|R_{FJ}|,
|R_{\rm seam}^{\rm pred}|
\right).
$$

定义：

$$
\tau_{\rm seam}
=
8192
\epsilon_{\rm mach}
S.
$$

要求：

$$
\left| R_{FJ}
-
R_{\rm seam}^{\rm pred}
\right|
\le
\tau_{\rm seam}.
$$

### 如果 PASS

结论记录为：

```text
B root cause confirmed:
periodic x seam / nonperiodic pairing-face
face-to-cell adjoint mismatch
```

进入 B2。

### 如果 FAIL

**立即停止。**

不得执行下面 B2–B7。

输出：

```text
B1 seam prediction failed.
Do not modify flux or pairing field.
Report E_left/E_right/J_first/J_last/R_FJ/R_seam_pred.
```

执行模型不得自行寻找另一个修复。

---

# 16. 阶段 B2：实现真正的 periodic x weighted adjoint

如果 B1 PASS，才执行本阶段。

本阶段核心原则：

> **不修改** **`charge_current_face`** **去追随 force current。**

而是：

> **使 u-force 使用的 cell field 成为当前 x current operator 在 Poisson face quadrature 下的严格离散转置。**

---

# 17. 为什么正确修复不是改 `J_charge`

当前：

$$
J^{C}_{i+1/2}
$$

直接来自实际 `x_flux_rate`。

它同时控制：

```text
mass_delta_x
charge continuity
Poisson charge evolution
```

因此它是**生产连续性 authority**。

不得写：

```cpp
charge_current_face =
    some_corrected_force_current;
```

不得：

```cpp
J *= correction_factor;
```

不得重新构造另一个 current。

否则就会重新回到历史根因：

```text
连续性电流
能量电流
```

各自独立。

---

# 18. 正确的 weighted transpose

定义：

$$
G:
J^F_{\rm cell}
\rightarrow
J^C_{\rm face}.
$$

当前 periodic centered operator 满足：

内部：

$$
J^C_f
=
\frac12
\left(
J^F_{f-1}+J^F_f
\right).
$$

周期接缝：

$$
J^C_0
=
J^C_N
=
\frac12
\left(
J^F_{N-1}+J^F_0
\right).
$$

Poisson face inner product 有：

```text
内部 face weight = 1
endpoint weight  = 1/2
```

因此 u-force field 必须使用：

$$
G^*E_{\rm pair}.
$$

---

# 19. 第一 cell 的正确 field

不能再使用：

$$
\frac12
(\mathcal E_0+\mathcal E_1).
$$

必须使用：

$$
E^{G^*}_{{\rm cell},0}
=
\frac12\mathcal E_1
+
\frac14
\left(
\mathcal E_0+\mathcal E_N
\right).
$$

---

# 20. 普通内部 cell

对于：

$$
1\le i\le N-2,
$$

继续使用：

$$
E^{G^*}_{{\rm cell},i}
=
\frac12
\left(
\mathcal E_i+\mathcal E_{i+1}
\right).
$$

---

# 21. 最后一个 cell

不能再使用：

$$
\frac12
(\mathcal E_{N-1}+\mathcal E_N).
$$

必须使用：

$$
E^{G^*}_{{\rm cell},N-1}
=
\frac12\mathcal E_{N-1}
+
\frac14
\left(
\mathcal E_0+\mathcal E_N
\right).
$$

---

# 22. 这个公式必须满足的恒等式

修改后的 cell field 必须保证：

$$
\Delta x \sum_{i=0}^{N-1} E^{G^*}_{{\rm cell},i} J_i^F
\Delta x
\left[
\frac12\mathcal E_0J^C_0
+
\sum_{f=1}^{N-1}
\mathcal E_fJ^C_f
+
\frac12\mathcal E_NJ^C_N
\right].
$$

即：

$$
W_F=W_J.
$$

这个关系是**代数恒等式**。

不得通过 tolerance 才“近似成立”。

---

# 23. 阶段 B3：新增唯一 helper，不允许两处手写公式

在：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
```

新增唯一 helper。

推荐函数名：

```cpp
build_periodic_x_adjoint_cell_field
```

接口：

```cpp
static bool build_periodic_x_adjoint_cell_field(
    const SpatialGrid& sg,
    const std::vector<double>& pairing_face,
    int mpi_rank,
    int mpi_size,
    std::vector<double>& pairing_cell);
```

---

# 24. B3.1 helper 输入检查

必须检查：

```cpp
sg.nx_local > 0
sg.nx_global >= 2

pairing_face.size()
    == sg.nx_local + 1

mpi_rank >= 0
mpi_rank < mpi_size
mpi_size >= 1
```

全部 pairing face 必须：

```cpp
std::isfinite(...)
```

失败：

```cpp
pairing_cell.clear();
return false;
```

不得 fallback。

---

# 25. B3.2 获取两个全局 Poisson endpoint pairing face

```cpp
double endpoint_local[2] = {0.0, 0.0};

if (mpi_rank == 0)
    endpoint_local[0] = pairing_face.front();

if (mpi_rank == mpi_size - 1)
    endpoint_local[1] = pairing_face.back();

double endpoint_global[2] = {0.0, 0.0};

MPI_Allreduce(
    endpoint_local,
    endpoint_global,
    2,
    MPI_DOUBLE,
    MPI_SUM,
    MPI_COMM_WORLD);
```

定义：

```cpp
const double e_left  = endpoint_global[0];
const double e_right = endpoint_global[1];
```

检查 finite。

---

# 26. B3.3 构造 local cell field

```cpp
pairing_cell.assign(
    static_cast<size_t>(sg.nx_local),
    0.0);
```

遍历：

```cpp
for (int ix = 0; ix < sg.nx_local; ++ix) {
    const int ig = sg.ix_start + ix;

    if (ig == 0) {
        pairing_cell[ix] =
            0.5 * pairing_face[ix + 1]
            + 0.25 * (e_left + e_right);
    }
    else if (ig == sg.nx_global - 1) {
        pairing_cell[ix] =
            0.5 * pairing_face[ix]
            + 0.25 * (e_left + e_right);
    }
    else {
        pairing_cell[ix] =
            0.5 *
            (pairing_face[ix]
             + pairing_face[ix + 1]);
    }
}
```

随后全部 finite-check。

---

# 27. helper 注释必须写明

必须加入含义等价于下面内容的注释：

```cpp
// J1 TEST TOPOLOGY ONLY.
//
// The J1 x operator is periodic: global faces 0 and Nx represent
// the same seam current, while the OpenElectrostaticSolver pairing
// field retains two distinct non-periodic physical endpoint faces
// with half-cell quadrature weights.
//
// Therefore the u-force cell field must be the exact weighted
// transpose G* of the periodic centered x-current map G.
//
// Do not replace the first/last-cell formulas with
// 0.5*(E_left_face + E_right_face).
//
// This is NOT the production open-boundary rule.
// J2 must replace the periodic seam operator with the real
// OpenBackgroundBoundary operator and derive its corresponding G*.
```

这一注释非常重要。

---

# 28. 阶段 B4：先增加 isolated adjoint 单元测试

**在生产 Newton 使用 helper 之前先测试 helper。**

修改：

```text
tests/joint_phase_space_midpoint_unit_test.cpp
```

不要新建 CMake target。

直接加入现有 J0 测试。

命名：

```text
J0-E2 periodic-seam weighted adjoint
```

---

# 29. B4.1 测试状态不得具有平凡对称性

构造 deterministic positive midpoint mass。

必须：

```text
first cell current != last cell current
```

并构造人工 pairing face，使：

```text
pairing_face[0] != pairing_face[N]
```

不得使用：

```text
E_left == E_right
J_first == J_last
```

否则 seam 错误会被自动隐藏。

---

# 30. B4.2 独立构造 force current

用 production：

```cpp
build_hamiltonian_velocity()
```

和测试的 midpoint mass 构造：

$$
J_i^F.
$$

---

# 31. B4.3 使用 production x bundle

调用：

```cpp
build_periodic_center_flux(...)
```

取得：

```text
bundle.x_flux_rate
bundle.charge_current_face
```

不得重新构造 `J_charge`。

---

# 32. B4.4 计算新的 adjoint cell field

调用刚新增的：

```cpp
build_periodic_x_adjoint_cell_field(...)
```

得到：

```text
e_adjoint_cell
```

---

# 33. B4.5 独立计算两边功

force side：

$$
W_F
=
\Delta x
\sum_i
E^{G^*}_i
J_i^F.
$$

face side：

$$
W_J
=
\Delta x
\left[
\frac12
\mathcal E_0J_0^C
+
\sum_{f=1}^{N-1}
\mathcal E_fJ_f^C
+
\frac12
\mathcal E_NJ_N^C
\right].
$$

必须分别独立求和。

不得调用同一个 helper 返回两边。

---

# 34. B4.6 必须同时验证旧 naive gather 会失败

额外构造：

$$
E_i^{\rm naive}
=
\frac12
(\mathcal E_i+\mathcal E_{i+1}).
$$

计算：

$$
W_F^{\rm naive}.
$$

然后验证：

$$
W_F^{\rm naive}-W_J
$$

等于：

$$
R_{\rm seam}^{\rm pred}.
$$

这一步用于证明测试真的覆盖了此前遗漏的 seam。

---

# 35. B4 阶段门

要求：

```text
periodic_seam_weighted_adjoint_pass=1
periodic_seam_prediction_pass=1
```

并且：

$$
\frac{
|W_F-W_J|
}{
\max(1,|W_F|,|W_J|)
}
\le
8192\epsilon_{\rm mach}.
$$

同时旧 naive mismatch 必须明显非零。

如果 helper 测试不通过：

> **不得进入 B5。**

---

# 36. 阶段 B5：在 Newton candidate 中真正使用 `G*`

当前 candidate evaluation 在 Poisson solve 后调用：

```cpp
build_potential_pairing_field(...)
```

这一部分不变。

当前随后执行的：

```cpp
for (int ix = 0; ix < grid_.nx_local; ++ix)
    e_pair_cell[ix] =
        0.5 *
        (pairing_face[ix] +
         pairing_face[ix + 1]);
```

必须删除。当前源码确实仍采用这一简单平均。

替换为：

```cpp
std::vector<double> e_pair_cell;

const bool adjoint_ok =
    JointPhaseSpaceMidpointOperator::
        build_periodic_x_adjoint_cell_field(
            grid_,
            pairing_face,
            mpi_rank,
            mpi_size,
            e_pair_cell);

local_ok = local_ok && adjoint_ok;
```

---

# 37. B5.1 无 fallback

如果：

```cpp
adjoint_ok == false
```

candidate evaluation 失败。

禁止：

```cpp
e_pair_cell = eval_fields.Ex;
```

禁止：

```cpp
e_pair_cell =
    0.5 * (pairing_face_left + pairing_face_right);
```

禁止：

```cpp
fallback to old averaging
```

---

# 38. B5.2 `evaluate_local_residual()` 不修改

仍然调用：

```cpp
JointPhaseSpaceMidpointOperator::
    evaluate_local_residual(
        ...,
        e_pair_cell,
        ...
    );
```

不得改：

```text
x_flux_rate
u_flux_rate
midpoint mass
vH
charge_current_face
```

本轮唯一改变的是：

> **u-force 接收到的 cell pairing field。**

---

# 39. B5.3 preconditioner 自动使用新 field

当前：

```cpp
accepted_e_local
```

来自 candidate evaluate。

因此只要：

```cpp
e_local_out = e_pair_cell;
```

保持不变，preconditioner 会自然接收到新的 `G*E`。

不要修改 preconditioner。

---

# 40. 阶段 B6：最终能量诊断也必须使用同一个 helper

当前最终能量诊断再次手工构造：

```cpp
final_pairing_cell[ix] =
    0.5 *
    (final_pairing_face[ix] +
     final_pairing_face[ix + 1]);
```

当前源码确实如此。

必须删除。

改为：

```cpp
std::vector<double> final_pairing_cell;

const bool final_adjoint_ok =
    JointPhaseSpaceMidpointOperator::
        build_periodic_x_adjoint_cell_field(
            grid_,
            final_pairing_face,
            mpi_rank,
            mpi_size,
            final_pairing_cell);

if (!final_adjoint_ok) {
    result.failure_code = 71;
    result.failure_stage =
        "joint_midpoint_final_x_adjoint_field";
    return result;
}
```

---

# 41. B6.1 `W_F` 必须代表实际生产 force field

修改之后：

```text
joint_midpoint_force_current_work
```

必须使用：

```text
final_pairing_cell
```

也就是实际 production：

$$
G^*E_{\rm pair}.
$$

这样才允许比较：

$$
W_u-W_F,
$$

以及：

$$
W_F-W_J.
$$

---

# 42. B6.2 可以保留旧 naive work 作为诊断

为了审计历史根因，可以额外构造：

```text
joint_midpoint_naive_force_current_work
```

只用于诊断。

它不得参与：

```text
residual
Newton
acceptance
energy gate
```

预期修改后仍有：

$$
W_F^{\rm naive}-W_J
\simeq
R_{\rm seam}^{\rm pred}.
$$

但真正 production 应为：

$$
W_F^{G^*}-W_J
\simeq0.
$$

这会构成非常清晰的修复证据。

---

# 43. 阶段 B7：重新运行 J0

必须重新编译：

```bash
cmake --build build \
  --target joint_phase_space_midpoint_unit_test \
  -j4
```

运行：

```bash
./build/joint_phase_space_midpoint_unit_test
```

旧六门全部必须继续：

```text
j0_a_cell_mass_pass=1
j0_b_u_flux_geometry_pass=1
j0_c_midpoint_consistency_pass=1
j0_d_u_work_force_adjoint_pass=1
j0_e_x_current_force_adjoint_pass=1
j0_f_poisson_pairing_pass=1
```

并新增：

```text
j0_e2_periodic_seam_weighted_adjoint_pass=1
```

不得把旧 J0-E 删除或替换。

---

# 44. 阶段 B8：重新运行 F10 单 rank

先运行：

```bash
cmake --build build \
  --target joint_phase_space_midpoint_energy_test \
  -j4
```

然后：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result output/b8_smooth-background_n1.result
```

再运行：

```bash
yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result output/b8_smooth-perturbed-background_n1.result
```

---

# 45. B8 必须检查的数值链

## 45.1 平衡态

必须继续：

```text
status=PASS
accepted=1
finite=1
gauss_ok=1
converged=1
failure_code=0
```

任何回归都立即停止。

---

## 45.2 非平凡场首先检查 u-work

要求：

$$
W_u-W_F
$$

仍为 roundoff。

如果这一项因为 B 修改明显变大：

> **立即停止。**

说明新的 field 没有被 u-flux 和 force diagnostic 一致使用。

不得继续调电流。

---

## 45.3 然后检查真正的 B 根因

要求：

$$
R_{FJ}
=
W_F-W_J
$$

达到 roundoff。

建议：

$$
\frac{
|W_F-W_J|
}{
\max(1,|W_F|,|W_J|)
}
\le
10^{-12}.
$$

理论上应远好于这个阈值。

---

## 45.4 再检查：

$$
R_{uJ}
=
W_u-W_J.
$$

因为：

$$
W_u\simeq W_F
$$

且：

$$
W_F\simeq W_J,
$$

所以必须得到：

$$
R_{uJ}\simeq0.
$$

---

# 46. B8 后严格重新进入第 11 节决策树

这一点非常重要。

B 修完后**不要因为代码修改完成就自动进入 F11。**

必须重新分类。

---

## B8 情形 1：F10 全部 PASS

如果：

```text
R_uJ PASS
R_PJ PASS
R_E PASS
accepted=1
failure_code=0
```

则：

> B 情形修复成功。

这时才允许进入原方案：

```text
J1 MPI gate
```

---

## B8 情形 2：`R_uJ` 已小，但 `R_PJ` 仍大

则新的失败已经变成：

> **第 11 节情形 A。**

立即停止 B。

不得继续修改 `G*`。

不得修改 Poisson。

必须另外执行 A 情形的 `charge continuity ↔ Poisson pairing` 审计。

尤其当前旧 F10 有：

$$
R_{PJ}
=
-5.2950189122213942\times10^{-6}
\ {\rm J/m^2}.
$$

它虽然比原 B 主残差小约六个数量级，但 B 修复以后必须重新判断它是否达到真正的 Poisson work tolerance。

不能预设它一定自动消失。

---

## B8 情形 3：所有能量 identity PASS，但出现 code 73、74、76

这时才属于原方案：

```text
情形 D
```

才允许进入 F11 positivity。

在此之前绝对不能改 Newton。

---

# 47. 为什么本方案不能直接把 J1 改成开放边界

不要采用这样的“修复”：

```text
既然 periodic seam 有问题，
那就直接把 J1 x boundary 改成 open。
```

原因是：

J1 当前的任务是先隔离验证共同 phase-space operator。

真正的：

```text
OpenBackgroundBoundary
reservoir inflow
physical outflow
boundary kinetic flux
boundary charge current
```

属于 J2。

如果现在直接把 J2 的开放边界引入 J1，那么一旦测试失败，将无法判断残差来自：

```text
joint x/u core
还是
open reservoir boundary
```

因此本方案修复的只是：

> **J1 manufactured periodic x topology 下的正确离散 adjoint。**

---

# 48. 更重要：这个 periodic `G*` 绝对不能直接变成最终生产边界

最终物理模型仍然是开放 background。

所以在 J2：

当前 periodic：

$$
G_{\rm periodic}
$$

必须被真正的：

$$
G_{\rm open}
$$

替代。

同时必须重新推导：

$$
G_{\rm open}^*.
$$

J2 中不能继续使用本方案的：

$$
\frac14
(\mathcal E_0+\mathcal E_N)
$$

periodic seam 项。

因为最终物理模型中：

```text
左 endpoint
右 endpoint
```

不是同一个拓扑 face。

这点必须在代码注释中明确写死，避免较弱模型以后把 J1 测试公式误带进生产模型。

---

# 49. 本轮修改后的完整离散链应该是什么

B 修复以后，J1 的离散链应成为：

```text
M_mid
  ↓
same vH
  ↓
x_flux_rate
  ↓
charge_current_face = G J_force
  ↓
Poisson pairing face
  ↓
E_force_cell = G* E_pair_face
  ↓
u_flux_rate
  ↓
kinetic work
```

于是第一条恒等式：

$$
W_u=W_F.
$$

第二条恒等式：

$$
W_F=W_J.
$$

第三条：

$$
\Delta U_E
-
W_{\rm electrode}
+
W_J
=0.
$$

因此：

$$
W_u + \Delta U_E
-
W_{\rm electrode}
=0.
$$

这才真正修复了根因文档所说的：

```text
连续性 current
kinetic-work current
Poisson field work
```

不是同一个离散对象的问题。

根因文档明确指出最终方向必须由同一候选 midpoint state、同一 x/u 通量和同一 Poisson pairing 建立共同离散结构，而不是继续局部补丁。

---

# 50. 为什么这不是“能量补丁”

本方案没有执行：

```text
根据 residual 修改 E
根据 residual 修改 f
根据 residual 修改 J
缩放 current
调整 Poisson
修改 kinetic energy
```

它只是把原本已经存在的 x current operator：

$$
G
$$

对应的 force-field gather 从错误的简单平均改为其数学上的严格 weighted transpose：

$$
G^*.
$$

所以它属于：

> **修复离散算子不共轭。**

而不是：

> **事后强迫能量守恒。**

这是两者最重要的区别。

---

# 51. 禁止执行模型做出的错误“简化”

执行 GPT-5.6 Luna / DeepSeek-v4-Flash 时，必须把下面内容原样作为禁止项。

```text
1. 不得修改 OpenElectrostaticSolver。

2. 不得修改 build_potential_pairing_field()。

3. 不得修改 endpoint half-weight。

4. 不得改 charge_current_face 数值来匹配 force current。

5. 不得缩放 J。

6. 不得缩放 E。

7. 不得修改 vH。

8. 不得修改 u_flux_rate。

9. 不得修改 midpoint mass。

10. 不得改变 dt。

11. 不得放宽 1e-8 energy gate。

12. 不得增加 Newton iteration 试图掩盖 code 75。

13. 不得进入 F11。

14. 不得把 negative mass clip 到 0。

15. 不得直接接入 J2 open boundary。

16. 不得恢复旧 PPM/Strang。

17. 不得恢复 Ampere field advance。

18. 不得把生产 Poisson 改成 periodic。

19. 不得把 periodic seam G* 当作最终开放边界公式。

20. B1 seam identity 未验证前，不得实施 B2 以后算法修改。
```

---

# 52. 每个阶段的执行报告格式

每阶段结束必须写入：

```text
docs/执行情况.md
```

格式固定：

````markdown
# 阶段：B?

## 当前阶段

B?

## baseline

commit=
working_tree_clean=

## 实际修改文件

- ...

## 修改内容

1. ...
2. ...

## 明确未修改

- OpenElectrostaticSolver
- Newton
- Beam
- Tail
- collision
- open boundary

## 编译命令

```text
...
```

## 测试命令

```text
...
```

## 数值结果

```text
...
```

## 阶段门

PASS / FAIL

## 判定

...

````

失败时：

> **只记录 FAIL，不得继续下一阶段。**

---

# 53. B 修复全部 PASS 后需要更新的根因文档

只有 B8 真正 PASS 后，才更新：

```text
docs/VPFP能量问题根因链_已确认结论.md
```

第 9 节补充：

```text
F10 情形 B 的 W_F-W_J 已进一步确认来自
J1 periodic centered x-current operator 与旧 naive
face-to-cell pairing-field gather 不互为 weighted adjoint。

旧 J0-E 只验证内部 face，因此未覆盖 periodic seam。

修复方式不是改变 charge current 或 Poisson，而是构造
periodic x-current map G 在 Poisson endpoint half-weight
face inner product 下的严格 transpose G*。
```

并记录实际修复后数值：

```text
W_u=
W_F=
W_J=
R_uF=
R_FJ=
R_uJ=
R_PJ=
R_E=
```

没有真实数值前不得提前写 PASS。

---

# 54. 对原修复方案第 11 节情形 B 的正式补充

原方案只写了：

```text
检查 u trace
检查 vH
检查 x flux
检查 midpoint mass
检查 geometry factors
```

在当前 F1–F10 的新证据下，这已经不够具体。

现在应把情形 B 的执行逻辑明确更新为：

```text
B1:
先验证 W_F-W_J 是否等于 periodic seam adjoint mismatch
的解析预测值。

B2:
若成立，不改变 x flux 和 charge current，
而构造 current operator G 的 weighted transpose G*。

B3:
增加唯一 face-to-cell adjoint helper。

B4:
用独立 J0-E2 测试证明 periodic seam 全局功恒等式。

B5:
Newton candidate 的 u-force 改用 G* E_pair。

B6:
最终 energy diagnostic 同样改用相同 G*。

B7:
全部 J0 回归。

B8:
重新运行单 rank F10。

之后重新进入 A/B/C/D 决策树。
```

---

# 55. 最终判断

基于当前仓库源码和 F10 最新证据，我认为**当前 B 情形不应该继续在 u-flux、****`vH`** **或 midpoint mass 上寻找问题**。

当前最关键的源码不一致是：

```text
x transport:
periodic centered seam

Poisson work:
non-periodic face topology
physical endpoint half weights

force field gather:
naive local 0.5*(left face + right face)
```

这三个对象在 first/last cell 上不是同一个离散伴随系统。

当前 J0 又只验证了内部 face：

$$
J^C_{i+1/2}
=
\frac12
(J_i^F+J_{i+1}^F),
$$

没有验证 periodic seam，因此 J0 全 PASS 与 F10 的：

$$
W_F-W_J
=
6.6676430276886904\ {\rm J/m^2}
$$

并不矛盾。

**因此本轮应首先执行 B1，用**

$$
R_{\rm seam}^{\rm pred}
=
\frac{
\Delta t\Delta x
}{4}
(\mathcal E_0-\mathcal E_N)
(J_0^F-J_{N-1}^F)
$$

**对当前 6.667643 J/m² 缺口做数值确认。**

如果 B1 在 roundoff 内完全解释 `W_F-W_J`，则后续修复就非常明确：

> **保持 Poisson 不变、保持 x flux 不变、保持 charge current 不变，唯一把 J1 u-force 的 face-to-cell gather 改成 periodic x current operator 的严格 weighted transpose。**

这比原方案第 11 节情形 B 的“检查 x flux / vH / midpoint”已经具体到了可以直接机械实现的程度，同时不改变 `开放Beam_非周期场_VPFP_PIC完整重构方案.md` 所规定的最终物理模型。

---

# 56. 执行命令与结果验收附录（本次整理新增）

本附录只补充原有 B0--B8 的执行命令，不更改其阶段门。所有命令在集群项目根目录执行：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
```

远程操作仍必须遵守 `docs/ssh远程连接协议.md`：无远程修改授权时只允许读取、构建已存在代码和运行调度器任务；本文不构成绕过该协议的授权。

## 56.1 共用构建与唯一输出目录

仅在首次构建、`CMakeLists.txt` 变更或新增/修改测试源文件后运行配置命令。不得删除 `build/`：

```bash
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON
```

每个阶段建立独立输出目录：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_${RUN_ID}"
mkdir -p "$OUT"
```

以下命令中的 `$OUT` 均指该独立目录。不得改成既有的 `output/b1_*`、`output/b4_*`、`output/b7_*` 或 `output/b8_*` 目录，也不得添加 `--overwrite-output`。

## 56.2 B0：冻结基线命令

对应原文 B0，只记录，不修改：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b0_${RUN_ID}"
mkdir -p "$OUT"

git status --short | tee "$OUT/git_status.txt"
git rev-parse HEAD | tee "$OUT/git_head.txt"
git diff --check | tee "$OUT/git_diff_check.txt"
```

验收：三个文件均已生成；工作树非空只代表必须记录已有修改，不代表可以 reset、clean 或覆盖。

## 56.3 B1：旧算法 seam 预测诊断命令

对应原文 B1。完成 B1.1--B1.5 的**仅诊断**代码修改后，先构建：

```bash
cmake --build build --target joint_phase_space_midpoint_energy_test -j4
```

运行旧算法的非平凡场 case。`code 75` 在此阶段可预期，因此必须保留退出码但不得因非零退出码跳过 result 检查：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b1_${RUN_ID}"
mkdir -p "$OUT"

set +e
yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/b1_seam_diagnostic.result" \
  > "$OUT/b1_seam_diagnostic.out" \
  2> "$OUT/b1_seam_diagnostic.err"
RC=$?
set -e
printf 'solver_exit_code=%s\n' "$RC" | tee "$OUT/b1_exit_code.txt"

test -s "$OUT/b1_seam_diagnostic.result" || exit 1
sed -n '1,320p' "$OUT/b1_seam_diagnostic.result"
```

必须从 result 读取原文 B1.5 的七项：

```text
joint_midpoint_pairing_face_left
joint_midpoint_pairing_face_right
joint_midpoint_force_current_first_cell
joint_midpoint_force_current_last_cell
joint_midpoint_force_charge_residual
joint_midpoint_seam_predicted_residual
joint_midpoint_seam_prediction_error
```

按原文 B1 阶段门计算 $S$ 和 $\tau_{\rm seam}$。只有 `abs(seam_prediction_error) <= tau_seam` 才可进入 B2；`solver_exit_code != 0` 本身不构成 B1 FAIL。

## 56.4 B2--B3：helper 实现后的编译检查

B2/B3 没有独立数值阶段门。完成原文 §16--§27 的代码修改后，只构建随后 B4 必须使用的既有 J0 target：

```bash
cmake --build build --target joint_phase_space_midpoint_unit_test -j4
```

构建失败时，停止并修复 B2/B3 白名单内的 helper 接口、实现或测试接线；不得运行 B5/B6 或修改生产 current。

## 56.5 B4：isolated periodic-seam weighted-adjoint 命令

对应原文 B4。B4 不新建 CMake target；J0-E2 必须加入既有 `joint_phase_space_midpoint_unit_test --case all`。

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b4_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/b4_periodic_seam_adjoint.result" \
  > "$OUT/b4_periodic_seam_adjoint.out" \
  2> "$OUT/b4_periodic_seam_adjoint.err" || exit 1

sed -n '1,420p' "$OUT/b4_periodic_seam_adjoint.result"
```

验收必须同时满足原文 B4.1--B4.6：

```text
j0_e2_periodic_seam_weighted_adjoint_pass=1
new_weighted_adjoint_relative_error <= 8192 * epsilon_machine
old_naive_mismatch 明显非零
first_cell_current != last_cell_current
pairing_face[0] != pairing_face[N]
```

B4 FAIL 时不得把 helper 接入 Newton。

## 56.6 B5--B6：production 接线后的构建检查

B5/B6 只是在 production J1 candidate 和最终诊断中使用同一个 `G^*` helper，数值门在 B7/B8。完成原文 §36--§42 后：

```bash
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
```

若构建失败，只允许检查 helper 调用、失败码 `71` 路径、`e_local_out` 与 final diagnostic 的同源性；不得修改 `evaluate_local_residual()`、`x_flux_rate`、`u_flux_rate`、`charge_current_face` 或 Poisson。

## 56.7 B7：全部 J0 回归命令

对应原文 B7。B5/B6 接线后必须先通过 J0，才可运行 B8：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b7_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_unit_test -j4

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/b7_j0_all.result" \
  > "$OUT/b7_j0_all.out" \
  2> "$OUT/b7_j0_all.err" || exit 1

sed -n '1,520p' "$OUT/b7_j0_all.result"
```

验收：原文 B7 的六个旧 J0 门均为 `1`，且：

```text
j0_e2_periodic_seam_weighted_adjoint_pass=1
```

任何旧门回归或 J0-E2 FAIL 均禁止进入 B8。

## 56.8 B8：F10 单 rank 回归命令

对应原文 B8。只有 B7 PASS 后运行。两个 case 必须分别保存结果：

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b8_${RUN_ID}"
mkdir -p "$OUT"

cmake --build build --target joint_phase_space_midpoint_energy_test -j4

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-background \
  --result "$OUT/b8_smooth-background_n1.result" \
  > "$OUT/b8_smooth-background_n1.out" \
  2> "$OUT/b8_smooth-background_n1.err" || exit 1

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result "$OUT/b8_smooth-perturbed-background_n1.result" \
  > "$OUT/b8_smooth-perturbed-background_n1.out" \
  2> "$OUT/b8_smooth-perturbed-background_n1.err" || exit 1

sed -n '1,360p' "$OUT/b8_smooth-background_n1.result"
sed -n '1,420p' "$OUT/b8_smooth-perturbed-background_n1.result"
```

验收严格使用原文 §45--§46：

```text
smooth-background: status=PASS, accepted=1, finite=1,
                   gauss_ok=1, converged=1, failure_code=0

smooth-perturbed-background:
  W_u-W_F 为 roundoff；
  abs(W_F-W_J) / max(1, abs(W_F), abs(W_J)) <= 1e-12；
  R_uJ 为 roundoff；
  之后才评估 R_PJ、R_E、accepted 与 failure_code。
```

若 B8 任何条件失败，严格按原文 §46 重新分流；不得直接进入 F11 或 J1 MPI gate。

## 56.9 结果文件与阶段关系速查

| 阶段 | 唯一关键 result | 允许的进展 |
|---|---|---|
| B0 | `git_status.txt`、`git_head.txt`、`git_diff_check.txt` | 记录后进入 B1 |
| B1 | `b1_seam_diagnostic.result` | 仅 seam prediction PASS 才进入 B2 |
| B4 | `b4_periodic_seam_adjoint.result` | J0-E2 PASS 才接入 B5/B6 |
| B7 | `b7_j0_all.result` | 旧 J0 六门和 J0-E2 均 PASS 才进入 B8 |
| B8 | 两个 `b8_*_n1.result` | 按 §46 决定 J1 MPI、情形 A、情形 D 或停止 |
