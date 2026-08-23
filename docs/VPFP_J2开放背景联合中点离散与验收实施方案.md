# VPFP J2：开放背景联合中点离散与验收实施方案

## 0. 文档用途与最高纪律

本文档是 J1 全部通过后，将联合中点 bulk 核心从 periodic manufactured topology 接入最终 open/reservoir background 的机械执行规格。

执行代理必须完整阅读：

```text
docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
docs/VPFP能量问题根因链_已确认结论.md
docs/情形A执行情况.md
docs/ssh远程连接协议.md
```

### 0.1 不可改变的生产物理

```text
background bulk boundary = reservoir/open
field boundary = DIRICHLET_PHI
phi_left = 0
phi_right = 0
Beam = OFF（J2）
Tail = OFF（J2）
conversion = OFF（J2）
return = OFF（J2）
collision = OFF（J2）
```

J2 结束前不得接入 J3a--J3e。

### 0.2 绝对禁止

1. 不得把生产 background 改成 periodic。
2. 不得在物理左/右端点 wrap 到另一端。
3. 不得直接调用 J1 的 `build_periodic_x_adjoint_cell_field()` 作为生产 open 规则。
4. 不得修改 `OpenElectrostaticSolver::solve()`、Poisson stencil 或零端点电势。
5. 不得从密度差反推并覆盖 open x current。
6. 不得把 reservoir 变成域内体源或每步重置边界 cell。
7. 不得用能量投影、场/电流缩放、分布裁剪或零模扣除闭账。
8. 不得放宽 J1 已验证的 phase、Poisson、pairing、总能量和 code-76 门。
9. 测试必须调用生产函数，不得复写公式后宣称生产实现通过。
10. 每次只执行一个阶段，当前阶段未 PASS 立即停止。
11. 所有修改在本地完成；集群测试由用户执行。没有结果时报告“待集群执行”。
12. 阶段报告统一写入 `docs/J2执行情况.md`。

---

## 1. J1 已完成与 J2 的边界

J1 已证明：

$$
W_u=W_F=W_J+O(\epsilon_{\mathrm{mach}})
$$

以及：

$$
\Delta K_{\mathrm{bulk}}+\Delta U_E-W_{\mathrm{electrode}}
=O(\epsilon_{\mathrm{solver}})
$$

但 J1 的 x operator 是 periodic manufactured topology。J2 必须保留 J1 的：

- cell-integrated mass；
- midpoint mass；
- u-face flux 与 Hamiltonian velocity；
- stable incremental charge assembly；
- stable Poisson work identity；
- phase/Poisson/pairing 三重收敛门；
- code-76 最终正性门。

J2 只替换 x boundary topology、对应的边界 ledger 和 open pairing 分解。

---

## 2. 开放 x-flux 的唯一离散定义

### 2.1 方向约定

所有 x face flux 均按正 x 方向为正。对速度槽：

$$
v^H_{j,k}
$$

必须继续使用 J1 已验证的 Hamiltonian velocity，不能改用解析 `vx` 或另一张速度表。

### 2.2 内部 face

对内部 face：

$$
T^x_{i+1/2,j,k}
=\frac{1}{2\Delta x}
\left(M_{i,j,k}^{n+1/2}+M_{i+1,j,k}^{n+1/2}\right)
$$

$$
F^x_{i+1/2,j,k}=v^H_{j,k}T^x_{i+1/2,j,k}
$$

内部 MPI face 必须由两侧相同 midpoint state 构造，并按单一 owner 计入全局 ledger。

### 2.3 左物理 face

若：

$$
v^H_{j,k}>0
$$

则为左侧 incoming characteristic，trace 必须来自：

```cpp
background_boundary.incoming_cell_average(
    PhysicalSide::LEFT, j, k, time_mid, electrons)
```

若：

$$
v^H_{j,k}<0
$$

则为左侧 outgoing characteristic，trace 来自第一个域内 cell 的 midpoint mass：

$$
T^x_{1/2,j,k}=M_{0,j,k}^{n+1/2}/\Delta x
$$

### 2.4 右物理 face

若：

$$
v^H_{j,k}<0
$$

则为右侧 incoming reservoir；若：

$$
v^H_{j,k}>0
$$

则为最后一个域内 cell 的 outgoing trace。

### 2.5 零速度

若 `vH==0`，flux 必须严格为零，不得随机选择 incoming/outgoing。

### 2.6 absorbing

`ABSORBING` incoming trace 为零；outgoing 仍来自域内状态。J2 主验收使用 reservoir/reservoir，同时必须保留 absorbing 单元测试。

---

## 3. 同源边界 ledger

所有 ledger 必须在构造 `F_x` 的同一个 face/velocity 循环内累加，禁止事后从密度或边界 cell 差值反推。

### 3.1 数目 ledger

定义正数 magnitude：

$$
N_{L,\mathrm{in}}
=\Delta t\sum_{v^H>0}F^x_{1/2,j,k}
$$

$$
N_{L,\mathrm{out}}
=-\Delta t\sum_{v^H<0}F^x_{1/2,j,k}
$$

$$
N_{R,\mathrm{in}}
=-\Delta t\sum_{v^H<0}F^x_{N+1/2,j,k}
$$

$$
N_{R,\mathrm{out}}
=\Delta t\sum_{v^H>0}F^x_{N+1/2,j,k}
$$

要求：

$$
N^{n+1}-N^n
=N_{L,\mathrm{in}}+N_{R,\mathrm{in}}
-N_{L,\mathrm{out}}-N_{R,\mathrm{out}}
$$

### 3.2 charge current

$$
J_f=q_e\sum_{j,k}F^x_{f,j,k}
$$

必须使用完整 signed flux，包括 incoming 与 outgoing。

### 3.3 kinetic ledger

使用同一速度槽的离散动能：

$$
F^K_f=\sum_{j,k}K_{j,k}F^x_{f,j,k}
$$

分别输出四个正数 magnitude 的 kinetic in/out。要求 x transport 满足：

$$
\Delta K_x
=K_{L,\mathrm{in}}+K_{R,\mathrm{in}}
-K_{L,\mathrm{out}}-K_{R,\mathrm{out}}
$$

### 3.4 动量 ledger

同一循环输出 longitudinal momentum in/out，用于验证 reservoir drift 和后续 J3 接入；J2 阶段不要求动量守恒，因为 reservoir 是外部动量源，但必须闭合账本。

---

## 4. open current map 的 affine 分解

开放边界 current 不是纯线性 periodic map。写成：

$$
J^{\mathrm{open}}
=G_{\mathrm{ref}}J^{\mathrm{cell}}
+J^{\mathrm{boundary}}
$$

其中 `G_ref` 是非周期 reference map：

- 内部 face：左右 cell current 各一半；
- 左 endpoint：第一个 cell current 的一半；
- 右 endpoint：最后一个 cell current 的一半。

实际 upwind/reservoir endpoint current 与 reference 的差为：

$$
\delta J_L=J^{\mathrm{open}}_L-\frac{1}{2}J^{\mathrm{cell}}_0
$$

$$
\delta J_R=J^{\mathrm{open}}_R-\frac{1}{2}J^{\mathrm{cell}}_{N-1}
$$

该差异是物理 boundary contribution，不得通过修改 endpoint current 消除。

---

## 5. 非周期 reference map 的加权转置

Poisson face inner product 使用内部 face 全权重、物理 endpoint 半权重。`G_ref` 的加权转置必须为：

第一个 cell：

$$
E_0^{\mathrm{ref}*}
=\frac{1}{4}E_{1/2}^{\mathrm{pair}}
+\frac{1}{2}E_{3/2}^{\mathrm{pair}}
$$

普通内部 cell：

$$
E_i^{\mathrm{ref}*}
=\frac{1}{2}
\left(E_{i-1/2}^{\mathrm{pair}}+E_{i+1/2}^{\mathrm{pair}}\right)
$$

最后一个 cell：

$$
E_{N-1}^{\mathrm{ref}*}
=\frac{1}{2}E_{N-3/2}^{\mathrm{pair}}
+\frac{1}{4}E_{N-1/2}^{\mathrm{pair}}
$$

必须先由 isolated weighted-adjoint test 证明该公式，再接入 Newton。禁止修改 J1 periodic helper；新增：

```cpp
build_open_x_reference_adjoint_cell_field(...)
```

---

## 6. 边界 pairing 功

定义：

$$
W_{\mathrm{boundary},E}
=\Delta t\frac{\Delta x}{2}
\left(
E_L^{\mathrm{pair}}\delta J_L
+E_R^{\mathrm{pair}}\delta J_R
\right)
$$

必须验证：

$$
\Delta t\langle E^{\mathrm{pair}},J^{\mathrm{open}}\rangle_f
=\Delta t\sum_iE_i^{\mathrm{ref}*}J_i^{\mathrm{cell}}\Delta x
+W_{\mathrm{boundary},E}
$$

该项是开放 current 与 cell force work 的边界配对余额，不是人为能量补丁。

---

## 7. J2 完整能量恒等式

定义净 reservoir kinetic 输入：

$$
K_{\mathrm{res,net}}
=K_{L,\mathrm{in}}+K_{R,\mathrm{in}}
-K_{L,\mathrm{out}}-K_{R,\mathrm{out}}
$$

J2 完整余额：

$$
R_{J2}
=\Delta K_{\mathrm{bulk}}
+\Delta U_E
-W_{\mathrm{electrode}}
-K_{\mathrm{res,net}}
+W_{\mathrm{boundary},E}
$$

应达到 solver tolerance。必须同时输出未加边界项的 residual，证明边界项非平凡时确实解释缺口。

当前零端点电势下 `W_electrode=0`，但字段和符号仍必须保留。

---

## 8. 文件级架构

### 8.1 `joint_phase_space_midpoint.h/.cpp`

新增：

```cpp
struct OpenXBoundaryLedger;
struct OpenJointPhaseSpaceFluxBundle;

static bool build_open_center_flux(...);
static bool evaluate_open_local_residual(...);
static bool build_open_x_reference_adjoint_cell_field(...);
static bool audit_open_x_pairing(...);
```

不得改变现有 periodic API 和 J0/J1 结果。

### 8.2 `vpfp_integrator.h/.cpp`

新增独立模式：

```text
JOINT_MIDPOINT_OPEN_BACKGROUND
```

不得让 `JOINT_MIDPOINT_ENERGY` 静默改变语义。新增独立推进函数：

```cpp
advance_joint_midpoint_open_background(...)
```

第一版可以复用 J1 Newton 框架，但必须调用 open residual、open adjoint 和 boundary ledger。

### 8.3 测试

新增：

```text
tests/joint_open_background_flux_test.cpp
tests/joint_open_background_energy_test.cpp
```

CMake target：

```text
joint_open_background_flux_test
joint_open_background_energy_test
```

---

## 9. 阶段顺序

```text
J2-0  冻结J1基线
J2-1  open flux与四类ledger
J2-2  open连续性/边界通量单元测试
J2-3  open reference G*与boundary pairing identity
J2-4  open local residual
J2-5  接入独立open-background integrator mode
J2-6  单rank物理case
J2-7  dt/dt2与10步
J2-8  1/2/5 rank
J2-9  全回归与完成门
```

任何阶段失败只允许修改当前白名单。

---

## 10. J2-0：冻结基线

记录 J1 A8 的结果路径、源码 hash 和当前 `pairing_tolerance`。运行 J0/J1 单 rank 和 1/2/5 rank，确认基线仍 PASS。不得修改代码。

阶段门：J1 全部保持 PASS。

---

## 11. J2-1：实现 open flux bundle

白名单：

```text
src/joint_phase_space_midpoint.h
src/joint_phase_space_midpoint.cpp
tests/joint_open_background_flux_test.cpp
CMakeLists.txt
docs/J2执行情况.md
```

实现第 2、3 节全部 flux 和 ledger。每个 physical face/velocity 槽必须输出 trace source：`RESERVOIR_IN`、`ABSORBING_IN` 或 `DOMAIN_OUT`。

阶段门仅编译，不接入 integrator。

---

## 12. J2-2：连续性与边界 ledger

测试 case：

```text
equilibrium-reservoir
left-inflow
right-inflow
left-outflow
right-outflow
absorbing-both
mpi-internal-face
```

每个 case 验证 number、charge current、kinetic、momentum ledger 同源；制造状态必须非对称，避免零值掩盖。

阶段门：cell continuity 与全局 number ledger 在 roundoff 内闭合，internal face 不重复计数。

---

## 13. J2-3：open adjoint 与边界功

先实现独立 helper，再在测试中构造非平凡 endpoint current。必须分别验证：

1. reference `G*` weighted identity；
2. `delta J_L/R` 非零；
3. 加上 `W_boundary,E` 后完整 face/cell identity 闭合；
4. 不加该项时 residual 明显非零；
5. 1/2/5 rank 结果只差舍入。

阶段门 PASS 前不得接入 Newton。

---

## 14. J2-4：open local residual

新增 `evaluate_open_local_residual()`，保持 J1 midpoint/u flux/正性语义，只替换 x flux topology。必须输出：

```text
phase_residual_linf
open_continuity_linf
boundary_number_residual
boundary_kinetic_residual
open_pairing_residual
```

测试 old=candidate、非零 candidate、signed accepted state 和左右 reservoir drift。

---

## 15. J2-5：接入独立 integrator mode

新增 mode 和 `advance_joint_midpoint_open_background()`。配置硬门：

```text
Beam OFF
Tail OFF
collision trivial
left/right boundary != PERIODIC
field DIRICHLET_PHI
```

不满足时明确 failure stage，不得 fallback 到 periodic J1 或 legacy。

Newton candidate 的 rho 仍由 stable mass increment 构造；open boundary 改变的总粒子数已经包含在 candidate mass 中，不得额外向 rho 添加 reservoir source。

---

## 16. J2-6：单 rank 物理验收

`joint_open_background_energy_test` 必须支持：

```text
equilibrium-reservoir
left-drift-inflow
right-drift-inflow
symmetric-outflow
smooth-perturbed-open
absorbing-open
```

正测试要求：finite、Gauss、phase、Poisson、pairing、boundary ledger、code-76 和完整 J2 energy 全部通过。

pairing relative 使用 J1 已验证门；总能量 relative 保持 `1e-8`。

---

## 17. J2-7：步长与多步

`smooth-perturbed-open` 运行：

```text
dt-scale=1.0, steps=1
dt-scale=0.5, steps=1
dt-scale=1.0, steps=10
```

累计余额必须包含 reservoir kinetic net 与 boundary pairing work。10 步每步都必须接受；100 步仍为非阻断可选测试。

---

## 18. J2-8：MPI

固定 `nx_global=20`，运行 1/2/5 rank。结构恒等式按 roundoff，nonlinear stopping 按 solver tolerance。必须输出实际 `pairing_tolerance`。

internal shared face 只计一次；global physical endpoint 只能由 rank 0 和 last rank 所有。

---

## 19. J2-9：完成定义

全部满足：

- J1 全回归不变；
- J2 flux/ledger/adjoint 单元测试通过；
- 单 rank 六 case 通过；
- dt、dt/2、10步通过；
- 1/2/5 rank 通过；
- 无 fallback、无能量补丁、无 periodic endpoint wrap。

才可写：

```text
J2 open-background joint midpoint PASS.
```

随后才允许 J3a Beam。

---

## 20. 编译与运行命令

所有命令由用户在集群根目录执行。

### 20.1 J2-1/J2-2

```bash
cmake --build build --target joint_open_background_flux_test -j4

RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/j2_flux_${RUN_ID}"
mkdir -p "$OUT"

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_open_background_flux_test \
  --case all --result "$OUT/flux_all.result" \
  > "$OUT/flux_all.out" 2> "$OUT/flux_all.err" || exit 1
```

### 20.2 J2-3

```bash
for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" --cpu-bind=cores \
    ./build/joint_open_background_flux_test \
    --case open-adjoint \
    --result "$OUT/open_adjoint_n${NP}.result" \
    > "$OUT/open_adjoint_n${NP}.out" \
    2> "$OUT/open_adjoint_n${NP}.err" || exit 1
done
```

### 20.3 J2-4/J2-5 构建

```bash
cmake --build build --target \
  joint_open_background_flux_test \
  joint_open_background_energy_test -j4
```

### 20.4 J2-6

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/j2_energy_${RUN_ID}"
mkdir -p "$OUT"

for CASE_NAME in equilibrium-reservoir left-drift-inflow right-drift-inflow symmetric-outflow smooth-perturbed-open absorbing-open; do
  yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
    ./build/joint_open_background_energy_test \
    --case "$CASE_NAME" --steps 1 --dt-scale 1.0 \
    --result "$OUT/${CASE_NAME}.result" \
    > "$OUT/${CASE_NAME}.out" 2> "$OUT/${CASE_NAME}.err" || exit 1
done
```

### 20.5 J2-7

```bash
for DT_SCALE in 1.0 0.5; do
  TAG=$(printf '%s' "$DT_SCALE" | tr '.' 'p')
  yhrun -N 1 -n 1 --cpu-bind=cores \
    ./build/joint_open_background_energy_test \
    --case smooth-perturbed-open --steps 1 --dt-scale "$DT_SCALE" \
    --result "$OUT/dt_${TAG}.result" \
    > "$OUT/dt_${TAG}.out" 2> "$OUT/dt_${TAG}.err" || exit 1
done

yhrun -N 1 -n 1 --cpu-bind=cores \
  ./build/joint_open_background_energy_test \
  --case smooth-perturbed-open --steps 10 --dt-scale 1.0 \
  --result "$OUT/cumulative_10steps.result" \
  > "$OUT/cumulative_10steps.out" \
  2> "$OUT/cumulative_10steps.err" || exit 1
```

### 20.6 J2-8

```bash
for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" --cpu-bind=cores \
    ./build/joint_open_background_energy_test \
    --case smooth-perturbed-open --steps 1 --dt-scale 1.0 \
    --result "$OUT/mpi_n${NP}.result" \
    > "$OUT/mpi_n${NP}.out" 2> "$OUT/mpi_n${NP}.err" || exit 1
done
```

---

## 21. 每阶段报告格式

每阶段在 `docs/J2执行情况.md` 追加：

````markdown
# 阶段：J2-?

## 当前阶段
`J2-?`

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
````
