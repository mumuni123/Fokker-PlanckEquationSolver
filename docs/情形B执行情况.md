# 阶段：B1 执行情况

状态：**已验收 —— PASS**（集群结果已读取并按 §15 阶段门判定通过；未修改任何 flux、field、residual、Newton、Poisson、energy gate、acceptance logic、dt 或阈值）

## 基线

```text
baseline_commit=42a15f60504ca3b55bc625a757bd307d8fb5214f
working_tree_clean=否（3 个未跟踪 docs 文件，均为方案文档）
集群与本地基线一致性（修改前核对）：
  src/vpfp_integrator.cpp                          md5 30ce9c0156502b016b2bc6d1d2cbecb4
  src/vpfp_integrator.h                            md5 214281da8d2f25a59fde7584b3bd29cd
  tests/joint_phase_space_midpoint_energy_test.cpp md5 e314746f0ae0254a4677f7af30b056c4
```

## 本地实际修改文件

- src/vpfp_integrator.h
- src/vpfp_integrator.cpp
- tests/joint_phase_space_midpoint_energy_test.cpp
- docs/情形B执行情况.md

## 修改内容（仅 B1.1--B1.5 诊断）

1. `src/vpfp_integrator.h`：`VpfpStepResult` 新增 7 个诊断字段并初始化 0.0：
   `joint_midpoint_pairing_face_left/right`、
   `joint_midpoint_force_current_first_cell/last_cell`、
   `joint_midpoint_naive_force_current_work`、
   `joint_midpoint_seam_predicted_residual`、
   `joint_midpoint_seam_prediction_error`。
2. `src/vpfp_integrator.cpp` `advance_joint_midpoint()`：
   - 入口处将 7 个新字段显式置零；
   - 在最终能量诊断处（`joint_midpoint_force_charge_residual` 计算之后、`failure_code=75` 早退之前）按 B1.2 用 production midpoint mass、production vH 独立构造 `force_current_cell`；按 B1.3 对 `final_pairing_face.front()/back()` 与 `force_current_cell.front()/back()` 做 4 元素 `MPI_Allreduce(SUM)`；另计算 naive gather 功诊断；按 B1.4 计算
     `R_seam_pred = 0.25*dt*grid_.dx*(E_left-E_right)*(J_first-J_last)`，
     `seam_prediction_error = joint_midpoint_force_charge_residual - R_seam_pred`。
   - 未修改任何生产算法对象。
3. `tests/joint_phase_space_midpoint_energy_test.cpp`：rank 0 result 输出按 B1.5 新增 7 行诊断输出。

## 明确未修改

- OpenElectrostaticSolver、build_potential_pairing_field、endpoint half-weight
- charge_current_face、x_flux_rate、u_flux_rate、midpoint mass、vH
- Newton、Poisson、energy gate、acceptance logic、dt、阈值
- Beam、Tail、collision、conversion、return、open boundary

## 编译命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
cmake --build build --target joint_phase_space_midpoint_energy_test -j4
```

## 测试命令（用户在集群执行）

```bash
yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_energy_test \
  --case smooth-perturbed-background \
  --result output/b1_seam_diagnostic.result \
  > output/b1_seam_diagnostic.out \
  2> output/b1_seam_diagnostic.err
```

注：本次运行使用文档 §14 原始路径 `output/b1_seam_diagnostic.*`（新建文件，未覆盖任何既有 output）。`.out` 为空文件；`.err` 仅含 yhrun 的退出码报告。

## 数值结果（来自 output/b1_seam_diagnostic.result，2026-08-22 只读审计取得）

```text
status=FAIL（旧算法预期失败：accepted=0, failure_code=75, joint_midpoint_energy_residual）
case=smooth-perturbed-background
finite=1 gauss_ok=1 converged=1 pairing_field_built=1 iterations=3
residual_linf=1.1480138226714858e-13
poisson_residual_linf=2.084279431871553e-16

W_u(joint_midpoint_u_face_work)=28.891863089349869
W_F(joint_midpoint_force_current_work)=28.891863089349886
W_J(joint_midpoint_charge_current_work)=22.224220061661196
R_uJ(current_pair_residual)=6.6676430276886727
R_FJ(force_charge_residual)=6.6676430276886904

E_left(joint_midpoint_pairing_face_left)=-51178802744.279716
E_right(joint_midpoint_pairing_face_right)=51178802743.112534
J_first(joint_midpoint_force_current_first_cell)=-130281340518476.31
J_last(joint_midpoint_force_current_last_cell)=130281340522669.23
naive_force_current_work=28.891863089349886（与生产 W_F 一致，因当前生产即 naive gather）

R_seam_pred(joint_midpoint_seam_predicted_residual)=6.6676430276867
seam_prediction_error=1.9904078385479806e-12
solver_exit_code=1（yhrun 报告 exit code 1，源于预期的 failure_code=75）
```

端点非对称性确认（B4 前置条件成立）：`E_left != E_right`，`J_first != J_last`。

## 阶段门（§15）

```text
S        = max(1, |R_FJ|, |R_seam_pred|) = 6.6676430276886904
tau_seam = 8192 * epsilon_mach * S = 1.212837201399214e-11
|seam_prediction_error| = 1.9904078385479806e-12 <= tau_seam
=> PASS
```

结论（§15 规定文本）：

```text
B root cause confirmed:
periodic x seam / nonperiodic pairing-face
face-to-cell adjoint mismatch
```

## 判定

- B1 PASS。R_FJ 与解析 seam 预测在 roundoff 内完全一致，情形 B 根因（periodic x current 算子与 naive face-to-cell pairing-field gather 不互为 weighted adjoint）已获数值确认。
- 下一步允许用户指定进入 B2；本报告不自行执行 B2。

---

# 阶段：B2 执行情况（根因确认）

状态：**PASS**（纯确认阶段；未修改任何源码、测试、CMake 或参数。集群测试由用户执行，本阶段直接读取已验收的 `output/b1_seam_diagnostic.result` 证据）

## 确认范围核对

1. **数值一致性**：
   `R_FJ = 6.6676430276886904` 与 `R_seam_pred = 6.6676430276867`
   之差 `seam_prediction_error = 1.9904078385479806e-12 <= tau_seam = 1.212837201399214e-11`
   （S = 6.6676430276886904，tau_seam = 8192 * eps_mach * S）。一致。
2. **根因确认**：B1 数值证据表明
   `W_F - W_J` 完全由
   periodic x seam current 算子 G（iface=0 与 iface=nx 同映射 last→first cell）
   与 endpoint half-weight 非周期 Poisson pairing face 的 naive gather
   `0.5*(E_i+E_{i+1})` 不互为 weighted adjoint 所解释；
   辅助证据：`joint_midpoint_naive_force_current_work = 28.891863089349886`
   与生产 `W_F` 完全一致（当前生产 W_F 即 naive gather），且
   `E_left != E_right`、`J_first != J_last`（seam 缺口未被对称性掩盖）。
   root_cause_confirmed=1。
3. **后续唯一允许变化**：构造 G 在 Poisson face quadrature 下的严格加权转置，
   即 u-force cell field 改用 `G* E_pair`（§18--§22 公式）。
4. **明确不修改**（后续阶段同样适用）：
   `x_flux_rate`、`charge_current_face`、Poisson、`build_potential_pairing_field`、
   endpoint half-weight、Newton、dt、vH、u_flux_rate、midpoint mass、
   energy gate、acceptance logic。

## 阶段门

PASS —— 下一步允许用户指定进入 B3（新增唯一 helper）；本报告不自行执行 B3。

---

# 阶段：B3 执行情况（helper 实现）

状态：**实现完成 —— 构建检查待集群执行**（本地 `-fsyntax-only` 通过；集群编译与后续数值验收由用户执行）

## 本地实际修改文件

- src/joint_phase_space_midpoint.h
- src/joint_phase_space_midpoint.cpp
- docs/情形B执行情况.md

## 修改内容

1. `JointPhaseSpaceMidpointOperator` 新增唯一静态 helper 声明与定义：
   `build_periodic_x_adjoint_cell_field(sg, pairing_face, mpi_rank, mpi_size, pairing_cell)`。
2. 输入检查（§24）：`sg.nx_local > 0`、`sg.nx_global >= 2`、
   `pairing_face.size() == sg.nx_local + 1`、
   `mpi_rank >= 0 && mpi_rank < mpi_size && mpi_size >= 1`、
   全部 `pairing_face` finite；任一失败即 `pairing_cell.clear(); return false;`，
   无任何 fallback。endpoint `E_left/E_right` 经 2 元素 `MPI_Allreduce(SUM)`
   获取并做 finite 检查。
3. cell field 构造严格按 §26/B3.3：
   - 全局首 cell：`0.5*pairing_face[ix+1] + 0.25*(E_left+E_right)`；
   - 全局末 cell：`0.5*pairing_face[ix] + 0.25*(E_left+E_right)`；
   - 内部 cell：`0.5*(pairing_face[ix]+pairing_face[ix+1])`；
   随后全部 finite-check，失败时 clear + false。
4. §27 要求的 `J1 TEST TOPOLOGY ONLY` 注释已原样写入 header 与 cpp 两处。

## 明确未修改

- vpfp_integrator、flux、charge_current_face、x_flux_rate、u_flux_rate、
  midpoint mass、vH、Poisson、OpenElectrostaticSolver、Newton、测试文件、CMake。

## 编译命令（由用户在集群执行）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
cmake --build build --target joint_phase_space_midpoint_unit_test -j4
```

## 数值结果

```text
build_exit_code=0
集群构建输出（用户提供）：
Linking CXX executable joint_phase_space_midpoint_unit_test
Built target joint_phase_space_midpoint_unit_test
```

## 阶段门

PASS —— B3 构建检查通过（§56.4）。B3 无独立数值阶段门；
下一步允许用户指定进入 B4（J0-E2 isolated periodic-seam weighted-adjoint 测试），
本报告不自行执行 B4。

---

# 阶段：B4 执行情况（J0-E2 isolated periodic-seam weighted adjoint）

状态：**已验收 —— PASS**（集群结果只读审计完成；`output/b4_periodic_seam_adjoint.result`，status=PASS，err 为空，退出码 0）

## 本地实际修改文件

- tests/joint_phase_space_midpoint_unit_test.cpp
- docs/情形B执行情况.md

## 修改内容

1. 在既有 `joint_phase_space_midpoint_unit_test` 的 `--case all` 中新增
   J0-E2 periodic-seam weighted adjoint（未新建 CMake target）：
   - B4.1：deterministic positive midpoint mass（单调 x tilt + u 非对称因子，
     无 Maxwellian 对称性），保证 `force_current_first_cell != last_cell`；
     人工 pairing face（含线性斜坡），保证
     `pairing_face[0] != pairing_face[N]`；
   - B4.2：用 production `build_hamiltonian_velocity()` 与测试 midpoint mass
     构造 `J_i^F`；
   - B4.3：调用 production `build_periodic_center_flux(...)` 取
     `bundle.charge_current_face`（未重建 J_charge）；
   - B4.4：调用 production `build_periodic_x_adjoint_cell_field(...)`
     得到 `e_adjoint_cell`；
   - B4.5：独立求和计算
     `W_F = dt*dx*sum(E_adj*J_F)` 与
     `W_J = dt*dx*[0.5*E_0*J_0 + sum E_f*J_f + 0.5*E_N*J_N]`
     （两个独立循环，无共享 helper）；
   - B4.6：额外构造 naive gather `E=0.5*(E_f+E_{f+1})` 计算 `W_F_naive`，
     验证 `W_F_naive - W_J` 明显非零且等于
     `R_seam_pred = 0.25*dt*dx*(E_left-E_right)*(J_first-J_last)`。
2. 输出与门控：
   `j0_e2_periodic_seam_weighted_adjoint_pass` 及 helper/asymmetry/
   weighted/prediction_nonzero/prediction_identity 五个子门、
   relative error、W_F/W_J/W_F_naive、mismatch、R_seam_pred、prediction_error、
   端点四值；E2 结果并入 `--case all` 的整体 status 与退出码；
   旧六门输出与语义未改动。

## 明确未修改

- production helper、vpfp_integrator、Poisson、flux、charge_current_face、
  Newton、物理参数、CMake、旧六门判定逻辑。

## 编译命令（由用户在集群执行）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
cmake --build build --target joint_phase_space_midpoint_unit_test -j4
```

## 测试命令（文档 §56.5，唯一输出目录）

```bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b4_${RUN_ID}"
mkdir -p "$OUT"

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/b4_periodic_seam_adjoint.result" \
  > "$OUT/b4_periodic_seam_adjoint.out" \
  2> "$OUT/b4_periodic_seam_adjoint.err" || exit 1

sed -n '1,420p' "$OUT/b4_periodic_seam_adjoint.result"
```

## 验收判据（收到 result 后）

```text
j0_a..j0_f 六个旧门全部 =1
j0_e2_periodic_seam_weighted_adjoint_pass=1
weighted_adjoint_relative_error <= 8192 * epsilon_machine
old_naive_mismatch 明显非零且 prediction_error <= 8192*eps*scale
first_cell_current != last_cell_current, pairing_face[0] != pairing_face[N]
```

## 数值结果（来自 output/b4_periodic_seam_adjoint.result，2026-08-22 只读审计取得）

```text
status=PASS
j0_a_cell_mass_pass=1
j0_b_u_flux_geometry_pass=1
j0_c_midpoint_consistency_pass=1
j0_d_u_work_force_adjoint_pass=1
j0_e_x_current_force_adjoint_pass=1
j0_f_poisson_pairing_pass=1
j0_all_six_pass=1

j0_e2_periodic_seam_weighted_adjoint_pass=1
j0_e2_helper_ok=1
j0_e2_endpoint_asymmetry_pass=1
j0_e2_weighted_adjoint_pass=1
j0_e2_prediction_nonzero_pass=1
j0_e2_prediction_identity_pass=1
weighted_adjoint_relative_error=4.163336342344337e-17   (<= 8192*eps = 1.819e-12)
W_F(G*)=-0.07460523480480262
W_J=-0.074605234804802661
W_F_naive=-0.074801909413251128
old_naive_mismatch=-0.00019667460844846607（明显非零，约为 |W_J| 的 0.26%）
seam_predicted_residual=-0.00019667460844851142
prediction_error=4.5346755864206223e-17                  (<= 8192*eps*scale)
pairing_face_left=200000000        pairing_face_right=249999999.99999997
force_current_first_cell=-76047515266758.438
force_current_last_cell=-91781483942639.359
```

## 阶段门

PASS —— §35 全部条件满足：

1. 旧六门全部保持 1，`j0_all_six_pass=1`，无回归；
2. `j0_e2_periodic_seam_weighted_adjoint_pass=1`，
   `|W_F-W_J|/max(1,|W_F|,|W_J|) = 4.16e-17 <= 8192*epsilon_machine`；
3. 旧 naive mismatch 明显非零，且与解析 seam 预测一致
   （prediction_error = 4.53e-17）；
4. 端点非对称性成立（E_left != E_right，J_first != J_last）。

结论：production helper `build_periodic_x_adjoint_cell_field`
满足 periodic seam 全局 weighted-adjoint 恒等式；
旧 naive gather 缺口被解析预测在 roundoff 内完全解释。
下一步允许用户指定进入 B5（Newton candidate 接线 G*）；本报告不自行执行 B5。

---

# 阶段：B5 执行情况（Newton candidate 使用 G*）

状态：**实现完成 —— 构建检查待集群执行**（本地 `-fsyntax-only` 通过；§56.6 构建检查由用户执行，无独立数值门，数值门在 B7/B8）

## 本地实际修改文件

- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

## 修改内容

1. `advance_joint_midpoint()` 内 candidate evaluate lambda 中：
   - `build_potential_pairing_field(...)` 调用保持不变；
   - 删除手写 naive 平均
     `e_pair_cell[ix] = 0.5*(pairing_face[ix]+pairing_face[ix+1])`；
   - 改为调用唯一 production helper
     `JointPhaseSpaceMidpointOperator::build_periodic_x_adjoint_cell_field(
         grid_, pairing_face, mpi_rank, mpi_size, e_pair_cell)`；
   - `adjoint_ok == false`（或 pairing face 尺寸非法）时
     `local_ok = false` → candidate evaluation 失败；无任何 fallback。
2. `e_local_out = e_pair_cell` 保持不变：preconditioner
   （`apply_local_block_diagonal_preconditioner`）与后续路径继续使用同一
   `G*E_pair` field。失败路径上 `accepted_e_local` 不会被消费
   （evaluate 失败即提前 return），行为安全。

## 明确未修改

- evaluate_local_residual、x_flux_rate、u_flux_rate、midpoint mass、vH、
  charge_current_face、OpenElectrostaticSolver、Poisson、Newton、GMRES、
  line search、energy gate、acceptance logic、dt、阈值。
- 最终能量诊断处的 naive gather（属 B6 范围，本阶段未触碰）。

## 编译命令（由用户在集群执行，文档 §56.6）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
```

## 数值结果

```text
build_exit_code=0
集群构建输出（用户提供）：
Linking CXX executable joint_phase_space_midpoint_unit_test
Built target joint_phase_space_midpoint_unit_test
Linking CXX executable joint_phase_space_midpoint_energy_test
Built target joint_phase_space_midpoint_energy_test
```

## 阶段门

PASS —— B5 构建检查通过（§56.6）。B5 无独立数值门（数值门在 B7/B8）。
下一步允许用户指定进入 B6（最终能量诊断接线同一 G* helper）；
本报告不自行执行 B6。

---

# 阶段：B6 执行情况（最终能量诊断使用 G*）

状态：**实现完成 —— 构建检查待集群执行**（本地 `-fsyntax-only` 通过；§56.6 构建检查由用户执行，数值门在 B7/B8）

## 本地实际修改文件

- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

## 修改内容

1. 最终能量诊断处：删除手写 naive 平均
   `final_pairing_cell[ix] = 0.5*(final_pairing_face[ix]+final_pairing_face[ix+1])`，
   改为调用唯一 production helper
   `build_periodic_x_adjoint_cell_field(grid_, final_pairing_face, mpi_rank,
   mpi_size, final_pairing_cell)`；
   `joint_midpoint_force_current_work` 继续由该
   `final_pairing_cell`（即实际生产 `G*E_pair`）构造，满足 §41。
2. §40 失败路径：`final_adjoint_ok == false` 时设置
   `failure_code = 71`、`failure_stage = "joint_midpoint_final_x_adjoint_field"`
   并返回；无任何 fallback。
3. B1 诊断字段 `joint_midpoint_naive_force_current_work`
   改为使用独立构造的 naive gather `0.5*(E_f+E_{f+1})`
   （仅诊断用途；不进入 residual、Newton、acceptance、energy gate，
   满足 §42）。

## 明确未修改

- candidate evaluate、evaluate_local_residual、x_flux_rate、u_flux_rate、
  charge_current_face、Poisson、OpenElectrostaticSolver、Newton、GMRES、
  line search、接受逻辑、energy gate、dt、阈值。

## 编译命令（由用户在集群执行，文档 §56.6）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
cmake --build build --target \
  joint_phase_space_midpoint_unit_test \
  joint_phase_space_midpoint_energy_test -j4
```

## 数值结果

```text
build_exit_code=0
集群构建输出（用户提供）：
Linking CXX executable joint_phase_space_midpoint_unit_test
Built target joint_phase_space_midpoint_unit_test
Linking CXX executable joint_phase_space_midpoint_energy_test
Built target joint_phase_space_midpoint_energy_test
```

## 阶段门

PASS —— B6 构建检查通过（§56.6）。B6 无独立数值门（数值门在 B7/B8）。
下一步允许用户指定进入 B7（全部 J0 回归）；
本报告不自行执行 B7。

---

# 阶段：B7 执行情况（全部 J0 回归）

状态：**已验收 —— PASS**（集群结果只读审计完成；`output/b7_j0_all.result`，status=PASS，err 为空，退出码 0。本阶段未修改任何源码、测试、CMake 或参数）

## 集群命令（文档 §56.7，唯一输出目录）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b7_${RUN_ID}"
mkdir -p "$OUT"

yhrun -N 1 -n 1 --cpu-bind=cores stdbuf -oL -eL \
  ./build/joint_phase_space_midpoint_unit_test \
  --case all \
  --result "$OUT/b7_j0_all.result" \
  > "$OUT/b7_j0_all.out" \
  2> "$OUT/b7_j0_all.err" || exit 1

sed -n '1,520p' "$OUT/b7_j0_all.result"
```

## PASS 判据（收到 result 后验收）

以下全部为 1：

```text
j0_a_cell_mass_pass
j0_b_u_flux_geometry_pass
j0_c_midpoint_consistency_pass
j0_d_u_work_force_adjoint_pass
j0_e_x_current_force_adjoint_pass
j0_f_poisson_pairing_pass
j0_e2_periodic_seam_weighted_adjoint_pass
```

任何旧门回归或 J0-E2 FAIL 均禁止进入 B8。

## 数值结果（来自 output/b7_j0_all.result，2026-08-22 只读审计取得）

```text
status=PASS
j0_a_cell_mass_pass=1
j0_b_u_flux_geometry_pass=1
j0_c_midpoint_consistency_pass=1
j0_d_u_work_force_adjoint_pass=1
j0_e_x_current_force_adjoint_pass=1
j0_f_poisson_pairing_pass=1
j0_all_six_pass=1
j0_e2_periodic_seam_weighted_adjoint_pass=1
j0_e2_helper_ok=1
j0_e2_endpoint_asymmetry_pass=1
j0_e2_weighted_adjoint_pass=1
j0_e2_prediction_nonzero_pass=1
j0_e2_prediction_identity_pass=1
j0_e2_weighted_adjoint_relative_error=4.163336342344337e-17
j0_e2_old_naive_mismatch=-0.00019667460844846607
j0_e2_seam_predicted_residual=-0.00019667460844851142
j0_e2_prediction_error=4.5346755864206223e-17
```

## 阶段门

PASS —— §43 验收满足：六个旧 J0 门全部保持 1（无回归），且

```text
j0_e2_periodic_seam_weighted_adjoint_pass=1
```

B5/B6 接线后 J0 全部回归通过。下一步允许用户指定进入 B8
（F10 两个单 rank case 回归）；本报告不自行执行 B8。

---

# 阶段：B8 执行情况（F10 两个单 rank case 回归）

状态：**已验收 —— 按 §46 命中情形 2（情形 A），B 阶段停止**
（集群结果只读审计完成：`output/b8_smooth-background_n1.result` 与
`output/b8_smooth-perturbed-background_n1.result`。本阶段未修改任何源码、测试、CMake 或参数）

## 数值结果（2026-08-22 只读审计取得）

### smooth-background（平衡态，§45.1）

```text
status=PASS
accepted=1
finite=1
gauss_ok=1
converged=1
failure_code=0
err 为空（退出码 0）
```

§45.1 全部满足，无回归。

### smooth-perturbed-background（非平凡场，§45.2--§45.4）

```text
W_u(joint_midpoint_u_face_work)      = 17.779409055027578
W_F(joint_midpoint_force_current_work) = 17.779409055027571
W_J(joint_midpoint_charge_current_work)= 17.779409055027124

W_u-W_F = 7.1e-15                     （roundoff，§45.2 通过）
R_FJ 相对误差 = 4.4764192352886312e-13 / 17.779... = 2.52e-14 <= 1e-12 （§45.3 通过）
R_uJ(current_pair_residual) = 4.5474735088646412e-13 （roundoff，§45.4 通过）

诊断自洽性：
joint_midpoint_naive_force_current_work = 22.224575607258192
joint_midpoint_seam_predicted_residual  = 4.4451665522306207
naive mismatch 与 seam 预测一致；
生产 W_F(G*)-W_J 已闭合到 roundoff。

后续门（§45 最后判断）：
delta_k_u                = 17.779409055028388
field_energy_change      = -17.779411022900604
electrode_work           = 0
poisson_potential_charge_work = -17.779411022868356
R_PJ(poisson_transport_residual) = -1.9678734801686915e-06   ← 仍大
energy_residual          = -1.9678730254213406e-06
combined_energy_residual / scale ≈ 1.107e-07 > 1e-8 能量门
accepted=0
failure_stage=joint_midpoint_energy_residual
failure_code=75
```

## 阶段判定（§46 决策树）

恒等式链验收：

1. `W_u-W_F` 为 roundoff —— 通过；
2. `abs(W_F-W_J)/max(1,|W_F|,|W_J|) = 2.52e-14 <= 1e-12` —— 通过；
3. `R_uJ` 为 roundoff —— 通过；
4. `R_PJ = -1.9679e-06` 仍大（相对尺度约 1.1e-07），
   且能量门触发 `failure_code=75`、`accepted=0`。

=> **§46 情形 2：R_uJ 已小而 R_PJ 仍大 —— 新失败为第 11 节情形 A。**

处置（严格按 §46）：

- 立即停止 B 阶段；
- 不得继续修改 `G*`；
- 不得修改 Poisson；
- 必须另行执行情形 A 的
  `charge continuity ↔ Poisson pairing` 审计；
- 不进入 J1 MPI gate，不进入 F11/J2/J3。

备注：旧 F10 的 `R_PJ = -5.295e-06`，修复后为 `-1.968e-06`
——与方案预警一致：它不会因 B 修复自动消失，
需按其自身的 Poisson work tolerance 重新判断。

## 集群命令（文档 §56.8，唯一输出目录，两个 case 分别保存）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
RUN_ID=$(date +%Y%m%d_%H%M%S)
OUT="./output/f10_case_b_b8_${RUN_ID}"
mkdir -p "$OUT"

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

## 验收顺序（§45--§46，已按此顺序执行）

1. `smooth-background`：status=PASS、accepted=1、finite=1、gauss_ok=1、
   converged=1、failure_code=0 —— 通过；
2. `smooth-perturbed-background`：W_u-W_F 为 roundoff —— 通过；
3. abs(W_F-W_J)/max(1,abs(W_F),abs(W_J)) = 2.52e-14 <= 1e-12 —— 通过；
4. R_uJ 为 roundoff —— 通过；
5. R_PJ = -1.9679e-06 仍大，能量门触发 failure_code=75、accepted=0
   —— 命中 §46 情形 2（情形 A）。

实际运行使用 `output/b8_*` 路径保存两份 result/out/err（未覆盖任何旧 output）。
数值结果见上方"数值结果（2026-08-22 只读审计取得）"小节。





