# B0

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守该文档第 1 节“不得改变的物理模型”、第 8 节 B0 的修改范围和禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 严格按 B0 -> B1 -> B2 -> B3 -> B4 -> B5 -> B6 -> B7 -> B8 的顺序执行。当前阶段未 PASS，立即停止，不得进入下一阶段。
3. 只修改当前阶段白名单中的文件；不得顺手重构、优化或清理无关代码。
4. 不得通过能量投影、电流/电场缩放、裁剪分布函数、放宽阈值或伪造 PASS 来处理失败。
5. 远程操作必须严格遵守 docs/ssh远程连接协议.md。
6. 测试必须直接调用生产实现；不得在测试中复写公式后宣称验证生产代码。
7. 不得删除、覆盖或清空旧 output 结果；每次测试使用唯一输出目录。
8. 如果文档与源码不一致、目标不存在、命令失败或结果不清楚，停止并报告事实，不要猜测或自行补方案。

本次只执行阶段：B0。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- docs/情形B执行情况.md

本阶段实施内容：
1. 不修改任何源码、测试、CMake、参数或已有 output。
2. 在项目根目录执行：
   git status --short
   git rev-parse HEAD
   git diff --check
3. 记录 baseline_commit、working_tree_clean 和 diff check。
4. 工作树非空不代表失败；不得执行 git reset、git clean、删除或覆盖文件。

执行流程：
1. 先阅读 B0 原文、SSH 协议和 docs/情形B执行情况.md。
2. 列出当前阶段、允许修改文件、计划执行命令和 PASS 判据。
3. 执行上述只读 Git 命令。
4. 读取本地 Git 输出，写入 docs/情形B执行情况.md。
5. 最终严格使用以下格式添加或更新本阶段输出：

# 阶段：B0 执行报告

## 当前阶段

B0

## 实际修改文件

- docs/情形B执行情况.md

## 修改内容

1. 未修改源码或测试。
2. 记录 Git 基线。

## 编译命令
~~~
无。
~~~

## 测试命令
~~~
git status --short
git rev-parse HEAD
git diff --check
~~~

数值结果：
~~~
baseline_commit=
working_tree_clean=
diff_check=
~~~

阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B1，不得自行进入。
- FAIL：列出失败命令和原因；停止，不得执行 B1。
~~~

# B1

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守该文档第 1 节物理模型、第 9--15 节 B1 和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B1。B1 未 PASS，立即停止，不得进入 B2。
3. 只修改当前阶段白名单中的文件；不得顺手重构、优化或清理无关代码。
4. 不得通过能量投影、电流/电场缩放、裁剪分布函数、放宽阈值或伪造 PASS 处理失败。
5. 远程操作必须严格遵守 docs/ssh远程连接协议.md。
6. 不得删除、覆盖或清空旧 output；必须使用唯一输出目录。
7. 文档与源码不一致或结果不清楚时停止，不得猜测。

本次只执行阶段：B1。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。


本阶段允许修改文件：
- src/vpfp_integrator.h
- src/vpfp_integrator.cpp
- tests/joint_phase_space_midpoint_energy_test.cpp
- docs/情形B执行情况.md

本阶段实施内容：
1. 只增加 B1.1--B1.5 所需诊断和结果输出：
   joint_midpoint_pairing_face_left/right
   joint_midpoint_force_current_first_cell/last_cell
   joint_midpoint_naive_force_current_work
   joint_midpoint_seam_predicted_residual
   joint_midpoint_seam_prediction_error
2. 从最终 accepted candidate 的 production midpoint mass、production vH、final_pairing_face 和 force current 构造：
   R_seam_pred = 0.25*dt*dx*(E_left-E_right)*(J_first-J_last)
3. 只增加诊断。不得修改 flux、field、residual、Newton、Poisson、energy gate、acceptance logic、dt 或测试阈值。

PASS 判据：
abs(seam_prediction_error) <= tau_seam

执行流程：
1. 先阅读 B1 原文、B0 报告和相关源码。
2. 列出计划修改的函数、字段、命令和 PASS 判据。
3. 再在本地修改代码，并向用户提供集群编译/运行命令。
4. 未收到用户提供的 result/out/err 前，在 docs/情形B执行情况.md 写“待集群执行”；收到结果后才验收。

最终严格使用以下格式：

# 阶段：B1 执行报告

## 当前阶段

B1

## 实际修改文件

- <实际文件>

## 修改内容

1. <新增诊断>
2. <新增输出>

数值结果：

E_left=
E_right=
J_first=
J_last=
R_FJ=
R_seam_pred=
seam_prediction_error=
tau_seam=
solver_exit_code=

阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B2，不得自行进入。
- FAIL：列出首个失败数值、对应源码位置和 B1 允许范围；停止，不得执行 B2--B8。
~~~

# B2

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守该文档第 16--22 节 B2 的根因确认范围。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B2。B2 未 PASS，停止，不得进入 B3。
3. 远程操作必须严格遵守 docs/ssh远程连接协议.md。
4. 不得删除、覆盖或清空旧 output。
5. 结果不清楚时停止，不得猜测。

本次只执行阶段：B2。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- docs/情形B执行情况.md

本阶段实施内容：
1. 不修改任何源码、测试、CMake 或参数。
2. 读取用户提供的 B1 result/out/err 与 docs/情形B执行情况.md。
3. 只确认：
   - R_FJ 与 R_seam_pred 在 tau_seam 内一致；
   - 根因是 periodic x seam current G 与 endpoint-half-weight pairing face 的 naive gather 不互为 weighted adjoint；
   - 后续唯一允许变化是构造 G*E_pair；
   - 不修改 x_flux_rate、charge_current_face、Poisson、Newton、dt、vH、u_flux_rate、midpoint mass。

最终严格使用以下格式：

# 阶段：B2 执行报告

## 当前阶段

B2

## 实际修改文件

- docs/情形B执行情况.md

## 修改内容

1. 未修改源码。
2. 核对 B1 根因数值。


数值结果：

R_FJ=
R_seam_pred=
seam_prediction_error=
tau_seam=
root_cause_confirmed=


阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B3，不得自行进入。
- FAIL：说明哪项 B1 证据不足；停止，不得实施 helper。
~~~

# B3

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 23--27 节 B3、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B3。B3 未 PASS，停止，不得进入 B4。
3. 只修改当前阶段白名单中的文件。
4. 远程操作必须严格遵守 docs/ssh远程连接协议.md。
5. 不得删除、覆盖或清空旧 output。

本次只执行阶段：B3。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- src/joint_phase_space_midpoint.h
- src/joint_phase_space_midpoint.cpp
- docs/情形B执行情况.md

本阶段实施内容：
1. 新增唯一静态 helper：
   build_periodic_x_adjoint_cell_field(
       const SpatialGrid& sg,
       const std::vector<double>& pairing_face,
       int mpi_rank,
       int mpi_size,
       std::vector<double>& pairing_cell)
2. 严格按 B3.3 构造：
   first cell = 0.5*pairing_face[ix+1] + 0.25*(E_left+E_right)
   last cell  = 0.5*pairing_face[ix]   + 0.25*(E_left+E_right)
   interior   = 0.5*(pairing_face[ix]+pairing_face[ix+1])
3. 检查输入尺寸、nx_global、MPI rank、MPI size 和 finite；用 MPI_Allreduce 获取 E_left/E_right；失败时 pairing_cell.clear(); return false;；不得 fallback。
4. 写入 B3.4 要求的 J1 TEST TOPOLOGY ONLY 注释。
5. 禁止修改 vpfp_integrator、flux、charge_current_face、Poisson、Newton、测试和 CMake。

最终严格使用以下格式：

# 阶段：B3 执行报告

## 当前阶段

B3

## 实际修改文件

- src/joint_phase_space_midpoint.h
- src/joint_phase_space_midpoint.cpp
- docs/情形B执行情况.md

## 修改内容

1. <helper 声明和定义>
2. <输入检查和注释>



数值结果：

build_exit_code=


阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B4，不得自行进入。
- FAIL：列出编译错误和允许修改范围；停止。
~~~

# B4

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 28--35 节 B4、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B4。B4 未 PASS，停止，不得进入 B5。
3. 测试必须调用 production helper 与 production x bundle；不得复写替代 current。
4. 不得新增 CMake target。
5. 远程操作必须严格遵守 SSH 协议；不得覆盖旧 output。

本次只执行阶段：B4。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- tests/joint_phase_space_midpoint_unit_test.cpp
- docs/情形B执行情况.md

本阶段实施内容：
1. 在既有 joint_phase_space_midpoint_unit_test 的 --case all 中新增 J0-E2 periodic-seam weighted adjoint。
2. 使用 production build_periodic_x_adjoint_cell_field 和 build_periodic_center_flux。
3. 构造 deterministic positive M_mid，并保证 J_first != J_last、pairing_face[0] != pairing_face[N]。
4. 验证 weighted-adjoint relative error <= 8192*epsilon_machine。
5. 验证 old naive gather mismatch 非零且与 B1 seam prediction 一致。
6. 输出 j0_e2_periodic_seam_weighted_adjoint_pass。
7. 禁止修改 production helper、vpfp_integrator、Poisson、flux、Newton 或物理参数。

最终严格使用以下格式：

# 阶段：B4 执行报告

## 当前阶段

B4

## 实际修改文件

- tests/joint_phase_space_midpoint_unit_test.cpp
- docs/情形B执行情况.md

## 修改内容

1. <J0-E2 实现>
2. <生产 helper/bundle 调用>

数值结果：
j0_a_cell_mass_pass=
j0_b_u_flux_geometry_pass=
j0_c_midpoint_consistency_pass=
j0_d_u_work_force_adjoint_pass=
j0_e_x_current_force_adjoint_pass=
j0_f_poisson_pairing_pass=
j0_e2_periodic_seam_weighted_adjoint_pass=
weighted_adjoint_relative_error=
naive_mismatch=

阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B5，不得自行进入。
- FAIL：列出首个失败数值和 B4 允许修改范围；停止。
~~~

# B5

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 36--39 节 B5、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B5。B5 未 PASS，停止，不得进入 B6。
3. 只修改当前阶段白名单中的文件。
4. 不得修改 evaluate_local_residual、x_flux_rate、u_flux_rate、midpoint mass、vH、charge_current_face、OpenElectrostaticSolver、Poisson、Newton、GMRES 或 line search。
5. 远程操作必须遵守 SSH 协议；不得覆盖旧 output。

本次只执行阶段：B5。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

本阶段实施内容：
1. 在 candidate evaluate 中保留 build_potential_pairing_field 不变。
2. 删除手写 naive average：
   e_pair_cell[ix] = 0.5*(pairing_face[ix]+pairing_face[ix+1])
3. 改为调用：
   JointPhaseSpaceMidpointOperator::build_periodic_x_adjoint_cell_field(
       grid_, pairing_face, mpi_rank, mpi_size, e_pair_cell)
4. adjoint_ok=false 时 candidate evaluation 失败；不得 fallback。
5. 保留 e_local_out=e_pair_cell，使现有 preconditioner 使用同一 field。

最终严格使用以下格式：

# 阶段：B5 执行报告

## 当前阶段

B5

## 实际修改文件

- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

## 修改内容

1. <candidate field 改为 G*E_pair>
2. <无 fallback 路径>

数值结果：

build_exit_code=

阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B6，不得自行进入。
- FAIL：列出错误和 B5 允许修改范围；停止。
~~~

# B6

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 40--42 节 B6、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B6。B6 未 PASS，停止，不得进入 B7。
3. 只修改当前阶段白名单中的文件。
4. 不得修改 candidate evaluate、evaluate_local_residual、x_flux_rate、u_flux_rate、charge_current_face、Poisson、Newton、接受逻辑或阈值。
5. 远程操作必须遵守 SSH 协议；不得覆盖旧 output。

本次只执行阶段：B6。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

本阶段实施内容：
1. 删除最终诊断的 naive average：
   final_pairing_cell[ix] = 0.5*(final_pairing_face[ix]+final_pairing_face[ix+1])
2. 改为调用 build_periodic_x_adjoint_cell_field。
3. final_adjoint_ok=false 时设置 failure_code=71、
   failure_stage=joint_midpoint_final_x_adjoint_field 并返回；不得 fallback。
4. joint_midpoint_force_current_work 必须使用 final_pairing_cell。
5. naive work 只可保留为诊断，不得进入 residual、Newton、acceptance 或 energy gate。


最终严格使用以下格式：

# 阶段：B6 执行报告

## 当前阶段

B6

## 实际修改文件

- src/vpfp_integrator.cpp
- docs/情形B执行情况.md

## 修改内容

1. <final diagnostic field 改为 G*E_pair>
2. <failure_code=71 路径>

数值结果：

build_exit_code=


阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B7，不得自行进入。
- FAIL：列出错误和 B6 允许修改范围；停止。
~~~

# B7

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 43 节 B7、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B7。B7 未 PASS，停止，不得进入 B8。
3. 本阶段禁止修改任何源码、测试、CMake 或参数。
4. 远程操作必须遵守 SSH 协议；不得覆盖旧 output。

本次只执行阶段：B7。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- docs/情形B执行情况.md

PASS 判据：以下全部为 1：
j0_a_cell_mass_pass
j0_b_u_flux_geometry_pass
j0_c_midpoint_consistency_pass
j0_d_u_work_force_adjoint_pass
j0_e_x_current_force_adjoint_pass
j0_f_poisson_pairing_pass
j0_e2_periodic_seam_weighted_adjoint_pass

最终严格使用以下格式：

# 阶段：B7 执行报告

## 当前阶段

B7

## 实际修改文件

- docs/情形B执行情况.md

## 修改内容

1. 未修改源码。
2. 向用户提供 J0 回归命令；收到用户提供的 result/out/err 后再验收。

数值结果：

j0_a_cell_mass_pass=
j0_b_u_flux_geometry_pass=
j0_c_midpoint_consistency_pass=
j0_d_u_work_force_adjoint_pass=
j0_e_x_current_force_adjoint_pass=
j0_f_poisson_pairing_pass=
j0_e2_periodic_seam_weighted_adjoint_pass=

阶段门：
PASS 或 FAIL

- PASS：仅说明下一步允许用户指定 B8，不得自行进入。
- FAIL：列出首个失败字段；停止。
~~~

# B8

~~~text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形B根因定位与严格修复实施方案.md

最高规则：
0. 绝对遵守第 44--46 节 B8、物理模型和第 51 节禁止项。
1. 本文档是本次唯一执行规格；不得自行修改算法、物理模型、边界条件、Poisson、Beam、Tail、碰撞模型或阈值。
2. 只执行 B8。不得自动进入 F11、J1 MPI、J2 或 J3。
3. 本阶段禁止修改任何源码、测试、CMake 或参数。
4. 远程操作必须遵守 SSH 协议；不得覆盖旧 output。

本次只执行阶段：B8。

补充约束：所有文件修改必须在本地工作区进行。集群测试由用户自行执行；代码代理不得通过 SSH 修改远程文件、编译或提交/运行测试作业。若尚未获得用户提供的集群 result/out/err，只能报告“待集群执行”，不得推测 PASS 或 FAIL。

本阶段允许修改文件：
- docs/情形B执行情况.md

验收顺序：
1. 平衡态必须 status=PASS、accepted=1、finite=1、gauss_ok=1、converged=1、failure_code=0。
2. 非平凡场先检查 W_u-W_F 是否为 roundoff。
3. 再检查 abs(W_F-W_J)/max(1,abs(W_F),abs(W_J)) <= 1e-12。
4. 再检查 R_uJ。
5. 最后判断 R_PJ、R_E、accepted 和 failure_code。

最终严格使用以下格式：

# 阶段：B8 执行报告

## 当前阶段

B8

## 实际修改文件

- docs/情形B执行情况.md

## 修改内容

1. 未修改源码。
2. 向用户提供两个单 rank J1 命令；收到用户提供的 result/out/err 后再重新分流。


数值结果：

smooth_background_status=
smooth_perturbed_status=
W_u=
W_F=
W_J=
R_uJ=
R_PJ=
R_E=
accepted=
failure_code=
failure_stage=

阶段门：
PASS 或 FAIL

- PASS：若全部 identity PASS 且 accepted=1/failure_code=0，仅说明下一步允许用户指定 J1 MPI gate。
- FAIL：若 R_uJ 小而 R_PJ 大，标记情形 A 并停止；若 identity PASS 但 code=73/74/76，仅说明允许用户另行指定 F11；其他失败记录首个失败项并停止。
~~~
