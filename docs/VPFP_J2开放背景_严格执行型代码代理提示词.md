# VPFP J2 开放背景：严格执行型代码代理提示词

一次只复制一个阶段代码块。不得同时执行多个阶段。

## J2-0

```text
你是严格执行型代码代理。完整阅读并严格执行：

docs/VPFP_J2开放背景联合中点离散与验收实施方案.md
docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
docs/VPFP能量问题根因链_已确认结论.md

本次只执行J2-0。

规则：
1. 绝对遵守J2方案第0节。
2. 不修改源码，只冻结J1基线和结果路径。
3. 本地只更新docs/J2执行情况.md；集群测试由用户执行。
4. 缺少集群结果时写“待集群执行”。
5. 不得进入J2-1。

执行：
1. 记录J1 A8、A6、A7结果路径和当前pairing_tolerance。
2. 列出J0/J1基线编译与测试命令交给用户。
3. 收到结果后确认全部PASS；否则停止。
4. 按方案第21节写报告。
```

## J2-1

```text
你是严格执行型代码代理。完整阅读：
docs/VPFP_J2开放背景联合中点离散与验收实施方案.md

前置：J2-0 stage gate=PASS。
本次只执行J2-1，不得进入J2-2。

白名单：
- src/joint_phase_space_midpoint.h
- src/joint_phase_space_midpoint.cpp
- tests/joint_open_background_flux_test.cpp
- CMakeLists.txt
- docs/J2执行情况.md

强制要求：
1. 实现方案第2、3节的open face flux和四类ledger。
2. 内部face保持midpoint centered trace；物理face按vH符号选择reservoir incoming或domain outgoing。
3. vH==0时flux严格为0。
4. 所有ledger在同一个face/velocity循环内累加。
5. 不修改periodic API，不接入integrator。
6. 不得从density差反推current/ledger。
7. 所有修改本地完成；用户运行集群编译。

执行：
1. 新增OpenXBoundaryLedger/OpenJointPhaseSpaceFluxBundle。
2. 新增build_open_center_flux生产函数。
3. 增加finite、尺寸、边界类型检查；失败不得fallback。
4. 新建测试骨架和CMake target，但本阶段只要求编译。
5. 本地静态检查和git diff --check。
6. 报告“待集群执行”并给出20.1编译命令。
```

## J2-2

```text
你是严格执行型代码代理。完整阅读J2方案第2、3、12、20.1节。

前置：J2-1集群编译PASS。
本次只执行J2-2，不修改integrator。

白名单：
- tests/joint_open_background_flux_test.cpp
- 必要时joint_phase_space_midpoint.h/.cpp，仅修测试发现的open flux错误
- docs/J2执行情况.md

执行：
1. 实现equilibrium-reservoir、left/right inflow、left/right outflow、absorbing-both、mpi-internal-face。
2. 每个case必须非对称并直接调用build_open_center_flux。
3. 验证cell continuity、全局number、charge current、kinetic和momentum ledger。
4. 验证internal shared face单owner，physical endpoint仅rank0/last拥有。
5. 不得为了PASS修改阈值或测试期待值。
6. 给用户20.1命令；结果缺失时待集群执行。
7. 任一case失败即停止，不进入J2-3。
```

## J2-3

```text
你是严格执行型代码代理。完整阅读J2方案第4--6、13、20.2节。

前置：J2-2 PASS。
本次只执行J2-3，不接入Newton。

白名单：
- src/joint_phase_space_midpoint.h
- src/joint_phase_space_midpoint.cpp
- tests/joint_open_background_flux_test.cpp
- docs/J2执行情况.md

执行：
1. 新增build_open_x_reference_adjoint_cell_field，逐字实现第5节首/中/末cell公式。
2. 保持build_periodic_x_adjoint_cell_field完全不变。
3. 从实际open endpoint current减去reference half-cell current，构造delta_J_L/R。
4. 构造W_boundary_E并验证第6节完整identity。
5. 测试必须让delta_J_L/R和未修正residual明显非零。
6. 运行1/2/5 rank；结构误差按roundoff门。
7. 任一门失败停止，不得进入J2-4。
```

## J2-4

```text
你是严格执行型代码代理。完整阅读J2方案第14节。

前置：J2-3 PASS。
本次只执行J2-4，不新增integrator mode。

白名单：
- src/joint_phase_space_midpoint.h/.cpp
- tests/joint_open_background_flux_test.cpp
- docs/J2执行情况.md

执行：
1. 新增evaluate_open_local_residual。
2. midpoint mass、u flux、vH、signed residual和code-76契约必须与J1一致。
3. 只替换x topology为open flux。
4. 输出phase/open continuity/boundary number/boundary kinetic/open pairing residual。
5. 测试old=candidate、非零candidate、signed accepted、左右reservoir drift。
6. 不得fallback到periodic evaluator。
7. 编译和单元结果PASS后才允许J2-5。
```

## J2-5

```text
你是严格执行型代码代理。完整阅读J2方案第7、8.2、15节。

前置：J2-4 PASS。
本次只执行J2-5，数值物理门留给J2-6。

白名单：
- src/vpfp_integrator.h
- src/vpfp_integrator.cpp
- tests/joint_open_background_energy_test.cpp
- CMakeLists.txt
- docs/J2执行情况.md

执行：
1. 新增JOINT_MIDPOINT_OPEN_BACKGROUND，禁止改变JOINT_MIDPOINT_ENERGY。
2. 新增advance_joint_midpoint_open_background。
3. 配置硬拒绝Beam/Tail/collision/periodic background/non-Dirichlet field。
4. 复用J1 Newton框架、stable rho、stable Poisson identity和三重门。
5. 调用open residual/open adjoint/boundary ledger。
6. 按第7节构造完整J2 energy residual。
7. 不得fallback到periodic或legacy。
8. 新建energy test target和CLI骨架。
9. 本阶段只验收编译；待集群执行。
```

## J2-6

```text
你是严格执行型代码代理。完整阅读J2方案第16、20.4节。

前置：J2-5集群编译PASS。
本次只执行J2-6。

白名单：
- tests/joint_open_background_energy_test.cpp
- 必要时当前阶段直接暴露的J2生产文件
- docs/J2执行情况.md

执行：
1. 实现六个指定case，不得减少。
2. 每个case直接调用JOINT_MIDPOINT_OPEN_BACKGROUND生产推进。
3. 输出finite/Gauss/三重门/code76、四类number/K/P ledger、W_boundary_E、W_electrode、完整R_J2。
4. equilibrium必须非平凡检查incoming/outgoing相消，不能全零蒙混。
5. drift case必须验证方向和符号。
6. 总能量relative<=1e-8，pairing使用结果中实际阈值。
7. 给用户20.4命令；任一case失败停止。
```

## J2-7

```text
你是严格执行型代码代理。完整阅读J2方案第17、20.5节。

前置：J2-6 PASS。
本次只执行J2-7。

规则：
1. 测试程序支持--dt-scale和--steps，默认1。
2. 运行dt、dt/2和10步；100步不是阻断门。
3. 每步真实连续推进接受态。
4. 累计量必须包含reservoir kinetic net、boundary pairing和electrode work。
5. 任一步失败仍写result并记录failed_step_index。
6. 不得清零累计量或只记录最后一步。
7. 修改本地，用户运行20.5命令。
8. 三档全部PASS才允许J2-8。
```

## J2-8

```text
你是严格执行型验收代理。完整阅读J2方案第18、20.6节。

前置：J2-7 PASS。
本次只执行J2-8，原则上不修改生产代码。

执行：
1. 固定nx_global=20和完全相同初态，给用户1/2/5 rank命令。
2. 结构恒等式按roundoff；nonlinear停止量按实际pairing_tolerance。
3. 迭代次数允许不同，不得要求R_J2跨rank逐位一致。
4. 检查internal face单owner、physical endpoint唯一owner、ledger不随rank倍增。
5. 缺少任一档结果则待集群执行。
6. ownership真失败才允许修改J2 face owner代码。
7. PASS只允许J2-9。
```

## J2-9

```text
你是严格执行型最终验收代理。完整阅读J2方案第19节。

前置：J2-0至J2-8全部stage gate=PASS。
本次只执行J2-9，不新增算法。

执行：
1. 重跑J1全回归，必须无回归。
2. 汇总J2 flux/ledger/adjoint/single-rank/dt/10step/MPI全部结果。
3. 检查无periodic endpoint wrap、无fallback、无能量补丁。
4. 只有全部满足才在docs/J2执行情况.md写：J2 open-background joint midpoint PASS。
5. PASS后只说明允许J3a Beam，不自行进入。
```
