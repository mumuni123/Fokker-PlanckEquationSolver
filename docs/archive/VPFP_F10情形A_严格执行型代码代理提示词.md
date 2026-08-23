# VPFP F10 情形 A：严格执行型代码代理提示词

一次只复制当前阶段对应的一个代码块，不得同时提交多个阶段。

## A0

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

同时阅读：

docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
docs/情形B执行情况.md

最高规则：
0. 绝对遵守主方案第0节和第16节。
1. 本次只执行A0，不得进入A1。
2. 本阶段只记录基线，不得修改源码。
3. 不得改变物理模型、边界、Poisson、Beam、Tail、碰撞、阈值或阶段顺序。
4. 所有文件修改在本地工作区完成。集群测试由用户执行；不得通过SSH修改远程文件、编译或提交作业。
5. 没有用户提供的集群结果时，报告写“待集群执行”，不得推测PASS或FAIL。
6. 不得删除、覆盖或清空旧output。

执行：
1. 核对情形B执行情况中B8与主方案第1节一致。
2. 在docs/情形A执行情况.md创建A0报告。
3. 写出主方案第14.1节命令给用户。
4. 按主方案第15节格式报告，未获得集群结果时停止。

报告路径固定为：docs/情形A执行情况.md
```

## A1

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：docs/情形A执行情况.md中的A0必须PASS，否则停止。

最高规则：
0. 绝对遵守主方案第0、5、16节。
1. 本次只执行A1，不得进入A2。
2. 只允许修改src/vpfp_integrator.h、src/vpfp_integrator.cpp、tests/joint_phase_space_midpoint_energy_test.cpp、docs/情形A执行情况.md。
3. 只能增加只读诊断。不得修改flux、Newton、line search、Poisson、G*、energy gate、acceptance、dt或容差。
4. 必须使用最终有限candidate的accepted_residual、accepted_bundle、fields和candidate_fields，不得重建电流或重新推进。
5. 严格实现第5.2至5.10节全部字段、符号、u边界项、paired potential、long double求和和roundoff bound。
6. 所有修改在本地完成；集群编译和运行由用户执行。
7. 没有集群结果时报告“待集群执行”。
8. 不得把failure_code=75改成PASS。

执行：
1. 阅读允许文件中的结构体、advance_joint_midpoint最终诊断块和测试输出。
2. 先列出将修改的函数、位置和字段，再修改。
3. 完成A1全部代码。
4. 本地做静态检查和git diff --check；无MPI环境不得声称集群编译通过。
5. 将第14.2节命令写入报告，等待用户执行。
6. 用户回传结果后按第5.11节验收；未回传则停止。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A2

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A1必须PASS，否则停止。

最高规则：
0. 绝对遵守主方案第0、6、16节。
1. 本次只执行A2，不得进入A3。
2. 只允许修改现有J0/J1测试文件；只有确实无法复用现有target时才允许最小修改CMakeLists.txt。
3. 不得修改任何生产算法。
4. 必须实现a2-zero-residual和a2-injected-residual两个case。
5. 两个case必须直接调用第6.2节列出的生产函数，禁止测试复写生产公式。
6. injected case必须让W_C显著非零，禁止用对称抵消或零场蒙混通过。
7. 修改在本地完成；用户执行集群测试。
8. 没有集群结果时报告“待集群执行”。

执行：
1. 阅读第2、3、6节和现有J0/J1测试。
2. 优先加入joint_phase_space_midpoint_unit_test的--case all。
3. 输出第6.3节全部字段和子门。
4. 本地静态检查和git diff --check。
5. 把第14.3节命令写入报告并等待用户。
6. 回传结果后按tau_C、tau_A和非平凡W_C验收。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A3

```text
你是严格执行型分析代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A1和A2必须PASS，且用户已提供完整result/out/err。否则停止。

最高规则：
0. 本次只执行A3，不修改源码，不进入A4。
1. 严格按第3节和第7节机械分类，不得创造新分类。
2. 最终只能输出一个：A-N_NONLINEAR_RESIDUAL_PROJECTION、A-T_TIME_LAYER_MISMATCH、A-O_FACE_OWNERSHIP_OR_ENDPOINT、A-S_SOURCE_OR_U_BOUNDARY_OR_ASSEMBLY、A-P_POISSON_SCALAR_IDENTITY、INCONCLUSIVE。
3. 多个分类同时满足时必须输出INCONCLUSIVE。
4. 文档修改在本地完成，不得SSH写入或运行集群内容。

执行：
1. 读取A1/A2全部数值。
2. 复算tau_C、tau_A、charge projection mismatch和prediction error。
3. 按第3节逐项排除。
4. 写docs/情形A执行情况.md的A3报告。
5. 若命中A-S，只说明允许进入A-S0，不得解释成真实物理source，不得进入A4。
6. 只有A-N时说明允许A4，但不要自行进入。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A-S0

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A3必须输出A-S_SOURCE_OR_U_BOUNDARY_OR_ASSEMBLY，否则禁止执行。

最高规则：
0. 绝对遵守主方案第0、7A、16节。
1. 本次只执行A-S0，不得进入A-S1、A-FS、A3R或A4。
2. 只允许修改src/vpfp_integrator.h、src/vpfp_integrator.cpp、tests/joint_phase_space_midpoint_energy_test.cpp、docs/情形A执行情况.md。
3. 只能增加只读诊断；不得修改rho装配、Poisson、flux、Newton、line search、G*、energy gate、acceptance、dt或容差。
4. 必须逐速度单元先计算candidate-m_old，再用long double求和；禁止分别求两个大总数再相减。
5. 必须独立输出assembly、transport、u-boundary和potential-weighted分解，严格实现第7A.3至7A.5节。
6. assembly roundoff bound必须包含离子/电子父电荷尺度和gamma_m；禁止继续用旧tau_C判断装配误差。
7. 不得仅凭first_bad_global_ix=0宣称边界问题。
8. 修改在本地完成；集群编译和运行由用户执行。
9. 没有集群结果时报告“待集群执行”。

执行：
1. 阅读A1诊断块、candidate rho装配、m_old/candidate布局和第7A节。
2. 增加第7A.5节全部字段并初始化。
3. 在code75前用同一个最终candidate计算分解。
4. 扩展J1 result输出，保留所有旧字段。
5. 本地静态检查和git diff --check。
6. 把第14.4A命令交给用户。
7. 收到结果后严格按第7A.6节输出A-S-CHARGE_ASSEMBLY_ROUNDOFF、A-S-PHYSICAL_SOURCE_OR_BOUNDARY或A-S-MIXED之一；不要自行进入下一阶段。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A-S1

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A-S0必须唯一分类为A-S-CHARGE_ASSEMBLY_ROUNDOFF，否则禁止执行。

最高规则：
0. 绝对遵守主方案第0、7B、16节。
1. 本次只执行A-S1，不得进入A-FS、A3R或A4。
2. 只允许修改src/vpfp_integrator.cpp、必要时的src/vpfp_integrator.h、现有J0/J1测试和docs/情形A执行情况.md。
3. 只修改J1 candidate evaluate的rho装配；禁止修改EMFields::set_charge_density通用生产路径、Poisson、flux、G*、Newton或阈值。
4. 必须逐速度单元用long double累加state-m_old，并按第7B.3节构造增量rho。
5. 禁止从current divergence构造rho；rho仍必须来自candidate mass。
6. 必须保留absolute form只读对照，不得用对照值覆盖增量结果。
7. 修改在本地完成；集群测试由用户执行。
8. 没有集群结果时报告“待集群执行”。

执行：
1. 精确定位advance_joint_midpoint的candidate rho装配循环。
2. 按第7B.3节替换为稳定增量形式，符号不得改动。
3. 增加第7B.4节只读对照字段。
4. 增加或更新单元测试：固定离子、近中性小扰动、正负扰动，验证增量与高精度参考一致。
5. 本地静态检查和git diff --check。
6. 把第14.4B命令交给用户。
7. 收到结果后按第7B.5节验收；PASS只允许进入A-FS，不得自行进入。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A-FS

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A-S1必须PASS，且 `docs/情形A执行情况.md` 已记录 A-FS 第二轮结果：stable sum被调用、large-baseline fixture输入不可解析、nonzero endpoint为当前模型外结构问题。否则停止。

最高规则：
0. 绝对遵守主方案第0、7C、16节。
1. 本次只执行A-FS-R1，不得进入A3R或A4。
2. 默认只允许修改tests/vpfp_poisson_work_identity_test.cpp、tests/joint_phase_space_midpoint_energy_test.cpp和docs/情形A执行情况.md。
3. 只有缺少两套Gate F诊断字段时，才允许最小修改open_electrostatic_solver.h/.cpp和vpfp_integrator.h/.cpp；只能新增或转发字段，不得再次修改稳定求和公式。
4. 不得修改OpenElectrostaticSolver::solve、reconstruct_phi、build_potential_pairing_field、boundary_energy_work、Poisson stencil、边界、场状态、J1 flux、Newton、G*、energy gate、dt。
5. Gate F按用户批准使用16384*epsilon*scale；禁止把absolute-term-sum加入门槛，禁止继续扩大16384。必须同时输出旧8192门。
6. 禁止通过修改production identity、降低1e4大基线比、降低128输入可解析性系数或删掉失败case来通过。
7. large-baseline的rho_delta必须由long double场增量经离散Gauss差分构造，禁止after.rho-before.rho。
8. 非零endpoint case必须保留真实identity布尔和残差；只允许通过production_scope=0将其从零端点总体门隔离，禁止伪造PASS。
9. 测试必须直接调用生产evaluate_work_identity，禁止复写生产稳定求和公式。
10. evaluate_work_identity调用前后全部场数组必须逐位不变。
11. F10 nontrivial case仍可能因W_transport返回code75；A-FS-R1只验收scalar identity和离散链，不得提前实施A4。
12. 修改在本地完成；集群编译和测试由用户执行。
13. 没有集群结果时报告“待集群执行”。

执行：
1. 先确认7C.3--7C.7稳定求和已经存在；本阶段不得重复重写。
2. 严格按7C.10.2构造确定性的base field和正负delta field，并验证零端点trapezoid积分兼容性。
3. 严格按7C.10.3从long double face增量构造rho_delta；before/after rho只作诊断。
4. 严格按7C.10.4用nextafter计算输入量化界并输出全部字段。
5. 严格按7C.10.5只在固定eta列表中选择首个同时满足1e4和128门的候选；无候选则fixture FAIL，不得改门。
6. 按7C.10.6输出8192/16384两套bound、ratio、pass，兼容alias必须映射16384。
7. 按7C.10.7保留nonzero endpoint真实结果，设置production_scope=0和known_limitation；不得强制identity_pass。
8. 按7C.10.8实现明确的zero_endpoint_production_gate_pass与status聚合，禁止把endpoint identity pass混入零端点门。
9. 扩展F10输出所需两套scalar字段；不改变F10 energy gate。
10. 本地静态检查和git diff --check。
11. 把第14.4C命令交给用户；不得自行在集群运行。
12. 收到全部结果后按7C.10.10逐项验收。smooth-perturbed即使code75，只要满足7C.10.9结构门，不因此判A-FS-R1失败。
13. 任一零端点制造解、状态只读门、J0/B4或scalar identity门失败，立即停止，不得进入A3R。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A3R

```text
你是严格执行型分析代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A-S0、A-S1和A-FS必须PASS并有完整集群结果，否则停止。

最高规则：
0. 本次只执行A3R，不修改源码，不自行进入A4/A5。
1. 严格按第7D节机械分类。
2. 只能输出A-ASSEMBLY_ONLY_RESOLVED、A-N_NONLINEAR_RESIDUAL_PROJECTION、A-S-PHYSICAL_SOURCE_OR_BOUNDARY、A-S-MIXED或INCONCLUSIVE。
3. 如果F10在稳定装配后已通过，必须选择assembly-only并明确跳过A4；禁止为了更漂亮能量账继续加solver门。
4. 文档修改在本地完成，不得运行集群测试。

执行：
1. 读取A-S0/A-S1/A-FS全部result/out/err。
2. 先确认A-FS所有制造解、Gate F、J0/B4均PASS，再比较assembly、transport、W_assembly、W_transport、R_P、R_PJ、R_uJ和F10状态。
3. 按第7D节输出唯一分类。
4. 写docs/情形A执行情况.md；只说明允许的下一阶段，不自行进入。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A4

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A3直接分类或完成A-FS后的A3R重新分类必须是A-N_NONLINEAR_RESIDUAL_PROJECTION，否则禁止执行。若A3R为A-ASSEMBLY_ONLY_RESOLVED，绝对禁止执行A4。

最高规则：
0. 绝对遵守主方案第0、8、16节。
1. 本次只执行A4，不得进入A5。
2. 只允许修改src/vpfp_integrator.h、src/vpfp_integrator.cpp、tests/joint_phase_space_midpoint_energy_test.cpp、docs/情形A执行情况.md。
3. 不得修改Poisson、flux、G*、energy gate 1e-8、max_iterations=20、Beam、Tail、collision或边界。
4. pairing_tolerance固定为1e-9，不得再次放宽；A6-R2已在该门下通过dt、dt/2和10步累计，总能量门仍为1e-8。
5. candidate自身的Poisson-current metric必须进入收敛条件；禁止使用旧candidate或后处理值。
6. phase未通过时保持旧line search；phase通过而pairing未通过时才用第8.4节第二规则。
7. 不得接受违反phase或Poisson门的状态。
8. 修改在本地完成；集群测试由用户执行。
9. 没有集群结果时报告“待集群执行”。

执行：
1. 阅读candidate evaluate、收敛判断、两个line-search路径、iteration record和失败返回。
2. 按第8.2节增加candidate metric。
3. 按第8.3节实现三重门。
4. 按第8.4节同步修改两个line-search路径。
5. 按第8.5节扩展日志。
6. 按第8.6节实现停滞失败；不得提高迭代上限。
7. 本地静态检查和git diff --check。
8. 把第14.5节构建命令写入报告并等待用户。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A5

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置条件二选一：
- A3R为A-ASSEMBLY_ONLY_RESOLVED，执行A5A；
- A4集群构建PASS，执行A5N。
两者均不满足则停止。

最高规则：
0. 本次只执行A5，不得进入A6。
1. A5A只做回归，不得增加pairing-gate-negative或solver hook。
2. A5N才允许修改J1 energy测试和必要的测试专用hook，不得修改生产公式或阈值。
3. A5A覆盖全部J0/B4/A2与两个F10 case；A5N覆盖两个F10 case和pairing-gate-negative。
4. 负测试hook不得成为生产CLI。
5. 非平凡正测试必须真实通过，禁止覆盖结果字段。
5. 修改在本地完成；用户执行集群测试。
6. 没有集群结果时报告“待集群执行”。

执行：
1. 先读取A3R/A4报告，机械选择A5A或A5N，禁止混跑两条路径。
2. A5A不修改代码，直接给用户第14.6节A5A命令。
3. A5N实现三个case及结果字段，并给用户第14.6节A5N命令。
4. 有代码修改时本地静态检查和git diff --check。
5. 回传后严格按第9节验收。
6. 任一必需测试失败，立即停止。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A6

```text
你是严格执行型代码代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A5A或A5N中实际选择的路径必须PASS，否则停止。

最高规则：
0. 本次只执行A6，不得进入A7。
1. 只允许给joint_phase_space_midpoint_energy_test增加--dt-scale及输出字段。
2. 默认值必须是1.0；不得修改生产dt或生产CLI。
3. 参数必须positive且finite，非法值退出码2。
4. dt和dt/2都必须独立PASS，不能用dt/2掩盖dt失败。
5. 必须实现--steps，默认1；运行10步时真实累计每一步signed/absolute energy residual，不得清零或只取最后一步。
6. 10步必须满足主方案第10.1节累计漂移预算；任一步失败则整个case FAIL。100步为A8后的可选非阻断测试，不得要求用户现在运行。
7. 多步初始candidate必须允许上一接受态的signed residual domain；初始evaluate使用allow_negative_probe=true。
8. 最终code-76正性门和negative_tolerance必须保持不变；禁止裁剪上一接受态。
9. 修改在本地完成；用户执行集群测试。
10. 没有结果时报告“待集群执行”。

执行：
1. 实现--dt-scale和--steps并输出actual_dt、dt_scale、accepted_step_count及全部累计字段。
2. 本地静态检查和git diff --check。
3. 给用户第14.7节命令。
4. 回传后按第10节验收；任一档FAIL即停止。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A7

```text
你是严格执行型验收代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A6必须PASS，否则停止。

最高规则：
0. 本次只执行A7，不得进入A8。
1. 原则上不修改生产代码；先给用户第14.8节命令。
2. 集群测试由用户执行；不得SSH写入、编译或提交作业。
3. 缺少1/2/5 rank任一result/out/err时报告“待集群执行”。
4. 若仅多rank失败，只允许在A-O范围检查ownership和endpoint重复计数。
5. 结构恒等式按tau_A/continuity roundoff验收；candidate R_PJ及迭代停止按实际pairing_tolerance验收，禁止要求跨rank R_PJ达到tau_A。
6. 迭代次数可以不同；只要全部在上限内满足phase/Poisson/pairing三重门，不得据此判FAIL。

执行：
1. 核对三档全局网格、初态、dt、case一致。
2. 读取第11节全部字段。
3. 复算tau_A、每档normalized pairing residual、normalized spread和absolute spread允许界。
4. 严格按第11.1--11.3节判定，不得混用舍入门与solver门。
5. PASS时只说明允许A8，不得自行进入。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```

## A8

```text
你是严格执行型最终验收代理。请完整阅读并严格执行：

docs/VPFP_F10情形A_连续性Poisson功配对严格修复实施方案.md

前置：A0、A1、A2、A3、A-S0、A-S1、A-FS、A3R、实际选择的A5路径、A6、A7均PASS；只有A3R为A-N时才额外要求A4和A5N PASS。

最高规则：
0. 本次只执行A8。
1. 不新增算法，只运行并验收第12节和第14.9节全回归。
2. 集群测试由用户执行；不得SSH写入、编译或提交作业。
3. 缺少完整result/out/err时报告“待集群执行”。
4. 任一J0、B4、A2、A-S0、A-S1、A-FS、实际A5路径、A6、A7或F10门失败，A8立即FAIL。
5. 不得把J1 periodic G*直接接入生产开放边界J2。

执行：
1. 给用户第14.9节命令，并列出需保留的A5负测试、A6、A7结果路径。
2. 收到结果后逐项验收第12节。
3. 只有全部PASS才写第12节固定结论。
4. PASS后只说明允许更新根因总结和进入J1 MPI最终门，不自行开始J2/J3。

报告路径固定为：docs/情形A执行情况.md
报告格式严格使用主方案第15节。
```
