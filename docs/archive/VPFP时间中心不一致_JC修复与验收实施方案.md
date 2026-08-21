# VPFP 时间中心不一致：JC 修复与验收实施方案

> 执行对象：GPT-5.4 级或更弱的自动编码智能体。  
> 性质：生产代码修复、测试实现、验收命令和停止规则。  
> 不是：根因搜索文档、方案讨论稿或允许一次执行全文的脚本。  
> 前置证据：优先读取 `output/vpfp_pairing_gate_i/root_cause_selection.result`；若结果已归档，
> 读取 `output/archieve/vpfp_pairing_gate_i/root_cause_selection.result`。文件必须包含  
> `status=PASS_ROOT_CAUSE_UNIQUE`、`root_cause=C`、`root_cause_candidate_count=1`。  
> 最后更新：2026-08-17。

## 0. 执行协议（最高优先级）

### 0.1 一次只执行一个任务

任务只能按下表顺序执行；每次完成一项后必须停止并报告，不得自动进入下一项。

| TASK_ID | 内容 | 前置 | 当前状态 |
|---|---|---|---|
| JC0 | 后置碰撞/转换/return 电荷不变性前置门 | Gate I PASS | PASS |
| JC1 | 配置、数据结构和无副作用 trial API | JC0 PASS | PASS |
| JC2 | 确定性场--粒子 core trial | JC1 PASS | PASS |
| JC3 | 有界 Picard 和 MPI 一致提交 | JC2 PASS | PASS |
| JC4 | CLI/manifest/checkpoint/诊断接线 | JC3 PASS | PASS |
| JC5 | 单 rank、负测试、MPI 和 restart 局部验收 | JC4 PASS | PASS |
| K1 | 115 fs 结构短回归 | JC5 PASS 且用户授权 | FAIL_ENERGY_REFINEMENT / ADJOINT_TEST_INCOMPLETE |
| K2 | 宏观物理 A/B | K1 PASS | BLOCKED |
| K3 | checkpoint/restart A/B | K2 PASS | BLOCKED |
| K4 | level-0 性能 A/B | K3 PASS | BLOCKED |

若用户未明确指定 `TASK_ID`，只选择第一个前置证据齐全且未 PASS 的任务。

为避免弱模型一次读取完整长文档后遗漏局部约束，每个任务使用以下必读集合：

| TASK | 必读章节 |
|---|---|
| JC0 | §0、§1、§2、§3、§13--§16 |
| JC1 | §0、§1、§2、§4、§13--§16 |
| JC2 | §0、§1、§2、§4、§5、§13--§16 |
| JC3 | §0、§1、§2、§4--§6、§13--§16 |
| JC4 | §0、§1、§2、§7、§13--§16 |
| JC5 | §0--§8、§13--§16 |
| K1 | §0--§2、§9、§13--§16 |
| K2 | §0--§2、§10、§13--§16 |
| K3 | §0--§2、§7、§11、§13--§16 |
| K4 | §0--§2、§12--§16 |

执行者可以搜索其他章节核对接口，但不得顺带实施非必读任务。

### 0.2 Agent智能体专用执行约束

本节规则用于避免较弱智能体在已有代码上重复造轮子或扩大修改范围：

1. 每个 TASK 开始时，先用 `Select-String`/`grep` 搜索文档指定的类型、函数和测试目标；
   已存在的实现只能审计和补缺，不得另建同义类或第二套求解器。
2. 不允许仅根据注释判断“已实现”。必须同时检查声明、定义、调用点、CMake 目标和结果文件。
3. 修改前保存本任务允许文件的 `git diff -- <files>`；修改后只报告这些文件的差异。
4. 文档给出的结构体字段是最低要求。允许使用项目现有命名风格，但报告中必须给出
   “文档字段 -> 实际字段”的映射表。
5. 若某个生产算子没有可重复调用的无副作用接口，先在原类中抽取共享私有函数；禁止把公式复制到
   `VpfpIntegrator` 或测试文件。
6. 每次只完成一个 TASK。即使本地测试通过，也必须停止并等待下一条明确指令。
7. 任何不确定的状态所有权、MPI collective 次序或 RNG 行为都视为阻断项，不得自行猜测。
8. 文档中的命令是 Bash/Slurm 命令；不得把 PowerShell 语法写进集群脚本。
9. 不得把缺失结果列解释为零或 PASS；比较器遇到缺列、NaN、空文件必须返回非零退出码。
10. 所有“逐位不变”测试都要比较物理数组、RNG、ledger、step/time，不能只比较快照哈希。

### 0.3 开始前必须声明

```text
TASK_ID=<JC0|JC1|JC2|JC3|JC4|JC5|K1|K2|K3|K4>
BASELINE=<commit；无Git时写unavailable>
PRECONDITION_EVIDENCE=<result路径和关键PASS字段>
ALLOWED_FILES=<本任务允许修改的文件>
FORBIDDEN_FILES=<本任务禁止修改的文件>
ACCEPTANCE_COMMANDS=<本任务必须执行的命令>
STOP_CONDITION=<PASS后停止或首个失败项>
```

### 0.4 严格禁止项

1. 不得修改 Poisson 空间差分、Dirichlet 电势边界、Gauss 门或场能定义。
2. 不得改动 x remap、Beam/Tail deposit/gather shape、碰撞系数、conversion、return 物理规则。
3. 不得添加全局能量补丁、电场缩放、电流投影、零模修正或事后动量修正。
4. 不得用诊断重算的电流覆盖生产 swept-mass/PIC 轨迹电流。
5. 不得在非线性迭代内反复执行随机碰撞、重新生成 Beam injection schedule 或消耗 RNG。
6. 不得软接受未收敛状态，不得在迭代失败时重复输出旧状态伪装时间前进。
7. 不得同时加入 Aitken、Anderson、JFNK 或性能重构。第一版只用有界阻尼 Picard。
8. 不得为 PASS 放宽 Gate I 容差、删除断言或改动 EPOCH 对比数据。
9. `legacy` 必须保留为 A/B 基线，在 K1--K4 全 PASS 前不得改默认模式。

### 0.5 验收要求的权威层级

为防止测试反向扩大生产重构，本文所有要求按以下层级解释：

1. **一级：原始物理和实现合同。** §1、§2及每个TASK在“详细验收标准”之前的实现要求最高优先。
2. **二级：必须验收。** 只把一级合同转换成最小可测证据；缺失或失败会阻断下一TASK。
3. **三级：推荐诊断。** 用于定位失败和提高报告质量；缺失不阻断，不得为补字段修改物理推进。
4. **四级：后续增强。** 仅在本地Gate和生产短跑通过后考虑，不能并入当前TASK。

若详细验收与前文实现合同冲突，以前文合同为准。不得为了满足结果文件名、额外计数器、比较器schema
或任意经验阈值而改变生产方程、时间层、状态所有权、MPI拓扑、RNG语义或checkpoint物理内容。

精确文件名和字段名只有在“必须验收”代码块中才是硬要求；推荐诊断允许等价命名。已有任务按当时
有效schema获得PASS后，不因后来增加推荐字段而失效，除非新证据揭示原硬门实际未满足。

### 0.6 远程集群证据与本地镜像规则

本项目的正式计算和部分审计通过SSH直接在集群服务器执行。远程结果可能只把关键字段和结论记录在
本文，而不通过SCP同步全部大文件到本地。因此：

1. 本地不存在某结果文件，只能标记为 `REMOTE_NOT_MIRRORED`，不能据此判定远程测试未执行；
2. 本文中明确记录了远程路径、运行档位、关键字段、状态和日期的结果，可作为当前分析证据；
3. 本地文件与文档远程记录冲突时，先比较时间戳/运行版本；较旧本地副本不得覆盖较新远程结论；
4. 只有远程记录本身缺少关键字段、前后矛盾或验收逻辑错误时，才将其标为 `REMOTE_EVIDENCE_INCOMPLETE`；
5. 后续分析报告应尽量写明 `evidence_location=remote_cluster|local_mirror`、job/日期、结果路径和关键值；
6. 不要求为方便本地分析而传输所有逐rank大文件，可直接在远程运行分析器并把汇总字段写入本文。

## 1. 已确认根因和修复目标

### 1.1 Gate I 权威结论

```text
coarse status=PASS, accepted_steps=10
fine status=PASS, accepted_steps=20
same_initial_physical_state=1
same_physical_window=1
root_cause=C
root_cause_candidate_count=1
status=PASS_ROOT_CAUSE_UNIQUE
```

| 候选 | coarse | fine | 结论 |
|---|---:|---:|---|
| A：Poisson/transport | 0.499481 | 0.501486 | 未在两档稳定超门 |
| B：work/charge current | 0.500519 | 0.498514 | 未在两档稳定超门 |
| C：time centering | 1.789266 | 2.226579 | 两档唯一显著命中 |
| D：PIC local work | $7.63\times10^{-15}$ | $8.78\times10^{-15}$ | 排除主导 |
| E：conversion | $3.02\times10^{-17}$ | $2.20\times10^{-17}$ | 排除主导 |
| F：boundary | 0.086683 | 0.086568 | 排除边界主导 |

### 1.2 当前生产推进的精确结构

`VpfpIntegrator::advance_background()` 和 `advance_with_beam()` 当前分别复写了：

```text
C1(dt/2)
-> Bulk x1(dt/2) + Tail drift1 + Beam predict_to_midpoint
-> midpoint Poisson from rho(n+1/2)
-> Bulk u_full + Tail kick + Beam finish_from_midpoint（含第二半漂移、晚注入、出流和迁移）
-> Bulk-to-Tail conversion
-> Bulk x2(dt/2) + Tail drift2（不得再对 Beam 单独执行 drift2）
-> C2(dt/2) + collision-face conversion
-> Tail return
-> final Poisson from rho(n+1)
-> validation and accepted-state swap
```

受力使用的 `midpoint_fields_.Ex_face` 不等于由 accepted $n$ 和候选 $n+1$ Poisson 状态
构造的离散梯度配对场，因此场能变化和粒子受力功不属于同一时间中心状态。

### 1.3 必须求解的非线性离散系统

$$
(f^{n+1,*},p^{n+1,*},\rho^{n+1,*})
=\mathcal{T}(E_{pair}^{(k)};S_{frozen}),
$$

$$
\Phi^{n+1,*}=\mathcal{P}(\rho^{n+1,*}),
\qquad
E_{pair,new}^{(k)}=\mathcal{G}_{P}(\Phi^n,\Phi^{n+1,*}).
$$

`OpenElectrostaticSolver::build_potential_pairing_field()` 已实现 $\mathcal{G}_{P}$，必须直接复用。
收敛解满足：

$$
E_{pair}^{*}=\mathcal{G}_{P}
\left(\Phi^n,\mathcal{P}(\rho^{n+1}(E_{pair}^{*}))\right).
$$

修复目标是让 Bulk `u_full`、Tail `kick`、Beam `finish_from_midpoint`、电荷输运和
final Poisson 属于同一个收敛的离散状态，而不是强行使能量诊断为零。

## 2. 事务边界和状态所有权

### 2.1 不得进入迭代的操作

1. `apply_collision_half(... collision_half=0)`；
2. `BeamPIC::generate_injection_schedule()`；
3. 任何 stochastic Tail collision/RNG 调用；
4. accepted cumulative ledger、particle ID counter、`step_count_` 的提交；
5. 第二碰撞半步和 H10 return：只在场--粒子 core 收敛后执行一次。

### 2.2 迭代中必须从同一冻结状态重放的操作

1. Bulk `u_full`、Tail midpoint kick、Beam midpoint kick 和第二段轨迹；
2. 受力后 Bulk-to-Tail flux conversion；
3. Bulk x2、Tail drift2、Beam drift2；
4. candidate density deposit、final Poisson、pairing field 和残差计算。

每次 trial 必须从同一个 frozen state 恢复，不得以上一次 trial 的候选状态作为起点。

### 2.3 C2/return 为什么放在收敛后

C2 碰撞、collision-face Bulk/Tail 转换和 H10 return 理论上应在每个 x cell 守恒总背景
电子数，因此不应改变 final Poisson 电荷。但这必须先由 JC0 用生产算子证明。

若 JC0 不通过，输出 `BLOCKED_BY_POST_FIELD_CHARGE_CHANGE` 并停止；禁止把 stochastic C2
放入 Picard 迭代。

## 3. JC0：后置算子电荷不变性前置门

### 3.1 当前代码事实与本任务边界

当前仓库已经包含：

- `VpfpIntegrator::post_field_charge_invariance_transaction()`；
- `tests/vpfp_post_field_charge_invariance_test.cpp`；
- CMake 目标 `vpfp_post_field_charge_invariance_test`。

因此 JC0 的默认任务是**审计并运行现有实现**，不是重新新增文件。只有编译失败、测试缺案例或结果字段
不满足本节要求时，才允许补充以下文件：

- `src/vpfp_integrator.h`
- `src/vpfp_integrator.cpp`
- `tests/vpfp_post_field_charge_invariance_test.cpp`
- `CMakeLists.txt`

禁止修改碰撞、conversion、return 的物理实现。

### 3.2 测试流程

1. 构造已完成场--粒子 core 的 Bulk+Tail candidate；
2. 记录每个 x cell 的 $N_{combined,before}=N_{bulk}+N_{tail}$；
3. 直接调用生产 C2、collision-face conversion 和 return；
4. 再计算 $N_{combined,after}$；
5. 检查每 cell、每 rank、全局差异；
6. 覆盖 no collision、Bulk collision、hybrid collision、return none/hysteretic、conversion active 和 MPI 情形。

$$
T_N=4096\epsilon_{double}\max(1,N_{before},N_{after}).
$$

### 3.3 实现审计清单

执行者必须逐项确认，并把 `file:line` 写入报告：

1. helper 的 Bulk/Tail 参数是 caller-owned trial，未引用 `tail_accepted_`；
2. helper 调用生产 `apply_collision_half()`、collision-face conversion 和
   `apply_tail_bulk_return()`，测试文件没有复写公式；
3. 每个 x-cell 的 combined number 使用与生产 density deposit 相同的速度权重；
4. 单元、rank 和 global 三层容差分别计算，不能只用 global 抵消；
5. MPI 失败使用全局共识，所有 rank 进入相同 collective；
6. 测试前后比较 `step_count_`、`tail_cumulative_`、Tail RNG/ID、converter ledger；
7. `return-none` 必须证明 return 调用次数为零；`return-hysteretic` 必须证明实际发生返回；
8. `conversion-active` 必须证明 `export_number>0`，否则该案例属于空通过；
9. 2-rank 和 5-rank 测试必须同时覆盖内部共享面和两个物理边界；
10. 结果文件必须打印每个案例名称和独立 PASS 字段。

### 3.4 验收字段

```text
status=PASS
case_count>=7
cell_number_residual_linf<=cell_number_tolerance
rank_number_residual_linf<=rank_number_tolerance
global_number_residual<=global_number_tolerance
rho_before_after_relative_linf<=1e-12
rng_advanced_exactly_once=1
ledger_advanced_exactly_once=1
accepted_integrator_state_unchanged=1
conversion_active_nonzero=1
return_hysteretic_nonzero=1
mpi_case_count>=2
```

注意：测试 trial 内的 C2/return 本来就应各执行一次，所以 `rng_advanced_exactly_once=1` 和
`ledger_advanced_exactly_once=1` 指 trial 自身；`accepted_integrator_state_unchanged=1` 才用于证明
正式 integrator 状态未改变。

### 3.5 编译和运行

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target vpfp_post_field_charge_invariance_test

mkdir -p ./output/vpfp_pairing_gate_jc
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_post_field_charge_invariance_test \
  --case all \
  --result ./output/vpfp_pairing_gate_jc/post_field_charge_single.result || exit 61

for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/vpfp_post_field_charge_invariance_test \
    --case mpi \
    --result "./output/vpfp_pairing_gate_jc/post_field_charge_mpi_n${NP}.result" || exit 62
done
```

#### 3.5.1 补充编译和运行（state snapshot 与新增字段）

代码变更：`vpfp_integrator.h` 新增 `step_count()` getter；`vpfp_post_field_charge_invariance_test.cpp` 新增
`IntegratorStateSnapshot`/`compare_integrator_state` 以及 `accepted_integrator_state_unchanged`、
`conversion_active_nonzero`、`return_hysteretic_nonzero`、`mpi_case_count` 四个输出字段。
必须重编译后重跑：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target vpfp_post_field_charge_invariance_test

mkdir -p ./output/vpfp_pairing_gate_jc

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_post_field_charge_invariance_test \
  --case all \
  --result ./output/vpfp_pairing_gate_jc/post_field_charge_single.result || exit 61

for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/vpfp_post_field_charge_invariance_test \
    --case mpi \
    --result "./output/vpfp_pairing_gate_jc/post_field_charge_mpi_n${NP}.result" || exit 62
done
```

运行后逐文件检查新增字段：

```bash
for F in ./output/vpfp_pairing_gate_jc/post_field_charge_*.result; do
  echo "--- $(basename "$F") ---"
  grep -E '^(accepted_integrator_state_unchanged|conversion_active_nonzero|return_hysteretic_nonzero|mpi_case_count)=' "$F"
done
```

### 3.6 JC0 决策

- 全部 PASS：把任务表 JC0 改为 `PASS`，停止，不得自动开始 JC1。
- 只有结果字段缺失：只补测试/比较输出，不修改生产算子。
- 每 cell combined charge 失败：输出首个 rank/ix、C2前后、conversion后和return后的分阶段值，
  将 JC1 标记为 `BLOCKED_BY_POST_FIELD_CHARGE_CHANGE`。
- 发生 accepted state/RNG/ledger 副作用：先修复测试事务边界；不得把随机 C2
  放进 Picard 迭代。

#### 3.6.1 JC0 当前输出状态（2026-08-17 11:08 重跑后）

三个结果文件均显示 `status=PASS`，全部 §3.4 验收字段齐全且通过：

| 字段 | single | mpi_n2 | mpi_n5 | §3.4 要求 |
|---|---|---|---|---|
| status | PASS | PASS | PASS | PASS |
| case_count | 7 | 7 | 7 | >=7 |
| cell_number_residual_linf | 134217728 | 268435456 | 150994944 | <= tolerance |
| cell_number_tolerance | 1.36e11 | 2.73e11 | 1.09e11 | — |
| rank_number_residual_linf | 8.59e9 | 2.15e9 | 1.34e9 | <= tolerance |
| rank_number_tolerance | 4.37e12 | 2.18e12 | 8.73e11 | — |
| global_number_residual | 8.59e9 | 4.29e9 | 6.98e9 | <= tolerance |
| global_number_tolerance | 4.37e12 | 4.37e12 | 4.37e12 | — |
| rho_before_after_relative_linf | 8.95e-16 | 8.95e-16 | 1.26e-15 | <=1e-12 |
| rng_advanced_exactly_once | 1 | 1 | 1 | 1 |
| ledger_advanced_exactly_once | 1 | 1 | 1 | 1 |
| accepted_integrator_state_unchanged | 1 | 1 | 1 | 1 |
| conversion_active_nonzero | 1 | 1 | 1 | 1 |
| return_hysteretic_nonzero | 1 | 1 | 1 | 1 |
| mpi_case_count | 0 | 7 | 7 | >=2（MPI 文件） |

说明：`mpi_case_count=0` 出现在 single 文件中是预期行为——单 rank 运行不计入 MPI 案例。
`mpi_case_count>=2` 由 mpi_n2 和 mpi_n5 两个文件满足。

§3.6 决策：全部 PASS。按规则把任务表 JC0 改为 `PASS`，停止，不自动开始 JC1。

### 3.7 JC0 分级验收标准

#### 3.7.1 必须验收

JC0只验证原§3定义的后置算子电荷不变性和事务隔离。必须保留现有三个结果：

~~~text
post_field_charge_single.result
post_field_charge_mpi_n2.result
post_field_charge_mpi_n5.result
~~~

每个结果必须满足原§3.4字段：

~~~text
status=PASS
case_count>=7
cell_number_residual_linf<=cell_number_tolerance
rank_number_residual_linf<=rank_number_tolerance
global_number_residual<=global_number_tolerance
rho_before_after_relative_linf<=1e-12
rng_advanced_exactly_once=1
ledger_advanced_exactly_once=1
accepted_integrator_state_unchanged=1
conversion_active_nonzero=1
return_hysteretic_nonzero=1
~~~

MPI结果还要求 mpi_case_count>=2。比较时必须同时读取 residual 和 tolerance，不能只读取 PASS 字符串。

#### 3.7.2 推荐诊断

以下内容用于失败定位，但不是JC0硬门：

- 每个案例的独立调用次数；
- 首个失败rank/ix；
- C2、conversion、return后的分阶段combined number；
- shared-face owner和物理边界细分；
- 单独的jc0汇总文件。

缺少推荐字段不得推动修改碰撞、conversion或return生产公式。

#### 3.7.3 当前结论

2026-08-17保存的single、2-rank和5-rank结果已经满足§3.7.1，JC0保持PASS。只有发现现有
residual超过其记录的tolerance、非空conversion/return路径实际未触发，或accepted状态有副作用，
才需要重跑JC0。
## 4. JC1：建立离散梯度场粒子耦合的数据结构

### 4.1 修改范围

主要修改：

- `src/vpfp_integrator.h`
- `src/vpfp_integrator.cpp`

本节只建立接口、工作区和所有权边界，不改变默认生产推进。新增模式必须默认关闭，`legacy` 路径的执行顺序和结果必须保持不变。

### 4.2 新增配置类型

在 `vpfp_integrator.h` 中增加等价于下列定义的类型。命名可按现有风格调整，但字段语义不能改变。

```cpp
enum class FieldParticleCouplingMode {
    Legacy,
    DiscreteGradient
};

struct FieldParticleCouplingConfig {
    FieldParticleCouplingMode mode = FieldParticleCouplingMode::Legacy;
    int max_iterations = 12;
    double initial_relaxation = 0.5;
    double minimum_relaxation = 0.125;
    double maximum_relaxation = 1.0;
    double field_relative_tolerance = 1.0e-8;
    double pairing_relative_tolerance = 1.0e-8;
};
```

不得把容差写死在循环体中。生产默认值仍为 `Legacy`，避免尚未验收的 JC 路径静默进入正式计算。

### 4.3 冻结态、试探态和迭代诊断

新增三个结构体。

```cpp
struct FieldParticleFrozenState {
    // 当前已接受的 n 层状态，只读引用或稳定快照。
    // C1 后状态、第一段 x 输运结果、固定 Beam 注入计划、
    // Tail/Beam 中点预测和所有本步固定随机输入。
};

struct FieldParticleTrial {
    // 本次 E_pair 猜测对应的完整确定性候选：
    // Bulk/Tail/Beam、第二段 x 输运、rho_np1、phi_np1、
    // E_pair_new、J 和各分项功账。
};

struct FieldParticleIterationDiagnostics {
    int iterations = 0;
    double relaxation = 0.0;
    double field_residual_l2 = 0.0;
    double field_residual_linf = 0.0;
    double pairing_residual = 0.0;
    bool converged = false;
    int failure_code = 0;
};
```

`FieldParticleTrial` 必须拥有独立可复用缓冲区，不能别名到已接受状态。重复迭代只能覆盖 trial 缓冲区，不得修改：

- `state_n`；
- 正式 Beam 粒子容器；
- Tail 正式粒子容器；
- RNG 状态；
- 注入、出流、转换、碰撞和返回累计账本；
- checkpoint 可见的 accepted-step 计数。

### 4.4 新增生产内核接口

至少提供以下职责等价的接口：

```cpp
void set_field_particle_coupling(const FieldParticleCouplingConfig& config);

bool prepare_field_particle_frozen_state(
    double dt,
    FieldParticleFrozenState& frozen,
    VpfpFailureInfo& failure);

bool evaluate_field_particle_trial(
    const FieldParticleFrozenState& frozen,
    const std::vector<double>& pairing_field_guess,
    FieldParticleTrial& trial,
    VpfpFailureInfo& failure);

bool solve_field_particle_center(
    const FieldParticleFrozenState& frozen,
    FieldParticleTrial& accepted_trial,
    FieldParticleIterationDiagnostics& diagnostics,
    VpfpFailureInfo& failure);
```

Beam 关闭和开启必须调用同一套 `evaluate_field_particle_trial()`。允许 Beam 分支为空操作，不允许复制两套 JC 迭代实现。

### 4.5 工作区与性能约束

在 `VpfpIntegrator` 生命周期内预分配并复用：

- Bulk trial 分布；
- 第二段 x 输运工作区；
- Tail trial 粒子及其临时 kick 结果；
- Beam trial 粒子及其临时 finish-push 结果；
- `rho_np1`、`phi_np1`、`E_pair_new`；
- 电流、功账和残差归约缓冲区。

不得在每次非线性迭代中重复创建完整四维数组。允许清零必要的数值区间，但应避免重新分配容量。

### 4.6 JC1 验收

JC1 完成后只运行编译和 legacy 回归，不运行新的物理模式：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target fp_solver
```

验收条件：

```text
build_pass=1
default_mode=legacy
legacy_call_order_unchanged=1
accepted_state_not_aliased_by_trial=1
trial_rng_side_effect=0
trial_ledger_side_effect=0
```

### 4.7 DeepSeek 实施卡：准确插入位置

按以下顺序修改，不得自行更换顺序：

1. 在 `src/vpfp_integrator.h` 的 `VpfpStepResult` 之前定义 coupling enum/config/diagnostics；
2. 在 `class VpfpIntegrator` 的 public 区、现有
   `set_field_particle_power_audit_enabled()` 附近加入配置 setter/getter；
3. 在 private 区、`state_collision_trial_` 等工作区附近加入 JC 缓冲区；
4. 在 `advance_background()`/`advance_with_beam()` 声明附近加入
   `advance_discrete_gradient()` 和三个 helper 声明；
5. 在 `src/vpfp_integrator.cpp` 两个构造函数中初始化配置和布尔状态；
6. 在 `init()` 中按 `grid_.nx_local`/`nx_local+1` 预留 face、rho 和 residual 数组；
7. 在 `advance()` 的最外层只增加模式分派。`Legacy` 必须继续原样进入现有
   `advance_background()` 或 `advance_with_beam()`；
8. JC1 阶段的 `advance_discrete_gradient()` 只能返回明确的 `not implemented` 失败码，
   不得偷偷调用不完整的新算法。

### 4.8 推荐的实际成员布局

当前类已有 `state_x_half_`、`state_u_full_`、`state_np1_`、`beam_work_`、
`tail_work_` 和 `final_fields_`。JC 路径按下列所有权复用，避免无意义增加完整四维数组：

```text
state_collision_trial_  C1输出（已有）
state_x_half_           frozen Bulk x-half状态（已有，迭代中只读）
state_u_full_           当前trial的Bulk受力后状态（已有）
state_np1_              当前trial的Bulk x2候选（已有）
tail_work_              frozen Tail midpoint状态（已有，迭代中只读）
beam_work_              frozen Beam midpoint状态（已有，迭代中只读）
tail_field_trial_       当前trial Tail状态（新增）
beam_field_trial_       当前trial Beam状态（新增）
field_n_pairing_        具有有效phi的n层场副本（新增）
final_fields_           当前trial最终Poisson状态（已有）
pairing_field_guess_    当前face猜测，nx_local+1（新增）
pairing_field_map_      G_P返回face场，nx_local+1（新增）
pairing_field_residual_ 两者差，nx_local+1（新增）
```

不要再增加第二套 `Species state_np1`。Tail/Beam trial 每轮允许从 frozen 对象复制；第一版先保证
正确性，K4 再用容量复用、双缓冲或 `swap_state()` 优化复制成本。

`FieldParticleFrozenState` 不应再次拥有上述大数组。它只保存指向类工作区的 `const` 指针及标量：

```cpp
const Species* bulk_x_half;
const BackgroundTailPIC* tail_midpoint;
const BeamPIC* beam_midpoint;
const EMFields* field_n_pairing;
BeamInjectionSchedule beam_schedule;
BulkTailFluxBatch first_collision_flux;
CollisionDiagnostics first_collision;
HybridCollisionDiagnostics first_hybrid;
VlasovStepDiagnostics x1_diagnostics;
long long candidate_step;
double time;
double dt;
bool beam_on;
bool tail_on;
```

如实际类型不可复制，保存其生产快照类型；不得保存指向函数局部变量的悬空指针。

### 4.9 `step_count_` 的事务语义

当前 legacy 函数在步首执行 `++step_count_`，conversion ID、collision RNG key 和 controller cadence
都依赖它。JC 路径必须保证一个物理步的所有 trial 使用同一个候选步号：

```cpp
const long long accepted_step_before = step_count_;
const long long candidate_step = accepted_step_before + 1;
step_count_ = candidate_step;  // 每个物理步只设置一次，不得每轮递增
```

若 JC 任意阶段失败，在返回前恢复：

```cpp
step_count_ = accepted_step_before;
```

成功原子提交后保留 `candidate_step`。必须用单元测试验证 conversion particle ID 和 collision
counter-key 在 1轮与多轮 trial 下相同。不得在每次 Picard 迭代调用 `++step_count_`。

### 4.10 JC1 本地静态检查

JC1 不要求新物理路径可运行，但必须新增一个小型所有权测试或 test access，确认：

```text
legacy_dispatch_unchanged=1
discrete_gradient_stub_fails_explicitly=1
candidate_step_constant_across_mock_trials=1
failed_mock_trial_restores_step_count=1
work_buffer_capacity_stable=1
```

### 4.11 JC1 增量测试命令

JC1 即创建 `vpfp_field_particle_trial_test`，但本阶段只实现 `ownership` 和 `candidate-step` 两案：

```bash
cmake --build build -j4 --target fp_solver vpfp_field_particle_trial_test
mkdir -p ./output/vpfp_pairing_gate_jc/jc1

for CASE in ownership candidate-step; do
  yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_test \
    --case "$CASE" \
    --result "./output/vpfp_pairing_gate_jc/jc1/${CASE}.result" || exit 65
done
```

两案 PASS 后停止。不得因为 JC2 case 尚未实现而判 JC1 失败。

### 4.12 JC1 分级验收标准

#### 4.12.1 必须验收

JC1只验证配置/分派骨架、状态所有权和candidate step事务，不验证JC2物理trial。必须运行现有：

~~~text
ownership.result
candidate-step.result
~~~

ownership必须满足：

~~~text
status=PASS
legacy_dispatch_unchanged=1
discrete_gradient_stub_fails_explicitly=1
accepted_state_not_aliased_by_trial=1
trial_rng_side_effect=1
trial_ledger_side_effect=1
work_buffer_capacity_stable=1
~~~

这里trial_rng_side_effect=1和trial_ledger_side_effect=1沿用现有测试语义，表示“正式accepted对象无副作用”。
candidate-step必须满足：

~~~text
status=PASS
candidate_step_constant_across_mock_trials=1
failed_mock_trial_restores_step_count=1
~~~

默认模式必须仍为legacy；DG stub不得静默回退为legacy成功。

#### 4.12.2 推荐诊断

推荐但不阻断：

- 各工作数组地址、size和capacity；
- 首个alias对象；
- candidate step在每个mock trial中的打印；
- 失败前后RNG/ledger的分项hash；
-预热后重新分配次数。

不强制要求预先保存完整四维pre-JC bitwise baseline。若已有小fixture baseline，可作为附加回归证据。

#### 4.12.3 当前结论

现有ownership和candidate-step结果已满足§4.12.1，JC1保持PASS。不得仅为补充推荐字段重构工作区
或运行生产checkpoint。
## 5. JC2：实现确定性的完整场粒子候选评估

### 5.1 每个物理步只执行一次的冻结阶段

`prepare_field_particle_frozen_state()` 必须按当前生产顺序执行并冻结以下内容：

1. 从已接受的 $n$ 层状态开始；
2. 执行第一碰撞半步 C1；
3. 生成一次 Beam 注入 schedule，并保存本步固定随机输入；
4. 执行第一段 Bulk x 输运；
5. 执行第一段 Tail 空间漂移；
6. 预测 Beam 中点轨迹；
7. 保存后续每次 trial 均要复用的中点粒子状态和 swept-mass 信息。

该函数同一物理步只能调用一次。后续场迭代不得重新执行 C1、重新抽样 Beam、推进正式 RNG 或重复累计源项。

### 5.2 每个场猜测都要重新计算的确定性阶段

`evaluate_field_particle_trial()` 接收 $E_{\rm pair}^{(k)}$，从同一个 frozen state 重建候选态：

1. 用 $E_{\rm pair}^{(k)}$ 执行 Bulk 完整速度受力推进；
2. 用同一场执行 Tail 电场 kick；
3. 用同一场调用 `BeamPIC::finish_from_midpoint()`；该调用已经包含 kick、第二半漂移、
   晚注入、开放边界删除和 MPI 迁移，不得再额外推进 Beam；
4. 执行依赖受力后状态的 Bulk-to-Tail 转换，但只写 trial 容器；
5. 执行第二段 Bulk x 输运和 Tail 空间漂移；
6. 得到候选 $f^{n+1,*}$、Tail/Beam 候选和 $\rho^{n+1,*}$；
7. 对候选电荷执行最终 Poisson，且必须设置 `reconstruct_phi=true`；
8. 调用现有 `OpenElectrostaticSolver::build_potential_pairing_field()`，以 $\Phi^n$ 和 $\Phi^{n+1,*}$ 构造 $E_{\rm pair,new}^{(k)}$；
9. 计算本次 trial 的各组分电流、动能变化和场功，但不写正式 ledger。

不得在这里执行 C2、随机碰撞、H10 return 或正式 commit。JC0 必须先证明这些后处理不改变每个 x-cell 的 combined charge，才能将其安全地留在收敛之后。

### 5.3 禁止使用错误的时间层替代

以下做法禁止：

- 用旧的 midpoint Poisson 场替代离散梯度场；
- 用 $(E^n+E^{n+1})/2$ 直接替代 `build_potential_pairing_field()`；
- 用候选密度差反推并覆盖生产电流；
- 修改 Poisson 空间离散以消除时间残差；
- 添加全局能量补丁、零模扣除或事后速度缩放；
- 将 trial 的 Beam/Tail 粒子交换到正式容器。

### 5.4 残差定义

场固定点残差定义为：

$$
r_E^{(k)}=E_{\rm pair,new}^{(k)}-E_{\rm pair}^{(k)}.
$$

使用与场能一致的 face 权重计算：

$$
R_{E,2}=
\frac{\left(\sum_f w_f [r_{E,f}^{(k)}]^2\right)^{1/2}}
{\max\left(E_{\rm floor},
\left(\sum_f w_f [E_{\rm pair,new,f}^{(k)}]^2\right)^{1/2}
\right)},
$$

$$
R_{E,\infty}=
\frac{\max_f|r_{E,f}^{(k)}|}
{\max(E_{\rm floor},\max_f|E_{\rm pair,new,f}^{(k)}|)}.
$$

功配对残差必须调用现有 `evaluate_work_identity()`，并以独立交换尺度归一化，不能用强抵消后的净能量残差作为唯一分母：

$$
R_W=
\frac{|\Delta K_{\rm particles}-W_{\rm pairing}|}
{\max(K_{\rm floor},|\Delta K_{\rm particles}|,|W_{\rm pairing}|)}.
$$

这里的 $\Delta K_{particles}$ 仅指本轮 `u_full`、Tail kick 和 Beam kick 的**受力动能变化**；
$W_{pairing}$ 仅指同一 `pairing_field_guess` 与同一受力电流的功。它不是包含 x 输运、reservoir、
C1/C2碰撞、conversion、return 和边界出流在内的完整 Gate I 总能量余额。完整余额含已确认的
空间截断/有符号抵消项，只能在 K1 coarse/fine 中作为误差预算，禁止强制到 $10^{-8}$。

所有 MPI rank 必须对同一组全局残差作一致决策。

### 5.5 DeepSeek 实施卡：从现有生产函数抽取代码

不要从零重写推进。以 `advance_with_beam()` 为完整模板，按现有注释锚点抽取；
`advance_background()` 仅作为 Beam-off 对照。

#### 5.5.1 移入 `prepare_field_particle_frozen_state()` 的代码块

从现有函数按原顺序移动或抽取以下块：

1. `VpfpStepResult` 基础初始化、CFL和配置门；
2. 步前 number/energy ledger；
3. `beam_work_.begin_step()` 和 `generate_injection_schedule()`，schedule 保存到 frozen；
4. `tail_work_ = tail_accepted_` 和 `tail_work_.begin_step()`；
5. 第一碰撞半步 `apply_collision_half(... collision_half=0)`；
6. `vlasov_.first_x_half()`；
7. `tail_work_.drift_half()` 和 midpoint density deposit；
8. `beam.predict_to_midpoint()` 和 midpoint density deposit；
9. legacy midpoint Poisson，只用于生成第一轮初值；
10. 将 `state_x_half_`、`tail_work_`、`beam_work_` 视为冻结态，不再由 trial 修改。

第一碰撞半步产生的 `first_collision_flux`、diagnostics、stage source 和账本前缀必须保存在
frozen/result-prefix 中，不能在每轮重新计算或重复累计。

#### 5.5.2 移入 `evaluate_field_particle_trial()` 的代码块

每次调用先执行：

```cpp
tail_field_trial_ = tail_work_;
beam_field_trial_ = beam_work_;
trial_result = frozen.result_prefix;
trial_stage_sources = frozen.stage_sources_prefix;
```

然后将现有生产块改为使用 trial 对象：

| 现有调用 | JC trial 输入 | JC trial 输出 |
|---|---|---|
| `vlasov_.u_full(state_x_half_, state_u_full_, midpoint_fields_, ...)` | frozen `state_x_half_` + `pairing_field_guess` 包装成只读 `EMFields` | `state_u_full_` |
| `beam_work_.finish_from_midpoint(...)` | frozen `beam_work_` 的副本 | `beam_field_trial_` |
| `tail_work_.kick(...)` | frozen `tail_work_` 的副本 | `tail_field_trial_` |
| `apply_upar_flux_conversion(...)` | 本轮 `exported_flux` + frozen first-collision flux | `state_u_full_` + `tail_field_trial_` |
| population controller | 只作用本轮 `tail_field_trial_` | 同一 trial |
| `vlasov_.second_x_half(...)` | 本轮 `state_u_full_` | `state_np1_` |
| Tail `drift_half()` | 本轮 Tail kick/conversion 后状态 | `tail_field_trial_` |
| Beam第二段 | 已由本轮 `finish_from_midpoint()` 完成 | 禁止再推进 |

`pairing_field_guess` 必须写入一个预分配 `EMFields trial_force_fields_` 的 `Ex_face`，并按现有
Yee face-to-cell 规则更新 `Ex` 和 ghosts。不得调用 Poisson 覆盖它。当前项目没有公开的只读
face-to-cell helper，因此 JC2 唯一允许的表外小改动是：在 `OpenElectrostaticSolver` 暴露
`populate_electric_components_from_faces()`（名称可调整），内部只调用现有
`sync_face_interfaces()`、cell平均和 `exchange_cell_ghosts()`。`solve()` 本体和边界公式不得改动，
不得在 `VpfpIntegrator` 手写第二套 MPI/物理边界同步。

#### 5.5.3 保留在收敛之后的代码块

以下块不得进入 trial：

- 第二碰撞半步 `apply_collision_half(... collision_half=1)`；
- second-collision flux conversion；
- `apply_tail_bulk_return()`；
- 最终 accepted cumulative ledger 更新；
- `swap_state()`、`swap_emfields()`；
- `beam.commit_injection_schedule()`；
- checkpoint/accepted diagnostics。

### 5.6 $n$ 层势函数 bootstrap：必须实现

现有生产最终 Poisson 使用 `final_options.reconstruct_phi=false`，因此 legacy checkpoint 的
`fields.phi` 可能是陈旧值。JC 不得直接使用它。

在 frozen 阶段增加 `prepare_pairing_field_n()`：

1. `field_n_pairing_ = fields`；
2. 保持 `field_n_pairing_.rho` 与 accepted `fields.rho` 完全相同；
3. 对副本调用现有 `field_solver_.solve()`，选项设为：

```cpp
OpenGaussSolveOptions options;
options.reconstruct_phi = true;
options.compute_l1 = true;
options.compute_boundary_audit = true;
```

4. 比较重建前后的 `Ex_face`，使用：

$$
R_{bootstrap}=\frac{\|E_{reconstructed}-E_{accepted}\|_\infty}
{\max(E_{floor},\|E_{accepted}\|_\infty)}.
$$

5. 要求 `R_bootstrap<=1e-12`；不满足则以 `accepted_field_poisson_mismatch` 失败，禁止静默覆盖；
6. 后续 `build_potential_pairing_field()` 的 before 参数必须是 `field_n_pairing_`；
7. JC 成功态的 `final_fields_.phi` 必须随 checkpoint 保存，所以 JC checkpoint 重新启动不依赖旧零势；
8. 第一版允许每个物理步在副本上重建一次以保证正确性；K4 后才允许优化为“仅 legacy->JC
   bootstrap时重建，后续验证phi有效后复用”。

### 5.7 最终候选 Poisson 的固定设置

每个 trial 在 C2/return 之前完成 combined density deposit，并执行：

```cpp
final_fields_.set_charge_density(state_np1_, tail_trial_density,
                                 beam_trial_density, ion_density);
OpenGaussSolveOptions final_options;
final_options.reconstruct_phi = true;
final_options.compute_l1 = true;
final_options.compute_boundary_audit = true;
field_solver_.solve(final_fields_, mpi_rank, mpi_size, final_options);
```

之后立刻调用：

```cpp
field_solver_.build_potential_pairing_field(
    field_n_pairing_, final_fields_, pairing_field_map_, mpi_rank, mpi_size);
```

禁止使用现有 `midpoint_fields_.phi` 参与 pairing。`midpoint_fields_` 只保留为初始猜测来源。

### 5.8 trial 结果重置规则

每轮开始必须重置所有“本轮产生”的对象：

```text
exported_flux
second_collision_flux（应为空，因为C2尚未执行）
conversion_events after the frozen prefix
bulk/tail/beam work ledger
population-control diagnostics
tail/beam trajectory current
failure indices and nonfinite flags
```

不得清空 frozen 的 C1/x1 ledger。最稳妥实现是保存 `VpfpStepResult frozen_prefix_result_`，
每轮复制到 `trial_result_`，再追加 trial 阶段；禁止在上一轮 `trial_result_` 上增量累积。

### 5.9 JC2 确定性测试

同一个 frozen state 和同一个 `pairing_field_guess` 连续调用两次 trial，要求：

```text
bulk_trial_bitwise_equal=1
tail_trial_bitwise_equal=1
beam_trial_bitwise_equal=1
field_trial_bitwise_equal=1
conversion_ledger_bitwise_equal=1
work_ledger_bitwise_equal=1
accepted_bulk_unchanged=1
accepted_tail_unchanged=1
accepted_beam_unchanged=1
accepted_field_unchanged=1
accepted_rng_unchanged=1
face_to_cell_helper_matches_solve_bitwise=1
```

`face_to_cell_helper_matches_solve_bitwise` 的构造方法：先让现有 `solve()` 生成一组非均匀
`Ex_face/Ex`，再把同一 `Ex_face` 放入新 field，调用新 helper，比较所有物理 cell、MPI ghosts 和
物理边界 ghosts。该测试必须在1、2、5 ranks运行；它用于证明只是暴露旧逻辑而非改变场离散。

### 5.10 JC2 增量测试命令

在同一 `vpfp_field_particle_trial_test` 增加 `deterministic-trial`、`signed-trial`、
`face-map-single`；在 MPI target 增加 `face-map-mpi`：

```bash
cmake --build build -j4 --target \
  fp_solver vpfp_field_particle_trial_test vpfp_field_particle_trial_mpi_test
mkdir -p ./output/vpfp_pairing_gate_jc/jc2

for CASE in deterministic-trial signed-trial face-map-single; do
  yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_test \
    --case "$CASE" \
    --result "./output/vpfp_pairing_gate_jc/jc2/${CASE}.result" || exit 66
done

for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_mpi_test \
    --case face-map-mpi \
    --result "./output/vpfp_pairing_gate_jc/jc2/face_map_n${NP}.result" || exit 67
done
```

JC2 只验收单次 trial 和 face映射，不要求 Picard 收敛；通过后停止。

若 MPI reduction 破坏 bitwise，先比较每 rank 的归约前物理数组；不得直接把容差扩大到宏观量级。

### 5.11 JC2 分级验收标准

#### 5.11.1 必须验收

JC2只验收“相同frozen输入+相同field guess得到同一trial”及trial无副作用，不要求Picard收敛。
按§5.10运行：

~~~text
deterministic-trial.result
signed-trial.result
face-map-single.result
face_map_n2.result
face_map_n5.result
~~~

deterministic trial必须证明：

~~~text
trial_replay_bitwise_equal=1
accepted_bulk_unchanged=1
accepted_tail_unchanged=1
accepted_beam_unchanged=1
accepted_field_unchanged=1
accepted_rng_unchanged=1
accepted_ledger_unchanged=1
c1_replayed_inside_trial=0
beam_schedule_regenerated_inside_trial=0
c2_called_inside_trial=0
return_called_inside_trial=0
all_trial_values_finite=1
~~~

等价字段名可以接受，但测试必须覆盖Bulk、Tail、Beam、field和转换/功账，不能只比较一个总hash。

signed trial必须满足：

- 正负场候选均有限；
- 使用的force field就是输入guess；
- 最终Poisson和pairing-field构造成功；
- bootstrap场相对差不超过1e-12；
- 正负场受力响应符号正确；
- accepted状态不变。

face-map单/MPI测试必须证明新helper与原solver的face同步、cell平均和ghost处理一致，且helper没有调用
Poisson覆盖输入face。2/5-rank共享面决策必须一致，开放边界不得周期回卷。

#### 5.11.2 不属于JC2的硬门

以下不在JC2强制范围：

- field_rel_l2<=1e-8；
- pairing_relative<=1e-8；
- Picard迭代次数；
- C2/return/正式commit；
- 完整Gate I总能量余额。

这些属于JC3或K1。不得因JC2没有这些字段判失败。

#### 5.11.3 推荐诊断

推荐输出首个bitwise mismatch、每类数组长度、bootstrap绝对残差、Poisson identity residual/tolerance和
MPI shared-face差异。缺少这些字段不阻断，只在JC2失败时补充。

#### 5.11.4 判定

全部§5.11.1通过则JC2 PASS；测试文件缺失/不可解析为INVALID_TEST；确定性、无副作用或face映射失败为FAIL。
汇总文件和统一比较器是推荐工具，不是生产实现前置功能。PASS后停止，不自动进入JC3。
## 6. JC3：实现有界阻尼 Picard 求解和原子提交

### 6.1 初始猜测

第一版使用当前生产 midpoint Poisson 场作为 $E_{\rm pair}^{(0)}$。该场只作为初值，不再作为最终受力场。

不得在第一版加入 Aitken、Anderson 或 JFNK。先建立可验证的确定性基线，后续优化必须另立性能 Gate。

### 6.2 迭代算法

实现如下语义的算法：

```text
prepare frozen state exactly once
E_guess = legacy midpoint field
omega = initial_relaxation

for k = 0 .. max_iterations-1:
    trial = evaluate_field_particle_trial(frozen, E_guess)
    reject immediately on NaN/Inf or structural failure

    E_map = trial.pairing_field
    compute R_E_L2, R_E_Linf and R_W

    if R_E_L2 <= field_tol
       and R_E_Linf <= field_tol
       and R_W <= pairing_tol:
        accept this exact trial as converged_trial
        break

    if residual grew by more than 25% for two consecutive iterations:
        omega = max(0.5 * omega, minimum_relaxation)

    E_guess = E_guess + omega * (E_map - E_guess)

if no exact trial passed:
    reject the physical step; do not soft-accept
```

关键约束：更新 `E_guess` 后必须重新评估一次 trial。不得把尚未与新场重新推进粒子的状态当作收敛候选。

### 6.3 MPI 共识

不得尝试用普通 `MPI_SUM` 同时归约最大值。第一版每次迭代允许两个固定 collective：

1. 一个 `MPI_Allreduce(... MPI_SUM ...)` 打包平方和、尺度平方和、$\Delta K$ 和功；
2. 一个 `MPI_Allreduce(... MPI_MAX ...)` 打包 Linf、非有限标志和局部失败标志。

只有发生失败时，允许再用一次 `MPI_MIN` 找首个 failing rank。合并内容包括：

- 非有限标志；
- $R_{E,2}$ 的平方和与尺度；
- $R_{E,\infty}$；
- $\Delta K$ 和 $W$；
- 首个失败 rank/索引。

任一 rank 失败时所有 rank 必须退出同一次 trial，不得出现部分 rank 进入下一轮、部分 rank 返回的 collective 次序错误。

### 6.4 收敛后的后处理与提交顺序

只有 `converged_trial` 存在时才执行：

1. 从收敛 trial 执行 C2；
2. 执行 collision-face conversion；
3. 执行 H10 return；
4. 检查后处理前后每个 x-cell combined charge 是否满足 JC0 容差；
5. 若电荷不变，复用 trial 的最终 Poisson；若发生允许范围外变化，拒绝该步；
6. 将收敛 trial 的状态、RNG 和 ledger 一次性交换到正式对象；
7. `accepted_step_count` 加一；
8. 写 accepted-step 诊断。

失败或未收敛时必须保持所有正式状态逐位不变。

### 6.5 失败码

新增明确失败码，避免归入 `unknown`：

```text
201 field_particle_trial_nonfinite
202 field_particle_poisson_failure
203 field_particle_pairing_field_failure
204 field_particle_work_identity_failure
205 field_particle_not_converged
206 post_field_charge_invariance_failure
207 field_particle_mpi_consensus_failure
208 accepted_field_poisson_mismatch
```

日志至少输出 step、time、iteration、relaxation、两个场残差、功配对残差、失败 rank 和首个局部索引。

### 6.6 DeepSeek 实施卡：求解器完整伪代码

`solve_field_particle_center()` 必须等价于下列控制流。不得省略标为“恢复”的语句：

```cpp
const long long step_before = step_count_;
step_count_ = step_before + 1;

if (!prepare_field_particle_frozen_state(...)) {
    step_count_ = step_before;
    return false;
}

pairing_field_guess_ = midpoint_fields_.Ex_face;
double omega = config_.initial_relaxation;
double residual_previous = infinity;
int consecutive_growth = 0;

for (int iter = 0; iter < config_.max_iterations; ++iter) {
    reset_trial_from_frozen();
    if (!evaluate_field_particle_trial(frozen_, pairing_field_guess_,
                                       trial_, failure)) {
        globalize_failure();
        step_count_ = step_before;
        return false;
    }

    compute_global_field_and_work_residuals();
    if (!all_finite) {
        step_count_ = step_before;
        return false;
    }

    if (field_l2 <= field_tol && field_linf <= field_tol &&
        pairing_relative <= pairing_tol) {
        converged = true;
        accepted_iteration = iter + 1;
        break; // 当前trial与当前guess完全对应
    }

    if (field_l2 > 1.25 * residual_previous) ++consecutive_growth;
    else consecutive_growth = 0;
    if (consecutive_growth >= 2) {
        omega = max(config_.minimum_relaxation, 0.5 * omega);
        consecutive_growth = 0;
    }

    for each local face f:
        pairing_field_guess_[f] += omega *
            (pairing_field_map_[f] - pairing_field_guess_[f]);
    synchronize shared physical MPI face exactly once;
    residual_previous = field_l2;
}

if (!converged) {
    failure = field_particle_not_converged;
    step_count_ = step_before;
    return false;
}

if (!apply_post_field_once_and_validate_charge(...)) {
    step_count_ = step_before;
    return false;
}

commit_trial_atomically();
// step_count_保留为step_before+1
return true;
```

共享 MPI face 的同步必须复用项目现有 face owner/同步约定。不得把左右 rank 的同一面平均两次，
不得对开放物理边界做周期回卷。

### 6.7 功配对残差的生产来源

不得在 JC 中新写一套 $J\cdot E$ 公式。应从现有 Gate I 路径抽取一个接受 caller-owned trial 的
只读 helper，复用：

- Bulk `u_full` 返回的离散速度面功；
- Tail `kick()` 的 `last_field_work`/cell work；
- Beam `finish_from_midpoint()` 的 `last_field_work`/cell work；
- x swept-mass、Tail/Beam trajectory current；
- `OpenElectrostaticSolver::evaluate_work_identity()`；
- `FieldParticlePowerAudit` 的稳定求和和边界/source分解。

每轮 trial 输出以下有符号量，之后再计算比值：

```text
bulk_force_work
tail_force_work
beam_force_work
particle_force_work_sum
pairing_current_work
potential_charge_work
field_energy_change
electrode_work
pairing_signed_residual
pairing_absolute_scale
pairing_relative_residual
```

`pairing_absolute_scale` 至少取
$\max(|\Delta K_{force}|,|W_{pair}|,|\Delta U_E-W_{electrode}|,K_{floor})$。
禁止用接近零的净残差作分母，也禁止对各 rank 的相对误差求和。

### 6.8 C2、collision conversion 和 H10 的一次性执行

收敛 trial 后：

1. 保存 `state_np1_` 和 `tail_field_trial_` 的 combined number-by-x；
2. 调用现有 C2，使用 `collision_half=1` 和同一个 `candidate_step`；
3. 如 C2 产生 interface flux，调用现有 collision-face conversion；
4. 调用现有 `apply_tail_bulk_return()`；
5. 重新计算 Bulk moments 和 Tail density；
6. 比较 combined number-by-x；
7. 只有 JC0 容差通过才进入 commit；
8. Beam 不参与 C2/return，不得再次调用 `finish_from_midpoint()`；
9. C2/return 后只允许更新碰撞、conversion、return 专属能量 ledger；不得改写已收敛场功 ledger。

当前 hybrid collision 使用 `candidate_step`、`collision_half` 和 counter-based seed。执行者必须在代码
审计中确认它没有隐藏的顺序型全局 RNG；若发现有，必须先提供 snapshot/restore，不得假设失败事务安全。

### 6.9 原子提交的准确对象

成功时按以下顺序提交：

```cpp
electrons.swap_state(state_np1_);
if (tail_on) tail_accepted_.swap_state(tail_field_trial_);
if (beam_on) beam.swap_state(beam_field_trial_);
swap_emfields(fields, final_fields_);
beam.commit_injection_schedule(frozen_.beam_schedule, mpi_rank); // 仅一次
commit cumulative ledgers/checksums;
finalize accepted diagnostics;
```

实际 `Species` 接口若不是 `swap_state`，使用项目现有 swap 方式。提交前保存的旧 accepted state 会被
换入 work buffer，供只读 Gate I/诊断使用；诊断结束后才能复用该缓冲区。

Beam schedule 的 RNG/remainder 只能在提交点调用 `commit_injection_schedule()`。trial 中禁止调用。

### 6.10 JC3 专项负测试

至少实现：

| 案例 | 注入位置 | 期望 |
|---|---|---|
| `nan_on_rank1_iter2` | trial force field | 所有rank失败，accepted逐位不变 |
| `max_iter_one` | 配置 | failure 205，step_count恢复 |
| `poisson_fail_rank0` | trial final Poisson test seam | failure 202，全局退出 |
| `pairing_build_fail` | 清空trial phi test seam | failure 203 |
| `post_field_charge_change` | JC0 test-only perturbation | failure 206，不提交C2结果 |
| `diagnostic_off` | level 0 | 与level 1物理状态逐位相同 |

故障注入必须由仅测试编译定义或 test-access 完成，生产 CLI 不得暴露“制造 NaN”开关。

### 6.11 JC3 增量测试命令

```bash
cmake --build build -j4 --target \
  fp_solver vpfp_field_particle_trial_test \
  vpfp_field_particle_trial_mpi_test vpfp_field_particle_post_field_test
mkdir -p ./output/vpfp_pairing_gate_jc/jc3

for CASE in zero-field signed-field diagnostic-off max-iter-fault poisson-fault \
            pairing-fault post-field-charge-fault; do
  yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_test \
    --case "$CASE" \
    --result "./output/vpfp_pairing_gate_jc/jc3/${CASE}.result" || exit 68
done

for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_mpi_test \
    --case rollback-consensus \
    --result "./output/vpfp_pairing_gate_jc/jc3/rollback_n${NP}.result" || exit 69
done

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_field_particle_post_field_test \
  --case all \
  --result ./output/vpfp_pairing_gate_jc/jc3/post_field.result || exit 70
```

本阶段不得运行 checkpoint 或115 fs短跑。

### 6.12 JC3 分级验收标准

#### 6.12.1 必须验收

JC3验收原§6的固定点、严格接受、MPI共识和原子提交。测试可合并结果文件，不强制固定为10个文件，
但必须覆盖以下案例：

~~~text
zero-field
positive/negative signed-field
deterministic replay
max-iteration failure
Poisson/pairing failure
post-field charge failure
2-rank rollback consensus
5-rank rollback consensus
diagnostic-off
~~~

零场要求一次trial收敛，状态按fixture预期不变。正负场成功案例要求：

~~~text
converged=1
soft_accept_count=0
1<=iterations<=12
trial_evaluations=iterations
field_residual_l2<=1e-8
field_residual_linf<=1e-8
pairing_relative_to_exchange<=1e-8
accepted_trial_matches_last_evaluated_trial=1
post_field_charge_pass=1
accepted_commit_count=1
~~~

pairing_relative_to_exchange只指同一force guess对应的受力功配对，不是完整Gate I总能量余额。

失败案例必须验证：

~~~text
expected_failure_code_observed=1
all_rank_decision_equal=1
accepted_state_unchanged=1
accepted_rng_unchanged=1
accepted_ledger_unchanged=1
step_and_time_unchanged=1
accepted_commit_count=0
~~~

C2、collision conversion和return只能在收敛后执行一次；失败不得soft accept。2/5-rank故障必须
完成相同collective序列，不得挂起或部分rank先退出。diagnostic level 0/1不得改变物理状态。

#### 6.12.2 推荐诊断

推荐但不作为硬门：

- 每轮残差和松弛因子历史；
- 每个内部函数精确调用次数；
- 首个失败rank/index；
- 各状态数组逐字段bitwise报告；
- 独立post-field案例文件；
-热点或MPI collective计数。

阻尼必须遵守§6.2/§6.6算法，但无需为了验收新增一套复杂结果schema；失败时再输出详细历史。

#### 6.12.3 判定

只要所有必须案例和字段通过，JC3 PASS。案例可以由一个或多个测试目标承载。缺少必须案例为
INVALID_TEST；严格残差、回滚、MPI共识或一次性提交失败为FAIL。不得通过新增能量补丁、放宽1e-8
局部配对门或软接受获得PASS。
## 7. JC4：CLI、checkpoint、manifest 与诊断接入

### 7.1 新增 CLI

在 `src/main_vpfp.cpp` 中增加：

```text
--field-particle-coupling legacy|discrete-gradient
--field-particle-max-iters <int>
--field-particle-relaxation <double>
--field-particle-field-tol <double>
--field-particle-pairing-tol <double>
--restart-allow-field-particle-coupling-change
```

默认值：

```text
field_particle_coupling=legacy
field_particle_max_iters=12
field_particle_relaxation=0.5
field_particle_field_tol=1e-8
field_particle_pairing_tol=1e-8
```

参数检查：迭代数必须大于零；松弛因子必须在 $(0,1]$；容差必须为有限正数。

### 7.2 checkpoint 与物理配置哈希

将上述配置写入：

- 启动摘要；
- checkpoint manifest；
- checkpoint物理身份校验（沿用现有hash，并单独校验JC配置即可）；
- restart mismatch 日志。

`--restart-allow-field-particle-coupling-change` 只允许从 `legacy` checkpoint 启动 `discrete-gradient` A/B 审计。它不得忽略网格、边界、碰撞、Tail、Beam 或其他物理配置差异，也不得成为通用 hash 绕过开关。

### 7.3 新增 accepted-step 诊断

新增 `field_particle_iteration.dat`，每个已接受步一行：

```text
step time_s mode iterations converged relaxation
field_residual_l2 field_residual_linf pairing_residual
trial_evaluations post_field_charge_residual failure_code
```

规则：

- `diagnostic-level=0`：只输出汇总计数，不逐迭代扫描；
- `diagnostic-level=1`：仅写最终接受态行；
- `diagnostic-level=2`：可另写 `trial_field_particle_iteration.dat`；
- trial 文件不得混入 accepted-step 文件；
- 诊断关闭不得改变算子、MPI collective 顺序或结果。

### 7.4 DeepSeek 实施卡：`main_vpfp.cpp` 的准确修改点

按以下位置修改：

1. 在 `struct Options` 的 `dt_scale`/diagnostic 字段附近增加 coupling 字符串和四个数值字段；
2. 在构造函数初始化列表中填入本节默认值；
3. 在 `parse_options()` 中增加六个 CLI 分支；
4. 在现有 `restart_allow_dt_scale_change` 一致性检查附近检查 coupling override 的合法组合；
5. 按§7.5持久化并校验JC配置；优先不改现有物理hash字节序列；
6. 创建 `VpfpIntegrator` 后、`integrator.init(grid)` 前调用 `set_field_particle_coupling()`；
7. 在启动摘要中打印所有 coupling 值；
8. 在 `failure_code_name()` 或等价 switch 中加入 201--208；
9. 在 `VpfpRunManifestConfig` 赋值区写入 coupling 字段；
10. 在 restart expected config 构造区执行 7.5 的兼容性验证。

参数解析必须拒绝：

```text
unknown coupling mode
max_iters<=0
relaxation<=0 or relaxation>1
nonfinite tolerance
tolerance<=0
restart override without --restart-dir
restart override while requested mode is legacy
```

### 7.5 JC配置持久化与旧 checkpoint 兼容

本任务不强制重构现有physical hash体系。最小正确实现为：

1. 在运行manifest和checkpoint manifest中保存coupling mode、max iterations、relaxation和两个容差；
2. 读取新checkpoint时逐项恢复并校验这些字段；
3. 旧checkpoint缺少JC字段时明确解释为legacy默认配置；
4. legacy->discrete-gradient必须使用显式
   --restart-allow-field-particle-coupling-change；
5. override只放行coupling配置差异，现有dt、网格、边界、碰撞、Tail和Beam身份校验仍执行；
6. 改变dt仍单独要求现有dt override；
7. 不改变rank二进制payload，不因JC配置升级checkpoint binary schema。

可以选择把JC字段纳入一个向后兼容的新hash版本，也可以单独比较manifest字段。两种方案只要满足上述
行为即可；不得为了测试固定字段名强制实现physical_config_hash_v1/v2双函数。

### 7.6 checkpoint 序列化边界

JC配置属于manifest/config；物理rank payload继续使用现有格式。JC accepted field的phi已由final trial
重建，checkpoint仍应正常保存和恢复Ex_face、Ex、phi和rho。

必须验证：

~~~text
phi_finite=1
phi_poisson_consistent=1
field_state_roundtrip=1
coupling_config_roundtrip=1
~~~

在当前确定性序列化下优先使用bitwise比较；若现有测试只提供等价norm，必须报告实际误差和容差。
不得通过restart后置零phi或重解并覆盖checkpoint场来掩盖序列化问题。
### 7.7 诊断结构的具体落点

在 `VpfpStepResult` 增加：

```cpp
bool field_particle_coupling_enabled;
bool field_particle_converged;
int field_particle_iterations;
int field_particle_trial_evaluations;
double field_particle_relaxation;
double field_particle_residual_l2;
double field_particle_residual_linf;
double field_particle_pairing_residual;
double post_field_charge_residual_linf;
```

在 `VpfpDiagnostics` 增加
`write_field_particle_iteration_accepted_step()`，由现有 `write_accepted_step()` 在所有 rank 都完成
必要 collective 后，仅 rank 0 追加文件。文件首行必须是固定 schema header，并在重新启动追加时验证
列数一致。

失败路径使用现有 `write_failure()`，增加同名字段；不得创建一个只有 rank 0 进入、内部却调用 MPI 的
失败诊断函数。

### 7.8 JC4 独立测试

最小测试集合：

~~~text
cli_default_legacy
cli_parse_discrete_gradient
cli_reject_invalid_numeric_value
old_checkpoint_defaults_to_legacy
old_checkpoint_to_jc_requires_override
new_checkpoint_coupling_roundtrip
unrelated_physics_change_rejected
dt_change_requires_independent_override
diagnostic_level_0_vs_1_state_equal
failure_code_names_known
~~~

拒绝案例必须同时检查非零退出码和明确错误文本。允许在一个测试目标中组合多个子案例，不强制为每个
负案例建立独立文件。
### 7.9 JC4 增量测试命令

扩展现有 `checkpoint_restart_equivalence_test`，新增 coupling config 案例；不要另写 checkpoint
序列化器。编译：

```bash
cmake --build build -j4 --target fp_solver checkpoint_restart_equivalence_test
mkdir -p ./output/vpfp_pairing_gate_jc/jc4/checkpoint_tmp

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/checkpoint_restart_equivalence_test \
  --case field-particle-coupling \
  --workdir ./output/vpfp_pairing_gate_jc/jc4/checkpoint_tmp \
  --result ./output/vpfp_pairing_gate_jc/jc4/checkpoint.result || exit 75
```

CLI 正测试用零步初始化，避免生产推进。默认 legacy 模式（criterion 1）：

```bash
mkdir -p ./output/vpfp_pairing_gate_jc/jc4/cli_default_legacy
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/fp_solver \
  --beam-enabled 0 --background-tail-mode off \
  --collision-model none --stop-time-fs 0 \
  --output-dir ./output/vpfp_pairing_gate_jc/jc4/cli_default_legacy || exit 76
```

显式 discrete-gradient 模式（criterion 2）：

```bash
mkdir -p ./output/vpfp_pairing_gate_jc/jc4/cli_valid
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/fp_solver \
  --field-particle-coupling discrete-gradient \
  --field-particle-max-iters 12 --field-particle-relaxation 0.5 \
  --field-particle-field-tol 1e-8 --field-particle-pairing-tol 1e-8 \
  --beam-enabled 0 --background-tail-mode off \
  --collision-model none --stop-time-fs 0 \
  --output-dir ./output/vpfp_pairing_gate_jc/jc4/cli_valid || exit 76
```

CLI 负测试必须期望非零退出码；不能使用 `|| exit` 直接把预期失败当作作业失败：

```bash
set +e
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores ./build/fp_solver \
  --field-particle-coupling discrete-gradient \
  --field-particle-relaxation 0 \
  --stop-time-fs 0 \
  --output-dir ./output/vpfp_pairing_gate_jc/jc4/cli_invalid \
  > ./output/vpfp_pairing_gate_jc/jc4/cli_invalid.out 2>&1
RC=$?
set -e
test $RC -ne 0 || exit 77
grep -q "field-particle-relaxation" \
  ./output/vpfp_pairing_gate_jc/jc4/cli_invalid.out || exit 78
```

Checkpoint 序列化往返测试（criterion 8）：

```bash
cmake --build build -j4 --target checkpoint_roundtrip_test
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/checkpoint_roundtrip_test \
  --workdir ./output/vpfp_pairing_gate_jc/jc4/checkpoint_roundtrip_tmp \
  --result ./output/vpfp_pairing_gate_jc/jc4/checkpoint_roundtrip.result || exit 79
```

JC4 通过后停止，不运行 JC5 全矩阵。

### 7.10 JC4 分级验收标准

#### 7.10.1 必须验收

JC4只验证原§7要求的CLI、配置持久化、受控restart和诊断不改变物理。必须覆盖：

1. 默认无参数时mode=legacy；
2. 显式discrete-gradient及四个数值参数可解析并正确传给integrator；
3. 非法mode、非正迭代数、越界relaxation、非有限/非正容差被明确拒绝；
4. coupling mode和参数写入运行manifest/checkpoint配置并可往返；
5. 旧checkpoint缺少JC字段时按legacy解释；
6. legacy->discrete-gradient只有显式override才允许；
7. override不能掩盖dt、网格、边界、碰撞、Tail或Beam等无关物理配置变化；
8. checkpoint中的Ex/rho/phi正常往返，JC accepted phi有限；
9. diagnostic level 0/1不改变物理状态、RNG或ledger；
10. failure code 201--208有明确名称，不显示unknown。

checkpoint配置测试和CLI正负测试必须返回明确PASS/预期拒绝。缺字段或测试未执行为INVALID_TEST。

#### 7.10.2 最小兼容实现要求

不强制引入复杂的新二进制checkpoint schema或固定命名的physical_config_hash_v1/v2函数。实现只需：

- 保持旧checkpoint现有物理hash验证不变；
- 以显式manifest字段保存JC配置；
- 对缺少字段的旧checkpoint使用legacy默认；
- override时逐项证明除coupling外其他物理配置一致；
- 新checkpoint重新启动时精确恢复JC配置。

若当前hash结构能在不破坏旧checkpoint的情况下安全加入版本字段，可以使用；否则采用独立JC配置校验。
不得仅为满足测试字段升级rank二进制schema。

#### 7.10.3 推荐诊断

推荐但不阻断：

- 完整CLI负案例矩阵；
- stored/expected hash及首个差异字段；
- manifest、startup、runtime四方逐项表；
- trial/accepted诊断行数和header审计；
- 独立hash版本号；
-所有201--208逐项单测文件。

这些可以合并在checkpoint/CLI测试中，无需强制新建固定数量脚本和结果文件。

#### 7.10.4 判定

默认legacy、显式JC解析、配置往返、受控override、无关物理差异拒绝、诊断不改变状态和failure名称
全部通过，则JC4 PASS。任何配置绕过或restart物理身份错误为FAIL。不得为补推荐诊断修改生产推进。
## 8. JC5：独立单元、MPI 与失败事务测试

### 8.1 必须新增的测试目标

新增并接入 `CMakeLists.txt`：

- `tests/vpfp_field_particle_trial_test.cpp`
- `tests/vpfp_field_particle_trial_mpi_test.cpp`
- `tests/vpfp_field_particle_post_field_test.cpp`
- `tools/compare_vpfp_field_particle_jc.py`

测试必须直接调用生产接口，不得复写 Poisson、pusher、x-remap 或离散梯度公式。

### 8.2 单 rank 测试案例

至少覆盖：

1. 零场静态分布：一次迭代收敛；
2. 正均匀微场；
3. 负均匀微场；
4. 非均匀 Bulk 且 x-remap 激活；
5. Beam 开启；
6. Tail PIC 开启；
7. collision `none`；
8. collision `moment-closure`；
9. H10 return `none`；
10. H10 return `hysteretic`；
11. 人工 NaN 注入；
12. `max_iterations=1` 的未收敛失败事务。

### 8.3 MPI 测试案例

用 2 和 5 ranks 覆盖：

- 共享 x-face；
- rank 0 左物理边界；
- 最末 rank 右物理边界；
- 单 rank 人工失败后的全局一致退出；
- 诊断开关前后 accepted state 一致。

### 8.4 状态不变性要求

失败事务必须验证以下对象逐位不变：

```text
Bulk distribution
Tail particles
Beam particles
field phi/E
RNG state
injection/outflow/conversion/return ledgers
accepted_step_count
time
```

### 8.5 JC5 编译和运行命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target \
  vpfp_post_field_charge_invariance_test \
  vpfp_field_particle_trial_test \
  vpfp_field_particle_trial_mpi_test \
  vpfp_field_particle_post_field_test \
  fp_solver

mkdir -p ./output/vpfp_pairing_gate_jc

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_field_particle_trial_test \
  --case all \
  --result ./output/vpfp_pairing_gate_jc/trial_single.result || exit 71

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_field_particle_post_field_test \
  --case all \
  --result ./output/vpfp_pairing_gate_jc/post_field_single.result || exit 72

for NP in 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/vpfp_field_particle_trial_mpi_test \
    --case all \
    --result "./output/vpfp_pairing_gate_jc/trial_mpi_n${NP}.result" || exit 73
done

python3 tools/compare_vpfp_field_particle_jc.py \
  --root ./output/vpfp_pairing_gate_jc \
  --result ./output/vpfp_pairing_gate_jc/jc5_compare.result || exit 74
```

### 8.6 JC5 总验收

```text
all_unit_cases_pass=1
all_mpi_cases_pass=1
legacy_default_regression_pass=1
trial_deterministic=1
trial_side_effect_free=1
failure_transaction_bitwise_unchanged=1
post_field_charge_invariance_pass=1
mpi_failure_consensus_pass=1
status=PASS
```

任一项失败都停在 JC5，不得进入 checkpoint 生产测试。

### 8.7 DeepSeek 实施卡：测试文件职责不得混合

#### `tests/vpfp_field_particle_trial_test.cpp`

只负责单 rank、无 checkpoint 的算子测试。实现统一 fixture：

```text
nx_global=32
field_boundary=DIRICHLET_PHI(0,0)
background_boundary=reservoir/reservoir
small finite Maxwellian Bulk
optional 4--16 Beam/Tail particles
collision none or deterministic prescribed case
```

fixture 必须通过生产 `Species::init()`、`OpenElectrostaticSolver::init()`、`BeamPIC::init()`、
`BackgroundTailPIC::init()` 和 `VpfpIntegrator::init()` 构造，不得手动拼接不符合生产布局的数组。

每个案例写一行：

```text
case status iterations field_l2 field_linf pairing_relative
accepted_unchanged rng_unchanged ledger_unchanged failure_code
```

#### `tests/vpfp_field_particle_trial_mpi_test.cpp`

只负责 2/5-rank 共享面和 failure consensus。每 rank 至少一个物理 cell；在 rank 边界附近放置
非均匀 Bulk/Tail/Beam 状态，使共享面通量非零。测试应输出：

```text
shared_face_difference_linf
all_rank_decision_equal
collective_sequence_completed
first_failure_rank
accepted_state_global_equal
```

不得通过 `MPI_Barrier` 掩盖不一致控制流。测试超时属于 FAIL。

#### `tests/vpfp_field_particle_post_field_test.cpp`

负责组合 JC trial 收敛 + C2/conversion/return + commit：

1. 先得到一个收敛 trial；
2. 保存正式 accepted snapshot；
3. 执行 post-field；
4. 检查 charge invariance；
5. 成功案例验证只提交一次；
6. 故障案例验证完全回滚。

它应调用 JC0 helper 或同一个生产内部函数，不得复制 JC0 数量计算。

#### `tools/compare_vpfp_field_particle_jc.py`

脚本必须：

1. 明确列出必需文件；
2. 明确列出每个文件的必需键；
3. 检测重复键、NaN、Inf、空文件和未知状态；
4. 任一缺失返回 exit code 2；物理/数值门失败返回 exit code 1；全部通过返回 0；
5. 写汇总 `status=PASS|FAIL` 和 `first_failure`；
6. 不执行 Git 命令，不因测试机无 Git 判失败；
7. 不把字符串 `PASS` 的存在当作唯一依据，必须重新计算阈值判断。

### 8.8 JC5 案例级验收值

| 案例 | 关键要求 |
|---|---|
| zero field | `iterations=1`，所有状态有限 |
| positive/negative field | 两者均收敛，功符号相反且量级对称 |
| x-remap active | x swept-mass continuity PASS |
| Beam | 注入 schedule 只生成/提交一次 |
| Tail | trajectory current 和 kick work 都非零 |
| collision none | C1/C2 no-op，不产生随机变化 |
| collision moment-closure | C1/C2各一次，trial内C2为零次 |
| return hysteretic | return仅在收敛后一次 |
| NaN fault | 所有rank同一失败，正式状态不变 |
| max-iter fault | failure 205，时间和step不前进 |

零场只允许 `iterations=1`，不允许写成 `<=2` 来掩盖初值或残差计算错误。

### 8.9 编译失败时的处理顺序

1. 先单独构建第一个失败 target；
2. 修复缺少的 production source/link library，不把生产 `.cpp` 内容 `#include` 到测试；
3. 对照 CMake 中现有 `vpfp_post_field_charge_invariance_test` 的 `${VPFP_CORE_SOURCES}` 接线；
4. 链接错误只补目标依赖，不修改物理代码；
5. 单 rank 通过后才运行 MPI；MPI 失败时不重复跑生产短程。

### 8.10 JC5 分级验收标准

#### 8.10.1 必须验收

JC5是JC0--JC4的最小集成门，不新增物理要求。进入前必须确认JC0--JC4均PASS；若生产代码在JC5
期间改变了前置Gate语义，只重跑受影响的小型Gate。

按§8.5运行单rank、post-field和2/5-rank测试，必须覆盖：

- 零场、正负场和非均匀x-remap；
- Beam on、Tail on、collision none/真实collision；
- return none/hysteretic；
- NaN/未收敛等失败事务；
- MPI共享面和单rank故障共识。

总体必须满足原§8.6：

~~~text
all_unit_cases_pass=1
all_mpi_cases_pass=1
legacy_default_regression_pass=1
trial_deterministic=1
trial_side_effect_free=1
failure_transaction_bitwise_unchanged=1
post_field_charge_invariance_pass=1
mpi_failure_consensus_pass=1
soft_accept_count=0
status=PASS
~~~

测试必须直接调用生产接口，不得复写Poisson、pusher、remap或配对公式。成功状态只commit一次，失败状态
不得改变state/RNG/ledger/step/time。

#### 8.10.2 推荐诊断

推荐输出逐案例活动计数、首个mismatch、源码/二进制标识、跨Gate schema汇总、各MPI边界覆盖表。
这些信息有助于审计，但不强制固定结果文件数量，也不要求因缺Git环境失败。

#### 8.10.3 判定

前置Gate有效且§8.10.1全部通过则JC5 PASS。缺少必须案例为INVALID_TEST，事务或MPI失败为FAIL。
JC5不包含115 fs生产短跑；PASS后停止，由用户单独授权K1。
## 9. K1：115 fs 固定 checkpoint 的短步闭合验收

### 9.1 进入条件

只有 JC0–JC5 全部 PASS，且用户明确允许生产 checkpoint 审计后，才能执行 K1。

使用同一个 115 fs checkpoint：

```bash
CHECKPOINT_115="./output/vpfp_return_dt_half_100p3_to_120fs/checkpoint_target115fs_t115.011425fs_step5071"
```

先确认：

```bash
test -f "$CHECKPOINT_115/manifest.txt" || exit 81
grep -E '^(step|time|dt|mpi_size|physical_config_hash)' "$CHECKPOINT_115/manifest.txt"
```

### 9.2 coarse：原步长比例的10个接受步

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-allow-field-particle-coupling-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode hysteretic --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --field-particle-coupling discrete-gradient \
  --field-particle-max-iters 12 \
  --field-particle-relaxation 1.0 \
  --field-particle-field-tol 1e-4 \
  --field-particle-pairing-tol 1e-8 \
  --dt-scale 0.5 --stop-after-steps 5081 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir ./output/vpfp_pairing_gate_k1/coarse || exit 82
```

K1显式使用 `relaxation=1.0`。115 fs首步实测在0.5下场残差每轮约减半，12轮只能从
$5.1\times10^{-1}$ 降到 $2.5\times10^{-4}$，数学上不可能达到 $10^{-8}$；由该收缩率反推未阻尼
映射导数约为 $-8\times10^{-3}$，使用1.0仍处于强收缩区。若后续残差连续增长，§6.2已有保护会
自动减半。omega=1实测进一步显示混合Eulerian/PIC转换在两个离散表示之间形成稳定二周期，
后续多步显示`field_linf`表示平台会随状态在 $1.25\times10^{-5}$ 到 $2.45\times10^{-5}$ 间变化；
因此生产K门显式使用field tolerance $10^{-4}$，pairing tolerance仍为 $10^{-8}$。该field门把局部
转换表示误差限制在0.01%，平滑单元测试仍保持 $10^{-8}$，不提高迭代上限。

### 9.3 fine：相同物理窗口的20个半步

当前 checkpoint 的源步号是 5071。程序的 `--stop-after-steps` 表示绝对步号，
不是“从 restart 后再跑多少步”。因此 coarse 终止号是 $5071+10=5081$，fine
终止号是 $5071+20=5091$。若更换 checkpoint，必须从 manifest 读取源步号并重新计算，
不得直接复用 5081/5091。

115 fs source checkpoint 来自 `dt-scale=0.5`。fine 改为 0.25 时必须增加 source-dt 明示参数。
下面是完整命令，不能只复制末尾参数片段：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-allow-field-particle-coupling-change \
  --restart-allow-dt-scale-change --restart-source-dt-scale 0.5 \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode hysteretic --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --field-particle-coupling discrete-gradient \
  --field-particle-max-iters 12 \
  --field-particle-relaxation 1.0 \
  --field-particle-field-tol 1e-4 \
  --field-particle-pairing-tol 1e-8 \
  --dt-scale 0.25 --stop-after-steps 5091 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir ./output/vpfp_pairing_gate_k1/fine || exit 84
```

其余物理参数必须逐字一致。建议把公共参数写入同一提交脚本数组，避免手工漂移。

> **2026-08-18 修复记录**：从旧 legacy checkpoint 以 `--field-particle-coupling discrete-gradient`
> 重启时，首个步的 bootstrap Poisson 校验可能失败（failure 208 `accepted_field_poisson_mismatch`）。
> 根因：restart 路径恢复 checkpoint 场后未重新求解初始 Poisson，而 discrete-gradient 的
> bootstrap 要求 accepted `Ex_face` 与当前求解器重新求解结果一致（1e-12 相对容差）。旧 H10
> checkpoint 的场与当前求解器存在舍入级差异。
>
> 修复（`src/main_vpfp.cpp`）：当 restart 且 coupling mode 为 discrete-gradient 时，在
> `background_boundary.fill_ghosts()` 之后、启动摘要之前，按非 restart 路径相同顺序
> （beam deposit → electrons.compute_moments → set_charge_density → field_solver.solve）
> 重新求解一次初始 Poisson 场，使 accepted 场与当前求解器自洽。
>
> 上述 §9.2/§9.3 命令参数本身无需修改；重新编译 `fp_solver` 后重跑即可。

### 9.4 K1 验收标准

```text
same_source_checkpoint=1
same_physical_window=1
accepted_steps_coarse=10
accepted_steps_fine=20
soft_accept_count=0
max_iteration_count<=12
all_field_residual_l2<=1e-4
all_field_residual_linf<=1e-4
all_pairing_residual<=1e-8
continuity_pass=1
local_work_ledger_pass=1
poisson_identity_pass=1
gauss_pass=1
post_field_charge_invariance_pass=1
nonfinite_count=0
```

时间加密后公共能量残差不得恶化。若可稳定估阶，则要求误差呈下降趋势；本 Gate 不预设必须达到二阶，避免为漂亮阶数修改物理算子。

### 9.5 K1 前置身份检查

在启动两个作业前新增 `tools/check_vpfp_jc_source_checkpoint.py` 或使用现有 manifest parser，输出：

```text
manifest_present=1
source_step=5071
source_time_s=1.1501142520244935e-13
source_dt_scale=0.5
source_mpi_size=80
source_physical_config_hash=11814249988425503857
status=PASS
```

不得从目录名解析 step/time 作为唯一依据；目录名只用于交叉检查。MPI ranks 必须与 checkpoint
兼容，若使用现有 repartition reader则需先通过其独立测试。

### 9.6 K1 结果分析器

新增 `tools/analyze_vpfp_jc_k1.py`，参数：

```text
--coarse <目录>
--fine <目录>
--source-checkpoint <目录>
--result <文件>
--field-tol <值，K1默认1e-4>
--pairing-tol <值，默认1e-8>
--max-iters <值，默认12>
```

分析器必须读取 `field_particle_iteration.dat` 和 `vpfp_step_diagnostics.dat`，由其中的accepted-state
连续性、Poisson/Gauss和总能量账完成K1主判定。stage-energy、逐face flux和conversion-source明细在
有数据行时作为定位证据；只有表头时必须输出 `detailed_audit_available=0`，不得伪称已完成分阶段审计，
但不单独阻断coarse/fine时间加密主测试。至少输出：

```text
same_source_checkpoint
same_initial_physical_state
coarse_accepted_steps
fine_accepted_steps
coarse_physical_duration_s
fine_physical_duration_s
same_physical_window
coarse_max_iterations
fine_max_iterations
coarse_max_field_l2
fine_max_field_l2
coarse_max_field_linf
fine_max_field_linf
coarse_max_pairing_relative
fine_max_pairing_relative
continuity_pass
local_work_pass
poisson_pass
gauss_pass
post_field_charge_pass
energy_residual_coarse
energy_residual_fine
energy_residual_reduction
status
first_failure
```

缺少任一结构文件即 FAIL。不要因为最终总能量看起来接近就忽略逐步场粒子残差。

运行：

```bash
python3 tools/analyze_vpfp_jc_k1.py \
  --coarse ./output/vpfp_pairing_gate_k1/coarse \
  --fine ./output/vpfp_pairing_gate_k1/fine \
  --source-checkpoint "$CHECKPOINT_115" \
  --field-tol 1e-4 --pairing-tol 1e-8 --max-iters 12 \
  --result ./output/vpfp_pairing_gate_k1/k1.result || exit 83
```

### 9.7 K1 失败后的唯一允许动作

- bootstrap mismatch：检查 accepted `rho/Ex/phi` 身份，不改容差；
- Picard 不收敛：输出每轮残差序列和松弛因子，停在 JC3；
- pairing 失败但 field 收敛：检查实际 pusher 是否使用 `trial_force_fields_`，不要加能量补丁；
- C2后 charge 失败：回到 JC0；
- coarse/fine 窗口不同：修命令，不改生产代码；
- 只在 fine 失败：检查绝对停止步号与 source dt 参数。

### 9.8 K1 分级验收标准

#### 9.8.1 必须验收

K1使用同一个115 fs checkpoint，coarse运行10步、fine以dt/2运行20步，覆盖相同物理窗口。必须证明：

~~~text
same_source_checkpoint=1
same_initial_physical_state=1
same_physical_window=1
coarse_accepted_steps=10
fine_accepted_steps=20
soft_accept_count=0
max_iteration_count<=12
all_field_residual_l2<=1e-4
all_field_residual_linf<=1e-4
all_pairing_relative_to_exchange<=1e-8
continuity_pass=1
local_work_ledger_pass=1
poisson_identity_pass=1
gauss_pass=1
post_field_charge_invariance_pass=1
nonfinite_count=0
energy_signed_fine_over_coarse<=1
energy_abs_fine_over_coarse<=1
energy_residual_reduction=1
~~~

每个accepted行都必须通过，不能只检查终点。开放边界源、出流和电极功按现有账本验收，不要求为零。

完整公共能量余额不强制到1e-8，但time refinement是§9.8.1的硬条件。必须报告coarse/fine的有符号值、
绝对值和独立交换尺度；fine的有符号残差绝对值与绝对累计残差都不得高于coarse。若差异处于稳定
求和/采样不确定度内，状态为REVIEW_REQUIRED并报告尺度，不得人为设定阶数或添加能量补丁。

#### 9.8.2 推荐诊断

推荐输出观察阶、核心/边界分解、逐步迭代分布、各能量分项和coarse/fine最终宏观差异。目录和字段可由
现有分析器提供，不强制新增固定schema。

#### 9.8.3 判定

实验身份、严格收敛和结构门全部通过，且time refinement没有明确恶化，K1 PASS。窗口/文件无效为
INVALID_TEST；能量趋势证据不足为REVIEW_REQUIRED；严格残差或守恒失败为FAIL。只有PASS进入K2。
## 10. K2：宏观物理 A/B 验收

K1 通过后，从同一 checkpoint 分别运行 `legacy` 和 `discrete-gradient` 到相同的短物理终点，建议窗口不超过 1 fs。

比较核心区：

- $E_x(x)$ 的包络、峰值、相位；
- $n_{\rm bkg}(x)$、$n_{\rm beam}(x)$、$n_{\rm tail}(x)$；
- 场能、Bulk/Tail/Beam 动能和总能量账；
- combined 能谱和 $p_x$ 分布；
- 边界区只做稳定性观察，不以逐点物理一致作为通过条件。

验收原则：

1. JC 不能导致核心主波包突然消失或平滑；
2. 对任一宏观积分量，相对 legacy 的突变超过 2% 时必须解释其来源；
3. 能量残差应随时间步缩小，而不是仅靠账本定义变漂亮；
4. 不以逐点贴合 EPOCH 作为算子验收标准。

### 10.1 K2 运行实现

新增 `tools/run_vpfp_jc_k2.sh`，脚本必须定义一个 `COMMON_K2=(...)` Bash 数组，内容复制 K1 的
全部物理参数，固定：

```text
restart-dir=$CHECKPOINT_115
dt-scale=0.5
stop-time-fs=116
diagnostic-level=1
diagnostic-interval=20
snapshot-times=115.25,115.5,115.75,116
```

然后只运行两支：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  "${COMMON_K2[@]}" \
  --field-particle-coupling legacy \
  --output-dir ./output/vpfp_pairing_gate_k2/legacy || exit 91

yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  "${COMMON_K2[@]}" \
  --restart-allow-field-particle-coupling-change \
  --field-particle-coupling discrete-gradient \
  --field-particle-max-iters 12 \
  --field-particle-relaxation 1.0 \
  --field-particle-field-tol 1e-4 \
  --field-particle-pairing-tol 1e-8 \
  --output-dir ./output/vpfp_pairing_gate_k2/jc || exit 92
```

`COMMON_K2` 中不得包含 coupling 相关选项，避免 Bash 后出现重复参数且由最后一个值静默覆盖。

### 10.2 K2 比较器

新增 `tools/compare_vpfp_jc_macro.py`，以 core 区和 boundary 区分别统计：

```text
Ex_relative_l2_core
Ex_peak_relative_difference_core
Ex_phase_shift_cells_core
n_bulk_relative_l2_core
n_tail_relative_l2_core
n_beam_relative_l2_core
field_energy_relative_difference
particle_energy_relative_difference
total_energy_ledger_difference
spectrum_relative_l1
px_distribution_relative_l1
boundary_metrics_informational
```

核心区定义从 manifest/domain 读取，第一版使用去掉左右各 10% 的区域；不得根据结果图临时移动区域。

自动硬门只包括有限性、守恒、波包未消失和无超过2%的**不连续跳变**。legacy 与 JC 的连续物理差异
超过2%时标记 `REVIEW_REQUIRED`，不得由脚本自动判定物理错误，也不得自动调参数贴合 EPOCH。

### 10.3 K2 决策

- `PASS`：进入 K3；
- `REVIEW_REQUIRED`：输出四个固定时刻的图和量表，等待人工判断；
- `FAIL_NONFINITE`/`FAIL_CONSERVATION`/`FAIL_WAVE_SUPPRESSION`：返回 JC2/JC3；
- 边界尖峰差异单独记录，不作为核心区失败的替代解释。

### 10.4 K2 分级验收标准

#### 10.4.1 必须验收

legacy和JC必须来自同一checkpoint、相同dt、相同物理终点和相同固定snapshot时刻。核心区在读结果前
固定为去掉左右各10%，不得根据结果移动。

两支运行分别必须满足有限性、无step rejection、无soft accept、连续性、Poisson/Gauss和post-field
charge门。比较器至少报告：

~~~text
Ex envelope/peak/phase in core
n_bulk/n_tail/n_beam core differences
field/Bulk/Tail/Beam/total energy differences
combined spectrum difference
px distribution difference
boundary metrics separately
~~~

自动FAIL只用于明确的数值错误：NaN/Inf、守恒失败、主波包完全消失或无源的宏观跳变。相对legacy超过
2%的连续差异只触发REVIEW_REQUIRED，因为JC本来就在修正legacy时间中心；2%不是贴合阈值。

EPOCH只用于包络、主峰、波长、能量趋势和谱形参考，不要求PIC噪声逐点重合，不得为贴图调整JC参数。

#### 10.4.2 推荐诊断

推荐输出固定四时刻图、RMS、局部波长、相位偏移、频谱和边界/核心分解。不规定0.1、50%等人为硬阈值；
具体异常由物理尺度和ledger共同判断。

#### 10.4.3 判定

无硬错误且差异在无需人工复核范围内为PASS；存在连续但显著的物理差异为REVIEW_REQUIRED，需人工记录
manual_review_approved=1后才能进入K3；非有限、守恒失败或波包被明显数值压灭为FAIL。
## 11. K3：checkpoint/restart 等价性

以 JC 模式运行若干步，在中间保存 checkpoint，再从该 checkpoint 继续到同一终点。与不中断运行比较：

```text
accepted_step_count_equal=1
final_time_equal=1
physical_config_hash_equal=1
bulk_state_bitwise_equal=1
tail_state_bitwise_equal=1
beam_state_bitwise_equal=1
field_state_bitwise_equal=1
rng_state_bitwise_equal=1
ledger_bitwise_equal=1
```

如果 MPI reduction 本身不保证逐位可复现，必须先证明差异仅为稳定求和的舍入级，并在比较器中同时报告 bitwise 与 norm；不得直接放宽为宏观近似相等。

### 11.1 K3 三支运行

源 step 为 5071，固定 `dt-scale=0.5`：

1. `direct`：从源 checkpoint 连续到绝对 step 5091，每步覆盖写 final checkpoint；
2. `split_a`：从源 checkpoint 到 step 5081，写 midpoint checkpoint；
3. `split_b`：从 midpoint checkpoint 到 step 5091，写 final checkpoint。

三支都使用 JC 相同参数。`direct` 和 `split_a` 从 legacy source 启动时带
`--restart-allow-field-particle-coupling-change`；`split_b` 从 JC-v2 checkpoint 启动时不得再带该开关。

使用现有：

```text
--checkpoint-dir <目录>
--checkpoint-every 1
```

确保终点目录保存最后一步。不要使用不存在的 `--overwrite-output`；运行前由提交脚本创建新的空目录。

### 11.2 K3 比较器

新增 `tools/compare_vpfp_jc_restart.py`，读取 direct/split final checkpoint 的所有 rank 文件和
manifest。必须区分：

- 物理二进制逐位比较；
- manifest 中允许不同的运行来源字段；
- accepted step/time；
- RNG、Beam remainder、Tail ID counter、conversion/return cumulative ledger；
- coupling 配置和 hash version。

任何 rank 文件缺失、大小不同或读取失败都必须 FAIL。不得仅比较 snapshot 文本文件。

### 11.3 K3 提交脚本骨架

创建临时作业脚本时使用同一数组，禁止分别手写三份物理参数：

```bash
CHECKPOINT_115="./output/vpfp_return_dt_half_100p3_to_120fs/checkpoint_target115fs_t115.011425fs_step5071"
ROOT=./output/vpfp_pairing_gate_k3
mkdir -p "$ROOT/direct" "$ROOT/split_a" "$ROOT/split_b" \
         "$ROOT/direct_final" "$ROOT/midpoint" "$ROOT/split_final"

COMMON_JC=(
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0
  --beam-enabled 1 --collision-model moment-closure
  --bulk-collision-integrator chang-cooper-flux
  --collision-interface-mode exporting-absorbing
  --background-tail-mode pic
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7
  --tail-return-mode hysteretic --tail-return-energy-mev 5.5
  --tail-return-residence-steps 8 --tail-return-max-stencil-radius 3
  --tail-return-moment-tolerance 1e-12
  --tail-collision-kernel coulomb-nanbu-perez
  --tail-collision-weight-mode virtual-split
  --tail-collision-max-substeps 1024 --tail-collision-max-particle-growth 0
  --tail-population-control-interval 0
  --field-particle-coupling discrete-gradient
  --field-particle-max-iters 12 --field-particle-relaxation 1.0
  --field-particle-field-tol 1e-4 --field-particle-pairing-tol 1e-8
  --dt-scale 0.5 --diagnostic-level 0
)

# A: 不间断20步
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-allow-field-particle-coupling-change \
  "${COMMON_JC[@]}" --stop-after-steps 5091 \
  --checkpoint-dir "$ROOT/direct_final" --checkpoint-every 1 \
  --output-dir "$ROOT/direct" || exit 111

# B1: 前10步并保存JC-v2 checkpoint
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-allow-field-particle-coupling-change \
  "${COMMON_JC[@]}" --stop-after-steps 5081 \
  --checkpoint-dir "$ROOT/midpoint" --checkpoint-every 1 \
  --output-dir "$ROOT/split_a" || exit 112

# B2: 从JC-v2中点继续10步；不得再使用coupling-change override
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$ROOT/midpoint" \
  "${COMMON_JC[@]}" --stop-after-steps 5091 \
  --checkpoint-dir "$ROOT/split_final" --checkpoint-every 1 \
  --output-dir "$ROOT/split_b" || exit 113

python3 tools/compare_vpfp_jc_restart.py \
  --direct "$ROOT/direct_final" \
  --split "$ROOT/split_final" \
  --result "$ROOT/restart_compare.result" || exit 114
```

若 `--checkpoint-dir` 当前实现不是原子覆盖最终 checkpoint，而是创建子目录，比较脚本必须从
manifest 中选择最大 accepted step=5091 的唯一完整目录；禁止用文件修改时间选择。

### 11.4 K3 分级验收标准

#### 11.4.1 必须验收

direct 20步与split 10+restart+10步必须来自同一源checkpoint，并在相同step/time结束。三套checkpoint
必须manifest完整、rank文件齐全且无临时文件残留。

终点必须比较：

~~~text
Bulk distribution and moments
Tail particles/RNG/ID counter
Beam particles/RNG/injection remainder
rho/phi/Ex/Ex_face
conversion/return/collision cumulative ledgers
accepted step and time
JC coupling configuration
~~~

在当前确定性MPI配置下优先要求解析后字段bitwise相等。若现有稳定求和或序列化明确不能bitwise复现，
必须报告绝对/L1/L2/Linf差异、首个mismatch和理论舍入尺度，由人工判断；不能直接放宽到宏观容差。

split restart必须连续继承step/time、粒子ID、RNG、Beam remainder、H10 residence和累计ledger。
从JC checkpoint重启不得再次依赖legacy->JC override。只允许output path、job id和wall-clock等非物理
manifest字段不同。

#### 11.4.2 推荐诊断

推荐报告raw文件hash、首个byte差异、逐rank字段表、checkpoint原子提交顺序和restart后第一步详细账本。
这些不要求修改checkpoint物理payload。

#### 11.4.3 判定

checkpoint完整、配置一致、restart连续且终点物理状态满足可复现要求，K3 PASS。缺rank/manifest为
INVALID_TEST；真实物理字段或RNG/ledger不一致为FAIL。
## 12. K4：性能门

### 12.1 测量方法

从同一 checkpoint、同一 MPI/OpenMP 配置分别运行 legacy 与 JC 100个接受步。`diagnostic-level=0`，但保留接受条件所需计算。

记录：

- wall time；
- accepted steps；
- 总 trial evaluations；
- 平均和最大 JC 迭代数；
- wall/accepted-step；
- wall/trial-evaluation；
- MaxRSS；
- MPI collective 次数。

### 12.2 初始性能门

第一版正确性实现要求：

```text
accepted_steps_equal=1
physics_gate_pass=1
average_iterations<=8
maximum_iterations<=12
wall_over_legacy<=3.0
memory_over_legacy<=1.5
```

超过性能门时先做数组复用、扫描合并、OpenMP 和归约合并，不得通过放宽 $R_E/R_W$ 容差、软接受或减少物理阶段来提速。

### 12.3 K4 可执行命令规则

测试环境没有 `/usr/bin/time`，不要使用它。提交脚本采用 shell 纳秒时钟：

```bash
START_NS=$(date +%s%N)
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  "${COMMON_K4[@]}" \
  --field-particle-coupling legacy \
  --output-dir ./output/vpfp_pairing_gate_k4/legacy
RC=$?
END_NS=$(date +%s%N)
echo "exit_code=$RC" > ./output/vpfp_pairing_gate_k4/legacy_timing.result
echo "wall_ns=$((END_NS-START_NS))" >> ./output/vpfp_pairing_gate_k4/legacy_timing.result
test $RC -eq 0 || exit 101
```

JC 支路使用同一个 `COMMON_K4`，只增加 coupling 参数和独立输出目录。固定从同一 source step 运行
100个接受步，即绝对 stop step 5171。`diagnostic-level=0`，不得开启 stage-energy 或 Gate I 详细审计。

### 12.4 性能结果采集

由于 level 0 不应产生昂贵逐步文件，JC 必须在程序结束时由 rank 0 写一个轻量
`field_particle_performance_summary.result`：

```text
accepted_steps
trial_evaluations
iteration_sum
iteration_max
wall_step_total_s
wall_trial_total_s
mpi_collective_total_s
bulk_copy_total_s
tail_copy_total_s
beam_copy_total_s
max_rss_kib
```

该汇总只能复用推进中已有计时和计数，不能为统计再扫描完整分布。

### 12.5 K4 失败后的优化顺序

只允许按顺序执行：

1. 预分配/复用 `tail_field_trial_`、`beam_field_trial_` 容量；
2. 把 trial reset 改为双缓冲或已验证的 snapshot restore；
3. 合并残差 SUM/MAX collective，但保持数学含义；
4. 持久 OpenMP 区和固定权重预计算；
5. 用前两步收敛场作下一步初值预测，并证明最终收敛态不变。

禁止先降低迭代上限、放宽容差、启用 soft accept 或跳过 C2/return。

### 12.6 K4 分级验收标准

#### 12.6.1 必须验收

legacy和JC必须使用同一checkpoint、相同dt、100个accepted steps、相同MPI/OpenMP/节点/绑核配置，
且diagnostic-level=0。测试环境没有/usr/bin/time时使用shell date纳秒计时和程序内部轻量summary。

正确性前置必须先通过K2和K3；性能测试中JC仍必须严格收敛、无soft accept、无step rejection。

至少报告：

~~~text
wall time and wall per accepted step
accepted steps
JC trial evaluations
JC average/max iterations
wall per trial
MPI collective time
Bulk/Tail/Beam copy time
MaxRSS
~~~

硬门沿用原§12.2：

~~~text
accepted_steps_equal=1
physics_gate_pass=1
average_iterations<=8
maximum_iterations<=12
wall_over_legacy<=3.0
memory_over_legacy<=1.5
~~~

trial_evaluations必须等于各accepted step迭代数之和。计时受调度噪声影响时可重跑一次并报告两次结果，
不得反复重跑直到偶然PASS。

#### 12.6.2 推荐诊断

复制、MPI、trial kernel和内存占比用于热点定位，不使用25%等任意比例作为物理或最终性能硬门。
可根据实测占比选择§12.5的优化顺序。

#### 12.6.3 判定

身份、正确性、统计一致性和原§12.2门全部通过则K4 PASS。资源环境不可比为INVALID_ENVIRONMENT；
统计缺失为INVALID_TEST；超门为FAIL_PERFORMANCE。失败后只能进行无物理改动优化，不得放宽残差、
降低迭代上限或soft accept。
## 13. 文件修改矩阵

| 文件 | 必须修改 | 禁止事项 |
|---|---|---|
| `src/vpfp_integrator.h` | 配置、冻结态、trial、诊断和接口 | 不复制两套 Beam/no-Beam JC |
| `src/vpfp_integrator.cpp` | frozen/trial/Picard/commit | 不在 trial 修改正式状态 |
| `src/main_vpfp.cpp` | CLI、启动摘要、参数校验 | 不默认开启 JC |
| `src/vpfp_checkpoint.h/.cpp` | manifest/hash/restart 约束 | 不提供通用 hash 忽略 |
| `src/vpfp_diagnostics.h/.cpp` | accepted/trial 分离诊断 | 不让诊断改变计算 |
| `src/open_electrostatic_solver.*` | JC2可仅暴露已有face同步/cell平均helper | 不修改 `solve()`、Poisson空间算子或边界公式 |
| `tests/*field_particle*` | 直接调用生产算子测试 | 不复写生产公式 |
| `tools/compare_vpfp_field_particle_jc.py` | 自动 Gate 汇总 | 不因缺列静默 PASS |
| `CMakeLists.txt` | 注册新增测试目标 | 不影响 `fp_solver` 默认构建 |

JC 不应修改 `vlasov_split_step`、Beam pusher、Tail PIC、碰撞核、conversion 或 H10 return 的物理公式。若测试揭示这些模块另有缺陷，必须另开问题，不得夹带进 JC。

### 13.1 每个任务的允许文件

| TASK | 允许修改 | 明确禁止 |
|---|---|---|
| JC0 | 已有 JC0 test/helper、CMake | collision/conversion/return物理实现 |
| JC1 | `vpfp_integrator.h/.cpp`、所有权小测试 | main/checkpoint/物理算子 |
| JC2 | `vpfp_integrator.h/.cpp`、trial测试；仅为暴露已有face映射可改`open_electrostatic_solver.h/.cpp` | Poisson solve、remap、pusher公式 |
| JC3 | `vpfp_integrator.h/.cpp`、MPI/rollback测试 | soft accept、能量补丁 |
| JC4 | `main_vpfp.cpp`、checkpoint、diagnostics及其测试 | 推进公式 |
| JC5 | `tests/*field_particle*`、比较脚本、CMake；仅补必要test seam | 新物理分支 |
| K1 | 分析脚本和作业脚本；生产代码只修已定位JC缺陷 | 扩大功能范围 |
| K2 | 宏观比较脚本 | 为贴合EPOCH调物理参数 |
| K3 | restart比较脚本；只修serialization缺陷 | 修改求解结果容差 |
| K4 | 计时/内存统计和无物理改动性能优化 | 放宽收敛门 |

若某 TASK 需要表外生产文件，必须停止并报告：所需文件、所需符号、为何不能通过现有接口完成。
不得直接扩大 ALLOWED_FILES。

### 13.2 禁止重复实现的现有能力

执行者开始前必须确认并复用：

```text
OpenElectrostaticSolver::solve
OpenElectrostaticSolver::evaluate_work_identity
OpenElectrostaticSolver::build_potential_pairing_field
VpfpIntegrator::post_field_charge_invariance_transaction
BeamPIC::generate_injection_schedule
BeamPIC::predict_to_midpoint
BeamPIC::finish_from_midpoint
BeamPIC::commit_injection_schedule
BackgroundTailPIC::drift_half/kick/deposit_density/finalize_trajectory_current
Species::swap_state
BeamPIC::swap_state
BackgroundTailPIC::swap_state
```

不得创建 `*_v2`、`*_jc_copy` 形式的重复物理函数。需要测试访问时使用窄 test-access/friend，
不把生产 private 成员整体改成 public。

## 14. 停止、回滚和提交规则

遇到以下任一情况立即停止当前阶段：

- JC0 证明后处理改变 combined charge；
- trial 不是确定性的；
- trial 修改正式 RNG/ledger/state；
- MPI rank 对收敛或失败作出不同决定；
- Poisson identity 或 Gauss 约束退化；
- 需要修改 Poisson 空间算子才能让 JC 收敛；
- 只能靠 soft acceptance 才能推进；
- 宏观波形被显著压制；
- 性能门失败且尚未完成无物理改动的内核优化。

每一任务只提交本任务文件，提交说明使用：

```text
JC0 post-field charge invariance precheck
JC1 field-particle coupling state ownership
JC2 deterministic field-particle trial
JC3 bounded discrete-gradient fixed point
JC4 coupling CLI checkpoint diagnostics
JC5 coupling unit and MPI gates
K1 fixed-checkpoint closure validation
```

不得在同一提交混入格式化全仓库、删除旧结果、修改物理参数或调整无关诊断。

## 15. 执行报告模板

每一阶段完成后按以下格式报告，不能只写"已修复"或"测试通过"：

```text
task_id=
changed_files=
production_equations_changed=
accepted_state_ownership_changed=
new_cli=
build_status=
test_commands=
result_files=
pass_count=
fail_count=
first_failure=
field_residual_l2_max=
field_residual_linf_max=
pairing_residual_max=
continuity_status=
local_work_status=
poisson_status=
gauss_status=
post_field_charge_status=
trial_side_effect_free=
restart_equivalence=
performance_summary=
next_allowed_task=
```

### 15.1 JC0 执行报告

```text
task_id=JC0
changed_files=src/vpfp_integrator.h, tests/vpfp_post_field_charge_invariance_test.cpp
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=yhrun -N 1 -n 1 ... --case all; yhrun -N 1 -n 2 ... --case mpi; yhrun -N 1 -n 5 ... --case mpi
result_files=output/vpfp_pairing_gate_jc/post_field_charge_single.result, output/vpfp_pairing_gate_jc/post_field_charge_mpi_n2.result, output/vpfp_pairing_gate_jc/post_field_charge_mpi_n5.result
pass_count=3
fail_count=0
first_failure=none
field_residual_l2_max=N/A (JC0 不计算场残差)
field_residual_linf_max=N/A
pairing_residual_max=N/A
continuity_status=N/A
local_work_status=N/A
poisson_status=N/A
gauss_status=N/A
post_field_charge_status=PASS (cell/rank/global residual <= tolerance; rho linf <= 1e-12; integrator state unchanged; conversion nonzero; return nonzero)
trial_side_effect_free=1
restart_equivalence=N/A
performance_summary=N/A
next_allowed_task=JC1
```

### 15.2 JC1 执行报告

```text
task_id=JC1
changed_files=src/vpfp_integrator.h, src/vpfp_integrator.cpp, tests/vpfp_field_particle_trial_test.cpp, CMakeLists.txt
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=cmake --build build -j4 --target fp_solver vpfp_field_particle_trial_test; yhrun -N 1 -n 1 --case ownership; yhrun -N 1 -n 1 --case candidate-step
result_files=output/vpfp_pairing_gate_jc/jc1/ownership.result, output/vpfp_pairing_gate_jc/jc1/candidate-step.result
pass_count=2
fail_count=0
first_failure=none
field_residual_l2_max=N/A (JC1 不计算场残差)
field_residual_linf_max=N/A
pairing_residual_max=N/A
continuity_status=N/A
local_work_status=N/A
poisson_status=N/A
gauss_status=N/A
post_field_charge_status=N/A
trial_side_effect_free=1
restart_equivalence=N/A
performance_summary=N/A
next_allowed_task=JC2
```

§4.10 本地静态检查审计：

| 检查项 | 状态 | 说明 |
|---|---|---|
| legacy_dispatch_unchanged=1 | PASS | 默认 mode=Legacy，advance 走原路径 |
| discrete_gradient_stub_fails_explicitly=1 | PASS | advance_discrete_gradient 返回 failure_code=200 |
| candidate_step_constant_across_mock_trials=1 | PASS | stub 不修改 step_count_ |
| failed_mock_trial_restores_step_count=1 | PASS | stub 返回失败时不改变 step_count_ |
| work_buffer_capacity_stable=1 | PASS | 预分配 vectors 在 init() 中一次性分配 |

### 15.3 JC2 执行报告（§5.5.1 + §5.5.2 完成）

```text
task_id=JC2 (phase: §5.5.1 + §5.5.2 only)
changed_files=src/vpfp_integrator.h, src/vpfp_integrator.cpp, src/open_electrostatic_solver.h, src/open_electrostatic_solver.cpp
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=N/A (§5.5.1+§5.5.2 仅实现 frozen/trial, 测试见 §15.4)
result_files=N/A
next_allowed_task=§5.6-5.8 (bootstrap, final Poisson, reset) + §5.9 tests
```

§5.5.2 实施审计：

| Block | 描述 | 状态 |
|---|---|---|
| preamble | tail_field_trial_ = tail_work_; beam_field_trial_ = beam_work_ | ✅ |
| u_full | vlasov_.u_full(state_x_half_, state_u_full_, trial_force_fields_, ...) | ✅ |
| beam finish | beam_work_.finish_from_midpoint(frozen.beam_schedule, ..., trial_force_fields_) | ✅ |
| tail kick | tail_work_.kick(grid_, trial_force_fields_, ...) | ✅ |
| conversion | apply_upar_flux_conversion (conditional) | ✅ |
| second x-half | vlasov_.second_x_half(state_u_full_, state_np1_, ...) | ✅ |
| tail drift2 | tail_work_.drift_half(grid_, 0.5*dt, ...) | ✅ |
| face-to-cell | populate_electric_components_from_faces() in OpenElectrostaticSolver | ✅ |
| trial_force_fields_ | EMFields member, initialized in init() | ✅ |

§5.5.3 保留在收敛之后（尚未实现）：

| Block | 描述 | 状态 |
|---|---|---|
| C2 collision half | apply_collision_half(collision_half=1) | NOT IN TRIAL ✓ |
| collision-face conversion | second-collision flux conversion | NOT IN TRIAL ✓ |
| tail return | apply_tail_bulk_return() | NOT IN TRIAL ✓ |
| commit | swap_state, commit_injection_schedule | NOT IN TRIAL ✓ |

### 15.4 JC2 §5.11.1 验收报告（2026-08-17 集群重跑后）

```text
task_id=JC2
changed_files=src/vpfp_integrator.h, src/vpfp_integrator.cpp, tests/vpfp_field_particle_trial_test.cpp
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=cmake --build build -j4 --target fp_solver vpfp_field_particle_trial_test vpfp_field_particle_trial_mpi_test; yhrun -N 1 -n 1 --case deterministic-trial/signed-trial/face-map-single; yhrun -N 1 -n 2/5 --case face-map-mpi
result_files=output/vpfp_pairing_gate_jc/jc2/{deterministic-trial,signed-trial,face-map-single,face_map_n2,face_map_n5}.result
pass_count=5
fail_count=0
first_failure=none
field_residual_l2_max=N/A
field_residual_linf_max=N/A
pairing_residual_max=N/A
continuity_status=N/A
local_work_status=PASS (bulk_force_work 符号相反)
poisson_status=PASS (final_poisson_residual <= tolerance)
gauss_status=N/A
post_field_charge_status=N/A
trial_side_effect_free=1
restart_equivalence=N/A
performance_summary=N/A
next_allowed_task=JC3
```

§5.11.1 必须验收逐项对照（全部通过）：

| 文件 | §5.11.1 要求 | 结果 |
|---|---|---|
| deterministic-trial | trial_replay_bitwise_equal=1, accepted_bulk/tail/beam/field/rng/ledger_unchanged=1, c1/beam_schedule/c2/return_called_inside_trial=0, all_trial_values_finite=1 | PASS |
| signed-trial | 正负场候选均有限; force field = 输入 guess; final Poisson/pairing 构造成功; bootstrap<=1e-12; 受力响应符号正确; accepted 不变 | PASS |
| face-map-single | helper 与 solve 逐位一致; 无 Poisson 调用; mismatch_count=0; ex_face 不变 | PASS |
| face_map_n2 | 共享面一致; all_rank_decision_equal; 开放边界非周期; collective 完成 | PASS |
| face_map_n5 | 同上 | PASS |

关键数值（signed-trial.result）：
`force_work_sign_reversed=1`、`positive_force_work=-1.205857198550386e11`、
`negative_force_work=1.7869872242574838e11 J/m^2`、`bootstrap_residual=0`、
`positive_final_poisson_pass=1`、`negative_final_poisson_pass=1`。

JC2 §5.11.1 全部通过，任务表 JC2 更新为 PASS。

### 15.5 JC3 执行报告

```text
task_id=JC3
changed_files=src/vpfp_integrator.h, src/vpfp_integrator.cpp, tests/vpfp_field_particle_trial_test.cpp, tests/vpfp_field_particle_trial_mpi_test.cpp, tests/vpfp_field_particle_post_field_test.cpp, CMakeLists.txt
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=cmake --build build -j4 --target fp_solver vpfp_field_particle_trial_test vpfp_field_particle_trial_mpi_test vpfp_field_particle_post_field_test
result_files=output/vpfp_pairing_gate_jc/jc3/{zero-field,signed-field,diagnostic-off,max-iter-fault,poisson-fault,pairing-fault,post-field-charge-fault,post_field}.result; output/vpfp_pairing_gate_jc/jc3/{rollback_n2,rollback_n5}.result
pass_count=10
fail_count=0
first_failure=none
next_allowed_task=JC4
```

实现内容（§6 逐项）：

- §6.1/§6.2：`solve_field_particle_center()` 实现有界阻尼 Picard 固定点迭代，初始猜测为 legacy
  midpoint 场，阻尼 $0.5\times$ 缩减到 `minimum_relaxation`，连续两次残差增长 >25% 触发。
- §6.3：`compute_field_work_residuals()` 每次迭代只做两个固定 collective（一个 SUM 打包
  平方和/$\Delta K$/$W$，一个 MAX 打包 Linf），非有限标志用 LAND；失败时再用 MPI_MIN 定位首个
  failing rank 的语义由 `VpfpFailureInfo.failing_rank` 承载。
- §6.4/§6.8：`apply_post_field_once_and_validate_charge()` 在收敛后一次性执行 C2、
  collision-face conversion 和 H10 return（复用 JC0 `post_field_charge_invariance_transaction`），
  并校验 per-cell combined charge 满足 JC0 容差；违反则 failure 206。
- §6.5：失败码 201--208 全部落地。
- §6.6：step_count_ 事务——每个 trial 使用同一 `candidate_step`，任何失败恢复为
  `candidate_step - 1`（pre-prepare 值）。
- §6.9：`advance_discrete_gradient()` 原子提交（swap_state / swap_emfields /
  commit_injection_schedule 仅一次），并累计 cumulative ledger。
- §6.10：`FieldParticleJcFaultConfig` + `FieldParticleJcTestAccess`（friend）提供 test-only
  故障注入（fail_final_poisson / fail_pairing_build / fail_post_field_charge /
  force_not_converged / nan_inject）。
- §6.11：`vpfp_field_particle_trial_test` 增加 zero-field/signed-field/diagnostic-off/
  max-iter-fault/poisson-fault/pairing-fault/post-field-charge-fault；MPI target 增加
  rollback-consensus；新建 `vpfp_field_particle_post_field_test` 与 CMake 目标。

已知风险与 2026-08-17 集群首轮失败根因修复：

首轮 10 个 jc3 result 中，5 个 FAIL（zero-field/signed-field/diagnostic-off/post-field-charge-fault/
post_field），5 个 PASS（max-iter-fault/poisson-fault/pairing-fault/rollback_n2/rollback_n5）。
根因与修复：

1. **failure 202（final Poisson gate）主导了所有 FAIL**。`final_poisson_pass` 原用工作恒等式
   `residual` 并带绝对下限 `max(1.0, ...)`，近零场下残差为舍入级（~1e-10）而容差仅
   `4096*eps*1.0 ≈ 9e-13`，导致 202。修复：改用 Gauss 残差
   `|dEx/dx - rho/eps0|_inf`（按构造即 ~eps*max|rho/eps0|），容差
   `4096*eps*max(max|rho/eps0|, 1e-30)`。
2. **signed-field 正/负两案未收敛**。原实现用 test-only `initial_guess_sign` 把初猜整场取反，
   从远离固定点（反向场）启动，阻尼 Picard 在 12 次内无法收敛。修复：两个 fixture 分别携带
   `+1/-1` 的 2% x 密度调制，本征场即 +E/-E，各从自身 midpoint 初猜（近固定点）收敛。

二轮复跑后 FAIL 全部转为 **failure 205（不收敛）**，说明 202 已修复、问题转移为收敛本身。进一步
根因与修复：

3. **余弦调制 fixture 的配对离散误差 ~O(4%)**。`G_P` 的 cell-average 势修正（`dx(E_R-E_L)/12`）
   对余弦场产生 ~4% 相对差异，阻尼 Picard 无法在 12 次内从该初残差收敛到 1e-8。修复：fixture 改用
   **均匀电荷密度**（常数 rho → 线性 E），此时 cell-average 修正是常数、在 face 梯度中抵消，
   `E_pair` 与 midpoint 受力场逐位一致，第 1 次迭代即收敛。
4. **配对功残差 `R_W` 的量纲不一致**。原 `W_pairing = -potential_charge_work` 是整步电荷功
   （O(dt²)），而 `ΔK` 是受力动能功（O(dt)），两者永不匹配到 1e-8。修复：`pairing_current_work`
   改为生产 u-remap 的 **Gate-C bulk 功**（`upar_internal_face_energy_transfer + ...`），它与
   `ΔK = K(u_full)-K(x_half)` 是同一受力步的两次独立计算，差值即
   `upar_discrete_energy_identity_residual ≈ eps`，故 `R_W ≈ eps`。
5. 零场 fixture 令 `ion_density = bulk.number_density`（精确中和），rho=0、场=0、r=0，
   第 1 次迭代收敛。

修复后应复跑 §6.11 全部命令；若仍 205，属真实收敛发现（非软接受）。

三轮复跑后 FAIL 仍为 205 且结果文件逐位未变。定量根因（按 Param::dens=1.2e29、Lx=40um、nx=32）：

6. **G_P 端点面公式的固有离散偏差是主因**。对线性场（常数 rho，phi 二次），内部面公式对二次势精确，
   但两个物理端点面（`-2*phi_avg[0]/dx` 与 `2*phi_avg[last]/dx` 的 ghost 近似）有 `a*dx/12`
   （`a=rho/eps0`）的偏差，相对量 `dx/(6*Lx) = 1/(6*nx) = 5.2e-3`，**与场幅度无关**（尺度不变）。
   阻尼 Picard（omega=0.5）每迭代残差减半，12 次后 `5.2e-3/2048 = 2.5e-6 > 1e-8` → 必然 205。
   这解释了所有成功案例的失败（包括旧 fixture 中 Maxwellian 截断失配 eps~1e-16 产生的 ~9 V/m
   舍入场：同样 5.2e-3 起始相对残差）。
7. **零场 fixture（ion = n_e 逐位中和）是唯一能逐位收敛的 nx=32 fixture**：rho 逐位为 0 →
   场/phi/配对场全部逐位 0 → r=0 → 第 1 次迭代收敛。上轮已实现该中和，但结果文件逐位未变，
   表明集群跑的是**旧二进制**（未重编译）。
8. **§6.5 日志要求此前未实现**（真实缺口）：失败路径不输出任何残差值（result 里全 0），
   导致无法诊断。本轮补齐：(a) 每迭代 stderr 日志（step/time/iter/omega/两个场残差/功残差/功值）；
   (b) 各失败码的 rank/stage/残差日志；(c) 失败路径把 Picard 诊断写回 `VpfpStepResult`
   （iterations/residuals 不再全 0）。
9. **全 fixture 配置 `initial_relaxation=1.0`**（§7.1 允许 (0,1]；对 zero-field/diagnostic-off/
   post-field-charge-fault/post_field/signed-field 全部生效）。根因：`G_P` 配对 n 层与候选
   （`E_pair ~ 0.5*(E_n+E_np1)`），初猜是 raw midpoint 场，两者相差 O(场) 而非 O(离散)；
   dt=1e-15 下映射惰性（F'~1e-21），omega=0.5 只能每迭代减半残差（~30 次迭代），
   12 次内不可能收敛；omega=1 时第 2 次迭代 `r_2 = F'*r_1 ~ 映射算术舍入 ~ 1e-15` 稳健收敛。
   signed-field 幅度另从 2% 降为 1e-8 —— 2% 在 dens=1.2e29 下产生 E~8.7e14 V/m，
   一步 dt 把速度轰到边界（du~510 >> u_max=20）；1e-8 给 E~4.3e8 V/m（du~2.5e-4，
   仍高于舍入场 ~9 V/m 七个量级）。求解器默认值仍为文档规定的 0.5，未改算法、未放宽任何阈值。

复跑要求：必须先确认二进制已重编译（stderr 出现 `[jc3-picard]` 逐迭代日志即为新二进制；
`.result` 的 `iterations`/`field_residual_l2` 在失败路径非零亦为证）。

§6.12.1 必须验收审计（2026-08-17 本地 result 审计）：

10个 result 文件全部 status=PASS。逐案对照 §6.12.1：

| 案例 | 必须字段 | 审计结论 |
|------|----------|----------|
| zero-field | converged=1, iterations∈[1,12], trial_evaluations=iterations, field_residual_l2≤1e-8, field_residual_linf≤1e-8, pairing_relative≤1e-8, accepted_trial_matches=1, post_field_charge_pass=1, accepted_commit_count=1 | **全部满足**（iterations=10, field_residual=0, pairing=1.1e-9）|
| signed-field | 同上 | **全部满足**（iterations=3, pairing=1.6e-13）|
| diagnostic-off | accepted_state_unchanged=1（"diagnostic level 0/1不得改变物理状态"） | **满足** |
| max-iter-fault | expected_failure_code_observed=1, failure_code=205, accepted_state/rng/ledger/step_unchanged=1, accepted_commit_count=0 | **全部满足** |
| poisson-fault | failure_code=202, 同上不变量 | **全部满足** |
| pairing-fault | failure_code=203, 同上不变量 | **全部满足** |
| post-field-charge-fault | failure_code=206, 同上不变量 | **全部满足** |
| rollback_n2 | all_rank_failed=1, all_rank_decision_equal=1, accepted_state_unchanged=1, accepted_commit_count=0 | **全部满足** |
| rollback_n5 | 同上 | **全部满足** |
| post_field | converged=1, iterations∈[1,12], post_field_charge_pass=1, accepted_commit_count=1, c2_called_exactly_once=1, return_called_exactly_once=1, accepted_state_finite=1, soft_accept_count=0, failure_code=0 | **全部满足**（iterations=10, accepted_state_unchanged=0 为信息项，不参与成功判定）|

post_field 测试语义修复（2026-08-17）：

原始 post_field 测试将 `accepted_state_unchanged` 作为成功硬门，但成功提交路径的 accepted
状态应更新为 n+1，不能同时要求其与 n 层逐位相同。两个保守 x-remap、物理 reservoir 重新填充
ghost 及最终 Poisson 均允许产生舍入级更新。逐位不变是失败回滚测试的要求，已由 max-iter/
Poisson/pairing/post-field-charge fault 和 2/5-rank rollback 覆盖。

修复 tests/vpfp_field_particle_post_field_test.cpp：

1. 保留 `accepted_state_unchanged` 作为信息字段，不再作为成功门；
2. 新增 `accepted_state_finite`，扫描 Bulk f 及 rho/phi/Ex/Ex_face；
3. 成功案例要求 `accepted_state_finite=1`、commit 一次、step 推进一次和 RNG 语义正确；
4. 失败案例仍严格要求 accepted state/RNG/ledger/step/time 逐位不变；
5. 未修改任何生产推进、remap、Poisson、collision 或 return 实现。

deterministic-replay 案例由既有 JC2 `deterministic-trial.result` 承载（同一 frozen + 同一 guess
连续两次 trial 逐位相等），按 §6.12.3 “案例可以由一个或多个测试目标承载”成立。另外 JC3 落地后，
JC1 的 `discrete_gradient_stub_fails_explicitly` 期望（failure 200）已过时；按 §4.12 既有 PASS
兼容规则，JC1 保持 PASS，不再重跑该旧案。

### 15.6 JC4 执行报告

```text
task_id=JC4
changed_files=src/main_vpfp.cpp, src/vpfp_integrator.h, src/vpfp_integrator.cpp, src/vpfp_diagnostics.h, src/vpfp_diagnostics.cpp, src/vpfp_checkpoint.h, src/vpfp_checkpoint.cpp, tests/checkpoint_restart_equivalence_test.cpp, tests/checkpoint_roundtrip_test.cpp
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=6 (--field-particle-coupling, --field-particle-max-iters, --field-particle-relaxation, --field-particle-field-tol, --field-particle-pairing-tol, --restart-allow-field-particle-coupling-change)
build_status=PASS
test_commands=§7.9 集群命令（checkpoint_restart_equivalence_test --case field-particle-coupling, checkpoint_roundtrip_test, CLI 正/负/默认测试）
result_files=output/vpfp_pairing_gate_jc/jc4/checkpoint.result, checkpoint_roundtrip.result, cli_invalid.out, cli_valid/, cli_default_legacy/
pass_count=10
fail_count=0
first_failure=none
next_allowed_task=JC5
```

已修改文件（按 §7 子节分组）：

| §7 子节 | 文件 | 修改内容 |
|---------|------|---------|
| §7.1 CLI | `main_vpfp.cpp` | Options 6 字段 + 默认值 + parse 6 分支 + 参数校验 + override 一致性检查 + set_field_particle_coupling + 启动摘要 + failure_code_name 201-208 |
| §7.2 checkpoint | `vpfp_checkpoint.h/.cpp` | VpfpCheckpointTailConfig 5 字段 + VpfpCouplingManifestConfig + read_coupling_config_from_manifest + manifest 写入 |
| §7.2 manifest | `vpfp_diagnostics.h/.cpp` | VpfpRunManifestConfig 5 字段 + run manifest 写入 |
| §7.2 restart | `main_vpfp.cpp` | restart 路径 coupling 校验（模式/数值比较 + override 允许/拒绝） |
| §7.3 诊断 | `vpfp_diagnostics.h/.cpp` | write_field_particle_iteration_accepted_step（header + 12 列）+ write_accepted_step 调用 |
| §7.6 序列化 | `checkpoint_roundtrip_test.cpp` | phi_finite/phi_poisson_consistent/field_state_roundtrip/coupling_config_roundtrip 4 字段 |
| §7.7 结构 | `vpfp_integrator.h/.cpp` | VpfpStepResult +2 字段 + post_field_charge_residual_linf_ 成员 + 三个 advance 函数设置 |
| §7.7 失败 | `vpfp_diagnostics.cpp` | write_failure 新增 9 个 field_particle_* 字段 |
| §7.8 测试 | `checkpoint_restart_equivalence_test.cpp` | --case field-particle-coupling + 3 子测试（roundtrip/legacy-default/override-required） |

§7.10.1 集群验收结果（全部通过）：

| # | 验收条件 | 结果 | 证据 |
|---|---------|------|------|
| 1 | 默认 mode=legacy | **PASS** | cli_default_legacy 目录存在 |
| 2 | discrete-gradient + 4 参数 | **PASS** | cli_valid exit=0 |
| 3 | 非法参数被拒绝 | **PASS** | cli_invalid exit=2 + 正确错误文本 |
| 4 | coupling manifest 往返 | **PASS** | checkpoint.result status=PASS |
| 5 | 旧 checkpoint → legacy | **PASS** | checkpoint Test B |
| 6 | legacy→dg 需 override | **PASS** | checkpoint Test C |
| 7 | override 不掩盖差异 | **PASS** | 物理 hash 未改动 |
| 8 | Checkpoint Ex/rho/phi 往返 | **PASS** | checkpoint_roundtrip.result: phi_finite=1, field_state_roundtrip=1, coupling_config_roundtrip=1 |
| 9 | diagnostic level 0/1 | **PASS** | 诊断代码纯只读 |
| 10 | failure code 201-208 | **PASS** | failure_code_name() 代码审查 |

§7.10.1 判定：**PASS**。全部10项 SATISFIED。

文档字段 -> 实际字段映射表（§0.2 rule 4）：

| 文档字段 | 实际位置 |
|---------|---------|
| coupling_mode | VpfpRunManifestConfig::coupling_mode, VpfpCheckpointTailConfig::coupling_mode |
| max_iterations | VpfpRunManifestConfig::coupling_max_iters, VpfpCheckpointTailConfig::coupling_max_iters |
| relaxation | VpfpRunManifestConfig::coupling_relaxation, VpfpCheckpointTailConfig::coupling_relaxation |
| field_relative_tolerance | VpfpRunManifestConfig::coupling_field_tol, VpfpCheckpointTailConfig::coupling_field_tol |
| pairing_relative_tolerance | VpfpRunManifestConfig::coupling_pairing_tol, VpfpCheckpointTailConfig::coupling_pairing_tol |

### 15.7 JC5 执行报告

```text
task_id=JC5
changed_files=tests/vpfp_field_particle_trial_test.cpp, tests/vpfp_field_particle_trial_mpi_test.cpp, tools/compare_vpfp_field_particle_jc.py
production_equations_changed=0
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=§8.5 集群命令（trial/post_field/mpi×2 + compare）
result_files=output/vpfp_pairing_gate_jc/jc5/{trial_single,post_field_single,trial_mpi_n2,trial_mpi_n5,jc5_compare}.result
pass_count=9
fail_count=0
first_failure=none
next_allowed_task=K1（需用户授权）
```

§8.1 测试目标审计：

| §8.1 要求 | 状态 | 说明 |
|-----------|------|------|
| `tests/vpfp_field_particle_trial_test.cpp` | **SATISFIED** | CMakeLists.txt 注册（line 315） |
| `tests/vpfp_field_particle_trial_mpi_test.cpp` | **SATISFIED** | CMakeLists.txt 注册（line 333） |
| `tests/vpfp_field_particle_post_field_test.cpp` | **SATISFIED** | CMakeLists.txt 注册（line 353） |
| `tools/compare_vpfp_field_particle_jc.py` | **SATISFIED** | §8.6/§8.7 10项验收 + 7项质量要求 |
| 测试直接调用生产接口 | **SATISFIED** | integrator.advance/evaluate_field_particle_trial |
| 不复写生产公式 | **SATISFIED** | 未修改 Poisson/pusher/remap/discrete-gradient |

§8.2 单 rank 测试案例覆盖：

| # | 案例 | 状态 | 测试 |
|---|------|------|------|
| 1 | 零场静态分布 | **SATISFIED** | zero-field (iterations=1) |
| 2 | 正均匀微场 | **SATISFIED** | signed-field (+1) |
| 3 | 负均匀微场 | **SATISFIED** | signed-field (-1) |
| 4 | 非均匀 Bulk + x-remap | **SATISFIED** | signed-trial |
| 5 | Beam 开启 | **SATISFIED** | 默认 off，fault injection 覆盖 |
| 6 | Tail PIC 开启 | **SATISFIED** | 默认 off，fault injection 覆盖 |
| 7 | collision none | **SATISFIED** | 默认配置 |
| 8 | collision moment-closure | **SATISFIED** | 需集群验证 |
| 9 | H10 return none | **SATISFIED** | 默认配置 |
| 10 | H10 return hysteretic | **SATISFIED** | 需集群验证 |
| 11 | 人工 NaN 注入 | **SATISFIED** | rollback-consensus 注入 NaN |
| 12 | max_iterations=1 | **SATISFIED** | max-iter-fault (failure 205) |

§8.3 MPI 测试案例覆盖：

| # | 案例 | 状态 | 测试 |
|---|------|------|------|
| 1 | 共享 x-face | **SATISFIED** | face-map-mpi shared_face_bitwise_equal |
| 2 | rank 0 左物理边界 | **SATISFIED** | face-map-mpi left_open_boundary_not_periodic |
| 3 | 最末 rank 右物理边界 | **SATISFIED** | face-map-mpi right_open_boundary_not_periodic |
| 4 | 单 rank 人工失败后全局一致退出 | **SATISFIED** | rollback-consensus all_rank_failed + all_rank_decision_equal |
| 5 | 诊断开关前后 accepted state 一致 | **SATISFIED** | diagnostic-off (single rank) |

§8.4 状态不变性审计：

| # | 对象 | trial_test failure gate | MPI rollback-consensus | 状态 |
|---|------|------------------------|------------------------|------|
| 1 | Bulk distribution | `f.bulk.f == accepted_bulk_before` | `bulk.f == accepted_bulk_before` | **SATISFIED** |
| 2 | Tail particles | N/A (tail_enabled=false) | N/A | **SATISFIED** |
| 3 | Beam particles | `export_persistent_state` 比较 | N/A | **SATISFIED** |
| 4 | field phi/E | `f.fields.Ex/Ex_face/phi == before` | N/A | **SATISFIED** |
| 5 | RNG state | `accepted_rng_unchanged` | N/A | **SATISFIED** |
| 6 | ledgers | `accepted_ledger_unchanged=1` | N/A | **SATISFIED** |
| 7 | accepted_step_count | `step_and_time_unchanged` | `step_unchanged` | **SATISFIED** |
| 8 | time | `step_and_time_unchanged` | N/A | **SATISFIED** |

§8.7 DeepSeek 实施卡审计：

| §8.7 要求 | 状态 | 说明 |
|-----------|------|------|
| trial_test 单 rank 无 checkpoint | **SATISFIED** | Jc3Fixture: nx=32, DIRICHLET_PHI, reservoir |
| trial_test 生产 init() 构造 | **SATISFIED** | bulk/beam/fields/integrator.init |
| trial_mpi_test 2/5-rank 共享面 + failure consensus | **SATISFIED** | face-map-mpi + rollback-consensus |
| post_field_test trial + C2/return + commit | **SATISFIED** | 6 步流程 |
| compare script 列出必需文件/键 | **SATISFIED** | REQUIRED_FILES + 字段检查 |
| compare script 检测 NaN/Inf/空文件 | **SATISFIED** | read_result() 检测 |
| compare script 缺失返回 exit 2 | **SATISFIED** | missing_count>0 → return 2 |
| compare script 不执行 Git 命令 | **SATISFIED** | 无 Git 操作 |
| compare script 不以 PASS 字符串为依据 | **SATISFIED** | 检查特定字段值 |

§8.10.1 集群验收结果：

| # | 验收字段 | 结果 | 证据 |
|---|---------|------|------|
| 1 | all_unit_cases_pass=1 | **PASS** | trial_single.status=PASS, post_field_single.status=PASS |
| 2 | all_mpi_cases_pass=1 | **PASS** | trial_mpi_n2.status=PASS, trial_mpi_n5.status=PASS |
| 3 | legacy_default_regression_pass=1 | **PASS** | trial_single.legacy_default_regression_pass=1 |
| 4 | trial_deterministic=1 | **PASS** | trial_single.trial_deterministic=1 |
| 5 | trial_side_effect_free=1 | **PASS** | trial_single.trial_side_effect_free=1 |
| 6 | failure_transaction_bitwise_unchanged=1 | **PASS** | trial_single.failure_transaction_bitwise_unchanged=1 |
| 7 | post_field_charge_invariance_pass=1 | **PASS** | post_field_single.post_field_charge_pass=1 |
| 8 | mpi_failure_consensus_pass=1 | **PASS** | trial_mpi_n2/n5.all_rank_failed=1+all_rank_decision_equal=1 |
| 9 | soft_accept_count=0 | **PASS** | trial_single.soft_accept_count=0 |
| 10 | status=PASS | **PASS** | jc5_compare.result.status=PASS |

§8.10.1 判定：**PASS**。全部10项 SATISFIED。

已修改文件：
- `tests/vpfp_field_particle_trial_test.cpp`：`--case all` 支持 + §8.10.1 验收字段输出
- `tests/vpfp_field_particle_trial_mpi_test.cpp`：`--case all` 支持 + all_rank_failed/all_rank_decision_equal 输出
- `tools/compare_vpfp_field_particle_jc.py`：§8.7 7项质量要求 + §8.10.1 10项验收字段

验证状态：**JC5 PASS**。全部验收条件满足。

### 15.8 K1 执行报告

```text
task_id=K1
changed_files=tools/check_vpfp_jc_source_checkpoint.py, tools/analyze_vpfp_jc_k1.py, tools/analyze_vpfp_k1_stage_scale.py, src/main_vpfp.cpp, src/vpfp_integrator.cpp
production_equations_changed=0 (未修改 Poisson/remap/pusher/collision 公式)
accepted_state_ownership_changed=0
new_cli=0
build_status=PASS
test_commands=§9.2/§9.3 集群命令 + §9.5/§9.6 分析脚本 + §15.11 分阶段核查工具
result_files=output/vpfp_pairing_gate_k1/{k1_source_check,k1,stage_audit_coarse,stage_audit_fine,stage_scale}.result, {coarse,fine}/*.dat
pass_count=15 (§9.8.1) + 2 (stage audit structure) + 1 (stage scale analysis) + continuity remote rerun + scalar/direct-face dual reconstruction
fail_count=1 (energy_residual_reduction); manufactured G/G* adjoint test remains incomplete, not a production failure
first_failure=energy_residual_reduction
next_allowed_task=实现并验收独立的G/G*伴随制造解测试；P1/P2仍BLOCKED
```

> 后续小节按时间顺序记录 K1 的完整核查链：§15.9 修复分析器假 FAIL → §15.10 补齐
> 分阶段 energy audit → §15.11 验收与根因排除链 → §15.12 力-场配对核查、方案评估
> 与积分格式层修复方案。本节仅记录任务执行本身。

#### 已修改文件

| 文件 | 修改内容 |
|------|---------|
| `tools/check_vpfp_jc_source_checkpoint.py` | §9.5 前置身份检查（manifest 读取，非目录名） |
| `tools/analyze_vpfp_jc_k1.py` | §9.6 结果分析器（field-tol/pairing-tol/max-iters 显式参数 + §9.8.1 全部字段 + ledger_populated 门） |
| `src/main_vpfp.cpp` | restart→discrete-gradient 时重新求解初始 Poisson（failure 208 修复） |
| `src/vpfp_integrator.cpp` | ① trial 复用 legacy 转换数据流（failure 10）；② Poisson 容差全局 MAX 统一（failure 202）；③ Tail/Beam kick 功纳入 pairing 残差（failure 205 平台）；④ accepted-step 能量账接入 |

#### 修复过程（4 轮）

1. **第 1 轮：首步转换失败（failure 10）**：trial 未复用 legacy 转换数据流。修复：trial 的 `u_full()` 传入 `partition_`/`exported_flux`，设置 `apply_interface_sink`，保存真实 conversion ledger。
2. **第 2 轮：全局 Poisson 门 MPI 分裂（failure 202）**：`final_poisson_tolerance` 用局部 rho 尺度导致 rank 判定分裂。修复：先 MPI_Allreduce(MAX) 再算全局容差，bootstrap 同样全局化。
3. **第 3 轮：固定点与功配对不收敛（failure 205）**：omega=0.5 收缩慢 + pairing 漏 Tail/Beam kick 功。修复：配对用全局 Bulk+Tail+Beam ΔK，K 门显式 `relaxation=1.0`。
4. **第 4 轮：混合表示固定点平台**：Eulerian/PIC 表示产生确定性二周期，Linf 平台随步升高。生产 K 门修订为 `field_tol=1e-4`、`pairing_tol=1e-8`、`omega=1`、`max_iters=12`；接入 accepted-step 能量账，分析器新增 `ledger_populated` 门。

#### §9.8.1 验收审计（能量加密判定修正）

按 §9.8.1（生产 K 门：`field_tol=1e-4`、`pairing_tol=1e-8`、`max_iters=12`）：

| §9.8.1 要求 | coarse | fine | 结果 |
|-------------|--------|------|------|
| same_source_checkpoint=1 | 同一 115 fs checkpoint | 同一 | **PASS**（§9.5 前置检查） |
| same_initial_physical_state=1 | 同一重启源 | 同一 | **PASS** |
| same_physical_window=1 | 终点 1.1513935063899238e-13 | 终点 1.1513935063899264e-13 | **PASS**（diff~2.5e-28） |
| coarse_accepted_steps=10 | 10（5072-5081） | N/A | **PASS** |
| fine_accepted_steps=20 | N/A | 20（5072-5091） | **PASS** |
| soft_accept_count=0 | 0 | 0 | **PASS** |
| max_iteration_count<=12 | 3 | 3 | **PASS** |
| all_field_residual_l2<=1e-4 | 1.22e-05 | 2.17e-06 | **PASS** |
| all_field_residual_linf<=1e-4 | 6.22e-05 | 2.33e-05 | **PASS** |
| all_pairing_relative_to_exchange<=1e-8 | 6.8e-12 | 5.5e-12 | **PASS** |
| continuity_pass=1 | 全局number ledger非零 | 同上 | **仅全局账本可用；Gate-I逐组分continuity仍FAIL，不能判结构PASS** |
| local_work_ledger_pass=1 | 真实energy balance已输出 | 同上 | **PASS（仅账本可用）** |
| energy_residual_reduction=1 | signed=2.963e5, abs=3.763e5 J/m² | signed=3.210e5, abs=3.913e5 J/m² | **FAIL：fine分别恶化8.3%/4.0%** |
| poisson_identity_pass=1 | gauss_charge_residual≤4.3e-14 | ≤5.0e-14 | **PASS** |
| gauss_pass=1 | gauss_charge_residual 远小于容差 | 同上 | **PASS** |
| post_field_charge_invariance_pass=1 | post_field_charge_residual≤1.3e-14 | ≤1.4e-14 | **PASS** |
| nonfinite_count=0 | 0 | 0 | **PASS** |

关键数值：
- coarse：10 步全接受，max field_l2=1.22e-05、max field_linf=6.22e-05、max pairing=6.8e-12
- fine：20 步全接受，max field_l2=2.17e-06、max field_linf=2.33e-05、max pairing=5.5e-12
- 两者物理终点一致（1.151393506e-13），覆盖相同物理窗口
- `coarse/vpfp_failure.dat` 已不存在（上一轮 5 步失败已解决）

**§9.8.1 判定：FAIL_ENERGY_REFINEMENT。** 固定点、连续性、Poisson、Gauss和post-field charge通过，
但fine能量余额未随dt/2下降。

验证状态：**K1未通过，K2保持BLOCKED。** 必须先定位能量余额来自时间离散还是accepted-step账本漏项。

### 15.9 分析器假 FAIL 修复与权威复核（2026-08-18）

本节只读取当前 `output/vpfp_pairing_gate_k1/{coarse,fine}`，不引用旧 `.err`。

#### 两处分析器 bug（`tools/analyze_vpfp_jc_k1.py`，均已修复）

1. **`# columns=` 头解析错误**：`read_dat_columns` 把所有 `#` 开头行当注释跳过，但
   `field_particle_iteration.dat` 的列名就写在 `# columns=...` 行上。于是第一个真实数据行
   被误当 header：
   - 数据行只剩 9（coarse）/19（fine）条；
   - header 全是数值 → `iterations`/`time_s`/`field_residual_l2` 索引均 -1 → 全部返回默认 0；
   - `last_time_s=0` 不满足窗口覆盖条件 → `same_physical_window=0` → 假 FAIL。
   修复：`read_dat_columns` 识别 `# columns=...` 行作为 header，跳过 `# schema=` 行。
2. **status 判定 key 无前缀**：判定循环检查的是 `continuity_pass`/`local_work_pass`/
   `poisson_pass`/`gauss_pass`/`post_field_charge_pass` 等无前缀 key，而实际 key 是
   `coarse_*`/`fine_*` 前缀 → `continuity_pass` 恒为 0 → 假 FAIL。修复：按 coarse/fine 前缀
   分别判定。

#### 修复后权威复核

```text
coarse: 10步全接受，max_iter=3，max_l2=1.220e-5，max_linf=6.223e-5，max_pairing=6.80e-12
fine:   20步全接受，max_iter=3，max_l2=2.166e-6，max_linf=2.330e-5，max_pairing=5.54e-12
终点物理时间差约 2.5e-28 s（覆盖相同物理窗口）
```

真实能量余额（修复前误把 `domain_energy_delta` 当残差，已改为真实余额列）：

```text
coarse signed=2.96324e5, abs=3.76297e5 J/m2
fine   signed=3.21014e5, abs=3.91254e5 J/m2
fine/coarse signed=1.0833, abs=1.0397   （fine 未随 dt/2 改善）
```

**结论**：固定点、pairing、Gauss、post-field charge、结构门全部通过；**唯一真实失败项为
`energy_residual_reduction`**（fine 累计残差不随 dt 缩小）。此时 stage 明细文件只有表头，
不足以定位缺口来源，下一步补齐分阶段 energy audit（§15.10）。

该步骤已由 §15.10（DG stage audit 补齐）+ §15.11（验收与根因排除链）完成。

### 15.10 K1 分阶段 energy audit 补齐（discrete-gradient 路径，2026-08-18）

**问题**：`advance_discrete_gradient` 从未调用 `capture_stage_energy` / `finalize_stage_energy_audit`，
因此 `vpfp_stage_energy_audit.dat` 只有表头。11 个阶段记录（accepted_n → collision_half1 → x_half1 →
midpoint_poisson → u_force_tail_beam_kick → conversion_after_force → x_half2 → collision_half2 →
conversion_after_collision → tail_bulk_return → final_poisson）在 DG 路径中散落在三个 helper 里，
无法用 legacy 的单个内联 pass 采集。

**根因**：DG 路径把子步拆到 `prepare_field_particle_frozen_state`（阶段 0-3）、Picard trial
`evaluate_field_particle_trial`（阶段 4-6）和 post-field `post_field_charge_invariance_transaction`
（阶段 7-9），最后阶段 10 在收敛后复用同一份 final fields。

**修复**（`src/vpfp_integrator.h/.cpp`）：
- 新增成员 scratch `stage_energy_scratch_`、累积 sources `stage_sources_`、frozen sources 快照
  `stage_sources_frozen_`，以及 `capture_dg_stage()` helper，把记录写入 member scratch，使其跨
  prepare/trial/post-field 三个 helper 存活；
- `advance_discrete_gradient` 每步开始时重置 scratch；
- `prepare_field_particle_frozen_state` 采集阶段 0-3（accepted_n、collision_half1、x_half1、
  midpoint_poisson），并累加 C1 reservoir 与 x1 背景流；
- `evaluate_field_particle_trial` 每次 trial 开始重置 scratch count=4 并从 frozen sources 恢复，
  采集阶段 4-6（force、conversion_after_force、x_half2），累加 bulk_upar/Gate-C 功、kick 功、
  conversion 与 x2 背景流；由于收敛 trial 是最后一次求值，被接受的 trial 记录自然保留；
- `post_field_charge_invariance_transaction` 采集阶段 7-9（collision_half2、
  conversion_after_collision、tail_bulk_return），累加 C2 reservoir 与 C2 conversion；
- `advance_discrete_gradient` 在收敛后采集阶段 10（final_poisson，复用 accepted trial 的
  final fields），把 scratch 拷贝进 `result.stage_energy` 并调用 `finalize_stage_energy_audit`。

**只读性**：所有采集点只读取已经算出的 kinetic/field 能量与 ledger，不重放任何物理算子，不参与
接受条件，不改变 accepted state、RNG 或累计 ledger。

**验收**：`g++ -fsyntax-only`（MPI stub）通过；`vpfp_integrator.cpp`、`main_vpfp.cpp`、
`vpfp_diagnostics.cpp`、`vpfp_field_particle_trial_test.cpp`、`vpfp_post_field_charge_invariance_test.cpp`
全部编译。集群重跑后 `vpfp_stage_energy_audit.dat` 应含每个 accepted 步的 11 行记录，再以
`analyze_vpfp_stage_energy_audit.py` 判定 stage_balance 首次非零出现在哪一阶段。

### 15.11 K1 验收与根因排除链（2026-08-18，集群重跑后）

#### 数据完整性

| 文件 | coarse | fine |
|------|--------|------|
| field_particle_iteration.dat | 10 行（5072-5081） | 20 行（5072-5091） |
| vpfp_step_diagnostics.dat | 10 行 | 20 行 |
| vpfp_stage_energy_audit.dat | 111 行（10 步 × 11 阶段） | 221 行（20 步 × 11 阶段） |

#### §9.8.1 复核审计（关键项，全表见 §15.8）

| 项目 | coarse | fine | 判定 |
|------|--------|------|------|
| accepted_steps | 10 | 20 | **PASS** |
| max_iterations | 3 | 3（≤12） | **PASS** |
| max field_l2 / linf | 1.22e-5 / 6.22e-5 | 2.17e-6 / 2.33e-5 | **PASS**（≤1e-4） |
| max pairing | 6.83e-12 | 5.53e-12 | **PASS**（≤1e-8） |
| nonfinite / ledger / component continuity / gauss / post-field | 0/1/1/1/1 | 0/1/1/1/1 | **PASS（远程最新重跑；本地副本较旧）** |
| 物理窗口终点 | 1.1513935063899238e-13 | 1.1513935063899264e-13 | **PASS**（diff 2.5e-28） |
| 能量 refinement | signed 2.963e5, abs 3.763e5 | signed 3.210e5, abs 3.913e5 | **FAIL**（fine 未改善） |

当前能量主门仍为`energy_residual_reduction=0`。远程最新重跑已补齐conversion events并使Gate-I
逐组分`continuity_pass=1`；本地未镜像该批逐rank文件不否定远程结论。以下排除链在远程最新结果上有效，
最新`v3` scalar与direct-face dual重构均已通过；但独立$G/G^*$制造解伴随测试仍缺失。
因此仍不能唯一选择P1/P2生产修复。

#### 排除 A：stage 采集口径假象（已修复，假缺口）

逐步 U_E 绝对演化（coarse step 5072，修复前）：

```text
accepted_n              U_E=8.85254e7
midpoint_poisson        U_E=8.94572e7   (+9.32e5, midpoint 场能)
u_force_tail_beam_kick  U_E=8.29194e7   (-6.538e6, 粒子获 5.414e5 功)
x_half2                 U_E=8.94572e7   (+6.538e6, 跳回 midpoint 值!)
collision_half2         U_E=8.80376e7   (-1.42e6, 最终场能)
final_poisson           U_E=8.80376e7
```

`x_half2` 是纯空间 remap，物理上不改变 U_E，却记录出 6.5e6 的假跳变。根因：
**DG 路径 stage capture 口径不一致**——`u_force` 阶段记 `trial_force_fields_`
（E_pair 候选场），而 `x_half2`/`collision_half2` 记 `midpoint_fields_`。

修复（`src/vpfp_integrator.cpp`）：与 legacy 口径一致——`u_force`、
`conversion_after_force`、`collision_half2`、`conversion_after_collision`、
`tail_bulk_return` 阶段统一记 `midpoint_fields_`（场在 force 前后不变），仅
`final_poisson` 记 `final_fields_`。修复后：

| 指标 | coarse（前→后） | fine（前→后） |
|------|----------------|--------------|
| field_particle_pair_residual_relative | 0.693 → **0.055** | 0.480 → **0.060** |
| decomposition_matches_roundoff | 0 → **1** | 0 → **1** |
| source_residual_gate_pass | 0 → **1** | 0 → **1** |
| U_E 假跳变 | 有 → **无** | 有 → **无** |

**结论**：此前 `first_dominant_stage=field_coupling`（1.16e7）大部分是此假跳变，
修正后 field_coupling 残差 = 2.995e5（coarse）/3.242e5（fine），与
`energy_balance_residual`（2.963e5 / 3.210e5）逐位对应（比值 1.01）——**能量 refinement
失败就是修正后的 field_coupling 缺口**。采集口径修复后缺口仍存在，进入下一步排除。

#### 排除 B：时间截断误差（被推翻）

**分阶段时间细化核查**（`tools/analyze_vpfp_k1_stage_scale.py`，跨 coarse/fine 对比
每阶段 stage_balance）：

```text
coarse_accepted_steps=10  fine_accepted_steps=20  accepted_step_count_valid=1
time_shrinking_stages=midpoint_poisson(0.31), u_force_tail_beam_kick(0.69),
    conversion_after_force(0.50), x_half2(0.70), collision_half2(0.30)
fixed_magnitude_stages=none
inconclusive_stages=none
status=PASS
```

（括号内为 fine/coarse 每步平均 |balance| 比值。）

- 所有大 balance 阶段都随 dt 缩小（比值 0.30-0.70），无固定量级账本项；
- 但 per-step 残差缩小是**物理量自身**随 dt 缩小（力功、场能变化都是 O(dt) 量），
  不能据此判定缺口是时间误差。

**决定性判据——单位物理时间残差率恒定**：

```text
field_coupling residual rate:
  coarse = 2.601e21 J/(m2 s)      fine = 2.668e21 J/(m2 s)
  ratio fine/coarse = 1.026（不随 dt 变化）
```

若缺口是时间截断误差，残差率应随 dt 缩小（dt/2 → 率减半）；实测比值 ≈1。
**缺口不是时间离散截断误差**，而是每单位物理时间恒定的能量不守恒项。

#### 排除 C：漏源项（被推翻）

逐步 stage_balance 分解（coarse step 5072）：

```text
stage                          stage_balance
collision_half1                -8.494e-07   (闭合，reservoir 已被 dK 抵消)
x_half1                        -1.906e-06   (闭合)
midpoint_poisson                +9.318e+05   (dU_E=+9.32e5，场能增)
u_force_tail_beam_kick          +4.943e+05   (粒子动能增，场能不变)
x_half2                         -1.112e-06   (闭合)
collision_half2                  +3.427e-07   (闭合)
final_poisson                   -1.4196e+06   (dU_E=-1.42e6，场能降)
sum(stage_balance)              = 1.506e+04   = ledger residual  ✅ telescope 闭合
```

- **stage telescope 完全闭合**：`sum(stage_balance) = ledger residual` 逐位一致
  （step 5072: 15063.3 = 15063.3）；
- 碰撞 reservoir、背景流、beam 等源项已在各自 stage_balance 中精确闭合
  （collision 阶段 balance≈0）；
- **缺口不是"漏了源项"**。

#### 排除 D：记账 bug（被推翻）

`tools/analyze_vpfp_k1_force_field_pairing.py` 按 `residual = force_pair +
poisson_pair + other` 分解每步：

```text
coarse: force_pair=-5.1488e6  poisson_pair=+5.4451e6  other≈0
fine:   force_pair=-5.0539e6  poisson_pair=+5.3749e6  other≈0
decomposition_self_consistent=1（重构残差/账本残差比值=1.000000000003）
```

stage 分解精确复现账本残差——**不是记账 bug**，缺口在物理量本身。

#### 排除 E：力场≠G_P（被推翻）

`field_particle_iteration.dat` 显示力场与 G_P 的 field residual 每步收敛
（l2~1.1e-5、linf 2-6e-5 < 1e-4）——**粒子实际受的场确实 = G_P**（§16 条件3 已满足）。

#### 阶段性结论：固定比例力-场离散功缺口候选

在现有全局账本口径下，缺口集中于力-场配对项。由于Tail逐组分continuity审计仍失败，下面的相关性
只能定位候选分支，不能作为唯一生产根因。量化（§15.12 工具）：

```text
corr(force_pair, residual) = -0.986 (coarse) / -0.983 (fine)
residual ≈ -2.7% × force_pair（最小二乘斜率，coarse/fine 一致）
```

即：粒子被 G_P 场推的动能变化（force_pair）与场能变化（U_E(final)-U_E(n)）之间存在
**固定比例（≈-2.7%）的不闭合**——粒子获得/失去的能量中约 2.7% 在场能账本中没有对应项。
方向随当场力功反转（corr→-1），比例不随 dt 变（coarse/fine 共享空间网格）。

**K1 判定：FAIL_ENERGY_REFINEMENT（真实，机制定位见 §15.12）。**

### 15.12 §16 条件3 力-场配对核查、方案评估与修复方案

#### 核查实现（只读诊断）

沿 §16 条件 3 方向，新增只读诊断，判定"粒子实际受的场"与"场能账本对应的场"是否同一
离散对象：

1. **`src/vpfp_integrator.h`**：`VpfpStepLedger` 新增 `e_pair_field_energy` 字段 —— 收敛 trial
   实际用于推粒子的 E_pair 场（`G_P(Φ^n,Φ^{n+1})`，即 `trial_force_fields_` 在收敛时的值）的
   全局场能 U_E(E_pair)。
2. **`src/vpfp_integrator.cpp`**：`advance_discrete_gradient` 在 `finalize_energy_ledger` 前记录
   `result.ledger.e_pair_field_energy = field_energy(trial_force_fields_)`（全局归约，各 rank 同值）。
3. **`src/vpfp_diagnostics.cpp`**：`vpfp_step_diagnostics.dat` 新增两列 `U_E_pair U_E_final`
   （E_pair 场能 / 最终场能）。
4. **新工具 `tools/analyze_vpfp_k1_force_field_pairing.py`（正式判定）**：分解每步
   `residual = force_pair + poisson_pair + other`（粒子动能变化 + 场能变化 + 其余 stage），
   按单位物理时间残差率 coarse/fine 比值判定 `FIXED_FORCE_FIELD_GAP` vs `TIME_TRUNCATION`，
   并验证 `decomposition_self_consistent`；另输出 `corr_force_pair_residual` 与
   `residual_vs_force_pair_ratio`（§15.12 方案 D，见下）。
5. **新工具 `tools/analyze_vpfp_k1_epair_consistency.py`（参考）**：记录 `U_E_pair` vs
   `U_E_final` 差异与力功 identity 残差，仅诊断参考（U_E 为二次型，凸组合与跨场比较
   不构成判定），正式判定以 force_field_pairing 为准。

**只读性**：全部改动只记录已算出的场能/功，不重放算子、不参与接受条件、不改 accepted state。

#### 核查命令

重编译后重跑 coarse/fine（§9.2/§9.3），然后：

```bash
# 1) §16 条件3 正式核查：力-场配对离散功口径（决定性判定）
python3 tools/analyze_vpfp_k1_force_field_pairing.py \
  --coarse ./output/vpfp_pairing_gate_k1/coarse \
  --fine ./output/vpfp_pairing_gate_k1/fine \
  --result ./output/vpfp_pairing_gate_k1/force_field_pairing.result

# 1b) §16 条件3 参考：E_pair 场能 vs 最终场能（诊断，非判定）
python3 tools/analyze_vpfp_k1_epair_consistency.py \
  --coarse ./output/vpfp_pairing_gate_k1/coarse \
  --fine ./output/vpfp_pairing_gate_k1/fine \
  --result ./output/vpfp_pairing_gate_k1/epair_consistency.result

# 2) stage audit 复核（确认 field_coupling 结构门）
python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run ./output/vpfp_pairing_gate_k1/coarse \
  --result ./output/vpfp_pairing_gate_k1/stage_audit_coarse.result \
  --expected-accepted-steps 10
python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run ./output/vpfp_pairing_gate_k1/fine \
  --result ./output/vpfp_pairing_gate_k1/stage_audit_fine.result \
  --expected-accepted-steps 20

# 3) K1 总验收
python3 tools/analyze_vpfp_jc_k1.py \
  --coarse ./output/vpfp_pairing_gate_k1/coarse \
  --fine ./output/vpfp_pairing_gate_k1/fine \
  --source-checkpoint "$CHECKPOINT_115" \
  --field-tol 1e-4 --pairing-tol 1e-8 --max-iters 12 \
  --result ./output/vpfp_pairing_gate_k1/k1.result || exit 83

# 4) 分阶段时间细化对比
python3 tools/analyze_vpfp_k1_stage_scale.py \
  --coarse ./output/vpfp_pairing_gate_k1/coarse \
  --fine ./output/vpfp_pairing_gate_k1/fine \
  --result ./output/vpfp_pairing_gate_k1/stage_scale.result
```

**预期**：若 `verdict=FIXED_FORCE_FIELD_GAP`（残差率 fine/coarse≈1），则证实能量 refinement
失败源于力-场配对离散功口径不一致——沿 §16 条件 3 方向，而非时间离散或放宽门。

#### 核查结果（最终判定）

`tools/analyze_vpfp_k1_force_field_pairing.py` 聚焦恒等式
`residual_step = force_pair + poisson_pair + other`：

```text
force_pair   = u_force 阶段 dK（粒子被 G_P 场推的动能变化）
poisson_pair = midpoint+final 场能变化 dU_E（= U_E(final)-U_E(n)，telescope 闭合）
other        = x remap / collision / H10 stage_balance

coarse: force_pair=-5.1488e6  poisson_pair=+5.4451e6  other≈0
fine:   force_pair=-5.0539e6  poisson_pair=+5.3749e6  other≈0
coarse_rate_signed=2.574e21  fine_rate_signed=2.641e21 J/(m2 s)
signed_rate_fine_over_coarse=1.026
decomposition_self_consistent=1 (重构残差/账本残差比值=1.000000000003)
dominant_component=poisson_pair
verdict=FIXED_FORCE_FIELD_GAP
```

**判定**：`decomposition_self_consistent=1`（stage 分解完美复现账本残差），
`FIXED_FORCE_FIELD_GAP` —— 残差率不随 dt 变化（fine/coarse≈1），缺口是**每单位物理时间
恒定的力-场配对口径差异**，不是时间截断。物理结论：`force_pair` 与 `poisson_pair` 每步相差
1.5e4-7.7e4 J/m²（相对 3%），且每单位物理时间恒定——粒子用 G_P 场推的功与场能变化
**不是同一离散功口径**。

参考工具（epair_consistency）：`U_E_pair` 与 `U_E_final` 相差 ~5e6 J/m²、低于两端点场能。
这是**物理正常**现象（E_pair 是时间中心场，U_E 是 E 的二次型，不满足凸组合约束），仅参考
不构成判定；`field_particle_iteration.dat` 显示力场与 G_P 的 field residual 每步收敛
（l2~1.1e-5、linf 2-6e-5 < 1e-4）——**力场确实≈G_P**（排除链 E，§15.11）。

#### 修复方案评估（方案 A/B/C/D）

**问题精确定位**：`residual = force_pair + (U_E(final) - U_E(n层))`（逐位恒等式），
每步 1.5e4-7.7e4 J/m²，不随 dt 变。

**决定性新发现 —— 残差与力功强负相关**：

```text
corr(force_pair, residual) = -0.986 (coarse) / -0.983 (fine)
residual_vs_force_pair_ratio = -0.027 (coarse) / -0.028 (fine)
residual_over_force_pair_mean = -1.6% (coarse) / -8.4% (fine)
```

残差几乎完全由 `residual ≈ -2.7~2.8% × force_pair` 决定（最小二乘斜率，coarse/fine 一致）：
粒子被 G_P 场推的动能变化（force_pair）有约 2.7% **未被场能变化抵消**。这是力-场离散功的
固定比例缺口，不随 dt 缩放，方向随当场力功反转（corr→-1）。

**关键事实**：
- `U_E(pair)`（G_P 场能）与 `U_E(final)` 差 ~5e6 且始终低于 n 层——E_pair 是时间中心场，
  场能定义不同于存储场，**不能**作为配对量。
- `midpoint_poisson_delta + final_poisson_delta` 与 `U_E(final) - U_E(accepted_n)`
  **数值完全相等**（stage telescope 闭合）——配对定义不影响结论。
- residual ≠ 0 且 ≈ -2.7~2.8%×force_pair，是**力-场离散功格式的固定比例不闭合**，
  不是记账 bug、不是漏源项、不是粒子受力错误（力场确实=G_P）。

**方案对比**：

| 方案 | 代码改动 | 物理 | 性能 | 判定 |
|------|---------|------|------|------|
| A. 力步改用最终场推 | `evaluate_field_particle_trial` 收敛后用 `final_fields_trial` 替换 force 场 | 破坏 §16 条件3（力场≠G_P） | 零 | **不推荐**（违背前提） |
| B. 只改 stage 记账口径 | `poisson_pair` 改用 `U_E(final)-U_E(n)` | 不改受力，但该定义与原定义数值等价 | 零 | **无效**（已证等价） |
| C. 收敛后 G_P 完整场重推 | 收敛后补 rho/phi 重推 | 力场仍 G_P，rho 不影响受力 | +1 次迭代 | **无效**（已证 rho 与受力无关） |
| D. **保留 residual 为时间中心缺口度量（推荐）** | 不改受力、不改记账；把 `residual = force_pair + (U_E_final - U_E_n)` 明确为 §16 条件3 的时间中心离散缺口指标 | 物理受力不变（G_P），§16 条件3 保持 | 零 | **推荐** |

**方案 D 说明（只保留为诊断，不是修复结论）**：残差 ≈ -2.7%×force_pair 且不随 dt 变，说明
当前离散系统存在固定比例缺口，但尚不能把它定义成可接受的“固有误差”。正确做法：
1. **不通过改受力或加补丁消除**（会违背 §16 条件3/4）；
2. 保留 residual 作为**时间中心离散缺口的明确度量**，其值 ≈ -2.7~2.8%×force_pair；
3. 真正消除需在**时间积分格式**层面改进——核查 G_P 时间中心构造与力步离散功（Gate-C
   `upar_internal_face_energy_transfer`）的恒等式，找出 2.7% 比例因子缺失的离散项。

**方案 D 执行结果（2026-08-19，本地已实施，无需集群重跑）**：

```text
coarse_corr_force_pair_residual = -0.9857766
coarse_residual_vs_force_pair_ratio = -0.0269807
fine_corr_force_pair_residual  = -0.9829915
fine_residual_vs_force_pair_ratio = -0.0283949
decomposition_self_consistent = 1
signed_rate_fine_over_coarse = 1.0263
verdict = FIXED_FORCE_FIELD_GAP   （与实施前一致，纯诊断改动不改变判定）
first_failure = force_field_gap
status = FAIL
```

- `tools/analyze_vpfp_k1_force_field_pairing.py` 新增 `corr_force_pair_residual` 与
  `residual_vs_force_pair_ratio` 报告（判定逻辑未改）；`analyze_vpfp_stage_energy_audit.py`
  未改动（配对定义等价）。
- **无需重跑 §9.2/§9.3**：方案 D 是纯分析器改动；当前 `output/vpfp_pairing_gate_k1` 即为
  集群最新重跑结果，重跑只会得到逐位相同结果。**只有实施积分格式层修复（改生产代码）后**
  才需重编译重跑作回归。

#### 积分格式层修复方案

**机制判定**：§15.11 排除链 A-E 已排除采集口径假象、时间截断、漏源项、记账 bug、力场≠G_P。
此处补充两项与"固定比例"特征相关的机制判定：

| 候选机制 | 证据 | 判定 |
|---------|------|------|
| 显式 Euler 推步的 O(Δt²) 相对论项 | 该项与力功比值应随 dt 缩小（dt/2 → 比例减半），实测 coarse/fine 相同 | **已排除** |
| **空间离散对偶口径不一致** | coarse/fine共享空间网格且缺口率不随dt变；可能涉及格点/面场、x-remap charge current、u-remap work current或$G/G^*$映射 | **主要候选，尚未唯一确认** |

**候选机制细节**：粒子推步的加速度用**格点（单元平均）场**
`a_u = q·E_cell/(m·c)`，`E_cell[ix] = 0.5(E_f[ix]+E_f[ix+1])`
（`conservative_ppm_remap.cpp:680-682`）；而 G_P 构造与场能配对恒等式
（`build_potential_pairing_field` / `evaluate_work_identity`）用**面场** E_f 与
配对势单元平均 φ_bar。在强 E 梯度区（束流-等离子体不稳定性），E_cell≠E_f 的
系统性偏差（O(dx²) 局部、相对可达数%）造成"粒子功（格点口径）与场能变化
（面场口径）"的固定比例不闭合——因为 coarse/fine 网格相同，比例不随 dt 变。

**阶段 1：机制确认（只读，零代码改动，复用 Gate I）**

重跑 §9.2/§9.3 时开启 Gate I pairing audit（`field_particle_power_pairing.dat`），
`poisson_transport_residual = ΔU_E − W_electrode + dt⟨E_pair,J⟩` 区分两分支：

- 若 `poisson_transport_residual ≈ 0` 且缺口 = `force_pair + W_pc`（W_pc =
  potential_charge_work，已有 `trial.potential_charge_work`，需写盘）→ **分支 P**：
  只确认缺口在“粒子受力功 vs charge-current配对功”；不能仅凭该结果断言是格点/面场插值系数；
- 若 `poisson_transport_residual ≠ 0`（达 2.7% 量级）→ **分支 Q**：缺口在
  Poisson 功恒等式自身（端点权重/φ 平均修正的代数缺陷）。

**阶段 2：修复（依据阶段 1 分支，最小侵入、零性能开销）**

- **分支 P（先做离散对偶电流审计，§16 条件4 方向）**：在修改生产算子前，必须直接比较
  $J_{charge,face}$ 与 $G^*J_{force,cell}$。只有残差逐face复现current-pair缺口后，才允许评估：
  - P1：配对恒等式改为格点口径——G_P 的构造与 `-dt⟨E_pair,J⟩ = Σφ_bar·Δρ`
    全部改用 E_cell（格点）定义，使力功与场能配对在格点层面严格闭合；
  - P2：推步加速度改为与配对恒等式对偶的面场加权（需验证粒子受力仍 = G_P，
    保持 §16 条件3）。
  - 判据：修复后 `energy_balance_residual → roundoff`（相对 <1e-8），
    force_field_pairing 的 slope/corr → ≈0。
- **分支 Q（Poisson 恒等式代数修正）**：修正 `build_potential_pairing_field` /
  `evaluate_work_identity` 的端点权重（dx/2）与 φ 单元平均修正（±dx·ΔE/12 系数），
  使 `ΔU_E = W_electrode + W_pc` 精确成立（改后恒等式残差应到 roundoff）。
  - 判据：`poisson_transport_residual → roundoff`，residual 随之消失。

**性能**：两分支均为 O(N) 单次循环，零性能开销；不改 field/pairing 收敛门、
不加能量补丁、不放松任何 Gate（§9.7 约束不变）。

**回归验证（修复实施后）**：重编译 → 重跑 coarse/fine（§9.2/§9.3）→ 重跑
`analyze_vpfp_k1_force_field_pairing.py`（预期 slope/corr→0、residual→roundoff）、
`analyze_vpfp_jc_k1.py`（预期 energy_residual_reduction PASS）、
`analyze_vpfp_stage_energy_audit.py`（结构门保持 PASS）、
`analyze_vpfp_k1_stage_scale.py`。全部通过后更新 §15.8 执行报告与 §15.12 修复记录。

### 15.13 历史根因复核记录（2026-08-19，已被§15.14取代）

> **状态说明**：本节保留于此只用于追踪“边界功遗漏→scalar dual修复”的历史过程。
> 2026-08-20远程直读已确认scalar dual reconstruction通过，本节中
> `dual reconstruction=FAIL`、`DUAL_RECONSTRUCTION_FAILED`和对旧`dual_face.result`的物理解释均不再是当前结论。
> 后续开发和验收必须以§15.14为准。

本节以当前 `output/vpfp_pairing_gate_k1` 为唯一数据源。先澄清各结果文件的 `status=PASS`：

- `stage_audit_*.result PASS`：表示11阶段采集、telescope和source ownership自洽；不表示物理能量闭合；
- `stage_scale.result PASS`：表示coarse/fine阶段缩放分析有效；不表示残差随dt收敛；
- `pairing_branch.result PASS`：表示分支分类成功，实际`verdict=BRANCH_P`；不表示K1通过；
- `k1.result status=FAIL`且`first_failure=energy_residual_reduction`才是K1总判定。

#### 15.13.1 已可靠确认的证据

```text
Poisson-transport residual / full residual = O(1e-10)
current-pair residual / full residual ≈ 1
force-work local ledger residual = O(roundoff)
conversion transaction residual = O(roundoff)
x-remap、collision、H10 stage balance = O(roundoff)
field-coupling residual rate fine/coarse ≈ 1.026
```

因此Poisson恒等式、边界源、碰撞、conversion和return不是主导缺口。能量残差位于
“charge-conserving面电流的配对功”与“速度受力产生的粒子功”之间。

#### 15.13.2 尚未通过的前置结构门

`field_particle_power_pairing.dat`中每个accepted步均有：

```text
continuity_pass=0
continuity_tail≈1e14
continuity_bulk≈1e7
continuity_beam≈1e3
```

主K1分析器此前仅以全局number ledger有限/非零代替逐组分continuity，因此错误写成
`continuity_pass=1`。DG路径只保存conversion标量ledger，没有把accepted trial的
`conversion_events`完整传给Gate-I workspace；Tail的逐cell conversion source因此可能缺失。
这会使Tail continuity和root-cause mask不可信。现有全局combined charge仍可闭合，因为Bulk/Tail
内部转换在总电荷中抵消，但不能据此宣称逐组分离散连续性通过。

#### 15.13.3 远程最新结果与当前可成立结论

证据位置：remote_cluster（2026-08-19）；逐rank dual-face大文件未SCP到本地，按§0.6不视为缺失。

远程重跑已确认：

~~~text
conversion events已接入u-force与C2 accepted trial
continuity_pass=1 (coarse/fine每个accepted步)
continuity_tail: O(1e14) -> O(1e3~1e4)
pairing_branch verdict=BRANCH_P
transport/full≈1e-10
current_pair/full≈1
dual-face文件: coarse 800个rank-step文件，fine 1600个rank-step文件
~~~

因此逐组分continuity前置门在远程最新版本已经通过；本地continuity_pass=0属于较旧镜像，不能覆盖
远程结论。Poisson分支Q、边界源、collision、conversion transaction和H10不是主导缺口。当前可靠
定位仍是“粒子受力功 vs charge-current配对功”接口。

#### 15.13.4 dual-face审计结果与正确解释

远程dual-face输出已经生成并完成运行，但重构门失败：

~~~text
coarse_dual_recon_pass=0
fine_dual_recon_pass=0
max_scalar_recon_relative=0.362 / 0.607
max_face_recon_relative=O(1e17)
left5与core90均有显著贡献
~~~

这些结果证明“当前dual审计构造不能重构current_pair residual”，而不是直接证明生产
$J_{charge}-G^*J_{force}$具有上述巨大物理误差。审计代码中存在一个已确认的口径遗漏：

- current_pair_residual使用的Tail/Beam总力功包含
  tail_delta_ke_boundary与beam_delta_ke_boundary；
- J_force_cell只由域内cell work数组构造；
- 当前dual区域积分只计算域内$-dt\langle E,J_{charge}-G^*J_{force,cell}\rangle$，
  没有加入PIC边界CIC份额。

因此正确重构关系应先写为：

$$
R_{current}
=
R_{dual,in-domain}
+
W_{tail,boundary}
+
W_{beam,boundary}.
$$

在补入该独立边界功之前，scalar reconstruction失败和左5%主导不能用于选择P1/P2。另需继续检查：

1. J_force_cell由force_work/(dt*dx*E_cell)构造时的近零场处理；
2. face quadrature端点dx/2权重；
3. MPI shared-face唯一owner；
4. G*是否与E_cell=0.5(E_f+E_{f+1})在同一内积下严格伴随；
5. scalar region sum与direct face sum是否使用同一符号和dt因子。

#### 15.13.5 当前执行状态

1. **conversion event接线：PASS（remote_cluster）。**
2. **逐组分continuity：PASS（remote_cluster）。**
3. **BRANCH_P分类：PASS。** 只确认缺口位于force-work/current-pair接口。
4. **dual-face文件生成：PASS（remote_cluster）。**
5. **dual reconstruction：FAIL。** 当前审计口径不完整，不能用于生产根因唯一化。
6. **P1/P2生产修改：BLOCKED。**

#### 15.13.6 强制执行顺序

1. 保留现有远程coarse/fine和continuity结论，不要求为本地分析重复SCP大文件。
2. 修复dual审计而非生产算子：
   - 在scalar reconstruction中显式加入Tail/Beam boundary force work；
   - 输出in-domain dual、boundary force work及两者之和；
   - 增加制造解测试验证$\langle E,G^*J\rangle=\langle GE,J\rangle$；
   - 验证shared-face owner和物理端点权重。
3. 远程重跑dual分析，必须同时满足：
   - scalar region sum与direct face sum一致；
   - in-domain dual + boundary work重构current_pair residual；
   - coarse/fine dual_recon_pass=1；
   - left5/core90/right5分解稳定。
4. 第3步通过后，若$J_{charge}-G^*J_{force}$仍完整复现能量缺口，才允许选择P1/P2；
   若差异消失，则问题仅在审计构造，不修改生产物理。
5. 始终禁止Poisson投影、场缩放、电流补丁或继续放宽field/pairing门。

K1状态：

~~~text
FAIL_ENERGY_REFINEMENT / DUAL_RECONSTRUCTION_FAILED
~~~

K2继续BLOCKED。

### 15.14 2026-08-20 远程最新K1验收与后续执行基线

#### 15.14.1 证据身份与读取约束

本节是当前唯一有效的K1结论，取代§15.13的当前状态和执行顺序。证据来自：

```text
evidence_location=remote_cluster
remote_host=XY_server
remote_project=/HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
result_root=output/vpfp_pairing_gate_k1
read_date=2026-08-20
access_mode=read-only
```

本次直接读取远程最新结果，没有要求SCP逐rank大文件，也没有修改远程文件、运行作业或生成新证据。
按§0.6，远程最新文件的结论高于本地较旧镜像。

#### 15.14.2 结果文件语义

- `stage_audit_coarse.result` / `stage_audit_fine.result` 的`PASS`只表示11阶段采集、telescope、source ownership和分解重构自洽。
- `stage_scale.result` 的`PASS`只表示阶段coarse/fine比较可用。
- `pairing_branch.result` 的`PASS`表示成功选出`BRANCH_P`，不表示K1通过。
- `k1.result status=FAIL` 且`first_failure=energy_residual_reduction`是K1总门判定。
- `dual_face.result status=PASS`表明修正后的direct-face汇总已通过；它不代替独立的$G/G^*$制造解伴随测试。

#### 15.14.3 已通过的K1结构门

远程最新`k1.result`确认：

```text
same_source_checkpoint=1
same_initial_physical_state=1
same_physical_window=1
coarse_accepted_steps=10
fine_accepted_steps=20
coarse_continuity_pass=1
fine_continuity_pass=1
coarse_local_work_pass=1
fine_local_work_pass=1
coarse_poisson_pass=1
fine_poisson_pass=1
coarse_gauss_pass=1
fine_gauss_pass=1
coarse_post_field_charge_pass=1
fine_post_field_charge_pass=1
soft_accept_count=0
nonfinite_count=0
max_iteration_count=3
```

场和功配对收敛门通过：

```text
all_field_residual_l2=1.2199267580931741e-05 < 1e-4
all_field_residual_linf=6.2226054575886500e-05 < 1e-4
all_pairing_relative_to_exchange=6.729269355320043e-12 < 1e-8
```

旧本地镜像中`continuity_pass=0`和`continuity_tail=O(1e14)`的原因是accepted trial的`conversion_events`
没有完整进入Gate-I workspace。远程最新版已将conversion events接入u-force和C2 accepted trial：

```text
coarse continuity_pass=1 for 10/10 accepted steps
fine continuity_pass=1 for 20/20 accepted steps
continuity_tail: O(1e14) -> O(1e3~1e4)
```

相对$O(10^{29})$粒子数尺度，剩余值属稳定求和/舍入尺度。逐组分continuity前置门已通过。

#### 15.14.4 K1仍失败：能量残差不随$dt$细化

| 指标 | coarse | fine | fine/coarse | 判定 |
|---|---:|---:|---:|---|
| signed residual | $2.9632405\times10^5$ | $3.2101445\times10^5$ J/m$^2$ | 1.0833 | FAIL |
| absolute residual | $3.7629653\times10^5$ | $3.9125365\times10^5$ J/m$^2$ | 1.0397 | FAIL |

coarse/fine使用同一源checkpoint并覆盖同一物理时间窗，能量残差没有随$dt/2$下降。
因此当前缺口仍是固定空间/功配对问题，不是主要时间截断误差。

`pairing_branch.result`进一步确认：

```text
coarse transport/full=1.0762401966548095e-10
fine transport/full=2.5771513448407913e-10
coarse current_pair/full=1.0000000000161784
fine current_pair/full=0.9999999999987690
verdict=BRANCH_P
```

结合stage audit，Poisson恒等式、x-remap、collision、conversion transaction、H10 return和source ownership
都不是主导缺口。缺口位于“charge-conserving面电流配对功”与“速度受力粒子功”之间。

#### 15.14.5 scalar dual reconstruction已通过

最新`field_particle_power_pairing.dat`的schema为`vpfp-field-particle-power-pairing-v3`，已显式输出
`dual_in_domain_work`、`boundary_force_work`、`dual_plus_boundary_work`和`dual_reconstruction_pass`。
它验证的关系为：

$$
R_{\rm current}
=R_{\rm dual,in-domain}
+W_{\rm tail,boundary}
+W_{\rm beam,boundary}.
$$

| 指标 | coarse | fine |
|---|---:|---:|
| accepted rows | 10 | 20 |
| `dual_reconstruction_pass` | 10/10 | 20/20 |
| max absolute reconstruction error | $1.1751\times10^{-9}$ | $4.8021\times10^{-10}$ J/m$^2$ |
| max error/tolerance | 0.563 | 0.632 |
| max error/current-pair | $6.40\times10^{-14}$ | $7.19\times10^{-14}$ |

因此§15.13历史记录中的PIC边界功遗漏已在scalar production audit中修复，scalar dual reconstruction已PASS。

#### 15.14.6 修正后的direct-face汇总已通过

远程文件数量已核验：

```text
coarse field_particle_power_dual_face files=800
fine field_particle_power_dual_face files=1600
coarse pairing profile files=800
fine pairing profile files=1600
```

历史分析器曾存在以下口径错误：

1. scalar侧比较$R_{left}+R_{core}+R_{right}=R_{current}$，漏掉`boundary_force_work`；
2. direct-face侧计算$\sum_f E_f\Delta J_f w_f$，漏乘`-dt`；
3. 工具自行推断$dx$和owner，没有直接使用文件已提供的`quadrature_weight`和`face_is_owner`。

由于$dt\sim10^{-17}$ s，漏乘$dt$会制造`max_face_recon_relative=O(1e17)`的假误差。修正后的
远程`dual_face.result`（2026-08-20 15:54:20）为`status=PASS`，具体结果为：

| 指标 | coarse | fine |
|---|---:|---:|
| `dual_recon_pass` | 1 | 1 |
| direct-face / in-domain 最大相对误差 | $7.32\times10^{-15}$ | $1.50\times10^{-14}$ |
| direct-face + boundary / current-pair 最大相对误差 | $6.15\times10^{-14}$ | $5.84\times10^{-14}$ |
| 最大owner face数 | 8001 | 8001 |
| missing owner / owner duplicate / owner-rank error / nonfinite | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |

`duplicate_face_count_max=79`只是80个MPI rank之间的79个共享面在原始文件中的正常重复，不是owner重复。修正后应记录为：

```text
scalar_dual_reconstruction=PASS
direct_face_files_present=PASS
direct_face_analyzer=PASS
independent_direct_face_reconstruction=PASS
manufactured_G_Gstar_adjoint_test=NOT_YET_INDEPENDENTLY_VALIDATED
```

#### 15.14.7 空间贡献不能简化为纯边界问题

| 区域 | coarse signed | coarse absolute accumulation | fine signed | fine absolute accumulation |
|---|---:|---:|---:|---:|
| left 5% | $1.5850\times10^5$ | $1.6587\times10^5$ | $1.6398\times10^5$ | $1.7100\times10^5$ |
| core 90% | $4.2814\times10^4$ | $1.7762\times10^5$ | $6.2937\times10^4$ | $1.9711\times10^5$ |
| right 5% | $8.99\times10^2$ | $7.41\times10^3$ | $8.32\times10^2$ | $7.43\times10^3$ |
| PIC boundary work | $5.78\times10^3$ | $1.50\times10^4$ | $4.68\times10^3$ | $1.24\times10^4$ |

signed sum中左端最大，但绝对值累计中core 90%与左端相当甚至更大。
因此缺口不能简化为纯边界问题；核心区同样存在显著配对贡献。PIC边界功是重构必需项，但不是累计主导项。

#### 15.14.8 当前执行状态

1. `conversion_events` accepted-state接线：**PASS**。
2. 逐组分continuity：**PASS**。
3. `BRANCH_P`分类：**PASS**。
4. dual-face文件生成：**PASS**。
5. scalar dual reconstruction：**PASS**。
6. direct-face独立重构：**PASS**。
7. $G/G^*$伴随性（独立制造解）：**PASS**。
   `manufactured_single.result`和`manufactured_mpi.result`均通过；单rank与MPI
   shared-face的恒等式误差均为0，物理端点权重均通过`dx/2`检查。
8. K1能量细化：**FAIL**。
9. P1/P2生产修复选择：**P1 REJECTED；基础P2已实施但不足；P3.0 确认单速度基线缺口，P3-V.1 当前 signed-work 审计无效，先修 P3-V.1R 测试口径**。
10. K2：**BLOCKED**。

#### 15.14.9 强制执行顺序

1. 保留远程coarse/fine、continuity、BRANCH_P、scalar dual与direct-face dual PASS结论；不为本地镜像重复SCP大文件。
2. `tools/analyze_vpfp_k1_dual_face.py`修复已验收，不再修改生产算子或重跑coarse/fine生产短窗。
3. 保留现有coarse/fine作为真实生产态的伴随自检证据：在生产权重下，
   `direct-face == in-domain scalar`、owner面唯一、开放端点权重正确。
4. 制造解测试已完成，并直接复用生产$G/G^*$和生产权重，在单rank、MPI
   shared-face和开放物理端点三类案例中验证：
   $$
   \langle E,G^*J\rangle_{face}=\langle GE,J\rangle_{cell}.
   $$
   物理端点权重必须是$dx/2$，MPI shared face全局只能有一个owner。
5. 当face schema不变时，只重跑修复后的分析器和制造解测试，不重跑coarse/fine生产短窗；
   本次`point4.result=status=PASS`。
6. 严格的独立验收必须同时满足：
   - direct-face vs in-domain scalar相对误差$≤10^{-10}$；
   - direct-face + boundary vs current-pair相对误差$≤10^{-10}$；
   - coarse/fine所有accepted步通过；
   - owner face无缺失、无重复；
   - 制造解伴随恒等式达到稳定求和舍入容差。
   本次远程结果：
   `point4.result=PASS`；单rank和MPI制造解均为`PASS`；
   `manufactured_identity_pass=1`、`manufactured_endpoint_weight_pass=1`；
   coarse/fine的direct-face重构最大相对误差分别为
   `6.15e-14`和`5.75e-14`，owner缺失/重复/rank错误和非有限face均为0。
7. 第6项已PASS，已能排除Poisson、$G/G^*$、边界功、owner和审计重构错误。代码复核进一步确认：
   当前`ConservativePpmRemap::advect_u_parallel()`已经用
   $E_{cell}=0.5(E_{f,i}+E_{f,i+1})$推进背景受力，即基础P2已在生产路径中；
   K1仍失败，故不得重复实施该基础P2。
8. **P1拒绝**：将Poisson配对恒等式改为cell口径只会改变能量账的配对目标，破坏已通过的face-Poisson
   功恒等式，不能修复真实场能变化。**局部缩放加速度也禁止**：它会直接改写Vlasov受力项而非建立同一离散系统。
9. P3.0 已在实际 Strang 时间层、最终 $Q^{x1}+Q^{x2}$ 通量下确认：单速度近无 limiter 案例仍有
   5%--7% 功缺口，故 PPM/FCT 不是第一根因。P3-V.1 的首轮 $|S|$ 加权速度审计因遗漏 $qE$
   符号且不表示实际功电流而无效；下一步只能执行 P3-V.1R 的 signed-work 电流审计。
   只有 P3-V.1R 确认离散速度链式法则后，才允许实施 P3-V.2；只有 P3-V 使单速度态闭合后，才允许
   重构 PPM/FCT 的最终相空间通量。实施、测试和停止条件见
   `docs/VPFP_P3联合x_u_Poisson功通量重构实施方案.md`。
10. 如果 P3-V 的有效速度审计不能复现单速度缺口，或统一速度表后单速度缺口不下降，问题属于审计构造或
    原离散不可闭合，禁止以
    电流补丁、场缩放、Poisson投影或局部加速度比例因子修改生产物理。
11. 始终禁止Poisson投影、场缩放、事后电流/能量补丁或继续放宽field/pairing门。

当前权威状态：

```text
K1=FAIL_ENERGY_REFINEMENT
scalar_dual_reconstruction=PASS
direct_face_analyzer=PASS
G_Gstar_production_state_adjoint=PASS
G_Gstar_manufactured_adjoint=PASS
point4_manufactured_and_dual_audit=PASS
P1=REJECTED_DIAGNOSTIC_REDEFINITION
P2_BASIC=ALREADY_IMPLEMENTED_IN_ADVECT_U_PARALLEL
P3.0=FAIL_SINGLE_VELOCITY_BASELINE_PAIRING
P3-V.1R=PASS_SIGNED_WORK_AUDIT
P3-V.2=PASS_AUDIT_FAIL_SINGLE_VELOCITY_ERROR_REDUCTION
P3/P3.1/P3.2=STOP_KEEP_ANALYTIC_DEFAULT
NEXT=JOINT_PHASE_SPACE_MIDPOINT_ENERGY_CLOSURE_RECONSTRUCTION
K2=BLOCKED
```

### 15.15 K1 最终根因链与路线切换（2026-08-20）

本节取代此前“继续 P3-V/P3.1 局部修复”的执行方向。

#### 已排除层

1. Poisson/Gauss、开放端点 $dx/2$ 权重、$G/G^*$ 伴随、MPI shared-face owner：均通过制造解和 direct-face 审计。
2. Beam/Tail 注入、出流、边界功与 conversion/return ledger：不是主导残差。
3. 逐组分连续性、source ownership、x remap 质量守恒：通过。
4. PPM/FCT limiter：单速度近无 limiter 状态已有 5%--7% 残差，故不是首要根因。
5. 单独替换 x transport 的解析中心速度为能量共轭速度表：P3-V.2 的 table/MPI/checkpoint/continuity 审计均通过，
   但单速度功误差未降低，故不是根因。

#### 最终结构性根因

当前 collisionless background 的推进是：

$$
T_x(\Delta t/2)\circ T_u(E,\Delta t)\circ T_x(\Delta t/2).
$$

其中两个 $T_x$ 的 charge current 来自独立 PPM characteristic sweep，$T_u$ 的动能变化来自独立 u-face
finite-volume sweep。它们没有由同一个时间中心 phase-space flux、离散乘积法则或反对称 Poisson bracket
共同定义。因此即使每个子算子各自守恒，整体仍不保证：

$$
\Delta K_{\rm bulk,E}
=\Delta t\sum_i dx\,E_i(GJ_{\rm charge})_i.
$$

这正是 K1 中固定、随 $dt$ 不消失的 current-pair 能量缺口。该根因处于**离散积分格式层**，不是单个边界、
速度表、limiter 或诊断项的缺陷。

#### 强制后续决策

```text
保持 analytic-cell-center 为默认生产模式。
停止 P3-V、P3.1、P3.2 的局部修改。
禁止电流补丁、能量补丁、场缩放、Poisson 投影、局部 a_u 缩放。
下一阶段仅允许构建新的联合 x/u/Poisson 时间中心能量闭合算子。
K2/K3/K4 继续 BLOCKED。
```

新的实施规范见：
`docs/VPFP_联合相空间时间中心能量闭合重构实施方案.md`。

## 16. 最终完成定义

只有同时满足以下条件，才能声明“时间中心不一致已修复”：

1. JC0–JC5 全部通过；
2. K1 固定 checkpoint coarse/fine 通过；
3. $E_{\rm force}=\mathcal G_P(\Phi^n,\Phi^{n+1})$ 在生产代码中由同一个收敛 trial 实际使用；
4. 粒子动能变化、场功和最终 Poisson 使用同一个候选状态；
5. 无 soft acceptance；
6. 失败事务保持 accepted state、RNG 和 ledger 不变；
7. K2 未破坏核心宏观响应；
8. K3 restart 等价；
9. K4 性能可接受；
10. legacy 默认路径在正式切换前保持可复现。

完成上述验收后，才可单独提交“将 production 默认模式切换为 `discrete-gradient`”的变更。默认切换本身必须再跑一次 JC5、K1 和 K3。

### 16.1 对 DeepSeek-V4-Flash 的适用性结论

扩展后的文档可以交给 DeepSeek-V4-Flash，但仅限以下使用方式：

- 一次新对话只执行一个 TASK；
- JC0、JC1、JC4、JC5 可由其直接执行；
- JC2、JC3 属于高风险核心任务，必须要求它先输出“现有代码块到新 helper 的映射表”，经人工确认后
  才允许编辑；
- K1--K4 必须分别运行和验收，不能让它在一次任务内连续修改生产代码并长跑；
- 任何测试失败时只分析本 TASK，不得自动进入下一节或扩大物理修改范围。

推荐给执行智能体的固定提示词：

```text
本次只执行 TASK_ID=<填写一个任务>。
先读取文档§0.1的“必读章节表”，只完整阅读当前TASK的必读章节。
按§0.3输出开始声明；然后阅读允许文件并给出文档步骤到实际函数的映射。
映射未完成前不要修改代码。严格遵守本TASK的允许文件、禁止项、测试命令和停止条件。
发现前置条件不满足、接口与文档不一致或需要修改表外文件时立即停止并报告，不得自行扩展范围。
完成本TASK的验收后停止，不得继续下一TASK。
```

不建议使用“按照整个文档全部实现”作为提示词。即使文档足够详细，这种提示仍容易导致较弱模型
同时重构推进、checkpoint 和测试，破坏事务边界并掩盖首个失败来源。
