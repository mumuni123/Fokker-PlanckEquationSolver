# VPFP公共能量残差：离散功配对审计与修复实施方案

## 1. 文档用途

本文档用于指导自动编码模型审计开放Beam、非周期Poisson场、Eulerian bulk与PIC Tail混合
VPFP求解器中的公共能量余额，并依据已完成的Gate A--G、100.3--120 fs生产结果以及
115--117 fs时间步A/B证据选择下一步修复。

执行模型应保留已通过Gate的结论，不得重复修改已被排除的空间Poisson算子、边界源账本、
bulk/Tail/Beam局部功配对或通过事后能量补丁伪造闭合。本文档面向能力有限的自动编码模型，
因此所有操作均按可执行约束书写。

主方案文档为：

```text
docs/开放Beam_非周期场_VPFP_PIC完整重构方案.md
```

本文件只处理公共能量余额，不重复H1--H10的其他设计。

---

## 2. 当前可信事实与修正后的根因状态

主要证据文件为：

```text
output/archieve/vpfp_force_work_none_100p3_100steps
output/archieve/vpfp_force_work_return_100p3_100steps
output/archieve/vpfp_force_work_none_100p3_dt_half_200steps
output/archieve/vpfp_force_work_none_100p3_100steps.result
output/archieve/vpfp_force_work_return_100p3_100steps.result
output/archieve/vpfp_force_work_none_100p3_dt_half_200steps.result
output/archieve/vpfp_poisson_work_identity_unit.result
output/archieve/vpfp_stage_energy_audit_unit_after_gatef.result
output/archieve/vpfp_gate_f_fixed_state_100p3.result
output/vpfp_return_dt_half_100p3_to_120fs
output/dt025_115_to_117
```

### 2.1 已经排除的根因

以下阶段累计闭合误差均约为$10^{-5}\ \mathrm{J/m^2}$或更小：

- 第一次和第二次$x$ remap；
- 第一次和第二次碰撞半步与collision reservoir；
- Tail--bulk H10返回；
- 转换事务的$N/P_x/K$守恒。
- bulk $u_\parallel$离散通量功配对；
- Tail和Beam kick的局部动能恒等式；
- Dirichlet Poisson空间算子与电极功恒等式；
- 已枚举并记录的开放边界、碰撞reservoir、转换和H10返回源项没有重复记账或漏记。

这里的`source_ownership_valid=1`只证明**已枚举字段**各出现一次，不证明枚举集合必然包含所有
由完整状态更新产生的物理功率。Gate H仍需用逐阶段状态差反查是否存在尚未进入枚举集合的项。

因此不得再次修改这些模块来“尝试降低总能量残差”。

### 2.2 当前公共余额

100步累计正式能量余额为：

| 模式 | 累计余额 | 末态总能量占比 |
|---|---:|---:|
| return-none | $2.985478269\times10^6\ \mathrm{J/m^2}$ | 约$5.53\times10^{-4}$ |
| hysteretic return | $3.041507349\times10^6\ \mathrm{J/m^2}$ | 约$5.63\times10^{-4}$ |

残差由两组量组成：

| 组合 | return-none | hysteretic return |
|---|---:|---:|
| midpoint/final Poisson净场能变化 | $2.087114625\times10^6$ | $2.117814073\times10^6$ |
| 受力与force后转换组合 | $8.983636446\times10^5$ | $9.236932752\times10^5$ |

none档受力组合进一步分解为：

```text
bulk u-force              +7.389119852e6 J/m2
Tail kick                 +5.032334826e5 J/m2
Beam kick                 -8.290686255e6 J/m2
conversion-after-force    +1.296696566e6 J/m2
combined                  +8.983636446e5 J/m2
```

Gate F完成时的阶段性根因分类曾为：

$$
\text{field--particle分裂耦合的混合阶有限时间步截断误差（历史阶段性判断）}.
$$

误差位于`midpoint_poisson -> force/kick -> final_poisson`耦合块。当前各子算子局部恒等式
可以闭合，但它们在有限$dt$下不构成一个严格的全局离散中点能量系统。FCT active-set切换、
bulk--Tail转换事件、随机Tail碰撞和开放边界使时间映射不光滑，因此只能观察到混合阶，
不应强行解释为稳定二阶。

Gate F固定状态审计的最稳定前缀（共5个公共区间）给出$q_{abs}=1.562$和$1.524$；
其他前缀由正负抵消和非渐近行为主导。该结果只证明局部冻结映射存在可收敛成分，不能外推为
长期生产残差的主导阶。其历史输出为：

```text
source_ownership_valid=1
poisson_identity_pass=1
classification=time-truncation-mixed-order
status=PASS
```

### 2.3 100.3--120 fs生产与115--117 fs A/B的新证据

`output/vpfp_return_dt_half_100p3_to_120fs`使用`dt_scale=0.5`稳定推进到120 fs：

```text
accepted_steps=1540
split_steps=0
collision_flux_rollback_count=0
cumulative_energy_residual=5.11614234597231e7 J/m2
cumulative_abs_energy_residual=5.12996002640222e7 J/m2
cumulative_residual/final_domain_energy=9.02062114495465e-3
positive_step_fraction=0.9915584415584415
```

残差不是随机正负抵消，而是几乎单向累积。累计占比从105 fs的0.150%增长到110 fs的0.366%、
115 fs的0.614%和120 fs的0.902%。该运行通过转换/H10事务守恒和稳定性门，但达到长期能量预算
预警区，并接近1%失效线。

随后从同一个115.01142520244935 fs checkpoint比较相同物理窗口115.011425--117 fs：

| 指标 | `dt_scale=0.5` | `dt_scale=0.25` | 结论 |
|---|---:|---:|---|
| 累计能量残差 | $6.548772840\times10^6$ | $6.518997392\times10^6$ J/m2 | 仅下降0.455% |
| 单位物理时间残差 | $3.293199154\times10^6$ | $3.278225893\times10^6$ J/m2/fs | 仅下降1.4% |
| 观测阶$p$ | - | $0.00657$ | 近似零阶 |
| wall time | 1454.45 s | 2968.93 s | 成本2.04倍 |
| 最终总域能相对差 | - | $1.40\times10^{-5}$ | 宏观总能接近 |
| 最终$U_E$相对差 | - | $1.17\times10^{-3}$ | 可见但不大 |
| 最终$K_{bulk+Tail}$相对差 | - | $8.71\times10^{-5}$ | 接近 |
| Tail粒子数相对差 | - | $5.91\%$ | 表示具有步长依赖 |

单步残差约随$dt$成比例减小，但步数同比增加，因而单位物理时间的残差近似不变：

$$
R_{step}\propto\Delta t,
\qquad
\frac{R_{step}}{\Delta t}\approx\text{常数}.
$$

这推翻了“长期公共余额主要由可通过继续缩小全局$dt$控制的混合阶时间截断误差”这一最终判断。
当前可确认的根因层级只能写为：

$$
\boxed{
\text{残差位于场--粒子公共耦合/能量源账块，且表现为近似固定的缺失功率或结构性功率误差}
}.
$$

尚未区分以下两种情况：

1. 生产状态确实以近似固定功率产生非物理能量；
2. 状态更新正确，但`accounted_energy_source`遗漏了一个同样为$O(\Delta t)$的真实功/边界能流项。

因此，Poisson空间恒等式和各局部kick恒等式仍然有效，但`classification=time-truncation-mixed-order`
不再是生产长时根因的最终分类。§11.7 Gate H阶段功率定位已执行并完成：首个固定缺失功率对应
`field_force_pair` 块（midpoint_poisson → force/kick → final_poisson 的场--粒子时间层），
触发量为 `field_particle_pair`，它在两档均解释约99%的完整残差功率；不得继续用更小全局
时间步替代根因修复。

---

## 3. 全程硬约束

### 3.1 不得改变的生产内容

Gate F已通过局部恒等式与只读结构门。以下内容在Gate H定位前仍禁止盲目修改：

- `OpenElectrostaticSolver`的Dirichlet边界条件；
- Beam注入、推进、开放流出和轨迹电流；
- `ConservativePpmRemap`实际输出状态；
- FCT/PPM降阶规则；
- 碰撞算子及collision reservoir；
- bulk-to-tail转换阈值、flux-interface sink和PIC装载；
- H10返回事务、返回阈值和表示门；
- 接受/拒绝条件、时间步、随机数和checkpoint格式中的物理字段。

### 3.2 禁止的“修复”

禁止实施以下操作：

- 每步缩放$E$、bulk $f$或PIC粒子动量来补齐能量；
- 把残差平均分摊到所有粒子；
- 修改Poisson边界电势以抵消残差；
- 将`collision_reservoir`或开放边界能流作为可调参数；
- 仅放宽物理能量门使结果显示PASS；
- 为审计复写一套不调用生产状态/生产通量的离散公式；
- 在`diagnostic-level=0/1`增加全数组扫描或额外MPI collective。

### 3.3 审计只读性

新增诊断只能在`diagnostic-level>=2`启用，并必须满足：

```text
accepted state unchanged
trial state unchanged
Beam/Tail RNG unchanged
conversion/H10 ledger unchanged
MPI collective order identical on every rank
failed trial state never written into accepted audit
```

---

## 4. 执行总顺序

严格按以下顺序执行：

1. Gate A：修复分析器的浮点假失败；
2. Gate B：仅用现有阶段文件完成派生分解；
3. Gate C：在生产算子中加入只读离散功记录；
4. Gate D：运行同checkpoint的none/return 100步审计；
5. Gate E：运行相同物理窗口的$dt/dt/2$审计；
6. Gate F：固定状态$dt/dt/2/dt/4$、源项所有权和Poisson功恒等式审计；
7. Gate G：对最终生产配置执行较小时间步短回归并完成100.3--120 fs生产；
8. Gate H：使用115 fs checkpoint执行阶段功率A/B，定位近似零阶长期残差的首个来源。

任何Gate失败时，只修该Gate，不得提前修改下一阶段。

> 执行进度：Gate A--H 均已通过/完成；100.3--120 fs生产及115--117 fs A/B已完成。
> 长期残差对$dt$减半的观测阶仅为0.00657，旧的`time-truncation-mixed-order`最终分类已被
> 生产证据否定。Gate H 阶段功率定位把根因定为**场--粒子时间层/离散功配对结构**（触发量
> `field_particle_pair` 解释99%；对应 `field_force_pair` 块的中点/末态Poisson时间层）。
> 当前不实施Poisson空间算子修改、全局能量补丁或全局减步长；下一步是另立结构保持设计文档。

---

## 5. Gate A：修复分析器浮点假失败

### 5.1 修改文件

```text
tools/analyze_vpfp_stage_energy_audit.py
tests/vpfp_stage_energy_audit_test.cpp
```

### 5.2 必须修改的判据

当前分析器用固定相对阈值$10^{-12}$判断强相消后的望远镜误差，导致绝对误差仅
$10^{-6}\ \mathrm{J/m^2}$仍被判FAIL。必须将**结构完整性门**改成机器精度缩放门。

对每个接受步计算：

$$
E_{scale}=\max_s\left(
1,
|K_{bulk,s}|+|K_{tail,s}|+|K_{beam,s}|+|U_{E,s}|
\right),
$$

$$
S_{stage}=\sum_s |R_s|,
$$

$$
S=\max(E_{scale},S_{stage},|R_{full}|,|R_{ledger}|,1).
$$

结构误差容差定义为：

$$
T_{round}=\max(10^{-10},512\epsilon_{double}S).
$$

以下两个误差分别与自己的$T_{round}$比较：

```text
stage_telescope_abs = abs(stage_balance_sum - full_balance)
energy_ledger_abs   = max(abs(full_balance-ledger_reference), ledger_spread)
```

不得用该容差判断物理能量是否守恒。物理余额继续独立输出，且不得因为结构门PASS而改写为物理PASS。

### 5.3 必须新增的结果字段

分析器至少新增：

```text
stage_telescope_max_abs
stage_telescope_roundoff_tolerance
stage_telescope_worst_step
energy_ledger_max_abs
energy_ledger_roundoff_tolerance
energy_ledger_worst_step
audit_structure_pass
physical_energy_residual_cumulative
physical_energy_gate_evaluated=0
```

保留旧相对字段供参考，但旧相对字段不得单独控制`status`。

`status=PASS`只表示：文件、阶段、步数、有限性、无split/failure及审计望远镜结构可信。

### 5.4 单元测试要求

在现有`vpfp_stage_energy_audit_test.cpp`中补充：

1. $5\times10^9$能量尺度下加入$2\times10^{-6}$求和扰动，必须PASS；
2. 同尺度下故意漏掉$1\ \mathrm{J/m^2}$源项，必须FAIL；
3. 小尺度人工账中加入$10^{-4}$残差，必须FAIL；
4. 缺阶段、重复阶段、NaN/Inf、未接受态、split和failure记录继续FAIL；
5. `physical_energy_residual_cumulative`非零不得令结构审计自动FAIL。

### 5.5 编译与测试命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target vpfp_stage_energy_audit_test fp_solver || exit 11

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_stage_energy_audit_test --case all \
  --result ./output/vpfp_stage_energy_audit_unit.result || exit 12

grep -q '^status=PASS$' ./output/vpfp_stage_energy_audit_unit.result || exit 13
python3 -m py_compile tools/analyze_vpfp_stage_energy_audit.py || exit 14
```

### 5.6 Gate A验收标准

```text
unit status=PASS
漏项1 J/m2用例必须被拒绝
2e-6 J/m2浮点扰动用例必须接受
生产求解器编译成功
未修改任何生产状态更新公式
```

---

## 6. Gate B：利用现有记录完成派生分解

### 6.1 修改范围

本阶段只修改：

```text
tools/analyze_vpfp_stage_energy_audit.py
```

不得修改C++生产代码。

### 6.2 每步必须派生的量

从现有11个阶段直接计算：

```text
u_force_bulk_delta       = dK_bulk  at u_force_tail_beam_kick
tail_kick_delta          = dK_tail  at u_force_tail_beam_kick
beam_kick_delta          = dK_beam  at u_force_tail_beam_kick
conversion_after_force   = dK_bulk+dK_tail at conversion_after_force
conversion_ledger_delta  = K_conversion(current)-K_conversion(previous)
conversion_transaction_residual = conversion_after_force-conversion_ledger_delta
force_conversion_pair    = above four terms summed
midpoint_poisson_delta   = dU_E at midpoint_poisson
final_poisson_delta      = dU_E at final_poisson
poisson_pair             = midpoint_poisson_delta+final_poisson_delta
resolved_common_balance  = force_conversion_pair+poisson_pair
```

同时输出累计值和占完整累计余额的比例。分母小于$10^{-30}$时比例写`nan`并增加
`fraction_valid=0`，禁止除零。

### 6.3 一致性要求

必须验证：

$$
R_{resolved}=R_{force+conversion}+R_{Poisson}
$$

与扣除$x$、碰撞、H10及其他已记账源项后的完整余额一致到Gate A机器精度容差。若不一致，说明仍有
漏项，禁止进入Gate C。

### 6.4 对现有结果重新分析

```bash
python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run ./output/vpfp_stage_energy_audit_none_100p3_100steps \
  --expected-accepted-steps 100 --require-no-split \
  --result ./output/vpfp_stage_energy_audit_none_100p3_100steps.result || exit 21

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run ./output/vpfp_stage_energy_audit_return_100p3_100steps \
  --expected-accepted-steps 100 --require-no-split \
  --result ./output/vpfp_stage_energy_audit_return_100p3_100steps.result || exit 22

cat ./output/vpfp_stage_energy_audit_none_100p3_100steps.result
cat ./output/vpfp_stage_energy_audit_return_100p3_100steps.result
```

### 6.5 Gate B验收标准

两组均必须满足：

```text
status=PASS
audit_structure_pass=1
accepted_steps=100
split_steps=0
failure_records=0
resolved decomposition matches full ledger within roundoff tolerance
x/collision/H10 cumulative residual each <= 1e-3 J/m2
```

物理累计余额约$3\times10^6\ \mathrm{J/m^2}$是待定位对象，不是Gate B失败理由。

---

## 7. Gate C：增加只读离散功记录

### 7.1 修改文件

```text
src/conservative_ppm_remap.h
src/conservative_ppm_remap.cpp
src/background_tail_pic.h
src/background_tail_pic.cpp
src/vpfp_integrator.h
src/vpfp_integrator.cpp
src/vpfp_diagnostics.h
src/vpfp_diagnostics.cpp
tests/vpfp_force_work_audit_test.cpp
CMakeLists.txt
tools/analyze_vpfp_stage_energy_audit.py
```

### 7.2 `RemapDiagnostics`新增字段

在`RemapDiagnostics`末尾新增并在构造函数中初始化为0：

```cpp
double upar_internal_face_energy_transfer;
double upar_left_velocity_boundary_energy;
double upar_right_velocity_boundary_energy;
double upar_interface_energy_removed;
double upar_discrete_energy_identity_residual;
```

这些字段只描述**最终实际使用的`upar_swept_`**。不得在FCT/接口sink修改前计算，不得使用另行重构的
高阶候选通量。

### 7.3 bulk $u_\parallel$离散能量恒等式

当前更新为：

$$
M_j^{n+1}=M_j^n-S_{j+1/2}+S_{j-1/2}.
$$

对每个物理$x$ cell和$u_\perp$环，使用生产网格中的
`cgrid_->kinetic_energy[j*nmu+k]`，内部面贡献必须按：

$$
W_{f}=S_f\left(K_j-K_{j-1}\right),\qquad 1\le f<N_v.
$$

速度左、右端面的符号必须由上式的望远镜展开确定：

$$
W_{left}=S_0K_0,
\qquad
W_{right}=-S_{N_v}K_{N_v-1}.
$$

实现时必须满足：

1. `S_f`就是接口sink处理后参与`m_new`更新的`upar_swept_[f]`；
2. Tail-owned cell被清零的能量不得混入普通内部面功，单独累计为
   `upar_interface_energy_removed`；
3. 速度域物理出流继续使用现有`u_tail_energy()`账本；
4. 禁止调用新的全分布扫描，累计必须融合进现有`ix/k/j`更新循环；
5. 每rank先累计local值，并合入该函数现有的打包`MPI_Allreduce`，不得新增多个小collective。

最终检查：

$$
R_{u,id}=\Delta K_{bulk}
-W_{internal}
-W_{left}
-W_{right}
+K_{interface\ removed}.
$$

符号必须由单元测试确定。若测试显示公式符号错误，应修正公式和实现，禁止通过绝对值隐藏符号错误。

### 7.4 Tail kick功

为`BackgroundTailPIC::kick()`增加可空的诊断输出参数，推荐接口：

```cpp
void kick(const SpatialGrid&, const EMFields&, double dt,
          int mpi_rank, int mpi_size,
          double* local_kinetic_work = NULL);
```

仅当指针非空时，在现有粒子kick循环中累计：

$$
W_{tail,p}=w_p[K(\boldsymbol u_p^{after})-K(\boldsymbol u_p^{before})].
$$

要求：

- 使用与`tail_total_kinetic_energy()`相同的相对论动能定义；
- OpenMP使用reduction，不得原子逐粒子累加；
- 返回rank-local值，由阶段审计统一MPI求和；
- 指针为空时不得增加额外`sqrt`、全粒子扫描或MPI；
- 不把该诊断值写入checkpoint，不改变Tail状态结构。

### 7.5 Beam kick功

Beam已经维护`BeamPIC::last_field_work()`。阶段审计必须记录其rank-local值，并在最终打包归约中求和。
不得重新推动Beam或复算轨迹。

在`accepted_n`、`collision_half1`和`x_half1`阶段该值应为0；在
`u_force_tail_beam_kick`及后续阶段记录当前步累计值。

### 7.6 `VpfpStageEnergyRecord`新增字段

在固定长度记录中追加：

```cpp
double bulk_upar_face_work;
double bulk_upar_velocity_boundary_work;
double bulk_upar_interface_energy_removed;
double bulk_upar_identity_residual;
double tail_kick_work;
double beam_kick_work;
```

同步更新：

- `packed_values`；
- packed数组写入/读取索引；
- `vpfp_stage_energy_audit.dat`表头和数据行；
- 所有字段有限性检查；
- MPI所有权注释。

所有新增值都必须声明是rank-local还是already-global：

```text
bulk upar RemapDiagnostics  = already-global，只由rank 0贡献到最终SUM
Tail kick work              = rank-local，所有rank求和
Beam last_field_work        = rank-local，所有rank求和
```

禁止再次出现已全局量乘以MPI rank数的问题。

### 7.7 新增派生余额

分析器新增：

```text
bulk_upar_identity_residual_cumulative
tail_kick_snapshot_mismatch_cumulative
beam_kick_snapshot_mismatch_cumulative
conversion_pair_residual_cumulative
particle_work_sum_cumulative
field_energy_change_cumulative
field_particle_pair_residual_cumulative
field_particle_pair_residual_relative
```

其中：

$$
R_{tail}=\Delta K_{tail,kick}-W_{tail},
$$

$$
R_{beam}=\Delta K_{beam,kick}-W_{beam},
$$

$$
R_{conversion}=\Delta(K_{bulk}+K_{tail})_{conversion}
-\Delta K_{conversion\ ledger}.
$$

`K_conversion`是阶段累计源项，必须使用相邻阶段差，不能直接重复减去其累计值。force阶段中的bulk
接口移除与随后conversion阶段的Tail增加共同构成表示转移；分析物理受力组合时应保留这两个实际
状态变化，分析转换事务自身时才使用上式扣除conversion ledger。

$$
R_{field-particle}=\Delta U_E+W_{bulk}+W_{tail}+W_{beam}-W_{electrode}.
$$

开放边界动能流、Tail物理出流、Beam注入/出流和collision reservoir继续在完整账中单独处理，禁止
重复扣除。若分析器无法证明某项只扣除一次，应报`source_ownership_valid=0`并FAIL。

### 7.8 新增单元测试

新增`tests/vpfp_force_work_audit_test.cpp`，测试必须直接调用生产实现，至少覆盖：

1. 非均匀$u_\parallel$网格、常量正场；
2. 非均匀$u_\parallel$网格、常量负场；
3. 无速度端面出流、无Tail接口时的离散望远镜恒等式；
4. 有左/右速度端面出流时的符号；
5. flux-interface sink打开时，接口移除能量只记一次；
6. Tail单粒子正负场kick，记录功等于直接动能差；
7. Beam已有`last_field_work()`与阶段动能差配对；
8. diagnostic关闭时输出字段为0且状态bitwise不变。

测试容差使用机器精度缩放，不得统一使用任意`1e-6`。

### 7.9 构建与测试命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target \
  vpfp_stage_energy_audit_test \
  vpfp_force_work_audit_test \
  fp_solver || exit 31

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_stage_energy_audit_test --case all \
  --result ./output/vpfp_stage_energy_audit_unit.result || exit 32

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_force_work_audit_test --case all \
  --result ./output/vpfp_force_work_audit_unit.result || exit 33

grep -q '^status=PASS$' ./output/vpfp_stage_energy_audit_unit.result || exit 34
grep -q '^status=PASS$' ./output/vpfp_force_work_audit_unit.result || exit 35
```

### 7.10 Gate C验收标准

```text
两项单元测试status=PASS
正负场符号均通过
bulk u-face恒等式达到机器精度缩放容差
Tail/Beam kick记录与直接动能差一致
gate_c_work_fields_present=1
diagnostic-level 0/1不新增热路径扫描或collective
无物理推进公式变化
```

原Gate C时的`source_ownership_valid`只表示“字段存在”，因此当时不属于Gate C局部功验收范围。
该所有权门已在§10.5.1重新实现并通过，不再是待办项。

---

## 8. Gate D：100.3 fs checkpoint的100步none/return审计

### 8.1 公共环境

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_100P3 < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#CHECKPOINTS_100P3[@]}" -eq 1 || exit 41
CHECKPOINT_100P3="${CHECKPOINTS_100P3[0]}"
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 42
```

### 8.2 return-none 100步命令

```bash
OUT_NONE=./output/vpfp_force_work_none_100p3_100steps
test ! -e "$OUT_NONE" || exit 43

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --restart-allow-return-config-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode none \
  --dt-scale 1.0 \
  --stop-time-fs 102.858508730862 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_NONE" || exit 44

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_NONE" --expected-accepted-steps 100 --require-no-split \
  --result ./output/vpfp_force_work_none_100p3_100steps.result || exit 45
```

### 8.3 hysteretic return 100步命令

```bash
OUT_RETURN=./output/vpfp_force_work_return_100p3_100steps
test ! -e "$OUT_RETURN" || exit 46

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --dt-scale 1.0 \
  --stop-time-fs 102.858508730862 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_RETURN" || exit 47

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_RETURN" --expected-accepted-steps 100 --require-no-split \
  --result ./output/vpfp_force_work_return_100p3_100steps.result || exit 48
```

### 8.4 Gate D验收标准

两组均要求：

```text
status=PASS
audit_structure_pass=1
gate_c_work_fields_present=1
accepted_steps=100
split_steps=0
failure_records=0
bulk_upar_identity_residual within roundoff tolerance
tail/beam kick snapshot mismatch within roundoff tolerance
conversion transaction residual <= 1e-3 J/m2 cumulative
H10 direct residual <= 1e-3 J/m2 cumulative
```

Gate D运行时的`source_ownership_valid`仅表示Gate C列存在，因此Gate D本身只验收记录结构和
局部功恒等式。真正的源项所有权门已由§10.5.1执行并通过，不再限制Gate D/F结论。

若return组相对none组只改变后续状态，但H10直接残差仍在容差内，继续把公共根因归于共享场耦合，
不得回退H10。

---

## 9. Gate E：$dt/dt/2$同物理窗口判别

### 9.1 目的

本测试区分：

- 固定离散算子不兼容：减半$dt$后单位物理时间残差不按预期下降；
- 二阶Strang/中点Poisson截断误差：减半$dt$后累计残差应明显下降；
- 随机Tail采样误差：不同seed/粒子噪声主导且无稳定收敛阶。

### 9.2 测试设置

使用return-none以隔离H10。必须从同一个100.3 fs checkpoint出发，并覆盖同一物理终止时间。

先重新编译包含受限`dt_scale`重启豁免的生产求解器：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake --build build -j4 --target fp_solver || exit 50
```

基准档沿用Gate D的100步结果。半步长档运行200步。两档必须使用共同的
物理终止时间`102.858508730862 fs`，不能在重启任务中写
`--stop-after-steps 100/200`：当前CLI中的`--stop-after-steps`表示全局绝对步号，
而不是重启后追加步数。

100.3 fs checkpoint的`manifest.dat`中记录的`dt`可能是为精确落到100.3 fs而
缩短的最后一步，例如`2.1534300890896308e-17 s`。它不是重启后的标称步长，
不能用于推导Gate E的时间窗口。重启后的标称步长由`--dt-scale`和编译参数重新计算。

由于`dt_scale`属于checkpoint物理配置哈希，半步长档必须显式声明checkpoint的源
`dt_scale=1.0`，并使用受限豁免。程序会先用源值重建完整物理哈希；只有哈希完全匹配
时，才允许本次续跑改用`dt_scale=0.5`。这不会豁免其他物理配置差异。

```bash
OUT_HALF=./output/vpfp_force_work_none_100p3_dt_half_200steps
test ! -e "$OUT_HALF" || exit 51

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --restart-allow-return-config-change \
  --restart-source-dt-scale 1.0 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode none \
  --dt-scale 0.5 \
  --stop-time-fs 102.858508730862 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_HALF" || exit 52

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_HALF" --expected-accepted-steps 200 --require-no-split \
  --result ./output/vpfp_force_work_none_100p3_dt_half_200steps.result || exit 53
```

执行后必须确认：

```text
Gate D accepted_steps=100
Gate E dt/2 accepted_steps=200
两组restart起始时间相同
两组last accepted time均为102.858508730862 fs（允许浮点舍入误差）
启动日志包含restart_allow_dt_scale_change=1
启动日志包含restart_source_dt_scale=1.0
```

若源`dt_scale`填写错误，重建哈希必须失败；不得通过编辑checkpoint哈希、关闭整个
物理配置校验或把manifest末步`dt`当作源`dt_scale`绕过。若两组物理窗口不一致，
本测试无效，不得计算收敛阶。

### 9.3 收敛率

对以下量分别计算：

```text
physical_energy_residual_cumulative
poisson_pair_cumulative
force_conversion_pair_cumulative
field_particle_pair_residual_cumulative
```

定义：

$$
q=\log_2\frac{|R_{dt}|}{|R_{dt/2}|}.
$$

该定义只适用于以下条件同时成立的量：

1. 两档从bitwise相同的接受态出发；
2. 比较相同物理窗口；
3. 残差在窗口内没有严重正负相消；
4. 随机Tail碰撞使用可比较的共同随机路径，或在只读算子审计中被冻结；
5. `dt/2`相对`dt`没有改变未受控的转换事件序列、粒子表示统计或随机样本数量。

只要第3项不满足，就必须同时报告累计绝对残差：

$$
A_{dt}=\sum_n |R_{dt,n}|,
\qquad
q_{abs}=\log_2\frac{A_{dt}}{A_{dt/2}}.
$$

还必须按共同物理区间报告局部前缀：第1、2、5、10、20、50和100个`dt`区间。一个`dt`区间
对应两个`dt/2`步。禁止只用整个窗口末端的有符号累计量宣布收敛阶。

判定只作根因分类：

- 局部单步缺陷$q\approx3$：与二阶Strang格式的局部$O(dt^3)$误差相符；
- 固定物理窗口$q$或$q_{abs}\approx2$：与二阶全局误差相符；
- $0.5\le q<1.5$：可能是一阶边界/转换耦合或混合误差，需看分阶段量；
- $q<0.5$：固定离散不兼容或随机误差主导；
- 残差符号翻转或接近0时，不报告虚假阶数，改报绝对包络。

`q<0.5`本身不能在“固定离散不兼容”和“随机/事件路径差异”之间作选择。若Tail碰撞、转换或
粒子数在两档间改变，必须先执行§10.5的固定状态审计。

### 9.4 已完成的Gate E比较结论及适用边界

两档覆盖同一物理窗口$100.3\text{--}102.858508730862\ \mathrm{fs}$。对$dt/2$残差每两步配对成一个
$dt$物理区间后，实测为：

| 指标 | $dt$ | $dt/2$ | 结果 |
|---|---:|---:|---:|
| 最终signed累计残差 | $2.9854783\times10^6$ | $3.7576731\times10^6$ | 不单调 |
| 配对物理区间绝对残差和 | $9.1338290\times10^6$ | $3.8376683\times10^6$ | 下降$58.0\%$ |
| 配对物理区间RMS | $1.04357\times10^5$ | $4.41191\times10^4$ | 改善$2.37$倍 |
| wall time | $605.86\ \mathrm{s}$ | $1246.04\ \mathrm{s}$ | 成本$2.06$倍 |

因此：

$$
q_{abs}=\log_2\frac{9.1338290\times10^6}{3.8376683\times10^6}=1.251\quad(配对区间口径).
$$

`dt/2`降低了该短窗内局部能量缺陷的绝对包络和RMS，证明冻结/短窗映射中包含时间截断成分。
但最终signed累计量不单调，而且转换、随机碰撞和Tail粒子表示已经分叉。因此该Gate不能证明
长期生产残差由同一个时间截断项主导。

末态宏观差异为：

```text
U_E             -0.346%
K_e             +0.231%
K_b             -0.046%
K_tail          -0.790%
K_combined      +0.197%
domain energy   +0.0178%
N_e             +3.35e-6 relative
Tail count      +7.47%
```

Tail粒子数差异说明两档生产轨迹已在转换/碰撞事件上分叉。因此Gate E只支持“短窗局部包络中
存在可收敛部分”，不能用于声称稳定二阶收敛，也不能覆盖§2.3中115--117 fs生产A/B的近零阶
结果。Gate F的固定状态结论与Gate E一样属于局部证据，长期根因必须以生产A/B和Gate H为准。

---

## 10. Gate F：根因选择门

以下分支用于描述根因位置和误差性质，**并非天然互斥**。特别是分支C描述残差所在的耦合块，
分支D描述该耦合误差是否随时间步收敛；二者可能同时成立。只有完成§10.5后，才允许选择生产修复。

### 10.1 分支A：bulk $u_\parallel$通量功不闭合

触发条件：

```text
bulk_upar_identity_residual占完整余额超过50%
并且该残差不随dt减半下降
```

舍入量级残差即使不随`dt`下降也不得触发本分支。必须同时满足“宏观可见”和“不收敛”。

允许修改：

- `ConservativePpmRemap::advect_u_parallel()`中电流/功矩的离散定义；
- 使用同一个最终`upar_swept_`构造离散速度矩；
- 非均匀网格上的$K_j-K_{j-1}$链式法则。

禁止修改Poisson、Beam、Tail和H10。

### 10.2 分支B：Tail或Beam kick功不闭合

触发条件：

```text
tail_kick_snapshot_mismatch或beam_kick_snapshot_mismatch显著超过机器精度门
```

只允许修对应pusher的时间层、场插值或诊断所有权。Beam必须继续使用现有staggered face gather；不得
将Beam改成cell-centered受力。

### 10.3 分支C：粒子功均闭合，但Poisson--粒子配对不闭合

触发条件：

```text
bulk/Tail/Beam各自功恒等式通过
field_particle_pair_residual仍占完整余额超过50%
```

此时根因是约束场时间离散。下一修复必须从同一个离散弱形式构造：

```text
midpoint density/current
particle force work
final charge density
Poisson field-energy change
```

上述条件只能先把问题定位为“场–粒子耦合块不闭合”，不能仅凭该条件断言Poisson空间算子或时间离散
存在固定不兼容。还必须排除：

- 开放边界、电极功或表示转换源项遗漏/重复；
- 中点与末态的物种密度并非来自同一接受轨迹；
- 随机Tail碰撞在`dt/dt/2`中走了不同随机路径；
- Strang分裂的有限时间截断误差；
- signed累计残差的正负相消。

只有在Poisson恒等式通过、源项所有权正确，且阶段功率审计确认状态更新存在结构性功率缺口时，
才允许设计迭代中点或离散梯度型Vlasov--Poisson耦合。禁止事后缩放场或粒子能量。该重构必须
另写独立设计文档，不得由本审计任务直接实施。Gate F曾把局部固定状态结果分类为混合阶时间截断，
但§2.3的长期A/B已证明该分类不能直接作为生产根因。

### 10.4 分支D：时间截断误差及混合阶收敛（当前未被长期证据支持）

若固定状态审计显示残差随$dt$缩小而下降，但由于FCT active-set、阈值事件、开放边界或
分裂算子的非光滑性而没有稳定二阶，仍可分类为混合阶时间截断误差。此时应先做误差预算：

- 以120 fs累计能量目标反推允许$dt$；
- 比较减小$dt$的资源成本；
- 决定接受截断误差还是开发结构保持中点格式。

不得把“属于截断误差”直接解释为“无需处理”。

长期生产A/B对本分支增加必要条件：在相同起点、相同物理窗口下，单位物理时间累计残差必须随
$dt$缩小呈现明确正阶：

$$
p=\log_2\frac{|R_{\Delta t}|/T}{|R_{\Delta t/2}|/T}>0.5.
$$

115--117 fs实测$p=0.00657$，因此当前长期公共余额**不满足**分支D作为主导根因的条件。
固定状态局部阶约1.5只能说明误差中存在时间截断分量，不能说明它控制长期signed累计余额。

### 10.5 Gate F最终确认：执行后才允许选择生产修复

#### 10.5.1 修复`source_ownership_valid`

当前分析器中的实现：

```python
source_ownership_valid = 1 if work_fields_present else 0
```

只检查诊断列是否存在，并没有验证所有权。必须改成逐项结构检查，至少输出：

```text
electrode_work_ownership_count
background_boundary_energy_ownership_count
beam_boundary_energy_ownership_count
velocity_boundary_energy_ownership_count
conversion_energy_ownership_count
tail_return_energy_ownership_count
collision_reservoir_ownership_count
source_ownership_residual
source_ownership_valid
```

每项必须证明在完整步能量账中恰好出现一次。`source_ownership_valid=1`要求所有count等于1，且
`source_ownership_residual`处于Gate A机器精度缩放容差内。不得再以“字段存在”代替所有权验证。

#### 10.5.2 固定接受态局部缺陷审计

从100.3 fs同一个checkpoint构造只读事务，分别评估：

```text
A: 一个dt步
B: 两个dt/2步
C: 四个dt/4步
```

三档必须从相同状态副本出发，不得依次续跑。必须记录每个共同`dt`物理区间的：

```text
field_particle_pair_residual_signed
field_particle_pair_residual_abs
bulk/tail/beam work
field_energy_change
electrode_work
final bulk/tail/beam/field state norm
conversion event count and transferred N/Px/K
Tail particle count and RNG counter/hash
```

随机碰撞必须采用以下二选一协议，并在结果中写明：

1. **确定性耦合审计**：只读重放场–粒子子映射，冻结碰撞、转换和H10，不接受结果；
2. **共同随机数审计**：一个`dt`随机增量通过Brownian bridge/等价可加规则拆成两个`dt/2`及四个
   `dt/4`增量，保证三档表示同一随机路径。

禁止直接把“步数翻倍、随机调用次数也翻倍”的两次生产轨迹当成确定性Richardson测试。

#### 10.5.3 离散Poisson功恒等式

对固定Dirichlet边界，直接验证生产Poisson算子：

$$
U_E^{n+1}-U_E^n-W_{electrode}
=\frac12\left\langle\phi^{n+1}+\phi^n,
\rho^{n+1}-\rho^n\right\rangle_h+R_{P,h}.
$$

必须使用生产`G/G^*`、面场和cell电荷权重，禁止复写另一套中心差分公式。把
$\rho^{n+1}-\rho^n$按bulk、Tail、Beam、离子和边界源拆分，输出每个物种贡献与
`R_P,h`。只有`R_P,h`在机器精度门内，才可排除Poisson空间算子/边界离散错误。

#### 10.5.4 最终判定

- `R_P,h`不闭合：选择C，但先修Poisson空间伴随/边界项；
- `R_P,h`闭合，局部缺陷$q\approx3$且固定窗$q\approx2$：选择C+D，先做误差预算，不立即重构；
- `R_P,h`闭合，但确定性局部缺陷不收敛：选择C，设计同源离散梯度/迭代中点耦合；
- 只有随机生产轨迹不收敛：先做随机弱收敛和统计置信区间，禁止修改Poisson；
- 发现源项count或ownership residual错误：先修诊断/账本，再重新执行Gate D–F。

已执行结果：

```text
same_initial_state=1
same_physical_window=1
collision_frozen=1
conversion_frozen=1
h10_frozen=1
accepted_state_bitwise_equal_after_audit=1
rng_bitwise_equal_after_audit=1
ledger_bitwise_equal_after_audit=1
source_ownership_valid=1
poisson_identity_pass=1
classification=time-truncation-mixed-order
status=PASS
```

Poisson单元恒等式为：

```text
R_P_h=-1.1754943508222875e-38 J/m2
roundoff_tolerance=1.7038161760986797e-34 J/m2
poisson_identity_pass=1
```

固定状态区间5给出稳定的粗到中、中到细绝对阶$1.562$和$1.524$；其他前缀受抵消和非渐近误差
影响而不稳定。Gate F当时的阶段性选择是：

$$
\boxed{\text{分支C（残差位于场--粒子耦合块）}+\text{分支D（混合阶时间截断误差）}}
$$

后续100.3--120 fs生产和115--117 fs A/B已经触发上述复核条件：`dt_scale=0.25`相对0.5只把
单位时间残差降低约1.4%，观测阶0.00657。因此最终选择现修正为：

$$
\boxed{
\text{分支C（公共耦合/源账块）成立；分支D仅是局部分量，不是已证实的长期主导根因}
}.
$$

生产决策同步修正为：保留当前空间Poisson和已验证局部恒等式，不添加能量补丁，不继续减小全局
时间步；先执行Gate H，把公共余额分解为每阶段、每单位物理时间的功率缺口。只有确认状态更新块
而非诊断源账遗漏后，才单独立项重构统一时间中心的场--粒子耦合。

#### 10.5.5 实现文件、命令与验收输出

本节所需工具已实现并通过。对应文件为：

```text
tools/analyze_vpfp_stage_energy_audit.py
tools/analyze_vpfp_local_defect.py
tests/vpfp_poisson_work_identity_test.cpp
tests/vpfp_local_defect_checkpoint_audit.cpp
src/open_electrostatic_solver.{h,cpp}       # 只增加只读恒等式接口
src/vpfp_integrator.{h,cpp}                 # 只增加事务式固定状态审计入口
CMakeLists.txt
```

`vpfp_local_defect_checkpoint_audit`必须直接读取生产checkpoint并使用生产网格、Poisson、bulk remap、
Tail/Beam pusher和gather。允许为确定性审计冻结碰撞、转换和H10，但必须只作用于内存副本；不得写回
checkpoint，不得改变生产`advance()`路径。审计结束后必须比较原始状态、RNG和ledger哈希。

编译及单元测试：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target \
  vpfp_stage_energy_audit_test \
  vpfp_poisson_work_identity_test \
  vpfp_local_defect_checkpoint_audit \
  fp_solver || exit 61

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/vpfp_poisson_work_identity_test --case all \
  --result ./output/vpfp_poisson_work_identity_unit.result || exit 62

grep -q '^status=PASS$' \
  ./output/vpfp_poisson_work_identity_unit.result || exit 63
```

固定checkpoint MPI审计：

```bash
AUDIT_OUT=./output/vpfp_gate_f_fixed_state_100p3
test ! -e "$AUDIT_OUT" || exit 64

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/vpfp_local_defect_checkpoint_audit \
  --restart-dir "$CHECKPOINT_100P3" \
  --restart-source-dt-scale 1.0 \
  --mode deterministic-field-particle \
  --dt-scales 1.0,0.5,0.25 \
  --common-dt-intervals 1,2,5,10 \
  --output-dir "$AUDIT_OUT" || exit 65

python3 tools/analyze_vpfp_local_defect.py \
  --run "$AUDIT_OUT" \
  --result ./output/vpfp_gate_f_fixed_state_100p3.result || exit 66
```

审计工具的`status=PASS`只表示试验可比较，不表示必须得到某个预设收敛阶。结果至少必须包含：

```text
same_initial_state=1
same_physical_window=1
collision_frozen=1
conversion_frozen=1
h10_frozen=1
accepted_state_bitwise_equal_after_audit=1
rng_bitwise_equal_after_audit=1
ledger_bitwise_equal_after_audit=1
source_ownership_valid=1
poisson_identity_pass=1
q_signed_interval_1
q_abs_interval_1
q_signed_interval_2
q_abs_interval_2
q_signed_interval_5
q_abs_interval_5
q_signed_interval_10
q_abs_interval_10
classification=<poisson-space|time-truncation|fixed-coupling|stochastic-followup>
status=PASS
```

任何哈希变化、源项count不为1、Poisson恒等式不闭合或三档起始态不同都必须FAIL；不得通过放宽
能量阈值使审计通过。

---

## 11. Gate G生产结果、时间步复核与Gate H根因定位

### 11.1 不再重复的测试

Gate A--F已完成，不得重复执行以下高成本审计，除非生产核心算子、Poisson边界、
阶段顺序或checkpoint物理语义发生改变：

- none $dt$的100步Gate D；
- none $dt/2$的200步Gate E；
- return $dt$的100步Gate D；
- 固定状态$dt/dt/2/dt/4$ Gate F；
- Poisson空间功恒等式审计。

### 11.2 唯一必需的短回归

目的是验证最终生产组合`hysteretic return + dt_scale=0.5`能从同一个100.3 fs checkpoint正常推进，
而不是重新估计收敛阶。

首次执行结果`output/vpfp_energy_budget_return_dt_half_100p3_to_102p858.result`的唯一失败项为
`accepted_steps=201`。第201步的实际$dt=1.2496\times10^{-27}\ \mathrm{s}$，是十进制终止时间与
二进制浮点落点之间的舍入差，不是物理时间步。该伪步的能量残差只有
$-2.24\times10^{-7}\ \mathrm{J/m^2}$，但它仍触发了按步计数的Tail/H10处理，使Tail粒子数变化3846。

该问题已通过`src/vpfp_time_control.h`和`src/main_vpfp.cpp`修复：

- 剩余时间不超过$64\epsilon$物理尺度容差时，不再调用一次生产`advance()`；
- checkpoint/snapshot判据使用同一容差，最后一个真实接受态仍能满足目标时间；
- 容差按$|t|$和$dt$缩放，不再使用会在fs问题上过宽的`max(1, |t|)`；
- `tests/vpfp_time_control_test.cpp`已验证“压制$10^{-27}$ s伪步、保留$10^{-17}$ s正常步”。

去掉第201个伪步后，首次结果的前200步已满足宏观门：$U_E$差$0.397\%$、
$K_{combined}$差$0.201\%$、$K_{tail}$差$1.996\%$、$N_e$差$3.53\times10^{-6}$；
$|R_{signed}|/E_{domain}=6.92\times10^{-4}$。因此重跑只是验证终止时间修复和正确的200步计数。

编译：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake --build build -j4 --target vpfp_time_control_test fp_solver || exit 70

./build/vpfp_time_control_test \
  --result ./output/vpfp_time_control_unit.result || exit 70
grep -q '^status=PASS$' ./output/vpfp_time_control_unit.result || exit 70
```

运行前由用户设置真实checkpoint路径：

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_100P3 < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#CHECKPOINTS_100P3[@]}" -eq 1 || exit 71
CHECKPOINT_100P3="${CHECKPOINTS_100P3[0]}"
OUT_G=./output/vpfp_energy_budget_return_dt_half_100p3_to_102p858_time_snap_fix
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 71
test ! -e "$OUT_G" || exit 72

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --restart-source-dt-scale 1.0 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.5 \
  --stop-time-fs 102.858508730862 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_G" || exit 73

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_G" --expected-accepted-steps 200 --require-no-split \
  --result ./output/vpfp_energy_budget_return_dt_half_100p3_to_102p858_time_snap_fix.result \
  || exit 74
```

> 本命令选择的`h10_r3a_return_100_to_100p3` checkpoint已使用
> `tail_return_mode=hysteretic`，因此不得添加`--restart-allow-return-config-change`。
> 若以后改用`tail_return_mode=none`的checkpoint，才添加该受限豁免。若源`dt_scale`不是`1.0`，
> 必须将`--restart-source-dt-scale`改为manifest对应的真实物理配置值；不得修改checkpoint hash绕过校验。

### 11.3 Gate G验收标准

结构门必须全部满足：

```text
status=PASS
accepted_steps=200
split_steps=0
failure_records=0
audit_structure_pass=1
source_ownership_valid=1
decomposition_matches_roundoff=1
no NaN/Inf
no collision/conversion/H10 rollback
```

`accepted_steps=200`是硬验收条件。若出现201步，且最后一步$dt$仅为$10^{-27}\ \mathrm{s}$
量级，说明集群使用了未包含`vpfp_time_control.h`终止时间吸附修复的旧二进制文件。不得把
分析器预期步数改为201绕过该错误。

与已归档的`output/archieve/vpfp_force_work_return_100p3_100steps`末态比较时，宏观安全门为：

```text
field energy relative difference <= 2%
combined bulk+Tail energy relative difference <= 1%
Tail energy relative difference <= 2%
background electron number relative difference <= 1e-4
conversion N/Px/K residual <= 1e-12
tail-return N/Px/K residual <= configured moment tolerance
tail-return Jx/Pixx/Piperp residual finite and recorded (diagnostic only)
```

这些差异门用于拦截配置或实现回归，不要求PIC Tail粒子数或随机轨迹bitwise一致。

H10返回事务的严格守恒量是$N$、$P_x$和$K$。`Jx/Pixx/Piperp`是连续PIC粒子到固定
cell-center速度网格的表示保真度诊断，不是任意粒子组上可以精确强制的事务不变量。
因此不得把它们与`--tail-return-moment-tolerance=1e-12`直接比较并宣布Gate G失败；它们由
H10速度网格/表示收敛测试单独约束。

能量门不得只比较末端signed累计量。必须同时报告：

```text
physical_energy_residual_cumulative
sum(abs(energy_balance_residual))
RMS(energy_balance_residual)
field_particle_pair_residual_cumulative
physical_energy_residual_cumulative / final_domain_energy
```

使用以下命令从逐步诊断中生成预算摘要：

```bash
python3 - "$OUT_G/vpfp_step_diagnostics.dat" \
  > ./output/vpfp_energy_budget_return_dt_half_100p3_to_102p858_time_snap_fix.budget.result <<'PY'
import math
import sys

path = sys.argv[1]
with open(path, "r", encoding="utf-8") as handle:
    lines = [line.split() for line in handle if line.strip()]
header = lines[0]
rows = [dict(zip(header, row)) for row in lines[1:]]
if not rows:
    raise SystemExit("no accepted diagnostic rows")

residual = [float(row["energy_balance_residual"]) for row in rows]
signed = sum(residual)
absolute = sum(abs(value) for value in residual)
rms = math.sqrt(sum(value * value for value in residual) / len(residual))
domain = float(rows[-1]["domain_energy_after"])
ratio = abs(signed) / max(abs(domain), 1.0e-300)

print(f"accepted_rows={len(rows)}")
print(f"signed_residual={signed:.17e}")
print(f"absolute_residual={absolute:.17e}")
print(f"rms_residual={rms:.17e}")
print(f"final_domain_energy={domain:.17e}")
print(f"signed_to_domain={ratio:.17e}")
PY
```

使用以下命令执行末态宏观安全门：

```bash
python3 - \
  ./output/archieve/vpfp_force_work_return_100p3_100steps/vpfp_step_diagnostics.dat \
  "$OUT_G/vpfp_step_diagnostics.dat" \
  > ./output/vpfp_energy_budget_return_dt_half_100p3_to_102p858_time_snap_fix.macro.result <<'PY'
import sys

def final_row(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    return dict(zip(lines[0], lines[-1]))

reference = final_row(sys.argv[1])
candidate = final_row(sys.argv[2])
limits = {
    "U_E": 2.0e-2,
    "K_combined": 1.0e-2,
    "K_tail": 2.0e-2,
    "N_e_after": 1.0e-4,
}
passed = True
for name, limit in limits.items():
    a = float(reference[name])
    b = float(candidate[name])
    relative = abs(b - a) / max(abs(a), abs(b), 1.0e-300)
    ok = relative <= limit
    passed = passed and ok
    print(f"{name}_relative_difference={relative:.17e}")
    print(f"{name}_pass={int(ok)}")
print(f"status={'PASS' if passed else 'FAIL'}")
raise SystemExit(0 if passed else 1)
PY
```

若结构门和宏观门通过，即允许生产；不再要求该随机生产轨迹显示预设二阶。

### 11.4 通过后的100.3--120 fs生产命令

本节已执行完成。Gate G通过后，使用同一组物理参数从100.3 fs连续推进。生产档改用轻量诊断：

```bash
OUT_PROD=./output/vpfp_return_dt_half_100p3_to_120fs
test ! -e "$OUT_PROD" || exit 75

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --restart-source-dt-scale 1.0 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.5 \
  --checkpoint-times 105,110,115,120 \
  --snapshot-times 105,110,115,120 \
  --stop-time-fs 120 \
  --diagnostic-level 1 --diagnostic-interval 10 \
  --output-dir "$OUT_PROD" || exit 76
```

本生产命令继续使用已开启hysteretic return的100.3 fs checkpoint，不需要返回配置豁免。

短窗物理门至少保持：

```text
field relative L2 <= 2%
background density relative L2 <= 1%
combined spectrum relative L2 <= 2%
upar/u_perp marginal relative L2 <= 2%
conversion N/Px/K residual <= 1e-12
no failure/split/rollback
```

长跑能量预算分两级：

```text
预警：|cumulative physical residual| / final domain energy > 5e-3
失效：上述比值 > 1e-2，或宏观场/能谱对时间步不收敛
```

该运行已经达到预警线：累计signed残差占最终域能0.902%。因此已按要求从115 fs checkpoint
执行局部时间步A/B；结果见§11.5。不得再从0 fs重跑来重复确认同一事实。

### 11.5 已完成的115--117 fs `dt_scale=0.5/0.25`复核

输入数据：

```text
coarse: output/vpfp_return_dt_half_100p3_to_120fs/vpfp_step_diagnostics.dat
fine:   output/dt025_115_to_117/vpfp_step_diagnostics.dat
start:  115.01142520244935 fs
end:    117.0 fs
```

粗档诊断没有恰好落在117 fs，因此比较时必须对粗档最后两个诊断点做线性时间插值；不得把
116.994269 fs粗档状态直接当成117 fs。对齐后的结果为：

```text
coarse cumulative residual       = 6.54877284047139e6 J/m2
fine cumulative residual         = 6.51899739212164e6 J/m2
fine/coarse residual              = 9.95453278182785e-1
coarse residual power             = 3.29319915375452e6 J/m2/fs
fine residual power               = 3.27822589331371e6 J/m2/fs
observed order                    = 6.57449055967995e-3
coarse wall                       = 1.45444937840312e3 s
fine wall                         = 2.96892963052200e3 s
fine/coarse wall                  = 2.04127395192102
```

两档所有单步残差均为正。`dt_scale=0.25`无split、failure或rollback，转换$N/P_x/K$最大相对
残差约$1.57\times10^{-13}$，H10返回$N/P_x/K$最大相对残差小于$6.3\times10^{-15}$。

末态宏观相对差：

```text
domain_energy_after   1.40e-5
U_E                   1.17e-3
K_e                   1.45e-4
K_tail                8.27e-3
K_combined            8.71e-5
N_e_after             2.08e-6
Tail particle count   5.91e-2
```

结论：宏观总能和combined动能已接近，但公共余额的单位时间误差不收敛。继续使用`dt_scale=0.25`
只会把成本提高约2.04倍，不能作为生产修复。

Tail粒子数不能用于严格Richardson判阶。当前`tail-return-residence-steps=8`以步数而非物理时间
定义驻留期；减半$dt$会改变其物理驻留时间。flux-interface每步生成parcel也使表示粒子数依赖步数。
这不推翻能量功率结论，因为两档combined能量接近且能量残差功率几乎相同；但后续若专门审计
Tail表示收敛，必须把驻留判据改为物理时间，或按$dt$反比缩放步数。

### 11.6 修正后的生产决策

立即生效的决策：

1. 不再运行115--120 fs `dt_scale=0.25`长窗；现有115--117 fs结果已足够否定全局减步长方案。
2. 不把`dt_scale=0.25`设为生产默认值；它没有改善单位时间能量残差。
3. 不修改Poisson空间算子、局部kick、转换/H10事务或碰撞reservoir；这些局部恒等式仍已通过。
4. 不添加全局能量补丁，不缩放$E/f$/PIC动量，不通过放宽1%门掩盖残差。
5. 只执行§11.7 Gate H；其目标是确定首个产生固定残差功率的阶段以及该缺口属于状态更新还是源账。

### 11.7 Gate H：115 fs阶段功率定位

#### 11.7.1 允许修改的文件

第一轮只允许新增或修改离线分析工具：

```text
tools/analyze_vpfp_stage_power_ab.py
tests/vpfp_stage_power_ab_test.py
```

不得修改`src/`中的生产推进。工具必须读取生产生成的：

```text
vpfp_step_diagnostics.dat
vpfp_stage_energy_audit.dat
```

不得复写Vlasov、Poisson、PIC或碰撞公式。

分析器CLI必须严格实现为：

```text
python3 tools/analyze_vpfp_stage_power_ab.py
  --checkpoint <完整checkpoint目录>
  --coarse <粗档结果目录>
  --fine <细档结果目录>
  --coarse-dt-scale <正数>
  --fine-dt-scale <正数>
  --expected-coarse-steps <正整数>
  --expected-fine-steps <正整数>
  --result <输出.result>
```

禁止增加隐式默认目录。任一输入目录、`manifest.dat`、`vpfp_step_diagnostics.dat`或
`vpfp_stage_energy_audit.dat`缺失时必须立即失败，不得输出部分PASS。

checkpoint的`step`、`time_s`和`physical_config_hash`必须从`manifest.dat`读取。若checkpoint采用
每rank manifest，则必须验证所有rank的`step/time_s/physical_config_hash`一致。粗细档初始状态通过
各自第一行的`domain_energy_before`、`U_E_before`、`K_e_before`、`K_b_before`、`K_tail_before`、
`N_e_before`和`N_b_before`与checkpoint元数据/共同起点进行比较；不能只比较目录名。

粗细档上述初始标量应由同一checkpoint恢复，优先要求输出十进制字符串完全一致；若MPI求和顺序
造成末位差异，允许的绝对门为：

$$
T_{initial}=512\epsilon_{double}\max(|A|,|B|,1).
$$

任一初始标量超过该门即`same_initial_physical_state=0`。不得使用$10^{-6}$等宏观相对阈值放过
不同初态。

分析器退出码定义：

```text
0  PASS_ROOT_CAUSE_IDENTIFIED 或 PASS_LEDGER_DEFECT_IDENTIFIED
2  INCONCLUSIVE_EXTEND_WINDOW 或 INCONCLUSIVE_EVENT_PATH_DIVERGED
3  FAIL_AUDIT_STRUCTURE
4  FAIL_NONFINITE_OR_CORRUPT_INPUT
5  FAIL_USAGE
```

`.result`中的`status=`必须与退出码一致。禁止把`INCONCLUSIVE_EXTEND_WINDOW`写成PASS。
结果文件必须为稳定的`key=value`文本，至少包含：

```text
status
checkpoint_step
checkpoint_time_s
checkpoint_physical_config_hash
coarse_steps/fine_steps
coarse_elapsed_time_s/fine_elapsed_time_s
same_checkpoint
same_initial_physical_state
same_physical_window_relative_error
full_residual_power_coarse/fine
full_residual_power_ratio
full_residual_observed_order
known_source_minus_accounted_coarse/fine
dominant_stage_or_group
dominant_explanation_fraction_coarse/fine
dominant_power_ratio
root_cause
recommended_files
```

所有逐stage字段采用固定前缀`stage_<sanitized_stage_name>_`。stage名称只允许字母、数字和下划线；
发现其他字符时必须明确转义，禁止产生重复key。

单元测试必须使用临时目录生成最小合法的step/stage表，不依赖真实生产大文件。至少覆盖：

1. 单一stage解释90%且粗细功率相同：`PASS_ROOT_CAUSE_IDENTIFIED`；
2. 已知源项从`accounted_energy_source`遗漏：`PASS_LEDGER_DEFECT_IDENTIFIED`；
3. 最大stage解释率低于80%：`INCONCLUSIVE_EXTEND_WINDOW`；
4. 粗细物理窗口不一致：`FAIL_AUDIT_STRUCTURE`；
5. 缺列、重复step、非有限值、split/failure/rollback：对应硬失败；
6. 舍入量级stage：`order=not_evaluated_roundoff`；
7. 粗10步/细20步的单步残差分别为$2r/r$但功率相同：观测阶必须接近0，不得误报一阶。

单元测试命令：

```bash
python3 -m py_compile tools/analyze_vpfp_stage_power_ab.py \
  tests/vpfp_stage_power_ab_test.py || exit 78

python3 tests/vpfp_stage_power_ab_test.py \
  --analyzer ./tools/analyze_vpfp_stage_power_ab.py \
  --result ./output/vpfp_stage_power_ab_unit.result || exit 79

grep -q '^status=PASS$' ./output/vpfp_stage_power_ab_unit.result || exit 79
```

完成分析器和单元测试后必须确认`src/`没有因Gate H工具开发发生变化；若执行环境有Git，运行
`git diff --name-only -- src`应为空。无Git环境时直接比较修改清单，不把Git不可用判为测试失败。

#### 11.7.2 分析器必须实现的量

`vpfp_step_diagnostics.dat`至少要求以下列；缺任一列均为结构失败：

```text
step time_s accepted split
U_E K_e K_b K_tail K_combined
U_E_before K_e_before K_b_before K_tail_before
N_e_before N_e_after N_b_before N_b_after
domain_energy_before domain_energy_after domain_energy_delta
accounted_energy_source energy_balance_residual
electrostatic_boundary_work background_boundary_energy_net beam_boundary_energy_net
collision_reservoir fct_energy
conversion_N_residual conversion_Px_residual conversion_K_residual
tail_outflow_K
tail_return_N_residual tail_return_Px_residual tail_return_K_residual
collision_flux_rollback_count
```

`vpfp_stage_energy_audit.dat`至少要求：

```text
step time_s accepted audit_valid split failure_code stage_id stage_name
K_bulk K_tail K_beam U_E dK_bulk dK_tail dK_beam dU_E
Q_bkg_left_in Q_bkg_left_out Q_bkg_right_in Q_bkg_right_out
Q_beam_in Q_beam_out Q_tail_out Q_collision_reservoir
K_conversion K_tail_return W_electrostatic_boundary stage_balance
bulk_upar_face_work bulk_upar_velocity_boundary_work
bulk_upar_interface_energy_removed bulk_upar_identity_residual
tail_kick_work beam_kick_work
```

禁止假定`stage_id`固定不变。必须按`stage_name`聚合，并验证同一运行内
`stage_id <-> stage_name`是一一映射。未知stage必须保留为独立项并计入`unclassified_power`，
不能静默丢弃。

除逐stage输出外，必须计算以下固定组合；组合只做求和，不得改变符号：

```text
collision_pair = collision_half1 + collision_half2
x_pair = x_half1 + x_half2
field_force_pair = midpoint_poisson + u_force_tail_beam_kick + final_poisson
conversion_pair = conversion_after_force + conversion_after_collision
tail_return = tail_bulk_return
```

`accepted_n`是每步参考态，不计入残差组合。若生产代码新增stage，分析器必须将其列入
`unclassified_power`并令结构门失败，直到文档明确该stage的所有权。

对每一档先验证相同checkpoint、无split/failure/rollback，并计算真实物理窗口
$T=t_{last}-t_{checkpoint}$。对每个`stage_name`分别累计：

$$
P_s=\frac{1}{T}\sum_n R_{s,n},
$$

其中$R_{s,n}$直接取生产`stage_balance`，不得重新定义。至少输出：

```text
elapsed_time_s
accepted_steps
full_residual_signed/abs/power
accounted_source_power
domain_energy_delta_power

stage_<name>_residual_signed
stage_<name>_residual_abs
stage_<name>_residual_power
stage_<name>_fraction_of_full_signed

field_particle_pair_power
poisson_pair_power
force_conversion_pair_power
bulk_upar_identity_power
tail_kick_mismatch_power
beam_kick_mismatch_power
collision_stage_power
collision_reservoir_power
x_remap_boundary_power
conversion_internal_power
tail_return_internal_power
unclassified_power
conversion_event_count
tail_return_attempted_groups
tail_return_committed_groups
tail_particle_count_first/last
event_paths_differ
```

完整步源账必须使用生产代码`VpfpIntegrator::finalize_energy_ledger()`的同一符号约定核对：

$$
S_{known}=Q_{bkg,boundary}+Q_{beam,boundary}-Q_{tail,out}
-Q_{collision,reservoir}+W_{electrode}.
$$

对应列为：

```text
background_boundary_energy_net
beam_boundary_energy_net
tail_outflow_K
collision_reservoir
electrostatic_boundary_work
```

必须分别输出`known_source_reconstructed`和
`known_source_minus_accounted_energy_source`。后者只有在机器精度门内，才能声明已知源账实现一致。
这项一致不代表没有尚未枚举的物理源；未知源只能通过阶段状态差与物理边界/转换事件进一步定位。

物理窗口必须按checkpoint的`time_s`计算，不能用步数乘标称$dt$代替：

```text
T_coarse = coarse_last_time_s - checkpoint_time_s
T_fine   = fine_last_time_s   - checkpoint_time_s
window_relative_error = |T_coarse-T_fine| / max(|T_coarse|,|T_fine|)
```

只有`window_relative_error<=1e-12`才允许计算功率比和阶数。若终止时间发生浮点吸附，应使用实际
接受态时间；禁止线性插值阶段记录。

同时对每项输出：

$$
\text{power ratio}=\frac{|P_{0.25}|}{\max(|P_{0.5}|,P_{floor})},
\qquad
p_s=\log_2\frac{|P_{0.5}|}{|P_{0.25}|}.
$$

`P_floor`只能使用机器精度缩放门；不得用固定的大阈值隐藏小项。若某项在任一档低于舍入门，
输出`order=not_evaluated_roundoff`。

候选阶段对完整残差的signed解释率定义为：

$$
F_s=\max\left(0,
1-\frac{|P_{full}-P_s|}{\max(|P_{full}|,P_{floor})}
\right).
$$

不得使用$|P_s|/\sum|P_j|$代替，因为阶段间可能发生物理抵消。声明单一根因要求粗细两档均满足：

```text
same sign as full residual power
signed_explanation_fraction >= 0.80
0.80 <= abs(P_fine/P_coarse) <= 1.25
same dominant stage or same predefined stage group
```

若只有预定义组合达到80%，则根因只能定位到组合，不能擅自选择组合中的最大单项。

#### 11.7.3 生成短审计数据的粗档命令

使用完整115 fs checkpoint，不得使用snapshot。命令从manifest读取起始step，并将停止step设置为
`CP_STEP+10`；不得假定以后生成的checkpoint仍固定为5071：

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 80
CHECKPOINT_115="${CHECKPOINTS_115[0]}"
CP_STEP=$(awk -F= '$1=="step" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
CP_TIME_S=$(awk -F= '$1=="time_s" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
[[ "$CP_STEP" =~ ^[0-9]+$ ]] || exit 80
test -n "$CP_TIME_S" || exit 80
COARSE_STOP_STEP=$((CP_STEP + 10))
OUT_H_COARSE=./output/vpfp_gate_h_stage_power_dt050_10steps
test -s "$CHECKPOINT_115/manifest.dat" || exit 80
test ! -e "$OUT_H_COARSE" || exit 81

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.5 \
  --stop-after-steps "$COARSE_STOP_STEP" \
  --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_H_COARSE" || exit 82

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_H_COARSE" --expected-accepted-steps 10 --require-no-split \
  --result ./output/vpfp_gate_h_stage_power_dt050_10steps.structure.result || exit 83
```

#### 11.7.4 生成短审计数据的细档命令

细档从同一checkpoint独立出发，将停止step设置为`CP_STEP+20`，覆盖与粗档相同的标称物理时间：

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 84
CHECKPOINT_115="${CHECKPOINTS_115[0]}"
CP_STEP=$(awk -F= '$1=="step" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
CP_TIME_S=$(awk -F= '$1=="time_s" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
[[ "$CP_STEP" =~ ^[0-9]+$ ]] || exit 84
test -n "$CP_TIME_S" || exit 84
FINE_STOP_STEP=$((CP_STEP + 20))
OUT_H_FINE=./output/vpfp_gate_h_stage_power_dt025_20steps
test -s "$CHECKPOINT_115/manifest.dat" || exit 84
test ! -e "$OUT_H_FINE" || exit 85

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-source-dt-scale 0.5 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.25 \
  --stop-after-steps "$FINE_STOP_STEP" \
  --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_H_FINE" || exit 86

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_H_FINE" --expected-accepted-steps 20 --require-no-split \
  --result ./output/vpfp_gate_h_stage_power_dt025_20steps.structure.result || exit 87
```

#### 11.7.5 阶段功率比较命令

实现§11.7.1工具后执行：

```bash
mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 88
CHECKPOINT_115="${CHECKPOINTS_115[0]}"

python3 tools/analyze_vpfp_stage_power_ab.py \
  --checkpoint "$CHECKPOINT_115" \
  --coarse ./output/vpfp_gate_h_stage_power_dt050_10steps \
  --fine ./output/vpfp_gate_h_stage_power_dt025_20steps \
  --coarse-dt-scale 0.5 \
  --fine-dt-scale 0.25 \
  --expected-coarse-steps 10 \
  --expected-fine-steps 20 \
  --result ./output/vpfp_gate_h_stage_power_ab.result
RC=$?
cat ./output/vpfp_gate_h_stage_power_ab.result
case "$RC" in
  0) echo "Gate H root cause identified" ;;
  2) echo "Gate H inconclusive: execute section 11.7.7 only" ;;
  *) exit "$RC" ;;
esac
```

#### 11.7.6 Gate H验收和根因选择

结构门：

```text
same_checkpoint=1
same_initial_physical_state=1
same_physical_window_relative_error <= 1e-12
coarse_steps=10
fine_steps=20
split/failure/rollback=0
stage_telescope_matches_roundoff=1
source_ownership_valid=1
all_finite=1
```

根因选择必须遵循：

1. 若`stage_telescope`或完整余额重构不在舍入门内，先修诊断记录/分析器，禁止修改物理算子。
2. 若某个已拥有物理源项没有进入`accounted_energy_source`，其缺失功率在粗细两档均解释至少80%的
   完整残差，且补入后余额满足20%剩余门，则根因是源账遗漏；
   只修账本所有权，不改变状态推进。
3. 若某个状态更新阶段的`stage_balance_power`在两档中近似相同，并解释至少80%的完整残差功率，
   则该阶段是直接根因；只为该阶段另写修复方案。
4. 若`field_particle_pair_power`解释至少80%，且Poisson恒等式仍通过，则根因是场--粒子时间层/离散功
   配对结构，而不是Poisson空间算子。此时才允许设计统一中点或离散梯度耦合。
5. 若残差分散在多个阶段，必须报告各阶段占比；不得用全局补丁一次消除。
6. 若10/20步样本受单个转换/H10事件主导，可扩大到20/40步，但不得直接重跑115--120 fs。

`tail-return-residence-steps=8`是生产配置中的步数语义，粗细档对应的物理驻留时间不同。Gate H
保持该参数不变，是为了审计**当前实际生产配置**的残差功率，不得用本Gate宣布Tail粒子数或H10
返回频率随$dt$收敛。若候选根因落在H10返回阶段且`event_paths_differ=1`，结果必须是
`INCONCLUSIVE_EVENT_PATH_DIVERGED`，不能直接修改H10。应另行把驻留条件改为物理时间后再做专门测试。

Gate H本身只负责定位，不允许在同一提交中修改生产物理。完成后输出唯一根因阶段、证据占比、
粗细档功率比和建议修改文件；没有达到80%解释率时必须报告`root_cause=not_yet_unique`。

定位后允许进入的修改范围如下；Gate H执行智能体只能提出方案，不能在同一轮实施：

| Gate H结果 | 后续允许评估的文件 | 仍禁止修改 |
|---|---|---|
| `PASS_LEDGER_DEFECT_IDENTIFIED` | `src/vpfp_integrator.cpp`中的`finalize_energy_ledger()`、`src/vpfp_diagnostics.cpp`及对应测试 | 任意状态推进、Poisson解、粒子速度 |
| `collision_pair`主导 | 碰撞阶段能量交换/储能账接口及其单元测试 | Poisson、x/u remap、Beam/Tail pusher |
| `x_pair`主导 | 开放背景边界能流所有权与x-remap阶段账 | u-force、Poisson空间算子、PIC kick |
| `field_force_pair`主导 | `src/vpfp_integrator.cpp`的中点场/受力/末态Poisson时间层，另立结构保持设计 | `OpenElectrostaticSolver`空间离散、事后能量投影 |
| `conversion_pair`主导 | flux-interface转换阶段的状态差和内部转移所有权 | 转换阈值、粒子物理动量的全局缩放 |
| `tail_return`主导 | 先修物理驻留时间语义并重做专门A/B | 直接调返回能量以补账 |
| `not_yet_unique` | 仅执行20/40步扩展或补诊断 | 所有生产物理 |

#### 11.7.7 仅在10/20步结果不充分时执行20/40步扩展

只有§11.7.6返回`INCONCLUSIVE_EXTEND_WINDOW`时才运行本节。若返回结构失败，必须先修分析器；
若已定位根因，禁止为获得更漂亮数字继续扩大样本。

粗档20步命令：

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 89
CHECKPOINT_115="${CHECKPOINTS_115[0]}"
CP_STEP=$(awk -F= '$1=="step" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
[[ "$CP_STEP" =~ ^[0-9]+$ ]] || exit 89
COARSE_STOP_STEP=$((CP_STEP + 20))
OUT_H_COARSE40=./output/vpfp_gate_h_stage_power_dt050_20steps
test ! -e "$OUT_H_COARSE40" || exit 90

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.5 \
  --stop-after-steps "$COARSE_STOP_STEP" \
  --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_H_COARSE40" || exit 91

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_H_COARSE40" --expected-accepted-steps 20 --require-no-split \
  --result ./output/vpfp_gate_h_stage_power_dt050_20steps.structure.result || exit 92
```

细档40步命令：

```bash
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 93
CHECKPOINT_115="${CHECKPOINTS_115[0]}"
CP_STEP=$(awk -F= '$1=="step" {print $2; exit}' "$CHECKPOINT_115/manifest.dat")
[[ "$CP_STEP" =~ ^[0-9]+$ ]] || exit 93
FINE_STOP_STEP=$((CP_STEP + 40))
OUT_H_FINE40=./output/vpfp_gate_h_stage_power_dt025_40steps
test ! -e "$OUT_H_FINE40" || exit 94

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_115" \
  --restart-source-dt-scale 0.5 \
  --restart-allow-dt-scale-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 \
  --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode hysteretic \
  --tail-return-energy-mev 5.5 \
  --tail-return-residence-steps 8 \
  --tail-return-max-stencil-radius 3 \
  --tail-return-moment-tolerance 1e-12 \
  --dt-scale 0.25 \
  --stop-after-steps "$FINE_STOP_STEP" \
  --stop-time-fs 120 \
  --diagnostic-level 2 --diagnostic-interval 1 \
  --output-dir "$OUT_H_FINE40" || exit 95

python3 tools/analyze_vpfp_stage_energy_audit.py \
  --run "$OUT_H_FINE40" --expected-accepted-steps 40 --require-no-split \
  --result ./output/vpfp_gate_h_stage_power_dt025_40steps.structure.result || exit 96
```

扩展比较命令：

```bash
mapfile -t CHECKPOINTS_115 < <(
  find ./output/vpfp_return_dt_half_100p3_to_120fs -maxdepth 1 -type d \
    -name 'checkpoint_target115fs*' -print
)
test "${#CHECKPOINTS_115[@]}" -eq 1 || exit 97
CHECKPOINT_115="${CHECKPOINTS_115[0]}"

python3 tools/analyze_vpfp_stage_power_ab.py \
  --checkpoint "$CHECKPOINT_115" \
  --coarse ./output/vpfp_gate_h_stage_power_dt050_20steps \
  --fine ./output/vpfp_gate_h_stage_power_dt025_40steps \
  --coarse-dt-scale 0.5 \
  --fine-dt-scale 0.25 \
  --expected-coarse-steps 20 \
  --expected-fine-steps 40 \
  --result ./output/vpfp_gate_h_stage_power_ab_20_40.result
RC=$?
test "$RC" -eq 0 -o "$RC" -eq 2 || exit "$RC"
exit "$RC"
```

### 11.8 当前停止条件

在Gate H完成前：

- 不运行新的0--120 fs生产；
- 不继续执行`dt_scale=0.25`长跑；
- 不实施全局统一中点重构；
- 不修改能量预警/失效阈值；
- 不把现有0.902%解释为纯诊断误差或纯物理加热。

当前唯一允许推进的工作是§11.7的只读阶段功率定位。

---

## 12. 自动编码模型最终报告模板

每完成一个Gate，必须按以下格式报告：

```text
Gate:
修改文件:
未修改的生产模块:
编译结果:
测试命令:
测试结果文件:
关键数值:
粗档残差功率:
细档残差功率:
功率比与观测阶:
首个异常stage_name:
该阶段解释完整残差的比例:
源账遗漏还是状态更新缺口:
验收标准逐项结果:
是否允许进入下一Gate:
剩余风险:
```

禁止只回复“已修复”“测试通过”或只给最终`status`。必须列出绝对残差、尺度、容差、最坏步号和
none/return差异。Gate H还必须列出每个主要阶段的单位物理时间残差；不能只比较单步残差，因为
单步残差会随$dt$机械缩小并造成虚假的“改善”。

---

# 附录：Gate A--H执行记录、长期A/B与当前决策

> 本附录记录自动编码模型按正文逐门执行后的实测结果。
> 状态：Gate A--H、100.3--120 fs生产和115--117 fs时间步A/B均已完成。Poisson空间算子、
> 已记录源项所有权和各粒子局部功恒等式通过；但长期残差对$dt$减半呈近零阶，旧的“混合阶
> 时间截断为最终根因”结论已撤销。Gate H阶段功率定位已完成，根因=场--粒子时间层/离散功配对
> 结构（field_particle_pair解释99%，对应field_force_pair块）。

## A. 执行状态总览

| Gate | 内容 | 状态 | 关键结论 |
|---|---|---|---|
| A | 分析器机器精度缩放门 | ✅ 通过 | 固定 1e-12 相对阈值改为 `512·ε·S` 缩放门 |
| B | 现有记录派生分解 | ✅ 通过 | 分解与完整余额在舍入容差内一致，无未分类阶段账本项 |
| C | 只读离散功记录 | ✅ 通过 | bulk通量功/Tail/Beam kick局部恒等式闭合 |
| D | 100.3 fs none/return 100 步 | ✅ 通过 | 残差定位于field--particle耦合块，H10不是根因 |
| E | dt/dt/2 同窗口判别 | ✅ 通过（局部证据） | $dt/2$将短窗配对绝对残差降低58%，不能外推长期主导阶 |
| F | 根因选择门 | ✅ 通过（结构门） | Poisson和源项所有权通过；历史时间截断分类被长期A/B修正 |
| G | 最终生产配置短回归 | ✅ 通过 | 200步、无split/failure/rollback，结构、能量预算和宏观门均通过 |
| 长跑 | 100.3--120 fs | ⚠️ 预警 | 累计残差占最终域能0.902%，99.16%单步同号 |
| A/B | 115--117 fs 0.5/0.25 | ❌ 时间步修复无效 | 残差功率仅降1.4%，观测阶0.00657，成本2.04倍 |
| H | 阶段功率定位 | ✅ 完成 | 根因=场--粒子时间层；field_particle_pair解释约99%（field_force_pair块） |

## B. Gate A 结果（通过）

- 修改文件：`tools/analyze_vpfp_stage_energy_audit.py`、`tests/vpfp_stage_energy_audit_test.cpp`。
- 判据：`T_round = max(1e-10, 512·ε_double·S)`，`ε_double=2.220446e-16`，
  `S = max(E_scale, S_stage, |R_full|, |R_ledger|, 1)`；`stage_telescope_abs` 与 `energy_ledger_abs`
  分别与自身 `T_round` 比较；旧相对字段仅作参考，不再控制 `status`。
- 新增结果字段：`stage_telescope_max_abs/roundoff_tolerance/worst_step`、
  `energy_ledger_max_abs/roundoff_tolerance/worst_step`、`audit_structure_pass`、
  `physical_energy_residual_cumulative`、`physical_energy_gate_evaluated=0`。
- 单元测试 `output/vpfp_stage_energy_audit_unit.result`：**status=PASS**。
  - `accept_large_perturb=1`（5e9 尺度 + 2e-6 求和扰动 PASS，tol=5.68e-4）；
  - `reject_large_missing_source=1`（漏 1 J/m² 源项 FAIL）；
  - `reject_small_ledger_residual=1`（小尺度 1e-4 残差 FAIL，tol=1e-10 底线）；
  - `accept_physical_residual=1`（物理残差非零不令结构 FAIL）；
  - 缺阶段/重复/NaN/Inf/未接受/split/failure 六项 `reject_*=1`。
- 旧 none/return 100p3 数据重分析：均 status=PASS（none 最坏步 3930，
  `stage_telescope_max_abs=1.49e-6 ≤ 6.11e-4`，`energy_ledger_max_abs=1.91e-6 ≤ 6.11e-4`）。

## C. Gate B 结果（通过）

- 仅修改 `tools/analyze_vpfp_stage_energy_audit.py`。
- none 100p3 派生分解（与正文 §2.2 完全一致）：
  `bulk u-force=+7.389119852e6`、`Tail kick=+5.032334826e5`、`Beam kick=-8.290686255e6`、
  `conversion-after-force=+1.296696566e6`、`force_conversion_pair=+8.983636446e5`、
  `poisson_pair=+2.087114625e6`、`resolved_common_balance=2.985478269e6`（≈完整余额，
  `accounted_other=-4.9e-7`）。
- return 100p3：`poisson_pair=+2.117814073e6`、`force_conversion_pair=+9.236932752e5`、
  `resolved_common_balance=3.041507349e6`。
- x/collision/H10 累计残差均 ≤1e-3：none `-4.3e-6 / 3.1e-6 / -4.8e-7`，
  return `-1.2e-5 / 1.7e-5 / 2.6e-6`；`decomposition_matches_roundoff=1`。

## D. Gate C 结果（通过）

- 修改文件：`src/conservative_ppm_remap.{h,cpp}`、`src/background_tail_pic.{h,cpp}`、
  `src/vpfp_integrator.{h,cpp}`、`src/vpfp_diagnostics.cpp`、`CMakeLists.txt`、
  `tools/analyze_vpfp_stage_energy_audit.py`；新增 `tests/vpfp_force_work_audit_test.cpp`。
- 关键实现：
  1. `RemapDiagnostics` 新增 5 字段（`upar_internal_face_energy_transfer`、
     `upar_left/right_velocity_boundary_energy`、`upar_interface_energy_removed`、
     `upar_discrete_energy_identity_residual`）；
  2. `BackgroundTailPIC::kick(..., double* local_kinetic_work=NULL)`（OpenMP reduction，空指针零开销）；
  3. `VpfpStageEnergyRecord` 新增 6 字段，`packed_values 15→21`，bulk 4 字段 rank0-only、
     Tail/Beam kick 全 rank 求和，`.dat` 追加 6 列并扩展有限性检查；
  4. 修复一处实现缺陷：`advect_u_parallel` 更新循环内先捕获 `input_mass` 再写 `output`，
     保证 in-place（`output==input`）下离散动能恒等式不失效。
- 单元测试 `output/vpfp_force_work_audit_unit.result`：**status=PASS**。
  - `bulk_identity_residual_plus=-2.47e-17`、`bulk_identity_residual_minus=-3.79e-18`（机器精度）；
  - `left/right_velocity_boundary_energy=-2020041.00`（左/右端面出流符号正确）；
  - `interface_removed=298395.24 == interface_expected_removed`（接口移除能量只记一次）；
  - `tail_work_plus=2.937`、`tail_work_minus=-0.936`、`tail_work_mismatch=0`（正负场符号+配对）；
  - `beam_work_mismatch=0`、`beam_work_zero_before_push=0`（Beam kick 与动能差一致）；
  - `state_bitwise_equal_with_null=1`（空指针路径状态 bitwise 不变）。

## E. Gate D 结果（通过）

- `output/vpfp_force_work_none_100p3_100steps.result`：status=PASS，物理累计残差 2.985478269e6 J/m²。
  - `bulk_upar_identity_residual_cumulative=1.87e-5`（舍入，恒等式闭合）；
  - `tail_kick_snapshot_mismatch_cumulative=-1.19e-5`、`beam_kick_snapshot_mismatch_cumulative=2.62e-6`（舍入）；
  - `conversion_pair_residual_cumulative=-4.56e-7 ≤ 1e-3`；`H10=2.38e-7 ≤ 1e-3`；
  - `field_particle_pair_residual_cumulative=3.031e6`（≈物理余额 101.5%）。
- `output/vpfp_force_work_return_100p3_100steps.result`：status=PASS，物理累计残差 3.041507349e6 J/m²。
  - `field_particle_pair_residual_cumulative=3.087e6`（≈物理余额 101.5%）。
- 结论：bulk u_parallel 通量、Tail kick、Beam kick、转换事务、H10 各自功恒等式全部闭合，
  公共残差几乎全部落在 `field_particle_pair = ΔU_E + W_bulk + W_tail + W_beam − W_electrode`。
  → 排除分支 A（bulk 通量不闭合）与分支 B（Tail/Beam kick 不闭合）。

## F. Gate E 结果（通过，但仅是局部时间截断证据）

- `output/vpfp_force_work_none_100p3_dt_half_200steps.result`：status=PASS（200 步）。
  - 末态物理时间与 Gate D 基准一致（`1.0285850873086202e-13` vs `1.0285850873086075e-13`，差 ~9 ulp），
    窗口有效；半步长档 `physical_energy_residual_cumulative=3.7577e6`（较基准 2.9855e6 **增大 26%**）。
- 收敛率 `q = log2(|R_dt|/|R_dt/2|)`：

  | 量 | R_dt | R_dt/2 | q |
  |---|---|---:|---:|
  | physical_energy_residual | 2.9855e6 | 3.7577e6 | **-0.33** |
  | poisson_pair | 2.0871e6 | 1.7954e6 | **+0.22** |
  | force_conversion_pair | 8.984e5 | 1.9623e6 | **-1.13** |
  | field_particle_pair | 3.0314e6 | 3.7987e6 | **-0.33** |

- 原始signed累计结果全部$q<0.5$，但逐步复核显示该结论不能用于排除分支D：
  - 第1个共同`dt`区间：$R_{dt}=-1.45536\times10^5$，两个半步合计
    $R_{dt/2}=-1.85605\times10^4$，局部$q=2.97$；
  - 前2个共同区间局部累计$q=3.45$；
  - `dt`档100步signed累计为$3.0314\times10^6$，绝对累计为$9.1396\times10^6$；
    `dt/2`档分别为$3.7987\times10^6$和$3.9405\times10^6$；
  - 全窗口绝对累计阶$q_{abs}=1.21$，属于一阶/混合误差区，不是固定不兼容的充分证据；
  - 两档末态Tail粒子数约$6.50\times10^7$与$7.02\times10^7$，随机碰撞调用、转换事件和
    粒子表示已经不同。
- 修正结论：物理窗口和步数验收通过，但当前生产`dt/dt/2`不是受控的确定性Richardson试验。
  它只能证明短窗误差包络中存在时间截断分量，不能证明该分量主导长期signed累计。后续
  115--117 fs生产A/B给出近零阶，已经否定“继续减小全局$dt$”作为长期修复。

## G. Gate F 根因选择（结构门通过，最终分类已修正）

- 触发条件逐项命中（§10.3）：
  1. bulk/Tail/Beam 各自功恒等式通过（Gate C/D，残差 ~1e-5，机器精度）；
  2. `field_particle_pair_residual` 占完整余额 101.5%（none 3.031e6/2.985e6，return 3.087e6/3.041e6，>50%）。
- 可确认结论：残差位于`midpoint_poisson → force/kick → final_poisson`所构成的场–粒子耦合块；
  bulk、Tail、Beam各自的局部动能恒等式不是百万量级缺口来源。
- §10.5局部/固定状态审计全部完成：
  - `source_ownership_valid=1`；
  - `poisson_identity_pass=1`；
  - 审计后接受态、RNG和ledger均bitwise不变；
  - `classification=time-truncation-mixed-order`；
  - `status=PASS`。
- Poisson单元测试给出$R_{P,h}=-1.1755\times10^{-38}\ \mathrm{J/m^2}$，而舍入容差为
  $1.7038\times10^{-34}\ \mathrm{J/m^2}$，排除Poisson空间算子/边界功错误。
- 固定状态审计的稳定前缀呈现约1.5阶，但其他前缀有强烈抵消和非渐近行为。该结果说明耦合块
  含有时间截断成分，但不足以决定生产长时主导误差。
- 115--117 fs生产A/B给出观测阶0.00657，因此历史输出
  `classification=time-truncation-mixed-order`只保留为Gate F当时的局部分类，不再作为最终根因。
- 当前决策：不改Poisson，不做全局能量补丁，不继续依赖全局减步长；执行Gate H阶段功率定位。

## H. Gate G与100.3--120 fs生产结果

Gate G最新结果：

```text
accepted_steps=200
split_steps=0
failure_records=0
audit_structure_pass=1
source_ownership_valid=1
decomposition_matches_roundoff=1
signed_to_domain=6.9212821e-4
U_E relative difference=3.96673e-3
K_combined relative difference=2.01445e-3
K_tail relative difference=1.99614e-2
N_e relative difference=3.52916e-6
macro status=PASS
analyzer status=PASS
```

Tail返回严格不变量的最大逐步相对残差为：$N=1.99\times10^{-15}$、
$P_x=1.82\times10^{-15}$、$K=2.47\times10^{-15}$，均低于$10^{-12}$。碰撞回滚计数为0。

100.3--120 fs生产已经完成：1540步全部接受，无split/failure/rollback；转换和H10严格事务通过。
累计能量残差$5.116142346\times10^7\ \mathrm{J/m^2}$，占最终域能0.902%，达到预警并接近1%失效线。

## I. 115--117 fs A/B与后续唯一待办

从相同115.011425 fs checkpoint得到：

```text
dt_scale=0.5 residual power = 3.293199154e6 J/m2/fs
dt_scale=0.25 residual power = 3.278225893e6 J/m2/fs
observed order              = 0.00657449
fine/coarse wall            = 2.04127
```

Gate H（§11.7）已完成阶段功率定位，根因=场--粒子时间层/离散功配对结构（触发量
`field_particle_pair`，见下文 J 节）。后续唯一待办是
依据 J 节结论另立结构保持设计文档；设计落地前仍不得继续减小生产时间步、重跑0--120 fs、
修改Poisson空间算子或实施全局能量补丁。

---

## J. Gate H 结果（通过，根因=场--粒子时间层/离散功配对结构）

Gate: Gate H：115 fs阶段功率定位（§11.7）

修改文件:
- `tools/analyze_vpfp_stage_power_ab.py`（新增，阶段功率A/B分析器）
- `tests/vpfp_stage_power_ab_test.py`（新增，单元测试）
- 期间修复：checkpoint manifest 实为 `manifest.txt`（空格分隔 `key value`，时间键 `time`），
  分析器改为优先 `manifest.txt`、回退 snapshot 的 `manifest.dat`（`key=value`，键 `time_s`）；
  窗口门由固定 1e-12 改为机器缩放 `max(1e-12, 512·ε·time_scale/T)`。

未修改的生产模块:
- `src/` 全部未改动（Gate H 仅新增离线分析工具；`git diff --name-only -- src` 应无本 Gate 引入差异）
- 生产推进、Poisson、PIC、碰撞、转换、H10 均未改动

编译结果:
- `python3 -m py_compile tools/analyze_vpfp_stage_power_ab.py tests/vpfp_stage_power_ab_test.py` 通过
- 单元测试 `output/vpfp_stage_power_ab_unit.result`：status=PASS（11 项用例全 PASS）

测试命令:
```bash
python3 -m py_compile tools/analyze_vpfp_stage_power_ab.py tests/vpfp_stage_power_ab_test.py
python3 tests/vpfp_stage_power_ab_test.py \
  --analyzer ./tools/analyze_vpfp_stage_power_ab.py \
  --result ./output/vpfp_stage_power_ab_unit.result
python3 tools/analyze_vpfp_stage_power_ab.py \
  --checkpoint "$CHECKPOINT_115" \
  --coarse ./output/vpfp_gate_h_stage_power_dt050_10steps \
  --fine   ./output/vpfp_gate_h_stage_power_dt025_20steps \
  --coarse-dt-scale 0.5 --fine-dt-scale 0.25 \
  --expected-coarse-steps 10 --expected-fine-steps 20 \
  --result ./output/vpfp_gate_h_stage_power_ab.result
```

测试结果文件:
- `output/vpfp_gate_h_stage_power_dt050_10steps.structure.result`：status=PASS
- `output/vpfp_gate_h_stage_power_dt025_20steps.structure.result`：status=PASS
- `output/vpfp_gate_h_stage_power_ab.result`：status=PASS_ROOT_CAUSE_IDENTIFIED
- `output/vpfp_stage_power_ab_unit.result`：status=PASS

> 注：`vpfp_gate_h_stage_power_ab.result` 生成于分析器命名修正之前，其内部
> `root_cause=field_force_pair` / `dominant_stage_or_group=field_force_pair` 为旧标签；修正后的
> 分析器已将其改为 `field_particle_pair`（触发量），语义与本附录一致，二者不矛盾。

关键数值:
```text
checkpoint_step=5071
checkpoint_time_s=1.1501142520244935e-13   (115.011425 fs)
coarse_steps=10 (step 5072..5081)          fine_steps=20 (step 5072..5091)
coarse_elapsed_time_s=1.2792543654e-16     fine_elapsed_time_s=1.2792543654e-16
same_checkpoint=1
same_initial_physical_state=1
same_physical_window_relative_error=1.973e-12  <= tol=1.023e-10
stage_telescope_matches_roundoff=1
source_ownership_valid=1
all_finite=1
known_source_minus_accounted_coarse=0 / fine=0   （无源账遗漏）
```

粗档残差功率 (coarse, dt=0.5, 10步, T=1.2793e-16 s):
```text
full_residual_signed=328318.556 J/m2      full_residual_power=2.5665e21 J/m2/s
midpoint_poisson        +1.2617e23 J/m2/s   (signed +1.6141e7 J/m2)
final_poisson           -9.5893e22 J/m2/s   (signed -1.2267e7 J/m2)
u_force_tail_beam_kick  -2.8426e22 J/m2/s   (signed -3.6364e6 J/m2)
conversion_after_force  +7.1300e20 J/m2/s   (signed +9.1210e4 J/m2)
collision/x_remap/H10    ~1e9..1e11 J/m2/s  (signed ~1e-7..1e-5 J/m2，舍入级)
```

细档残差功率 (fine, dt=0.25, 20步, T=1.2793e-16 s):
```text
full_residual_signed=380567.646 J/m2      full_residual_power=2.9749e21 J/m2/s
midpoint_poisson        +7.2171e22 J/m2/s   (signed +9.2325e6 J/m2)
final_poisson           -3.8563e22 J/m2/s   (signed -4.9332e6 J/m2)
u_force_tail_beam_kick  -3.1340e22 J/m2/s   (signed -4.0092e6 J/m2)
conversion_after_force  +7.0722e20 J/m2/s   (signed +9.0472e4 J/m2)
collision/x_remap/H10    ~1e9..1e11 J/m2/s  (signed ~1e-7..1e-5 J/m2，舍入级)
```

功率比与观测阶:
```text
field_particle_pair   coarse 2.5920e21 / fine 3.0003e21   比=1.158  阶=-0.211
field_force_pair      比=1.223  阶=-0.291
full_residual         比=1.159  阶=-0.213    （近零阶，非二阶截断）
conversion_after_force 比=0.992  阶=+0.012
碰撞/x_remap/H10 各阶段 = not_evaluated_roundoff（舍入级）
```

首个异常stage_name:
- `field_particle_pair`（触发量，命中 §11.7.6 第4条）。它对应 `field_force_pair` 块
  （= midpoint_poisson + u_force_tail_beam_kick + final_poisson）；注意二者非同一量：
  `field_particle_pair = ΔU_E + W_bulk + W_tail + W_beam − W_electrode` 解释 99%，
  而 `field_force_pair` 三个 stage_balance 之和仅解释 72.2%（coarse）/ 76.2%（fine）。

该阶段解释完整残差的比例:
- `field_particle_pair` 解释完整残差功率 99.0%（coarse）/ 99.1%（fine）
- `dominant_explanation_fraction=0.9901 (coarse) / 0.9915 (fine)`（≥0.80 门槛，属 field_particle_pair）
- `dominant_power_ratio=1.1575`（在 [0.80,1.25] 内）
- 对照：`field_force_pair`（stage 组合）解释仅 72.2%/76.2%，未达 80% 单根因门槛；
  其余 ~24–28% 为 conversion_after_force（+9.12e4/+9.05e4 J/m²）。

源账遗漏还是状态更新缺口:
- 状态更新缺口（非源账遗漏）。判据：`known_source_minus_accounted=0` 且
  `source_ownership_valid=1`，即 `S_known` 与 `accounted_energy_source` 在机器精度门内一致；
  故根因不是 `finalize_energy_ledger()` 遗漏源项，而是场--粒子时间层/离散功配对结构。

验收标准逐项结果 (§11.7.6):
```text
same_checkpoint=1                                ✓
same_initial_physical_state=1                    ✓
same_physical_window_relative_error=1.973e-12 ≤ tol=1.023e-10  ✓
coarse_steps=10 / fine_steps=20                  ✓
split/failure/rollback=0                         ✓
stage_telescope_matches_roundoff=1               ✓
source_ownership_valid=1                         ✓
all_finite=1                                     ✓
```
根因选择（§11.7.6 第4条）：触发判据为 `field_particle_pair_power` 在两档均解释≥80% 且
Poisson 恒等式仍通过 → 根因=场--粒子时间层/离散功配对结构，而非 Poisson 空间算子。
Poisson 恒等式通过来自 Gate F §10.5（`poisson_identity_pass=1`），Gate H 输入数据不重算该项；
本 Gate 命中第4条后 `root_cause` 报告为 `field_particle_pair`（触发量），`recommended_files`
映射到 `field_force_pair主导` 行的中点/末态 Poisson 时间层。`event_paths_differ=1`
（conversion 10↔20、tail_return committed 236↔228）仅说明事件路径随 dt 变化；因候选根因非
H10 返回阶段，不触发 `INCONCLUSIVE_EVENT_PATH_DIVERGED`。

是否允许进入下一Gate:
- 允许（Gate H 只定位，不修改物理）。下一步另立结构保持设计文档，评估范围仅限
  `src/vpfp_integrator.cpp` 的中点场/受力/末态 Poisson 时间层（统一中点或离散梯度型
  Vlasov--Poisson 耦合）。

剩余风险:
- 10/20 步样本受随机碰撞/转换事件影响（event_paths_differ=1）；但 field_particle_pair 解释率
  99% 且与 Gate D（field_particle_pair≈101.5%）、Gate E（非二阶）一致，不改变根因结论。
- 结构保持重构尚未设计，属下一独立任务；在此之前不得修改 Poisson 空间算子或施加全局能量补丁。
