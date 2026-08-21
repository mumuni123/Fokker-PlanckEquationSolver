# dual-u 当前结论、性能瓶颈与后续执行顺序

## 1. 文档状态

本文档以当前源码和以下结果目录为准：

```text
output/dual_u_perf_validation
output/dual_u_prod_reference/24_to_30_after_beam_fix
output/beam_ledger_io_ab
output/midpoint_accel_ab/24fs
```

已经完成且不再重复的工作：

- dual-u 阶段1--4的离散闭合、端点和MPI一致性测试；
- checkpoint/restart基础功能测试；
- 24 fs原时间步的40/80次迭代上限测试；
- `dt/2`的10步和500步严格收敛测试；
- `dt/2`的24--30 fs严格参考运行；
- Beam ledger的`full/summary`分离；
- 第一版场端点Aitken和深度3 Anderson；
- 24 fs中点加速A/B/C/D测试。

当前结论已经发生变化：

1. Beam ledger清理实现正确，但只有约`0.8%`性能收益；
2. 第一版Aitken和Anderson没有恢复原`dt`严格收敛，反而增加算子评估和wall time；
3. 当前不能继续30 fs中点加速测试，也不能启动原`dt`加速生产；
4. 当前可信生产基线仍是`dt/2 + midpoint-acceleration=none`；
5. 下一步优先寻找“仍可严格收敛的最大稳定时间步”，而不是继续盲调Aitken/Anderson。

---

## 2. 已确认的数值基线

### 2.1 原时间步的不收敛不是迭代上限不足

24 fs checkpoint结果：

| 运行 | 步数 | 总迭代 | strict/soft | 代表性E残差 | 代表性J_bkg残差 |
|---|---:|---:|---:|---:|---:|
| 原`dt`, max40 | 5 | 200 | `0/5` | `6.91e-6` | `8.86e-4` |
| 原`dt`, max80 | 5 | 400 | `0/5` | `5.81e-6` | `6.58e-4` |

残差在迭代过程中有界振荡，而不是稳定单调收缩。将上限从40提高到80近似翻倍成本，
但不能恢复严格收敛。因此不得把提高`midpoint-max-iters`作为生产修复。

### 2.2 `dt/2`恢复严格收敛

10步结果：

| 运行 | 总迭代 | 平均迭代 | strict/soft |
|---|---:|---:|---:|
| `dt/2`, max40 | 137 | 13.7 | `10/0` |
| `dt/2`, max80 | 137 | 13.7 | `10/0` |

500步结果：

```text
accepted_steps                 500
strict_accepted_steps          500
soft_accepted_steps            0
total_nonlinear_iterations     6622
mean_iterations_per_step       13.244
final_residual_E               6.98e-7
final_residual_J_bkg           1.75e-6
```

同时满足：

- 背景粒子数相对漂移约`-2.71e-8`；
- 新增能量账误差约占总能量`0.00663%`；
- Beam连续性保持在`1e-11--1e-13`量级；
- 已保存接受态没有负分布；
- 没有NaN/Inf、MPI错配或冻结状态。

所以`dt/2`是当前可信的严格参考和保底生产路径。

### 2.3 “严格收敛”的准确含义

当前生产strict判据是：

```text
residual_E     < 1e-6
residual_J_bkg < 1e-5
```

`residual_f`只做诊断，不参与接受。因此本文中的严格收敛指：

> 场和背景电流满足当前生产中点接受条件。

它不能单独证明完整四维分布函数Picard残差也达到同样精度。后续测试必须继续记录
`residual_f`，但暂时不能增加未经标定的`f`硬阈值。

### 2.4 24--30 fs严格参考段

`output/dual_u_prod_reference/24_to_30_after_beam_fix`显示：

| 项目 | 结果 |
|---|---:|
| `soft_unconverged` | 全部为`0` |
| `coupled_iter_global` | `13--18`，平均约`13.55` |
| 最大E残差 | `9.82e-7` |
| 最大J_bkg残差 | `8.83e-6` |
| FCT active fraction | `0.428--0.438` |
| `x_limiter_min_alpha` | `0` |
| 背景粒子数相对变化 | `-7.40e-7` |
| 30 fs能量账误差/总能量 | `0.3477%` |

因此24--30 fs参考段通过当前数值验收。FCT覆盖率高，但在该时间段没有阻止严格收敛。

30 fs checkpoint已经在服务器生成，本机未复制。无需重新运行24--30 fs来获取checkpoint。

---

## 3. Beam ledger I/O A/B验收

### 3.1 实测结果

同一24 fs checkpoint、相同80 MPI rank和4线程、500个`dt/2`半步：

| 指标 | `full` | `summary` |
|---|---:|---:|
| accepted steps | 500 | 500 |
| strict/soft | `500/0` | `500/0` |
| 总非线性迭代 | 6622 | 6622 |
| 平均迭代/步 | 13.244 | 13.244 |
| 内部wall | 483.155 s | 479.431 s |
| 外部wall | 487.267 s | 483.445 s |
| wall/迭代 | 0.072962 s | 0.072400 s |

性能变化：

```text
内部wall降低 0.77%
外部wall降低 0.78%
speedup约为 1.008
```

`full`生成：

```text
beam_half_step_ledger.dat
501行
125158 bytes
```

`summary`不生成该逐步文件。

### 3.2 物理一致性

以下最终输出逐字节一致：

```text
fields_00002.dat
fields_face_00002.dat
density_00002.dat
current_00002.dat
fv_bkg_e_00002.dat
scalars.dat
step_diagnostics.dat
```

两组`beam_ledger_summary.result`除模式名外完全一致：

```text
accepted_substeps                 500
N_in_total                        1.0409842921512564e20
N_out_total                       0
injected_energy_total             6.5265741900813684e7 J/m2
injected_current_impulse_total   -1.6678442131989001e1
```

background和field hash一致。Beam微观hash不一致，但Beam密度、网格电流、宏观量和ledger
均一致。按当前高效率原则，Beam微观hash差异记为非阻断项，不再单独追查；只有后续出现
Beam密度、电流、能量或连续性差异时才重新打开审计。

### 3.3 验收结论

Beam ledger修改通过：

- `summary`没有改变生产物理解；
- `summary`取消逐步ledger文件和对应长期存储；
- `full`保留短测试审计能力；
- `summary`应作为`diagnostic-level=0/1`的生产默认；
- `full`只用于两半步参考和短程定点审计。

但是该项不是主要性能瓶颈。500步只写约125 KB，Vlasov-Ampere算子成本远大于ledger
I/O，因此不再继续优化Beam ledger。

---

## 4. 24 fs中点加速A/B结果

### 4.1 测试定义

所有测试从同一个24 fs checkpoint开始，并覆盖相同物理时间：

| 组 | 时间步 | 加速 | 步数 |
|---|---|---|---:|
| A | `dt/2` | none | 10 |
| B | 原`dt` | none | 5 |
| C | 原`dt` | Aitken | 5 |
| D | 原`dt` | Anderson depth 3 | 5 |

### 4.2 收敛和性能

| 组 | strict/soft | 平均迭代 | 算子评估 | 内部wall | wall/物理fs |
|---|---:|---:|---:|---:|---:|
| A | `10/0` | 13.7 | 137 | 10.522 s | 1577.2 s/fs |
| B | `0/5` | 40 | 200 | 24.731 s | 3707.1 s/fs |
| C | `0/5` | 40 | 262 | 32.379 s | 4853.5 s/fs |
| D | `0/5` | 40 | 248 | 30.787 s | 4614.9 s/fs |

相对原`dt`无加速B：

```text
Aitken算子评估增加 31%
Aitken内部wall增加约 31%
Anderson算子评估增加 24%
Anderson内部wall增加约 24%
```

两种加速都没有减少Picard迭代，仍然每步做到40次后软接受。

### 4.3 加速候选统计

| 指标 | Aitken | Anderson3 |
|---|---:|---:|
| attempts | 89 | 90 |
| accepted | 27 | 42 |
| fallback evaluations | 62 | 48 |
| residual rejection | 62 | 48 |
| coefficient rejection | 19 | 0 |
| history resets | 81 | 48 |

当前保护逻辑在候选失败后执行额外普通路径算子评估。由于大量候选被拒绝，保护机制本身
显著增加了计算成本。

### 4.4 最终残差

| 组 | `residual_E` | `residual_J_bkg` | `residual_f` |
|---|---:|---:|---:|
| A：`dt/2-none` | `5.15e-7` | `3.21e-6` | `1.164e-3` |
| B：`dt-none` | `5.45e-6` | `5.07e-4` | `3.668e-1` |
| C：Aitken | `1.52e-5` | `1.10e-3` | `3.905e-1` |
| D：Anderson3 | `1.44e-5` | `1.32e-3` | `3.938e-1` |

Aitken和Anderson不仅没有恢复严格收敛，还同时恶化E、J和f残差。

### 4.5 宏观量比较

相对A组`dt/2`严格参考，原`dt`三组在该极短物理区间尚未出现宏观爆炸：

- `Ex_face`相对L2约`7.3e-5--7.5e-5`；
- 背景密度相对L2约`4e-5`；
- 主要背景/总电流相对L2约`0.9%`；
- 个别诊断电流列约`4.7%`；
- Beam ledger累计量在四组覆盖相同物理时间时一致。

但B/C/D全部是软接受状态，因此这些只是相近的未收敛候选，不能作为正确物理解。

加速组相对B组的宏观差异很小：

```text
Aitken vs B:
  Ex_face L2约 3.09e-6
  n_bkg L2约 1.36e-6
  current最大列L2约 8.01e-5

Anderson vs B:
  Ex_face L2约 4.00e-6
  n_bkg L2约 1.64e-6
  current最大列L2约 1.41e-4
```

这不能证明加速正确，只说明三种原`dt`路径停留在彼此接近的软接受状态。

### 4.6 第6节验收结论

当前第一版中点加速明确不通过：

- C/D没有做到`5/5`严格接受；
- 平均迭代数没有低于20；
- 算子评估数没有下降；
- wall time高于无加速B，也远高于`dt/2`参考A；
- E/J/f最终残差均恶化；
- 当前参数不能进入30 fs压力测试；
- 当前参数不能进入原`dt`生产。

---

## 5. 当前性能瓶颈的重新定位

### 5.1 不是Beam ledger

ledger切换只能节省约`0.8%`，因此长期生产慢的主因不是逐步Beam文件写出。

### 5.2 主要成本是非线性算子评估

每次完整评估包含：

- 四维背景分布的x方向高阶通量和FCT；
- 速度空间u方向推进和FCT；
- 背景矩、连续性电流和dual-u功电流；
- Yee face Ampere更新；
- halo交换和MPI collective。

`dt/2`每步约需13次完整评估；原`dt`每步做满40次仍不收敛。当前Aitken/Anderson又为
被拒候选增加额外评估，所以其性能必然更差。

### 5.3 当前固定点映射过刚且不够光滑

主要根因链为：

```text
时间步增大
  -> E/f/J耦合刚性增强
  -> FCT active set随迭代切换
  -> 固定点残差有界振荡而非单调收缩
  -> 场端点外推候选频繁被拒
  -> 回退引入额外完整算子评估
  -> 仍达到40次并软接受
```

第一版只加速场端点，背景`guess=work`仍是未松弛Picard更新。实测说明，仅对场做标量
Aitken或低深度Anderson不足以控制原`dt`下的完整E-f-J耦合映射。

不得因此直接对完整四维`f`做Anderson线性组合。这样可能制造负分布、改变FCT active set，
并显著增加内存流量。

---

## 6. 当前执行决策

### 6.1 立即保留

生产和后续测试继续保留：

- `background-coupling-mode=dual-u`；
- 当前x/u高阶通量和FCT；
- Beam轨迹沉积和开放边界；
- Yee face Ampere；
- 原E/J严格阈值；
- `beam-ledger-mode=summary`；
- `diagnostic-level=1`；
- 当前已验证的算子内核性能优化。

### 6.2 立即冻结

在新方案通过24 fs测试之前：

- 不执行30 fs Aitken/Anderson压力测试；
- 不执行原`dt`加速版30--36 fs；
- 不将当前Aitken或Anderson用于生产；
- 不继续只调`accept-ratio`、history depth或最大系数；
- 不提高`midpoint-max-iters`到80；
- 不放宽E/J接受阈值；
- 不把软接受重新作为长期生产常态。

### 6.3 当前可信生产方式

当前唯一已经验证的生产方式是：

```text
checkpoint保存dt为原dt时:
  --dt-scale 0.5
  --midpoint-acceleration none

checkpoint保存dt为dt/2时:
  --dt-scale 1.0
  --midpoint-acceleration none
```

如果现有24--36 fs `dt/2`参考作业仍在运行，可以继续完成，不要重复提交同一区间。

---

## 7. 下一步：寻找最大严格稳定时间步

### 7.1 为什么先做时间步包络测试

目前只知道：

```text
0.5 * dt_original：稳定，平均约13--14次
1.0 * dt_original：40次仍软接受
```

中间可能存在一个明显优于`dt/2`、同时仍严格收敛的时间步。相比继续修改物理算子或对
四维分布做非线性外推，寻找稳定时间步包络风险更低，且不会改变离散方程。

第一轮测试：

```text
0.60 * dt_original
0.70 * dt_original
0.80 * dt_original
0.90 * dt_original
```

每档先跑10步。判据：

1. `10/10`严格接受；
2. 无NaN/Inf和hard failure；
3. 平均迭代数目标不超过25；
4. `residual_f`不相对`dt/2`出现数量级恶化；
5. 按`wall_seconds_per_physical_fs`比较真实效率；
6. Beam连续性、背景质量和能量不得出现新退化。

若相邻两档一档通过、一档失败，再在两者之间二分一次。不要一次创建大量参数扫描。

### 7.2 编译命令

只要源码没有继续修改，不需要重新编译。需要重编译时使用：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

rm -rf build
mkdir build
cd build

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake .. -DCMAKE_CXX_COMPILER=mpicxx
make -j4 fp_solver
```

### 7.3 24 fs，`0.60 dt`测试

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J dt060
#SBATCH -o dt060_%j.out
#SBATCH -e dt060_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
CHECKPOINT_24="/ABSOLUTE/PATH/TO/24FS_CHECKPOINT"
OUT="./output/dt_envelope_24fs/dt_060"
mkdir -p "$OUT"

TIMEFORMAT=$'real_s=%R\nuser_s=%U\nsys_s=%S'
{ time yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.60 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 10 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir "$OUT" \
  --overwrite-output \
  > "$OUT/run.out" 2> "$OUT/run.err"; } 2> "$OUT/timing.txt"
```

### 7.4 24 fs，`0.70 dt`测试

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J dt070
#SBATCH -o dt070_%j.out
#SBATCH -e dt070_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
CHECKPOINT_24="/ABSOLUTE/PATH/TO/24FS_CHECKPOINT"
OUT="./output/dt_envelope_24fs/dt_070"
mkdir -p "$OUT"

TIMEFORMAT=$'real_s=%R\nuser_s=%U\nsys_s=%S'
{ time yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.70 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 10 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir "$OUT" \
  --overwrite-output \
  > "$OUT/run.out" 2> "$OUT/run.err"; } 2> "$OUT/timing.txt"
```

### 7.5 24 fs，`0.80 dt`测试

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J dt080
#SBATCH -o dt080_%j.out
#SBATCH -e dt080_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
CHECKPOINT_24="/ABSOLUTE/PATH/TO/24FS_CHECKPOINT"
OUT="./output/dt_envelope_24fs/dt_080"
mkdir -p "$OUT"

TIMEFORMAT=$'real_s=%R\nuser_s=%U\nsys_s=%S'
{ time yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.80 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 10 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir "$OUT" \
  --overwrite-output \
  > "$OUT/run.out" 2> "$OUT/run.err"; } 2> "$OUT/timing.txt"
```

### 7.6 24 fs，`0.90 dt`测试

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J dt090
#SBATCH -o dt090_%j.out
#SBATCH -e dt090_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
CHECKPOINT_24="/ABSOLUTE/PATH/TO/24FS_CHECKPOINT"
OUT="./output/dt_envelope_24fs/dt_090"
mkdir -p "$OUT"

TIMEFORMAT=$'real_s=%R\nuser_s=%U\nsys_s=%S'
{ time yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.90 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 10 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir "$OUT" \
  --overwrite-output \
  > "$OUT/run.out" 2> "$OUT/run.err"; } 2> "$OUT/timing.txt"
```

四档必须分别提交为独立作业，不能在同一个作业内串行运行后直接比较wall time。

### 7.7 选择规则

先按以下顺序筛选：

1. 删除任何出现软接受的档位；
2. 删除平均迭代数大于25的档位；
3. 在剩余档位中比较`wall_seconds_per_physical_fs`；
4. 选择效率最高者，不一定选择数值最大的`dt_scale`；
5. 对选中档位再做500步稳定性测试；
6. 500步全部严格接受后，才允许到30 fs checkpoint复验。

如果`0.60--0.90`全部失败，则追加`0.55`测试；如果`0.90`通过，则不自动采用原`dt`，
先测试`0.95`并仍要求全部严格接受。

---

## 8. 第一版Aitken/Anderson的后续处理

### 8.1 当前状态

源码中的：

```text
--midpoint-acceleration aitken
--midpoint-acceleration anderson
```

可以保留用于研究和短测试，但生产必须显式使用：

```text
--midpoint-acceleration none
```

不需要立即删除实现，因为`none`路径已经用于可信参考；但不得让默认值变成Aitken或Anderson。

### 8.2 如果以后重新修复加速

必须先解决“拒绝候选导致额外完整算子评估”的结构性问题：

1. 每个循环最多执行一次生产算子；
2. 被拒候选不得在同一循环立即再执行一次完整fallback算子；
3. 拒绝后恢复保存的普通场候选和背景输入，在下一循环自然评估；
4. `total_operator_evaluations`必须与真实评估次数一致；
5. 短程诊断模式可同时评估plain/trial来标定，但生产模式不能永久双评估；
6. 候选必须相对“普通路径预计收缩”获得净收益，而不只是略小于上一残差；
7. 修复后重新从24 fs执行A/B/C，不得直接从30 fs开始；
8. 若场端点加速仍失败，不升级为完整`f`的Anderson组合。

只有新的24 fs测试达到以下要求，才重新开放30 fs测试：

```text
strict/soft = 5/0
mean iterations < 20
operator evaluations < 原dt无加速
wall/physical-fs < dt/2参考
E/J/f残差不退化
```

---

## 9. 30 fs checkpoint验收

30 fs checkpoint已在服务器，不需要传回本机。使用保存时相同的80 MPI rank重启1步：

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J ckpt30_check
#SBATCH -o ckpt30_check_%j.out
#SBATCH -e ckpt30_check_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
CHECKPOINT_30="/ABSOLUTE/PATH/TO/30FS_CHECKPOINT"
OUT="./output/checkpoint30_restart_check"
mkdir -p "$OUT"

yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --dt-scale 1.0 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir "$OUT" \
  --overwrite-output \
  > "$OUT/run.out" 2> "$OUT/run.err"
```

必须满足：

1. manifest、MPI size和网格元数据兼容；
2. checkpoint保存态与读回态hash一致；
3. 接受步数为1，物理时间真实增加一个checkpoint保存的`dt/2`；
4. 没有NaN/Inf、MPI错误或冻结状态。

这项restart验收可以与第7节时间步测试独立进行，不需要等待加速修复。

---

## 10. 30 fs后期复验

只有第7节选出的时间步通过24 fs的500步测试后才执行。

30 fs checkpoint保存的是`dt/2`。若24 fs选出的比例是：

```text
dt_test = s * dt_original
```

则从30 fs checkpoint恢复时应传：

```text
--dt-scale 2*s
```

例如24 fs选出`0.70*dt_original`，30 fs复验使用：

```text
--dt-scale 1.40
```

30 fs先运行10步，要求：

- `10/10`严格接受；
- 平均迭代不超过24 fs同档结果的合理波动范围；
- 不出现新的质量、能量、Beam连续性或`residual_f`退化；
- `wall_seconds_per_physical_fs`仍优于保存的`dt/2`参考。

若30 fs失败，生产继续使用checkpoint保存的`dt/2`，不得只凭24 fs测试提高时间步。

---

## 11. 后续生产策略

### 11.1 当前可执行路径

若现有`dt/2`参考作业已经生成36 fs checkpoint，优先沿该checkpoint继续，不重复24--36 fs。

分段建议：

```text
36 -> 45 fs
45 -> 60 fs
60 -> 90 fs
90 -> 120 fs
```

每段：

- 保存段末checkpoint；
- 长段增加中间checkpoint；
- 使用`diagnostic-level=1`；
- 使用`beam-ledger-mode=summary`；
- 使用`midpoint-acceleration=none`；
- 记录strict/soft、平均/最大迭代和wall/physical-fs；
- 检查核心Ex包络、相位、n_bkg、J、能量和质量；
- 边界尖峰只做隔离监测，不用单点幅值否定核心结果；
- 一旦出现持续软接受或迭代数持续增长，停止扩展下一段。

### 11.2 时间步包络通过后的路径

只有24 fs和30 fs都验证通过，才能用选出的时间步从最近checkpoint开始新的生产分支。
新分支必须与`dt/2`参考在相同物理时刻比较：

- 核心区Ex相对L2、峰值、主峰位置和相位；
- n_bkg和J_bkg；
- 能量账和背景质量；
- Beam连续性和累计注入/出流；
- 接受态负分布诊断。

不能只因为计算更快就替换严格参考。

---

## 12. 总体验收条件

任何新时间步或新加速方案进入生产必须同时满足：

1. 不改变dual-u、FCT、Beam和Yee Ampere离散；
2. 不放宽E/J阈值；
3. 短测试全部严格接受；
4. 500步稳定性测试全部严格接受；
5. 24 fs和30 fs两个状态点均通过；
6. 每物理fs的wall time优于`dt/2`；
7. `residual_f`不出现数量级退化；
8. 无NaN/Inf、MPI不一致或冻结状态；
9. 背景质量、能量和Beam连续性不退化；
10. 核心波形与`dt/2`时间收敛参考相容。

必须停止的情况：

- 依靠软接受维持推进；
- 只降低日志残差而没有减少真实算子评估；
- 候选回退使wall time增加；
- 通过关闭物理算子、放宽阈值或减少真实推进换取速度；
- 24 fs通过但30 fs失败；
- 宏观波形、质量或能量出现系统性退化。

---

## 13. 当前执行清单

```text
[x] dual-u离散闭合与MPI回归
[x] 生产算子基础性能优化
[x] 24 fs原dt的40/80次判别
[x] 24 fs dt/2的10步严格收敛
[x] 24 fs dt/2的500步稳定性
[x] 24--30 fs严格参考段
[x] 服务器生成30 fs checkpoint
[ ] 服务器端验收30 fs checkpoint restart/hash
[x] Beam ledger summary/full实现
[x] Beam ledger 500步A/B：物理一致
[x] Beam ledger性能结论：仅约0.8%，不再深挖
[x] 第一版Aitken实现
[x] 第一版Anderson3实现
[x] 24 fs加速A/B/C/D
[x] 第一版加速判定失败并冻结生产使用
[ ] 24 fs时间步包络：0.60/0.70/0.80/0.90
[ ] 选择最优严格档并运行500步
[ ] 在30 fs checkpoint复验最优档
[ ] 通过后决定是否建立更大时间步生产分支
[ ] 分段扩展到45/60/90/120 fs
```

核心原则：

> 当前首要目标不是让原`dt`在日志上看起来收敛，而是在不改变Vlasov-Ampere离散和
> 不依赖软接受的条件下，最小化每单位物理时间的真实算子评估成本。`dt/2 + none`
> 是可信基线；ledger优化已经完成但不是主瓶颈；第一版Aitken/Anderson实测失败，
> 在重新设计并通过24 fs验收前不得进入后期测试或生产。
