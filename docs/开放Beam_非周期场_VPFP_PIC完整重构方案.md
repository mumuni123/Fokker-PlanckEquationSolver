# 开放 Beam、非周期场、Eulerian bulk + PIC 高能尾部 VPFP-PIC 完整重构方案

## 0. 文档地位与执行原则

### 0.1 文档用途

本文档是当前仓库后续重构的主实施规格，面向直接执行代码修改、测试和验收的自动编码模型。执行者应把本文档视为接口契约和阶段门，而不是建议列表。

本文档基于当前代码状态重新编写：

- 阶段 1 至阶段 4 已经建立开放空间边界、非周期静电场、背景 Eulerian Vlasov 主体、Beam PIC 和无碰撞 Strang 推进。这些内容继续保留。
- 当前阶段 5 通过不断扩大 $u_{\parallel,\max}$ 和 $N_{u_\parallel}$ 容纳高能尾部的路线，不再作为最终生产路线。
- `docs/阶段5无碰撞短生产_问题修复实施方案.md` 中的“有限解析区 + guard-tail”方案仅保留历史参考，不再执行。
- 新阶段 5 起采用同一背景电子物种的混合表示：低能主体由 Eulerian Vlasov/FP 网格表示，高能尾部由独立 PIC 粒子表示。
- tail PIC 设计已对照本机 EPOCH `v4.20.1-0-gf294c484-clean` 的 1D pusher、沉积、边界、MPI 分箱和碰撞源码完成审计；审计结论集中在第 3.4、6、10 节。

旧周期 Vlasov-Ampere 求解器已经另库存档。本仓库不要求保留旧求解器、旧 checkpoint 兼容层或双生产路径。

### 0.2 状态标记

- **[已完成-保留]**：当前实现方向正确，禁止无关重写。
- **[已完成-需接口扩展]**：现有算法保留，只增加混合尾部所需接口。
- **[待实现]**：当前没有生产实现。
- **[待替换]**：已有占位实现，但不足以求解目标物理模型。
- **[验收后启用]**：实现后必须通过本节指定测试才能接入下一阶段。
- **[条件阶段]**：只有诊断满足触发条件时才实施。

### 0.3 不可违反的原则

1. `bulk` 与 `tail` 是同一背景电子物种的两种数值表示，不能被当作两个独立电子物种。
2. bulk-to-tail 或 tail-to-bulk 转换是内部表示变换，不能产生或删除物理电子，也不能成为外部能量源。
3. 场始终由非周期 Gauss/Poisson 约束求解，不恢复旧 Ampere 推进，不做周期零模扣除。
4. Beam 继续采用开放注入、开放流出，不周期回卷。
5. 不允许通过全局能量补丁、人工缩放电流、强制 $\langle E\rangle=0$、删除高能尾部或无记账裁剪得到“漂亮结果”。
6. 不以逐点贴合 EPOCH 为验收标准。验收对象是守恒关系、收敛阶、宏观包络、能谱、能量分配和参数收敛性。
7. 所有试探推进必须是事务式的。失败状态不得污染已接受状态、粒子 ID、随机数流、累计账本或 checkpoint。
8. 不允许一次性大改后直接长跑。必须按阶段门逐层接入。
9. 本文对 EPOCH 的引用只服务于 background-tail PIC 子系统。EPOCH 不是本项目 Poisson、Eulerian bulk、bulk--tail 转换或混合碰撞反作用的模板。
10. “参考 EPOCH”表示复核其离散思想和源代码行为，不表示逐行移植 GPL Fortran 源码。生产实现必须使用本仓库自己的接口、数据所有权、测试和许可证边界。

### 0.4 最近修改记录：不等权 tail--tail 人口爆炸根因修复

**状态：代码侧完成，并已通过本地回归和 H9 多 rank Beam 12 fs 人口增长阻断门。**

修改日期：2026-08-06。

问题证据：最新 H9 Beam 12 fs 在近轴 SDE 修复后跨过旧 step 413 失败点，
但 `conversion_N=0` 时 tail 宏粒子仍从 10 增至 380272350；单 rank 最大
295290241，collision wall time 增至约 282 s/步。该现象是旧不等权
`virtual-split` 永久实体化 residual weight 导致的表示人口近指数增长，
不是 MPI 死锁。

生产代码修改：

- `src/background_tail_nanbu_perez.cpp`：用有界 Sentoku--Kemp 不等权修正
  替代永久 residual 粒子生成；每个 pair 固定输出两个原 ID/weight；
- `src/background_tail_collision.h`：默认宏粒子增长预算改为 0，并明确
  `virtual-split` 只是兼容 token；
- `src/vpfp_integrator.cpp`：CLI 的增长预算 0 不再重解释为 1.0；
- `src/vpfp_checkpoint.cpp`、`src/vpfp_diagnostics.cpp`、
  `src/vpfp_diagnostics.h`、`src/main_vpfp.cpp`：checkpoint/snapshot manifest
  和启动日志增加 `tail_collision_weight_algorithm=sentoku-kemp-bounded-v1`
  与实际 `tail_collision_max_particle_growth`；
- `tests/background_tail_collision_weight_test.cpp`：增加零增长和 1024-pair
  集合统计守恒测试；
- `tests/checkpoint_roundtrip_test.cpp`：验收 manifest 算法版本。

当前数值语义：等权 pair 逐事件严格保持加权三动量和相对论能量；不等权
Sentoku--Kemp pair 逐事件保持粒子数、weight 和相对论能量，动量按多
cell/多 seed 集合统计闭合。不得再把不等权模式描述为单事件严格动量守恒。

本地验收：单 pair 始终 2 粒子且能量误差为 0；1024-pair 集合动量相对
误差 $4.36\times10^{-4}$、能量误差 $7.82\times10^{-15}$；近轴 SDE、
双温热化、各向同性化、bulk--tail 反作用、Beam 四通道生产路径、checkpoint
往返和完整生产目标编译均通过。

失效数据：旧 `hybrid_h9_beam_12fs` 目录和任何已保存数亿 tail 粒子的
checkpoint 不得续跑。新的 H8/H9 结果统一写入带 `sk_bounded_v1` 的独立
目录；具体命令和执行顺序见 §17.9--17.10。

### 0.5 历史 H9 多 rank 结论：碰撞人口爆炸已修复，旧静态阈值接口未通过

**状态：H9 部分通过；禁止进入 120 fs 生产。**

本节记录2026-08-06旧 static-cell conversion 的历史结果，已被 §0.9 的
flux-interface 最新生产回归更新；保留本节只用于解释为何停止旧静态转换路线。

2026-08-06 使用 80 rank、`DIRICHLET_PHI(0,0)`、
`sentoku-kemp-bounded-v1` 和零碰撞人口增长预算完成了新的 H9 回归。结论必须
分层表述，不能再把 H9 写成“尚未重跑”，也不能因 40 fs 作业跑完而判定整体
通过。

已经通过的部分：

- 不等权碰撞单元/集合测试通过：单 pair 固定输出 2 粒子且能量逐事件闭合；
  1024-pair 集合动量相对误差 $4.36\times10^{-4}$、能量误差
  $7.82\times10^{-15}$；
- checkpoint 往返、Beam 四碰撞通道和无 Beam 回归通过；
- Beam 12 fs 阻断档跨过旧 step 413，首次 tail 约在 10.46 fs 出现；最终
  tail 宏粒子数 2880，且 `conversion_N=0` 时没有碰撞诱发的人口增长；
- no-Beam 40 fs 共 1564 步全部接受，背景动能基本不变，Gauss 残差约
  $1.16\times10^{-21}$；
- Beam 40 fs controller-off 共 1564 步全部接受，旧版数亿宏粒子的指数增长
  已消失；tail--bulk 和 bulk 反作用的全局累计账本在当前求和/采样误差内闭合。

尚未通过的部分：

1. **bulk--tail 阈值接口存在显著人工谱谷。** controller-off 的 combined
   $dN/dK$ 在 5.8--6.0、6.0--6.2、6.2--6.4 MeV 三个 bin 约为
   $1.04\times10^{19}$、$3.82\times10^{16}$、$3.00\times10^{18}$，即跨
   $K_{out}=6$ MeV 先下降约 271 倍、再回升约 78.5 倍。这不满足 §19.2 的
   阈值连续性要求。
2. **controller-off 资源预算仍偏高。** 40 fs 时全局 tail 宏粒子数
   10466420，最大 rank 1678802；最后 100 步平均 wall time 约 2.21 s，其中
   collision 约 1.72 s。人口不再指数爆炸，但负载不平衡和绝对粒子数仍不适合
   直接外推到 120 fs。
3. **现有人口控制器不能进入生产。** controller-on 虽把 40 fs 粒子数降至
   1365113，但相对 controller-off 已使场能改变约 3.47%、combined 背景动能
   改变约 4.96%、tail 数和能量改变约 24%--26%；全程平均 wall time 反而约
   增加 2.2 倍，控制步平均约 29 s、最慢约 189 s。它同时未消除阈值谱谷。
4. 当前输出中没有保存 MPI collision failure-consensus 的正式结果，因此该门
   仍是“缺少证据”，不能写成已通过。

当前优先级固定为：先按 §7.11 独立定位并修复 bulk-to-tail 阈值接口；生产
基线保持 population controller 关闭；不得通过降低阈值、扩大谱 bin、平滑
输出或强制矩补丁掩盖谱谷。Sentoku--Kemp 有界碰撞修复已经通过，不应在本轮
阈值修复中改动。

### 0.6 2026-08-08 阈值接口独立测试与 12 fs 回归复核

**最新状态：转换器基础守恒通过，阈值感知粗分组通过；真实近轴阈值接口仍未
闭环，H9 继续判定为 partial pass。**

本轮复核覆盖 `output/h9_threshold_interface`、
`output/hybrid_h9_threshold_fixed_beam_12fs` 和
`output/no_dir_multibin.result`，结论如下。

已经确认的正确部分：

- 单 cell 转换六矩残差为 0；Poisson 转换前后密度相对 $L_2$ 误差约
  $3.72\times10^{-13}$、场相对 $L_2$ 误差约 $4.16\times10^{-13}$、能量误差
  约 $6.67\times10^{-15}$；事务测试全部通过；
- multibin 的 1/2/5 rank 回归均通过，最大转换矩残差约
  $5.5\times10^{-15}$--$7.4\times10^{-15}$，rank 一致性约
  $3\times10^{-12}$；
- 黄金 quartet 不压缩路径对测试谱的误差约 $10^{-14}$；当前旧压缩路径会产生
  $0.22$--$1.27$ 的归一化谱 $L_1$ 误差；阈值感知分组在显式粗能量组上将
  $N/K$ 和 edge spectrum 保持到约 $10^{-14}$；
- 新 12 fs Beam 回归 470/470 步接受，无 split；首次 tail 在 step 409、
  $t\approx10.4643$ fs 出现，末态 tail 宏粒子数 2892，最大 rank 1600；
  `unexplained_growth_steps=0`，Gauss 电荷残差最大约
  $3.11\times10^{-15}$，转换 $N/K$ 残差分别不超过
  $1.42\times10^{-15}$ 和 $1.30\times10^{-18}$。这说明有界不等权碰撞修复
  没有回归。

尚未闭环且决定下一步方向的证据：

1. 当前支撑审计把全部 $u_\perp$ 环混合后报告 6.0--6.2 MeV 内有 22 个中心、
   最大空隙约 0.0533 MeV，但这不能代表 Beam 驱动的近轴高能尾部。按实际生产
   网格，在最低 $u_\perp\approx0.002529$ 环上，相邻正向 $u_\parallel$ 中心的
   能量约为 5.809279、6.228272、6.675382 MeV；真实活跃路径在
   6.0--6.2 MeV 内没有 cell center。现有全局支撑审计因此会给出误导性的
   “有支撑”结论。
2. 12 fs 生产结果从首次转换开始就存在硬空 bin。10.515 fs 时
   5.8--6.0、6.0--6.2、6.2--6.4 MeV 的 combined 数量约为
   $1.2892\times10^{10}$、0、$5.6126\times10^7$；12 fs 时约为
   $2.48799\times10^{11}$、0、$4.71087\times10^{10}$。该空洞不是长期碰撞
   或人口控制器造成，因为 population controller 关闭，且它在首次转换快照中
   已经存在。
3. 阈值感知策略只保证显式粗组不跨组搬运质量。使用更细直方图时，其谱
   $L_1$ 误差仍约为：0.05 MeV bin 下 $0.63$--$0.84$，0.1 MeV 下
   $0.48$--$0.54$，0.2 MeV 下 $0.19$--$0.335$。因此“edge group 守恒通过”
   不能写成“局部谱保真通过”。
4. 当前 MPI 测试的 `status=PASS` 判据不完整：
   `tests/bulk_tail_threshold_interface_mpi_test.cpp` 只对 golden 路强制检查
   `spectrum_L1_rel<=1e-10`，threshold-aware 路即使谱误差很大也可通过；现有
   输出也只有不完整/重复追加的 MPI 结果，不能替代独立的 1/2/5 rank 验收。
5. 当时的 `conversion_source_bins 0` 说明旧版转换源谱诊断没有形成全局有效
   数据。该缺陷现已修复：最新 12 fs 回归已写出带显式边界的全局四阶段源谱，
   并证明空 bin 在进入转换器前已经存在。该旧结果只保留为诊断演进记录。

当前 `BulkTailConverter` 已将默认策略提前切换为
`THRESHOLD_AWARE_COMPRESSION`。鉴于上述门槛尚未完成，该实现只能标记为
**生产候选**，不能标记为已验收生产默认。后续不得先跑 40 fs/120 fs，也不得
继续扩大压缩 support 掩盖近轴离散支撑问题。已经完成的定位过程见
§7.11.2--§7.11.12；当前实施入口为 §7.11.13--§7.11.15 与 §17.10.1。

### 0.7 2026-08-08 全局源谱阻断回归：根因已定位到转换前离散表示

**状态：12 fs 阻断回归通过运行、守恒和诊断门；阈值物理接口仍未通过。**

最新目录 `output/hybrid_h9_threshold_source_beam_12fs` 完成 470/470 步，未用
split。最终 tail 宏粒子数 2764、最大 rank 1600；Gauss 电荷残差最大约
$3.55\times10^{-15}$，转换数目和能量相对残差最大分别约
$7.89\times10^{-16}$、$1.08\times10^{-18}$。碰撞人口、Gauss、转换总矩和
事务语义均无回归。

全局 `conversion_source_accepted_steps.dat` 已正确写出 30080 条事件-bin 数据，
不再出现 `conversion_source_bins 0`。首次非零转换位于 step 409、
$t=10.4643007$ fs，唯一非零源 bin 为约 6.215085--6.270048 MeV，并严格满足：

$$
N_{\rm pre\ bulk}=N_{\rm removed\ bulk}=N_{\rm created\ tail}
=9.1922729\times10^6.
$$

截至 12 fs 共 54 个非零转换事件，**全部只来自上述同一个源 bin**。累计创建
物理粒子数为 $4.8121025\times10^{10}$，最终接受态
`accepted_tail_total_N` 与其一致；逐 bin 的 `pre-removed` 最大绝对误差为 0，
`removed-created` 最大绝对误差约 $1.67\times10^{-6}$，相对累计量约
$3.5\times10^{-17}$。因此转换器没有删除阈值质量，后续推进、MPI 迁移和碰撞
也没有丢失 tail 数目。

80 rank 快照聚合给出：

| 时间 | 5.8--6.0 MeV | 6.0--6.2 MeV | 6.2--6.4 MeV |
|---|---:|---:|---:|
| 10.413 fs，首次转换前 | $1.3331\times10^{10}$ | 0 | $2.7599\times10^7$ |
| 10.515 fs | $1.2892\times10^{10}$ | 0 | $5.6126\times10^7$ |
| 10.618 fs | $1.1046\times10^{10}$ | 0 | $8.4024\times10^7$ |
| 12 fs | $2.3824\times10^{11}$ | 0 | $4.6294\times10^{10}$ |

最强证据是 10.413 fs 尚未转换时空 bin 已存在。结合近轴 `k=0..7` 的相邻中心
约为 5.809 和 6.228 MeV，可确定：当前谱谷首先来自把有限体积单元的全部质量
按 cell-center 能量当作点质量，并按 cell-center 判定 conversion。旧均匀压缩、
threshold-aware 压缩、碰撞核和 population controller 都不是该硬空洞的根因。

这里不能立即得出“必须扩大 `Nv/Nmu`”或“必须静态切掉跨阈值 cell”的结论。
有限体积 cell 表示一段 $(u_\parallel,u_\perp)$ 区域；应先建立单元体积积分的
能谱参考和单元内守恒 PIC 加载，判断连续单元表示是否已经覆盖 6.0--6.2 MeV。
静态 cut-cell 若每步对剩余 cell mass 重复按比例抽取，会产生非物理持续泄漏，
明确禁止。已完成结论见 §7.11.8--§7.11.12，下一步详见
§7.11.13--§7.11.15。

### 0.8 前一阶段决策：停止局部固定速度网格路线

§7.11.13--§7.11.16 已完成。真实80-rank、12 fs审计与固定网格离线设计共同给出：

- 阈值以下转换质量只有 $1.11\times10^{-16}$，转换事务和阈值判定不是主因；
- 当前网格的 $N/J_x/K/\Pi_{xx}$ 表示差较小，但 $\Pi_\perp$ 累计相对差为
  $3.26\times10^{-2}$，单cell最大约11.3%；
- 真实转换cell的中心六矩对非负subcell支撑不可行；
- 16B扫描中没有GREEN候选。最佳有效候选仍有
  $R_{L1}^{\Pi_\perp}=1.43\times10^{-2}$；误差更小的候选违反预先声明的最小网格宽度
  约束，不能用于生产。

因此局部固定速度网格路线已经停止。§7.11.17 的“最终有限体积面通量 -> 正权扫掠
质量包 -> tail PIC粒子”通量式 representation conversion 已完成 17A--17F 的全部独立
与阶段验收；其后的H9生产级回归也已在§0.9完成。以下纪律仍适用于生产路径：

1. 当前subcell PIC加载继续关闭；不得扩大16B参数扫描，也不得增加全局速度cell数；
2. 不得只修改converter后继续静态清空高能cell；bulk sink和PIC source必须来自同一份
   最终受限面通量；
3. 17A--17F 均已通过；后续改动不得破坏其接口、只读、事务、MPI、checkpoint 或碰撞面通量门；
4. 17B必须是只读审计，且与static-cell基线状态逐bit一致；
5. 17C/17D通过后才能进入17E无碰撞生产A/B；
6. `chang-cooper-flux + exporting-absorbing` 已提供碰撞最终保守面通量；H9生产推进门
   已通过，120 fs仍须等待最新flux-interface阈值能谱和资源外推完成；
7. 不允许用全局能量补丁、负权粒子、速度事后缩放或放宽残差门掩盖转换误差。

### 0.9 2026-08-11 最新状态：flux-interface H9 40 fs 生产门通过

**状态：H9 controller-off 生产推进门通过；120 fs 最终生产仍需阈值能谱和资源外推验收。**

最新 `output/h9_flux_interface_17f_beam_25p5fs_tail_owned_fix` 已完成25.5 fs：

- 最后 accepted step 为997，80个rank最终均为 `step_accepted`；
- 无 rejected step、failure文件或转换控制流分叉；
- `support_limit_violation_count=0`、`duplicate_id_count=0`、
  `face_ledger_mismatch_count=0`；
- 六矩转换残差约为 $10^{-15}$，tail-owned余额残差约为 $5.19\times10^{-24}$。

随后 `output/h9_flux_interface_17f_beam_40fs` 完成生产级阻断回归，
`h9_production.result` 为 `status=PASS`。主要结果为：

| 指标 | 结果 | 既定门槛 |
|---|---:|---:|
| accepted/rejected steps | 1564 / 0 | 1564 / 0 |
| 最终时间 | 40 fs | 40 fs |
| 最大Gauss电荷归一残差 | $1.14\times10^{-13}$ | validator通过 |
| 最大tail数量账误差 | $3.75\times10^{-14}$ | validator通过 |
| 最大conversion残差 | $3.62\times10^{-14}$ | validator通过 |
| 最终/最大全局tail宏粒子 | 7,829,024 | $<12,000,000$ |
| 最大单rank tail宏粒子 | 900,184 | $<2,000,000$ |
| 后100步最大wall time | 2.446 s/step | $<5$ s/step |

80个rank最后均处于 step 1564 的 `step_accepted`；所有validator gate均为1。
这证明tail-owned账本、秩亏压缩、转换失败MPI共识以及40 fs controller-off生产路径
已经闭环。它尚未证明：40 fs阈值附近combined $dN/dK$ 已满足 §19.2，也不能仅按
40 fs粒子数线性外推就宣称120 fs内存可接受。下一步只做两项：

1. 使用25 fs和40 fs快照生成cell-volume bulk + PIC tail的统一能谱，执行6 MeV阈值
   $pm0.4$ MeV连续性与分箱收敛验收；
2. 根据10.5--40 fs真实tail增长率、每rank最大值和wall time给出120 fs资源上界。

两项通过后才进入 §17.11；不得重新启用旧static-cell conversion，也不得为降低粒子数
直接开启尚未通过物理A/B的population controller。

### 0.10 2026-08-12 最新权威状态：controller-off 已连续推进到 120 fs

**状态：Gate C 连续推进与守恒链通过；H11 最终物理验收尚未完成。**

`output/hybrid_flux_interface_controller_off_40_to_120fs` 已从约 40 fs 连续推进到
120 fs，共 3127 个接受步。全部时间步满足 `accepted=1`，无 split、collision rollback、
NaN/Inf 或速度边界损失。主要离散闭合结果为：

- Gauss 电荷残差中位数约 $1.42\times10^{-14}$，最大约 $1.56\times10^{-13}$；
- Tail 数目平衡残差中位数约 $7.78\times10^{-15}$；
- conversion 的 $N/P_x/K$ 残差通常约 $10^{-14}$；
- 40--120 fs 累计 `conversion_N` 与 Tail 物理数增量的相对差约
  $3.41\times10^{-14}$；
- 6 MeV 附近 combined 能谱没有重新出现 static-cell 路线的人工硬谷。

因此，flux-interface 转换、事务接受、碰撞推进和长期 MPI 控制流已经通过 Gate C 的
**数值连续推进门**。但真实长跑推翻了 40 fs 时关于 H10 和资源增长的两个前瞻判断：

1. Tail 宏粒子由约 $7.85\times10^6$ 增长到 $1.4317\times10^8$；最大 rank 由
   $8.95\times10^5$ 增长到 $5.30\times10^6$。实际全局粒子数约为原
   $4.83\times10^7$ 外推的 2.97 倍。
2. `tail_return_mode=none` 使碰撞减速后的低能 Tail 永久留在 PIC 表示中。低于
   6 MeV 的 Tail 数目占比从 50 fs 的约 4.96% 增至 80 fs 的 50.98%，
   120 fs 完整checkpoint审计值为58.80%；Tail 平均动能由约 6.91 MeV 降至 5.17 MeV。
3. 平均单步 wall time 从 40--50 fs 的 2.54 s 增至 110--120 fs 的 11.65 s；
   40--120 fs 实测约 5.99 h。增长主要来自 collision、conversion 和 Tail push，
   不是文件 I/O 或 MPI collective。

Tail 在 120 fs 的物理数分数仍仅约 $1.95\times10^{-5}$，尚未主导背景密度；问题首先是
**表示归属和计算复杂度失控**，而不是 bulk--tail 数目不守恒。原 Gate B 在 40 fs
checkpoint 上的 `h10_required=0` 只能保留为当时的准入决定，不能继续作为 120 fs 的最终
H10 结论。后续H10 Tail-to-bulk返回现已完成R0--R4验证，最终状态见§0.11；未经验证的
Population Controller 仍不得代替返回机制。

能量方面，40--120 fs 域内 $U_E+K_{\rm bulk+tail}+K_{\rm beam}$ 增加约
$9.30\times10^8\ \mathrm{J/m^2}$，同期逐步累计 `collision_reservoir` 约为
$-8.93\times10^8\ \mathrm{J/m^2}$。计入后仍余约
$3.70\times10^7\ \mathrm{J/m^2}$，约为 40 fs 域内能量的 0.79%。现有输出没有完整的
背景 reservoir 开放边界能量通量累计，因此该余额不能直接判为数值生能，也不能判为
完全闭合。下一步必须先补齐只读能量账，再决定是否修改物理推进。

本节是全文当前状态的权威入口。§0.5--§0.9 中“禁止/尚未进入 120 fs”的句子均为历史
阶段记录，不得覆盖本节。

### 0.11 2026-08-13 最新权威状态：H10 R0--R4已通过

H10迟滞Tail-to-bulk返回的独立回归、checkpoint可行性审计和动态短程A/B已完成。最终R3采用
三路验收：R3a验证首次返还；R3b的同一100.3 fs预热态none/hysteretic分支只验证物理扰动；
原始100 fs controller-off分支只提供性能基线。最终结果为：

- R3a：`status=PASS`，11/11步接受；
- R3b：`status=PASS`，同预热态两侧均28/28步接受；
- 返回$N/P_x/K$最大残差$2.18\times10^{-15}$，MPI请求残差为0；
- 密度相对$L_2=3.94\times10^{-5}$，场相对$L_2=3.18\times10^{-4}$；
- 场能代理相对差$5.93\times10^{-6}$，场包络相对$L_2=3.82\times10^{-5}$；
- combined能谱相对$L_2=8.86\times10^{-5}$；
- 独立性能基线已接入，Tail内核时间下降41.66%，Tail宏粒子由99,194,116降至60,119,066。

因此R3的守恒、物理不变和性能门均已关闭。随后完成的R4五档参数扫描也已通过：全部运行具有
相同accepted窗口、完整快照且无失败；返回不变量最大残差为$4.67\times10^{-15}$，MPI请求残差为0；
相邻档电场相关系数为$0.999994$--$0.999997$，场相对$L_2$为$0.228\%$--$0.336\%$，场能差不超过
$0.0229\%$，密度差不超过$0.0381\%$，combined谱差不超过$0.1102\%$。故$K_{in}=5.5$ MeV、
$N_{res}=8$登记为准生产中央档。R3/R4历史失败记录只用于解释测试设计演进，不得覆盖本节。

### 0.12 2026-08-14 最新权威状态：H10长期资源门通过，长期物理门未通过

`output/h10_quasiprod_100p3_to_120` 已从预热后的100.3 fs hysteretic checkpoint连续推进到
120 fs。770个时间步全部接受，未发生split、rollback或未解决failure；返回$N/P_x/K$最大
相对残差为$4.98\times10^{-15}$，MPI请求残差为0。因此H10事务、硬守恒和连续推进保持通过。

相对同窗口controller-off历史基线，H10取得了明确资源收益：末态Tail宏粒子由
$1.4317\times10^8$降至$8.0115\times10^7$，减少44.0%；平均wall time由10.918 s/step降至
6.372 s/step，减少41.6%；末100步减少45.6%，单rank最大粒子数约减半。H10资源门通过。

但是，长期物理等价门未通过。120 fs末态相对controller-off参考的电场相对$L_2$为16.29%，
combined能谱为5.18%，$u_\parallel$边缘分布为2.44%，bulk $u_\perp$边缘分布为4.38%；combined
背景动能高3.13%。该差异从约102 fs的0.68%动能偏差单调累计到120 fs，不是单个末态噪声。
密度相对$L_2$仅0.595%，说明主要问题不是电荷守恒，而是PIC粒子返回Eulerian bulk后，虽然
$N/P_x/K$严格闭合，$J_x/\Pi_{xx}/\Pi_\perp$及更高阶分布矩并未逐事件严格保持，随后又切换到
Chang--Cooper bulk碰撞算子，长期产生表示和算子切换偏差。

旧输出还不能关闭能量门。100.3--120 fs按已有列重建的余额约
$4.13\times10^7\ \mathrm{J/m^2}$，但旧诊断没有逐步输出固定电势边界功和完整accepted-step能量
恒等式。现已增加：

- `electrostatic_boundary_work`及完整单步域能量源/余额；
- Tail返回前后$N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$六矩的有符号差；
- `tools/analyze_h10_quasiproduction.py`分离数值、守恒、性能、能量和长期物理五类Gate。

因此，R0--R4结论仍然有效，但只证明短程安全和参数平台；H10不能据此直接签字为最终生产
算法。当前禁止再次跑0--120 fs。下一步只执行§17.15的100.3--105 fs同checkpoint短A/B，先
确定能量余额来自漏记账还是推进误差，并量化每步高阶矩偏差与宏观漂移的相关性。不得使用
全局能量补丁或事后缩放返回分布。

**2026-08-14短A/B结论更新：** 100.3--105 fs两档均完成184个接受步，且无failure、split或
rollback。返回档相对关闭档的场、背景密度、combined能谱、$u_\parallel$和bulk $u_\perp$
相对$L_2$误差分别为1.42%、0.146%、0.206%、0.108%和0.184%，短程物理门通过。返回
$N/P_x/K$事务门也通过。能量余额在关闭档和返回档分别为约$6.218\times10^6$和
$6.255\times10^6\ \mathrm{J/m^2}$，两者仅差约0.6%，故该余额不是H10返回造成，禁止把它
并入返回投影补偿。H10可确认的独立缺陷是局部返回路径跳过表示优化后，$\Pi_\perp$差具有
稳定符号并累计至$-13.61$；修复见§17.15.5。短窗只有4.67 fs，性能数据只作报告，不构成
20%资源门。

---

## 1. 最终物理模型

### 1.1 背景电子总分布

背景电子的物理分布为

$$
f_e=f_{\rm bulk}+f_{\rm tail}.
$$

其中：

- $f_{\rm bulk}$ 使用轴对称 Eulerian 网格表示；
- $f_{\rm tail}$ 使用 1D3V PIC 宏粒子表示；
- 两者共享相同的电子质量 $m_e$、电荷 $-e$、空间边界和电场；
- 两者的密度、动量和动能必须共同进入物理诊断。

Eulerian 主体使用

$$
u_\parallel=\frac{p_x}{m_ec},\qquad
u_\perp=\frac{p_\perp}{m_ec},\qquad
\gamma=\sqrt{1+u_\parallel^2+u_\perp^2},
$$

$$
v_x=c\frac{u_\parallel}{\gamma},\qquad
\dot u_\parallel=-\frac{eE_x}{m_ec},
$$

并采用轴对称测度

$$
d^3u=2\pi u_\perp\,du_\parallel du_\perp.
$$

因此 bulk 是 1D2V 数值表示，代表陀螺对称的 1D3V 分布。tail 粒子则显式保存 $(u_x,u_y,u_z)$，为以后加入俯仰角散射保留完整 3V 自由度。

### 1.2 无碰撞背景动力学

在内部表示转换之外，背景电子满足

$$
\frac{\partial f_e}{\partial t}
+v_x\frac{\partial f_e}{\partial x}
-\frac{eE_x}{m_ec}\frac{\partial f_e}{\partial u_\parallel}=0.
$$

在混合离散中可写为

$$
\frac{\partial f_{\rm bulk}}{\partial t}
+v_x\frac{\partial f_{\rm bulk}}{\partial x}
-\frac{eE_x}{m_ec}\frac{\partial f_{\rm bulk}}{\partial u_\parallel}
=-S_{b\to t}+S_{t\to b},
$$

tail 粒子接收完全相反的内部源。对整个背景电子系统，$S_{b\to t}$ 与 $S_{t\to b}$ 必须严格相消。

### 1.3 碰撞动力学

最终目标包含

$$
C_{FP}[f]
=-\frac{\partial}{\partial\boldsymbol v}\cdot
\left(f\langle\Delta\boldsymbol v\rangle\right)
+\frac{1}{2}
\frac{\partial^2}{\partial\boldsymbol v\partial\boldsymbol v}:
\left(f\langle\Delta\boldsymbol v\Delta\boldsymbol v\rangle\right).
$$

bulk 使用守恒的圆柱速度空间 FP 离散。tail 的碰撞后端由 $\omega$ 的物理类型决定：一般 Kramers--Moyal 系数使用与同一漂移、扩散张量对应的随机微分方程，Coulomb/Landau 核的 tail--tail 碰撞优先使用相对论 Nanbu--Perez 二体 Monte Carlo：

$$
d\boldsymbol u=\boldsymbol A_u\,dt+\boldsymbol B_u\,d\boldsymbol W,
\qquad
\boldsymbol B_u\boldsymbol B_u^T=\boldsymbol D_u.
$$

SDE 后端从 $\boldsymbol v$ 到相对论归一化动量 $\boldsymbol u$ 的变换必须包含 Ito 修正和坐标 Jacobian。禁止把圆柱网格系数直接复制到 Cartesian 粒子增量中。二体 Monte Carlo 后端必须在相对论质心系采样散射，并处理不等宏粒子权重；两种后端不能对同一个 pair contribution 重复生效。

### 1.4 Beam

Beam 是独立的外部注入电子物种：

$$
\frac{dp_{x,p}}{dt}=-eE_x(x_p,t),\qquad
\frac{dx_p}{dt}=c\frac{u_{x,p}}{\gamma_p}.
$$

Beam 在 $0$ 至 $25\ \mathrm{fs}$ 从左边界按通量注入，从任一物理边界离开后删除。Beam 与 background-tail PIC 必须使用不同类、不同状态、不同账本和不同 checkpoint 段。

### 1.5 非周期静电场

场由

$$
\frac{\partial E_x}{\partial x}
=\frac{e}{\varepsilon_0}
\left(Zn_i-n_{\rm bulk}-n_{\rm tail}-n_b\right)
$$

求解。生产默认采用固定电势的非周期静电边界：

$$
\phi(x_{\min})=0,\qquad \phi(x_{\max})=0.
$$

即 `DIRICHLET_PHI`，两端均为零电势。离散 Poisson 解允许两端表面电荷承担域内净电荷对应的法向场跳变，不要求 $E_x$ 在任一端为零。该选择表示计算域两端连接到同一参考电势的静电边界，而不是周期闭合，也不是每步强制净电荷为零。

如果以后切换有限靶 absorbing 模式，必须显式给出电极电势或外部真空闭合；不能只改粒子边界而沿用 reservoir 的场边界解释。

**当前代码状态**：`main_vpfp.cpp::Options` 已默认使用 `DIRICHLET_PHI`，且 `phi_left=phi_right=0`，与本节生产模型一致。H9 和生产命令应显式写出 `--field-boundary dirichlet-phi --phi-left 0 --phi-right 0`，以避免脚本或未来默认值变化造成歧义；manifest 必须记录实际值。`LEFT_E` 只保留为 manufactured-solution/敏感性对照，不是生产默认。

---

## 2. 空间域与边界条件

### 2.1 生产空间网格 [已完成-保留]

- $L_x=40\ \mu\mathrm m$；
- $N_x=8000$；
- $\Delta x=0.005\ \mu\mathrm m$；
- MPI 只沿 $x$ 分区；
- 场和背景电子空间边界均为非周期。

不得把 `dx` 改回 $0.01\ \mu\mathrm m$ 以换取速度。空间收敛测试可通过编译宏修改 $N_x$，但物理长度保持 $40\ \mu\mathrm m$。

### 2.2 背景电子空间边界 [已完成-保留]

当前主模式为 reservoir：

- 域内初始完整充满背景等离子体；
- 左右计算边界外假定存在同温度、同参考密度的背景储库；
- 对流入特征，从储库 Maxwellian 提供入流；
- 对流出特征，使用域内分布自然出流；
- reservoir 交换必须进入粒子数、动量和能量账本。

`reservoir` 不表示在域内额外添加体源，也不表示每步把边界单元重置为 Maxwellian。

### 2.3 tail 粒子空间边界 [已完成-保留]

背景高能 tail 粒子采用开放流出：

- 从左或右边界离开后删除；
- 删除前的最后一段域内轨迹仍参与密度、电流和场功诊断；
- 删除的粒子数、平行动量、动能按左右边界分别累计；
- 不从 reservoir 自动补回对应的高能粒子。

低能 reservoir 入流仍由 bulk 表示。只有当 reservoir 的物理分布本身在转换阈值以上具有不可忽略质量时，初始化阶段才生成对应 tail 入流；100 eV Maxwellian 在 MeV 阈值上的质量通常低于数值可见度，必须通过离散积分确认，不能口头假设。

**当前实现记录**（阶段 H1/H3）：`BackgroundTailPIC` 已实现左右物理边界开放删除、出流前轨迹截断沉积、outside-shape ledger、左右出流数目/动量/能量账本及 MPI 多跳迁移。H1 的 1/2/5 rank 轨迹、连续性和开放出流测试已通过。该状态标记仅说明粒子边界算法完成，不代表 H9 的资源、噪声和碰撞物理已验收。

### 2.4 Beam 空间边界 [已完成-保留]

继续使用当前注入事件表、D-K-D 推进、MPI 迁移和开放删除。重构 tail PIC 时不得改变 Beam 的注入统计、随机数序列和连续性结果。

### 2.5 场边界的物理有效性门

`DIRICHLET_PHI` 是固定电势的物理建模选择，不是数值吸收层。静电 Poisson 系统也不存在通过吸收电磁波消除边界影响的机制。reservoir 粒子边界与固定电势场边界可以共存，但边界粒子携带的电荷、能量通量以及电极静电功必须分别记账。

进入生产前必须完成：

1. manufactured charge 下生产 `DIRICHLET_PHI(0,0)` 的解析/离散 Gauss 测试，并保留 `LEFT_E` 作为非生产对照；
2. reservoir 平衡态无 Beam 测试，确认边界不自发产生场；
3. 一次较短的域长敏感性测试，确认目标核心区在边界到达时间之前对 $L_x$ 稳定；
4. 分别报告边界附近和核心区的场能、净电荷与密度误差。

若核心区对场边界选择敏感，不能用更宽的绘图裁剪区隐藏，必须重新定义物理边界闭合。

---

## 3. 当前代码状态与保留边界

### 3.1 已完成并保留

以下文件的核心算法继续使用：

- `src/open_electrostatic_solver.*`：非周期 Gauss/Poisson 求解；
- `src/open_boundary.*`：reservoir/absorbing 背景空间边界；
- `src/species.*`：背景 cell-integrated mass 存储与矩计算；
- `src/conservative_ppm_remap.*`：背景 $x$ 和 $u_\parallel$ 保守 remap；
- `src/vlasov_split_step.*`：无碰撞 Strang 子步；
- `src/beam_pic.*`：Beam 注入、D-K-D、沉积、迁移和开放删除；
- `src/maxwell.*`：电荷密度容器和场数组；
- `src/vpfp_integrator.*`：当前事务推进框架，可作为混合调度器的基础；
- `src/vpfp_diagnostics.*`：已接受态诊断框架；
- `src/vpfp_checkpoint.*`：checkpoint 框架，但格式必须升级。

### 3.2 只扩展接口，不重写算法

1. `EMFields::set_charge_density()` 增加 tail 密度参数。
2. `VpfpIntegrator` 增加 tail 状态和转换器，不改变 Poisson 算法。
3. `VpfpDiagnostics` 增加 tail 与转换诊断，不改变只写已接受态的原则。
4. `VpfpCheckpoint` 增加 tail 粒子、ID、转换参数和账本。

### 3.3 必须替换或废止

- 以持续扩大 `FP_VELOCITY_GRID_UPAR_CORE_MAX` 和 `FP_VELOCITY_GRID_NV` 追逐高能尾部的生产策略；
- 只凭单步 `umax_loss_abort_fraction` 判定物理解失败；
- 把高能尾部删除后仅记作数值损失的处理；
- 当前仅支持零系数或常数对角系数的碰撞生产路线；
- 碰撞系数只在单一速度点计算并用于整个速度网格的占位实现。

### 3.4 EPOCH 4.20.1 源码审计结论与适用边界

本节依据本机只读参考源码：

```text
E:/ScientificComputation/epoch-4.20.1/epoch_release-4.20.1/epoch1d/src
```

审计对象和结论如下。

参考树的 `epoch1d/src/COMMIT` 为 `v4.20.1-0-gf294c484-clean`。后续若更换 EPOCH 版本，必须重新核对本节行号和默认碰撞后端，不能把本节结论无版本地外推。

| 功能 | EPOCH 源码位置 | 对本项目的结论 |
|---|---|---|
| 相对论粒子推进 | `particles.F90:28-575` | EPOCH 先半漂移到场采样位置，再用 stagger-aware 权重插值场、更新三维动量并完成后半漂移。纯静电 $B=0$ 时，本项目现有 D-K-D 是该时间层结构的专化形式，无需引入完整 Boris 旋转。 |
| staggered $E_x$ gather | `particles.F90:284-338`、`include/triangle/e_part.inc` | EPOCH 对 $E_x$ 使用相对 cell-centered 形函数平移半格的 `hx`，而不是先把面场平均成 cell 场。tail 必须复用本项目已经验证的 shifted face gather。 |
| 电荷守恒电流 | `particles.F90:443-507` | EPOCH 用末态形函数与初态形函数之差 `hx-gx` 累积 $J_x$，实质上离散求解连续性方程，不使用瞬时 $qv$ 矩。本项目保留真实轨迹分段沉积；它允许多 cell 轨迹，比 EPOCH 该实现所假定的“单步最多跨一格”更一般。 |
| 粒子开放边界 | `boundary.F90:661-877` | EPOCH 在完成域内轨迹电流后处理粒子边界；默认非周期越界粒子被删除，内部 MPI 边界粒子发送给相邻 rank。tail 也必须先完成截断轨迹与出流账本，再删除粒子。 |
| 粒子电流边界 | `boundary.F90:372-458,881-890` | 周期物种会合并周期 guard；非周期物种不回卷并清除 guard。本项目为非周期 Poisson 场，tail 轨迹电流只作连续性/功审计，绝不能周期闭合。 |
| MPI 粒子迁移 | `boundary.F90:849-869` | EPOCH 使用相邻 rank 点对点交换。其正确性依赖粒子 CFL；本项目必须显式保证一次子步最多跨一个 rank，或实现确定终止的多跳迁移循环。 |
| cell-local 碰撞配对 | `epoch1d.F90:220-256`、`housekeeping/secondary_list.F90` | 碰撞前按粒子当前位置重排到空间 cell，碰撞完成后再接回主粒子表并重新做边界处理。tail 碰撞也必须按物理 cell/碰撞 bin 分组，不能跨远距离 cell 随机配对。 |
| 默认二体碰撞 | `physics_packages/collisions.F90`、`deck/deck_collision_block.F90:118-137` | 4.20.1 默认 `use_nanbu=T`，采用相对论 Nanbu--Perez；Sentoku--Kemp 是可选旧后端。它是 Coulomb/Landau 碰撞的 Monte Carlo 近似，不是任意给定 Kramers--Moyal 系数的通用实现。 |
| 快粒子--背景碰撞 | `physics_packages/background_collisions.F90` | 该路径假定慢背景静止，只改变快粒子，不提供完整背景反作用。它可作 trace-tail 或高能极限基准，不能直接作为本项目完整 tail--bulk 生产碰撞。 |
| 注入 | `physics_packages/injectors.F90:132-323` | EPOCH injector 与外部 Beam 有参考价值，但 background-tail 由域内表示转换产生，不是边界注入；tail 类不得继承 injector 状态。 |

明确保留的改进：

1. EPOCH 默认随机数依赖粒子列表顺序和 shuffle；本项目继续使用基于 particle ID/accepted step 的 counter-based RNG，以获得 MPI/OpenMP 可重启性。
2. EPOCH 默认三点抛物线粒子形函数不应在 H1 阶段盲目替换本项目已验证的 CIC。H1 先复用 Beam 的 CIC/shifted gather；只有独立 TSC 电荷、场插值、边界和 MPI 测试通过后，才能进行 CIC/TSC A/B。
3. EPOCH 的碰撞密度估计与场沉积形函数不是同一个接口。tail 碰撞率应使用明确记录的 cell-local number/volume 估计，不能误把 Poisson 的平滑密度数组当成宏粒子配对权重。
4. EPOCH 的全 PIC 二体碰撞天然拥有两个粒子化碰撞对象；本项目的 tail--bulk 是跨表示耦合，必须额外构造 Eulerian 反作用，不能据 EPOCH 结果省略。
5. EPOCH 使用链表和临时 `secondary_list` 按 cell 重排，这是算法语义参考而非 C++ 性能模板。本项目应保留连续粒子存储，使用预分配 cell offsets/index permutation 完成碰撞分箱，避免每步粒子节点分配。

---

## 4. 混合背景电子的离散不变量

### 4.1 统一物理矩

所有物理诊断使用

$$
N_e=N_{\rm bulk}+N_{\rm tail},
$$

$$
P_{e,x}=P_{\rm bulk,x}+P_{\rm tail,x},
$$

$$
K_e=K_{\rm bulk}+K_{\rm tail}.
$$

Poisson 密度必须使用

$$
n_e(x)=n_{\rm bulk}(x)+n_{\rm tail}(x).
$$

禁止只输出 bulk 平均动能却标记为“背景电子平均动能”。应同时输出 bulk、tail 和 combined 三套量。

统一单位契约如下：

- `Species::f` 的 cell-integrated mass 是每单位横向面积的电子数，单位为 $\mathrm{m}^{-2}$；
- `BackgroundTailParticle::weight` 与 `BeamParticle::weight` 使用同一单位 $\mathrm{m}^{-2}$；
- 粒子向空间 cell 沉积时除以 $\Delta x$，得到 $\mathrm{m}^{-3}$；
- 全域动能使用 $\mathrm{J/m^2}$；
- 电流密度使用 $\mathrm{A/m^2}$。

实现时在类型注释、checkpoint schema 和单元测试中都写明该单位，禁止让 converter 自行引入另一套“宏粒子数”归一化。

### 4.2 内部转换守恒

每次 bulk-to-tail 转换必须满足

$$
\Delta N_{\rm bulk}^{\rm conv}+\Delta N_{\rm tail}^{\rm conv}=0,
$$

$$
\Delta P_{\rm bulk,x}^{\rm conv}+\Delta P_{\rm tail,x}^{\rm conv}=0,
$$

$$
\Delta K_{\rm bulk}^{\rm conv}+\Delta K_{\rm tail}^{\rm conv}=0.
$$

转换不是注入、碰撞、场功或边界损失，不得进入外部能量源项。

### 4.3 Poisson 电荷等价

转换前后的离散空间电荷应满足

$$
\left\|\rho_{\rm before}-\rho_{\rm after}\right\|_2
\le \epsilon_\rho
\max(\|\rho_{\rm before}\|_2,\rho_{\rm floor}).
$$

这一要求比全局电子数守恒更严格。即使总电子数相同，若宏粒子位置和形函数改变了局部密度，Poisson 场也会立即出现人工脉冲。

### 4.4 事务性

一次时间步必须维护：

- `bulk_n` 与 `bulk_trial`；
- `tail_n` 与 `tail_trial`；
- `beam_n` 与 `beam_trial`；
- `fields_n` 与 `fields_trial`；
- tail 粒子 ID 计数器、随机数状态和转换账本的 accepted/trial 副本。

只有全部验收通过后才统一 `swap`。任何失败分支不得消耗永久粒子 ID、推进永久 RNG 或累加永久账本。

### 4.5 combined 离散连续性

转换必须作为两个表示之间符号相反的内部源出现。以数通量 $\Gamma$ 表示：

$$
\frac{n_{\rm bulk}^{n+1}-n_{\rm bulk}^n}{\Delta t}
+D_x\Gamma_{\rm bulk}
=S_{\rm reservoir}-S_{b\to t}+S_{t\to b},
$$

$$
\frac{n_{\rm tail}^{n+1}-n_{\rm tail}^n}{\Delta t}
+D_x\Gamma_{\rm tail}
=S_{b\to t}-S_{t\to b}-S_{\rm tail,out}.
$$

相加后，内部转换源必须逐 cell 抵消。`BackgroundTailPIC` 必须保存完整轨迹电流和 conversion source density；不能只比较时间步前后的总粒子数。

---

## 5. 有限 bulk 速度域与转换阈值

### 5.1 设计目标

混合方案不是把当前 $u_{\parallel,\max}=20$ 的大网格原样保留后再增加 PIC。bulk 网格只负责热核、波粒共振主体和阈值前缓冲区；tail PIC 负责稀疏高能电子。

必须同时满足：

1. bulk 网格解析热核和主要共振区；
2. 转换阈值位于 bulk 的可信解析区内；
3. 阈值与数值速度边界之间保留缓冲；
4. 任何到达最外速度面的 bulk 质量都属于实现失败，而不是可接受尾损失。

### 5.2 阈值定义

阈值使用离散相对论动能：

$$
K_{jk}=m_ec^2\left(\sqrt{1+u_{\parallel,j}^2+u_{\perp,k}^2}-1\right).
$$

运行参数使用物理单位：

```text
--background-tail-mode pic
--tail-convert-energy-mev <K_out>
```

首轮候选值可以覆盖 $4$、$6$、$8\ \mathrm{MeV}$，但不得把其中一个值直接硬编码成最终物理参数。最终阈值由以下三项共同决定：

- 阈值以下 bulk 结果的网格收敛；
- 阈值变化时 combined 密度、能谱、场能和能量转移的稳定性；
- tail 粒子数和噪声的可接受性。

一个生产运行中的 $K_{out}$ 必须固定。禁止根据瞬时粒子数、场幅或与 EPOCH 的差异自动移动阈值。动态阈值只有在另行证明移动界面本身严格守恒且不产生谱形滞后后才能引入。

### 5.3 网格要求

第一版实现可以暂时使用当前网格验证转换器，但不能直接据此进行 120 fs 生产。进入生产前应建立有限的 `resolved + buffer` 网格：

- 中心热核保持当前已验证分辨率；
- 阈值位于 resolved 区，不位于最后一个可信单元；
- 阈值外至少保留 8 至 16 个平滑变宽的 $u_\parallel$ 缓冲单元；
- $u_\perp$ 外边界也必须有占据率和通量诊断。

禁止恢复 U40R 一类全局高分辨率大域作为默认生产配置。

### 5.4 阈值掩码

`SpatialGrid` 或新的 `HybridVelocityPartition` 必须预计算：

- 每个速度单元的 $K_{jk}$；
- `is_bulk_resolved[jk]`；
- `is_conversion_cell[jk]`；
- `is_buffer_cell[jk]`；
- 速度网格和阈值配置哈希。

这些数组初始化一次后只读共享，禁止在每个时间步重复计算。

首版以速度 cell center 的 $K_{jk}$ 对完整 cell 分类，不在一个速度 cell 内按连续阈值切割质量。这样转换对象仍是明确的 cell-integrated mass。阈值界面因此是离散阶梯面，其误差必须通过阈值和速度网格收敛测试量化。禁止在没有子 cell 重构和守恒积分的情况下按几何面积比例切 cell。

### 5.5 初始尾部

初始化 Maxwellian 后应积分阈值以上离散质量：

- 若小于 `initial_tail_quadrature_tolerance`，记录为初始化离散截断，不生成大量极小权重粒子；
- 若不可忽略，必须通过转换器生成初始 tail 粒子；
- 该容差必须相对总背景电子数定义，不得用绝对浮点数随网格改变。

---

## 6. BackgroundTailPIC 设计

### 6.1 新文件 [已完成-保留]

新增：

```text
src/background_tail_pic.h
src/background_tail_pic.cpp
tests/background_tail_pusher_test.cpp
tests/background_tail_deposition_test.cpp
tests/background_tail_shape_difference_test.cpp
tests/background_tail_stagger_test.cpp
tests/background_tail_open_boundary_mpi_test.cpp
```

不要继承 `BeamPIC`。可以复用经过测试的形函数、MPI 迁移和开放边界辅助函数，但 Beam 的注入状态不能进入 tail 类。

### 6.2 粒子结构

```cpp
struct BackgroundTailParticle {
    double x;
    double ux;
    double uy;
    double uz;
    double weight;
    std::uint64_t id;
};
```

`BackgroundTailPIC` 至少保存：

- 本 rank 粒子数组；
- cell-centered density；
- face current 或轨迹电流诊断；
- 左右出流数目、动量和能量；
- 当前最大动量和能量；
- 唯一 ID 生成状态；
- 用于碰撞的 counter-based RNG key；
- accepted/trial 双缓冲。

tail 的 `weight` 与现有 Beam 权重相同，表示每单位横向面积的真实电子数。不得保存“还需要乘 $\Delta x$”的隐含权重。

tail 内部动量统一保存为无量纲 $\boldsymbol u=\boldsymbol p/(m_ec)$。与参考 EPOCH 的 SI 动量 `part_p` 交换公式时，必须在一个命名明确的转换函数内完成单位换算；禁止在 pusher、碰撞和诊断中混用 `p` 与 `u`。

### 6.3 无碰撞推进

tail 使用相对论 D-K-D：

$$
x^{n+1/2}=x^n+\frac{\Delta t}{2}c\frac{u_x^n}{\gamma^n},
$$

$$
u_x^{n+1}=u_x^n-\frac{e\Delta t}{m_ec}E_x^{n+1/2}(x^{n+1/2}),
$$

$$
x^{n+1}=x^{n+1/2}+\frac{\Delta t}{2}c\frac{u_x^{n+1}}{\gamma^{n+1}}.
$$

$u_y$、$u_z$ 在无碰撞电场推进中保持不变。

场插值与密度沉积必须使用同阶成对形函数。首版优先使用已验证、低通信成本的 CIC；若采用 TSC，必须同时实现 guard 宽度、MPI 求和和边界截断测试。

轨迹电流必须沿真实粒子路径逐段沉积。若单步跨越多个 cell 或 MPI 子域，应在每个穿越面分段，不能只用端点速度矩。该电流不用于推进 Poisson 场，但用于 combined 连续性和场功审计。

实现时把 EPOCH 的离散思想映射为以下本项目接口：

```cpp
struct ParticleShape1D {
    static FaceGatherWeights shifted_face_weights(double x,
                                                   const SpatialGrid& grid);
    static CellDepositWeights cell_weights(double x,
                                           const SpatialGrid& grid);
};

struct ChargeConservingTrajectory1D {
    static void deposit_segment(double x0,
                                double x1,
                                double charge_weight,
                                double dt,
                                FaceCurrentAccumulator& current);
};
```

约束如下：

1. `shifted_face_weights()` 必须直接作用于 `Ex_face`；禁止先构造临时 cell-centered $E_x$。
2. `cell_weights()` 同时供 tail midpoint/final density、转换前后电荷等价测试和连续性重构使用。
3. `deposit_segment()` 以形函数变化满足逐 cell 离散连续性；多 cell 路径先在真实穿越面切段。
4. 开放出流段只沉积到首次物理边界穿越时刻；边界外路径不沉积、不回卷。
5. 内部 MPI 边界不是物理边界。轨迹段和粒子所有权都必须传给邻居，不能把 rank 边界误记为出流。
6. 若选用 CIC，H1 必须证明该实现与现有 Beam 的 shifted gather 和 density shape 完全一致；若抽取公共 helper，Beam 回归 hash/矩结果必须不变。
7. 物理边界附近不得把被截断的粒子形函数重新归一化到域内。域内 shape 份额进入 Poisson 密度，guard/outside 份额进入边界 shape ledger；二者之和才等于粒子总权重。否则粒子靠近边界但尚未出流时会产生人工电荷源。

### 6.4 粒子 ID 与随机数

ID 必须在 MPI 和 OpenMP 下确定性唯一。推荐：

```text
high bits: creation MPI rank or global creator block
low bits : accepted conversion counter
```

碰撞随机数使用 `(particle_id, accepted_step, collision_half, component)` 作为 counter-based key。禁止依赖粒子 vector 顺序，否则 MPI 迁移和线程调度会改变物理解。

### 6.5 粒子权重与噪声

tail 不能使用任意宽权重分布。必须输出：

- `weight_min/max/mean/std`；
- 每个空间单元的有效宏粒子数；
- 密度噪声估计；
- 最大单粒子电荷占局部背景电荷的比例。

首版允许确定性 splitting，禁止直接实现会改变数目、动量或能量的启发式 merging。需要 merging 时单独设计守恒测试。

### 6.6 EPOCH 参考下的 H1 附加验收

H1 除原有测试外必须新增：

1. **stagger test**：给定非线性 `Ex_face`，逐粒子比较公共 shifted gather 与手算 face 权重，防止退化为 cell 平均。
2. **shape-difference continuity test**：随机生成不跨格、跨一格和跨多格轨迹，验证 $\Delta\rho/\Delta t+D_xJ=0$ 达到求和误差。
3. **open truncation test**：粒子从边界内半个 shape 宽度逐步移动到域外，验证域内密度、outside-shape ledger、截断轨迹与最终出流账本连续；出流粒子从末态密度消失且无周期回卷，边界 shape 不做域内重归一化。
4. **MPI ownership test**：同一轨迹分别以 1、2、5 rank 运行，最终粒子、密度、轨迹电流和出流账本一致。
5. **CFL contract test**：若生产迁移只支持相邻 rank，启动时验证 $c\Delta t<L_{x,\mathrm{rank}}$；否则测试多跳迁移并设置最大跳数失败保护。
6. **electrostatic specialization test**：$B=0$ 时 D-K-D 结果与完整相对论 Boris 的静电专化结果一致，证明无需把 EPOCH 的磁旋转代码移植进来。

---

## 7. bulk-to-tail 保守转换器

### 7.1 新文件 [已完成-保留]

```text
src/bulk_tail_converter.h
src/bulk_tail_converter.cpp
tests/bulk_tail_single_cell_test.cpp
tests/bulk_tail_multibin_test.cpp
tests/bulk_tail_poisson_invariance_test.cpp
tests/bulk_tail_transaction_test.cpp
```

建议接口：

```cpp
struct BulkTailConversionDiagnostics {
    double number_removed;
    double number_created;
    double px_removed;
    double px_created;
    double energy_removed;
    double energy_created;
    double rho_l2_before_after;
    double rho_linf_before_after;
    std::uint64_t particles_created;
    bool finite;
    bool conservative;
};

class BulkTailConverter {
public:
    BulkTailConversionDiagnostics extract_after_substep(
        Species& bulk_trial,
        BackgroundTailPIC& tail_trial,
        const SpatialGrid& grid,
        const HybridVelocityPartition& partition,
        int accepted_step,
        ConversionLocation location);
};
```

### 7.2 转换时间层

无碰撞首版只在完整 $u_\parallel$ 受力子步结束后转换：

1. bulk 已经使用 $E^{n+1/2}$ 完成 $u$ remap；
2. tail 中原有粒子已完成同一个场的 kick；
3. 扫描 `bulk_trial` 的转换单元；
4. 从 bulk cell-mass 中扣除转换质量；
5. 在 $x^{n+1/2}$ 时间层生成新 tail 粒子；
6. 新粒子不再重复 kick，只参加第二个空间半漂移。

这样新粒子在前半步仍由 bulk 表示，在后半步由 PIC 表示，不会重复受力或漏掉受力。

碰撞启用后，在每个 bulk 碰撞半步结束后也必须执行转换，因为碰撞扩散可能把质量送入转换区。

### 7.3 不允许逐 cell 随机采样

每个非零速度单元直接随机生成粒子会导致：

- 粒子数不可控；
- 小权重跨度巨大；
- Poisson 噪声过强；
- 重启不可复现。

应先在每个本地空间单元内按 `u_parallel` 符号和能量区间聚合转换质量。聚合 bin 的分辨率通过运行参数控制：

```text
--tail-conversion-upar-bins <N>
--tail-conversion-energy-bins <N>
```

禁止把整个高能尾部压成一个粒子 bin。转换能谱必须随 bin 数收敛。

### 7.4 正确性参考加载

第一版正确性参考不能直接把一个聚合 bin 压成 4 个粒子。参考加载对每个被移除的 $(x,j,k)$ cell-mass $M_{xjk}$ 生成一个方位角 quartet：

$$
u_x=u_{\parallel,j},\qquad
(u_y,u_z)=u_{\perp,k}(\cos\phi_p,\sin\phi_p),
$$

$$
\phi_p\in\{0,\pi/2,\pi,3\pi/2\},\qquad w_p=M_{xjk}/4.
$$

该参考加载不仅保留 $N/P_x/K$，还保留离散 $v_x$ 电流矩、平行应力和横向应力。它可能产生较多粒子，但必须先作为 H2 的黄金参考。

### 7.5 生产用矩约束压缩

直接使用原文档中的

$$
\bar u_x=P_x/(m_ecN),\qquad
u_{\perp,\mathrm{eff}}^2=\bar\gamma^2-1-\bar u_x^2
$$

并让所有粒子具有相同 $u_x$，虽然能保留 $N/P_x/K$，却会消除 bin 内平行速度方差并把部分平行能量错误地转成横向能量。该算法不得作为生产转换器。

生产转换采用两层流程：

1. 先构造第 7.4 节的可行非负参考粒子集；
2. 再在同一个空间沉积 stencil、$u_\parallel$ 符号和能量 bin 内做确定性的稀疏矩压缩。

每个压缩组至少保持以下 cell-integrated 离散矩：

$$
\mathcal{M}=
\left(
N,
P_x,
\mathcal{J}_x,
K,
\mathcal{\Pi}_{xx},
\mathcal{\Pi}_{\perp}
\right),
$$

其中

$$
\mathcal{J}_x=J_x\Delta x
=-e\sum_p w_p c\frac{u_{x,p}}{\gamma_p},
$$

$$
\mathcal{\Pi}_{xx}=\Pi_{xx}\Delta x
=m_ec^2\sum_p w_p\frac{u_{x,p}^2}{\gamma_p},
$$

$$
\mathcal{\Pi}_{\perp}=\Pi_{\perp}\Delta x
=m_ec^2\sum_p w_p\frac{u_{y,p}^2+u_{z,p}^2}{\gamma_p}.
$$

需要输出物理密度量时再除以 $\Delta x$。禁止在 converter 一侧使用积分矩、在 tail 诊断一侧使用密度矩却直接比较。

推荐使用候选参考速度上的确定性 Carathéodory/null-space 消元。每个候选支撑点是一个完整 gyrophase quartet，而不是 quartet 中的单个粒子，因此横向净动量始终为零。

具体算法：

1. 用六个尺度归一化后的约束矩构造 $A\in\mathbb R^{6\times m}$，初始非负权重为 $w$，满足 $Aw=\mathcal{M}_{\rm ref}$。
2. 当活动支撑数大于约束秩加一时，对活动列使用固定 pivot 规则的 QR/SVD 求非零 null vector $z$，满足 $Az=0$。
3. 统一 $z$ 的符号后取
   $$
   \theta=\min_{z_i>0}\frac{w_i}{z_i}.
   $$
4. 更新 $w\leftarrow w-\theta z$。至少一个权重变为零，$Aw$ 保持不变，且其他权重仍非负。
5. 删除舍入级零权重，重复至活动支撑数不超过 7 个 ring support。
6. 用未缩放物理矩重新计算全部残差，不允许只检查缩放后的线性系统。

常见组最多保留 7 个 quartet，即 28 个粒子。若矩矩阵数值秩异常、残差超限或出现实质负权重，则保留未压缩参考 quartet，并记录 `conversion_compression_fallback=1`；不能降低矩约束换取更少粒子。

### 7.6 离散矩和压缩验收

对一个转换 bin，依据被移除的离散 cell-mass 计算：

$$
N=\sum_a M_a,
$$

$$
P_x=m_ec\sum_a M_a u_{\parallel,a},
$$

$$
K=m_ec^2\sum_a M_a(\gamma_a-1).
$$

转换器必须分别报告 $N/P_x/\mathcal{J}_x/K/\mathcal{\Pi}_{xx}/\mathcal{\Pi}_\perp$ 的 requested、created 和 residual。$N/P_x/K$ 是硬守恒门；$\mathcal{J}_x$ 和两个应力矩是物理保真门。后者阈值由参考加载和压缩求解精度确定，不能只用“bin 足够窄”代替测量。

### 7.7 空间沉积等价

若所有转换粒子放在 cell center，必须验证所用沉积形函数能把总权重严格返回原空间单元。若形函数会分配到相邻单元，则使用成对确定性位置，使转换后的离散空间密度等于转换前被移除的密度。

不得在转换后通过 Poisson 解的结果反推或修补粒子权重。电荷等价必须由加载构造本身保证。

对于 CIC，同一个空间 cell 内的粒子压缩还必须保持总权重和加权位置 $\sum_p w_p x_p$。因为 CIC 权重在一个 cell 内对 $x$ 是线性的，这两项可以保持相邻两个 cell 的离散沉积。跨越不同 CIC stencil 的粒子不得放入同一个压缩组。

### 7.8 转换后的 bulk 清理

只有已经转为粒子的 mass 才能从 `bulk_trial.f` 扣除。禁止：

- 逐点 `max(f,0)`；
- 把转换区直接清零但按另一个积分生成粒子；
- 扣除后重新归一化整个分布；
- 将矩残差分摊到热核。

PIC 粒子权重必须非负。若待转换 cell-mass 出现负值：

- 舍入级负值必须由现有保守 remap 的局部 roundoff 处理路径在进入转换器前解决，并保留修正账本；
- 超过舍入容差的负值说明 bulk 推进本身失败，转换器必须返回失败；
- 禁止创建负权重 tail 粒子，也禁止用相邻正质量在转换器内部静默抵消。

转换完成后重新计算 bulk moments，并独立计算 conversion residual。

### 7.9 首版返回策略

无碰撞阶段首版采用 **单向永久转换**：一旦进入 tail PIC，即使后续被减速，也暂不投回 Eulerian bulk。

原因是未经验证的 PIC-to-grid 返回更容易造成负分布、能量修补和阈值抖动。首版必须先验证单向转换的正确性和尾部规模。

但碰撞生产前必须检查：

- tail 数目占背景总数的比例；
- tail 中低于候选返回阈值的比例；
- tail 粒子在热核速度区的停留时间；
- tail PIC 对核心密度噪声的贡献。

`tail number fraction=1\%` 和 `thermalized-tail fraction=10\%` 只能作为首轮资源/有效性预警值，不是已验证物理常数。H4/H9 必须根据 PIC 噪声、内存、碰撞反作用和参数收敛确定最终门槛。超过当前配置门槛时必须实施第 16 节的双向返回或停止该配置，不能继续永久 PIC 化。

### 7.10 tail 粒子人口控制

只压缩“本步新转换粒子”不足以保证 120 fs 粒子数受控。新增 `TailPopulationController`，但只有在 H4 已证明基础转换正确后接入。

控制器只能在同一物理空间 cell、相同 CIC stencil 和相邻 tail 相空间 bin 内工作。其输入是现有粒子加新粒子，输出继续满足第 7.5 节的非负矩约束，并额外保持 $\sum_p w_p x_p$。

新增：

```text
src/tail_population_controller.h
src/tail_population_controller.cpp
tests/tail_population_controller_test.cpp
tests/tail_population_controller_mpi_test.cpp
```

控制策略：

- 低于 `target_particles_per_phase_bin` 时不操作；
- 超过上限时做确定性守恒压缩；
- 权重比过大时做等权或近等权 splitting，子粒子总权重、位置矩和速度矩不变；
- 边界 cell、MPI 迁移中的粒子和本步即将出流粒子不参与合并；
- 每次压缩输出 $N/P_x/\mathcal{J}_x/K/\mathcal{\Pi}/\rho$ 残差；
- 失败时保留原粒子，不得接受低阶矩匹配结果。

粒子人口控制是表示压缩，不是碰撞，不能改变能谱或产生熵增。

### 7.11 bulk-to-tail 阈值接口专项修复 [当前最高优先级]

#### 7.11.1 已知缺陷与禁止操作

H9 的 $K_{out}=6$ MeV、controller-off 结果在 6.0--6.2 MeV 出现硬空 bin，
而 6.2--6.4 MeV 又明显回升。现有转换器虽然保持
$N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$ 六个组积分矩，但当前
`HybridVelocityPartition::energy_bin()` 只把整个
$[K_{out},K_{max}]$ 均匀分成少量组，随后每组最多保留 7 个 ring support。
因此，“总矩闭合”并不能证明阈值附近的局部能谱被保留。

最新证据要求把问题拆成两个互相独立的层次：

- **压缩分组误差**：旧均匀大组会跨阈值诊断 bin 搬运谱质量；阈值感知粗组已
  修复这一层，并在组积分意义下通过；
- **活跃近轴路径的离散支撑误差**：生产高能尾部主要沿低 $u_\perp$ 环进入
  转换区，而最低环的相邻中心能量从约 5.809 MeV 直接跳到 6.228 MeV。
  这一层不是增加压缩 support 就能修复，必须先确认转换前离散占据和源谱。

本轮禁止：

- 直接平滑最终能谱、向空 bin 人工补粒子或跨 bin 搬运权重；
- 修改 $K_{out}$ 以绕开 6 MeV 谱谷；
- 通过扩大绘图 bin 掩盖谱谷；
- 对转换后粒子做全局能量/动量补丁；
- 同时修改碰撞核、Sentoku--Kemp 修正或 population controller；
- 未定位前实现复杂 cut-cell 转换或全局增加 `Nv/Nmu`。

#### 7.11.2 第一步：新增只读速度支撑审计

新增：

```text
tests/hybrid_threshold_support_audit_test.cpp
```

并在 `CMakeLists.txt` 增加同名可执行目标。该测试必须直接构造生产使用的
`CylindricalVelocityGrid` 和 `HybridVelocityPartition`，不得复写网格公式。

现有全局支撑审计已实现，但它把全部 $u_\perp$ 环合并，不能作为真实 Beam
活跃路径的验收。必须扩展为“逐 $u_\perp$ 环 + 占据加权”的只读审计，并对
$K_{out}\pm1$ MeV 范围输出：

```text
bin_low_mev bin_high_mev
bulk_center_count conversion_center_count
min_center_energy_mev max_center_energy_mev
max_uncovered_energy_gap_mev
phase_volume_sum
uperp_index uperp_center
nearest_below_Kout_mev nearest_above_Kout_mev crossing_gap_mev
pre_conversion_number pre_conversion_energy_J
```

同时输出全部 conversion cell-center 能量的有序唯一值及相邻间距。至少使用
0.05、0.1、0.2 MeV 三档只读诊断 bin。测试不推进方程、不创建粒子、不改变
生产参数。逐环表至少覆盖 `k=0..7`，并额外输出按实际转换前质量加权的全局
统计；不得再用“所有环合并后的中心计数”替代近轴支撑判断。

判定：若活跃近轴环的 6.0--6.2 MeV 内本来就没有或几乎没有 conversion
cell-center
支撑，而 6.2--6.4 MeV 有大量支撑，则现有 0.2 MeV 直方图不能被当作连续谱
物理量，首要问题是局部速度格点/阈值单元表示；若提取前该区有非零占据而
创建后消失，则问题位于 extraction、聚合压缩或转换时间序列。

#### 7.11.3 第二步：建立转换器三路独立 A/B/C

新增：

```text
tests/bulk_tail_threshold_interface_test.cpp
tests/bulk_tail_threshold_interface_mpi_test.cpp
```

在 `BulkTailConverter` 增加仅供测试选择的加载策略枚举，不通过全局变量或
编译宏切换：

```cpp
enum class BulkTailLoadingPolicy {
    GOLDEN_QUARTETS_NO_COMPRESSION,
    CURRENT_PRODUCTION_COMPRESSION,
    THRESHOLD_AWARE_COMPRESSION
};
```

阶段门要求：专项测试通过前，生产基线应保持
`CURRENT_PRODUCTION_COMPRESSION`，或将 threshold-aware 明确标记为 candidate
并禁止进入正式长跑。当前代码已提前使用 candidate，不能据此反向宣称验收完成。策略只改变
候选支撑压缩，不得改变 bulk 质量提取、转换时间层、空间位置、粒子 ID、沉积
或事务语义。

测试输入不得使用“所有 conversion cell 等质量”的不连续构造。应在多个 x cell
上构造至少四类解析光滑、严格非负的 cell-integrated 分布：

1. $f(K)\propto\exp(-K/T)$ 的单调谱；
2. 跨 $K_{out}$ 的宽高斯谱包；
3. 带 $u_\parallel$ 漂移且 $u_\perp\ne0$ 的各向异性谱。
4. `near-axis-narrow`：质量集中在最低若干 $u_\perp$ 环，$u_\parallel$ 谱包横跨
   $K_{out}$，其宽度和漂移应覆盖生产 10.4--12 fs 首次转换时的活跃单元。

每个 case 保存转换前的离散参考谱，并分别执行：

- A：不压缩黄金 quartet；
- B：当前生产压缩；
- C：阈值感知压缩。

所有谱必须使用同一组显式边界，并同时报告 raw count、bin width 和
$dN/dK$。输出至少包括：

```text
N/Px/Jx/K/Pixx/Piperp residual
rho_l2/rho_linf
spectrum_L1_relative
spectrum_Linf_relative
near_threshold_bin_relative[...]
near_threshold_log_curvature_max
particles_created
compression_fallback_count
```

MPI 测试以同一个全局解析状态分别使用 1、2、5 rank，比较全局谱、六矩和粒子
计数。随机数不得参与该测试。

两个测试程序必须提供稳定 CLI，便于集群脚本逐档执行：

```text
bulk_tail_threshold_interface_test
  --case <smooth-exp|broad-gaussian|anisotropic-drift|near-axis-narrow|all>
  --policy <golden|current|threshold-aware|all>
  --bin-width-mev <value>
  --result <path>

bulk_tail_threshold_interface_mpi_test
  --case all
  --policy <golden|threshold-aware>
  --bin-width-mev <value>
  --result <path>
```

所有测试只在 rank 0 写 `.result`，末行统一为 `status=PASS|FAIL`；任一 rank 的
有限性、守恒或事务检查失败，必须通过一次集体归约使全部 rank 返回非零，禁止
让部分 rank 提前退出。结果文件必须使用截断写入，禁止追加旧结果。PASS 必须
逐策略显式检查对应门槛，禁止使用“非 golden 路无条件忽略谱误差”的逻辑。

#### 7.11.4 第三步：根据 A/B/C 结果选择修复分支

**分支 A：黄金 quartet 已经出现同样谱谷。**

这说明压缩不是主因。先比较“转换前离散中心谱”和“连续 cell-volume 积分
参考谱”：

- 若离散中心谱与转换后黄金谱一致，而两者都相对连续参考出现谱谷，则属于
  当前速度网格对 0.2 MeV 分箱的欠分辨。此时只允许采用局部阈值网格设计：
  在 $K_{out}\pm\Delta K$ 对应的 $(u_\parallel,u_\perp)$ 区域增加有限支撑，
  保持热核网格、外速度边界和总内存基本不变；不得全局扩大 `Nv/Nmu`。
- 若转换前离散谱连续、黄金转换后不连续，则检查
  `is_conversion_cell`、`extract_conversion_masses()` 和粒子能量构造是否使用
  完全相同的 `kinetic_energy[j,k]`。不得进入压缩修复分支。

对真正跨越 $K_{out}$ 的 cut cell，只有在局部网格仍不能达到收敛时，才实现
基于非负重构和守恒子单元积分的 cut-cell converter。禁止按几何面积比例直接
切质量。该项属于后续独立设计，不是本轮默认方案。

**分支 B：黄金 quartet 连续，而当前生产压缩产生谱谷。**

这是当前最可能且改动最小的分支。修改：

- `src/grid.h`：为 `HybridVelocityPartition` 增加初始化后只读的
  `conversion_energy_edges`。边界必须显式包含
  $K_{out}$、$K_{out}+0.2$ MeV、$K_{out}+0.4$ MeV、
  $K_{out}+0.8$ MeV；更高能区可使用对数增宽。不得继续只用 4 个覆盖
  $[K_{out},K_{max}]$ 的等宽大组。
- `energy_bin()` 改为对上述边界执行 `upper_bound`，并给出端点单元测试。
- `src/bulk_tail_converter.cpp`：`GroupKey::energy_bin` 使用新边界；任何压缩组
  不得跨越阈值审计 bin。继续在每组内调用现有六矩非负压缩，不改变粒子能量
  或做谱后处理。
- `src/tail_population_controller.cpp`：本轮不启用控制器，但其分组函数必须
  复用同一 `partition.energy_bin()`，禁止另建不一致能量边界。
- checkpoint/snapshot manifest 写出边界数量、边界值和边界哈希；restart
  必须校验，不允许用旧分组配置静默续跑。
- 最后一条显式能量边界必须 clamp 到实际 conversion 网格最大能量；当前 manifest
  中约 13.6 MeV 的末边界超过实际可表示上限，需在构造和端点测试中消除该
  配置/支撑不一致。

若 C 仍因每组 7 个 support 产生可见谱失真，应先增加“每个阈值细 bin 的
$N$ 与 $K$ 不跨组”这一结构约束，而不是增加全局矩补丁。只有专项性能测试证明
资源可接受时，才允许把阈值附近组的 `max_support` 从 7 局部提高到 8--10；高能
远尾保持原上限。

**当前分支判定（2026-08-08）**：旧均匀分组确实属于分支 B，阈值感知边界已
修复粗组跨 bin 搬运；但真实生产的 6.0--6.2 MeV 空洞还同时满足分支 A 的
近轴支撑不足特征。因此不得继续把问题仅归因于 7-support 压缩。下一步必须先
完成全局 conversion source 和逐环占据审计：

- 若 `pre_extraction_bulk` 在 6.0--6.2 MeV 已为零，而相邻区非零，则设计
  $K_{out}$ 附近的**局部 $u_\parallel$ 加密**或守恒 cut-cell/subcell 表示；
  不全局增加 `Nv/Nmu`，不改变热核与远尾网格；
- 若提取前非零而 `created_tail` 为零，则修复 `is_conversion_cell()`、
  `extract_conversion_masses()`、能量分组或 quartet 压缩；
- 若创建后非零而接受态 combined 变零，则审计转换后的 drift/kick、MPI 迁移、
  碰撞和快照沉积顺序。

#### 7.11.5 诊断同步

修改 `BulkTailConversionDiagnostics`、`VpfpStepResult` 和 `VpfpDiagnostics`：

- conversion source spectrum 必须携带显式 `bin_edges`，不得只写 64 个无边界
  数值；至少同时支持 0.05、0.1、0.2 MeV 三档固定阈值窗口；
- 每个 rank 先形成 local source，再由 `MPI_Allreduce` 求 global source，只有
  global source 可以交给 rank 0 写文件；禁止 rank 0 直接保存自己的 local
  conversion 结果；
- 同一接受步若发生多个 conversion location，必须逐位置记录并另行累计 accepted
  step 总源，禁止后一次覆盖前一次；
- 对每个显式 bin 同时输出 `pre_extraction_bulk_N/K`、`removed_bulk_N/K`、
  `created_tail_N/K` 和 `accepted_tail_N/K`，并写明
  `conversion_location`、`conversion_source_stage`、step 和 time；
- `tail_threshold_interface_rank*.dat` 同时输出 conversion 前、黄金参考（仅
  专项测试）、最终 tail 和 combined 的相同 bin 谱；
- MPI 汇总后再计算全局 $dN/dK$，不能对各 rank 的相对量取平均；
- 报告每个 energy group 的输入 support、输出 support、fallback 和六矩残差；
- 这些详细量只在独立测试或 `diagnostic-level=2` 的快照计算，不进入每步生产
  热路径。

#### 7.11.6 独立验收标准

专项修复必须同时满足：

1. A 路黄金参考相对转换前离散参考的六矩残差均不超过 $10^{-10}$，空间密度
   残差保持现有舍入误差级；
2. C 路在每个有非零参考质量的显式 edge group 内，$N$ 与 $K$ 相对 A 路误差
   不超过 $10^{-10}$；这只验收粗组守恒。0.05/0.1/0.2 MeV 子 bin 的谱
   $L_1/L_\infty$ 必须单独报告并随加密收敛，不得用 edge-group $L_1$ 替代；
3. 1/2/5 rank 的全局六矩和逐 bin 谱差异不超过确定性求和误差，建议门槛
   $10^{-12}$；
4. 不允许出现负权重、NaN/Inf、conversion bulk 残留或事务失败后状态污染；
5. 新分组下创建粒子数不得超过当前生产压缩的 2 倍，转换器单次 wall time
   不得超过 1.5 倍；超过时先优化分组和工作区，不能删掉谱保真约束；
6. `near-axis-narrow` 在 1/2/5 rank 下必须通过六矩、逐 bin 谱和确定性一致性
   门；生产 conversion source 必须满足逐 bin
   `pre_extraction -> removed -> created -> accepted` 可追踪；
7. 只有独立测试通过后才跑 Beam 12 fs；12 fs 的
   `unexplained_growth_steps` 必须仍为 0；
8. 最终 40 fs controller-off 中，阈值三 bin 不得再出现“中间 bin 比两侧同时
   低两个数量级”的孤立人工谷。对同一个checkpoint做0.1/0.2/0.4 MeV后处理
   重新分箱，只验收窗口积分守恒和人工谷是否存在，不得把均匀重构后的逐细bin
   $L_1$当成生产离散收敛。只有采用不同生产转换离散参数得到的独立状态之间，才要求
   阈值$\pm0.4$ MeV区域归一化$L_1$差不超过5%、各非稀疏bin差不超过10%。

第 8 项只验收生产积分结果，不能替代前七项。禁止把“更接近 EPOCH”作为接口
通过条件。

#### 7.11.7 下一步详细实施顺序

按以下顺序执行；任一步失败即停止，不得跳到长跑：

1. **修复诊断数据流**：在 `src/bulk_tail_converter.*` 中让转换诊断携带显式
   bin edges、位置和四阶段谱；在 `src/vpfp_integrator.*` 中累计同一接受步的
   多次转换源；在 `src/vpfp_diagnostics.*` 中对谱做 MPI 全局求和后由 rank 0
   写出。该修改只增加 `diagnostic-level=2` 的审计，不改变生产转换结果。
2. **强化支撑审计**：修改
   `tests/hybrid_threshold_support_audit_test.cpp`，增加逐 $u_\perp$ 环、近轴
   `k=0..7`、阈值上下最近中心和占据加权统计。测试直接调用生产网格对象。
3. **强化 A/B/C 测试**：修改两个 `bulk_tail_threshold_interface*` 测试，加入
   `near-axis-narrow`，删除非 golden 路的无条件 PASS，结果改为截断写入，并
   分开报告 edge-group 守恒和细 bin 谱保真。
4. **恢复阶段门纪律**：在上述测试完成前，`THRESHOLD_AWARE_COMPRESSION` 只作
   candidate。生产运行必须显式记录策略；不得根据旧的、不完整的 PASS 记录
   宣称它已验收。
5. **运行独立门**：清空本轮测试目录，依次运行单 rank 0.05/0.1/0.2 MeV、
   原转换器回归、MPI 1/2/5 rank 的 golden 与 threshold-aware，以及近轴案例。
6. **根据源谱选择唯一修复分支**：输入已空则做局部网格/cut-cell 设计；输入
   非空而创建后为空则修转换器；创建后非空而接受态为空则修推进/迁移/沉积。
   在证据出现前不改碰撞、人口控制器或全局速度网格。
7. **只跑一次 12 fs 阻断回归**：在独立门全部通过后，以
   `diagnostic-level=2` 保存首次转换前后和 12 fs 快照。验收源谱链、Gauss、
   转换守恒、人口增长和阈值谱。
8. **最后才决定 40 fs**：12 fs 的 6.0--6.2 MeV 空洞得到物理/离散解释并通过
   设定门槛后，才允许 controller-off 40 fs。120 fs 仍不允许启动。

#### 7.11.8 源谱定位结论 [已完成]

§7.11.7 第 1--7 项已经由 `h9_threshold_interface_v2` 和
`hybrid_h9_threshold_source_beam_12fs` 完成。结论不再是三选一的
“尚未定位”，而是：

```text
pre-extraction bulk source  仅在 6.215--6.270 MeV 非零
removed bulk source         与 pre-extraction 逐 bin 相同
created tail source         与 removed 逐 bin 相同
accepted tail total         与累计 created 数目相同
6.0--6.2 MeV empty bin      首次转换前已经存在
```

因此关闭以下错误方向：

- 不继续调整 threshold-aware 压缩组或 `max_support` 解决硬空洞；
- 不修改 Sentoku--Kemp、SDE、碰撞 pair 或人口控制器；
- 不通过绘图平滑、扩大 bin 或向空 bin 补粒子；
- 不全局增加 `Nv/Nmu`；
- 不对跨阈值 cell 每步重复执行固定几何比例抽取。

#### 7.11.9 第四步：建立有限体积单元能谱参考

新增只读生产模块：

```text
src/tail_subcell_quadrature.h
src/tail_subcell_quadrature.cpp
tests/bulk_tail_cell_volume_spectrum_test.cpp
```

该模块接收 `CylindricalVelocityGrid` 的真实面坐标
`upar_faces[j:j+1]`、`uperp_faces[k:k+1]` 和 cell-integrated mass，不得把
cell mass 再放到中心能量。第一版使用确定性的张量 Gauss 求积；$u_\perp$ 方向
权重必须包含圆柱测度 $2\pi u_\perp\,du_\perp$。每个 cell 的子点权重先归一，
必须严格满足：

$$
\sum_q w_q=M_{jk},\qquad w_q\ge0.
$$

对能谱 bin $[K_a,K_b)$，直接累计满足
$K(u_{\parallel,q},u_{\perp,q})\in[K_a,K_b)$ 的子点。至少输出：

```text
cell_center_histogram_N
cell_volume_histogram_N
cell_volume_histogram_K
cell_volume_N_residual
cell_volume_K_residual
straddling_cell_count
straddling_cell_mass
```

独立测试必须包括最低八个 $u_\perp$ 环、5.809/6.228 MeV 相邻中心及
0.05/0.1/0.2 MeV 三档 bin。测试目标是回答：虽然没有 6.0--6.2 MeV 的 cell
center，相关有限体积 cell 的几何范围和非负单元内表示是否覆盖该区。

第一版只用 cell-constant 重构建立基线。若需要斜率，只允许使用与生产 remap
同源的单调受限重构；禁止为得到平滑图使用无约束高阶插值。单元体积参考仅用于
诊断，不改变 `Species::f`、转换质量或场推进。

#### 7.11.10 第五步：单元内守恒 PIC 加载候选 [失败，停止接入]

只有 §7.11.9 证明有限体积参考在 6.0--6.2 MeV 非零时，才修改
`src/bulk_tail_converter.*`。转换时仍提取当前生产判定的完整 conversion cell
mass，不改变提取时间层和事务语义；改变的仅是该 mass 在 PIC 速度空间中的
加载方式。

具体实现：

1. 对每个被提取的 $(x,j,k)$ cell，调用 `TailSubcellQuadrature` 在真实速度面内
   构造确定性候选节点；候选节点必须满足 $K\ge K_{out}$，并包含足够的
   $u_\parallel$ 和 $u_\perp$ 支撑；
2. 每个非零 $u_\perp$ 节点仍生成方位对称 quartet，确保横向一阶矩严格抵消；
3. 使用现有 `TailMomentConstraint` 对候选权重做非负约束求解，目标矩直接取本次
   被扣除 bulk mass 的
   $N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$，不得先生成谱再做全局矩补丁；
4. 所有权重必须非负。若局部约束不可行，保留当前 golden quartet 作为显式
   fallback，并记录 cell、原因、目标矩和 fallback 粒子数；禁止接受负权重；
5. 在阈值窗口内按 0.05/0.1/0.2 MeV 细 bin 报告有限体积参考与加载后 PIC 谱，
   同时保留 edge-group 六矩账；
6. 相同 cell、step、rank 和配置必须产生逐位确定的候选节点、权重和粒子 ID；
   不引入随机数；
7. 预计算每种 $(j,k)$ 的几何节点、能量和矩列，生产转换时只缩放/求权重，禁止
   每个 $x$ cell 重建几何工作区。

fallback 只用于保证失败事务可诊断，不代表新方案通过。若 12 fs 中阈值相关 cell
出现任何 fallback，或独立测试中 fallback 比例非零，则单元内加载门判定失败并
转入 §7.11.11；不得用 golden fallback 混合出一条表面平滑的生产谱。

**[2026-08-08 第10步首次失败与三路可行性审计]**：首次
`near-axis-narrow` 结果为 `subcell_cells_loaded=0`、`fallback reason=2`，六矩
残差为0只来自中心quartet fallback；0.05/0.1/0.2 MeV细谱误差分别约
2.21、2.21、1.42，不能进入后续回归。为区分求解器错误与离散目标不兼容，新增
`tests/bulk_tail_subcell_feasibility_test.cpp`，必须分别执行：

对应的历史测试命令已统一迁移到 §17.10.1 的“subcell 可行性审计 [已完成、路线已停止]”；
本节不再保留命令副本。

三项输出的含义固定为：

1. `self-moments` 必须通过，否则先修非负矩求解器；
2. `center-augmented_center_seed` 必须通过；若展开先验失败或最终非中心权重趋近0，
   说明中心目标只能退化回中心支撑；
3. `constraint-hierarchy` 依次加入 $N,P_x,K,J_x,\Pi_{xx},\Pi_\perp$，首次失败
   层就是当前候选表示与中心目标的冲突起点。

本地MPI桩和集群正式结果一致：self-moments和精确中心控制均以零残差通过；仅约束
$N,P_x$ 时subcell候选可行，加入 $K$ 后首次失败，完整六矩残差约
$8.11\times10^{-5}$。这支持“cell-center目标与展开的非负subcell表示不兼容”，
而非求解器完全失效。第10步正式判定失败，生产默认必须保持
`subcell_loading_enabled=false`；不得提高迭代次数、放宽六矩门或依靠fallback
接入生产。

该方案的目的不是凭空生成连续谱，而是把有限体积 cell 已代表的子单元信息以
守恒的 PIC 支撑表达出来。若六矩约束要求部分节点高于 6.228 MeV，同时允许部分
节点位于 6.0--6.2 MeV，这是正常的非负矩配对，不是人工补 bin。

#### 7.11.11 候选失败时的备选分支

若 §7.11.9 的有限体积参考在 6.0--6.2 MeV 仍为零，或 §7.11.10 的非负六矩
问题在真实近轴 cell 中系统性不可行，不继续增加求积节点碰运气。此时比较两个
独立方案：

1. **固定总 `Nv` 的局部 $u_\parallel$ 重分配**：在正负阈值速度附近对称增加
   face/center，将同等数量的极远尾低占据 cell 移出；保持总内存不变。必须重跑
   remap、Maxwellian 矩、碰撞平衡和阈值支撑测试；
2. **通量式表示转换**：在 Vlasov 速度通量实际跨越 resolved/tail 接口时，按
   该面通量创建 tail，并从同一通量更新 bulk。它最符合连续性，但会改变生产
   remap/转换耦合，只在局部重分配仍不收敛时实施。

禁止的“静态 cut-cell”是：每步用同一几何比例从一个尚未跨界的剩余 cell mass
中继续抽取。若以后实现 cut-cell，必须保存子单元状态或由真实跨界通量驱动，
使同一份质量不会被重复转换。

#### 7.11.12 已完成路线的结论

1. cell-center谱在6.0--6.2 MeV为空，但有限体积几何参考非零；
2. 非负矩求解器能精确恢复subcell自身矩；
3. 中心目标只允许精确中心支撑；展开先验从加入 $K$ 开始失败；
4. 因此“保持旧中心六矩不变，同时展开成细谱”的要求本身不兼容；
5. 旧第11、12步暂停，下一步转入真实转换质量的矩表示差审计。

#### 7.11.13 第六步：真实转换事件的矩表示差审计 [已完成]

该审计不能从10.5 fs已接受checkpoint直接计算，因为该状态中的conversion cells
已经被converter清空。必须在 `BulkTailConverter::extract_after_substep()` 完成
正质量请求扫描之后、修改 `bulk_trial.f` 之前计算；结果随
`BulkTailConversionDiagnostics` 返回，只有完整VPFP步被接受后才全局归约和落盘。
拒绝的trial不得进入审计文件。

新增文件：

```text
src/bulk_tail_moment_audit.h
src/bulk_tail_moment_audit.cpp
tests/bulk_tail_real_moment_audit_test.cpp
tests/bulk_tail_real_moment_audit_mpi_test.cpp
```

现有文件的修改位置固定为：

1. `src/bulk_tail_converter.h`：在 `BulkTailConversionDiagnostics` 中增加
   `MomentRepresentationAudit`，包含六矩的center/full-volume/eligible-raw/
   eligible-normalized、L1增量、signed增量、最大cell误差、计数和top-cell；增加
   `set_moment_audit_enabled(bool)`，默认关闭；
2. `src/bulk_tail_converter.cpp`：在 `requests` 完整形成后、任何
   `bulk_trial.f=0` 之前调用只读审计。不得改变request、group、粒子创建或fallback
   决策；
3. `src/vpfp_integrator.*`：把trial诊断保存在step result中。仅当
   `state_advanced=1` 时把该诊断交给全局化和输出；失败/split/retry trial全部丢弃；
4. `src/vpfp_diagnostics.*`：增加
   `write_bulk_tail_moment_audit_accepted_step()`。rank-local top-cell先携带global
   x index和rank，rank 0合并排序；固定长度sum/max数组合并为一次Allreduce；
5. `src/main_vpfp.cpp`：解析两个新CLI，将开关传给integrator/converter并写入
   manifest。审计开启且subcell加载开启时，在step 0前失败；
6. `CMakeLists.txt`：加入两个独立测试target及 `ctest` 注册；测试不得依赖12 fs
   生产输出。

实现时优先复用 `mass_cell_moments()`、`TailSubcellQuadrature` 和现有conversion
globalization协议。禁止在diagnostics中复制第四套矩公式，也禁止为审计重新扫描
完整 $(x,u_\parallel,u_\perp)$ 数组；只遍历本次已有的正质量conversion requests。

新增CLI：

```text
--tail-cell-moment-audit
--tail-cell-moment-audit-top-cells <N>   # 默认64
```

不开启CLI时不得执行子点扫描、额外MPI归约或文件输出。manifest必须写出
`tail_cell_moment_audit=0|1` 和 `tail_subcell_loading=0|1`；审计模式若发现
subcell加载开启，立即报错，避免审计对象被改变。

对每个真实请求cell $(x,j,k,M)$ 同时构造三套矩：

1. **中心矩 $C$**：继续调用当前生产 `mass_cell_moments(M,u_j,u_k)`；
2. **完整体积矩 $V$**：使用 `TailSubcellQuadrature` 的全部非负节点和原始归一
   权重；
3. **tail可表示矩 $R$**：只保留 $K_q\ge K_{out}$ 的节点，并分别输出未归一
   `R_raw` 和重新归一到总数 $M$ 的 `R_norm`。

六个分量固定为：

$$
N,\quad P_x,\quad J_x,\quad K,\quad\Pi_{xx},\quad\Pi_\perp.
$$

每个分量 $X$ 必须累计：

$$
R_{L1}^{X}=\frac{\sum_c|X_{V,c}-X_{C,c}|}
{\max(\sum_c|X_{C,c}|,X_{floor})},
$$

$$
R_{signed}^{X}=\frac{|\sum_c(X_{V,c}-X_{C,c})|}
{\max(\sum_c|X_{C,c}|,X_{floor})}.
$$

不得用全局净 $P_x/J_x$ 作唯一分母，因为正负速度支撑会相互抵消。还必须输出：

```text
accepted_step time_fs conversion_location
request_cell_count positive_request_cell_count
center_* volume_* eligible_raw_* eligible_norm_*
delta_* rel_l1_* rel_signed_* max_cell_rel_*
eligible_number_fraction below_threshold_number_fraction
threshold_window_6p0_6p2_N
center_target_feasible_count center_target_failed_count
volume_self_feasible_count volume_self_failed_count
eligible_target_feasible_count eligible_target_failed_count
```

rank 0写 `bulk_tail_moment_audit_accepted_steps.dat`；另外只输出按
`|delta_K|+|delta_Jx|+|delta_Pixx|+|delta_Piperp|` 排序最大的64个cell。所有sum
使用 `long double` 本地累加，MPI只对最终固定长度数组做一次 `MPI_Allreduce`。

#### 7.11.14 审计实现的独立验收 [已完成]

`bulk_tail_real_moment_audit_test` 至少覆盖：

1. `single-cell-known`：中心/体积/eligible三套矩与直接求积一致；
2. `symmetric-pair`：正负 $u_\parallel$ 的signed动量抵消，但L1尺度非零；
3. `threshold-face`：阈值位于公共面时允许 `straddling=0`，eligible质量仍正确；
4. `clipped-cell`：存在阈值以下子点时，`R_raw`质量亏损和`R_norm`重归一可见；
5. `accepted-only`：失败trial不写文件，接受步恰好写一次。

MPI测试必须用1/2/5 rank产生独立文件，并验证全局sum、L1、max及top-cell归属与
单rank结果一致。任何测试都不得推进生产物理状态。

**实现记录（2026-08-08）：**

- 两个测试目标均已实现 `--case` 和 `--result` 参数；结果文件无法创建时以退出码3
  明确失败，不再静默忽略路径；
- 单rank测试覆盖 `single-cell-known`、`symmetric-pair`、`threshold-face`、
  `clipped-cell` 和 accepted-only 输出门，并使用非单位质量暴露量纲错误；
- MPI测试使用20个全局x cell和固定的4个全局转换请求。每个请求只由拥有其
  `ix_global` 的rank处理，因此1/2/5 rank表示同一物理输入；验收全局七组矩数组、
  六组计数、逐分量max、每个top-cell的唯一owner及全局排序；
- 已修复审计内核中被 `mass=1` 掩盖的二次乘质量错误：`volume` 已是积分矩，
  正确差值为 `volume - mass*center_unit`；可行性目标也统一为积分中心矩；
- 生产accepted-step写出增加共享门函数，拒绝态或禁用态不能进入审计文件。

本地MPI stub下全部单rank案例通过；集群1/2/5 rank固定全局输入测试也已全部
`status=PASS`。旧的 `global_N=1/2` 假阳性结果已废弃，不得用于验收。

#### 7.11.15 真实结果的决策门

先检查硬正确性：$N$闭合不超过 $10^{-12}$；无NaN/Inf；MPI 1/2/5 rank一致；
accepted-only语义正确。随后按真实转换质量加权结果分类：

**绿色：可继续评估统一有限体积矩。** 同时满足：

- `below_threshold_number_fraction <= 1e-6`；
- $J_x,K,\Pi_{xx},\Pi_\perp$ 的事件累计 $R_{L1}\le1e-3$；
- 累计 $|\Delta K|$ 不超过累计removed-tail能量的 $10^{-3}$；
- 不存在占总转换质量 $10^{-6}$ 以上、单cell相对误差大于 $10^{-2}$ 的区域；
- volume-target非负可行性失败数为0。

绿色只允许进入“统一矩定义”的独立原型，不允许直接把 $V$ 塞进converter。原型
必须预计算cell-volume的 $v_x/K/\Pi_{xx}/\Pi_\perp$ 权重，并同步用于
`Species::compute_moments()`、背景能量/压力诊断、碰撞矩和converter目标；密度权重
保持严格为1。不得只改converter制造隐蔽的场功或能量跳变。

**红色：拒绝统一有限体积矩。** 任一条件成立即进入固定总 $N_v$ 的局部
$u_\parallel$ 网格重分布原型：

- 任一主要矩的事件累计 $R_{L1}\ge1e-2$；
- `below_threshold_number_fraction >= 1e-4`；
- 累计 $|\Delta K|$ 达removed-tail能量的1%；
- volume-target在真实占据cell中系统性不可行。

**灰色：$10^{-3}<R_{L1}<10^{-2}$。** 不修改生产算子；先对固定总 $N_v$ 的局部
重分布做只读几何A/B，比较阈值中心间距、Maxwellian矩误差和预计内存/时间。只有
局部重分布不能把上述误差压到绿色门内时，才设计通量式representation
conversion。

#### 7.11.16 真实12 fs结论与红色分支实施方案

##### 7.11.16A 真实结果 [已完成]

目录 `output/hybrid_h9_real_moment_audit_12fs` 使用80 rank运行470个接受步并到达
12 fs；无split、拒绝步或NaN/Inf。manifest确认：

```text
ranks=80
field_boundary=dirichlet-phi
collision_model=moment-closure
tail_cell_moment_audit=1
tail_subcell_loading=0
```

首次真实转换为step 409、10.464300709 fs。到12 fs共有54个含转换的接受步、743个
正质量转换cell，累计转换数目为
$5.068692683154256\times10^{10}$。全局结果为：

| 分量 | 累计 $R_{L1}$ | 最大已记录cell真实相对差 | 结论 |
|---|---:|---:|---|
| $N$ | $1.10\times10^{-16}$ | $5.90\times10^{-16}$ | 通过 |
| $P_x$ | $1.81\times10^{-16}$ | $4.76\times10^{-16}$ | 通过 |
| $J_x$ | $3.00\times10^{-6}$ | $3.00\times10^{-6}$ | 绿色 |
| $K$ | $1.10\times10^{-6}$ | $1.10\times10^{-6}$ | 绿色 |
| $\Pi_{xx}$ | $1.00\times10^{-6}$ | $1.00\times10^{-6}$ | 绿色 |
| $\Pi_\perp$ | $3.2617\times10^{-2}$ | $1.1275\times10^{-1}$ | **红色** |

累计阈值以下数目分数为 $1.11\times10^{-16}$，所以错误不是由低于6 MeV的质量
越界转换造成。转换事务的最大相对残差为
$R_N=7.19\times10^{-16}$、$R_{P_x}=4.04\times10^{-27}$、
$R_K=1.52\times10^{-18}$，说明提取和创建账本闭合。

`volume_target_feasible_count=0`、`volume_target_failed_count=743`，eligible目标同样
0/743。失败集中在 `iv=185`、低 `imu=1..4`；这是cell-center $u_\perp$ 与包含
$2\pi u_\perp\,du_\perp$ 权重的cell-volume矩不一致。该证据同时否决当前subcell
展开和“只改成有限体积矩”两条捷径。

##### 7.11.16B 先修审计显示，不改物理推进 [下一步第1项]

当前主表有470行，其中416行是无真实转换的全零占位；另外 `rel_l1_*` 和
`max_cell_rel_*` 使用了有量纲的 `max(1,scale)`，当 $J/K/\Pi<1$ 时实际输出绝对
误差而不是相对误差。必须先做以下小修：

1. 修改 `src/vpfp_diagnostics.cpp`：对每个event先全局归约
   `positive_request_cell_count`；全局值为0时所有rank一致跳过后续归约和文件写入；
2. 修改 `src/bulk_tail_moment_audit.cpp`：每个cell的相对误差使用
   $|\Delta X|/|X_C|$；若二者均为0则记0，若 $X_C=0$ 但 $\Delta X\ne0$ 则写
   `inf`并将该cell标记为不可比较，禁止使用有量纲常数1作分母；
3. 修改 `src/vpfp_diagnostics.cpp`：事件累计分母直接使用归约后的
   `center_l1_*`。分母为0且分子为0写0；分母为0且分子非零写 `inf`；
4. 在主表中显式写出 `center_l1_*` 和 `delta_l1_*` 原始量，避免只能从
   `rel_l1_*` 反推；
5. 扩展 `bulk_tail_real_moment_audit_test` 的 `accepted-only`：启用但零请求的
   接受event不得写行；非零接受event恰好写一行；
6. 扩展MPI测试：1/2/5 rank下空event都不写，非空event的原始L1和相对值一致；
7. 增加 `bulk_tail_moment_audit_velocity_histogram.dat`：只对真实接受event按
   `(iv,imu)` 写全局request质量和cell计数，作为后续离线回放输入。该文件不得包含
   x数组、trial数据或零质量记录。

**16A必须按以下接口实施，禁止自由发挥：**

1. 在 `src/bulk_tail_moment_audit.h` 新增：

   ```cpp
   struct BulkTailVelocityBinAudit {
       int iv;
       int imu;
       std::uint64_t request_cell_count;
       double request_number;
   };
   ```

   并在 `MomentRepresentationAudit` 中加入
   `std::vector<BulkTailVelocityBinAudit> velocity_bins`。这里的
   `request_number` 是request质量之和，不是PIC宏粒子数。

2. `bulk_tail_audit_conversion_requests()` 只遍历已有request。先按
   `key=iv*Param::Nmu+imu` 排序或映射聚合，再按key升序写入 `velocity_bins`。
   不得扫描完整背景分布，不得按x展开输出。
3. 将全局归约从 `VpfpDiagnostics` 中拆成可独立测试的生产函数，建议放入：

   ```text
   src/bulk_tail_moment_audit_io.h
   src/bulk_tail_moment_audit_io.cpp
   ```

   接口至少返回 `has_positive_requests`、全局sum/L1/max/count和全局速度bin。
   `VpfpDiagnostics::write_bulk_tail_moment_audit_accepted_step()` 必须直接调用该函数；
   测试也调用同一函数，禁止在测试中复写MPI公式。
4. 每个event的MPI顺序固定为：
   - 在进入逐event循环前，对本rank的 `conversion_events.size()` 分别做 `MIN/MAX`
     归约；二者不等时返回 `INVALID_EVENT_LAYOUT`，所有rank走同一失败路径，禁止部分
     rank继续进入collective；
   - 对每个event的 `conversion_location` 也做 `MIN/MAX` 一致性检查；不一致时同样
     返回 `INVALID_EVENT_LAYOUT`。不能假设各rank恰好具有相同event布局；
   - `MPI_Allreduce(local_positive_count, SUM)`；
   - 若全局为0，所有rank立即 `continue`，不得进入后续collective；
   - 一次固定长度 `MPI_Allreduce(..., SUM)` 合并原始矩和L1；
   - 一次 `MPI_Allreduce(..., MAX)` 合并逐分量cell最大值；
   - 速度bin可使用长度 `Param::Nvmu` 的number/count工作数组做两次归约，审计关闭时
     不得分配；
   - rank 0只输出全局count大于0的bin，并按 `(iv,imu)` 排序；
   - `top_cell_limit` 必须作为审计配置随event传入。每个rank可先保留本地前
     `top_cell_limit` 项以限制通信量（这足以构造精确的全局top-N），随后rank 0汇集
     这些局部候选，按确定的 `(relative_error desc, global_ix, iv, imu)` 排序并再次截断
     为全局前 `top_cell_limit` 项；禁止用“最大局部条数”替代第二次全局截断。
     1/2/5 rank测试必须逐项比较最终top-cell的 `global_ix/iv/imu` 和原始矩。
5. 相对量统一使用：

   $$
   r(a,b)=\begin{cases}
   a/b,&b>0,\\
   0,&a=0,\ b=0,\\
   +\infty,&a>0,\ b=0.
   \end{cases}
   $$

   同时输出 `relative_defined_0..5`。不得用 `DBL_MIN`、`1` 或其他有量纲常数
   偷换分母；`inf` 只表示诊断比例未定义，不得使物理步失败。
6. 主表每个分量的列顺序固定为：

   ```text
   center_X volume_X eligible_raw_X eligible_normalized_X
   delta_signed_X center_l1_X delta_l1_X
   rel_l1_X rel_signed_X max_cell_rel_X relative_defined_X
   ```

7. 速度直方图列固定为：

   ```text
   accepted_step time_fs conversion_location iv imu
   request_cell_count request_number
   ```

   同一event、同一 `(iv,imu)` 只能有一行。
8. 新增 `tools/analyze_bulk_tail_moment_audit.py`，输入主表和速度直方图，输出
   `audit_summary.result`。脚本必须从原始 `center_l1_X/delta_l1_X` 求累计比例，禁止
   对逐步比例做算术平均。输出54个真实event、743个cell、累计六矩、最大cell比例、
   阈值以下比例、可行性计数和 `decision=GREEN|GRAY|RED|INVALID`。脚本遇到
   `relative_defined_X=0` 时只增加 `undefined_relative_count_X`，不得把数学上未定义的
   比例本身判为物理失败；event布局不一致、缺列、重复event/bin或原始量非有限才判
   `INVALID`。
9. `tests/bulk_tail_real_moment_audit_test.cpp` 必须新增非单位质量和零分母案例；
   `tests/bulk_tail_real_moment_audit_mpi_test.cpp` 必须让固定物理请求跨rank重新分区，
   不得让每个rank各自新增一份相同请求。
10. MPI测试必须再加入故意构造的event条数不一致与location不一致案例，验收所有rank
    均返回 `INVALID_EVENT_LAYOUT` 且不死锁；加入每rank局部候选数不同的top-64案例，
    验收1/2/5 rank得到完全相同的全局top集合。

该项只允许修改审计数据结构、归约和写出，不得改 `bulk_trial.f`、conversion
request、PIC创建、碰撞或Poisson。完成后从首次转换前checkpoint重跑10.4--12 fs
即可，不需要从0 fs重跑；红色物理结论不依赖这次复跑才能成立。

**16A 实施记录（2026-08-08，本机 MPI 桩 + g++ -O2）**：

- 新增 `src/bulk_tail_moment_audit_io.h/.cpp`：可独立测试的生产归约
  `bulk_tail_moment_audit_check_event_layout`（event 数 MIN/MAX →
  INVALID_EVENT_LAYOUT）与 `bulk_tail_moment_audit_reduce_event`
  （location MIN/MAX → 正请求 SUM → 原始矩/L1 一次 SUM → 逐分量 cell-max
  一次 MAX → 速度 bin 两次 Nvmu 归约 → top-cell 全局二次截断）。
  `VpfpDiagnostics::write_bulk_tail_moment_audit_accepted_step` 与两个测试
  均直接调用该函数，禁止测试复写 MPI 公式。
- `bulk_tail_moment_audit.h/.cpp`：新增 `BulkTailVelocityBinAudit` 与
  `velocity_bins`（按 `iv*Nmu+imu` 聚合、仅正质量、key 升序）、事件级
  `relative_defined`；`max_cell_relative`/top-cell score 改用
  r(a,b)=a/b | 0 | +inf 规则，删除有量纲常数 1 分母；inf 只表示诊断比例
  未定义，不使 `finite` 失效。
- `vpfp_diagnostics.cpp`：主表按 §7.11.16B 第 6 项列序写出（含
  `center_l1_X`/`delta_l1_X` 原始量与 `relative_defined_X`）；全局正请求为
  0 的 event 不写行（去掉全零占位）；新增
  `bulk_tail_moment_audit_velocity_histogram.dat`（仅真实接受 event、
  `(iv,imu)` 一行、无零质量记录）；manifest 增加
  `tail_cell_moment_audit_top_cells`；审计关闭时零集体零分配。
- 测试：单 rank 新增 `zero-denominator` 案例并把 `accepted-only` 改为走 IO
  归约（空 event 不写行）；MPI 测试改用共享归约函数、固定物理请求跨 rank
  分区、新增 event 数与 location 不一致的 INVALID_EVENT_LAYOUT 案例（全体
  rank 同一失败路径、不死锁）与局部候选数不同的 top-64 一致性校验。
- 新增 `tools/analyze_bulk_tail_moment_audit.py`：从原始
  `center_l1/delta_l1` 求累计比例（不做逐步比例算术平均），输出真实
  event/cell 数、累计六矩、最大 cell 比例、阈值以下分数、可行性计数与
  `decision=GREEN|GRAY|RED|INVALID`；`relative_defined_X=0` 仅计入
  `undefined_relative_count_X`；缺列/重复 event 或 bin/非有限原始量判
  INVALID。阈值：GREEN≤1e-6、RED≥1e-2、GRAY 居中（对应 §7.11.16A 表格
  Jx/K/Pixx≈1e-6 绿色、Piperp 3.26e-2 红色）。
- 本机验证：两个审计测试 `status=PASS`（含新案例）；single_cell/multibin/
  poisson/transaction 回归 PASS；fp_solver 链语法通过；脚本合成数据验证
  （0.03 来自原始 3/100、RED）。集群 1/2/5 rank 与 10.4--12 fs 复跑待按
  §7.11 第 16A 步执行。
- 集群首跑修正（2026-08-08）：`bulk_tail_real_moment_audit_mpi_test` n=2/5
  FAIL 的根因是 `velocity_bins`/`top_cells` 只在 rank 0 由归约填充，而测试
  在全部 rank 上与 reference 比较（非 rank-0 的全局列为空）。已修复：这两
  个 rank-0 专属字段只在 rank 0 比较；归约本身（所有 rank 参与 collective）
  不变。本机 n=1 PASS。
- 16A 集群验收（2026-08-08，0--12 fs 全程重跑，非 checkpoint 续跑）：
  四档测试全 `status=PASS`；生产审计 470 步全 accepted，首次转换 step 409
  /10.464300709 fs、末步 470/12 fs 与原 §7.11.16A 记录一致；真实事件 54
  一致；`zero_event_rows=0`、主表与速度直方图非空、列格式含
  `center_l1_5`/`delta_l1_5`、manifest `tail_cell_moment_audit=1`/
  `tail_subcell_loading=0`、脚本 `input_valid=1`；
  `R_L1^Piperp=3.33e-2`（≈旧 3.26e-2，决策 RED 保持，未因诊断修改变绿）、
  可行性 0/611。**cell 数 611 vs 原记录 743**：按 0 fs 重跑分支验收，已比对
  首次转换与末步（均一致），差异归因于本次二进制含 7.11 分支 B 及 H9
  碰撞修复（threshold-aware 分组、BGK 后端、转换边界零通量墙等）后
  bulk 状态变化，叠加文档已记录的束流-等离子体不稳定性 ULP 发散；不属
  浮点终止误差，且不改变红色物理结论。据此 16A 判通过。

##### 7.11.16C 固定内存的局部速度网格离线原型 [下一步第2项]

新增独立模块，禁止首先修改 `CylindricalVelocityGrid::init()` 的生产默认：

```text
src/tail_interface_grid_design.h
src/tail_interface_grid_design.cpp
tests/tail_interface_grid_design_test.cpp
tests/tail_interface_grid_replay_test.cpp
tests/tail_interface_grid_replay_mpi_test.cpp
```

头文件至少提供以下独立接口；命名可按现有风格调整，但职责不得合并进测试：

```cpp
struct TailInterfaceGridDesignConfig {
    double ax;
    double aperp;
    double sigma_x_cells;
    double sigma_perp_cells;
    double min_width_ratio;
    double max_adjacent_width_ratio;
};

struct TailInterfaceGridCandidate {
    std::vector<double> upar_faces;
    std::vector<double> uperp_faces;
    bool valid;
    std::string invalid_reason;
};

TailInterfaceGridCandidate build_tail_interface_grid_candidate(
    const CylindricalVelocityGrid& baseline,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    const TailInterfaceGridDesignConfig& config);

TailInterfaceReplayResult replay_tail_interface_histogram(
    const CylindricalVelocityGrid& baseline,
    const TailInterfaceGridCandidate& candidate,
    const std::vector<BulkTailVelocityBinAudit>& histogram,
    double conversion_energy);
```

测试必须调用这两个生产函数。禁止在测试文件中复写monitor积分、face反演或圆柱交叠
公式。`tail_interface_grid_design.cpp` 不能引用 `main_vpfp.cpp` 的全局CLI状态，也不能
修改 `Param::Nv/Param::Nmu`。

原型必须固定：

```text
Nu_parallel = 192
Nu_perp = 64
u_parallel_min/max 不变
u_perp_min/max 不变
总 phase-space cell 数不变
```

构造四个候选：`G0`为当前网格；`Gx`只移动阈值附近的 $u_\parallel$ faces；
`Gp`只移动低 $u_\perp$ faces；`G2`同时移动两组faces。任何候选均不得新增cell、
扩大速度域或改变空间网格。face生成必须满足严格单调、镜像对称的
$u_\parallel$、$u_\perp=0$固定、端点固定；相邻cell宽度比默认不超过2，最小宽度
默认不小于当前对应区域宽度的0.5倍。参数必须由CLI提供，不能硬编码到生产网格。

**候选face必须由确定性monitor equidistribution生成：**

$$
w_x(u)=1+A_x\exp\left[-\left(\frac{|u|-u_*}{\sigma_x}\right)^2\right],
$$

$$
w_\perp(u_\perp)=1+A_\perp\exp\left[-\left(\frac{u_\perp}{\sigma_\perp}\right)^2\right].
$$

$u_*$ 由真实速度直方图质量加权的低 $u_\perp$ 转换cell计算，不得写死为某个iv。
在原速度域上高分辨率积分monitor累计函数，然后按等累计值反求新faces。参数第一轮
只扫描：

```text
Ax, Aperp = 0.5, 1.0, 2.0
sigma_x = 1, 2, 4 个当前阈值附近du_parallel
sigma_perp = 1, 2, 4 个当前低u_perp cell宽度
```

生成后应用相邻宽度比和最小宽度约束；违反约束的候选直接标记
`grid_valid=0`，不能自动裁剪后假装通过。`Gx`令 $A_\perp=0$，`Gp`令 $A_x=0$，
`G2`使用二者组合。相同参数和输入必须逐bit生成相同face数组。

真实请求回放需要新增审计文件
`bulk_tail_moment_audit_velocity_histogram.dat`。它按接受步和 `(iv,imu)` 聚合真实
request的数目质量，不输出全x数组。每个rank先本地稀疏聚合，再进行确定性的全局
归并；空event不写。离线测试读取该文件，把旧cell常数分布通过二维cell交叠体积
保守重映射到候选网格。$u_\perp$ 交叠权重必须积分
$\int 2\pi u_\perp du_\perp$，禁止按几何宽度线性分配。

旧cell $(j,k)$ 向新cell $(j',k')$ 的质量份额固定为：

$$
\theta_{jk\rightarrow j'k'}=
\frac{\Delta u_{\parallel}^{\rm overlap}}
     {u_{\parallel,j+1/2}-u_{\parallel,j-1/2}}
\frac{u_{\perp,hi}^{2}-u_{\perp,lo}^{2}}
     {u_{\perp,k+1/2}^{2}-u_{\perp,k-1/2}^{2}}.
$$

无交叠时 $\theta=0$。每个旧cell必须满足
$|\sum_{j'k'}\theta-1|\le10^{-14}$；新质量为
$M_{j'k'}=\sum_{jk}\theta M_{jk}$。不允许在重映射后clip负值或做全局质量补丁。

离线原型必须有以下控制案例：

1. `G0 identity`：faces、质量和六矩必须与输入审计重构一致；若失败，后续候选全部
   `INVALID`；
2. 单cell解析交叠：手算份额与实现误差不超过 $10^{-14}$；
3. 常数cell分布：重映射后数目严格闭合；
4. 对称正负 $u_\parallel$：$P_x/J_x$ signed和保持抵消，L1尺度非零；
5. Maxwellian和漂移Maxwellian：旧到新再回旧，无负质量，数目误差不超过
   $10^{-12}$，六矩误差报告而非强制修正；
6. 真实直方图：1/2/5 rank输入分区后所有候选的faces和汇总指标一致。

每个候选必须输出：

```text
grid_name
number_residual
negative_mass
below_threshold_number_fraction
R_L1_N R_L1_Px R_L1_Jx R_L1_K R_L1_Pixx R_L1_Piperp
max_cell_relative_*
min_dupar min_duperp max_adjacent_width_ratio
estimated_cell_count_ratio estimated_memory_ratio
estimated_velocity_dt_ratio estimated_operator_work_ratio
center_target_feasible_count center_target_failed_count
volume_self_feasible_count volume_self_failed_count
volume_self_sparse_failed_count max_sparse_support_count
estimated_created_macroparticles estimated_particle_ratio_to_center_quartet
```

`center_target_feasible_count` 只作信息项，不能作为硬门。对任意非零宽度cell，要求
一组展开节点同时严格复现位于中心点的动能和二阶矩，通常会因凸性而不可行；继续
要求失败数为0会把原型设计成必然失败。真正硬门是 `volume_self_failed_count=0`，即
候选cell-volume矩必须能由其自身非负求积节点复现。

但“存在16节点非负体积求积解”仍不足以进入生产。当前center quartet每个转换cell只创建
4个宏粒子；若直接保留全部16个体积求积节点，tail人口和碰撞成本可能立即增至约4倍。
因此离线原型必须在保持权重非负的前提下，对volume-self节点执行确定性的稀疏支撑约简：

1. 目标仍为同一组六矩，不允许删除 $\Pi_\perp$ 或用全局能量补丁；
2. 使用固定pivot顺序的NNLS/Caratheodory型约简，只删除线性相关支撑，不改变目标矩；
3. 每cell最多保留7个非零节点。残差不得除以可能接近0的signed矩；每个分量使用
   `abs(residual_X)/max(l1_target_X, scale_floor_X)`，其中 `scale_floor_X` 由该cell数目
   乘对应速度/能量特征尺度构造并显式写入结果。六矩尺度化残差均不超过
   $10^{-10}$，数目残差不超过 $10^{-12}$；
4. 输出失败cell数、最大支撑数、按真实请求估算的新增宏粒子总数，以及相对center
   quartet的粒子数比例；
5. 若存在失败cell，或估算粒子数比例超过2.0，即使矩误差通过，也因资源门判失败。

固定cell数也不能视为固定运行成本。原型必须用生产算子的实际速度输运/受力约束，分别
在当前网格和候选网格上计算同一状态允许的最大步长，输出
`estimated_velocity_dt_ratio=dt_candidate/dt_current`；同时用固定cell扫描成本、预计PIC
创建数和碰撞pair/substep预算估计 `estimated_operator_work_ratio`。禁止仅用最小cell宽度
代替生产约束，也禁止运行完整12 fs后才发现候选使时间步大幅缩小。

离线门为：数目残差不超过 $10^{-12}$；无负质量；$J_x/K/\Pi_{xx}/\Pi_\perp$
全部 $R_{L1}\le10^{-3}$；重要cell最大相对差不超过 $10^{-2}$；阈值以下分数不超过
$10^{-6}$；volume-self失败数为0；cell数和预计内存比必须严格为1。由于当前主导
误差来自 $\Pi_\perp$，`Gx`即使通过其他矩但未修复 $\Pi_\perp$，仍判失败。
此外必须满足 `volume_self_sparse_failed_count=0`、`max_sparse_support_count<=7` 和
`estimated_particle_ratio_to_center_quartet<=2.0`。这些是性能/内存硬门，不得留到
12 fs生产A/B后才检查。还必须满足 `estimated_velocity_dt_ratio>=0.8` 和
`estimated_operator_work_ratio<=1.5`；估计模型及其各项贡献必须写入result，不能只写
最终比例。

##### 7.11.16D 接入顺序 [后续条件分支]

1. 若 `G2` 或更简单候选通过离线门，先以
   `--velocity-grid-profile current|tail-interface-v1` 非默认CLI接入；
2. 对Maxwellian、漂移Maxwellian和阈值附近窄分布做旧网格到新网格再回映的守恒
   回归，检查六矩、非负性和二阶平滑分布误差；
3. 旧10.4 fs checkpoint不能直接由候选网格读取。必须新增独立
   `vpfp_checkpoint_velocity_remap` 工具，读取旧checkpoint，对每个x cell的bulk
   分布执行与离线原型完全相同的二维圆柱体积保守重映射，保持Beam、tail、场、
   RNG和累计账本不变，并写入新的网格profile/hash；禁止在restart中隐式remap；
4. checkpoint转换工具必须先通过1/2/5 rank往返测试：$N$严格闭合、无负质量，
   六矩误差满足离线绿色门，旧checkpoint不被覆盖，转换后checkpoint可独立重启；
5. `current`分支使用原10.4 fs checkpoint，候选分支使用显式转换后的对应checkpoint，
   分别运行到12 fs做A/B；不能把两种网格指向同一个checkpoint目录；
6. A/B通过后才讨论修改生产默认；
7. 若所有固定cell候选失败，停止调网格，进入通量式representation conversion：
   在Eulerian速度面通量穿越 $K_{out}$ 时创建tail，不再把整个中心cell静态切走。

候选网格离线回放只证明表示误差和资源预算，不证明非线性动力学正确。即使全部离线门
通过，也必须完成本节的checkpoint显式转换和10.4--12 fs生产A/B；不得直接替换默认网格。

#### 7.11.17 通量式 representation conversion 开发规范 [17A--17F 已通过]

本节替代继续调局部速度网格的路线。17A--17F 已按顺序完成并通过；不得因为阶段
验收已通过，就把“能够生成PIC粒子”误写为 120 fs 生产完成，也不得恢复当前
\`extract_after_substep()\` 的逐 cell 静态扫描。

**完成记录（2026-08-10）**：17A 的接口/parcel 纯函数、17B 的只读 audit、17C 的
sink+loader 事务闭合、17D 的 1/2/5 rank 与 checkpoint/restart、17E 的无碰撞
10.4--12 fs static/flux A/B，以及 17F 的碰撞保守面通量独立总门和 1/2/5 rank MPI
门均已通过。17F 的 `all.result` 九个单 rank 字段与三个 `mpi_n*.result` 的分区字段均为
PASS。该记录只覆盖本节的 staged tests；后续 H9 Beam/no-Beam 生产回归仍按 §17.10 执行。

##### 7.11.17.1 根本目标与禁止项

当前静态转换执行的是：速度推进后扫描所有 \`partition.is_conversion(j,k)\` cell，取走
整个cell质量，再用少量点粒子拟合cell-center六矩。16A/16B已经证明该目标与圆柱
cell-volume的 $\Pi_\perp$ 不兼容。新路径必须改成：

$$
\text{bulk有限体积最终面通量}
\longrightarrow
\text{正权扫掠相空间质量包}
\longrightarrow
\text{tail PIC粒子}.
$$

必须满足：

1. bulk只损失真正穿越离散bulk--tail接口的质量，不再逐步静态清空整个高能cell；
2. tail创建量来自同一个最终、已限制、实际用于更新bulk的通量；禁止从未限制高阶通量、
   cell终态差或另一次重算得到；
3. 转换对象不能只有标量面通量。必须保留PPM扫掠区间及其正权求积节点，否则把全部
   质量放在阈值面仍不能闭合离散动能和二阶矩；
4. 不做全局能量补丁、负权粒子、事后缩放速度、随机拒绝采样或重复静态抽取；
5. 第一版不做延迟创建/pending accumulator。每个接受trial内立即创建粒子，避免pending
   电荷在Poisson、碰撞和checkpoint中成为第四种未推进表示；
6. \`tail_return_mode=none\` 保持不变。tail粒子以后即使降到阈值以下仍由PIC表示；
7. 旧静态路径只在A/B期间通过非默认CLI保留，通量路径通过17E/17F后删除生产调用点。

##### 7.11.17.2 离散接口拓扑

修改 \`src/grid.h\` 中 \`HybridVelocityPartition\`，新增显式面接口，不允许在remap热循环中
反复用能量条件搜索：

\`\`\`cpp
enum class VelocityFaceDirection { U_PARALLEL, U_PERP };

struct BulkTailInterfaceFace {
    VelocityFaceDirection direction;
    int face_index;
    int transverse_index;
    int bulk_iv;
    int bulk_imu;
    int tail_iv;
    int tail_imu;
    int outward_sign;
};

std::vector<BulkTailInterfaceFace> upar_interface_faces;
std::vector<BulkTailInterfaceFace> uperp_interface_faces;
std::vector<unsigned char> bulk_owned_cell;
\`\`\`

接口构造规则固定为：

1. 第一版使用face-aligned staircase接口；不能同时实现cut-cell/embedded boundary；
2. 只有当一个速度cell的全部体积求积节点均满足 $K\ge K_{out}$ 时，该cell才标记为
   tail-owned；与阈值曲面相交的cell继续归bulk，避免转换阈值以下支撑；
3. bulk-owned与tail-owned相邻时建立唯一接口面；每个共享面只能出现一次；
4. $u_\parallel$正负两侧必须镜像；$u_\perp=0$轴不能成为tail接口；
5. tail-owned区域必须连接到速度域外侧，禁止内部孤立tail孔洞；
6. 初始化和每次接受步后必须满足
   \`bulk_mass_in_tail_owned_cells <= roundoff_floor\`；物质级残留直接失败；
7. 将接口版本、$K_{out}$、mask哈希和face列表哈希写入manifest/checkpoint配置哈希。

新增 \`tests/bulk_tail_flux_interface_test.cpp\`，覆盖正负 $u_\parallel$ 对称、低
$u_\perp$ 阈值相交cell、无孤立孔洞、无重复面、所有tail侧求积节点均不低于阈值、
mask/hash确定性。

##### 7.11.17.3 扫掠质量包数据结构

新增独立、无MPI、无PIC依赖的文件：

\`\`\`text
src/bulk_tail_flux_parcel.h
src/bulk_tail_flux_parcel.cpp
\`\`\`

固定接口：

\`\`\`cpp
struct FluxParcelNode {
    double upar;
    double uperp;
    double mass;       // electrons per transverse area, >= 0
};

struct BulkTailFluxParcel {
    int ix_local;
    int ix_global;
    VelocityFaceDirection direction;
    int face_index;
    int transverse_index;
    int operator_stage;
    std::vector<FluxParcelNode> nodes;
    double number;
    double px;
    double jx_dx;
    double kinetic_energy;
    double pixx_dx;
    double piperp_dx;
};

struct BulkTailFluxBatch {
    std::vector<BulkTailFluxParcel> parcels;
    bool finite;
    bool nonnegative;
    double quadrature_error_max;
};
\`\`\`

六矩必须由 \`nodes\` 调用现有统一 \`mass_cell_moments()\` 求和得到，禁止分别维护另一套
解析公式。\`mass\` 与 \`Species::f\` 一样是cell线密度数目，不再乘 \`dx\`。

PPM扫掠包不是“face速度乘标量质量”。对 \`u_parallel\` remap：

1. \`upar_swept_mass()\` 当前积分的每个donor完整cell/部分cell都要保留积分区间；
2. 在每个区间上用固定Gauss--Legendre求积积分最终速度
   $u_\parallel^{n+1}=u_\parallel^n+a_u\Delta t$；
3. $u_\perp$ 使用圆柱体积变量 $s=u_\perp^2$ 求积，不能按 $u_\perp$ 几何宽度均匀取点；
4. 初始使用 \`quadrature_order_upar=4\`、\`quadrature_order_uperp=4\`；测试另跑4/8阶。
   8阶相对4阶的六矩差必须小于 $10^{-11}$，否则生产不能固定4阶；
5. 节点权重必须来自实际已限制PPM抛物线且非负。舍入级小负权可归零并记入
   \`roundoff_weight_debt\`，物质级负权使trial失败；
6. 每个节点最终能量必须满足tail支撑条件。低于阈值超过 $64\epsilon$ 尺度的节点
   表示接口或扫掠区间错误，禁止静默删除。

新增 \`tests/upar_flux_parcel_test.cpp\`，包含常数、线性、受限抛物线、跨1格、跨多格、
正/负加速度镜像和阈值擦边案例。解析数目积分误差不超过 $10^{-13}$；4/8阶六矩差
不超过 $10^{-11}$；所有节点非负且无阈值以下泄漏。

##### 7.11.17.4 将PPM更新改成“接口sink + parcel输出”

修改：

\`\`\`text
src/conservative_ppm_remap.h
src/conservative_ppm_remap.cpp
src/vlasov_split_step.h
src/vlasov_split_step.cpp
\`\`\`

\`advect_u_parallel()\` 增加可选参数，不能把partition设成全局变量：

\`\`\`cpp
RemapDiagnostics advect_u_parallel(
    const Species& input,
    Species& output,
    const EMFields& midpoint_fields,
    double dt,
    double time,
    const HybridVelocityPartition* partition,
    BulkTailFluxBatch* exported_flux);
\`\`\`

实现顺序固定：

1. 先按现有PPM重构和限制器得到最终 \`upar_swept_\`；
2. bulk--bulk面继续执行原flux-difference更新；
3. bulk--tail接口面仅在通量从bulk指向tail时，从bulk扣除最终通量、不写入tail-owned
   Eulerian cell，并由同一扫掠区间生成parcel；
4. tail侧Eulerian cell始终写0；tail-to-bulk方向没有Eulerian donor，通量固定为0；
5. 多cell departure只能在遇到的第一个bulk--tail接口导出一次；
6. 最外端 \`tail_number_loss\` 只记录未被内部接口捕获的真实越界质量；
7. PPM在接口bulk侧使用确定的一侧重构，不能把tail零cell纳入四点高阶stencil；
   点数不足时局部退化为线性或常数；
8. \`partition==NULL\` 或 \`exported_flux==NULL\` 时必须逐bit保持当前remap结果。

\`RemapDiagnostics\` 新增：

\`\`\`text
interface_export_number
interface_export_energy
interface_parcel_count
interface_node_count
interface_duplicate_count
interface_below_threshold_number
tail_owned_bulk_residual
\`\`\`

新增 \`tests/upar_flux_sink_test.cpp\`。要求

$$
N_{bulk}^{after}+N_{parcel}=N_{bulk}^{before}
$$

相对误差不超过 $10^{-13}$；tail-owned Eulerian质量为舍入级；多格扫掠
\`interface_duplicate_count=0\`。原 \`conservative_upar_remap_test\` 必须在NULL接口模式
继续逐bit通过。

##### 7.11.17.5 parcel到tail PIC的确定性加载

重写 \`src/bulk_tail_converter.h/.cpp\` 的生产职责。保留诊断结构和共享压缩器，但新增：

\`\`\`cpp
BulkTailConversionDiagnostics convert_flux_batch(
    const BulkTailFluxBatch& batch,
    BackgroundTailPIC& tail_trial,
    const SpatialGrid& grid,
    const HybridVelocityPartition& partition,
    int accepted_step,
    ConversionLocation location,
    int mpi_rank);
\`\`\`

\`extract_after_substep()\` 在开发期只供 \`static-cell\` A/B调用；\`flux-interface\` 模式禁止
调用。加载规则：

1. 按 \`(ix_global, operator_stage, sign, energy_bin)\` 聚合parcel节点；
2. 每个 $(u_\parallel,u_\perp)$ 支撑创建
   $\phi=0,\pi/2,\pi,3\pi/2$ 的等权quartet；
3. 先由原始正权节点计算目标六矩，再调用 \`tail_compress_moment_supports()\` 压缩到最多
   7个quartet支撑；
4. 压缩失败时可保留原始求积支撑，但必须通过本步粒子数预算；不得退回cell-center
   quartet或负权最小二乘；
5. 粒子空间位置取Eulerian空间cell中心，保证CIC重新沉积与bulk source同cell；
6. ID只能从 \`tail_trial.next_particle_id(mpi_rank)\` 获取，失败trial不得推进已接受ID；
7. created六矩必须与parcel目标比较，不能与旧静态cell目标比较；
8. 创建后对trial做局部密度沉积审计，但正式轨迹电流仍由后半漂移产生，不能伪造
   representation conversion的x电流。

硬门：$N/P_x/K$相对L1残差不超过 $10^{-10}$；$J_x/\Pi_{xx}/\Pi_\perp$ 不超过
$10^{-9}$；负权、非有限值和重复ID均为0；每group最多7个压缩支撑。

新增：

\`\`\`text
tests/bulk_tail_flux_loader_test.cpp
tests/bulk_tail_flux_transaction_test.cpp
tests/bulk_tail_flux_loader_mpi_test.cpp
\`\`\`

覆盖非单位质量、正负 $u_\parallel$、低 $u_\perp$、多parcel同group、压缩成功、
压缩fallback、故意NaN/负权失败、失败后tail粒子/ID/账本不变，以及固定全局输入跨
1/2/5 rank重新分区后全局六矩一致。

##### 7.11.17.6 生产积分器接入与CLI

修改：

\`\`\`text
src/vpfp_integrator.h
src/vpfp_integrator.cpp
src/main_vpfp.cpp
src/vpfp_diagnostics.h
src/vpfp_diagnostics.cpp
src/vpfp_checkpoint.h
src/vpfp_checkpoint.cpp
\`\`\`

新增CLI：

\`\`\`text
--tail-conversion-mode static-cell|flux-audit|flux-interface
--tail-flux-quadrature-order 4
--tail-flux-max-supports 7
--tail-flux-max-created-particles-per-step <N>
\`\`\`

开发期间默认仍为 \`static-cell\`；17E通过后默认改为 \`flux-interface\`，随后删除生产中
\`extract_after_substep()\` 调用。三种模式：

- \`static-cell\`：当前基线，仅作A/B；
- \`flux-audit\`：PPM输出parcel并写审计，但bulk/tail仍走原静态路径，物理状态必须与
  static-cell逐bit一致；
- \`flux-interface\`：接口sink生效且只用parcel创建tail，禁止再扫描转换cell。

无Beam和有Beam两条 \`advance_*()\` 必须调用同一个私有函数
\`apply_upar_flux_conversion()\`。调用位置保持在tail kick之后、tail第二次空间半漂移
之前；新粒子不重复接受bulk remap已经包含的电场kick，但必须参加后半漂移。

背景质量账改为

$$
N_{bulk}^{after}=N_{bulk}^{before}+N_{x,in}-N_{x,out}
-N_{u,outer}-N_{interface}.
$$

manifest/checkpoint必须保存conversion mode、接口hash、求积阶数和最大支撑数。restart
配置不一致必须拒绝。第一版无pending状态，因此不得虚构pending checkpoint字段。

旧schema checkpoint的兼容规则必须显式实现：缺少conversion字段时只允许解释为
`static-cell`。`flux-audit` 是只读诊断覆盖，可从该旧checkpoint启动，但写出的新checkpoint
必须包含完整新字段；`flux-interface` 不得直接从缺少接口hash的旧checkpoint启动。
物理状态hash与诊断配置hash必须分开计算，避免只读audit选项使物理A/B无法比较。

新增接受态 \`bulk_tail_flux_accepted_steps.dat\`，至少写step/time/location、
parcel/node/support/particle数、六矩exported/created/residual、阈值以下数目、
tail-owned残留、重复导出数、fallback和wall time。trial失败只写
\`trial_bulk_tail_flux_failures.dat\`。

##### 7.11.17.7 碰撞路径的明确边界

当前 \`moment-closure\` 生产分支是局部BGK式更新，且 \`bulk_collision_mask_\` 在
bulk--tail边界设置零通量墙。它没有可用于转换的最终速度面通量。因此：

1. 17A--17E先验收电场 $u_\parallel$ remap产生的通量转换，不得宣称完整碰撞VPFP
   转换已完成；
2. \`flux-interface + collision-model!=none\` 在17F完成前只能以显式
   \`collision_interface_mode=zero-wall-validation\` 运行，并在manifest写
   \`collision_induced_conversion=0\`；该模式不能用于最终生产；
3. 不允许从BGK更新前后cell差反推伪面通量；
4. 最终碰撞生产必须让漂移、平行扩散、垂直扩散和交叉扩散都暴露实际最终保守面通量；
5. 当前mask接口face必须显式选择 \`zero-wall\` 或 \`exporting-absorbing\`；
6. $u_\perp$ parcel使用面半径和 $u_\parallel$ cell正权求积；交叉扩散必须先重构成
   共享面保守通量；
7. 每个碰撞substep可产生parcel，但只在完整collision half成功后提交；
8. 新增 \`CollisionDiagnostics::interface_flux_batch\` 或等价事务返回值，并分别标识
   \`AFTER_COLLISION_HALF_1/2\`。

具体重构步骤固定如下，实施者不得自行用cell差分简化：

1. 将“碰撞系数模型”和“碰撞离散积分器”解耦。`moment-closure` 继续表示
   `CollisionCoefficientProvider`，新增
   `--bulk-collision-integrator bgk-validation|chang-cooper-flux`；当前BGK路径只允许作为
   无接口基线，`flux-interface` 生产必须使用 `chang-cooper-flux`；
2. 在 `cylindrical_fp_collision.h` 新增 `CollisionFaceFluxes`，分别保存每个collision
   substep最终使用的 `upar_flux`、`uperp_flux` 和cross-flux；数组所有权属于一次trial，
   预分配后复用；
3. 对Backward-Euler/Chang--Cooper隐式扫掠，必须先求得新状态，再用新状态和组装矩阵时
   完全相同的face系数重算最终通量，并检查
   $M^{new}-M^{old}+\Delta t D_vF^{new}=0$。该残差超过线性求解容差即trial失败；
4. bulk--tail接口不能继续跳过矩阵face。将该face作为单向吸收边界写入bulk侧矩阵行；
   outward通量进入parcel，inward通量在 `tail_return_mode=none` 下固定为0。禁止先用
   zero-wall求解再从结果外推一个“导出通量”；
5. `u_parallel` 接口parcel的法向速度固定为face坐标，横向在对应
   $u_\perp^2$ cell内正权求积；`u_perp` 接口parcel的法向速度固定为face半径，横向在
   $u_\parallel$ cell内正权求积。扩散穿越只表示穿过该面，不使用电场PPM扫掠区间；
6. cross-diffusion当前为显式共享面通量。必须先形成唯一face数组，再同时用于两侧
   FV更新和接口导出；不得分别在两个方向重算。符号指向tail的一侧才导出；
7. 每个substep先在trial `CollisionFaceFluxes` 中累计，完整collision half通过有限性、
   正性、质量和线性残差门后，再合并为 `BulkTailFluxBatch` 并调用loader。任一substep
   失败时丢弃整批parcel、tail trial和ID变化；
8. 碰撞parcel六矩按“穿越面上的正质量”定义；若任一接口face得到物质级负导出质量，
   必须判为方向/边界离散错误，不能创建负权粒子，也不能对全局结果取绝对值；
9. `CollisionDiagnostics` 至少新增
   `interface_export_number/energy`、`implicit_flux_residual_linf`、
   `cross_flux_pair_residual_linf`、`interface_inward_clipped_number`、
   `interface_parcel_count` 和 `transaction_rollback_count`；
10. `interface_inward_clipped_number` 只能作为 `tail_return_mode=none` 的模型诊断。若其
    相对combined背景数目长期可见，说明单向表示分解不再适用，必须启动bulk-return设计，
    不能通过继续放宽门槛进入120 fs生产。

新增：

\`\`\`text
tests/collision_flux_interface_test.cpp
tests/collision_flux_interface_cross_test.cpp
tests/collision_flux_conversion_transaction_test.cpp
tests/collision_flux_interface_mpi_test.cpp
\`\`\`

验收零系数严格零parcel；纯漂移、纯平行扩散、纯垂直扩散解析案例质量闭合；交叉扩散
共享面全局和为0；两个碰撞半步分别闭合；1/2/5 rank结果一致。

##### 7.11.17.8 分阶段验收门

**17A：接口与parcel纯函数。** 只新增拓扑、数据结构和纯函数测试。要求所有测试PASS、
节点非负、4/8阶矩差 $\le10^{-11}$、接口无重复/孔洞/阈值以下支撑。

**17B：只读flux-audit。** 接入PPM parcel输出但保持静态物理路径。严格只读门改为在
同一次 $u_\parallel$ remap 中比较parcel观察前后的最终面通量，并在观察调用前后比较
背景/场/Beam/tail状态、RNG和物理账本摘要；三项必须逐bit一致。独立OpenMP/MPI运行的
checkpoint全文件hash仅作信息记录，因为已实测static-A/static-B在约$10^{-11}$相对量级
存在并行浮点非确定性，不能作为audit因果门。每个parcel的数目必须与产生它的最终PPM接口面通量逐面闭合到
$10^{-13}$，全局parcel数目等于所有接口导出通量之和到 $10^{-12}$；无重复、无阈值
以下质量；额外wall time不超过15%。旧静态整cell抽取量与接口通量不是同一个离散对象，
只能记录二者差异，禁止把二者相等设为通过条件。

**17C：sink + loader独立闭合。** 要求bulk损失=parcel=tail创建，六矩满足17.5硬门，
tail-owned残留为舍入级，失败事务不污染状态。

**17D：MPI与checkpoint。** 1/2/5 rank固定物理输入，checkpoint往返和一步restart。
要求六矩、接口hash、粒子ID集合、下一ID和累计账本一致。

**17E：无碰撞10.4--12 fs生产A/B。** A为static-cell，B为flux-interface，均使用
\`collision-model=none\`。B必须全部接受，无静态扫描调用；组合数目闭合
$\le10^{-10}$；转换六矩满足硬门；Gauss不退化；tail粒子数和step wall time不超过A的
2倍；旧约3.3%的 $\Pi_\perp$ conversion residual降至 $10^{-9}$ 内。

**17F：碰撞面通量。** 完成17.7后才重跑collision pair、反作用、Beam 12 fs和no-Beam
40 fs。未完成17F时禁止120 fs生产。

**状态：已通过。** 17F 的 aggregate `all.result` 已验证零通量、接口导出、交叉扩散、
纯漂移、平行/垂直扩散、两个 collision half 与失败事务；`mpi_n1/n2/n5.result` 已验证
MPI 分区一致性。其后的H9生产前置门、Beam 12 fs和40 fs controller-off也已通过；
当前状态和后续工作统一见§0.11、§15.15和§17.13；§17.11仅保留H9/H11历史命令。

##### 7.11.17.9 命令位置

本阶段的全部编译、运行和结果判定命令已统一迁移到
§17.10.2“17A–17F 通量式转换验收命令 [已完成归档]”。本节只保留算法规范、
验收定义和完成结论，禁止在算法章节继续维护命令副本。


##### 7.11.17.10 最终完成定义

**当前状态**：17A--17F及collision-enabled H9生产级阻断门均已通过。下列第1--9项
的代码与推进验收已经完成；启动120 fs前仍需完成独立的阈值能谱物理门和资源门，
执行位置为 §17.11，而不是重新运行 §17.10。

只有以下条件同时满足，第7项才算完成：

1. 生产积分器不再调用逐步静态 \`extract_after_substep()\`；
2. 电场remap和碰撞算子均从实际最终保守通量产生parcel；
3. bulk损失、parcel、tail创建六矩和组合背景账本闭合；
4. tail-owned Eulerian区域保持空且没有每步静态清扫；
5. 1/2/5/80 rank、checkpoint/restart和失败事务通过；
6. collision-none 12 fs与collision-enabled H9阻断门通过；
7. 不再出现旧静态路径约3.3%的 $\Pi_\perp$ representation residual；
8. 性能和粒子数通过本节门；
9. manifest写明 \`tail_conversion_mode=flux-interface\`、接口hash、求积阶数和
   \`collision_induced_conversion=1\`。

任一项缺失，只能标记为“flux conversion部分实现”，不能开始120 fs生产。

禁止把审计误差用全局能量补丁、负权PIC粒子、重复静态抽取、增大拟合迭代次数或
无限扩大 $N_v/u_{max}$ 隐藏。

---

## 8. 混合无碰撞时间推进

### 8.1 接受态

时间层 $n$ 的接受态为

```text
bulk^n, tail^n, beam^n, E^n,
conversion_ledger^n, boundary_ledger^n,
tail_id_state^n, tail_rng_state^n
```

### 8.2 一个完整步的顺序

无碰撞、Beam 开启时按以下顺序实现：

1. 复制或交换得到完整 trial 状态。
2. bulk 执行 $x$ 半步，得到 `bulk_x_half`。
3. 原有 tail 粒子执行第一个空间半漂移，得到 `tail_x_half`，完成 MPI 迁移和开放边界截断。
4. Beam 按当前事件表执行 `predict_to_midpoint()`，得到 `beam_x_half`。
5. 分别沉积 `bulk_x_half`、`tail_x_half`、`beam_x_half` 的 midpoint 密度。
6. 组装
   $$
   \rho^{n+1/2}=e(Zn_i-n_{\rm bulk}^{n+1/2}-n_{\rm tail}^{n+1/2}-n_b^{n+1/2}).
   $$
7. 用现有非周期求解器得到 $E^{n+1/2}$。
8. bulk 使用 $E^{n+1/2}$ 执行完整 $u_\parallel$ remap。
9. 原有 tail 粒子使用同一 $E^{n+1/2}$ 执行完整 kick。
10. Beam 使用同一 midpoint 场完成 kick。
11. 对 post-$u$ bulk 执行 bulk-to-tail 转换；新粒子创建在 midpoint 空间位置，不重复 kick。
12. 所有 tail 粒子执行第二个空间半漂移，完成 MPI 迁移和开放删除。
13. bulk 执行第二个 $x$ 半步。
14. Beam 完成第二个漂移和开放删除。
15. 计算最终 bulk、tail、Beam 密度和所有矩。
16. 求解 $E^{n+1}$ 并执行完整 Gauss 审计。
17. 验证转换、边界、粒子迁移、总能量和有限性。
18. 全部通过后一次性接受 trial 状态。

### 8.3 不能改变的时间层关系

- midpoint Poisson 必须使用三个电子分量同一时间层的密度；
- bulk 与 tail kick 必须使用同一个 $E^{n+1/2}$；
- 新转换粒子不能再次接受已经包含在 bulk $u$ remap 中的 kick；
- 新转换粒子必须参加后半空间漂移；
- 最终 Poisson 必须包含转换后的 tail 密度和剩余 bulk 密度。

### 8.4 Beam 关闭路径

`--beam-enabled 0` 仍要走同一个混合积分器，只令 Beam 密度、粒子和注入账本为空。禁止维护一套独立且逐渐分叉的 background-only 生产算法。

---

## 9. Poisson、电荷沉积与连续性

### 9.1 场接口修改

修改 `src/maxwell.h/.cpp`：

```cpp
void EMFields::set_charge_density(
    const Species& bulk,
    const std::vector<double>& tail_density,
    const std::vector<double>& beam_density,
    const std::vector<double>& ion_density);
```

不要在调用端先把 tail 和 Beam 密度混成一个无标签数组。诊断需要分别识别来源。

### 9.2 场方程不使用 PIC 电流推进

Poisson 只需要三个电子分量的密度。tail 和 Beam 的轨迹电流仍应沉积，用于：

- 离散连续性审计；
- 场功与粒子动能诊断；
- 开放边界粒子收支；
- 检查沉积与 pusher 的一致性。

禁止因为已有电流沉积而恢复 Ampere 作为主场推进。

tail 连续性审计必须包含：

- 步开始和结束密度；
- 左右开放出流源；
- bulk-to-tail 转换源；
- 以后可能存在的 tail-to-bulk 返回汇；
- 完整轨迹电流。

combined 连续性残差应把 bulk reservoir 源、tail 开放出流和 Beam 注入/出流分别输出，禁止把所有边界项合成一个无法定位的标量。

### 9.3 Gauss 验收

每个接受步至少检查：

$$
r_G=D_xE-\rho/\varepsilon_0.
$$

输出 `gauss_l2`、`gauss_linf` 和端点积分误差。转换测试还必须比较 conversion 前后的 $\rho$，以区分 Poisson 算子错误和转换沉积错误。

---

## 10. 碰撞重构

### 10.1 当前状态 [已完成-需接口扩展]

当前 `ZeroCollisionCoefficients` 可继续用于无碰撞验收。`PrescribedCollisionCoefficients` 和现有圆柱碰撞器只适合单元测试，不足以代表目标的速度依赖 FP 碰撞。

**[历史 H7 实现记录]**：bulk 碰撞算子曾按第 10.2 节重构——
系数在每个 $(x,u_\parallel,u_\perp)$ 速度点计算、支持非对角扩散
`d_parallel_perp`、圆柱几何项来自保守散度形式；新增第 10.2.1 节
`moment-closure` 模式（`MomentClosureCollisionCoefficients`，局部漂移
Maxwellian 闭合 + 速度依赖碰撞率 $\nu(u)=\nu_0/(1+(u/u_{th})^2)^{3/2}$，
与 Maxwellian 满足 Einstein 关系）。`self_consistent_landau`（全分布
积分）留待 H8。`tail collision=none` 保持（tail 后端属 H8）。

**当前实现定位**（H9 后复核）：

- `prescribed` 模式仍使用速度依赖的圆柱 FP/Chang--Cooper 通量算子，支持非对角扩散 `d_parallel_perp`；
- `moment-closure` 在 H9 稳定性修复后已改为局部、速率依赖的隐式 BGK 松弛，并用局部矩投影恢复 $N/P_\parallel/K$；
- 该 BGK 路径是 **moment-closure 近似后端**，不等价于目标 Kramers--Moyal/Landau FP 漂移-扩散算子，禁止在报告中将它标记为“完整 VPFP”或“self-consistent Landau”；
- `self_consistent_landau` 尚未实现；
- H7 通过证明了数值算子、保守通量和单元测试基础可用，不等于最终目标碰撞物理已完成。

### 10.2 bulk 碰撞系数

修改：

```text
src/collision_coefficients.h
src/collision_coefficients.cpp
src/cylindrical_fp_collision.h
src/cylindrical_fp_collision.cpp
```

要求：

1. 在每个 $(x,u_\parallel,u_\perp)$ 点计算漂移和扩散系数；
2. 支持必要的非对角扩散项；
3. 圆柱几何项必须来自保守散度形式；
4. 系数输入必须遵守第 10.2.1 节的模式契约，不能默认用少量 combined moments 代替完整分布；
5. 离散碰撞器分别验证数目守恒、动量交换、能量交换、Maxwellian 平衡和 H 定理趋势；
6. 碰撞引起的外部储库能量交换必须显式记账。

#### 10.2.1 碰撞系数模式契约

必须在 manifest 中明确选择以下一种模式：

1. `prescribed`：$\langle\Delta\boldsymbol v\rangle$ 和 $\langle\Delta\boldsymbol v\Delta\boldsymbol v\rangle$ 由外部模型直接给定。bulk 和 tail 在相同 $(x,\boldsymbol u,t)$ 上调用同一个 provider。
2. `moment_closure`：系数由密度、漂移、温度等有限矩近似得到。输入矩使用 bulk+tail combined 值，并明确这是闭合近似。
3. `self_consistent_landau`：系数来自总背景分布的 Rosenbluth 势或等价全分布积分。此时 tail 对系数的贡献必须通过速度空间沉积或粒子求和进入，不能只折算为局部温度。

`LocalCollisionMoments` 只能服务于第 2 种模式。它不能作为第 3 种模式的完整输入接口。

**当前限制**：H9 四档结果使用的是第 2 种 `moment_closure/BGK` 配置，它只能验证混合调度、数值稳定性和闭包账本，不能替代对第 1 种用户指定 $\omega$ 或第 3 种自洽 Landau 系数的验收。

#### 10.2.2 同一物种分解

背景电子被拆成两种表示后，自碰撞算子在概念上分解为

$$
C[f_e,f_e]
=C_{bb}+C_{bt}+C_{tb}+C_{tt},
$$

其中第一个下标表示被推进的 target，第二个下标表示产生碰撞系数的 field population。

- bulk FP 更新包含 $C_{bb}+C_{bt}$；
- tail PIC 碰撞更新包含 $C_{tb}+C_{tt}$，具体由 SDE 或 Nanbu--Perez 后端承担；
- 同一个 pair contribution 只能出现一次；
- 若省略 $C_{tt}$ 或 $C_{bt}$，必须记录稀疏近似和误差指标。

实现前先画出 provider 调用图和 pair mask。测试必须证明关闭某一 pair 后只有对应贡献消失，防止 combined provider 与显式反作用同时接入造成双计数。

#### 10.2.3 多物种 pair registry

新增显式 `CollisionPairRegistry`。每个条目至少包含：

```text
target_population
field_population
coefficient_model
discrete_operator
tail_collision_kernel
weight_mode
reaction_destination
conservation_contract
enabled
```

至少区分 `bulk-bulk`、`bulk-tail`、`tail-bulk`、`tail-tail`、`electron-ion`。Beam 默认保持当前 collisionless 行为；若以后启用 Beam-background 碰撞，必须新增 `beam-background` 与 `background-beam` pair，并把反作用写入 Beam 或 background，不能只让其中一方减速。

### 10.3 tail PIC 碰撞后端

**[统一生产调度和不等权有界修正已完成；H9 多 rank 基本路径已验收]**：统一接口与 Nanbu--Perez/SDE 后端、不等权重 `equal-strata/virtual-split`、冷等离子体截止及内部子循环已经实现。为兼容既有 CLI 和 checkpoint，`virtual-split` 名称暂时保留，但其生产语义已经改为**不生成永久 residual 宏粒子的有界 Sentoku--Kemp 修正**。`--tail-collision-weight-mode`、`--tail-collision-max-substeps` 和零增长预算从 CLI 贯穿 `VpfpIntegrator`、`HybridCollisionConfig`、具体碰撞请求、checkpoint 和 snapshot manifest；超过子步上限时事务失败，不接受部分更新。默认子步上限为 `1024`，默认宏粒子增长预算为 0。2026-08-06 已抽取 `VpfpIntegrator::apply_collision_half()`，`advance_background()` 与 `advance_with_beam()` 的两个 Strang 碰撞半步均只通过该入口选择 bulk-only 或完整 `HybridCollisionStep`，不再存在 Beam 分支绕过 tail 后端的问题。单 rank 和新 80-rank 生产路径均已观察到 `bb/tt/tb/bt=1/1/1/1`，且 Beam 12 fs 不再发生碰撞诱发的 residual 粒子增长。`TraceStationaryBackground` 仍是未实现的可选参考模式（§10.3.4）；H9 的剩余阻断项是 §7.11 的阈值接口和资源预算，不是该碰撞后端。

EPOCH 源码审计表明，tail PIC 碰撞必须按目标碰撞核选择后端，不能把一种算法宣称为所有 $\omega(\boldsymbol v,\Delta\boldsymbol v)$ 的通用实现。

新增统一接口：

```text
src/background_tail_collision.h
src/background_tail_collision.cpp
src/background_tail_nanbu_perez.h
src/background_tail_nanbu_perez.cpp
src/background_tail_collision_sde.h
src/background_tail_collision_sde.cpp
tests/background_tail_nanbu_perez_test.cpp
tests/background_tail_collision_moments_test.cpp
tests/background_tail_collision_equilibrium_test.cpp
tests/background_tail_collision_weight_test.cpp
tests/tail_collision_isotropisation_test.cpp
tests/tail_collision_two_population_test.cpp
```

```cpp
enum class TailCollisionKernel {
    None,
    CoulombLandauNanbuPerez,
    KramersMoyalSDE,
    TraceStationaryBackground
};

struct TailCollisionRequest {
    TailCollisionKernel kernel;
    double dt;
    int accepted_step;
    int collision_half;
    CollisionPairMask pairs;
};
```

`TraceStationaryBackground` 只能用于与 EPOCH `background_collisions.F90` 对照的受控近似，禁止作为完整 VPFP 默认生产模式。

#### 10.3.1 Coulomb/Landau：Nanbu--Perez tail--tail 后端

当 $\omega$ 对应库仑小角散射、目标连续算子为 Landau/Fokker--Planck 极限时，tail--tail 优先参考 EPOCH 4.20.1 默认 Nanbu--Perez 算法：

1. 每个碰撞子步开始时，按物理空间 cell 对 tail 粒子建立 index bins；禁止修改主粒子存储顺序。
2. 在每个 cell 内用 counter-based permutation 随机配对。permutation key 至少由 `(particle_id, accepted_step, collision_half, substep, pair_pass)` 构成，不能依赖 vector 顺序或 MPI rank；同一物种配对数按 $\lceil N_p/2\rceil$ 处理，奇数粒子通过确定性轮转避免长期固定自配对。
3. 把两粒子四动量 Lorentz 变换到质心系；根据局部密度、相对速度、Coulomb logarithm 和 $\Delta t_c$ 计算散射强度 $s_{12}$。
4. 按 Perez 分段反演采样 $\cos\theta$，均匀采样方位角，在质心系旋转相对动量，再变换回实验室系。
5. 对每个实际发生的等权重二体事件，必须在浮点求和误差内保持加权总三动量和相对论总能量。
6. `coulomb_log`、冷等离子体修正和碰撞频率上限必须由独立配置控制并写入 manifest；不得隐藏沿用 EPOCH 常数。
7. tail--tail 碰撞不产生 bulk 反作用，因为它是 tail 内部交换；其 combined tail $N/P/K$ 残差必须单独验收。

在本项目 1D、单位横向面积约定下，碰撞 cell 的物理数密度必须定义为

$$
n_{t,i}=\frac{1}{\Delta x_i}\sum_{p\in i}w_p,
$$

单位为 $\mathrm{m^{-3}}$。pair 数只决定 Monte Carlo 采样数，不能替代物理密度。所有由 EPOCH/Perez 公式移植的 sampling factor 必须重新按 `weight [m^-2]`、$\Delta x$ 和 $\Delta t_c$ 做量纲推导，并用均匀 cell 加倍 $w$、加倍 $\Delta x$、加倍粒子数的独立缩放测试验收。

这里参考的是 `collisions.F90:442-642,896-1111` 的相对论质心系散射结构，不逐行复制其 Fortran 实现。

#### 10.3.2 任意 Kramers--Moyal 核：SDE 后端

若用户提供的是一般

$$
\boldsymbol A_u=\frac{\langle\Delta\boldsymbol u\rangle}{\Delta t},
\qquad
\boldsymbol D_u=\frac{\langle\Delta\boldsymbol u\Delta\boldsymbol u\rangle}{\Delta t},
$$

而不是明确的 Coulomb/Landau 核，则使用 SDE 后端：

$$
d\boldsymbol u=\boldsymbol A_u\,dt+\boldsymbol B_u\,d\boldsymbol W,
\qquad \boldsymbol B_u\boldsymbol B_u^T=\boldsymbol D_u.
$$

它必须从与 bulk 同源的 coefficient provider 构造 $\boldsymbol A_u$、$\boldsymbol D_u$，并包含从 $\boldsymbol v$ 到相对论 $\boldsymbol u$ 的 Ito/Jacobian 修正。固定速度点大样本测试必须验证一、二阶弱矩收敛；不能用 EPOCH 二体散射替代一个物理核不同的用户给定 $\omega$。

#### 10.3.3 不等宏粒子权重

bulk-to-tail 转换与人口控制会产生不等 `weight`。EPOCH 4.20.1 的 `collisions.F90::weighted_particles_correction()` 使用 Sentoku--Kemp 不等权修正：轻粒子采用完整散射态；重粒子采用未散射态和散射态的权重混合，并添加随机横向动量恢复目标能量。该方法逐事件保持加权相对论能量，动量在统计期望上守恒，并保持固定宏粒子数。

生产实现按以下优先级处理：

1. converter 和 population controller 尽可能把粒子限制在少量离散权重层；碰撞只在同 cell、同权重层或权重比接近 1 的层间配对。
2. 等权重对使用严格二体质心系更新。
3. 不等权重对使用有界 Sentoku--Kemp 修正。轻粒子采用完整二体散射态；重粒子的混合比例为 $r=w_{\min}/w_{\max}$，先构造 $\boldsymbol p_m=(1-r)\boldsymbol p_{old}+r\boldsymbol p_{scat}$，再沿随机横向方向添加满足目标加权能量的动量。输入两个宏粒子，输出仍严格为两个，禁止把剩余权重实体化为永久新粒子。
4. 单事件验收要求粒子数和权重严格不变、能量达到求和精度；动量必须采用多 cell/多 seed 集合统计验收，不能错误要求单事件逐位守恒，也不能只检查均值而忽略方差和长期随机游走。
5. 禁止通过碰撞后全局缩放所有 tail 动量修正能量。

CLI token `virtual-split` 仅为 checkpoint/脚本兼容名，不再表示实体权重分裂。生产命令必须显式设置 `--tail-collision-max-particle-growth 0`；若该模式仍出现 `conversion_N=0` 而 tail 宏粒子数增加，应视为实现回归并立即失败，而不能依靠 population controller 事后压缩。

#### 10.3.4 快 tail--静止背景参考模式

EPOCH `background_collisions.F90` 的算法假定慢背景粒子静止，只旋转快粒子动量，不把反冲写回背景。它适用于验证以下极限：

- tail 能量远高于 bulk 热能；
- tail number/energy fraction 足够小；
- 只关心首阶俯仰角散射或 stopping benchmark。

本项目将其实现为可选 `TraceStationaryBackground` 测试后端。该后端必须在输出中标记：

```text
collision_approximation=trace_stationary_background
background_reaction_included=0
```

它不能通过第 19.3 节完整碰撞验收，也不能用于最终 120 fs 生产结论。

#### 10.3.5 碰撞子循环

所有后端都使用独立碰撞子步 $\Delta t_c$：

- SDE 根据局部漂移与扩散尺度限制典型 $|\Delta\boldsymbol u|$；
- Nanbu--Perez 根据最大 $s_{12}$、大角散射概率和权重分裂增长率限制；
- `coll_n_step>1` 式超循环只能在 $\Delta t_c$ 收敛测试通过后启用，不能仅因 EPOCH 支持就默认采用。

输出 `collision_substeps`、`max_s12`、`large_angle_fraction`、`weight_split_count`、`max_du` 和各后端 wall time。

#### 10.3.6 EPOCH 参考算例的本项目化

EPOCH 参考树包含：

```text
epoch1d/example_decks/electron_isotropisation.deck
epoch1d/example_decks/electron_ion_equilibration.deck
```

本项目不直接把这些周期全 PIC deck 当作生产测试，而是提取两个零维碰撞基准：

1. `tail_collision_isotropisation_test`：单 cell、关闭输运和场，初始化 $T_x\ne T_y=T_z$ 的等权重 tail，验证温度各向异性随粒子数增加和 $\Delta t_c$ 减小而收敛，同时保持总动量和总能量。
2. `tail_collision_two_population_test`：两个 PIC 测试种群初温不同，验证互碰后的温度趋同、总动量和总能量守恒；该测试只验证 PIC--PIC Nanbu--Perez 内核，不替代 bulk--tail 混合反作用测试。

两个基准都至少运行两组粒子数和 `dt/dt/2`，统计误差用多个 counter-based seed 估计。禁止只凭一条噪声曲线判定算法阶数。

`tail_collision_isotropisation_test` 按本条实现为 16 个 counter-based seed
的均值判据：同一物理阈值（mean_A<0.2·aniso_initial、|mean_A−mean_B|<0.02、
mean_C≤mean_A+1e-3、逐 seed P/K 守恒 <1e-12）只作用于 seed 均值，并输出逐
seed 扫描、首 seed 衰减曲线、子步数/max_s12/大角占比诊断；单条噪声曲线
不作为通过依据。

### 10.4 完整碰撞 Strang 顺序

新增：

```text
src/hybrid_collision_step.h
src/hybrid_collision_step.cpp
tests/hybrid_collision_pair_balance_test.cpp
tests/hybrid_collision_transaction_test.cpp
```

`HybridCollisionStep::advance()` 同时接收 bulk trial 和 tail trial，在同一个碰撞半步内完成 pair mask、系数构造、tail 后端选择、tail 随机/二体增量、bulk FP 更新和交叉反作用。只有 combined 碰撞账本通过后才返回成功。

最终一步为：

1. `HybridCollisionStep(bulk,tail,\Delta t/2)`；
2. collision 后 bulk-to-tail 转换；
3. 第 8 节的完整无碰撞混合步；
4. `HybridCollisionStep(bulk,tail,\Delta t/2)`；
5. collision 后 bulk-to-tail 转换；
6. 统一验收并接受。

同一碰撞半步中的 bulk 和 tail 必须使用相同时间层的系数输入。若采用隐式或迭代碰撞更新，combined 状态必须收敛后才能接受；不得让 bulk 使用 half-step moments 而 tail 使用 step-start moments。EPOCH 在主循环中于 particle push 后调用碰撞的顺序只作为全 PIC 参考，本项目必须保留此处与目标 VPFP 分裂一致的碰撞半步，不得为了“像 EPOCH”改变外层 Strang 顺序。

**[独立模块与两条生产路径均已接入；H9 多 rank applied 路径已验证]**：`CollisionPairRegistry` 已与算法后端解耦，独立 `HybridCollisionStep` 可以依次执行 $C_{bb}$ 的 bulk FP/BGK、$C_{tt}$ 的 Nanbu--Perez、$C_{tb}$ 的 tail--bulk SDE，以及与实际 tail 增量配对的 $C_{bt}$ bulk 局部反作用。SDE 单独请求仍可只启用 $C_{tb}+C_{bt}$；四个 pair 均可独立关闭，全部关闭时为精确 no-op。`advance_background()` 与 `advance_with_beam()` 现共用 `apply_collision_half()` 和 `finalize_collision_ledger()`；单 rank Beam+非零 tail 生产回归中四个 applied flag 均为 1，最大反作用相对残差为 $2.48\times10^{-15}$。新 80-rank Beam 12/40 fs 接受态也记录到四通道实际执行；局部相对残差在交换量接近零时可能病态，正式判定应同时使用带尺度底限的局部绝对残差和全局累计账本。

反作用使用离散模板实际施加后的 $\Delta P_x/\Delta K$ 记账，而不是请求值；事务测试保证任一后端或反作用失败时 bulk/tail 均不提交。该独立模块已通过本地 `hybrid_collision_pair_balance_test` 和 `hybrid_collision_transaction_test`，但只有 Beam 生产分支也调用同一模块后，这些结果才能支撑 H9。

### 10.5 tail-tail 与 bulk-tail 物理范围

不能默认声称 tail-tail 碰撞已包含。首版至少实现 tail 在 bulk 背景上的 test-particle 碰撞，并根据 tail number fraction 评估：

- 若 tail 稀疏，忽略 tail-tail 可作为受控近似；
- 若 tail 数目或能量占比不再稀疏，必须加入 field-particle 反作用或全耦合碰撞模型；Coulomb/Landau tail-tail 优先采用第 10.3.1 节的 pairwise Monte Carlo，一般 Kramers--Moyal 核才使用独立 SDE；
- 所有近似必须写入 manifest 和诊断，不能隐藏。

### 10.6 电子碰撞反作用

如果 tail 以 test-particle SDE 或 `TraceStationaryBackground` 形式实现 $C_{tb}$，而 bulk 的 $C_{bt}$ 没有由同一个混合 pair 算子自动产生，则 tail 的实际随机动量和能量变化必须在 bulk 方程中出现相反的 field-particle 反作用。否则 combined 背景电子会凭空丢失或获得动量和能量。tail--tail Nanbu--Perez 属于 tail 内部交换，不得重复向 bulk 写反作用。

每个空间 cell、每个碰撞半步至少累计：

$$
\Delta P_{t\leftarrow b},\qquad \Delta K_{t\leftarrow b}.
$$

bulk 的交叉碰撞源应满足

$$
\Delta P_{b\leftarrow t}=-\Delta P_{t\leftarrow b},
$$

$$
\Delta K_{b\leftarrow t}=-\Delta K_{t\leftarrow b}.
$$

这不是事后全局能量补丁，而是离散 $C_{bt}+C_{tb}$ 的必要组成。反作用必须局部施加，并通过非负、守恒的速度空间源离散实现。若 $C_{bt}$ 已由完整 field-population 算子显式计算，则不得再施加同一份反作用。

碰撞账本必须区分：

- electron-electron：combined 电子数、动量和能量内部守恒；
- electron-ion：动量和能量可进入 ion reservoir ledger；
- 外部规定散射：进入明确命名的 external collision reservoir。

若首版只采用 trace-tail test-particle 近似且忽略反作用，必须写入 `tail_collision_approximation=trace_test_particle`，并用 tail 能量占比给出误差上限。该模式不能称为完整 VPFP 生产解。

### 10.7 tail--bulk 生产实现的固定顺序

tail--bulk 不能直接调用 EPOCH 的 PIC--PIC `inter_collisions_np()`，因为 bulk 没有真实碰撞伙伴粒子；也不能直接调用其 stationary-background 路径，因为该路径没有背景反冲。生产实现按以下顺序：

1. 从碰撞半步开始的 combined 状态构造 bulk field-population 系数，得到作用于 tail 的 $\boldsymbol A_{t\leftarrow b}$、$\boldsymbol D_{t\leftarrow b}$。
2. 用 SDE 后端推进 tail 的 $C_{tb}$，逐 cell 累计所有实际宏粒子增量对应的 $\Delta P_{t\leftarrow b}$ 与 $\Delta K_{t\leftarrow b}$；累计使用真实 particle weight。
3. 独立计算 bulk 的 $C_{bt}$ 候选。若 `self_consistent_landau` 已从 tail 速度分布得到完整 field-particle contribution，则直接计算其离散 $\Delta P_{b\leftarrow t}$、$\Delta K_{b\leftarrow t}$。
4. 比较两侧局部交换矩。离散不匹配只能通过 **bulk 碰撞通量空间内的局部约束修正** 解决：修正保持 bulk 电子数不变，并使该 cell 的动量、能量交换等于 tail 的反量；不得修改无碰撞通量、Poisson 场或其他 cell。
5. 局部约束问题无非负可行解时，减小碰撞子步重算；超过最大子步仍不可行则事务失败并保留 accepted 状态。
6. 只有 $C_{tb}$、$C_{bt}$ 和局部反作用账本共同通过后，才接受该碰撞半步。

对于 `moment_closure`，第 3 步可以没有完整 $C_{bt}$ 形状，此时必须用明确命名的 `bulk_field_particle_reaction` 在 bulk 速度网格上分配反作用。分配基函数至少要能独立控制数目、平行动量和能量，且只在具有正质量预算的局部速度 stencil 上工作。该模式是闭合近似，manifest 必须标记，不能称为完整 Landau 离散。

可选研究后端可以从 bulk cell 分布确定性采样临时 partner markers，再做 EPOCH 风格 PIC--PIC 散射并把 partner 的前后差投影回 bulk；但在证明采样噪声、投影非负性和六矩误差优于上述 SDE+reaction 方案前，不作为默认生产路径。

---

## 11. 账本、保护与接受条件

### 11.1 电子数账本

接受步应满足

$$
N_{\rm bulk}^{n+1}+N_{\rm tail}^{n+1}
=N_{\rm bulk}^{n}+N_{\rm tail}^{n}
+N_{\rm reservoir,in}-N_{\rm reservoir,out}
-N_{\rm tail,leftout}-N_{\rm tail,rightout}.
$$

conversion 项不出现在右端，因为它在 combined 系统中相消。

### 11.2 能量账本

至少输出

$$
U_E+K_{\rm bulk}+K_{\rm tail}+K_{\rm beam},
$$

并累计：

- Beam 注入和左右出流能量；
- bulk reservoir 入流和出流能量；
- tail 左右出流能量；
- 碰撞储库能量；
- 数值 remap 误差；
- conversion 能量残差。

conversion 本身的净能量必须为零到离散求和误差。

开放非周期静电系统还可能存在由场边界条件和边界电荷通量共同产生的静电边界功。该项不能照搬周期系统公式，也不能凭连续公式单独估算。实现时必须从当前离散 Gauss/Poisson 算子、实际指定的左端场或电势条件以及离散连续性方程推导同一套边界能量项，并输出：

```text
electrostatic_boundary_work_step
electrostatic_boundary_work_cumulative
```

在该项完成前，总能量残差只能标为诊断量，不能用它触发物理状态修补。

当前 `OpenElectrostaticSolver::boundary_power()` 在生产 `DIRICHLET_PHI` 模式返回 $\phi_LI_L-\phi_RI_R$。当 $\phi_L=\phi_R=0$ 时该外加电极功为零，但仍必须通过离散 Poisson 能量恒等式验证符号、时间层和边界表面电荷项；不能仅凭公式中电势为零就省略边界审计。`LEFT_E` 返回零只适用于对照实现，不能外推为一般左场边界已能量闭合。

### 11.3 必须硬失败

- 任意接受候选中出现 NaN/Inf；
- Gauss 残差超过已验证阈值；
- conversion 的 $N/P_x/K$ 残差超过阈值；
- MPI 粒子迁移造成粒子数或权重不闭合；
- post-conversion bulk 在禁止区仍有超过容差的质量；
- bulk 质量到达最外速度边界并产生未表示出流；
- checkpoint 写入或读取不完整；
- 粒子 ID 重复；
- tail 出现非正或非有限宏粒子权重；
- 等权 tail--tail 二体事件的加权三动量或相对论能量残差超过求和阈值；
- 有界 Sentoku--Kemp 不等权事件改变宏粒子数/weight，或单事件加权相对论能量残差超过求和阈值；其动量按集合统计门验收，不对单事件随机横向修正使用严格动量硬门；
- 无 conversion/边界粒子源时，tail--tail 碰撞造成宏粒子数增长；
- 完整 electron--electron 模式下 tail--bulk 局部反作用不可行或 combined $P/K$ 不闭合；
- 内存或粒子数超过显式上限。

### 11.4 不应硬失败

- tail 粒子动量很高但仍有限；
- 边界场不为零；
- tail 粒子数发生物理增长但仍在资源预算内；
- EPOCH 与本程序后期逐点波形不同；
- 单个统计 realization 的 SDE 样本矩或概率权重碰撞不等于理论期望，但集合统计仍在置信区间内；
- 单个低权重粒子离开计算域。

### 11.5 资源与有效性门

新增：

```text
--tail-max-particle-count
--tail-max-number-fraction
--tail-max-weight-ratio
--tail-noise-warning-level
```

达到 warning level 时保存 checkpoint 和完整诊断；达到硬资源上限时有序停止，不允许操作系统 OOM 杀死。

---

## 12. checkpoint 与重启

### 12.1 格式升级 [已完成-保留]

升级 `src/vpfp_checkpoint.*` 的版本。每 rank 文件新增：

- 所有 `BackgroundTailParticle`；
- tail accepted ID counter；
- tail counter-based RNG key/version；
- tail 左右出流账本；
- conversion 累计 $N/P_x/K$；
- conversion 参数和聚合 bin 配置；
- population controller 配置、累计 splitting/merging/fallback 计数；
- bulk/tail 分区哈希；
- tail collision kernel、实现版本、weight mode、pair registry hash；
- collision subcycle 配置、Coulomb-log 模式和 counter-based RNG 算法版本；
- combined moment 校验值。

**[已完成-保留]**（2026-08-04 实施，阶段 H6）：

- schema 升到 v2（magic 不变、version=2），每 rank 二进制文件在
  f/beam/Ex_face/Ex/phi 之后追加 tail 段：全部 `BackgroundTailParticle`、
  accepted tail density、ID counter、counter-based collision RNG key、
  左右出流账本、截断/沉积 shape 账本、max_abs_u/max_kinetic_energy；
- config 段：partition 哈希、convert_energy_mev、buffer 宽、upar/energy
  bin、return_mode、collision kernel/weight mode/substeps/particle
  growth、population controller 全套配置与累计 groups/fallbacks、
  conversion 累计 N/Px/K 与 particles_created/outflow、combined
  N/K_kin/field 校验值；
- `BackgroundTailPIC` 新增 `export_accepted_state/import_accepted_state`
  （§14.4 扩展），checkpoint 写入/读出共用，逐位往返；
- manifest 新增 `background_representation`、`tail_return_mode`、
  `tail_convert_energy_mev`、`tail_conversion_bins`、
  `partition_config_hash`、tail 粒子/pusher/deposition/collision backend
  哈希、population control 参数、combined 校验值与
  `ix_start/nx_local` 空间分区；
- 写侧由 `main_vpfp.cpp` 从积分器 accepted tail + 累计账本组装
  `VpfpCheckpointTailState`；读侧校验 partition 哈希与 population
  control 配置，恢复 tail 状态、累计账本与 `step_count`。

manifest 新增：

- `background_representation=eulerian_bulk_plus_pic_tail`；
- `tail_return_mode`；
- `tail_convert_energy_mev`；
- 网格、阈值、形函数、粒子 schema 哈希；
- tail pusher/deposition/collision backend 哈希；
- MPI size 和空间分区。

### 12.2 兼容策略

旧的无 tail checkpoint 不自动升级为混合生产 checkpoint。允许提供一次性离线转换工具：

```text
vpfp_checkpoint_convert --input <old> --output <new> --tail-mode pic
```

转换工具必须调用同一个 `BulkTailConverter`，并写出 $N/P_x/K/\rho$ 残差。禁止在 `read_checkpoint()` 内静默修改物理状态。

**[已完成-保留]**（阶段 H6）：v1（无 tail）checkpoint 仅在 tail-off
求解器下可读；tail-on 求解器读 v1 直接报错（不自动升级）；含 tail 的
v2 checkpoint 被 tail-off 求解器拒绝。未提供 `vpfp_checkpoint_convert`
离线工具（§12.2 为可选项，当前无历史 v1 混合 checkpoint 需要转换）。

### 12.3 重启验收

必须比较连续两步与“一步 + checkpoint/restart + 一步”：

- bulk 全数组 hash 或相对 $L_2$；
- tail 粒子按 ID 排序后的状态；
- Beam 状态；
- 场；
- 所有账本；
- 下一次新建 tail 粒子 ID；
- 碰撞开启时 RNG 后继状态。

**[已完成-保留]**（阶段 H6，`tests/checkpoint_restart_equivalence_test`
单 rank）：Sim A 连续两步 vs Sim B 一步+checkpoint/restart+一步，最终
bulk f、beam 粒子、Ex_face/Ex、tail 粒子（按 ID 排序）、step 账本、
累计账本、combined 校验值与下一次新建 tail 粒子 ID 逐位一致。
说明：① 等价测试固定 OMP_NUM_THREADS=1（阶段 3/4 已记录的 OpenMP 归约
ULP 非确定性会使多线程下两步逐位比较不成立；生产多线程 restart 恢复的
是同一 accepted 状态）；② 比较排除 `fields.phi`——生产 final Poisson
`reconstruct_phi=false`，accepted phi 是 swap 携带的陈旧值、不参与任何
推进（checkpoint 仍逐位往返 phi，见 round-trip 测试）。

---

## 13. 诊断与输出

### 13.1 每个轻量接受步

只计算和输出：

- step、time、dt；
- $N/K/P_x$ 的 bulk、tail、combined、Beam 值；
- field energy；
- conversion 本步和累计残差；
- tail particle count（H5 起写入 `vpfp_step_diagnostics.dat` 末尾列
  `tail_particle_count`，全局宏粒子数；旧固定列号解析脚本会错位）；
- bulk forbidden-region mass；
- Gauss $L_2/L_\infty$；
- wall time 分解。

禁止在每步执行完整四维扫描以外的重复审计。

**[已完成-保留]**（阶段 H6）：`vpfp_step_diagnostics.dat` 追加
`wall_tail_push_s/wall_tail_deposit_s/wall_tail_migrate_s/
wall_conversion_s/wall_diagnostics_s/tail_particles_local_max` 六列
（§13.3；`wall_diagnostics_s` 为调用方填写的上一个接受步的诊断写墙钟，
其余在 advance 内计时）；固定列号解析脚本需同步。

**[66 列结构与代码侧生产语义已修复；新多 rank H9 已验证写出路径]**：`vpfp_step_diagnostics.dat` 已扩展为 66 列，且只在状态真实接受后写出。2026-08-06 已统一 Beam/background 碰撞路径，并在接受态写出前由 `globalize_conversion_ledger()` 对 conversion 的 $N/P_x/K$ 和创建粒子数做 MPI 全局归约，对相对残差取全局最大值；checkpoint 累计量仍在归约前按 rank-local 状态更新，避免破坏逐 rank restart 语义。阈值 forbidden-region 的四维扫描仍只在快照时执行，避免每步重复扫描拖慢生产。旧 H9 文件中的 pair 全零和 rank-local conversion 列仍是历史坏数据，不能用于验收新代码；新 H9 已验证四通道 applied 与全局 conversion 列，但 §7.11 所需的显式谱边界和压缩分组审计仍需补充。

### 13.2 快照时

输出：

- $E_x(x)$；
- $n_{\rm bulk}(x)$、$n_{\rm tail}(x)$、$n_b(x)$、combined density；
- bulk/tail/combined 平均动能；
- bulk Eulerian 能谱；
- tail 粒子能谱和 $u_x$ 分布；
- combined 能谱；
- tail 粒子每 cell 统计和噪声估计；
- conversion source spectrum；
- 左右边界通量。

还应输出阈值附近的表示接口诊断：

- conversion 前后的 combined 能谱；
- $K_{out}$ 上下各两个能量 bin 的粒子数和能量；
- bulk 与 tail 的 $J_x$、$\Pi_{xx}$、$\Pi_\perp$；
- 压缩前后六个约束矩的误差；
- `compression_fallback_count` 和 splitting/merging 次数。

碰撞开启时还要输出：

- 当前 tail collision kernel、pair mask 和 weight mode；
- 每个 cell 的候选粒子数、实际 pair 数和空 cell 比例；
- `max_s12`、大角散射比例、collision substeps；
- 等权重/不等权重 pair 数、virtual split 数和由此增加的粒子数；
- tail--tail 单独的 $\Delta P/\Delta K$ 残差；
- tail--bulk 的实际 tail 增量与 bulk 反作用残差；
- SDE 一、二阶样本矩或 Nanbu--Perez 散射角统计摘要。

combined 能谱必须使用统一的能量 bin，不得把两种表示的图简单叠加但使用不同归一化。

**[已完成-保留]**（阶段 H6，碰撞相关项待 H7/H8）：快照新增/补齐——
`energy_spectrum` 增加 `count_tail/count_combined` 与对应 per-eV 列
（统一 log 能量 bin）；`momentum_distribution_upar` 增加 tail 列；
新增 `tail_per_cell_stats_rank*.dat`（逐 cell 宏粒子数/权重和/密度 +
weight_stats 噪声估计）；新增 `tail_threshold_interface_rank*.dat`
（bulk/tail 的 Jx/Πxx/Π⊥、K_out±2 细能量 bin 的粒子数与能量、
最后接受步 conversion source 谱、最后控制步 groups/fallbacks）；
manifest 列示 tail 配置与新增文件。conversion source 谱由
`BulkTailConversionDiagnostics` 在转换扫描内填充（64 log bin，
[K_out, max_conversion_energy]），快照写最后接受步。

**[2026-08-08 最终复核]**：全局 MPI 归约、显式 bin 边界、多 location 累计和
`pre_extraction/removed/created/accepted` 四阶段谱已经在
`hybrid_h9_threshold_source_beam_12fs` 中通过。30080 条事件-bin 数据完整，
不再出现 `conversion_source_bins 0`；源数目链相对累计误差约
$3.5\times10^{-17}$。该诊断项恢复为 **PASS**。

**[配置溯源结构已完成；本地执行校验通过，多 rank 待验收]**：snapshot/checkpoint manifest 由请求配置生成，记录 field/background boundary、collision model、请求 kernel、tail--tail/tail--bulk 后端、四通道 pair mask、weight mode、weight algorithm version、子步上限和 population controller；已删除硬编码 `tail_collision_backend=none_v1`。其中兼容 token `virtual-split` 必须同时写出 `tail_collision_weight_algorithm=sentoku-kemp-bounded-v1`，否则视为旧实现或来源不明。阈值接口同时写出每个细 bin 的左右边界、宽度及归一化 $dN/dK$，并保留 forbidden bulk 数目/能量。旧 H9 曾出现 manifest pair mask 为 1 而 Beam applied flags 全零；统一碰撞入口后，本地 Beam+非零 tail 生产回归已要求并得到四个 applied flags 全为 1。manifest 仍只表示“请求配置”，H9 验收必须同时检查接受态 applied ledger，不能仅凭 manifest 判定碰撞执行。

### 13.3 性能诊断

在 `VpfpStepResult` 中增加：

```text
wall_bulk_x
wall_bulk_u
wall_tail_push
wall_tail_deposit
wall_tail_migrate
wall_conversion
wall_beam
wall_poisson
wall_collision_bulk
wall_collision_tail
wall_collision_tail_binning
wall_collision_tail_pairing
wall_collision_tail_scatter
wall_collision_reaction
wall_diagnostics
tail_particles_local_max
```

MPI 汇总只在接受步末一次完成。生产 `diagnostic-level=1` 不得在每个子步做全局粒子排序或全谱统计。

**[已完成-保留]**（阶段 H6）：`VpfpStepResult` 已增加
`wall_tail_push/deposit/migrate/conversion/diagnostics_seconds` 与
`tail_particles_local_max`（MPI 全局最大）；`wall_tail_migrate` 在
`BackgroundTailPIC::exchange_carriers` 内计时。碰撞相关计时
（wall_collision_tail_*）待 H8。

---

## 14. 文件级修改清单

### 14.1 `src/parameters.h`

1. 保留 $L_x=40\ \mu\mathrm m$、$N_x=8000$、$\Delta x=0.005\ \mu\mathrm m$。
2. 删除“生产必须无限扩大 `u_parallel` 域”的注释和默认策略。
3. 增加 `HybridTailParameters`，不要继续堆叠松散全局常量。
4. 参数包含阈值、转换 bin、粒子上限、返回模式和诊断门。
5. 所有运行时可变参数进入 checkpoint manifest。

### 14.2 `src/grid.h`

1. 新增 `HybridVelocityPartition`。
2. 初始化时预计算离散能量、resolved/conversion/buffer 掩码和矩权重。
3. 提供只读 accessor，禁止在内核中重复分支计算能量阈值。
4. 网格构造必须验证阈值位于 resolved 区且远离外边界。
5. 阈值支撑审计必须能逐 `u_perp` 环查询 $K_{out}$ 上下最近的
   `u_parallel` cell center；禁止只提供把所有环混合后的全局中心列表。

### 14.3 `src/species.h/.cpp`

1. 保留 cell-integrated mass 语义。
2. 提供转换器需要的安全局部访问或批量 extraction API。
3. extraction API 接收明确的 mass 数组并返回实际扣除矩，不允许外部直接修改私有数组后猜测账本。
4. 增加 forbidden-region mass 的单次融合扫描。

### 14.4 `src/background_tail_pic.*`

按第 6 节完整新增。不要混入 Beam injection 代码。

可抽取 Beam 中已验证的 shifted face-field gather、CIC deposit、轨迹分段和邻居 MPI 迁移为无注入语义的公共辅助模块，但必须先用 Beam 回归证明抽取前后结果一致。禁止通过继承 `BeamPIC` 共享状态。

### 14.5 `src/bulk_tail_converter.*`

按第 7 节完整新增。所有转换矩和加载逻辑只能有一份生产实现；测试直接调用该实现，不复写公式。

当前 `THRESHOLD_AWARE_COMPRESSION` 仅视为候选策略。需补充：

1. `BulkTailConversionDiagnostics` 携带显式谱边界、conversion location、
   pre-extraction/removed/created 四阶段 $N/K$；
2. 多次 conversion location 的诊断可加和，不能以后一次覆盖前一次；
3. `set_loading_policy()` 在测试和生产 CLI/manifest 中可追溯；独立门完成前不得
   仅靠构造函数默认值宣称生产策略已验收；
4. 压缩算法保持现有六矩闭合，不在本轮加入谱平滑或全局矩补丁。

最新源谱结论将文件级工作更新为：

5. 新增 `src/tail_subcell_quadrature.*`，统一提供 cell-volume 能谱求积和 converter
   候选支撑；诊断和生产加载不得复制两套子单元公式；
6. `BulkTailConverter` 第一候选只改变单元内 PIC 支撑，不改变 bulk extraction
   cell 集合、时间层、事务和全局账本；
7. 子单元候选权重继续调用 `TailMomentConstraint` 做非负六矩配对；fallback 必须
   显式记账；
8. 不实现每步重复比例抽取的静态 cut-cell；需要部分单元转换时必须使用保存的
   子单元状态或真实跨界通量。

#### 14.5.1 `src/tail_population_controller.*`

按第 7.10 节独立实现。该类只能接收 tail trial 粒子并返回新的 tail trial；不得访问 Poisson、Beam、bulk remap 或 accepted 状态。矩约束构造与 converter 共用一个小型 `TailMomentConstraint` 模块，不能复制两套公式。

### 14.6 `src/maxwell.*`

扩展电荷组装接口，算法保持非周期 Gauss/Poisson。

### 14.7 `src/vpfp_integrator.*`

1. 构造函数接收 `BackgroundTailPIC` 所需配置和 `BulkTailConverter`。
2. 增加 tail accepted/trial 工作区。
3. 按第 8 节重排无碰撞步骤。
4. 碰撞开启后按第 10.4 节调度。
5. 接受前统一计算 combined ledgers。
6. 失败路径不得改变传入的 accepted 对象。
7. 预分配所有 bulk 和 tail 工作数组，禁止每步重复大数组分配。
8. 当前真实碰撞接入时，禁止继续对传入的 accepted `Species& electrons` 原地调用 `collision_.apply()`；两个碰撞半步都必须作用于 trial bulk，失败后完整丢弃 trial。
9. 将一个接受步内各 conversion location 的 source spectrum 累计到
   `VpfpStepResult`；只累计实际接受态，失败 trial 必须丢弃。

### 14.8 `src/main_vpfp.cpp`

新增 CLI：

```text
--background-tail-mode <off|pic>
--tail-convert-energy-mev <value>
--tail-conversion-upar-bins <value>
--tail-conversion-energy-bins <value>
--tail-max-particle-count <value>
--tail-max-number-fraction <value>
--tail-target-particles-per-bin <value>
--tail-max-particles-per-bin <value>
--tail-population-control-interval <value>
--tail-max-weight-ratio <value>
--tail-return-mode <none|hysteretic>
--tail-collision-kernel <none|coulomb-nanbu-perez|kramers-moyal-sde|trace-stationary-background>
--tail-collision-weight-mode <equal-strata|virtual-split>
--tail-collision-max-substeps <value>
--tail-collision-max-particle-growth <value>
```

`virtual-split` 是向后兼容 token，当前唯一生产实现为
`sentoku-kemp-bounded-v1`；`probabilistic-ab` 尚未实现，不得写入生产命令。

启动时打印最终配置并写入 manifest。未知参数必须失败，不能静默使用默认值。
人口控制参数（`--tail-population-control-interval`/
`--tail-target-particles-per-bin`/`--tail-max-particles-per-bin`/
`--tail-max-weight-ratio`）在阶段 H5 起生效（interval>0 且 tail=pic 时
激活 TailPopulationController）；`--tail-max-weight-ratio` 为 H5 新增
（默认 8.0，近等权分裂的权重比上限，>1 才合法）。

### 14.9 `src/vpfp_checkpoint.*`

按第 12 节升级 schema。先完成 round-trip 测试，再允许生产保存混合 checkpoint。

**[已完成-保留]**（阶段 H6）：schema v2（tail 状态/config/账本/manifest，
见 §12.1）；新增 `VpfpCheckpointTailState/VpfpCheckpointTailConfig` 与
`tests/checkpoint_roundtrip_test`、`tests/checkpoint_restart_equivalence_test`
（§12.2/12.3 验收）；`main_vpfp.cpp` 的 checkpoint 写/读已接线 tail 状态
与配置校验。

### 14.10 `src/vpfp_diagnostics.*`

按第 13 节输出分表示和 combined 结果。生产模式只写已接受态。

**[已完成-保留]**（阶段 H6）：快照 tail 谱/逐 cell 统计/阈值界面/
conversion source 谱与 §13.3 时间列已实现（详见 §13 记录）；碰撞诊断
列待 H7/H8。

**[已完成-保留]**：conversion source 谱已在 rank 0 写出前完成全局归约，并
输出 `pre_extraction/removed/created/accepted` 四阶段谱和显式 bin 边界。
详细扫描仍只允许在 `diagnostic-level=2` 和指定快照执行，禁止进入每步生产热
路径。当前阻断项已转移到有限体积单元内部的阈值支撑与守恒 PIC 加载。

### 14.10.1 `tests/hybrid_threshold_support_audit_test.cpp`

扩展逐 `u_perp` 环和占据加权审计，至少覆盖 `k=0..7`。输出阈值上下最近中心、
跨越空隙和生产首次转换态的局部占据，避免全环聚合掩盖近轴空隙。

### 14.10.2 `tests/bulk_tail_threshold_interface*_test.cpp`

加入 `near-axis-narrow`；对 golden/current/threshold-aware 分别定义真实 PASS
门；结果使用截断写入；MPI 测试必须实际运行 1/2/5 rank 并比较逐 bin 全局谱。
当前 `threshold-aware` 不得因策略名不是 golden 而跳过谱误差检查。

### 14.10.3 单元体积谱与子单元加载测试

新增：

```text
tests/bulk_tail_cell_volume_spectrum_test.cpp
tests/bulk_tail_subcell_loading_test.cpp
```

前者只验证真实 cell faces、圆柱测度、0.05/0.1/0.2 MeV 能谱和 $N/K$ 求积；
后者调用生产 converter 的子单元加载路径，验证非负权重、六矩、Poisson、事务、
粒子数和 MPI 确定性。两个测试不得复写生产求积或矩约束公式。

### 14.11 `src/collision_coefficients.*` 与 `src/cylindrical_fp_collision.*`

无碰撞混合阶段不修改生产公式。待混合无碰撞验收后，再按第 10 节重构，避免同时调试表示转换和碰撞。

新增第 10.3 节列出的 `background_tail_collision.*`、`background_tail_nanbu_perez.*` 和 `background_tail_collision_sde.*`。`background_tail_nanbu_perez` 只拥有 tail PIC 内部或 PIC--PIC 二体散射；它不得直接访问 Eulerian bulk 数组。跨表示 pair 由 `HybridCollisionStep` 统一调度。

**[已完成-保留]**（阶段 H7）：`collision_coefficients.*` 增加
`CollisionCoefficientMode` 契约与 `MomentClosureCollisionCoefficients`
（moment-closure）；`cylindrical_fp_collision.cpp` 重构为逐速度点系数、
非对角扩散（显式守恒交叉步）、Chang-Cooper 上风格式修正（原 delta 权重
反置导致 pe~2 时负质量，已修正为平衡保持形式）、碰撞子步
（$\nu\Delta t\gg1$ 时按 §10.3.5 子循环，系数固定在同一时间层）、
储库能量/动量显式记账；积分器首个碰撞半步作用于 trial 副本（§14.7-8，
拒绝步不改 accepted）。tail 碰撞后端与 `HybridCollisionStep` 属 H8。

#### 14.11.1 `src/hybrid_collision_step.*`

H8 新增。它拥有 pair mask 和 collision trial 调度，不拥有 accepted 状态。`VpfpIntegrator` 只调用该统一入口，不能分别在外层调用 bulk collision 和 tail collision。

**[已完成-保留]**（阶段 H8）：`HybridCollisionStep::advance(bulk,tail,dt)`
在同一碰撞半步内完成 pair-masked bulk FP + tail 后端 + tail--bulk 反作用，
失败回滚双 trial；`VpfpIntegrator` 两个 Strang 碰撞半步在 tail 碰撞激活
时调用该入口（§15 H8 记录）。

### 14.12 `CMakeLists.txt`

新增源文件和测试 target。每一阶段只启用已完成的 target，禁止用 glob 自动把未完成源文件加入生产可执行文件。

---

## 15. 分阶段实施顺序

### 15.1 权威状态与使用规则

**权威状态（2026-08-13）**：基础阶段1--4和混合阶段H1--H9均已完成，Gate A/B曾作为
40 fs生产准入门通过；Gate C已连续推进到120 fs且数值稳定性/转换守恒通过。真实长跑
显示低能Tail长期滞留和宏粒子增长明显超过前期外推；H10的R0--R4现已通过，准生产中央档为
$K_{in}=5.5$ MeV、$N_{res}=8$。H11最终物理验收仍未签字。第15节只记录阶段目标、最终实现和验收结论；
全部命令统一见§17。

Gate A/B已经完成，不得重新打开H1--H9。当前生产路径固定为
`DIRICHLET_PHI(0,0)`、开放Beam、reservoir背景、`flux-interface`转换、
Chang--Cooper碰撞面通量和controller-off。

### 15.2 阶段状态总表

| 阶段 | 状态 | 最终结论 | 命令位置 |
|---|---|---|---|
| 1--4 | 已完成 | 开放拓扑、非周期Poisson和基础VPFP可用 | §17.1--17.2 |
| H1 | 已完成 | Tail PIC基础设施通过 | §17.3 |
| H2 | 已完成 | flux-interface转换与事务闭合通过 | §17.4、§17.10.2 |
| H3 | 已完成 | 无Beam混合与Gauss门通过 | §17.5 |
| H4 | 已完成 | Beam开放耦合通过 | §17.6 |
| H5 | 已完成 | 开发级阈值、分箱和转换门通过 | §17.7 |
| H6 | 已完成 | checkpoint、诊断和基础性能门通过 | §17.8 |
| H7 | 已完成 | bulk FP碰撞通过 | §17.9 |
| H8 | 已完成 | tail碰撞和反作用通过 | §17.9 |
| H9 | 已完成 | 40 fs controller-off生产门通过 | §17.10 |
| H10 | 短程通过/长期物理未通过 | R0--R4通过；100.3--120 fs资源收益通过，但场、能谱、动量和能量门未关闭 | §0.12、§16、§17.15 |
| H11 | 连续推进完成/物理验收待定 | 3127步全部接受；待能量边界账和Tail表示闭环 | §17.12、§19.4 |

H5已完成的是开发级验证；§17.11 Gate A是40 fs真实combined谱的生产准入审计。

### 15.3 基础阶段1--4 [已完成]

- 开放背景/Beam粒子拓扑和非周期静电Poisson场；
- 生产默认`DIRICHLET_PHI`且$\phi_L=\phi_R=0$；
- Eulerian圆柱速度坐标背景主体；
- 空间网格$\Delta x=0.005\ \mu\mathrm{m}$；
- reservoir背景边界和开放出流；
- 已验证的守恒有限体积、Poisson、Beam推进和checkpoint部件。

改变上述边界、网格或场模型属于新的物理模型变更，不能直接沿用H1--H9结论。

### 15.4 阶段H1：tail PIC基础设施 [已完成]

已完成Tail粒子状态、确定性ID/RNG、pusher、形函数沉积、开放边界、MPI迁移和矩诊断。
单粒子轨迹、沉积守恒、并行边界和staggered场插值通过。

### 15.5 阶段H2：保守bulk--tail转换 [已完成]

H2建立转换事务、六矩账本和失败回滚。真实审计否决`static-cell`/subcell路线，最终
17A--17F完成通量式转换：parcel与最终受限面通量同源，bulk sink、parcel source和
tail创建共用事务；非负权、tail-owned、MPI、checkpoint和碰撞面通量门均通过。

H2完成状态对应`flux-interface`，不代表旧static-cell方案仍受支持。

### 15.6 阶段H3：无Beam混合Vlasov--Poisson [已完成]

验证combined电荷进入同一个Poisson右端、tail-off退化、reservoir平衡和Gauss残差。
历史`build_phase4`已结束职责，不再要求重建。

### 15.7 阶段H4：Beam与无碰撞混合推进 [已完成]

Beam左侧注入、开放出流；bulk/tail/Beam电荷进入同一非周期Poisson系统；转换保持
接受步事务性。无碰撞短跑和组合电荷门通过。

### 15.8 阶段H5：阈值、分箱和有限bulk速度域 [已完成]

开发级验收覆盖阈值掩码、有限体积速度权重、分箱、近轴支撑、转换六矩和MPI一致性。
真实12 fs源谱审计推动flux-interface替代static-cell。

6 MeV附近0.1/0.2/0.4 MeV真实combined谱检查属于Gate A，不重新打开H5。

### 15.9 阶段H6：checkpoint、诊断和性能基础 [已完成]

checkpoint覆盖bulk、tail、Beam、场、ID/RNG和转换/碰撞账本；restart、一步等价、
失败事务和manifest往返通过。重启必须使用完整checkpoint，不能使用snapshot。

H6性能门不替代§18.3的120 fs资源预算。

### 15.10 阶段H7：bulk FP碰撞 [已完成]

bulk碰撞采用守恒有限体积/Chang--Cooper面通量；`exporting-absorbing`把越过表示
边界的最终碰撞通量交给tail。数目、动量、能量、零通量和平衡回归通过。

禁止恢复碰撞后静态扫描并清空高能cell的旧转换路径。

### 15.11 阶段H8：tail PIC碰撞与反作用 [已完成]

已完成tail--tail Nanbu--Perez、tail--bulk Kramers--Moyal SDE、bulk反作用、
有界Sentoku--Kemp不等权处理、近轴秩亏、子循环和全rank失败共识。生产永久宏粒子
增长预算为0，失败步骤事务回滚。旧人口爆炸checkpoint不得续跑。

已通过碰撞pair、反作用闭合、事务回滚、近轴秩亏、1/2/5 rank和Beam/no-Beam路径。
旧`virtual-split`名称只为CLI/checkpoint兼容保留，不代表仍创建永久residual粒子。

### 15.12 阶段H9：碰撞混合短跑 [已完成]

H9完成前置门、Beam 12 fs、restart、有源转换、Beam 3 fs、no-Beam 40 fs、25.5 fs
专项和Beam 40 fs controller-off。40 fs为1564/1564接受、0拒绝，
`h9_production.result`为`status=PASS`。

- conversion残差$\le3.62\times10^{-14}$；
- tail数量账误差$\le3.75\times10^{-14}$；
- Gauss残差$\le1.14\times10^{-13}$；
- 全局tail最大7,829,024，单rank最大900,184；
- 后100步最大2.446 s/step。

H9的碰撞接线、MPI控制流、tail-owned账本、转换事务和40 fs推进门均已关闭。
population controller继续关闭。

### 15.13 阶段H10：PIC-to-bulk双向返回 [短程通过/长期物理门未通过]

40 fs checkpoint 的 Gate B 指标曾支持 `h10_required=0`，该结论足以允许一次
controller-off Gate C，不是120 fs最终结论。实际120 fs中，Tail物理数分数仍仅
$1.95\times10^{-5}$，但Tail内部低于6 MeV的数目比例已达58.80%，能量比例达
42.21%，宏粒子达到$1.43165344\times10^8$，单步成本增长到约12 s。故低能Tail驻留已从可选机理诊断变为H10
设计依据。H10采用§16的迟滞和驻留条件；R0--R3已证明返回硬守恒、短程combined谱和场
扰动受控，并使Tail内核成本下降41.66%。R4进一步证明相邻阈值和驻留步数档保持同一宏观
波形平台，因此$K_{in}=5.5$ MeV、$N_{res}=8$已作为准生产中央档；Population Controller
仍保持关闭。但100.3--120 fs准生产显示资源成本下降41.6%的同时，末态场、能谱和动量分布
出现系统性漂移。故该中央档只保留为诊断候选，不能作为最终生产签字配置。详见§0.12和
§17.15。

### 15.14 阶段H11：120 fs生产 [连续推进完成/最终验收待定]

H11已沿用H9算法从40 fs推进到120 fs，3127步全部接受，未发生split或rollback。
转换守恒和阈值谱连续性通过；实际耗时5.99 h，Tail宏粒子达到1.4317亿。该结果是有效的
controller-off参考基线，但不能作为最终推荐生产配置。完成§19.4能量边界账审计和H10
长期表示误差定位后，才能决定是否需要生成新的最终生产结果。

### 15.15 当前唯一剩余执行顺序

1. 将现有40--120 fs结果冻结为controller-off参考基线，不重跑；
2. [已完成] 边界能量诊断接线和120 fs checkpoint单步验收；旧40--120 fs输出没有
   边界能量历史，不得事后伪造累计闭合；
3. [已完成] 对60/80/100/120 fs完整checkpoint执行低能Tail只读审计，确认低于
   6 MeV的Tail数目/能量占比到120 fs已达58.80%/42.21%，且主要位于内部移动波包；
4. [已完成] 带迟滞、可checkpoint驻留计数的Tail-to-bulk返回及R0--R3短程A/B；
5. [已完成] §16.6.5/§17.13.7的R4参数敏感性确认$K_{in}$和$N_{res}$存在稳定平台；
6. [已完成/部分失败] §17.14连续准生产的稳定性、守恒和资源门通过，长期物理等价与能量门
   未关闭；
7. [当前下一项] 完成能量恒等式和返回六矩有符号差诊断后，按§17.15从同一100.3 fs
   checkpoint执行none/hysteretic到105 fs的短A/B；
8. 只有短A/B能区分漏记账与真实推进误差，并给出可接受的宏观漂移，才允许扩展到111 fs。
   在此之前禁止新的0--120 fs生产；Population Controller继续关闭。

任一Gate失败即停止分析，不得通过修改阈值、放宽守恒门、恢复static-cell或开启
population controller绕过失败。

---
## 16. PIC-to-bulk 双向返回的条件方案

### 16.1 目标、边界和非目标

H10的目标是把已经碰撞减速、长时间停留在低能区的Tail电子从PIC表示返回
Eulerian bulk，同时保持combined物理状态。它的直接收益是降低Tail push、沉积、碰撞、
MPI迁移和checkpoint的成本；次要收益是降低已热化低能粒子的PIC采样噪声。

H10不是人口控制器，不允许合并后丢失物理数；不是能量补丁，不允许为闭合账本而改动
$E$或全局bulk；也不是边界清理器。尚未热化的束团、被俘获粒子、强各向异性粒子和
局部非负投影不可行的粒子必须继续留在PIC。

120 fs完整checkpoint审计已经触发H10开发门：$K<6$ MeV的Tail占Tail物理数的
58.80%、占Tail能量的42.21%，全部Tail达143,165,344个宏粒子，末期成本约12 s/step。

### 16.2 生产参数和迟滞状态机

新增参数：

```text
--tail-return-mode none|hysteretic
--tail-return-energy-mev <K_in>
--tail-return-residence-steps <N_res>
--tail-return-max-stencil-radius <R>
--tail-return-moment-tolerance <tol>
```

正向bulk-to-tail阈值仍是$K_{out}=$`--tail-convert-energy-mev`，必须满足
$0<K_{in}<K_{out}$。首轮A/B使用$K_{out}=6$ MeV、$K_{in}=5.5$ MeV、
$N_{res}=8$和$R=3$；这些只是候选值，不得在敏感性测试前标记为生产默认。

在`BackgroundTailParticle`中新增`std::uint32_t return_residence_steps`。它只在完整物理步的最终
试探态上更新一次：

- $K_p<K_{in}$：计数加1，在`UINT32_MAX`饱和；
- $K_p\ge K_{in}$：计数清0；
- 新的bulk-to-tail转换粒子：计数为0；
- 仅有计数$≥N_{res}$的粒子可进入投影候选集。

计数必须跟随粒子迁移、试探态复制、事务回滚和checkpoint/restart。失败步不得增加计数。

### 16.3 文件级实现要求

#### 16.3.1 粒子、checkpoint和CLI

1. 修改`src/background_tail_pic.h`，为`BackgroundTailParticle`增加驻留计数；所有粒子创建点必须
   显式初始化为0，不能依赖未初始化内存。
2. 修改`src/vpfp_checkpoint.h/.cpp`，把Tail粒子schema升级为
   `background_tail_particle_v2`，不再用`sizeof(BackgroundTailParticle)`直接写入未固定ABI的结构体；
   改用显式字段序列化，写入`x/ux/uy/uz/weight/id/return_residence_steps`。
3. 读取旧`background_tail_particle_v1`时将计数置0。旧checkpoint只能作为H10 A/B的共同
   初态；重启后必须经过$N_{res}$个新接受步才可返回，不得根据当前能量伪造历史驻留时间。
4. 修改`src/main_vpfp.cpp`，允许`--tail-return-mode hysteretic`，解析上述参数，写入启动日志、
   `physical_config_hash`和checkpoint manifest。新增
   `--restart-allow-return-config-change`，它只允许在H10的`none`与`hysteretic`配置之间建立受控
   双向A/B分支；除返回配置外的网格、场、Beam、碰撞、转换和MPI分区哈希仍必须严格相等。

#### 16.3.2 新增独立返回执行器

新增`src/tail_bulk_return.h/.cpp`，定义：

```cpp
struct TailBulkReturnConfig;
struct TailBulkReturnDiagnostics;
class TailBulkReturn;
```

`TailBulkReturn::apply()`的输入必须是`Species& bulk_trial`、`BackgroundTailPIC& tail_trial`、
`SpatialGrid`、`HybridVelocityPartition`、当前accepted step号和MPI信息。它不得访问或修改
`electrons`、`tail_accepted_`或已接受累计账本。

执行器按以下顺序实现：

1. 扫描当前rank的`tail_trial.particles`，更新驻留计数并建立候选索引；
2. 用现有`ParticleShape1D::cell_weights()`将每个候选粒子按CIC份额分到两个空间cell。
   这是表示转换，必须保持转换前Tail已有的cell密度沉积；不允许把全粒子权重归到
   单个最近cell；
3. 按共享同一`(cell0,cell1)`的候选粒子建立原子返回组。同组的两个CIC cell必须都能
   完成非负投影才能删除该组原粒子；任一cell不可行时整组保留。禁止只返回一个CIC份额
   却删除完整粒子，也禁止将剩余份额人为归一化；
4. 对MPI接缝处的第二个CIC cell，将目标矩和返回请求发送到该cell的所属rank。
   不得因为rank接缝而延迟、重复或遗漏返回；
5. 按原子组的两个`global_ix`分别构造局部目标矩，调用§16.4的非负投影；
6. 先将所有成功投影写入临时`ReturnMassRequest`列表，并建立待删粒子ID集合；
7. 执行本地有限性、非负性、重复ID和矩残差检查，然后用一次MPI全局共识决定是否提交；
8. 共识通过后，调用`Species::add_return_masses()`向`bulk_trial.f`加入非负cell-integrated mass，
   再从`tail_trial.particles`删除对应ID；共识失败时两者都不改变。

对某个cell的投影不可行是可预期的局部退化：该cell粒子原样保留并累计
`deferred_infeasible_groups`，不应使整步失败。只有NaN/Inf、负权重、重复删除、超预算加质量、
MPI请求不闭合或事务矩残差超限才是硬失败。

#### 16.3.3 接入生产Strang步

在`src/vpfp_integrator.h/.cpp`中持有`TailBulkReturn`和配置。两条生产路径
`advance_background()`与`advance_with_beam()`必须调用同一个私有辅助函数，不得复写算法。

返回只执行一次，位置固定为：

```text
第二个碰撞半步
-> 处理该半步的collision flux-interface bulk-to-tail导出
-> H10 tail-to-bulk返回
-> 最终Tail density/current沉积
-> 最终combined charge Poisson
-> 完整有限性/连续性/守恒门
-> accepted swap
```

不得在第一个碰撞半步后或Tail kick/drift中途返回，否则同一物理步会混用两套表示。
返回后必须重新计算bulk moments和Tail density，使最终Poisson使用返回后的combined电荷。
被返回粒子在本步Tail drift中已经产生的电荷守恒轨迹电流必须保留在本步
`tail_work_.current_face_x`中，不得因末态删除粒子而撤销或重算为零。返回只改变步末表示，
不改写本步已发生的轨迹历史。

### 16.4 非负局部速度投影

对每个空间cell，用CIC份额后的候选粒子构造完整审计矩：

$$
\mathbf b=(N,P_x,J_x,K,\Pi_{xx},\Pi_\perp)^T.
$$

速度支撑只能从`HybridVelocityPartition`标记为bulk-owned且$K_s<K_{out}$的Eulerian
$(u_\parallel,u_\perp)$ cell中选取。列向量必须复用当前`mass_cell_moments`和
`tail_particle_moments`的生产定义，禁止在H10里再写一套相对论速度/能量公式。

不得只以目标平均$u_\parallel$和$u_\perp$所在cell为中心构造单块stencil。真实checkpoint中
一个空间CIC组包含宽速度分布，单一均值附近小块只能表示均值，通常不能表示速度方差和
各向异性。执行器必须对每个CIC目标cell累计$u_\parallel/u_\perp$的加权均值、方差和实际范围，
在最小值、均值减标准差、均值、均值加标准差和最大值处建立去重锚点。每一档$R=1,2,3$
分别使用所有锚点附近$3\times3$、$5\times5$和$7\times7$小块的稀疏并集，而不是覆盖整个
速度矩形。支撑仍必须是bulk-owned且$K_s<K_{out}$的生产速度cell。

每组还要把粒子映射到最近生产速度cell，累计实际占用cell的候选质量。为控制成本，至多保留
24个确定性分桶锚点；少于24个时全部保留。跨MPI接缝请求必须携带相同的锚点索引和质量，
不能由接收rank重新猜测。非负求解器的prior优先使用该实际占用质量，并叠加仅用于数值启动的
极小正高斯权重，最后严格归一到目标$N$；prior仅用于降低迭代成本，不能作为最终答案。对每个
支撑并集调用`tail_solve_nonnegative_moment_weights()`。严格投影约束为$N$、$P_x$和$K$，
求得的$w_s$必须有限且$w_s\ge0$。这三个量分别保证转换瞬间电荷密度、平行动量和动能闭合。

$J_x$、$\Pi_{xx}$和$\Pi_\perp$仍使用同一最终权重重算并写入诊断，但不作为逐组$10^{-12}$
硬约束。连续PIC速度样本的完整六矩通常不位于固定cell-centre二维速度网格的同一个非负离散
凸包中，强制六矩精确会让真实组全部不可行。派生三矩必须通过速度网格收敛、combined能谱、
密度、Poisson场和短程A/B验收，不能被隐藏或用全局补丁抵消。全部验收使用未缩放物理矩。

`Species`新增`add_return_masses(const std::vector<ReturnMassRequest>&)`。它必须采用与
`extract_conversion_masses()`对称的预检—提交协议：先检查所有索引、mass有限非负和加和
溢出，全部有效后才修改`f`。

与正向转换不同，H10不需要为返回质量寻找“bulk正质量预算”：它是向bulk加入非负质量，
不会破坏非负性。真正的可行性条件是目标矩是否落在局部非负Eulerian支撑的凸包内。

### 16.5 账本、诊断和失败语义

新增`TailBulkReturnDiagnostics`并在`vpfp_step_diagnostics.dat`每个accepted step输出：

```text
tail_return_candidate_particles
tail_return_resident_particles
tail_return_attempted_groups
tail_return_committed_groups
tail_return_deferred_infeasible_groups
tail_return_deferred_rank_boundary_groups
tail_return_particles_removed
tail_return_number
tail_return_px
tail_return_jx_dx
tail_return_energy
tail_return_pixx_dx
tail_return_piperp_dx
tail_return_number_residual
tail_return_px_residual
tail_return_jx_residual
tail_return_energy_residual
tail_return_pixx_residual
tail_return_piperp_residual
tail_return_mpi_request_residual
tail_return_wall_seconds
```

`VpfpTailCumulativeLedger`增加返回的累计$N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$、删除宏粒子数和
退化组数，并写入checkpoint。正向conversion和反向return必须分列记账，不得仅写净值。

对每个成功返回组，定义相对残差

$$
r_q=\frac{|q_{bulk,added}-q_{tail,removed}|}
{\max(|q_{tail,removed}|,q_{L1},q_{floor})}.
$$

$q_{L1}$使用组内单粒子绝对贡献的和，避免$P_x$或$J_x$近零相消造成虚假大相对误差。
默认硬门为$N/P_x/K$最大$r_q\le10^{-12}$、MPI请求数目账$le10^{-13}$；
$J_x/\Pi_{xx}/\Pi_\perp$残差必须另记为`representation_residual_max`。若确定性求和精度包络
更宽，只能根据单元测试的网格收敛证据调整，不能为让短跑通过而放宽。

`tail-return-mode none`必须跳过所有候选扫描、MPI交换和诊断全局归约，保持现有二进制路径。

### 16.6 分阶段测试与验收标准

#### 16.6.1 阶段R0：`none`模式零回归

- 现有H1--H9单元测试全部保持PASS；
- 相同checkpoint、MPI rank数和参数下，`none`模式的bulk/Tail/Beam/field/RNG/ledger与
  H10接入前二进制状态一致；
- `tail_return_*`accepted量全部为0，且无新MPI collective。

#### 16.6.2 阶段R1：独立单元和MPI测试

必须新增以下真实调用生产实现的测试，不得在测试中复写投影公式：

- `tests/tail_bulk_return_hysteresis_test.cpp`：验证连续驻留、离开$K_{in}$清零、新粒子不立即返回；
- `tests/tail_bulk_return_parameter_switch_test.cpp`：验证从预热checkpoint分叉后，降低阈值会清零
  新阈值下不合格粒子的旧年龄、提高阈值时新增能带从1开始计数、改变$N_{res}$复用真实已保存年龄；
- `tests/tail_bulk_return_projection_test.cpp`：单cell、多cell和宽速度云的非负投影、不可行组保留；
- `tests/tail_bulk_return_transaction_test.cpp`：故障注入后bulk/Tail/RNG/ledger完全不变；
- `tests/tail_bulk_return_checkpoint_test.cpp`：驻留计数和累计账往返一致，v1粒子计数正确置0；
- `tests/tail_bulk_return_mpi_test.cpp`：1/2/5 rank下MPI接缝CIC返回、无重复ID、分区不变性。

所有测试必须输出`status=PASS`、$N/P_x/K$残差$le10^{-12}$、返回后bulk mass非负，且失败分支
的accepted状态哈希不变。

#### 16.6.3 阶段R2：完整checkpoint只读可行性审计

新增`tail_bulk_return_checkpoint_audit`，对60/80/100/120 fs状态在不修改checkpoint的前提下，分别统计：

- $K_{in}=5.25/5.5/5.75$ MeV下的候选物理数、能量和宏粒子数；
- $3\times3$、$5\times5$、$7\times7$ stencil的成功/不可行组数与可返回物理数；
- 每cell的$N/P_x/K$守恒残差、派生三矩表示残差、需要的最大stencil、核心/边界分布和预计删除宏粒子数；
- `audit_read_only_state_unchanged=1`。

只有在100和120 fs上，中位候选阈值下的可行物理数占候选数至少90%，且不可行组
没有系统集中在物理核心时，才进入动态A/B。若低于90%，先分析凸包不可行原因，不能扩大
全局stencil或放宽矩残差强行通过。

#### 16.6.4 阶段R3：固定checkpoint短程A/B

**状态：已通过。** 最终验收采用三路结构。R3a从同一个100 fs controller-off checkpoint运行
none/hysteretic到100.3 fs，验证首次返还；R3b物理A/B从同一个hysteretic预热后的100.3 fs
checkpoint分叉到101 fs；第三路从原始100 fs controller-off checkpoint到101 fs，只作为性能
基线。不得再用不同初态比较物理量，也不得要求粒子数近似相同的预热态A/B提供20%加速。

1. **硬守恒门**：全部步accepted，NaN/Inf、split、rollback和重复ID为0；每步返回$N/P_x/K$
   残差$≤10^{-12}$，Tail连续性、conversion、collision和Gauss保持现有硬门。
2. **物理不变门**：立即返回前后combined的$N/P_x/K$按各自L1尺度的相对差$≤10^{-12}$；
   $J_x/\Pi_{xx}/\Pi_\perp$差异必须显式报告并随速度网格加密收敛。动态R3a不是同一表示的
   bitwise replay：数千万Tail宏粒子返回bulk后，未约束高阶矩会使后续0.3 fs轨迹产生有限差异，
   因此不得继续用$10^{-9}$--$10^{-10}$逐点相等作为动态门。R3a采用：密度相对$L_2≤10^{-3}$、
   场相对$L_2≤2\times10^{-3}$、场能代理$\sum E^2$相对差$≤10^{-3}$、combined谱相对
   $L_2≤5\times10^{-3}$，以及按100个空间cell分块的场/密度RMS包络相对$L_2≤10^{-3}$。
   这些是短程宏观扰动上限；严格不变量仍由第一项硬门控制，不能以宏观门掩盖守恒失败。
   首轮动态A/B中，combined能谱在5--6.5 MeV不得产生新的单bin人工谷，电场能、combined动能和
   核心波包包络差异必须位于`none`重复运行的数值发散包络内。
3. **性能门**：B必须实际返回粒子，且末态Tail宏粒子数不高于独立controller-off性能基线。
   性能结论只用100.3--101 fs窗口判定；
   排除checkpoint和快照I/O后，Tail collision+push+deposit平均时间至少下降20%，否则H10虽可
   物理通过，但不值得进入新生产。

最终R3a和R3b均为`status=PASS`；R3b物理差异和硬守恒均低于上述门，独立性能门测得
`tail_kernel_wall_reduction=0.416586`。因此R3关闭，允许进入R4。

#### 16.6.5 阶段R4：参数敏感性和准生产门

**状态：已通过。**

只在R3通过后，比较$K_{in}=5.25/5.5/5.75$ MeV和$N_{res}=4/8/16$，每次只改一个参数。
所有档从R3a hysteretic分支的同一个100.3 fs预热checkpoint运行到102.3 fs，比较combined谱、
场能、压强矩、宏粒子数、返回率和wall time。禁止使用原始100 fs controller-off checkpoint，
避免重新引入驻留计数同步冷启动。

准生产参数必须处于敏感性平台：相邻参数档不得改变主波包的定性结构，combined能量和
5--6.5 MeV谱的差异不得呈单调的强参数驱动。若不存在平台，H10不能用于生产，保留
`tail_return_mode=none`基线。

最终五档均完整推进到102.3 fs，硬守恒门全部通过。原比较器错误复用了R3同参数A/B的场
相对$L_2\le2\times10^{-3}$，导致四个实际高度相关的参数分支被误判失败。R4是有意改变控制
参数的扫描，最终采用“有限L2包络+显式形状门”：场相对$L_2\le5\times10^{-3}$、相关系数
$\ge0.99998$、$\max|\Delta E|/\max|E|\le10^{-2}$；场能、密度、combined谱和硬守恒门保持原值。
现有结果的最差值分别为$3.357\times10^{-3}$、$0.9999944$、$8.078\times10^{-3}$、
$2.288\times10^{-4}$、$3.809\times10^{-4}$和$1.102\times10^{-3}$，均通过。中央档
$K_{in}=5.5$ MeV、$N_{res}=8$位于两组扫描内部，且平均Tail内核时间约3.893 s，登记为准生产配置。

---

## 17. 构建和运行命令

本节是全文唯一允许维护可执行命令的章节。其他章节只描述算法、验收标准和状态，并
引用本节。H1--H9 已全部完成；§17.3--§17.10.2 为完成记录和复现归档，默认不重跑。
§17.11.6的边界能量单步验收、低能Tail只读审计及H10 R0--R4均已完成；§17.14连续准生产
也已完成，但只通过稳定性、守恒和资源门，长期物理与能量门未通过。§17.13.7保留为R4复现
记录，不再默认重跑。Gate A/B和Gate C不再重复；下一实际入口是§17.15短A/B。

### 17.1 集群环境

工作目录：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
```

### 17.2 独立构建目录

重构期间不要覆盖已验证构建：

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4 --target fp_solver
```

旧 `build_phase4` 基线已完成历史职责，当前不再构建或作为生产前置条件。所有当前
测试和生产均使用同一 `./build/` 目录中的最新程序。

### 17.3 H1 测试 [已完成归档]

```bash
cmake --build build -j4 --target \
  background_tail_pusher_test \
  background_tail_deposition_test \
  background_tail_shape_difference_test \
  background_tail_open_boundary_mpi_test \
  background_tail_stagger_test

./build/background_tail_pusher_test
./build/background_tail_deposition_test
./build/background_tail_shape_difference_test
./build/background_tail_stagger_test
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/background_tail_open_boundary_mpi_test
```

多 rank 测试再分别使用 `-n 1` 和 `-n 5` 重跑，结果文件不得混写。

### 17.4 H2 测试 [已完成归档]

```bash
cmake --build build -j4 --target \
  bulk_tail_single_cell_test \
  bulk_tail_multibin_test \
  bulk_tail_poisson_invariance_test \
  bulk_tail_transaction_test

./build/bulk_tail_single_cell_test
./build/bulk_tail_multibin_test
./build/bulk_tail_poisson_invariance_test
./build/bulk_tail_transaction_test
```

`tail_population_controller_test` 不属于 H2：TailPopulationController 按
§7.10 在 H4 证明基础转换正确后才接入，已在 §15 阶段 H5 实现并加入显式
CMake target（`tail_population_controller_test` 与
`tail_population_controller_mpi_test`，构建/运行命令见 §17.7）。

### 17.5 H3 无 Beam 短跑 [已完成归档]

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --beam-enabled 0 \
  --collision-model none \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-return-mode none \
  --stop-time-fs 3 \
  --diagnostic-level 1 \
  --output-dir ./output/hybrid_h3_no_beam_3fs
```
**H3 验收状态：已完成。** 历史 `build_phase4` A/B 已完成其接线回归职责，
不再保留或重跑旧二进制命令。当前只保留上面的H3复现命令；它不属于生产前置门。



### 17.6 H4 Beam 短跑 [已完成归档]

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --beam-enabled 1 \
  --collision-model none \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 \
  --tail-return-mode none \
  --checkpoint-times 3,12,25,40 \
  --stop-time-fs 40 \
  --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_h4_beam_40fs
```

不要在 H1/H2 未通过时提交该作业。

单 rank 本地验收（H4 接线后）：

```bash
cmake --build build -j4 --target hybrid_beam_short_test
./build/hybrid_beam_short_test --case all
```

§15 阶段 H4 的 3/12/25/40 fs 验收按上述命令分别设置
`--stop-time-fs 3|12|25|40` 与对应 `--snapshot-times`（结果目录不混写）；
对照无转换热核时快照应与 `--background-tail-mode off` 跑一致到求和误差，
高能 tail 出现后检查 bulk 最外速度面占据率不再增长与 combined 能谱连续性。

### 17.7 阶段 H5 收敛矩阵与人口控制 A/B [已完成归档]

每种测试使用独立命令和目录，结果文件不得混写。先构建并运行 H5 测试
（单 rank 测试直接跑，MPI 一致性测试用 `yhrun -n 1 / -n 2 / -n 5`
分别重跑，`--result` 路径不混写）：

```bash
cmake --build build -j4 --target \
  tail_population_controller_test \
  tail_population_controller_mpi_test

./build/tail_population_controller_test --case all
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/tail_population_controller_mpi_test
```

**首轮 12 fs / 3 fs 矩阵结论（2026-08-04）**：全部步 accepted、无失败，
但物理时间落在转换前区（转换起点 ~step 820，t≈21 fs，见阶段 H4 记录），
`tail_particle_count` 全程为 0、人口控制从未触发、阈值/bin 差异落在同
配置运行间的 ULP 放大发散带内（U_E 相对差可达 ~1e-2），矩阵无效。以下
命令统一改为 40 fs，比较指标按“对运行间发散稳健的量”执行（见本节末尾
“比较内容”）。

注意：生产 `fp_solver` **不支持** `--overwrite-output`（未知参数直接
失败）；重跑前先清理目标目录：

```bash
rm -rf ./output/hybrid_threshold_4mev_40fs \
  ./output/hybrid_threshold_6mev_40fs \
  ./output/hybrid_threshold_8mev_40fs \
  ./output/hybrid_bin_1x1_6mev_40fs \
  ./output/hybrid_bin_4x4_6mev_40fs \
  ./output/hybrid_nv192_6mev_40fs \
  ./output/hybrid_nv384_6mev_40fs \
  ./output/hybrid_popctrl_off_6mev_40fs \
  ./output/hybrid_popctrl_on_6mev_40fs
```

**矩阵第 1 档：固定网格，$K_{out}$ 三档**（4/6/8 MeV，40 fs，转换 bin
4x4，人口控制 off）：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 4.0 --tail-return-mode none \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_threshold_4mev_40fs
```

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_threshold_6mev_40fs
```

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 8.0 --tail-return-mode none \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_threshold_8mev_40fs
```

**矩阵第 2 档：固定阈值（6 MeV），转换 bin 两档**（1x1 vs 4x4，40 fs，
人口控制 off）：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-upar-bins 1 --tail-conversion-energy-bins 1 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_bin_1x1_6mev_40fs
```

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-upar-bins 4 --tail-conversion-energy-bins 4 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_bin_4x4_6mev_40fs
```

**矩阵第 3 档：固定阈值（6 MeV）和 bin（4x4），bulk resolved 区分辨率
两档**（编译宏改 `FP_VELOCITY_GRID_NV`，40 fs）。基线即
`build`（Nv=192）；第二档在同一 `build` 目录串行重配置：

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DFP_VELOCITY_GRID_NV=384
cmake --build build -j4 --target fp_solver
```

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-upar-bins 4 --tail-conversion-energy-bins 4 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_nv384_6mev_40fs
```

（完成 Nv=384 后，须将同一 `build` 目录重新配置回 `FP_VELOCITY_GRID_NV=192`
并重新编译，再执行 Nv=192 对照跑；不得误用 Nv=384 二进制写入
`./output/hybrid_nv192_6mev_40fs`。）

**人口控制 off/on A/B**（固定 6 MeV、4x4 bin、40 fs；off 用
`--tail-population-control-interval 0`，on 用下面参数；阈值/谱/场/六个
约束矩的差异必须处于离散压缩误差内才能用于长跑）：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-upar-bins 4 --tail-conversion-energy-bins 4 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_popctrl_off_6mev_40fs
```

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --collision-model none --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-upar-bins 4 --tail-conversion-energy-bins 4 \
  --tail-population-control-interval 20 \
  --tail-target-particles-per-bin 64 \
  --tail-max-particles-per-bin 256 \
  --tail-max-weight-ratio 8.0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_popctrl_on_6mev_40fs
```

**比较内容**（不能以某档最接近 EPOCH 为唯一理由选择参数）：

- `vpfp_step_diagnostics.dat`：combined 密度相关列（N_e、U_E、K_e）、
  `tail_particle_count`（H5 新增列，宏粒子数）与 wall time；
- 快照目录：combined 能谱（`energy_spectrum`）、场与背景动能
  （`distribution_moments`）、tail 每 cell 统计；
- 六个约束矩 N/Px/Jx/K/Πxx/Π⊥ 与 Xw 由快照逐 rank 文件计算后 MPI
  汇总比较（人口控制 A/B 要求差异 ≤ 离散压缩误差 1e-10 量级）；
- `tail_population_control.dat`（on 档）核对每次控制的 groups/fallbacks
  与七矩残差；
- 40 fs 时间窗内转换约在 step 820 后开始，`tail_particle_count` 应在
  各阈值档出现可分辨差异；若某档仍为 0，说明该阈值在 40 fs 内未被触及，
  需延长或下调；
- 由于 Beam 开启时 OpenMP 归约 ULP 差异被束流-等离子体不稳定性指数
  放大（同配置运行间 U_E 相对差可达 ~1e-2，见阶段 H3/H4 记录），
  **禁止用单步场值或单个快照做档位判据**；改用对发散稳健的量：tail
  粒子数演化、时间窗平均场能/背景动能、late-time combined 能谱、总
  能量转移与性能。人口控制 off/on 还可对同一 accepted checkpoint 的
  tail 状态做离线同态 A/B（施加/不施加控制器后比较 combined 密度、
  场、谱与六矩），彻底消除双运行发散。

### 17.8 H6 checkpoint、诊断与 restart 测试 [已完成归档]

阶段 H6 的 schema-v2 checkpoint round-trip 与一步 restart 等价测试
（§12.2/12.3）。round-trip 测试单 rank 直接跑；restart 等价测试用
`yhrun -n 1 / -n 2 / -n 5` 分别重跑（`--workdir` 与 `--result` 各档
不混写）：

本节是H6通用回归，与17D的flux-conversion checkpoint测试相互补充但不能
互相替代。`checkpoint_restart_equivalence_test` 仅接受
`--case all --workdir <path> --result <path>`；尤其不接受 `--work-dir`。
若出现其usage后以exit code 2退出，应先检查调用的是否为错误可执行文件或
是否误用了17D专用命令的参数。

```bash
cmake --build build -j4 --target \
  checkpoint_roundtrip_test \
  checkpoint_restart_equivalence_test

./build/checkpoint_roundtrip_test
./build/checkpoint_restart_equivalence_test
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/checkpoint_restart_equivalence_test
```

生产混合跑 checkpoint 验证（tail=pic，短物理时间，写 checkpoint 后从
`--restart-dir` 继续，确认 step/time/累计账本/tail 粒子数与 next ID
衔接）：

```bash
rm -rf ./output/hybrid_h6_cp_check ./output/hybrid_h6_cp_restart

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --beam-enabled 1 --collision-model none \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --stop-time-fs 25 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 25 \
  --checkpoint-every 200 --checkpoint-dir ./output/hybrid_h6_cp_check \
  --output-dir ./output/hybrid_h6_cp_check

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --beam-enabled 1 --collision-model none \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 40 \
  --restart-dir ./output/hybrid_h6_cp_check \
  --output-dir ./output/hybrid_h6_cp_restart
```

注意：生产 `fp_solver` 不支持 `--overwrite-output`；重跑前先清理目录。

### 17.9 H7--H8 碰撞测试 [已完成归档]

在 H7/H8 前不得运行碰撞生产。完成后先执行：

**阶段 H7 范围**（tail collision=none，只验证 bulk 模块）：

```bash
cmake --build build -j4 --target cylindrical_fp_collision_test
./build/cylindrical_fp_collision_test --case all \
  --result ./output/cylindrical_fp_collision_test.txt
```

**阶段 H8+ 范围**（下述目标全部属于 H8，H7 不构建）：

```bash
cmake --build build -j4 --target \
  fp_solver \
  cylindrical_fp_collision_test \
  background_tail_nanbu_perez_test \
  background_tail_collision_moments_test \
  background_tail_collision_equilibrium_test \
  background_tail_collision_weight_test \
  tail_collision_isotropisation_test \
  tail_collision_two_population_test \
  hybrid_collision_pair_balance_test \
  hybrid_collision_transaction_test

rm -f ./output/h8_sk_bounded_v1_*.result

./build/background_tail_nanbu_perez_test \
  --result ./output/h8_sk_bounded_v1_nanbu_perez.result
./build/background_tail_collision_moments_test \
  --result ./output/h8_sk_bounded_v1_collision_moments.result
./build/background_tail_collision_equilibrium_test \
  --result ./output/h8_sk_bounded_v1_collision_equilibrium.result
./build/background_tail_collision_weight_test \
  --result ./output/h8_sk_bounded_v1_collision_weight.result
./build/tail_collision_isotropisation_test \
  --result ./output/h8_sk_bounded_v1_isotropisation.result
./build/tail_collision_two_population_test \
  --result ./output/h8_sk_bounded_v1_two_population.result
./build/hybrid_collision_pair_balance_test \
  --result ./output/h8_sk_bounded_v1_pair_balance.result
./build/hybrid_collision_transaction_test \
  --result ./output/h8_sk_bounded_v1_transaction.result
```

（`cylindrical_fp_collision_test` 是 H7 目标，已在 §17.9 上一段列出；
上面 8 个目标全部单 rank，`--result` 各自落盘、不混写。）

**生产 tail 碰撞短跑**（H8 验收，tail=pic + moment-closure + tail 后端，
短物理时间验证共存与账本）：

```bash
rm -rf ./output/hybrid_h8_sk_bounded_v1_tailcollision_3fs

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 3 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3 \
  --output-dir ./output/hybrid_h8_sk_bounded_v1_tailcollision_3fs
```

**生产 tail 碰撞短跑（SDE 后端）**（同样验收路径，验证
kramers-moyal-sde 与 bulk 反作用在生产共存）：

```bash
rm -rf ./output/hybrid_h8_sk_bounded_v1_tailcollision_sde_3fs

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/fp_solver \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel kramers-moyal-sde \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 3 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3 \
  --output-dir ./output/hybrid_h8_sk_bounded_v1_tailcollision_sde_3fs
```

碰撞测试分两组执行，不得用其中一组替代另一组：

1. `coulomb-nanbu-perez`：等权重二体事件的相对论 $P/K$ 闭合、散射角分布、Maxwellian 平衡、权重分层与有界 Sentoku--Kemp 不等权统计闭合；
2. `kramers-moyal-sde`：固定速度点的一、二阶增量矩、$dt/dt/2$ 弱收敛和 bulk--tail 反作用。

`trace-stationary-background` 只需与解析 stopping/pitch-angle benchmark 和 EPOCH 风格近似做对照，不计入完整生产通过标志。

**H8 集群验收结果**（2026-08-05，详见 §15 阶段 H8 记录）：

- 8 个 H8 单元测试全部 `pass=1`（output/ 根目录各 `--result` 落盘）：
  `background_tail_nanbu_perez_test`、`background_tail_collision_moments_test`、
  `background_tail_collision_equilibrium_test`、
  `background_tail_collision_weight_test`、`tail_collision_isotropisation_test`、
  `tail_collision_two_population_test`、`hybrid_collision_pair_balance_test`、
  `hybrid_collision_transaction_test`；
- `tail_collision_isotropisation_test` 采用 16 counter-based seed 均值判据
  （§10.3.6），实测 mean_A=0.0334、mean_B=0.0282、mean_C=0.0162、
  |A−B|=0.0052<0.02、C−A=−0.0172≤1e-3，P/K 守恒 1e-14；
- 两个生产短跑（`hybrid_h8_tailcollision_3fs` 与
  `hybrid_h8_tailcollision_sde_3fs`，80 rank）各 118/118 步全 accepted，
  collision_reservoir 30381→5425（累计 2.39e6）、K_e −1.92e6、U_E 1.48e5、
  tail 恒 0、电荷残差 ≤3.5e-15，指标与 §15 记录一致。

上述 2026-08-05 集群结果属于旧不等权实现。后续 `h8_sk_bounded_v1_*` 已用新
二进制完成重跑并通过，因此H8现为已完成；旧
`background_tail_collision_weight_test` 结果只作历史对照，不再形成待办。

### 17.10 H9 碰撞混合短跑 [已完成归档]

**归档说明**：本节记录17F之后已经完成的H9执行链。不要进入后文§17.10.1的
已暂停subcell路线，也不要重复提交这些长跑。历史执行顺序为：

1. 运行本节的生产前置门，重验不等权 pair、bulk--tail 反作用、事务与
   checkpoint 配置；
2. 提交 Beam 12 fs **无源稳定性档**；
3. 完成一个有明确 bulk-to-tail 转换源的独立触发档；
4. 仅在上述两档通过后，依次提交 Beam 3 fs 与 no-Beam 40 fs；
5. 再决定是否提交 Beam 40 fs controller-off。controller-on A/B 继续禁止。

**当前状态（2026-08-11）**：第1、2、2.5、2.6、3和4项均已通过。第4项的修复链和
最终证据如下：

1. tail-owned账本不再把MPI x-ghost副本作为独立物理质量；物理face转移、roundoff
   丢弃量和未解释余额分别记账；
2. 六矩支撑压缩改为固定宽度Caratheodory消元，秩亏重复求积点允许使用long-double
   守恒消元结果，但仍由原六矩残差硬门验收；
3. 所有conversion路径在进入tail drift collective前执行全rank结果共识，禁止
   rank-local失败造成控制流分叉；
4. 25.5 fs专项回归完成997步且80 rank最终均为 `step_accepted`；
5. Beam 40 fs controller-off完成1564步、0 rejected step，
   `h9_production.result` 的全部gate均为1且 `status=PASS`。

40 fs结果的最大conversion残差为 $3.62\times10^{-14}$，最大tail数量账误差为
$3.75\times10^{-14}$，最大Gauss电荷归一残差为 $1.14\times10^{-13}$；全局tail宏粒子
最大7,829,024，单rank最大900,184，后100步最大wall time为2.446 s/step。第4项不得再
标记为“待重跑”。population controller继续关闭；后续只验收统一能谱的阈值连续性和
120 fs资源外推，不再修改已经通过的转换守恒路径。

对所有碰撞开启的 `flux-interface` 生产命令，以下四项是当前程序的硬要求，
缺少其中任意一项会使运行被 CLI 拒绝或退回非生产语义：

```text
--bulk-collision-integrator chang-cooper-flux
--collision-interface-mode exporting-absorbing
--tail-conversion-mode flux-interface
--tail-flux-quadrature-order 4 --tail-flux-max-supports 7
```

#### 当前权威执行清单

以下命令是 17F 通过后的唯一执行顺序。工作目录固定为
`/HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver`；所有
可执行文件均位于 `./build/`，`fp_solver` 不支持 `--overwrite-output`。

**第 1 项：生产前置门 [已完成]。** 先编译，再依次运行全部求解器门与一项离线验收工具
fixture。任一命令返回非零或其 `.result` 不是 `status=PASS`，停止，不能提交 Beam
12 fs。

```bash
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target \
  fp_solver background_tail_collision_weight_test checkpoint_roundtrip_test \
  cylindrical_fp_collision_test bulk_tail_flux_loader_test \
  hybrid_collision_pair_balance_test hybrid_collision_transaction_test \
  hybrid_beam_short_test hybrid_no_beam_test \
  hybrid_collision_failure_consensus_mpi_test

rm -rf ./output/h9_flux_interface_17f_precheck
mkdir -p ./output/h9_flux_interface_17f_precheck

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/background_tail_collision_weight_test \
  --result ./output/h9_flux_interface_17f_precheck/weight.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/checkpoint_roundtrip_test \
  --result ./output/h9_flux_interface_17f_precheck/checkpoint.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/cylindrical_fp_collision_test --case roundoff-negative-cleanup \
  --result ./output/h9_flux_interface_17f_precheck/roundoff_negative.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_loader_test --case all \
  --result ./output/h9_flux_interface_17f_precheck/flux_loader.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/hybrid_collision_pair_balance_test \
  --result ./output/h9_flux_interface_17f_precheck/pair_balance.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/hybrid_collision_transaction_test \
  --result ./output/h9_flux_interface_17f_precheck/transaction.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/hybrid_beam_short_test --case beam-hybrid-collision-pairs \
  --result ./output/h9_flux_interface_17f_precheck/beam_pairs.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/hybrid_beam_short_test --case all \
  --result ./output/h9_flux_interface_17f_precheck/beam_all.result || exit 1
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/hybrid_no_beam_test --case all \
  --result ./output/h9_flux_interface_17f_precheck/nobeam_all.result || exit 1
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/hybrid_collision_failure_consensus_mpi_test \
  > ./output/h9_flux_interface_17f_precheck/failure_consensus.out 2>&1 || exit 1
python3 ./tests/h9_production_validator_test.py \
  > ./output/h9_flux_interface_17f_precheck/production_validator_fixture.result 2>&1 || exit 1

grep -H '^status=' ./output/h9_flux_interface_17f_precheck/*.result
grep -E 'accepted_min=0 accepted_max=0 failure_min=5 failure_max=5|status=PASS' \
  ./output/h9_flux_interface_17f_precheck/failure_consensus.out
grep '^status=PASS$' \
  ./output/h9_flux_interface_17f_precheck/production_validator_fixture.result
```

**第 2 项：Beam 12 fs 无源稳定性档 [已完成]。** 仅在第 1 项全部前置门通过后运行。当前
`DIRICHLET_PHI(0,0)`、6 MeV 阈值配置在 12 fs 的阈值快照中，$K\ge6\,\mathrm{MeV}$
bulk/tail 均为零，只有 $5.8$--$6.0\,\mathrm{MeV}$ 的约 $10^{-56}$ 舍入残量。因此本项
不得再要求非空 tail；它验证无真实源时不会由碰撞 roundoff、接口 loader 或 population
logic 凭空产生 tail。必须保存 10.5 fs 与 12 fs checkpoint 及 12 fs 阈值快照；10.5 fs
checkpoint 用于无源 restart A/B，12 fs checkpoint 用于最终状态保存。

```bash
rm -rf ./output/h9_flux_interface_17f_beam_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --checkpoint-times 10.5,12 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 12 \
  --output-dir ./output/h9_flux_interface_17f_beam_12fs
```

运行完成后先执行单次生产验收。该脚本只读取 accepted 输出；`status=PASS` 是进入
restart 对照和第 3 项的必要条件。

```bash
python3 ./tools/validate_h9_production.py \
  --run ./output/h9_flux_interface_17f_beam_12fs \
  --mode beam12 --min-accepted-steps 470 --require-tail no \
  --max-tail-particles 100000 --max-local-tail-particles 20000 \
  --max-last-wall-s 5 --require-threshold-snapshot \
  --result ./output/h9_flux_interface_17f_beam_12fs/h9_production.result
cat ./output/h9_flux_interface_17f_beam_12fs/h9_production.result
```

**第 2.5 项：无源 checkpoint/restart 对照 [已完成]。** 第 2 项通过后运行。它只验证完整
checkpoint 的物理配置、bulk 状态、空 tail、随机数和账本可恢复；不再把本项伪装为
“转换后非空 tail restart”。`CHECKPOINT_10P5` 必须指向第 2 项实际生成的 10.5 fs
checkpoint，这样 restart 的 10.5--12 fs 接受步能与直跑输出逐步重叠。checkpoint
配置文件名是 `manifest.txt`，不是快照的 `manifest.dat`。

```bash
CHECKPOINT_10P5=$(find ./output/h9_flux_interface_17f_beam_12fs -maxdepth 1 \
  -type d -name 'checkpoint_target10.5fs_*' | head -n 1)
test -n "$CHECKPOINT_10P5" || { echo 'missing 10.5 fs checkpoint'; exit 1; }
rm -rf ./output/h9_flux_interface_17f_beam_12fs_restart
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_10P5" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --checkpoint-times 12 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 12 \
  --output-dir ./output/h9_flux_interface_17f_beam_12fs_restart

python3 ./tools/validate_h9_production.py \
  --run ./output/h9_flux_interface_17f_beam_12fs_restart \
  --mode beam12 --min-accepted-steps 1 --require-tail no \
  --max-tail-particles 100000 --max-local-tail-particles 20000 \
  --max-last-wall-s 5 --require-threshold-snapshot \
  --result ./output/h9_flux_interface_17f_beam_12fs_restart/h9_production.result
# 两个变量都必须指向各自实际生成的 checkpoint 目录；同为 80 rank。
CHECKPOINT_12_DIRECT=$(find ./output/h9_flux_interface_17f_beam_12fs -maxdepth 1 \
  -type d -name 'checkpoint_target12fs_*' | head -n 1)
CHECKPOINT_12_RESTART=$(find ./output/h9_flux_interface_17f_beam_12fs_restart -maxdepth 1 -type d -name 'checkpoint_*' | head -n 1)
test -n "$CHECKPOINT_12_DIRECT" -a -n "$CHECKPOINT_12_RESTART" || \
  { echo 'missing 12 fs checkpoint'; exit 1; }
python3 ./tools/compare_h9_restart_production.py \
  --direct ./output/h9_flux_interface_17f_beam_12fs \
  --restart ./output/h9_flux_interface_17f_beam_12fs_restart \
  --min-common-steps 1 --rtol 1e-8 \
  --direct-checkpoint "$CHECKPOINT_12_DIRECT" \
  --restart-checkpoint "$CHECKPOINT_12_RESTART" --checkpoint-hash-informational \
  --result ./output/h9_flux_interface_17f_beam_12fs_restart/restart_compare.result
cat ./output/h9_flux_interface_17f_beam_12fs_restart/restart_compare.result
```

验收：第 2 和第 2.5 项的两个 `h9_production.result`、以及
`restart_compare.result` 均为 `status=PASS`。它覆盖 470/470 接受步、无
`vpfp_failure.dat`、`conversion_N=0` 时无碰撞诱发 tail 人口增长、Gauss、唯一界面、
`manifest.txt`、空 tail checkpoint 哈希诊断及 restart 后宏观状态一致性。二进制
checkpoint 哈希仅用于定位 OpenMP/MPI 归约或序列化路径的非逐字节差异，不是物理
restart 等价性的硬门；硬门使用守恒量、`gauss_charge_residual` 的绝对差和接受步序列。
若本项
出现非零 tail，或任一项不满足则停止。

**第 2.6 项：有源 bulk-to-tail 触发档 [已完成]。** 第 2 项的零源
结果不能验证真实转换、tail--bulk SDE 和其反作用。新增
`bulk_tail_flux_source_trigger_test`，它不改变生产 6 MeV 阈值、也不向生产 Beam
档人工注入粒子：测试在真实 `HybridVelocityPartition` 的阈值下 bulk donor cell 中
构造正、有限的四单元 PPM 平台，并施加常场使电子在一次生产
`ConservativePpmRemap::advect_u_parallel()` 中跨越真实的 outward
`u_parallel` interface face。该 remap 的 `apply_interface_sink=true` 会从 Eulerian
bulk 扣除面通量，再把同一 `BulkTailFluxBatch` 交给生产
`BulkTailConverter::convert_flux_batch()`。

该门的三个独立子契约分别覆盖：

- `bulk_tail_flux_source_trigger_test`：material face export、bulk 实际损失与
  face transfer 一致、非空 tail 创建，以及 $N,P_x,J_x,K,\Pi_{xx},\Pi_\perp$
  六矩闭合；
- `hybrid_beam_short_test --case beam-hybrid-collision-pairs`：非空 tail 进入
  生产 collision selector，四个 pair 与 bulk reaction 闭合；
- `bulk_tail_flux_checkpoint_test`：转换分区元数据、tail 粒子 ID 和累计账本的
  checkpoint round trip。

运行命令（单 rank，均为独立快速测试）：

```bash
mkdir -p ./output/h9_source_trigger_26
cmake --build ./build --target \
  bulk_tail_flux_source_trigger_test hybrid_beam_short_test \
  bulk_tail_flux_checkpoint_test -j 4 || exit 1

yhrun -n 1 --cpu-bind=cores ./build/bulk_tail_flux_source_trigger_test \
  --case all \
  --result ./output/h9_source_trigger_26/source_trigger.result || exit 1

yhrun -n 1 --cpu-bind=cores ./build/hybrid_beam_short_test \
  --case beam-hybrid-collision-pairs \
  --result ./output/h9_source_trigger_26/collision_pairs.result || exit 1

yhrun -n 1 --cpu-bind=cores ./build/bulk_tail_flux_checkpoint_test \
  --case all --workdir ./output/h9_source_trigger_26/checkpoint_roundtrip \
  --result ./output/h9_source_trigger_26/checkpoint.result || exit 1

python3 ./tools/validate_h9_source_trigger.py \
  --source ./output/h9_source_trigger_26/source_trigger.result \
  --collision ./output/h9_source_trigger_26/collision_pairs.result \
  --checkpoint ./output/h9_source_trigger_26/checkpoint.result \
  --result ./output/h9_source_trigger_26/h9_source_trigger.result
cat ./output/h9_source_trigger_26/h9_source_trigger.result
```

验收要求是 `h9_source_trigger.result` 的所有 `gate_*=1` 和 `status=PASS`。
该测试通过后，才允许把 `--require-tail yes` 用于具有已确认真实转换源的生产/检查点档。

**实际验收结果**：`output/h9_source_trigger_26/h9_source_trigger.result` 的全部
`gate_*=1` 且 `status=PASS`；material face export、bulk loss、parcel number、tail
number 的相对闭合误差为 $9.40\times10^{-15}$，六矩最大相对残差为
$3.41\times10^{-15}$。碰撞 pair/reaction 与 checkpoint round-trip 子门也均为
`status=PASS`。

**第 3 项：Beam 3 fs 与 no-Beam 40 fs [已完成]。** 两个命令都只能在第 2、2.5 和 2.6 项通过后
执行；可分别提交，但任一失败都不进入第 4 项。

```bash
rm -rf ./output/h9_flux_interface_17f_beam_3fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 3 --checkpoint-times 3 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3 \
  --output-dir ./output/h9_flux_interface_17f_beam_3fs

rm -rf ./output/h9_flux_interface_17f_nobeam_40fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 0 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --checkpoint-times 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 40 \
  --output-dir ./output/h9_flux_interface_17f_nobeam_40fs
```

```bash
python3 ./tools/validate_h9_production.py \
  --run ./output/h9_flux_interface_17f_beam_3fs \
  --mode beam3 --min-accepted-steps 118 \
  --result ./output/h9_flux_interface_17f_beam_3fs/h9_production.result
python3 ./tools/validate_h9_production.py \
  --run ./output/h9_flux_interface_17f_nobeam_40fs \
  --mode nobeam40 --min-accepted-steps 1564 \
  --max-tail-particles 0 --max-local-tail-particles 0 --max-last-wall-s 5 \
  --result ./output/h9_flux_interface_17f_nobeam_40fs/h9_production.result
grep -H '^status=' ./output/h9_flux_interface_17f_beam_3fs/h9_production.result \
  ./output/h9_flux_interface_17f_nobeam_40fs/h9_production.result
```

**实际验收结果**：两项均通过。Beam 3 fs 完成 118 个 accepted step、无 rejected
step，`max_gauss_charge_residual=4.44\times10^{-15}`；no-Beam 40 fs 完成 1564
个 accepted step、无 rejected step、tail 始终为零，
`max_gauss_charge_residual=5.42\times10^{-20}`。两份
`h9_production.result` 均为 `status=PASS`。

**第 4 项：Beam 40 fs controller-off [已通过]。** 第1--3项及本项均已通过。
早期复现曾定位并修复 tail-owned 账本重复统计 x ghost 的问题；该结论仍有效，但它不是
最新停滞的根因。2026-08-11 最新一批 25.5 fs trace 最后接受 step 为 954
($t=24.4081732924\,\mathrm{fs}$)。候选 step 955 的 80-rank trace 明确分裂为：

- rank 9 最后一行为 `conversion_begin`，本地 tail 粒子数为 57904；
- 其余 79 个 rank 最后一行为 `tail_drift2_begin`；
- 最新接受态的 `wall_conversion_s` 仅约 $0.017$--$0.019\,\mathrm{s}$，且 step 954
  全局只有 927 parcels、14832 quadrature nodes，因此不能把本次现象解释为正常的
  大分组压缩耗时；
- 目录中的 step 870 `vpfp_failure.dat` 和 `trial_bulk_tail_flux_failures.dat` 是前次作业
  遗留文件，不属于本次清空目录后的 step 955 候选态，不能继续作为本次根因证据。

该 trace 说明 rank 9 没有完成 conversion，而其他 rank 已进入含 MPI 粒子迁移的第二次
tail drift。结合前一 accepted step 的转换只需约 $0.02\,\mathrm{s}$，以及当前压缩器已
消除大分组超线性热点，最可能路径是 rank 9 遇到秩亏求积组后本地 conversion 返回失败，
直接离开 `advance()`；其他 rank 则在 drift collective 中等待。仅靠 trace 无法严格排除
rank 9 仍在 conversion 内部计算，因此修复同时覆盖两个分支：消除合法秩亏组的错误失败，
并在 conversion 后增加全 rank 结果共识。若仍失败，新共识日志将直接给出实际子原因，
不再以“长期卡住”隐藏失败。

生产修复包含两项，均不得放宽守恒或 fidelity 阈值：

1. `apply_upar_flux_conversion()` 以及第二 collision-half conversion 在进入后续 MPI
   drift 前必须调用统一的 `synchronize_conversion_outcome()`。任一 rank 失败时，全体
   rank 同步得到相同的 `failure_code`、`failing_rank` 和失败子原因，并在同一位置返回；
   禁止 rank-local early return。
2. `tail_moment_constraint.cpp` 的最后 full-row polish 允许秩亏输入回退到 long-double
   null-elimination 权重。重复 quadrature node 会使 full-row normal matrix 奇异，但这
   不等于物理矩不可守恒。回退后仍必须经过原有、未放宽的非缩放物理矩残差门；不得用
   负权、放宽阈值或增加 support 上限掩盖失败。

新增 `--tail-stage-trace` 只写每 rank 的本地 liveness marker，不加入任何物理计算、
接受条件或 MPI 归约。它只用于本次复现，正式生产必须关闭。

先编译并运行快速接口回归：

```bash
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target fp_solver upar_flux_sink_test \
  bulk_tail_flux_loader_test tail_moment_constraint_large_group_test \
  hybrid_collision_failure_consensus_mpi_test || exit 1
mkdir -p ./output/h9_flux_interface_17f_tail_owned_balance
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/upar_flux_sink_test \
  --result ./output/h9_flux_interface_17f_tail_owned_balance/upar_flux_sink.result || exit 1
cat ./output/h9_flux_interface_17f_tail_owned_balance/upar_flux_sink.result
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_loader_test --case all \
  --result ./output/h9_flux_interface_17f_tail_owned_balance/flux_loader.result || exit 1
cat ./output/h9_flux_interface_17f_tail_owned_balance/flux_loader.result
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_moment_constraint_large_group_test \
  > ./output/h9_flux_interface_17f_tail_owned_balance/large_group.result || exit 1
cat ./output/h9_flux_interface_17f_tail_owned_balance/large_group.result
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/hybrid_collision_failure_consensus_mpi_test \
  > ./output/h9_flux_interface_17f_tail_owned_balance/conversion_consensus.result || exit 1
cat ./output/h9_flux_interface_17f_tail_owned_balance/conversion_consensus.result
```

快速回归验收：`status=PASS`；
`guard_ledger_invariant=1`；
`tail_owned_expected_transfer_number` 与 `face_export_number` 相对差不超过
$10^{-12}$；`tail_owned_unexplained_relative\le10^{-12}$。该单元测试失败时不得提交
40 fs 作业。`flux_loader.result` 也必须为 `status=PASS`。大分组回归还必须给出
`status=PASS`、`active_supports<=7`、
`nonnegative=1` 且 `max_relative_residual<=1e-10`。新增的秩亏大分组还必须满足
`rank_deficient_compressed=1`、`rank_deficient_active_supports<=7` 和
`rank_deficient_max_relative_residual<=1e-10`。MPI 共识回归必须同时给出：

```text
conversion_ok_min=0 conversion_ok_max=0
conversion_failure_min=10 conversion_failure_max=10
conversion_rank_min=0 conversion_rank_max=0
status=PASS
```

它验证一个 rank 的 conversion failure 会在所有 rank 上一致拒绝，不会再让部分 rank
进入 tail drift collective。

2026-08-11 第二次短跑在约 $22.259\,\mathrm{fs}$ 的候选 step 870 再次触发
`upar_tail_interface`。该候选的 face export 和 expected transfer 均为
$4.0918\times10^7$，但旧统计另报 $1.0157\times10^7$ material residual。源码复核
确认：连通 tail 区域收支仍遍历了 MPI x-ghost；ghost 是相邻 rank 的通信副本，既没有
独立的物理 interface parcel，却被再次累计进 tail 余额。修复后只有
`ix in [nghost,nghost+nx_local)` 的物理单元进入 tail-owned 输入、roundoff 和残差
账本；ghost 输出仍按 one-way tail 语义清零，但不拥有独立守恒量。测试中的
`guard_ledger_invariant` 专门防止该错误回归。

快速回归通过后，先复现并跨过旧 step 958，而不是立即重跑到 40 fs：

```bash
rm -rf ./output/h9_flux_interface_17f_beam_25p5fs_tail_owned_fix
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-stage-trace \
  --stop-time-fs 25.5 --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h9_flux_interface_17f_beam_25p5fs_tail_owned_fix
```

验收：时间必须最终到达25.5 fs（当前固定步长下最后 accepted step 约为997），至少
跨过本次 step 955；不得出现所有 rank 长时间停在同一 accepted step。若 conversion
确实不满足物理硬门，程序必须快速、整齐地拒绝该步，并输出
`[bulk-tail-conversion-fail]`，其中包含准确的 `failing_rank`、`reason`、support/duplicate/
ledger 计数和六矩残差；不得再出现一个 rank 留在 `conversion_begin`、其余 rank 留在
`tail_drift2_begin` 的控制流分叉。成功跑到25.5 fs时，不得存在本次运行新产生的
`vpfp_failure.dat`；
`bulk_tail_flux_accepted_steps.dat` 中新的
`tail_owned_expected_transfer_number` 可为非零，而
`tail_owned_roundoff_discarded_number` 可在零源阶段为非零；
`tail_owned_bulk_residual` 必须保持在代码硬门内。若该短跑失败，先读取
`[tail-owned-balance-fail]` 的 `unexplained`、`expected_transfer` 与
`face_export`，禁止再次以全局 $N_{\rm bulk}$ 放宽阈值。

**实际结果：已通过。** 最后accepted step为997，时间为25.5 fs；80个rank最后均为
`step_accepted`，没有本次运行产生的failure文件。step 997的
`support_limit_violation_count`、`duplicate_id_count` 和
`face_ledger_mismatch_count` 均为0；六矩转换残差为 $O(10^{-15})$，
`tail_owned_bulk_residual` 约为 $5.19\times10^{-24}$。此前step 955的rank控制流分叉
没有复现。

阶段追踪文件为：

```text
output/h9_flux_interface_17f_beam_25p5fs_tail_owned_fix/tail_stage_trace_rank000000.dat
...
output/h9_flux_interface_17f_beam_25p5fs_tail_owned_fix/tail_stage_trace_rank000079.dat
```

列为：`step time_s rank stage local_tail_particles monotonic_s`。若某 rank 的最后一行是
`collision_half1_begin` 或 `collision_half2_begin`，该 rank 停在碰撞路径。若各 rank
最后阶段不同，必须先区分两种情况：同一阶段附近且时间戳持续更新属于负载不均衡；
单个 rank 停在 `conversion_begin`、其余 rank 已进入含 MPI 的 `tail_drift2_begin`，属于
转换失败后的控制流分叉，必须由转换共识处理，不能继续等待。25.5 fs复现已经通过，
随后完成了本项40 fs生产验收。

```bash
rm -rf ./output/h9_flux_interface_17f_beam_40fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-stage-trace \
  --stop-time-fs 40 --checkpoint-times 12,40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/h9_flux_interface_17f_beam_40fs
```

在作业停滞时，收集而不是等待无限久：

```bash
for f in ./output/h9_flux_interface_17f_beam_40fs/tail_stage_trace_rank*.dat; do
  tail -n 1 "$f"
done | sort -k4,4 -k1,1n
```

原运行命令如下，仅保留作历史对照：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 40 --checkpoint-times 12,40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/h9_flux_interface_17f_beam_40fs
```

```bash
python3 ./tools/validate_h9_production.py \
  --run ./output/h9_flux_interface_17f_beam_40fs \
  --mode beam40 --min-accepted-steps 1564 --require-tail yes \
  --max-tail-particles 12000000 --max-local-tail-particles 2000000 \
  --max-last-wall-s 5 --require-threshold-snapshot \
  --result ./output/h9_flux_interface_17f_beam_40fs/h9_production.result
cat ./output/h9_flux_interface_17f_beam_40fs/h9_production.result
```

**实际验收结果（2026-08-11）**：

- `accepted_step_count=1564`，`rejected_step_count=0`，最终时间40 fs；
- 80个rank最终均为step 1564的 `step_accepted`；
- `max_gauss_charge_residual=1.1368683772161603e-13`；
- `max_tail_number_balance_error=3.7537874674718047e-14`；
- `max_conversion_residual=3.6190277159198601e-14`；
- `tail_particles_final=tail_particles_max=7829024`；
- `tail_particles_local_max=900184`；
- `last100_wall_s_max=2.446282813`；
- 无collision诱发的零源tail增长，无static extractor调用，checkpoint和阈值快照齐全；
- validator全部gate为1，最终 `status=PASS`。

因此第4项生产推进、守恒、MPI控制流、既定资源预算和wall-time门均已通过。这里的
`gate_threshold_snapshot=1` 只表示阈值快照存在，不等于 §19.2 的combined能谱连续性
已经通过；该物理门仍需离线分析25 fs和40 fs快照。

H9的无Beam平衡、Beam 3 fs、12 fs和40 fs controller-off四档现均已完成并通过。
已执行顺序为：**本地/双rank前置门 → Beam 12 fs阻断门 → Beam 3 fs与无Beam
40 fs基线 → Beam 40 fs controller-off**。controller-on A/B不属于本次通过项，仍不得
为减少宏粒子而直接启用。四档全部使用
`moment-closure` bulk 碰撞 + tail PIC + `coulomb-nanbu-perez` tail 后端
（生产配置，与 H8 短跑同参数）；如需后端对比，可把
`--tail-collision-kernel` 换成 `kramers-moyal-sde` 后重跑同档目录对比。

**2026-08-06 H9 复核后的执行前置门**：旧 80-rank 四档均运行完成，但不能通过 H9。代码侧第 1--4 项现已完成，集群重跑前仍须先执行下方生产分支回归：

1. **[已完成]** 抽取唯一的 `apply_collision_half()`，使 `advance_background()` 与 `advance_with_beam()` 的两个 Strang 碰撞半步调用同一个 `HybridCollisionStep` 选择器；生产分支不再直接调用两次 bulk-only `collision_.apply()`；
2. **[已完成-本地通过]** `hybrid_beam_short_test --case beam-hybrid-collision-pairs` 以 Beam+非零 tail 直接驱动生产 `VpfpIntegrator::advance()`，要求 `pair_bb/pair_tt/pair_tb/pair_bt=1/1/1/1` 且 $C_{tb}+C_{bt}$ 反作用残差通过；
3. **[已完成-多 rank 待验收]** 接受态写出前对 conversion 的 $N/P_x/K$、创建粒子数及残差做 MPI 全局归约；tail outflow 与 combined 量沿用已有全局归约。rank-local checkpoint 累计量在全局化前更新，禁止把全局值重复写入每个 rank 的累计账本；
4. **[命令已统一]** 所有 H9 命令显式使用 `DIRICHLET_PHI`、$\phi_L=\phi_R=0$，并验收 manifest；
5. 明确本节验证的是 `moment_closure/BGK` 近似，不得将结果改名为 self-consistent Landau。

**2026-08-06 Beam 12 fs 卡顿根因与修复**：旧实现把不等权 pair 的剩余权重永久实体化为新宏粒子。在 `conversion_N=0` 后，tail 数仍从 step 413 的 10 增长到 step 444 的 380272350，单 rank 最大占 295290241；单步 collision wall time 同期从亚秒增长到约 282 s。这不是 MPI 死锁，也不是 SDE 子循环问题，而是每个碰撞半步最多翻倍的表示人口爆炸。当前代码已改为固定粒子数的有界 Sentoku--Kemp 修正，并把默认增长预算改为 0。旧输出目录以及已经保存了膨胀 tail 列表的 checkpoint 不得续跑；必须从首次 tail 转换前的 checkpoint 或从头以新二进制重跑 Beam 12 fs 门。

本轮代码侧预检结果（Windows 单 rank MPI stub，仅用于接口/事务/账本验收，不替代集群物理测试）：

- 完整 `fp_solver` 目标以 `-O2 -Wall -Wextra -fopenmp` 编译通过；
- `hybrid_collision_pair_balance_test`：四通道均被实际调用，$C_{tb}+C_{bt}$ 的 $P_x$ 相对余额为 $6.44\times10^{-16}$，能量余额为 0；
- `hybrid_collision_transaction_test`：失败态不提交，成功态反作用单元数非零；
- `background_tail_collision_weight_test`：有界 Sentoku--Kemp 修正保持粒子数固定，单事件能量误差为 0；1024 个独立不等权 pair 的集合动量相对误差为 $4.36\times10^{-4}$、集合能量误差为 $7.82\times10^{-15}$；零增长预算通过；
- `checkpoint_roundtrip_test`：weight mode、子步上限和碰撞配置往返一致；
- 一步 CLI 冒烟：snapshot manifest 与命令行一致，逐步诊断表头/数据均为 66 列。

统一后的新增本地结果如下：

- `hybrid_beam_short_test --case beam-hybrid-collision-pairs`：直接通过 Beam 生产分支推进一个非零 tail 状态，`bb/tt/tb/bt=1/1/1/1`，最大反作用相对残差 $2.48\times10^{-15}$；
- `hybrid_beam_short_test --case all`：无 tail 等价、Beam 转换连续性和碰撞 pair 回归同时通过；combined 数目相对误差 $2.80\times10^{-15}$，转换残差最大值 $2.34\times10^{-15}$；
- `hybrid_no_beam_test --case all`：共享碰撞入口重构后仍通过，最大账本差 $2.03\times10^{-16}$。

这些结果证明模块已经进入 Beam 生产分支，并且没有破坏无碰撞 Beam 或无 Beam 路径；仍需由下述 H9 四档验收 MPI 多 rank 的 applied flags、全局 conversion 账本、统计碰撞率、长期 tail 人口和阈值接口。

下方 `-n 80` 是当前 H9 配置。2026-08-06 旧四档 manifest 均为 80 rank，但使用了错误的 `left-E` 场边界且 Beam pair 未执行，因此必须使用新构建和新输出目录重跑，禁止在旧目录续写。

先构建生产程序和两条生产路径回归，并在登录/测试节点完成单 rank 门：

```bash
cmake --build build -j4 --target \
  fp_solver hybrid_beam_short_test hybrid_no_beam_test \
  hybrid_collision_failure_consensus_mpi_test \
  background_tail_collision_weight_test checkpoint_roundtrip_test \
  hybrid_collision_pair_balance_test hybrid_collision_transaction_test

rm -rf ./output/h9_flux_interface_17f_precheck
mkdir -p ./output/h9_flux_interface_17f_precheck

./build/background_tail_collision_weight_test \
  --result ./output/h9_flux_interface_17f_precheck/weight.result
./build/checkpoint_roundtrip_test \
  --result ./output/h9_flux_interface_17f_precheck/checkpoint.result
./build/hybrid_collision_pair_balance_test \
  --result ./output/h9_flux_interface_17f_precheck/pair_balance.result
./build/hybrid_collision_transaction_test \
  --result ./output/h9_flux_interface_17f_precheck/transaction.result

./build/hybrid_beam_short_test \
  --case beam-hybrid-collision-pairs \
  --result ./output/h9_flux_interface_17f_precheck/beam_pairs.result
./build/hybrid_beam_short_test --case all \
  --result ./output/h9_flux_interface_17f_precheck/beam_all.result
./build/hybrid_no_beam_test --case all \
  --result ./output/h9_flux_interface_17f_precheck/nobeam_all.result
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/hybrid_collision_failure_consensus_mpi_test
```

八项均须输出 `status=PASS`。其中 `hybrid_collision_pair_balance_test` 是 bulk--tail
反作用闭合门，`hybrid_collision_transaction_test` 是失败事务门；Beam pair 测试还必须包含
`collision_pairs_ok=1` 且 `collision_reaction_residual_max<=1e-10`。多 rank 失败共识测试必须
输出 `accepted_min=accepted_max=0` 和 `failure_min=failure_max=5`，证明单 rank 碰撞失败不会再
造成 MPI 控制流分叉。`background_tail_collision_weight_test` 还必须输出
`ensemble-statistical=1`、`particles=2048`；checkpoint manifest 必须包含
`tail_collision_weight_algorithm sentoku-kemp-bounded-v1`。否则不得提交 Beam 12 fs 阻断门。

生产 `fp_solver` 不支持 `--overwrite-output`，重跑前先清理目录：

```bash
rm -rf ./output/h9_flux_interface_17f_nobeam_40fs \
  ./output/h9_flux_interface_17f_beam_3fs \
  ./output/h9_flux_interface_17f_beam_12fs \
  ./output/h9_flux_interface_17f_beam_40fs
```

**基线档 A：无 Beam 平衡（仅在 Beam 12 fs 阻断档通过后提交）**：Beam
关闭，验证碰撞热化下 Maxwellian 平衡不自发产生场/噪声，40 fs 长跑。

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 0 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 40 \
  --output-dir ./output/h9_flux_interface_17f_nobeam_40fs
```

**基线档 B：Beam 3 fs（仅在 Beam 12 fs 阻断档通过后提交）**：对照 H8
短跑，验证碰撞开启后短跑账本闭合。

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 3 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3 \
  --output-dir ./output/h9_flux_interface_17f_beam_3fs
```

**优先阻断档：Beam 12 fs（四个生产命令中首先提交）**：覆盖当前碰撞配置
下约 10.18 fs 的首次转换，验证 tail 初始生成、短时跨表示账本和本次人口
爆炸修复。

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 12 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 12 \
  --output-dir ./output/h9_flux_interface_17f_beam_12fs
```

12 fs 作业结束后立即执行阻断验收；任何一项失败都不要提交 40 fs：

```bash
H9_12=./output/h9_flux_interface_17f_beam_12fs

grep -R "^tail_collision_weight_algorithm=sentoku-kemp-bounded-v1$" \
  "$H9_12"/snapshot_*/manifest.dat
grep -R "^tail_collision_max_particle_growth=0$" \
  "$H9_12"/snapshot_*/manifest.dat

awk '
NR==1 {
  for (i=1; i<=NF; ++i) col[$i]=i;
  next;
}
{
  tail=$(col["tail_particle_count"]);
  conv=$(col["conversion_N"]);
  if (have_previous && conv == 0 && tail > previous_tail) {
    ++unexplained_growth_steps;
    print "UNEXPLAINED_TAIL_GROWTH", $1, $2, previous_tail, tail;
  }
  previous_tail=tail;
  have_previous=1;
  if (tail > max_tail) max_tail=tail;
  if ($(col["collision_s"]) > max_collision_s)
    max_collision_s=$(col["collision_s"]);
}
END {
  print "unexplained_growth_steps=" unexplained_growth_steps+0;
  print "max_tail_particle_count=" max_tail+0;
  print "max_collision_s=" max_collision_s+0;
  exit(unexplained_growth_steps != 0);
}' "$H9_12/vpfp_step_diagnostics.dat" \
  | tee "$H9_12/tail_growth_gate.result"
```

验收要求：两个 `grep` 均命中；`unexplained_growth_steps=0`；step 413--444
不再出现 tail 数近翻倍及 collision wall time 指数增长。这里的 awk 门只用于
发现“无 conversion 时粒子数增加”的明确回归；最终仍需结合开放出流和
conversion 累计账本审查全局人口守恒。

**长期档：Beam 40 fs [controller-off 已通过]**：覆盖当前约
10.18 fs 的转换起点以及后续 tail 生成、碰撞和人口增长，快照
3/12/25/40 fs。

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/h9_flux_interface_17f_beam_40fs
```

本档已完成1564步且无拒绝步，正式validator结果为 `status=PASS`。全局/局部tail
宏粒子、Gauss、转换残差、tail数量账和后100步wall-time均通过既定门槛；详细数值见
§0.9和本节第4项实际验收结果。

四档结果目录不得混写。比较内容：combined 热化（bulk 温度与 tail 温度
趋同）、能量交换（`collision_reservoir` 与 K_e/K_b 账本闭合）、tail
population（`tail_particle_count` 演化与 combined 能谱连续性）；并对照
H8 短跑与无碰撞 H4/H5 短跑评估碰撞物理速率（s12 绝对定标与 moment-closure
ν₀ 合理性，§15 阶段 H8 未覆盖风险项）。

**人口控制 A/B [当前禁止执行，历史命令]**：阈值接口谱连续性与 controller-off 资源门通过前，
不得运行本段。保留以下内容仅用于未来独立 A/B，不是 17F 通过后的下一步：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 20 \
  --tail-target-particles-per-bin 64 \
  --tail-max-particles-per-bin 1024 \
  --tail-max-weight-ratio 8 \
  --stop-time-fs 40 --diagnostic-level 1 \
  --diagnostic-interval 1 --snapshot-times 3,12,25,40 \
  --output-dir ./output/hybrid_h9_sk_bounded_v1_beam_40fs_control_on
```

输出到独立目录 `output/hybrid_h9_sk_bounded_v1_beam_40fs_control_on`。验收不是只看宏粒子减少，而是同时比较 combined $N/P/K$、阈值两侧归一化能谱、$J_x/\Pi_{xx}/\Pi_\perp$、电场包络、噪声和 wall time。

**H9 历史结论（2026-08-08，全局源谱阻断回归后；已由2026-08-11结果更新）**：
Beam 12 fs 人口爆炸阻断门、无 Beam 40 fs 平衡和 Beam 40 fs controller-off
的基本碰撞/账本路径已通过；不等权碰撞不再在 `conversion_N=0` 时生成 residual
宏粒子。H9 整体仍为 **partial pass**，原因是 6 MeV bulk--tail 阈值处存在
显著人工谱谷。全局四阶段源谱已经证明该空洞在 converter 提取前就存在，
converter、碰撞和 population controller 均不是其直接来源。population
controller 必须保持关闭；单元内守恒加载候选已经失败，下一步只执行
§7.11.13--§7.11.15 的真实转换事件矩审计与分支决策。只有新表示方案完成独立
验收后才决定是否跑40 fs controller-off。该后续工作现已完成并通过；当前是否进入
120 fs应以§0.9所列的阈值统一能谱验收和资源外推为准，而不是继续重复本历史判断。

#### 17.10.1 阈值接口专项修复的编译与执行顺序

**subcell 可行性审计 [已完成、路线已停止]。** 以下命令从 §7.11.10 迁入，
仅用于复核历史否决结论，不属于当前生产前置门，也不得据此恢复 subcell 路线：

```bash
cmake --build build -j4 --target \
  bulk_tail_subcell_feasibility_test

mkdir -p ./output/h9_subcell_interface

./build/bulk_tail_subcell_feasibility_test \
  --case self-moments \
  --result ./output/h9_subcell_interface/feasibility_self_moments.result

./build/bulk_tail_subcell_feasibility_test \
  --case center-augmented \
  --result ./output/h9_subcell_interface/feasibility_center_augmented.result

./build/bulk_tail_subcell_feasibility_test \
  --case constraint-hierarchy \
  --result ./output/h9_subcell_interface/feasibility_constraint_hierarchy.result
```

三项历史结果已完成分析，并共同否决固定 subcell 支撑路线。除非离散表示设计发生
根本变化，否则不重跑。

以下命令针对 §7.11。第 0--8 步已经完成根因定位，现保留为历史验收记录；当前
执行入口是后面的第 9--13 步。每一步必须先分析 `.result`，不得一次提交到
40 fs。本节后面的旧作业号和实施记录只用于解释历史结果。

**第 0--8 步：全局源谱与近轴中心支撑定位 [已完成]。**

本阶段当时要求先完成 §7.11.7 第 1--4 项：全局 conversion-source 诊断、逐环
支撑审计、`near-axis-narrow` 案例、真实 PASS 判据和截断写文件。上述内容现已
完成，以下命令只作为定位证据保留。

**[2026-08-08 小修复已完成]**：

- `bulk_tail_threshold_interface_mpi_test` 现在始终检查 threshold-aware 的显式
  edge-group 谱；`near-axis-narrow` 和 golden 路额外强制细谱误差
  $\le10^{-10}$。普通光滑谱的细 bin 误差不再静默忽略，而是写出
  `fidelity_warning=1`、`fine_spectrum_required=0` 和独立 edge/fine PASS 标志；
- `VpfpDiagnostics::write_conversion_source_ledger()` 不再把
  `created_tail` 复制成 accepted 数据，而是从完整步最终接受的
  `BackgroundTailPIC` 重新分箱。列名明确改为
  `accepted_tail_total_N/K`，表示包含既有 tail 的最终接受态总谱；
- 本机 MPI 桩下 golden 和 threshold-aware 单 rank 新判据均通过，生产
  `fp_solver` 已完成全量编译/链接检查。真实 1/2/5 rank 集体语义仍以集群结果
  为准。

**第 1 步：编译独立测试和生产程序。**

```bash
cmake --build build -j4 --target \
  hybrid_threshold_support_audit_test \
  bulk_tail_threshold_interface_test \
  bulk_tail_threshold_interface_mpi_test \
  bulk_tail_multibin_test bulk_tail_single_cell_test \
  bulk_tail_poisson_invariance_test bulk_tail_transaction_test \
  fp_solver
```

**第 2 步：建立全新的结果目录。** 禁止在旧 `.result` 后追加新块。

```bash
rm -rf ./output/h9_threshold_interface_v2
mkdir -p ./output/h9_threshold_interface_v2
```

**第 3 步：逐环只读速度支撑审计。**

```bash
./build/hybrid_threshold_support_audit_test \
  --convert-energy-mev 6.0 \
  --uperp-rings 0,1,2,3,4,5,6,7 \
  --bin-widths-mev 0.05,0.1,0.2 \
  --result ./output/h9_threshold_interface_v2/support_by_uperp.result
```

验收时必须直接检查最低环的
`nearest_below_Kout_mev`、`nearest_above_Kout_mev` 和
`crossing_gap_mev`。若 CLI 尚不支持这些参数，说明 §7.11.2 尚未实现，不能用
旧 `support_grid.result` 替代。

**第 4 步：单 rank A/B/C/D 测试。** D 即 `near-axis-narrow`。

```bash
for BW in 0.05 0.1 0.2; do
  ./build/bulk_tail_threshold_interface_test \
    --case all --policy all --bin-width-mev "$BW" \
    --result "./output/h9_threshold_interface_v2/single_bw${BW}.result" || exit 1
done
```

必须分别查看 `edge_spectrum_L1` 与细 bin `spectrum_L1_relative`。前者通过只
证明粗组守恒，不能覆盖后者。

**第 5 步：原转换器守恒/Poisson/事务回归。**

```bash
./build/bulk_tail_single_cell_test \
  --result ./output/h9_threshold_interface_v2/single_cell.result || exit 1
./build/bulk_tail_poisson_invariance_test \
  --result ./output/h9_threshold_interface_v2/poisson_invariance.result || exit 1
./build/bulk_tail_transaction_test \
  --result ./output/h9_threshold_interface_v2/transaction.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_multibin_test \
    --result "./output/h9_threshold_interface_v2/multibin_n${NP}.result" || exit 1
done
```

**第 6 步：MPI 阈值接口双策略验收。** golden 与 threshold-aware 必须分别
运行，不能只跑 candidate。

```bash
for POLICY in golden threshold-aware; do
  for NP in 1 2 5; do
    yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
      ./build/bulk_tail_threshold_interface_mpi_test \
      --case all --policy "$POLICY" --bin-width-mev 0.1 \
      --result "./output/h9_threshold_interface_v2/mpi_${POLICY}_n${NP}.result" \
      || exit 1
  done
done

grep -H '^status=' ./output/h9_threshold_interface_v2/*.result
```

只有 1/2/5 rank 的两种策略都生成各自独立文件，且 `near-axis-narrow` 的六矩、
谱和 rank 一致性达到 §7.11.6，MPI 门才完成。

**第 7 步：12 fs 全局源谱阻断回归 [已完成]。** 当时仅在第 3--6 步通过后执行。

```bash
rm -rf ./output/hybrid_h9_threshold_source_beam_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.4,10.5,10.6,12 \
  --output-dir ./output/hybrid_h9_threshold_source_beam_12fs
```

该运行必须在快照中产生非零且带显式边界的 global conversion source，并能逐
bin 核对：

```text
pre_extraction_bulk -> removed_bulk -> created_tail -> accepted_tail
```

若仍写出 `conversion_source_bins 0`，测试直接失败，不分析 40 fs。

**第 8 步：分支决策 [已完成]。**

- 提取前已空：先建立有限体积单元积分谱，判断空洞是单元中心采样假象还是单元
  体积内确实无支撑；不得直接跳到静态 cut-cell；
- 提取前非空、创建后为空：只修 converter；
- 创建后非空、接受态为空：只修对应 drift/kick、MPI 迁移或沉积阶段。

最新结果属于“提取前中心采样已空”，因此执行下面第 9--13 步。完成对应修复并
重新通过第 3--7 步后，才允许 40 fs controller-off 回归。

**第 9 步：实现只读有限体积单元谱审计。**

实现 `src/tail_subcell_quadrature.h/.cpp` 和
`tests/bulk_tail_cell_volume_spectrum_test.cpp`。积分必须使用真实
$u_\parallel/u_\perp$ 面、圆柱测度 $2\pi u_\perp\,du_\perp$ 和确定性求和；
第一版按 cell 常值积分，不得擅自引入与生产 remap 不一致的高阶重构。实现完成后：

```bash
cmake --build build -j4 --target \
  bulk_tail_cell_volume_spectrum_test fp_solver

rm -rf ./output/h9_subcell_interface
mkdir -p ./output/h9_subcell_interface

./build/bulk_tail_cell_volume_spectrum_test \
  --convert-energy-mev 6.0 \
  --uperp-rings 0,1,2,3,4,5,6,7 \
  --bin-widths-mev 0.05,0.1,0.2 \
  --result ./output/h9_subcell_interface/cell_volume_spectrum.result
```

验收要求：总数积分相对误差不大于 $10^{-12}$；输出每个跨阈值 cell 的几何体积、
阈值上下体积分数、$N/P_x/K$；明确报告 6.0--6.2 MeV 在 cell-volume 参考中是否
非零。只允许读取状态，不得修改 bulk 或创建 tail。

**[2026-08-08 第一次测试假失败与测试修正]**：第一次结果虽然写出
`status=FAIL`，但三档分箱均已有
`threshold_window_6p0_6p2_volume_N=0.005204062676...`，失败仅来自旧判据强制
`straddling_cell_count>0`。当前6 MeV阈值位于相邻速度单元的公共面附近，因此
没有正体积单元同时跨越阈值是合法情况。测试已改为直接统计生产
`is_conversion_cell` 掩码内的子单元支撑，并补充窗口 $N/P_x/K$、实际贡献cell
数以及cell-center到cell-volume的能量表示差，不再把能量残差硬编码为0。

本地MPI桩回归得到：`selected_cell_count=112`、
`selected_threshold_support_cell_count=16`、窗口
$N=0.005204062676...$、$P_x\simeq-1.10\times10^{-39}$，三档分箱均
`status=PASS`。该结果说明第10步在数学上存在候选支撑；但集群必须用新二进制
重跑本节命令后，才把第9步正式标记为通过。

**第 10 步：仅在 cell-volume 参考显示阈值内有支撑时，实现单元内守恒 PIC
加载。**

实现 `tests/bulk_tail_subcell_loading_test.cpp`，并让 converter 仅对**现有中心
判据已经选中且本步将完整清零的 cell** 使用确定性子单元积分点。不得从中心仍在
bulk 侧的阈值相交 cell 每步抽取固定比例，否则会重复耗尽同一 cell。所有 PIC
support 必须满足 $K\ge K_{\rm out}$，采用
方位角四元组消除伪横向动量，并通过现有 `TailMomentConstraint` 求非负权重以匹配
$N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$。若非负问题不可行，必须显式返回失败；禁止负
权、全局能量补丁或把 support 移到错误能量 bin。

```bash
cmake --build build -j4 --target \
  bulk_tail_subcell_loading_test \
  bulk_tail_single_cell_test bulk_tail_poisson_invariance_test \
  bulk_tail_transaction_test bulk_tail_multibin_test \
  bulk_tail_threshold_interface_test \
  bulk_tail_threshold_interface_mpi_test fp_solver

for BW in 0.05 0.1 0.2; do
  ./build/bulk_tail_subcell_loading_test \
    --case near-axis-narrow --bin-width-mev "$BW" \
    --result "./output/h9_subcell_interface/subcell_bw${BW}.result" || exit 1
done
```

验收要求：六矩相对误差均不大于 $10^{-10}$；权重非负；阈值上下无非法 support；
相对 cell-volume 参考的细谱 $L_1$ 误差不大于 5%。几何节点和约束矩阵必须预计算，
生产转换不得逐 cell 重复分配大数组。

**第 11 步：重跑转换器基础门和 MPI 门 [暂停]。**

第10步已经失败，本步命令保留但当前禁止执行。只有 §7.11.15 选择并实现新的统一
矩定义或局部网格候选后，才重新启用。

```bash
./build/bulk_tail_single_cell_test \
  --result ./output/h9_subcell_interface/single_cell.result || exit 1
./build/bulk_tail_poisson_invariance_test \
  --result ./output/h9_subcell_interface/poisson_invariance.result || exit 1
./build/bulk_tail_transaction_test \
  --result ./output/h9_subcell_interface/transaction.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_multibin_test \
    --result "./output/h9_subcell_interface/multibin_n${NP}.result" || exit 1
done

for POLICY in golden threshold-aware; do
  for NP in 1 2 5; do
    yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
      ./build/bulk_tail_threshold_interface_mpi_test \
      --case all --policy "$POLICY" --bin-width-mev 0.1 \
      --result "./output/h9_subcell_interface/mpi_${POLICY}_n${NP}.result" \
      || exit 1
  done
done
```

**第 12 步：重跑 12 fs subcell阻断回归 [取消]。** 原候选没有通过第10步，
以下命令不得执行。

```bash
rm -rf ./output/hybrid_h9_subcell_beam_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.4,10.5,10.6,12 \
  --output-dir ./output/hybrid_h9_subcell_beam_12fs
```

12 fs 门要求：无固定硬空 bin；全局四阶段数目链继续闭合；转换 $N/P_x/K$ 误差
保持原量级；Gauss 不退化；tail 宏粒子数不超过当前基线的 2 倍；转换阶段 wall
time 不超过当前基线的 1.5 倍；不得出现 NaN、负权、重复 ID 或 rank 不一致。

**第 13 步：失败分支 [已触发]。** 第10步已证明当前中心六矩目标与展开subcell
支撑不兼容。禁止通过无限扩大 $u_{\max}/N_v$、开启population controller、
反复静态切割同一cell、提高非负拟合迭代数或全局能量修正强行通过。

**第 14 步：实现真实转换事件矩审计 [已完成]。** 单rank全部案例和真实MPI
1/2/5 rank均为 `status=PASS`。以下命令保留为回归入口。

```bash
cmake --build build -j4 --target \
  bulk_tail_real_moment_audit_test \
  bulk_tail_real_moment_audit_mpi_test \
  fp_solver

rm -rf ./output/h9_real_moment_audit_tests
mkdir -p ./output/h9_real_moment_audit_tests

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_real_moment_audit_test \
  --case all \
  --result ./output/h9_real_moment_audit_tests/single_rank.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_real_moment_audit_mpi_test \
    --case all \
    --result "./output/h9_real_moment_audit_tests/mpi_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_real_moment_audit_tests/*.result
```

只有四类单rank案例和1/2/5 rank归约全部通过，才能运行真实12 fs审计。

**第 15 步：真实12 fs接受事件审计 [已完成，判定为红色]。** 不使用subcell
加载；只增加只读统计。

从0 fs运行的标准命令：

```bash
rm -rf ./output/hybrid_h9_real_moment_audit_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-cell-moment-audit \
  --tail-cell-moment-audit-top-cells 64 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.4,10.5,10.6,12 \
  --checkpoint-times 10.4 \
  --output-dir ./output/hybrid_h9_real_moment_audit_12fs
```

上面是今后复跑使用的修正版命令，增加了 `--checkpoint-times 10.4`。已经完成的第15步
历史运行没有该参数，因此其10.4 fs目录只是snapshot，**不能供 `--restart-dir` 使用**。
时间checkpoint由程序写入 `output-dir` 下的 `checkpoint_target...` 完整目录；使用前必须
检查manifest及全部80个rank文件，并确认完整rank状态、RNG、Beam、tail、bulk和累计
账本均已保存。snapshot不能伪装成checkpoint；若不存在该完整目录，第16A必须从0 fs
重跑。

若已有**首次转换前**且配置一致的10.4 fs checkpoint，可节省前段运行：

```bash
CHECKPOINT_104=/absolute/path/to/checkpoint_target10.4fs_...

rm -rf ./output/hybrid_h9_real_moment_audit_10p4_to_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_104" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-cell-moment-audit \
  --tail-cell-moment-audit-top-cells 64 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.5,10.6,12 \
  --output-dir ./output/hybrid_h9_real_moment_audit_10p4_to_12fs
```

checkpoint必须使用80 rank、相同速度网格/边界/碰撞/tail配置，且累计
`conversion_N=0`。10.5 fs或更晚的已转换checkpoint不能用于替代提取前审计。

验收文件必须存在且非空：

```bash
AUDIT_DIR=./output/hybrid_h9_real_moment_audit_12fs
test -s "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" || exit 1
grep -R -E 'tail_cell_moment_audit|tail_subcell_loading' "$AUDIT_DIR"/manifest*.dat
tail -n 20 "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat"
```

manifest实际确认 `tail_cell_moment_audit=1`、`tail_subcell_loading=0`。54个真实事件
给出 $R_{L1}^{\Pi_\perp}=3.26\times10^{-2}$、单cell最大约11.3%，且volume目标
可行性为0/743，触发红色分支。另发现416行无转换全零占位和相对误差分母使用
有量纲常数1；两者属于诊断显示缺陷，不改变红色结论。

**第 16 步：执行红色分支。** 必须按16A、16B顺序执行；16A不得顺带修改物理，
16B不得先接入生产。

**第 16A 步：修复审计显示并生成真实速度cell直方图。** 实现§7.11.16B后构建：

```bash
cmake --build build -j4 --target \
  bulk_tail_real_moment_audit_test \
  bulk_tail_real_moment_audit_mpi_test \
  fp_solver

rm -rf ./output/h9_real_moment_audit_tests_v2
mkdir -p ./output/h9_real_moment_audit_tests_v2

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_real_moment_audit_test --case all \
  --result ./output/h9_real_moment_audit_tests_v2/single.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_real_moment_audit_mpi_test --case all \
    --result "./output/h9_real_moment_audit_tests_v2/mpi_n${NP}.result" \
    || exit 1
done
grep -H '^status=' ./output/h9_real_moment_audit_tests_v2/*.result
```

随后从首次转换前、`conversion_N=0`、80 rank且配置一致的10.4 fs**完整checkpoint**
重跑到12 fs；第15步命令生成的10.4 fs snapshot不满足此条件。若没有另行保存的完整
checkpoint，必须使用第15步的0--12 fs命令并改用新的输出目录。checkpoint复跑命令：

```bash
CHECKPOINT_104=/absolute/path/to/checkpoint_target10.4fs_...

rm -rf ./output/hybrid_h9_real_moment_audit_v2_10p4_to_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_104" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-cell-moment-audit \
  --tail-cell-moment-audit-top-cells 64 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.5,10.6,12 \
  --output-dir ./output/hybrid_h9_real_moment_audit_v2_10p4_to_12fs
```

验收：

```bash
AUDIT_DIR=./output/hybrid_h9_real_moment_audit_v2_10p4_to_12fs
test -s "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" || exit 1
test -s "$AUDIT_DIR/bulk_tail_moment_audit_velocity_histogram.dat" || exit 1

# 主表不得含正请求数为0的占位行；第5列是positive_request_cell_count。
awk 'NR>1 && $5==0 {bad++} END {print "zero_event_rows=" bad+0; exit(bad!=0)}' \
  "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" || exit 1

grep -q 'center_l1_5' "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" || exit 1
grep -q 'delta_l1_5' "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" || exit 1
find "$AUDIT_DIR" -name 'manifest_rank0.dat' -exec \
  grep -H -E 'tail_cell_moment_audit=1|tail_subcell_loading=0' {} \;

python tools/analyze_bulk_tail_moment_audit.py \
  --audit "$AUDIT_DIR/bulk_tail_moment_audit_accepted_steps.dat" \
  --velocity-histogram "$AUDIT_DIR/bulk_tail_moment_audit_velocity_histogram.dat" \
  --result "$AUDIT_DIR/audit_summary.result" || exit 1

cat "$AUDIT_DIR/audit_summary.result"
```

16A通过标准：四个独立测试结果均 `status=PASS`；`zero_event_rows=0`；主表和速度
直方图均非空；分析脚本给出 `input_valid=1`；在相同checkpoint、终止时间、时间步和
物理配置下，真实event/cell计数必须分别严格为54和743，不允许用浮点终止误差解释
计数变化；若从0 fs以不同步序重跑，则必须先比对首次转换时间和最后接受步，再单独
说明差异，不能直接判通过。新算得
$R_{L1}^{\Pi_\perp}$ 应约为旧数据重构的 $3.26\times10^{-2}$，不能因修改诊断而
突然降到绿色。若数值发生数量级变化，停止16B并审计归约/列定义。

**第 16B 步：固定 $(192,64)$ 网格离线原型。** 只有16A全部通过后实现
§7.11.16C；新增target并构建：

```bash
cmake --build build -j4 --target \
  tail_interface_grid_design_test \
  tail_interface_grid_replay_test \
  tail_interface_grid_replay_mpi_test

rm -rf ./output/h9_tail_interface_grid_design
mkdir -p ./output/h9_tail_interface_grid_design

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_interface_grid_design_test --case all \
  --result ./output/h9_tail_interface_grid_design/design.result || exit 1

HIST=./output/hybrid_h9_real_moment_audit_v2_10p4_to_12fs/bulk_tail_moment_audit_velocity_histogram.dat

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_interface_grid_replay_test \
  --input "$HIST" --profiles G0,Gx,Gp,G2 \
  --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0 \
  --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 \
  --result ./output/h9_tail_interface_grid_design/replay_n1.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/tail_interface_grid_replay_mpi_test \
    --input "$HIST" --profiles G0,Gx,Gp,G2 \
    --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0 \
    --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 \
    --result "./output/h9_tail_interface_grid_design/replay_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_tail_interface_grid_design/*.result

# G0必须先复现16A的当前网格基线；候选必须显式写出稀疏支撑、人口和步长预算。
grep -H -E '^(grid_name|candidate_status|R_L1_|volume_self_sparse_failed_count|max_sparse_support_count|estimated_created_macroparticles|estimated_particle_ratio_to_center_quartet|estimated_velocity_dt_ratio|estimated_operator_work_ratio)=' \
  ./output/h9_tail_interface_grid_design/*.result
```

此处文件末尾的 `status=PASS` 只代表测试程序、守恒交叠和MPI一致性正确执行；每个
候选必须另写 `candidate_status=GREEN|GRAY|RED|INVALID`。选择规则固定为：先排除
`INVALID/RED`；在GREEN中优先 `Gp`、再 `Gx`、最后 `G2`；同类候选选择最大
`min_du`、最小相邻宽度比者，不选择误差最小但网格最僵硬者。若没有GREEN，停止
局部网格路线并进入通量式转换设计，不得扩大参数扫描或增加速度cell。

在评价Gx/Gp/G2前，G0必须复现16A原始直方图的event/cell计数和六矩指标；特别是
$R_{L1}^{\Pi_\perp}$ 应约为 $3.26\times10^{-2}$。G0不一致说明回放输入、圆柱体积
权重或矩定义错误，此时所有候选结果均为 `INVALID`。任何GREEN候选还必须同时满足
`volume_self_sparse_failed_count=0`、`max_sparse_support_count<=7`、
`estimated_particle_ratio_to_center_quartet<=2.0`；缺少这些列不能判GREEN。
同时要求 `estimated_velocity_dt_ratio>=0.8`、`estimated_operator_work_ratio<=1.5`；否则
即使六矩通过，也因生产性能不可接受判为RED。

**[已完成-16B实施记录]** 新增独立模块并接入三个测试target（`CMakeLists.txt`）：

```text
src/tail_interface_grid_design.h
src/tail_interface_grid_design.cpp
tests/tail_interface_grid_design_test.cpp
tests/tail_interface_grid_replay_test.cpp
tests/tail_interface_grid_replay_mpi_test.cpp
tests/tail_interface_replay_common.h
```

生产模块不修改 `CylindricalVelocityGrid::init()` 的生产默认，不引用
`main_vpfp.cpp` 全局CLI状态，也不改动 `Param::Nv/Param::Nmu`。候选网格通过
`make_candidate_grid()` 复用公开的 `build_cell_geometry_and_moments` /
`build_moment_closure_table` 构造（先按 `init_grid` 的尺寸预分配公开字段），
六矩复用 `mass_cell_moments`，子单元节点复用 `TailSubcellQuadrature::nodes`，
稀疏支撑约简复用 `tail_compress_moment_supports`。接口、候选构造与约束、圆柱
交叠重映射、volume-self可行性、稀疏约简、步长/工作量预算和离线门均在
`tail_interface_grid_design.cpp` 单点实现；测试只调用生产函数，不复写monitor
积分、face反演或圆柱交叠公式。

关键实现约定（写入代码注释，集群结果必须按此解释）：

- monitor equidistribution 在基线网格坐标 $\xi$（基线face位于均匀
  $\xi=j/n_{cells}$）上做高分辨率累计积分并按等monitor质量反演，
  $w=1$ 时逐位复现基线faces（`Gx` 只移动阈值附近的faces，不会把 sinh 拉伸
  区变成均匀网格）；每基线cell细分512个子区间。$u_*$ 为直方图中低
  $u_\perp$（$\le 4\times$ 首个 $u_\perp$ cell宽度）转换cell的质量加权
  $|u_\parallel|$ 质心，不写死为某个iv。
- 圆柱交叠 $\theta$ 严格按文档公式（$u_\perp$ 用 $u^2$ 差加权），每个旧cell
  校验 $|\sum\theta-1|\le10^{-14}$（`max_partition_error`，纳入 `status=PASS`）。
- volume-self 稀疏约简目标为六矩；尺度化残差
  $|\mathrm{res}_X|/\max(|l1\_target_X|,\ scale\_floor_X)$，其中
  `scale_floor_X` 由该cell数目乘速度/能量特征尺度构造并显式写入结果
  （`scale_floor_max_*`）；数目残差 $\le10^{-12}$、其余五矩 $\le10^{-10}$。
- `center_target_feasible_count/failed_count` 仅信息项，不作硬门；
  `volume_self_failed_count=0` 与 `volume_self_sparse_failed_count=0` 为硬门。
- 步长预算使用生产算子约束：u_parallel PPM 受力输运（
  `semi_lagrangian_cfl=2.5` × 每内face的4-cell重构stencil半跨度 / $a_u$，
  名义 $E_{\rm rep}=10^{12}\ \mathrm{V/m}$，比值与场幅无关）与碰撞显式
  交叉扩散（`velocity_space_cfl=0.35$ × $\min(du_\parallel^2,du_\perp^2)/D$，
  系数用 moment-closure 公式按同一状态采样）。
- 工作量预算 = $0.7\times$固定cell扫描 + $0.15\times$预计PIC创建 +
  $0.15\times$碰撞pair×substep；碰撞pair数随粒子数线性（每粒子每substep配
  一对），权重0.7/0.15/0.15对应生产步wall-time中 bulk Vlasov 扫描主导的结构
  （诊断 `vlasov_s`~0.1 s vs `beam_s`/`collision_s`~$10^{-4}$--$10^{-3}$ s）。
  各分项均写入result（`scan_cost_ratio`、`pic_creation_ratio`、
  `collision_pair_ratio`、`collision_substep_ratio`）。
- 候选分类：`INVALID`（网格违反宽度约束或G0 identity失败）；`GREEN`（全部
  离线门通过）；`GRAY`（硬门通过但六矩/重要cell物理门落在灰带
  $10^{-3}<R_{L1}\le10^{-2}$ 或 $10^{-2}<$重要cell相对差$\le10^{-1}$）；
  其余 `RED`。`status=PASS` 只代表程序、守恒交叠与MPI一致性正确。

本地验证（g++ + MPI stub，单rank；与16A相同方式）：

```text
tail_interface_grid_design_test  --case all              -> status=PASS（6控制用例全过）
tail_interface_grid_replay_test  （真实直方图，100组参数）-> status=PASS
tail_interface_grid_replay_mpi_test （n=1，与串行结果逐块一致）-> status=PASS
```

控制用例（g0-identity、single-cell-overlap、constant-cells、symmetric-pair、
Maxwellian/漂移Maxwellian往返）全部PASS；两次replay运行结果逐bit一致；
串行与MPI(n=1)候选块逐块一致。G0在真实直方图上精确复现16A：
`histogram_events=54`、`histogram_cell_requests=611`、`histogram_bins=4`，
$R_{L1}^{\Pi_\perp}=3.3335474067264968\times10^{-2}$（16A为
$3.3335474067264885\times10^{-2}$）、`max_cell_rel_Piperp=0.11275207511488894`，
其余分量一致到机器精度。

集群执行（16B验收入口）：

```bash
cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release \
         -DBUILD_TESTING=ON
make -j4 tail_interface_grid_design_test \
         tail_interface_grid_replay_test \
         tail_interface_grid_replay_mpi_test

rm -rf ./output/h9_tail_interface_grid_design
mkdir -p ./output/h9_tail_interface_grid_design

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_interface_grid_design_test --case all \
  --result ./output/h9_tail_interface_grid_design/design.result || exit 1

HIST=./output/hybrid_h9_real_moment_audit_12fs/bulk_tail_moment_audit_velocity_histogram.dat

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_interface_grid_replay_test \
  --input "$HIST" --profiles G0,Gx,Gp,G2 \
  --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0 \
  --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 \
  --result ./output/h9_tail_interface_grid_design/replay_n1.result || exit 1

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/tail_interface_grid_replay_mpi_test \
    --input "$HIST" --profiles G0,Gx,Gp,G2 \
    --ax-values 0.5,1.0,2.0 --aperp-values 0.5,1.0,2.0 \
    --sigma-x-cells 1,2,4 --sigma-perp-cells 1,2,4 \
    --result "./output/h9_tail_interface_grid_design/replay_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_tail_interface_grid_design/*.result
grep -H -E '^(grid_name|candidate_status|R_L1_|volume_self_sparse_failed_count|max_sparse_support_count|estimated_created_macroparticles|estimated_particle_ratio_to_center_quartet|estimated_velocity_dt_ratio|estimated_operator_work_ratio)=' \
  ./output/h9_tail_interface_grid_design/*.result
```

注：文档正文此处 `HIST` 原指向 `hybrid_h9_real_moment_audit_v2_10p4_to_12fs`
（16A阶段目录名），该目录未实际创建；16A实际产物在
`hybrid_h9_real_moment_audit_12fs`，集群验收统一使用后者路径。本地预跑结论：
第一轮扫描（$A_x,A_\perp\in\{0.5,1,2\}$，$\sigma_x,\sigma_\perp\in\{1,2,4\}$）
共100个候选（G0×1、Gx×9、Gp×9、G2×81），无任何 `GREEN`：G0/Gx 的
$R_{L1}^{\Pi_\perp}$ 保持 $3.3\times10^{-2}$（Gx 不修 $\Pi_\perp$，按文档自动判
失败），最激进且仍满足宽度约束的 Gp/G2（$A_\perp=1$，$\sigma_\perp=4$）仅把
$\Pi_\perp$ 降到 $1.4\times10^{-2}$，仍超灰带上限；`Gx_ax1_sx1` 等因阈值下
分数 $9.5\times10^{-2}$ 超 $10^{-6}$ 门判RED，`ap2`/`ax2` 强细化组合因最小
宽度约束判INVALID。按选择规则，无GREEN即停止局部网格路线并进入通量式转换
设计（§7.11.16D 第7项），不得扩大参数扫描或增加速度cell；最终以集群
`grep -H '^status='` 与候选 `candidate_status=` 结果复核后为准。

**[已完成-16B集群验收与独立数据核算]** 集群执行上述命令后的产物：

```text
output/h9_tail_interface_grid_design/
  design.result        status=PASS
  replay_n1.result     status=PASS
  replay_n2.result     status=PASS
  replay_n5.result     status=PASS
```

四个结果文件 `status=PASS`（只代表测试程序、守恒交叠与MPI一致性正确）。在此
基础上按文档验收标准做了不依赖程序自报的独立数据核算，全部通过：

1. 直方图输入独立统计（Python 直接读
   `hybrid_h9_real_moment_audit_12fs/bulk_tail_moment_audit_velocity_histogram.dat`）：
   行数135（=16A `velocity_bin_rows`）、事件54（=16A `real_events`）、
   唯一$(iv,imu)$ bin 4 个（iv=185, imu=1..4）、request cell 总数611
   （=16A `real_cells`）、总质量 $4.24175828477356\times10^{10}$。
2. G0 六矩独立重算（Python 按文档公式重建192×64 sinh 网格 + 16节点
   Gauss-Legendre 体积矩）与集群 G0 块一致到机器精度：
   $R_{L1}^{J_x}=3.0016464\times10^{-6}$、$R_{L1}^{K}=1.0963442\times10^{-6}$、
   $R_{L1}^{\Pi_{xx}}=1.0017053\times10^{-6}$、
   $R_{L1}^{\Pi_\perp}=0.033335474067264961$（集群 0.033335474067264968，
   差 $2\times10^{-16}$）、`max_cell_rel_Piperp=0.11275207511488877`。
3. 候选关键数值独立重建（重建 monitor 面 + θ 重映射，不调用测试程序）：
   `Gx_ax1_sx1` 的 `below_threshold_number_fraction` 独立重算
   $=0.094565993494598252$ 与集群输出逐位一致；`Gp_ap1_sp4` 的
   $R_{L1}^{\Pi_\perp}=0.014324981765568267$ 与集群
   $0.014324981765568161$ 差 $7\times10^{-15}$（Jx/K/Pixx 差
   $\le1.3\times10^{-10}$）。u* 独立复核 $=13.1504565221269$、
   $j^*=185$、$du_{\rm threshold}=0.849021429707413$。
4. 1/2/5 rank 一致性：三文件候选状态分布完全相同
   （`GREEN=0 GRAY=0 RED=49 INVALID=51`）；全部物理量跨 rank 最大相对差
   $\le6\times10^{-10}$（MPI归并求和顺序的浮点末位差异；N/Px 的
   ~$10^{-16}$ 绝对噪声不构成差异）。同一运行内 faces 逐位一致、
   remap 相对差 $\le10^{-10}$ 由测试内部校验。
5. 模型自洽核算：`estimated_operator_work_ratio=1.225`，即
   $0.7\times1+0.15\times1.75+0.15\times1.75\times1$；
   Gp 不动 u∥ → `estimated_velocity_dt_ratio=1`、`min_dupar` 与 G0 相同
   （0.00524636）；Gx 的 $R_{L1}^{\Pi_\perp}=3.308\times10^{-2}\approx$ G0 的
   $3.334\times10^{-2}$，确认 u∥ 细化不能修复 u⊥ 主导误差；所有候选
   `volume_self_sparse_failed_count=0`、`max_sparse_support_count=7`、
   粒子比1.75。

候选失败分解（n1/n2/n5 一致）：`INVALID=51`（`u_perp cell 0` 最小宽度违反
21个=全部 $A_\perp=2$ 的 Gp/G2；`u_parallel cell 6/5` 最小宽度违反30个=
$A_x=2$ 的 Gx/G2）；`RED=49`（42个 `physics-moment gate violation beyond
gray band`：G0/Gx 全部及 Gp/G2 的 $A_\perp\le1$ 组合，Π⊥ 最低
$1.4\times10^{-2}$ 仍超灰带上限；7个 `hard offline gate failed`：σx=1 的
Gx/G2，阈值下分数 $9.46\times10^{-2}$ 超 $10^{-6}$ 门）。

注意：文档命令中 `for NP in 1 2 5` 会用 MPI n=1 输出覆盖串行 `replay_test`
写出的 `replay_n1.result`（集群产物含 `ranks=1` 字段即此原因）；两者候选块
逐块一致，不影响验收，但若需保留串行产物可把 MPI 循环改为 `NP in 2 5`。

**16B 结论（数据核算层面）**：G0 精确复现 16A 基线；候选生成、monitor 反演、
θ 重映射、六矩与预算指标均可独立复现；唯一能把
$R_{L1}^{\Pi_\perp}$ 修进 $10^{-2}$ 的组合（$A_\perp=2$）全部被最小宽度约束
判 INVALID，所有有效候选 $\Pi_\perp\ge1.4\times10^{-2}$。**无 GREEN 是数据
结论而非程序假象**。按 §7.11.16D 选择规则：停止局部网格路线，进入通量式
转换设计（第7项），不得扩大参数扫描或增加速度cell。

以下旧命令和记录保留为历史，不再作为当前执行入口。

**以下为历史执行记录，不构成当前验收步骤。**

**历史第 2 步：旧版只读速度支撑审计。**

```bash
mkdir -p ./output/h9_threshold_interface

./build/hybrid_threshold_support_audit_test \
  --convert-energy-mev 6.0 \
  --result ./output/h9_threshold_interface/support_grid.result
```

**历史第 3 步：单 rank 黄金/当前/候选三路测试。**

```bash
for BW in 0.05 0.1 0.2; do
  ./build/bulk_tail_threshold_interface_test \
    --case all --policy all --bin-width-mev "$BW" \
    --result "./output/h9_threshold_interface/single_bw${BW}.result" || exit 1
done
```

先比较 `golden` 和 `current`。若二者都继承相同谱谷，停止实现
threshold-aware compression，按 §7.11.4 分支 A 处理；只有黄金连续而 current
失真时，才执行分支 B。

**历史第 4 步：原转换器回归。** 候选实现不得破坏已经通过的守恒、Poisson 和
事务语义。

```bash
./build/bulk_tail_single_cell_test \
  --result ./output/h9_threshold_interface/single_cell.result
./build/bulk_tail_poisson_invariance_test \
  --result ./output/h9_threshold_interface/poisson_invariance.result
./build/bulk_tail_transaction_test \
  --result ./output/h9_threshold_interface/transaction.result
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_multibin_test \
  --result ./output/h9_threshold_interface/multibin_n1.result
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_multibin_test \
  --result ./output/h9_threshold_interface/multibin_n2.result
yhrun -N 1 -n 5 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_multibin_test \
  --result ./output/h9_threshold_interface/multibin_n5.result
```

**历史第 5 步：MPI 阈值接口验收。**

```bash
for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_threshold_interface_mpi_test \
    --case all --policy threshold-aware --bin-width-mev 0.1 \
    --result "./output/h9_threshold_interface/mpi_n${NP}.result" || exit 1
done

grep -H '^status=' ./output/h9_threshold_interface/*.result
```

只有所有相关 `.result` 为 `status=PASS` 且满足 §7.11.6 的数值门，才可接入
生产默认路径。

**第 1-2 步实施记录（2026-08-06，本机 MPI 桩 + g++ -O2）**：

- 新增 `tests/hybrid_threshold_support_audit_test.cpp`（§7.11.2）：直接构造
  生产 `CylindricalVelocityGrid`/`HybridVelocityPartition`，输出
  K_out±1 MeV 窗口在 0.05/0.1/0.2 MeV 诊断 bin 下的 bulk/conversion
  中心计数、能量范围、最大未覆盖空隙、相体积和，以及全部 conversion
  cell-center 能量的有序唯一值及相邻间距；不推进方程、不创建粒子。
- 新增 `tests/bulk_tail_threshold_interface_test.cpp` 与
  `tests/bulk_tail_threshold_interface_mpi_test.cpp`（§7.11.3）：
  `BulkTailConverter` 增加 `BulkTailLoadingPolicy`
  （golden/current/threshold-aware），经 `set_loading_policy()` 选择
  （per-converter 状态，非全局/宏）；当时设计要求生产默认保持 current，但后续
  实现已提前切换为 threshold-aware，当前按 §0.6 重新归类为未验收 candidate；
  `HybridVelocityPartition` 增加只读 `conversion_energy_edges`
  （K_out/+0.2/+0.4/+0.8 MeV 显式边界，更高能区对数增宽）与
  `energy_bin_threshold_aware()`（upper_bound），供 C 路使用；策略只改变
  候选支撑压缩，不改变质量提取、时间层、位置、ID、沉积或事务语义。
  三个解析光滑非负分布（exp 单调谱/跨阈值宽高斯/带漂移各向异性谱）分别
  走 A/B/C，报告六矩残差、rho_l2/linf、谱 L1/Linf、阈值细 bin 相对误差、
  对数曲率、谱谷比值、粒子数与 fallback。MPI 测试以 1/2/5 rank 集体归约
  全局六矩、逐 bin 谱与粒子数，失败经 LOR 归约使全部 rank 非零退出；
  无随机数。
- 本机结果：审计测试 `status=PASS`，揭示 6.0–6.2 MeV 仅 22 个 conversion
  支撑（最大空隙 0.053 MeV）而 6.2–6.25 MeV 有 82 个密集支撑；A/B/C 测试
  golden 谱 L1≈5e-15（精确），current 谱 L1≈0.4–1.0 且 anisotropic-drift
  阈值谷比值 1.3e-4（中间 bin 比两侧低 4 个数量级，复现 §7.11.1 谱谷），
  threshold-aware 谱 L1≈0.5（组内 7-support 压缩仍失真，指向 §7.11.4
  分支 B + 组内细 bin 结构约束）；MPI n=1 golden 谱 L1≈3e-14。
  全部 6 个 `.result` 为 `status=PASS`；`bulk_tail_single_cell_test`、
  `bulk_tail_multibin_test`、`bulk_tail_poisson_invariance_test`、
  `bulk_tail_transaction_test` 无回归。集群 1/2/5 rank 与多 bin 宽度结果
  待按本手册执行。

**历史第 3 步实施记录（2026-08-06，§7.11.4 分支 B）**：当时依据不完整的
集群 MPI/edge-group 结果提前实施分支 B；最新复核表明该记录不能等价于完整的
1/2/5 rank 细谱和近轴验收：

- `HybridVelocityPartition::energy_bin()` 改为对显式
  `conversion_energy_edges` 执行 upper_bound（K_out/+0.2/+0.4/+0.8 MeV +
  更高能区对数增宽），生产分组不再跨越阈值审计 bin；旧均匀分组保留为
  `energy_bin_uniform()`，仅作为 §7.11.3 A/B/C 测试的 CURRENT 参考策略；
- `BulkTailConverter` 默认加载策略切换为
  `THRESHOLD_AWARE_COMPRESSION`；该切换当前只能视为 candidate 接入，CURRENT 策略走
  旧均匀分组，golden 不压缩；
- `tail_population_controller` 无需改动（已复用 `partition.energy_bin()`，
  随方法自动同步）；
- `config_hash` 纳入能量边界（restart 校验自动拒绝旧分组配置）；
  checkpoint/snapshot manifest 写出边界数量、边界值与
  `conversion_energy_edges_hash`；
- 接口测试新增 edge 分辨率谱与 §7.11.6 第 2 项门：C 路在每个有非零参考
  质量的阈值细 bin 内 N/K 相对 A 误差 ≤1e-10，edge 分辨率全谱 L1 ≤1e-10。
- 本机验证：C 路 edge 谱 L1≈3.7e-15–5.2e-15、edge bin N/K max_rel≈
  1.4e-14–6.5e-14（全部 ≤1e-10）；B 路（旧均匀）edge bin N/K max_rel≈
  0.38–0.62——量化确认旧分组把阈值 bin 质量搬过边界是谱谷根源；
  `bulk_tail_single_cell_test`/`multibin`/`poisson_invariance`/`transaction`
  /审计的基础回归 PASS，`checkpoint_roundtrip_test` PASS；但 MPI threshold
  测试的旧 PASS 判据没有约束 threshold-aware 细谱，不能据此关闭 §7.11。
- 注意：旧分组配置写入的 checkpoint 与本二进制不兼容（restart 会拒绝），
  属 §7.11.4 要求的"不允许静默续跑"；0.1 MeV 细分辨率下 C 路组内仍有
  谱失真（L1≈0.5），对应 §7.11.4"组内 7-support 结构约束"或生产
  §7.11.6 第 7 项 0.1/0.2/0.4 MeV 收敛验收，需在第 6 步 Beam 12 fs/40 fs
  中评估。

**第 4 步集群回归分析（2026-08-07，作业 6865979）**：四个回归测试
`status=FAIL`，但物理指标全部在门内（single_cell 残差 0/1e-16、
poisson density_l2 3.7e-13/守恒 6.7e-15、transaction 全 flag=1、
multibin max_residual_rel 5e-15）。两因：
1. single_cell/poisson/transaction 的 FAIL 是 `write_result_file` 在
   `./output/h9_threshold_interface/` 不存在时返回 false（完美指标仍
   FAIL，本机已复现该机制）；需先 `mkdir -p ./output/h9_threshold_interface`
   再重跑，与运行手册第 2 步一致；
2. multibin 是真回归：`energy_bin()` 改为阈值边界后忽略 `energy_bins`，
   旧"coarse 单能量组 ≤12×2×7×4=672 粒子"假设失效（实测 coarse 2688）。
   已修复：粒子数上界按 `conversion_energy_edges` 组数计算
   （12×2×N_edge×7×4，本网格 ~6720），本机 PASS；矩 bin 无关性
   （moment_bin_convergence_max≈2.3e-12）保持。
3. 二次集群重跑（作业 6866299）：multibin `all_ok=1`（上界修复生效），
   但所有测试仍 `status=FAIL`，且 multibin n2/n5 出现 rank 分歧（rank 0
   FAIL、其余 PASS）——根因仍是集群 `./output/h9_threshold_interface/`
   目录缺失，rank-0 的 `write_result_file` 在集体归约后失败（其他 rank
   保持 PASS）。已在 multibin 测试中把落盘结果纳入集体归约，任一 rank
   （含 I/O 失败）全体非零退出；集群重跑前先
   `mkdir -p ./output/h9_threshold_interface`。
4. 第 4 步闭环（2026-08-07）：`mkdir -p` 后重跑，single_cell、
   poisson_invariance、transaction、multibin n1/n2/n5 全部 `status=PASS`，
   §17.10.1 第 4 步验收完成。

**历史第 6 步：Beam 12 fs 阻断回归。** 使用新目录，population controller 必须
关闭。

```bash
rm -rf ./output/hybrid_h9_threshold_fixed_beam_12fs
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 10.5,12 \
  --output-dir ./output/hybrid_h9_threshold_fixed_beam_12fs
```

本段是旧static-cell路线的历史计划。后续17A--17F通量式转换、25.5 fs专项复现和
40 fs controller-off均已完成并通过，不得再按本段重复第0--13步。controller-on算法
在重新设计并通过独立A/B前仍不得用于生产；最新下一步以§0.9为准。

#### 17.10.2 17A–17F 通量式转换验收命令 [已完成归档]

本节命令面向星逸集群，工作目录固定为：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
```

新增测试程序必须支持 `--case all --result <file>`；成功时结果文件末尾必须有
`status=PASS`，失败必须返回非零退出码。测试不得复写生产公式，必须直接调用本节新增的
生产类或生产函数。所有输出目录先显式创建；不得因为目录不存在而只在rank 0失败。
当前所有可执行文件统一从 `./build/` 运行，构建命令统一使用 `cmake --build build`。
生产 `fp_solver` 不支持 `--overwrite-output`；每个重跑命令必须先删除并重新创建其
目标输出目录，禁止在已有目录上追加或混写结果。

**17A：接口、parcel和remap sink纯函数。**

实现文件：

```text
src/grid.h
src/grid.cpp（若当前网格实现不在头文件中）
src/bulk_tail_flux_parcel.h
src/bulk_tail_flux_parcel.cpp
src/conservative_ppm_remap.h
src/conservative_ppm_remap.cpp
src/vlasov_split_step.h
src/vlasov_split_step.cpp
tests/bulk_tail_flux_interface_test.cpp
tests/upar_flux_parcel_test.cpp
tests/upar_flux_sink_test.cpp
CMakeLists.txt
```

编译和运行：

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON

cmake --build build -j4 --target \
  bulk_tail_flux_interface_test \
  upar_flux_parcel_test \
  upar_flux_sink_test \
  conservative_upar_remap_test

mkdir -p ./output/h9_flux_conversion_17A

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_interface_test \
  --case all \
  --result ./output/h9_flux_conversion_17A/interface.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/upar_flux_parcel_test \
  --case all \
  --result ./output/h9_flux_conversion_17A/parcel.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/upar_flux_sink_test \
  --case all \
  --result ./output/h9_flux_conversion_17A/sink.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/conservative_upar_remap_test \
  --case all \
  --result ./output/h9_flux_conversion_17A/remap_regression.result

grep -H '^status=' ./output/h9_flux_conversion_17A/*.result
```

17A只有在四个结果均PASS，并同时满足以下逐项字段时才能结束：

```text
interface_duplicate_count=0
interface_hole_count=0
axis_interface_count=0
below_threshold_number_relative<=1e-12
negative_node_count=0
quadrature_4_vs_8_relative_max<=1e-11
sink_number_relative_error<=1e-13
null_interface_bitwise_equal=1
```

若 `null_interface_bitwise_equal=0`，说明新增接口改变了原remap路径，必须先修复；不得以
“误差很小”替代逐bit回归门。

**17B：只读 `flux-audit` 接入。**

实现文件：

```text
src/vpfp_integrator.h/.cpp
src/main_vpfp.cpp
src/vpfp_diagnostics.h/.cpp
src/vpfp_checkpoint.h/.cpp
tools/compare_flux_audit_ab.py
```

`flux-audit` 只能生成parcel和诊断，不能改变bulk sink、tail粒子、场、ID、随机数或累计
账本。实现使用预分配的面通量副本核对观察前后字节，不重放完整$u_\parallel$算子；因此
不会把审计成本扩大为双倍速度推进。17B不要求checkpoint；若没有兼容checkpoint，直接
从0 fs分别运行static与audit。
两次运行必须使用同一可执行文件、相同80-rank布局和完全相同的物理参数；两次运行之间
禁止重新编译。当前程序没有独立Beam seed CLI，因此必须保留现有确定性默认初始化。
先编译生产求解器：

```bash
cmake --build build -j4 --target fp_solver
rm -rf ./output/h9_flux_conversion_17B/static ./output/h9_flux_conversion_17B/audit
mkdir -p ./output/h9_flux_conversion_17B/static ./output/h9_flux_conversion_17B/audit
```

从0 fs运行两次；除conversion mode、audit专用参数和输出目录外，其余参数必须相同：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode static-cell \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 12 \
  --output-dir ./output/h9_flux_conversion_17B/static

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-audit \
  --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 \
  --tail-population-control-interval 0 \
  --stop-time-fs 12 --diagnostic-level 2 \
  --diagnostic-interval 1 --snapshot-times 12 \
  --output-dir ./output/h9_flux_conversion_17B/audit

python3 ./tools/compare_flux_audit_ab.py \
  --static ./output/h9_flux_conversion_17B/static \
  --audit ./output/h9_flux_conversion_17B/audit \
  --result ./output/h9_flux_conversion_17B/compare.result

grep -H '^status=' ./output/h9_flux_conversion_17B/*.result
```

比较脚本将跨独立运行的完整接受态hash保留为信息项，严格验收采用同次调用原位门。验收字段：

```text
audit_inplace_state_bitwise_equal=1
audit_inplace_rng_equal=1
audit_inplace_ledger_equal=1
accepted_step_count_equal=1
parcel_vs_face_flux_relative_l1<=1e-13
global_parcel_vs_interface_flux_relative<=1e-12
interface_duplicate_count=0
below_threshold_number_relative<=1e-12
audit_wall_overhead_relative<=0.15
status=PASS
```

旧静态抽取数与parcel数仍需写入诊断，但二者不是同一个离散对象，差异不参与PASS/FAIL。
`state_hash_equal`不再单独决定PASS/FAIL；若三项原位门任一为0，再查audit路径对面通量、
状态、ID/RNG或累计账本的写副作用。

**17C：loader和失败事务独立闭合。**

```bash
cmake --build build -j4 --target \
  bulk_tail_flux_loader_test \
  bulk_tail_flux_transaction_test \
  bulk_tail_flux_loader_mpi_test

mkdir -p ./output/h9_flux_conversion_17C

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_loader_test \
  --case all --result ./output/h9_flux_conversion_17C/loader.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_transaction_test \
  --case all --result ./output/h9_flux_conversion_17C/transaction.result

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_flux_loader_mpi_test \
    --case all \
    --result "./output/h9_flux_conversion_17C/loader_mpi_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_flux_conversion_17C/*.result
```

必须报告并通过：

```text
number_relative_l1<=1e-10
px_relative_l1<=1e-10
kinetic_energy_relative_l1<=1e-10
jx_relative_l1<=1e-9
pixx_relative_l1<=1e-9
piperp_relative_l1<=1e-9
negative_weight_count=0
duplicate_id_count=0
max_compressed_supports<=7
failed_trial_state_unchanged=1
failed_trial_next_id_unchanged=1
mpi_global_moment_equal=1
```

**17D：checkpoint/restart。**

新增测试必须保存并恢复conversion mode、接口hash、求积阶数、最大支撑数、tail粒子、
下一ID和累计conversion账本：

17D使用的是下列两个 **flux-conversion 专用** 测试；不要以 H6 的
`checkpoint_restart_equivalence_test` 代替，也不要向后者传递17D的
`--work-dir`、`--restart-dir`、`--checkpoint-dir`或`--output-dir`参数。
所有17D命令统一使用拼写 `--workdir`。专用测试目前也兼容旧别名
`--work-dir`，但运行手册不再使用该别名，避免与通用测试的参数接口混淆。

```bash
cmake --build build -j4 --target \
  bulk_tail_flux_checkpoint_test \
  bulk_tail_flux_restart_mpi_test

mkdir -p ./output/h9_flux_conversion_17D

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/bulk_tail_flux_checkpoint_test \
  --case all \
  --workdir ./output/h9_flux_conversion_17D/checkpoint_work \
  --result ./output/h9_flux_conversion_17D/checkpoint.result

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/bulk_tail_flux_restart_mpi_test \
    --case all \
    --workdir "./output/h9_flux_conversion_17D/restart_n${NP}" \
    --result "./output/h9_flux_conversion_17D/restart_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_flux_conversion_17D/*.result
```

除所有 `status=PASS` 外，还要求配置不匹配重启被明确拒绝，且一步连续运行与restart运行的
完整状态hash一致。测试失败时禁止通过checkpoint读取后重新采样tail粒子来规避。

**17E：无碰撞生产A/B。**

不能把用 `collision-model=moment-closure` 生成的10.4 fs checkpoint改参数重启为
`collision-model=none`。必须先生成一批无碰撞checkpoint，并从诊断中选择“首次
conversion之前最后一个已接受checkpoint”。不要假定无碰撞首次转换仍恰好在10.4 fs。

选择器以 checkpoint manifest 的 `conversion_cumulative_number` 判定转换前状态。
因此生成 checkpoint 的二进制必须包含该 manifest 字段；旧 checkpoint 缺少该字段时会被
保守拒绝，不能手工改写 manifest。若 locator 已确认至某时间 $t_{\rm safe}$ conversion
累计数与 tail 粒子数均为零，只需用同一无碰撞配置重跑到 $t_{\rm safe}$ 并保存一个新的
checkpoint，无需重跑更晚的候选时间。

```bash
cmake --build build -j4 --target fp_solver
rm -rf ./output/h9_flux_conversion_17E/preconversion_clean
mkdir -p ./output/h9_flux_conversion_17E/preconversion_clean

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model none \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode static-cell \
  --tail-population-control-interval 0 \
  --checkpoint-times 9 \
  --stop-time-fs 9 --diagnostic-level 1 \
  --output-dir ./output/h9_flux_conversion_17E/preconversion_clean

python3 ./tools/select_preconversion_checkpoint.py \
  --input ./output/h9_flux_conversion_17E/preconversion_clean \
  --require-conversion-number-zero \
  --result ./output/h9_flux_conversion_17E/preconversion_selection.result

PRECONVERSION_CHECKPOINT=$(awk -F= \
  '$1=="selected_checkpoint" {print $2}' \
  ./output/h9_flux_conversion_17E/preconversion_selection.result)

test -n "$PRECONVERSION_CHECKPOINT" || { echo 'no valid preconversion checkpoint'; exit 2; }
```

17E 专用 restart 规则：当且仅当 checkpoint 的累计 conversion $N/P_x/K$、
已创建 tail 数、tail 出流、人口控制账本及全局 tail 粒子数均严格为零时，生产程序允许
`static-cell -> flux-interface` 切换。该例外只用于尚未发生表示转换的共同初态；任何
post-conversion checkpoint 的模式切换仍必须以 `tail conversion mode mismatch` 失败。
因此修改该门后只需重新编译 `fp_solver`，已有通过 selector 的 9 fs checkpoint 无需重跑。

然后从同一无碰撞checkpoint推进到12 fs：

```bash
rm -rf ./output/h9_flux_conversion_17E/static ./output/h9_flux_conversion_17E/flux
mkdir -p ./output/h9_flux_conversion_17E/static ./output/h9_flux_conversion_17E/flux

for MODE in static-cell flux-interface; do
  if [ "$MODE" = "static-cell" ]; then OUT=static; else OUT=flux; fi
  yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
    stdbuf -oL -eL ./build/fp_solver \
    --restart-dir "$PRECONVERSION_CHECKPOINT" \
    --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
    --beam-enabled 1 --collision-model none \
    --background-tail-mode pic \
    --tail-convert-energy-mev 6.0 --tail-return-mode none \
    --tail-conversion-mode "$MODE" \
    --tail-flux-quadrature-order 4 \
    --tail-flux-max-supports 7 \
    --tail-population-control-interval 0 \
    --stop-time-fs 12 --diagnostic-level 2 \
    --diagnostic-interval 1 --snapshot-times 12 \
    --output-dir "./output/h9_flux_conversion_17E/${OUT}" \
    || exit 1
done

python3 ./tools/compare_flux_conversion_production.py \
  --static ./output/h9_flux_conversion_17E/static \
  --flux ./output/h9_flux_conversion_17E/flux \
  --result ./output/h9_flux_conversion_17E/compare.result

grep -H '^status=' ./output/h9_flux_conversion_17E/*.result
```

17E验收不是要求A/B波形逐bit相同；转换数学定义已经改变。必须满足：

```text
flux_all_steps_accepted=1
flux_static_extractor_call_count=0
combined_number_relative_error<=1e-10
conversion_N_Px_K_relative_l1<=1e-10
conversion_Jx_Pixx_Piperp_relative_l1<=1e-9
conversion_piperp_old_residual_removed=1
below_threshold_number_relative<=1e-12
interface_duplicate_count=0
gauss_error_not_regressed=1
tail_particle_count_ratio_flux_over_static<=2
step_wall_time_ratio_flux_over_static<=2
```

若17E失败，只修复电场 `u_parallel` remap、parcel或loader，不得提前修改碰撞算子。

**17F：碰撞保守面通量及完整H9回归 [独立测试已通过]。**

当前碰撞路径已使用实际最终保守速度面通量；同一面通量同时更新 bulk 并生成 parcel，
未使用碰撞前后 cell 差构造伪 parcel。下方命令保留为可重复的回归，而非待开发项。

```bash
cmake --build build -j4 --target \
  cylindrical_collision_flux_test \
  collision_flux_interface_test \
  collision_flux_interface_cross_test \
  collision_flux_conversion_transaction_test \
  collision_flux_interface_mpi_test

mkdir -p ./output/h9_flux_conversion_17F

# This is the required single-rank aggregate gate.  Unlike the documented
# wrappers below, this target has no VPFP_COLLISION_FLUX_CASE definition, so
# --case all actually executes zero-flux, interface, cross, transaction,
# pure drift, parallel/perpendicular diffusion, and both collision halves.
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/cylindrical_collision_flux_test \
  --case all --result ./output/h9_flux_conversion_17F/all.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/collision_flux_interface_test \
  --case all --result ./output/h9_flux_conversion_17F/interface.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/collision_flux_interface_cross_test \
  --case all --result ./output/h9_flux_conversion_17F/cross.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/collision_flux_conversion_transaction_test \
  --case all --result ./output/h9_flux_conversion_17F/transaction.result

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/collision_flux_interface_mpi_test \
    --case all \
    --result "./output/h9_flux_conversion_17F/mpi_n${NP}.result" \
    || exit 1
done

grep -H '^status=' ./output/h9_flux_conversion_17F/*.result
```

独立碰撞门要求：零系数严格零parcel；纯漂移、平行扩散和垂直扩散质量闭合；交叉扩散
共享面全局和达到确定性求和精度；两个collision half分别闭合；失败事务不污染接受态；
1/2/5 rank全局六矩一致。`all.result` 必须同时满足以下字段均为 `1`，单个固定-case
wrapper 的 `status=PASS` 只证明其对应专项，不可替代该总门：

```text
zero_flux_pass
interface_export_pass
cross_flux_pass
pure_drift_pass
parallel_diffusion_pass
perp_diffusion_pass
collision_half_1_pass
collision_half_2_pass
transaction_rollback_pass
```

`mpi_n1.result`、`mpi_n2.result` 与 `mpi_n5.result` 的 `mpi_partition_pass=1` 仍是
独立 MPI 门。全部通过后，再按 §17.10 原H9命令重跑，但必须增加：

```text
--tail-conversion-mode flux-interface
--tail-flux-quadrature-order 4
--tail-flux-max-supports 7
```

并在manifest中验收：

```text
tail_conversion_mode=flux-interface
collision_induced_conversion=1
population_control_interval=0
```

最后依次重跑：collision pair单元测试、反作用闭合、Beam 12 fs、no-Beam 40 fs。
任何一步失败都停止；不得直接跳到Beam 40 fs或120 fs。

### 17.11 120 fs 后置验收与H10决策 [已完成]

本节已完成120 fs基线、边界能量单步接线和低能Tail只读审计，现保留为历史与决策记录。
§17.1--§17.10同样为已完成测试和历史追溯，不得重新从H1开始顺序执行。
当前状态链为：

```text
Gate A: 最新flux-interface阈值统一能谱 [已通过]
  -> Gate B: 40 fs准入资源外推 [已通过，但长期外推低估]
  -> Gate C: 40 -> 120 fs连续生产 [已完成]
  -> §17.11.6单步边界能量和低能Tail只读审计 [已完成]
  -> §16 / §17.13 H10迟滞返回R0--R3 [已通过]
  -> §16.6.5 / §17.13.7 H10 R4参数敏感性 [已通过]
  -> §17.14 H10连续准生产 [资源通过/长期物理未通过]
  -> §17.15完整能量账与返回表示误差短A/B [当前下一项]
```

与全文章节的对应关系如下：

| 顺序 | 文档任务 | 当前状态 | 执行位置 |
|---:|---|---|---|
| 1 | §7.11.6第8项 / §19.2：6 MeV附近combined能谱与0.1/0.2/0.4 MeV分箱审计 | 已通过 | §17.11.1--17.11.2 |
| 2 | §18.3：120 fs粒子数、单rank内存、负载和耗时外推 | 已通过 | §17.11.3 |
| 3 | §7.11.17.10：登记阈值物理门和资源门结果 | 已完成 | §17.11.4 |
| 4 | §17.12 / H11：40--120 fs连续生产 | 已完成 | §17.11.5 |
| 5 | §17.11.6：边界能量单步和低能Tail checkpoint审计 | 已通过 | §17.11.6 |
| 6 | §16 / H10：双向返回R0--R3 | 已通过 | §17.13.1--17.13.6 |
| 7 | H10 R4参数敏感性和准生产门 | 已通过 | §17.13.7 |
| 8 | H10中央档100.3--120 fs连续准生产 | 资源通过/长期物理未通过 | §17.14 |
| 9 | 完整能量恒等式与返回六矩有符号差短A/B | 当前下一项 | §17.15 |
| 10 | §19.4：能量、能谱、动量和宏观波包最终验收 | 阻塞于第9项 | §19.4 |

H1--H9、17A--17F及H9的Beam/no-Beam推进门均已完成，禁止把它们重新加入当前待办
清单。若只需继续项目，应进入§17.15的100.3--105 fs短A/B，而不是重跑R0--R4、
§17.3、§17.11.5或新的0--120 fs生产。

任一Gate失败即停止，不得跳过，也不得通过开启population controller、放宽转换残差、
修改6 MeV阈值或增加速度网格来使结果“通过”。

#### 17.11.1 一次性补充只读审计目标 [已完成]

现有 `energy_spectrum_rank*.dat` 只有256个全域对数bin；
`tail_threshold_interface_rank*.dat` 只有0.2 MeV阈值bin。它们可复核0.2/0.4 MeV，
但不能可靠构造0.1 MeV结果，也不包含低能Tail驻留时间。因此禁止用插值或把粗bin均分
伪造§7.11.6第8项。先新增一个只读目标：

```text
src/hybrid_checkpoint_gate_audit.cpp
build/hybrid_checkpoint_gate_audit
```

该目标必须复用生产checkpoint读取器、圆柱速度cell体积权重和Tail粒子能量定义；不得
复写另一套物理公式。它只读40 fs checkpoint以及25/40 fs snapshot，不得推进状态、修改
RNG或写回checkpoint。CLI必须为：

```text
--checkpoint <40fs完整checkpoint目录>
--snapshot25 <25fs快照目录>
--snapshot40 <40fs快照目录>
--kout-mev 6.0
--kin-mev <预声明的候选返回阈值>
--bin-widths-mev 0.1,0.2,0.4
--core-noise-limit 1e-3
--project-to-fs 120
--result <result文件>
```

输出至少包含：

```text
status
spectrum_convergence_available
spectrum_gate_pass
spectrum_rebin_metrics_informational
spectrum_rebin_shape_pass
bin_width_0p1_pass / bin_width_0p2_pass / bin_width_0p4_pass
threshold_l1_0p1 / threshold_l1_0p2 / threshold_l1_0p4
threshold_max_bin_relative
threshold_compared_non_sparse_bin_count
threshold_skipped_sparse_bin_count
threshold_reference_bin_width_mev
threshold_artificial_gap_count
threshold_deep_valley_count
threshold_max_adjacent_jump_ratio
threshold_artificial_gap_pass
threshold_spectrum_file
tail_physical_number_fraction
thermalized_tail_fraction
thermalized_tail_residence_available
core_tail_density_noise_relative
core_tail_density_noise_relative_to_tail_envelope
core_tail_density_rms_relative_to_background
core_tail_density_max_relative_to_background
core_tail_density_noise_limit
core_tail_density_noise_pass
tail_physical_number_projected_120fs
tail_particles_projected_120fs
tail_particles_local_projected_120fs
tail_particles_checkpoint_global
tail_particles_checkpoint_local_max
tail_particles_budget_global
tail_particles_budget_local
memory_projected_global_bytes
memory_projected_max_rank_bytes
wall_seconds_projected_40_to_120fs
wall_seconds_per_step_projected
performance_dt_fs
particle_projection_available
particle_budget_pass
performance_projection_available
resource_projection_available
resource_gate_pass
h10_required
population_controller_required
```

实现约束：细谱参考必须使用0.05 MeV，并与全部候选谱一样先做MPI全局归约；不得拿
rank-local细谱和global粗谱比较。Eulerian bulk是速度单元积分质量，构造能谱时必须
复用`TailSubcellQuadrature::energy_bin_fractions()`，解析计算圆柱速度单元与能量
壳层的分段体积交叠；禁止把整个单元质量集中到cell-center能量，也禁止用固定
$4\times4$点积生成细直方图。PIC Tail仍按粒子实际能量入bin。
0.1/0.2/0.4 MeV谱必须由0.05 MeV有限体积细谱守恒聚合，不能重复扫描完整分布。

必须输出`<result>.spectrum.dat`，其中分别包含`fine_bulk_0p05`、`fine_tail_0p05`和
`fine_combined_0p05`。人工结构检查必须同时覆盖：单bin孤立谷、由多个连续bin组成且
两侧恢复的深谷，以及相邻非稀疏bin超过100倍的异常跃变；不能只检查单个中间bin。
`threshold_max_bin_relative`只统计阈值窗口内的非稀疏bin。

同一checkpoint重新分成0.1/0.2/0.4 MeV直方图，不能证明生产离散随网格加密收敛；
它只衡量同一谱对后处理分箱的敏感性。因此这些$L_1$量必须标记为informational，
`spectrum_convergence_available=0`、`spectrum_gate_pass=-1`，不得因其超过5%直接判定
生产转换算法失败。Gate A不再要求重复执行昂贵的0--40 fs全程序A/B。生产
`u_parallel`转换算子的4/8阶求积收敛由`upar_flux_sink_test`直接调用同一生产remap
实现验收，真实40 fs checkpoint负责检验combined谱形。两类证据必须由
`tools/validate_hybrid_gate_a.py`联合汇总，不能用任意一类单独宣称Gate A通过。

`tail_particles_*`固定表示宏粒子数，Tail物理数必须写入独立的
`tail_physical_number_*`字段。核心噪声使用中央80%区域相对21-cell平滑包络的高频
残差，并以同区域背景密度RMS归一化作为物理门；相对Tail自身包络的比值只保留为
稀疏度/粗糙度诊断，不得据此判断核心区受污染。默认物理门为
`core_tail_density_noise_relative <= 1e-3`，命令中必须显式写出该阈值。

`tail_particles_budget_global/local`和`particle_budget_pass`是保守规划告警，不是物理门，
也不能替代真实硬件资源判定。Gate B必须使用`memory_projected_max_rank_bytes`、全局内存、
wall-time外推、作业实际内存上限和队列时限判断可运行性。固定的$1.2\times10^7$全局和
$2.0\times10^6$单rank宏粒子阈值不得单独令Gate B失败，也不得据此强制开启population
controller。`performance_projection_available`只表示wall-time外推是否可用；背景总粒子数
必须调用生产`Species::total_particle_number()`，不得在审计工具中复写速度空间积分公式。

若现有checkpoint没有逐粒子驻留计数，仍应输出
`thermalized_tail_residence_available=0`，但该字段只是可选机理诊断。Gate B的H10判据使用
当前Tail物理数分数、低能Tail瞬时比例和核心噪声；三者均通过时可判定
`h10_required=0`。审计工具旧版输出的`status=INCOMPLETE`和`h10_required=-1`只表示其
综合状态机尚未按本节更新，不得覆盖§17.11.3的正式Gate B判定。

编译命令：

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j4 --target \
  tail_energy_shell_overlap_test upar_flux_sink_test hybrid_checkpoint_gate_audit
mkdir -p ./output/h9_flux_interface_17f_post40_audit
./build/tail_energy_shell_overlap_test \
  --result ./output/h9_flux_interface_17f_post40_audit/shell_overlap.result
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/upar_flux_sink_test \
  --result ./output/h9_flux_interface_17f_post40_audit/conversion_q4_q8.result
```

#### 17.11.2 Gate A：阈值能谱命令 [已通过]

以下路径在服务器上按实际目录替换；`CHECKPOINT_40`必须是完整checkpoint目录，不是
snapshot目录。

```bash
ROOT=/HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
RUN40="$ROOT/output/h9_flux_interface_17f_beam_40fs"
CHECKPOINT_40="$RUN40/<40fs完整checkpoint目录>"
SNAP25="$RUN40/snapshot_t25.0222154fs_step978"
SNAP40="$RUN40/snapshot_t40fs_step1564"
AUDIT_OUT="$ROOT/output/h9_flux_interface_17f_post40_audit"
KIN_MEV=5.6  # 首轮迟滞候选，仅用于必要性审计，不是已验证物理常数

mkdir -p "$AUDIT_OUT"
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  "$ROOT/build/hybrid_checkpoint_gate_audit" \
  --checkpoint "$CHECKPOINT_40" \
  --snapshot25 "$SNAP25" --snapshot40 "$SNAP40" \
  --kout-mev 6.0 --kin-mev "$KIN_MEV" \
  --bin-widths-mev 0.1,0.2,0.4 \
  --core-noise-limit 1e-3 \
  --project-to-fs 120 \
  --result "$AUDIT_OUT/gates.result"

cat "$AUDIT_OUT/gates.result"

python3 "$ROOT/tools/validate_hybrid_gate_a.py" \
  --geometry-result "$AUDIT_OUT/shell_overlap.result" \
  --checkpoint-audit "$AUDIT_OUT/gates.result" \
  --conversion-result "$AUDIT_OUT/conversion_q4_q8.result" \
  --result "$AUDIT_OUT/gate_a.result"

cat "$AUDIT_OUT/gate_a.result"
```

第一条MPI命令完成固定checkpoint谱形审计；最后的汇总器才给出Gate A最终结论。
固定checkpoint审计必须满足：

```text
audit_schema=hybrid_checkpoint_gate_audit_v5
threshold_artificial_gap_pass=1
threshold_deep_valley_count=0
checkpoint_audit_finite=1
```

`spectrum_convergence_available=0`、`spectrum_gate_pass=-1`和总体`status=INCOMPLETE`
是当前单checkpoint调用的预期结果，不等于Gate A失败。0.1/0.2/0.4 MeV重新分箱的
$L_1$与最大bin差仅为信息项，不能作为生产算子FAIL条件。最终`gate_a.result`必须满足：

```text
gate_a_schema=hybrid_gate_a_v1
status=PASS
geometry_overlap_pass=1
geometry_mass_partition_pass=1
conversion_operator_pass=1
spectrum_shape_pass=1
spectrum_number_conservation_pass=1
core_noise_pass=1
```

其中4/8阶比较直接调用生产`ConservativePpmRemap::advect_u_parallel()`，检验转换接口
求积、sink守恒、重复面和负权；40 fs真实checkpoint检验6 MeV附近combined谱。该组合
替代重复0--40 fs的全程序参数A/B，但不替代Gate C的连续推进验收。

当前验收记录为：`gate_a_schema=hybrid_gate_a_v1`、`status=PASS`；4/8阶求积相对差
为$3.19516\times10^{-13}$，谱数目相对残差为$1.68517\times10^{-13}$，人工深谷数为0，
核心Tail噪声相对背景为$1.16656\times10^{-5}$。旧17A结果未包含`tiny_tail_ok`字段，
汇总器将其记为`tiny_tail_evidence_available=0`而不否决谱门；若使用新版
`upar_flux_sink_test`重跑并提供该字段，则字段存在时必须为1。

#### 17.11.3 Gate B：40 fs准入资源与H10判定 [准入通过/长期结论已被更新]

Gate B使用同一个 `gates.result`，不再提交一次生产作业。判定顺序为：

1. `tail_physical_number_fraction < 0.01`；
2. `thermalized_tail_fraction < 0.10`；
3. 低能Tail未在核心区造成不可接受的密度/电流噪声；
4. 120 fs全局内存、单rank内存和wall-time投影低于实际作业资源硬上限。

当前结果为：

```text
tail_physical_number_fraction=1.8045267899599466e-06       PASS
thermalized_tail_fraction=6.9459650846039642e-04          PASS
core_tail_density_noise_relative=1.1665626963395677e-05   PASS
memory_projected_global_bytes=3150268544                   PASS
memory_projected_max_rank_bytes=276769018.47706789         PASS
wall_seconds_projected_40_to_120fs=15997.605687648609      PASS
```

约3.15 GB全局内存、277 MB最大rank内存和4.44小时预计耗时均低于当前5节点作业的实际
资源上限。预计$4.83\times10^7$个全局宏粒子和$5.55\times10^6$个最大rank宏粒子会
造成负载不均衡并值得生产期间监控，但它们不是内存溢出或队列超时证据。审计工具内
硬编码的$1.2\times10^7/2.0\times10^6$规划阈值不属于已声明硬件限制，因此本次不采用
`particle_budget_pass=0`和`population_controller_required=1`作为阻断判据。

以下是启动Gate C时使用的前瞻分支，不是120 fs最终判定：

```text
h10_required=0, population_controller_required=0
  -> 当时将H10暂记为“不阻断Gate C”；120 fs后必须重新评估。

h10_required=1
  -> 停止；执行§16双向返回，未通过前禁止120 fs。

h10_required=0, population_controller_required=1
  -> 停止；只做§7.10人口控制物理A/B。宏粒子多但低能Tail少时，禁止误用H10。
```

当时的正式准入判定为：

```text
gate_b_status=PASS
h10_required=0
population_controller_required=0
resource_gate_pass=1
```

这里`resource_gate_pass=1`是按实际内存和wall-time上限作出的文档级正式判定；在审计
代码同步修改前，旧`gates.result`中的同名字段仍可能为0，不得混用。

120 fs真实运行没有否定“作业可以在现有资源内完成”，但否定了简单增长外推和
`h10_required=0`的长期有效性：全局Tail宏粒子实际达到$1.4317\times10^8$，低于6 MeV
的Tail数目占比达到58.80%。因此Gate B仍记为**生产准入通过**；该证据随后触发H10开发，
目前R0--R4已经通过，最终状态以§0.11、§15.13和§17.13为准。

#### 17.11.4 更新§7.11.17.10完成状态 [已完成]

该步骤当时将H9生产门、阈值物理门和资源准入门登记为通过，并据40 fs证据暂记H10
不需要。120 fs结果随后重新触发H10，现已完成R0--R4；该状态变化不否定已经通过的
flux-interface算子。

#### 17.11.5 Gate C：从40 fs checkpoint推进到120 fs [已完成]

只有以下三项均满足才运行；不得再要求单checkpoint的`gates.result`伪造
`spectrum_gate_pass=1`：

1. `gate_a.result`满足`status=PASS`；
2. §17.11.3的正式记录满足`gate_b_status=PASS`和`h10_required=0`；
3. 按实际资源上限判定`population_controller_required=0`。

```bash
ROOT=/HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
RUN40="$ROOT/output/h9_flux_interface_17f_beam_40fs"
CHECKPOINT_40="$RUN40/<40fs完整checkpoint目录>"
OUT120="$ROOT/output/hybrid_flux_interface_controller_off_40_to_120fs"

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

rm -rf "$OUT120"
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL "$ROOT/build/fp_solver" \
  --restart-dir "$CHECKPOINT_40" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 120 --checkpoint-times 60,80,100,120 \
  --snapshot-times 50,60,80,100,120 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --checkpoint-dir "$OUT120/checkpoints" \
  --output-dir "$OUT120"
```

正式生产不得启用 `--tail-stage-trace`。运行期间只在60、80、100 fs检查有限性、
Tail预算和wall-time趋势；未触发硬门时不要中途改变算法或参数。完成后执行§19.4，输出
电场能、combined背景动能/平均动能、Beam能量、总能量账、combined能谱、动量分布及
波包包络。EPOCH只比较宏观包络、能量分配和统计量，不要求PIC噪声逐点一致。

实际运行结果：step 1565--4691共3127步全部接受，0 split、0 rollback，最终到达
120 fs。转换和Tail数目守恒通过，实际耗时约5.99 h；Tail宏粒子最终为143,165,344，
最大rank为5,301,652。该命令保留为可复现实验记录，不得立即重复提交。

#### 17.11.6 120 fs 后置审计与下一开发项 [诊断接线与只读审计已完成]

当前不重跑120 fs，也不先修改Poisson、Beam、flux-interface转换或碰撞pair。按以下顺序：

1. 从 `vpfp_step_diagnostics.dat` 对每步累计 `collision_reservoir`，不得只读取末行；
2. 增加背景左右reservoir的数目、动量、动能流率与累计账。该诊断必须直接复用边界面
   remap通量，不得从相邻cell状态重新估算；
3. 验证
   $$
   \Delta(U_E+K_{\rm bulk}+K_{\rm tail}+K_{\rm beam})
   +\sum\Delta E_{\rm collision}
   -E_{\rm boundary,in}+E_{\rm boundary,out}
   $$
   是否收敛到时间离散、求和和PIC采样误差；
4. 用50/60/80/100/120 fs快照统计低于5.5、5.75和6 MeV的Tail数目、能量、空间分布、
   每cell有效粒子数和核心区密度/电流噪声；
5. 若低能Tail主要是碰撞热化并长期留在核心区，执行§16迟滞返回原型。推荐初始A/B为
   $K_{out}=6$ MeV，$K_{in}=5.5$ MeV，并要求连续驻留若干碰撞步后才返回；阈值必须通过
   5.25/5.5/5.75 MeV敏感性比较，不能直接固定为生产值；
6. 返回原型只在固定checkpoint短程运行中验收。必须比较combined的$N/P_x/K$、
   $\Pi_{xx}/\Pi_\perp$、6 MeV谱连续性、Poisson场、噪声、宏粒子下降率和wall time；
7. Population Controller继续关闭。它可合并宏粒子，但不能纠正“低能粒子属于哪种表示”
   的模型问题。

##### 17.11.6.1 已补充的只读/诊断代码（2026-08-12）

本轮只增加诊断，不改变Vlasov、Poisson、Beam、碰撞或flux-interface转换：

- `VpfpStepLedger`和`vpfp_step_diagnostics.dat`新增
  `background_left_inflow_energy`、`background_left_outflow_energy`、
  `background_right_inflow_energy`和`background_right_outflow_energy`；
- 四项能量直接取自两个$x$ remap半步已经计算的物理边界面能量通量，不从邻近cell
  重新估算，也不增加新的全分布扫描或MPI归约；
- 同一接受步ledger同时输出已有Beam状态中的`beam_injected_energy`和
  `beam_outflow_energy`，从而使0--25 fs注入阶段和后续开放出流也能进入完整能量账；
- 新增`low_energy_tail_checkpoint_audit`，只读完整checkpoint，默认对5.5、5.75和
  6 MeV三个阈值输出物理数、动能、宏粒子数、核心区占比、最大cell密度以及CIC逐cell
  空间分布；
- 审计输出`audit_read_only_state_unchanged=1`，并检查CIC域内份额与边界外形函数份额
  相加后重新得到原粒子数和能量。

旧40--120 fs诊断没有保存这四个边界能量分量，不能根据相邻cell或粒子数通量事后精确
重构累计边界能流。因此现有120 fs基线仍可用于低能Tail审计，但完整历史能量账必须由
未来连续运行记录，或从已有checkpoint分段重跑获得。不得伪造旧累计值。

##### 17.11.6.2 构建和最小验收命令

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx -DBUILD_TESTING=ON
cmake --build build --target \
  fp_solver low_energy_tail_audit_test low_energy_tail_checkpoint_audit -j 4

./build/low_energy_tail_audit_test
```

单元测试必须输出`status=PASS`。边界能量列采用一个完整checkpoint只推进1步验收；示例：

```bash
CHECKPOINT_120=<120fs完整checkpoint目录>
OUT_BOUNDARY=./output/post120_boundary_energy_one_step
rm -rf "$OUT_BOUNDARY"

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_120" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-return-mode none \
  --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 120.1 --stop-after-steps 4692 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir "$OUT_BOUNDARY"
```

这里`--stop-after-steps`是绝对接受步号，不是“再跑几步”；120 fs checkpoint为step 4691，
故设为4692。使用其他checkpoint时必须改为`checkpoint_step+1`。

验收要求：新一行诊断被接受；四个边界能量字段均为有限值且非负；数目连续性、Gauss、
conversion和碰撞门仍通过。零值允许，表示该方向本步没有对应能量通量。

**验收结果（2026-08-12）：通过。** 120 fs checkpoint从step 4691真实推进到
step 4692，`accepted=1`。本步四个背景边界能量为：

| 分量 | 能量 $(\mathrm{J/m^2})$ |
|---|---:|
| 左端流入 | 92.0520 |
| 左端流出 | 906.8624 |
| 右端流入 | 110.3139 |
| 右端流出 | 118.5901 |

因此本步背景边界净能量输入为$-823.09\ \mathrm{J/m^2}$。Beam注入/出流能量为零，
符合120 fs时的本步状态。`gauss_charge_residual=-1.42e-14`，Tail数目账残差
$4.20\times10^{-15}$，conversion的$N/P_x/K$相对残差分别约为
$1.25\times10^{-14}$、$9.42\times10^{-15}$和$1.64\times10^{-14}$；碰撞无rollback。
该结果验证了诊断接线，不代表旧40--120 fs历史能量账已经闭合。

##### 17.11.6.3 低能Tail完整checkpoint只读审计命令

每个checkpoint必须用保存时相同的80个MPI rank读取。分别运行60、80、100和120 fs：

```bash
AUDIT_ROOT=./output/low_energy_tail_checkpoint_audit
mkdir -p "$AUDIT_ROOT"

for T in 60 80 100 120; do
  CHECKPOINT=<对应${T}fs完整checkpoint目录>
  yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
    ./build/low_energy_tail_checkpoint_audit \
    --checkpoint "$CHECKPOINT" \
    --thresholds-mev 5.5,5.75,6.0 \
    --result "$AUDIT_ROOT/t${T}fs.result" || exit 1
done
```

每个`*.result`必须满足：

```text
status=PASS
finite=1
audit_read_only_state_unchanged=1
```

并要求每个阈值的`shape_number_residual`和`shape_energy_residual`处于MPI求和误差量级。
对应`*.result.cells.dat`是逐cell的低能Tail数目、密度、动能密度和CIC宏粒子支撑，供后续
判定低能Tail是否集中在核心区以及H10返回应在哪里生效。该审计不提供驻留时间；严格的
`tail_return_residence_steps`仍需H10为粒子增加checkpoint化状态后才能测量。

**验收结果（2026-08-12）：60/80/100/120 fs四个checkpoint全部通过。**
所有结果均为`status=PASS`、`finite=1`和`audit_read_only_state_unchanged=1`；
形函数数目残差最大$6.27\times10^{-14}$，能量残差最大
$2.55\times10^{-14}$，处于MPI求和误差量级，所有`outside_shape_number/energy=0`。

| 时间 | 全部Tail宏粒子数 | $K<6$ MeV数目占Tail | $K<6$ MeV能量占Tail |
|---:|---:|---:|---:|
| 60 fs | 27,897,968 | 25.19% | 16.57% |
| 80 fs | 58,919,504 | 50.98% | 34.13% |
| 100 fs | 97,105,848 | 58.56% | 40.90% |
| 120 fs | 143,165,344 | 58.80% | 42.21% |

$K<6$ MeV粒子的加权空间中位置从60 fs的$10.93\ \mu\mathrm{m}$移动到120 fs的
$29.13\ \mu\mathrm{m}$；120 fs的5%--95%加权区间约为$27.95$--$30.92\ \mu\mathrm{m}$。
这表明低能Tail主要位于域内移动波包，不是可以当作边界废粒子删除的成分。
审计结果已构成开发H10迟滞Tail-to-bulk返回的明确依据；实际返回仍必须经过§16的
守恒性、谱连续性和短程A/B验收。

### 17.12 120 fs 生产 [连续推进完成/最终物理验收待定]

本节表示阶段H11的最终生产状态。Gate C已经完成，不得重复提交§17.11.5。现有结果是
有效的controller-off长期参考基线。边界能量诊断的单步接线和60--120 fs低能Tail只读审计
已通过；但旧40--120 fs输出未保存边界能量历史，不能事后建立精确累计账。同时，低能Tail
占比和宏粒子数已证明H10返回需求。因此状态保持“连续推进完成/最终物理验收待定”。不得在本节维护第二份
生产命令。

### 17.13 H10迟滞Tail-to-bulk返回 [R0--R4已通过]

本节是§16的唯一命令入口。严格按R0、R1、R2、R3、R4顺序执行；前一阶段未通过时
不得提交后一阶段。下列target和工具已补齐并加入`CMakeLists.txt`。

R2早期的“投影可行比例为零”属于旧实现历史，随后已通过返回投影、事务和checkpoint修复关闭。
截至2026-08-13，R0--R4全部通过；最终状态和准生产参数以§0.11、§16.6.5和§17.13.7为准。

#### 17.13.1 构建H10生产程序和测试

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j4 --target \
  fp_solver \
  tail_bulk_return_hysteresis_test \
  tail_bulk_return_projection_test \
  tail_bulk_return_transaction_test \
  tail_bulk_return_checkpoint_test \
  tail_bulk_return_mpi_test \
  tail_bulk_return_none_regression_test \
  tail_bulk_return_checkpoint_audit
```

#### 17.13.2 R0：`none`模式零回归

**状态：已通过。** `output/h10_r0_none_regression.result`给出：

- `physical_state_bitwise_equal=1`；
- `rng_equal=1`；
- `ledger_equal=1`；
- `return_collective_count=0`；
- `status=PASS`。

```bash
./build/tail_bulk_return_none_regression_test \
  --case all \
  --result ./output/h10_r0_none_regression.result
```

验收：`status=PASS`、`physical_state_bitwise_equal=1`、`rng_equal=1`、
`ledger_equal=1`、`return_collective_count=0`。

#### 17.13.3 R1：单元、事务、checkpoint和MPI回归

**状态：已全部通过。** 当前结果包括：

- 迟滞驻留测试：阈值上方计数正确复位，满足驻留步数后正确删除，`status=PASS`；
- 六矩投影测试：`moment_residual_max=5.2083333333333335e-16`，bulk非负，
  不可行样例能够延迟而不是错误提交；
- 事务测试：故障状态被拒绝，`accepted_state_unchanged=1`；
- checkpoint测试：驻留计数、返回配置和累计账往返一致，旧v1兼容路径已覆盖；
- 真实MPICH的1/2/5 rank测试均`status=PASS`，`decomposition_invariant=1`、
  `duplicate_id_count=0`、`request_balance_error=0`，最大六矩残差不超过
  $1.735\times10^{-18}$。

这些结果证明H10的合成可行样例、事务协议、checkpoint持久化和MPI分解路径成立，但不证明
真实checkpoint中的候选矩一定存在非负离散投影。

```bash
./build/tail_bulk_return_hysteresis_test \
  --case all --result ./output/h10_r1_hysteresis.result
```

```bash
./build/tail_bulk_return_projection_test \
  --case all --result ./output/h10_r1_projection.result
```

```bash
./build/tail_bulk_return_transaction_test \
  --case all --result ./output/h10_r1_transaction.result
```

```bash
./build/tail_bulk_return_checkpoint_test \
  --case all --workdir ./output/h10_r1_checkpoint_work \
  --result ./output/h10_r1_checkpoint.result
```

```bash
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_mpi_test \
  --case all --result ./output/h10_r1_mpi_n1.result
```

```bash
yhrun -N 1 -n 2 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_mpi_test \
  --case all --result ./output/h10_r1_mpi_n2.result
```

```bash
yhrun -N 1 -n 5 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_mpi_test \
  --case all --result ./output/h10_r1_mpi_n5.result
```

每个`*.result`均必须`status=PASS`。投影守恒不变量$N/P_x/K$最大相对残差$\le10^{-12}$；
$J_x/\Pi_{xx}/\Pi_\perp$继续输出为速度网格表示误差，不作为逐组精确约束。MPI三档必须
`decomposition_invariant=1`、`duplicate_id_count=0`、`request_balance_error<=1e-13`；故障注入必须
`accepted_state_unchanged=1`。

#### 17.13.4 R2：完整checkpoint只读可行性审计

**状态：审计执行层通过，物理可行性门失败。** 四个结果文件虽然都写出`status=PASS`，该字段
当前只表示读取、有限性、生产投影调用和只读不变性通过，不能替代下面的
`feasible_candidate_number_fraction`验收。

| checkpoint | 5.5 MeV候选数（radius=3） | 返回数 | 可行比例 | 结论 |
|---|---:|---:|---:|---|
| 60 fs | $8.0447\times10^{18}$ | 0 | 0 | 不通过 |
| 80 fs | $3.2941\times10^{19}$ | 0 | 0 | 不通过 |
| 100 fs | $4.4084\times10^{19}$ | 0 | 0 | 不通过 |
| 120 fs | $4.7927\times10^{19}$ | 0 | 0 | 不通过 |

而且这不是中央参数的偶然现象：5.25、5.5、5.75 MeV和radius=1、2、3的全部36个组合均
`returned_number=0`、`feasible_number_fraction=0`。所有候选均留在Tail中；候选主要位于核心区，
不是空间边界样本导致的假失败。四次审计均`finite=1`且
`audit_read_only_state_unchanged=1`，所以该结论不是checkpoint被审计器修改所致。

根因定位为旧投影契约不可行：它把同一CIC cell pair内的候选作为不可拆事务组，并要求在固定
cell-centre二维速度网格上以$10^{-12}$同时重构六矩。R1人工样例位于该离散凸锥内，但真实
连续PIC样本的完整六矩通常不在同一个非负离散凸锥内；把radius从1扩大到3不能解决契约本身的
矛盾。因此问题不是“没有低能Tail”，而是旧返回表示约束过强，生产算子实际上不会提交。

**2026-08-13 R2 v2实测补充：** 门控语义已经按预期工作。60/80 fs因只作趋势审计仍
`status=PASS`；100/120 fs均为`feasibility_gate_required=1`、
`feasibility_gate_pass=0`、`status=FAIL`。中央档 attempted/deferred 组数分别为
3238/3238和4012/4012，committed为0，确认失败发生在每个真实空间组的投影可行性，而不是
候选扫描、MPI统计、checkpoint读取或边界样本。

上一轮四个结果缺少`projection_schema`、`invariant_residual_max`和
`representation_residual_max`，只能作为旧六矩投影契约的失败证据；它们现已被下面带实现
指纹的第二轮R2取代。

**2026-08-13第二轮R1/R2结果与修复：** 新R1七项全部`status=PASS`，包括1/2/5 rank
分解一致性；连续离网格宽云的$N/P_x/K$残差为$6.51\times10^{-14}$。新R2也已确认正确链接
`conservative_n_px_k_multianchor_v1`，但60/80/100/120 fs所有真实组仍全部defer，中央档
可行率均为0。100 fs为3238/3238组defer，120 fs为4012/4012组defer；候选数分别为
$4.4084\times10^{19}$和$4.7927\times10^{19}$。这排除了旧二进制、候选扫描、checkpoint读取、
MPI事务和六矩错误硬门，失败位于局部非负$N/P_x/K$求解本身。

根因是R1人工分布尚不足以代表真实checkpoint中的连续速度小组。原实现只使用固定12000步
投影梯度NNLS；真实稀疏多锚点矩阵条件数较差时，即使目标位于非负凸包内也可能无法在
$10^{-12}$达到停止门，而且失败前没有输出原因，因而结果中的残差为0只是“未生成请求”，
不是零误差。

现增加确定性几何后备解。把$N$约去后，$N/P_x/K$约束化为
$(u_\parallel,\gamma-1)$二维凸组合：先在目标附近选两条$u_\parallel$线，每条用上下能量cell
插值，最多四个cell精确重构；局部四点不可用时，对自适应支撑构造二维凸包并作确定性三角
分解。它不放宽$10^{-12}$门、不产生负质量，也不修改$J_x/\Pi_{xx}/\Pi_\perp$诊断。新实现指纹为
`projection_schema=conservative_n_px_k_local_histogram_v5`。

同时新增三类失败计数：`projection_invalid_input_cells`、
`projection_insufficient_support_cells`和`projection_infeasible_invariant_cells`。它们在R2中作
全MPI域归约，也进入生产接受态诊断。若v2仍有defer，可直接区分输入错误、候选支撑不足和目标
确实位于凸包外，禁止再用残差0猜测原因。

第一轮已修复投影支撑构造：由旧的“平均速度附近单一小块”改成覆盖组内均值、方差和极值的
稀疏多锚点支撑，并加入实际占用速度cell的确定性质量锚点；同时使用占用质量加微小高斯prior
控制NNLS成本。新增宽速度云回归，专门覆盖旧实现无法表达的分布。

第二轮根据真实R2仍为零可行率的证据，把数学契约修正为严格守恒$N/P_x/K$，并把
$J_x/\Pi_{xx}/\Pi_\perp$作为显式速度网格表示误差。非负性、$N/P_x/K$的$10^{-12}$门、
CIC原子组和MPI事务均保持不变。实现指纹为
`projection_schema=conservative_n_px_k_local_histogram_v5`。

本地单rank生产实现回归结果为：15个跨越多个远隔速度cell的宽云粒子全部返回，
`broad_cloud_removed=15`，bulk保持非负；该人工网格对齐样例的完整六矩最大残差为
$7.40\times10^{-14}$。这证明新契约和多锚点路径成立，但不能替代80-rank真实checkpoint R2。

下一步验证多锚点$N/P_x/K$守恒方案。新R2必须带上述`projection_schema`；缺少该字段说明运行
的仍是旧二进制，结果无效。如果可行率仍低于90%，再根据`invariant_residual_max`和
`representation_residual_max`区分主约束不可行与派生矩离散误差。不得仅放宽
`moment_tolerance`，也不得在R2未通过时继续做性能A/B。

60 fs：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_checkpoint_audit \
  --checkpoint <60fs完整checkpoint目录> \
  --return-thresholds-mev 5.25,5.5,5.75 \
  --stencil-radii 1,2,3 \
  --moment-tolerance 1e-12 \
  --minimum-feasible-fraction 0.90 \
  --result ./output/h10_r2_60fs.result
```

80 fs：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_checkpoint_audit \
  --checkpoint <80fs完整checkpoint目录> \
  --return-thresholds-mev 5.25,5.5,5.75 \
  --stencil-radii 1,2,3 \
  --moment-tolerance 1e-12 \
  --minimum-feasible-fraction 0.90 \
  --result ./output/h10_r2_80fs.result
```

100 fs：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_checkpoint_audit \
  --checkpoint <100fs完整checkpoint目录> \
  --return-thresholds-mev 5.25,5.5,5.75 \
  --stencil-radii 1,2,3 \
  --moment-tolerance 1e-12 \
  --minimum-feasible-fraction 0.90 \
  --result ./output/h10_r2_100fs.result
```

120 fs：

```bash
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_checkpoint_audit \
  --checkpoint <120fs完整checkpoint目录> \
  --return-thresholds-mev 5.25,5.5,5.75 \
  --stencil-radii 1,2,3 \
  --moment-tolerance 1e-12 \
  --minimum-feasible-fraction 0.90 \
  --result ./output/h10_r2_120fs.result
```

四个结果必须包含
`projection_schema=conservative_n_px_k_local_histogram_v5`，并满足`status=PASS`、`finite=1`、
`audit_read_only_state_unchanged=1`。缺少该投影指纹的R2结果来自旧二进制，必须判为无效，
不能用于评价当前修复。
100和120 fs在5.5 MeV、$7\times7$ stencil下必须
`feasible_candidate_number_fraction>=0.90`，否则停在R2分析不可行组。

审计器schema已升级为`tail_bulk_return_checkpoint_audit_v2`。60/80 fs只记录
`feasibility_gate_pass`供趋势分析；100/120 fs会把中央档90%门纳入进程退出状态。也就是说，
100或120 fs可行率不足时必须同时得到`feasibility_gate_required=1`、
`feasibility_gate_pass=0`、`status=FAIL`和非零退出码，不能再出现旧v1文件中“可行率为0但
status仍为PASS”的矛盾。`attempted_groups`、`committed_groups`、`deferred_groups`和
`particles_removed`现在是全MPI域总和，不再误用rank 0局部计数。

修改审计器后，只需重新编译并重跑R2，不要重跑60--120 fs生产checkpoint：

```bash
cmake --build build -j4 --target \
  tail_bulk_return_projection_test \
  tail_bulk_return_mpi_test \
  tail_bulk_return_checkpoint_audit

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_projection_test \
  --case all --result ./output/h10_r1_projection_after_multianchor.result

for NP in 1 2 5; do
  yhrun -N 1 -n "$NP" -c 4 --cpu-bind=cores \
    ./build/tail_bulk_return_mpi_test --case all \
    --result "./output/h10_r1_mpi_multianchor_n${NP}.result" || exit 1
done

grep -H '^status=' \
  ./output/h10_r1_projection_after_multianchor.result \
  ./output/h10_r1_mpi_multianchor_n*.result

# R1全部PASS后，按上面的四条命令重跑R2并统一检查。
grep -H -E '^(schema|projection_schema|feasible_candidate_number_fraction|feasibility_gate_required|feasibility_gate_pass|threshold_1_radius_3_projection_(invalid_input|insufficient_support|infeasible_invariant)_cells|status)=' \
  ./output/h10_r2_60fs.result \
  ./output/h10_r2_80fs.result \
  ./output/h10_r2_100fs.result \
  ./output/h10_r2_120fs.result
```

审计器内部复制checkpoint中的bulk/Tail状态，在副本中使低能候选满足一次审计所需的驻留条件，
并直接调用生产`TailBulkReturn::apply()`。因此它测试的是真实投影和MPI事务路径，而不是在
审计器中复写一套可行性公式；原checkpoint状态用显式字段hash复核为只读不变。

#### 17.13.5 R3a：100 fs checkpoint到100.3 fs物理A/B

**最新状态：A/B均已连续推进到100.3 fs并生成快照；背景质量账比较对象错误已修复，R3a短程物理门通过。**

快照提交协议已经更新，因此再次执行R3前必须先重编译。测试机已有有效CMake build目录时执行：

```bash
cmake --build build -j4 --target fp_solver
python3 -m py_compile tools/compare_tail_return_ab.py
```

必须使用本次编译生成的`./build/fp_solver`。旧二进制不会写`_COMPLETE`，旧R3快照也不能通过
新版完整性门。

2026-08-13重跑后，A/B均接受11步并生成100.3 fs快照，R3a比较器给出`status=PASS`：

- 返回$N/P_x/K$最大残差$4.67\times10^{-16}$，MPI请求余额为0；
- Tail宏粒子由97,695,024降至59,581,182，短档已观察到实际缩减；
- 密度相对$L_2=4.92\times10^{-5}$；场相对$L_2=5.47\times10^{-4}$；
- 场能代理相对差$1.57\times10^{-4}$，场包络相对$L_2=3.09\times10^{-4}$；
- combined能谱相对$L_2=1.33\times10^{-3}$；
- $J_x/\Pi_{xx}/\Pi_\perp$表示残差最大值0.991，必须保留报告；它不属于$N/P_x/K$硬守恒量，
  不能再由错误列名默认为0。该风险随后已由R3b/R4的动态宏观量、combined谱和参数平台联合验收关闭。

`tools/compare_tail_return_ab.py`已修复三处问题：使用诊断中的真实大小写列名；忽略已被后续
accepted step覆盖的旧`vpfp_failure.dat`；实际计算100-cell分块RMS包络和场能差，不再把
`wave_envelope_pass`简单等同于逐点场/密度门。R3a通过只允许进入R3b，不等于H10最终生产验收。

B侧已成功接受3911--3917步。期间每步约有$3.8\times10^7$个低能候选，但尚未满足8步驻留条件，故实际返回量为零。下一步首次批量返回时，旧实现用返回前的
`vlasov_diag.x_second.number_after`校验包含`tail_return.number`的期望背景数，必然触发
`failure_code=6`。现已改为比较返回后的
`global_sum(state_np1_.total_particle_number())`，并分别输出
`background_ledger_after_tail_return`和`tail_ledger_after_tail_return`失败阶段。该修复只更正最终守恒账本，不修改返回投影、粒子删除、bulk增量或场求解。

此前无有效步的旧记录如下，仅保留为历史说明。当前
`output/h10_r3a_compare.result`仅包含：

```text
baseline_accepted_step_count=0
baseline_snapshot_count=0
candidate_accepted_step_count=0
candidate_snapshot_count=0
comparison_error=baseline produced no accepted step
comparison_valid=0
status=FAIL
```

该失败的直接原因是比较输入不完整，不是密度、电场、能谱或守恒门超限：

1. `h10_r3a_compare.result`的时间戳早于A/B两个运行目录中诊断文件的最终时间戳，说明比较器在
   两个作业完成前已经执行；
2. 两个A/B目录均不存在`snapshot_*`；
3. 两侧`vpfp_step_diagnostics.dat`以及四个通量诊断文件都只有表头，没有任何accepted-step行，
   因而短跑没有留下可比较的已接受物理状态；
4. 比较器在读取基线快照时即退出，尚未计算`conservation_pass`、`spectrum_pass`、
   `field_pass`、`wave_envelope_pass`或粒子数缩减。

`tools/compare_tail_return_ab.py`已经修复为先验收两侧accepted-step和快照完整性，并在失败
结果中保留四个计数及`comparison_valid=0`。比较器的守恒硬门只检查
`tail_return_invariant_residual_max`（$N/P_x/K$）和MPI请求余额；另行输出
`tail_return_representation_residual_max`（$J_x/\Pi_{xx}/\Pi_\perp$），由速度网格收敛、
能谱、密度和场A/B共同验收，不能再误作逐组精确守恒门。这只修复测试判定与错误定位，
不改变生产物理状态。

**历史处置（已关闭）**：当时单独重跑比较脚本不能解决无accepted-step问题，必须先修复R2零
可行率再重跑A/B。后续R2投影、事务和checkpoint问题已经修复，R3与R4均通过；本段不得再作为
当前阻断条件。

A：返回关闭。

```bash
rm -rf ./output/h10_r3a_none_100_to_100p3
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir <100fs完整checkpoint目录> \
  --restart-allow-return-config-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 100.3 --snapshot-times 100.3 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_r3a_none_100_to_100p3

test "$(wc -l < ./output/h10_r3a_none_100_to_100p3/vpfp_step_diagnostics.dat)" -gt 1 || {
  echo 'R3a baseline produced no accepted step'; exit 3;
}
find ./output/h10_r3a_none_100_to_100p3 -maxdepth 1 -type d \
  -name 'snapshot_*' -print -quit | grep -q . || {
  echo 'R3a baseline produced no snapshot'; exit 4;
}
find ./output/h10_r3a_none_100_to_100p3 -maxdepth 2 -type f \
  -name '_COMPLETE' -print -quit | grep -q . || {
  echo 'R3a baseline snapshot was not committed'; exit 7;
}
```

这里A、B两侧都必须带`--restart-allow-return-config-change`。100 fs checkpoint由
`tail_return_mode=none`写出；该开关只允许读取checkpoint中保存的H10返回参数，并据此重建
预期物理哈希，仍会严格拒绝网格、时间步、场边界、背景边界、Beam、碰撞、转换阈值和
flux-interface配置的任何变化。A侧仍显式使用`--tail-return-mode none`，因此该开关不会启用
返回，也不会改变A侧物理推进；它只消除“关闭状态下未参与动力学的return默认参数变化”造成的
假性restart不兼容。

旧的`vpfp-open-v2` checkpoint物理哈希只包含`tail_return_mode`，尚未包含后来增加的
`return_energy/residence/radius/tolerance`四项。读取器以四个精确零表示这些缺失字段；兼容路径
必须复现“只哈希mode”的旧布局，不能把四个合成零值追加到哈希中。已知100 fs checkpoint的
旧哈希为`10068791302029124601`；当前实现只在上述四项均为零且mode为`none`时采用旧布局，
新checkpoint仍使用完整H10物理哈希。

B：迟滞返回开启。

```bash
rm -rf ./output/h10_r3a_return_100_to_100p3
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir <100fs完整checkpoint目录> \
  --restart-allow-return-config-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
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
  --checkpoint-times 100.3 \
  --stop-time-fs 100.3 --snapshot-times 100.3 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_r3a_return_100_to_100p3

test "$(wc -l < ./output/h10_r3a_return_100_to_100p3/vpfp_step_diagnostics.dat)" -gt 1 || {
  echo 'R3a candidate produced no accepted step'; exit 5;
}
find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
  -name 'snapshot_*' -print -quit | grep -q . || {
  echo 'R3a candidate produced no snapshot'; exit 6;
}
find ./output/h10_r3a_return_100_to_100p3 -maxdepth 2 -type f \
  -name '_COMPLETE' -print -quit | grep -q . || {
  echo 'R3a candidate snapshot was not committed'; exit 8;
}
find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
  -name 'checkpoint_target100.3fs*' -print -quit | grep -q . || {
  echo 'R3a candidate produced no 100.3 fs checkpoint'; exit 9;
}
```

比较：

```bash
python3 tools/compare_tail_return_ab.py \
  --baseline ./output/h10_r3a_none_100_to_100p3 \
  --candidate ./output/h10_r3a_return_100_to_100p3 \
  --result ./output/h10_r3a_compare.result

grep -E '^(comparison_error|conservation_pass|tail_return_invariant_residual_max|tail_return_representation_residual_max|spectrum_pass|field_pass|wave_envelope_pass|particle_reduction_observed|status)=' \
  ./output/h10_r3a_compare.result
grep -q '^status=PASS$' ./output/h10_r3a_compare.result || {
  echo 'R3a failed; do not run R3b'; exit 10;
}
```

`compare_tail_return_ab.py`必须实现§16.6.4的守恒门和物理不变门，并在结果中逐项输出
`conservation_pass`、`spectrum_pass`、`field_pass`、`wave_envelope_pass`、
`particle_reduction_observed`和最终`status`；不得只比较目录hash或文件大小。

#### 17.13.6 R3b：同一预热态物理A/B与独立controller-off性能基线

**2026-08-13最终根因：原R3b启动方法无效；必须从同一个已完成返还预热的checkpoint分叉。**

首轮v2结果中，A/B均接受39步，返回侧Tail宏粒子由99,194,116降至60,145,110，
Tail内核墙钟下降42.0%；$N/P_x/K$最大相对残差为$4.61\times10^{-16}$，因此失败不是
返回事务、硬守恒或性能问题。失败项为：场相对$L_2=4.40\times10^{-3}$、场能代理相对差
$1.20\times10^{-3}$、场包络相对$L_2=2.27\times10^{-3}$、combined能谱相对
$L_2=8.43\times10^{-3}$。它们均超过§16.6给出的短程扰动门。

根因是v2投影仅把$N/P_x/K$放入非负硬约束，观测速度格点超过24个时又按索引等数量分桶，
把整桶质量放到单个代表格点；该操作可保持粒子数，却可丢失大部分$\Pi_\perp$。首轮R3b中
$J_x/\Pi_{xx}/\Pi_\perp$表示残差最大值为0.991，1 fs窗口已经足以使该压力结构误差进入
能谱和场响应。因此不得通过放宽比较阈值把首轮结果改判为通过。

生产投影已升级为v5：本rank拥有的物理cell保留完整聚合速度直方图，并在非负
$N/P_x/K$精确可行集合内选择距该直方图最近的权重；只有跨MPI接缝的固定大小通信记录才做
六矩守恒压缩。旧v4虽然改进了最终权重，却仍在所有本地cell投影前把完整直方图压缩成少量
支撑，因此真实R3b中的combined能谱几乎没有改善。旧v3仅改善先验，
但随后的普通三矩求解会覆盖该先验，因此真实R3b中表示残差没有改善，现已被v4替代。
$J_x/\Pi_{xx}/\Pi_\perp$
没有被错误提升为可能不可行的硬约束。单rank生产投影回归中，宽连续速度云的硬不变量残差
为$6.46\times10^{-14}$、表示矩残差为$6.87\times10^{-4}$；单个离格粒子的表示矩残差为
$1.36\times10^{-3}$，且bulk保持非负。none、hysteresis、projection、transaction四项回归均
通过。以上只证明局部投影修复有效；仍须按本节顺序重跑R2、R3a和R3b，最终以真实checkpoint
动态结果验收。

**2026-08-13第五轮R3b根因收敛与v5修复：** 最新原目录结果仍为`status=FAIL`，但它证明
v4确实进入生产：表示矩最大残差由0.995降至0.0838；$N/P_x/K$残差为
$2.66\times10^{-15}$，Tail内核墙钟下降42.1%。物理失败仍为场相对$L_2=4.41\times10^{-3}$、
场能相对差$1.14\times10^{-3}$、场包络$2.16\times10^{-3}$和能谱$8.43\times10^{-3}$。
逐步对齐显示3911--3917步A/B完全一致；3918步驻留门同时释放38,018,745个宏粒子，返回权重
$4.40\times10^{19}$。该事务完成瞬间combined的$N/P_x/K$仍在$10^{-15}$量级一致，但随后
场、combined动量和能量逐步分离，说明根因是完整速度分布被少量支撑替换后，在不同
PIC/Eulerian推进算子下演化分叉，而不是命令、Poisson场或事务守恒错误。

v5因此取消本地物理cell的少支撑压缩：本地投影使用完整聚合速度直方图作为先验，只在跨rank
固定通信记录中保留六矩压缩；最终仍以非负最近点投影精确闭合$N/P_x/K$。新增45支撑宽云
回归，验证全部粒子返回、非负、硬不变量残差$1.11\times10^{-16}$、表示矩残差
$1.18\times10^{-3}$，并显式输出速度直方图相对$L_1=5.82\times10^{-2}$。新的生产指纹为
`conservative_n_px_k_local_histogram_v5`。下一轮不得再使用旧二进制或旧快照混合比较；启动日志
必须出现该指纹，且A/B两项结束后才能执行比较器。

第六轮结果进一步证明，旧R3b并非正常的稳态返还测试。100 fs checkpoint由
`tail_return_mode=none`生成，所有粒子的`return_residence_steps`均为0；从该checkpoint临时开启
`hysteretic`后，约$3.8\times10^7$个既有低能宏粒子在第8步同时达到驻留门并一次性迁移到
Eulerian bulk。迁移事务当步的combined $N/P_x/K$仍闭合到$10^{-15}$，但之后PIC Tail与
Eulerian bulk由不同输运/碰撞算子推进，1 fs内场相对$L_2$达到$4.39\times10^{-3}$、combined
能谱相对$L_2$达到$8.41\times10^{-3}$。这属于“从关闭状态热切换控制器”造成的同步启动冲击，
不能用来代表从计算开始持续启用H10时自然错开的返还队列。

不得通过按运行时长线性放宽物理阈值将该结果改判为PASS，也不得在生产算子中给新转换粒子
伪造驻留年龄。正确的R3b流程是：先在R3a的hysteretic分支完成100--100.3 fs预热并保存
checkpoint；随后让R3b的none和hysteretic两侧都从这个完全相同的预热后checkpoint出发。
这样初态、粒子位置/速度、随机数、bulk分布和既有驻留年龄严格相同，A侧只停止后续返还，B侧
继续正常返还，测试不会再包含一次性冷启动迁移。

快照写出协议已增加`_COMPLETE`提交标记：全部rank完成快照和rank manifest写出后，rank 0才
原子意义上提交该标记。比较器要求标记存在、rank数有效，并核对
`manifest_rank/fields_rank/density_background_rank/energy_spectrum_rank`文件数量和非空性。
因此更新后必须重编译`fp_solver`并重跑A/B；旧快照没有完成标记，不能继续作为最终验收输入。

只有R3a为`status=PASS`才运行。上述R3a候选命令已经包含`--checkpoint-times 100.3`。
先取得唯一的预热checkpoint；若得到0个或多个目录，都必须停止，不能猜测目录：

```bash
mapfile -t R3B_CHECKPOINTS < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#R3B_CHECKPOINTS[@]}" -eq 1 || {
  echo "expected one R3b warm checkpoint, found ${#R3B_CHECKPOINTS[@]}"; exit 20;
}
R3B_START="${R3B_CHECKPOINTS[0]}"
test -s "$R3B_START/manifest.dat" || {
  echo 'R3b warm checkpoint has no manifest'; exit 21;
}
```

R3b A关闭后续返还，但从同一个预热后状态开始：

```bash
rm -rf ./output/h10_r3b_none_100p3_to_101
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$R3B_START" \
  --restart-allow-return-config-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 101 --snapshot-times 101 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_r3b_none_100p3_to_101
```

R3b B继续hysteretic返还，并严格复用相同的`R3B_START`：

```bash
rm -rf ./output/h10_r3b_return_100p3_to_101
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$R3B_START" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
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
  --stop-time-fs 101 --snapshot-times 101 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_r3b_return_100p3_to_101
```

这里仅A侧需要`--restart-allow-return-config-change`，因为它把预热checkpoint中的
`hysteretic`切换为`none`；B侧与checkpoint的return配置完全相同，故使用普通严格重启。
该兼容开关现在允许`none -> hysteretic`和`hysteretic -> none`双向A/B，但实现仍会先用checkpoint
中的H10参数重建完整物理哈希，任何非return参数变化都会拒绝，不能将它当作通用的忽略哈希开关。

两个`yhrun`都正常返回后，检查步数、失败日志和快照提交标记。任何一项失败都不得运行比较器：

```bash
BASE=./output/h10_r3b_none_100p3_to_101
CAND=./output/h10_r3b_return_100p3_to_101
test "$(wc -l < "$BASE/vpfp_step_diagnostics.dat")" -gt 1 || exit 22
test "$(wc -l < "$CAND/vpfp_step_diagnostics.dat")" -gt 1 || exit 23
test ! -s "$BASE/vpfp_failure.dat" || exit 24
test ! -s "$CAND/vpfp_failure.dat" || exit 25
test -n "$(find "$BASE" -maxdepth 2 -name _COMPLETE -print -quit)" || exit 26
test -n "$(find "$CAND" -maxdepth 2 -name _COMPLETE -print -quit)" || exit 27
```

同一起点A/B只用于物理扰动验收。由于100.3 fs预热态已经删除了约$3.8\times10^7$个Tail宏粒子，
A、B在随后0.7 fs内的平均Tail数只相差约0.81%，该对照在数学上不可能产生20%的内核加速。
性能门必须另设第三路controller-off基线：从原始100 fs、`tail_return_mode=none` checkpoint推进到
101 fs，仅将100.3--101 fs的Tail计时和粒子数送入性能门。该第三路不参与密度、场、能谱或守恒
比较，因而不会用不同初态伪造物理等价性。

```bash
rm -rf ./output/h10_r3b_perf_none_100_to_101
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir <100fs完整checkpoint目录> \
  --restart-allow-return-config-change \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7 \
  --tail-return-mode none \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --stop-time-fs 101 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_r3b_perf_none_100_to_101

test "$(wc -l < ./output/h10_r3b_perf_none_100_to_101/vpfp_step_diagnostics.dat)" -gt 1 || exit 28
test ! -s ./output/h10_r3b_perf_none_100_to_101/vpfp_failure.dat || exit 29
```

最后单独执行：

```bash
python3 tools/compare_tail_return_ab.py \
  --baseline ./output/h10_r3b_none_100p3_to_101 \
  --candidate ./output/h10_r3b_return_100p3_to_101 \
  --performance-baseline ./output/h10_r3b_perf_none_100_to_101 \
  --performance-window-fs 100.3,101 \
  --result ./output/h10_r3b_compare.result

cat ./output/h10_r3b_compare.result
grep -q '^status=PASS$' ./output/h10_r3b_compare.result || exit 30
```

验收：物理门仅使用同预热起点的A/B并继续通过；性能门使用独立controller-off基线，要求
`tail_particle_count_candidate<=tail_particle_count_baseline`且排除I/O后
`tail_kernel_wall_reduction>=0.20`。比较结果必须输出
`performance_baseline_is_separate=1`，否则三路测试未正确接入。

**最终验收结果（2026-08-13）：R3整体通过。**

R3a：

```text
comparison_valid=1
conservation_pass=1
density_pass=1
field_pass=1
spectrum_pass=1
wave_envelope_pass=1
particle_reduction_observed=1
status=PASS
```

R3b三路比较：

```text
baseline_accepted_step_count=28
candidate_accepted_step_count=28
conservation_pass=1
density_relative_l2=3.942517827962118e-05
field_relative_l2=3.180202507159027e-04
field_energy_relative_difference=5.933669286621386e-06
field_envelope_relative_l2=3.823698352005829e-05
spectrum_relative_l2=8.855296998707048e-05
tail_return_invariant_residual_max=2.182938714185434e-15
tail_return_request_residual_max=0
performance_baseline_is_separate=1
tail_particle_count_baseline=99194116
tail_particle_count_candidate=60119066
tail_kernel_wall_reduction=0.4165862789469342
performance_pass=1
status=PASS
```

R3a/R3b四个物理分支均有非空accepted-step诊断和已提交`_COMPLETE`快照，未发现末态失败日志。
独立性能基线目录保留在集群侧、未同步到本机；本机保存的`h10_r3b_compare.result`已经记录其
粒子数、计时和`performance_baseline_is_separate=1`。因此结论可登记为通过，但归档时应连同
集群侧性能基线诊断一起保存，保证性能门可复算。

#### 17.13.7 R4：参数敏感性 [已通过，复现归档]

所有参数档从R3a hysteretic分支生成的同一个100.3 fs预热checkpoint分叉，并推进到102.3 fs，
形成严格2 fs的局部参数敏感性窗口。禁止回到原始100 fs controller-off checkpoint；否则所有
驻留计数重新从0开始，会再次制造数千万粒子同步返还的冷启动冲击。

分两组，不得同时改变阈值与驻留步数：

1. 固定`--tail-return-residence-steps 8`，分别运行
   `--tail-return-energy-mev 5.25|5.5|5.75`；
2. 固定`--tail-return-energy-mev 5.5`，分别运行
   `--tail-return-residence-steps 4|8|16`。

每档都必须复用§17.13.5的B命令全部非H10物理参数，仅替换上述单一H10参数、
`--stop-time-fs 102.3 --snapshot-times 102.3`和独立`--output-dir`。$K_{in}=5.5$ MeV、$N_{res}=8$
是两组共享的中央档，只运行一次；因此共五个不同参数作业，不得为了补齐矩阵重复运行中央档。

该测试是从同一已接受物理状态开始的**局部参数分叉**。降低阈值后，不再满足新阈值的粒子在
首个接受步正确清零驻留计数；提高阈值后，新纳入能带从0开始计数。不得修改checkpoint中的逐粒子
驻留年龄，也不得为不同参数伪造历史。

先重新编译生产程序、R4参数切换回归和比较器，并通过前置测试：

```bash
cmake --build build -j4 --target \
  fp_solver tail_bulk_return_parameter_switch_test

rm -f ./output/h10_r4_parameter_switch_precheck.result
yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_parameter_switch_test \
  --case all \
  --result ./output/h10_r4_parameter_switch_precheck.result

cat ./output/h10_r4_parameter_switch_precheck.result
grep -q '^status=PASS$' ./output/h10_r4_parameter_switch_precheck.result || exit 30
python3 -m py_compile tools/compare_tail_return_sensitivity.py || exit 31
```

前置结果必须同时满足：

```text
lower_threshold_resets_ineligible_age=1
raised_threshold_starts_age_at_one=1
residence_switch_uses_checkpointed_age=1
status=PASS
```

先定义与R3完全相同的物理参数数组：

```bash
mapfile -t R4_CHECKPOINTS < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#R4_CHECKPOINTS[@]}" -eq 1 || {
  echo "expected one R4 warm checkpoint, found ${#R4_CHECKPOINTS[@]}"; exit 32;
}
CHECKPOINT_100P3="${R4_CHECKPOINTS[0]}"
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 33

COMMON_H10=(
  --restart-dir "$CHECKPOINT_100P3"
  --restart-allow-return-config-change
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0
  --beam-enabled 1 --collision-model moment-closure
  --bulk-collision-integrator chang-cooper-flux
  --collision-interface-mode exporting-absorbing
  --background-tail-mode pic
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface
  --tail-flux-quadrature-order 4 --tail-flux-max-supports 7
  --tail-return-mode hysteretic
  --tail-return-max-stencil-radius 3
  --tail-return-moment-tolerance 1e-12
  --tail-collision-kernel coulomb-nanbu-perez
  --tail-collision-weight-mode virtual-split
  --tail-collision-max-substeps 1024
  --tail-collision-max-particle-growth 0
  --tail-population-control-interval 0
  --stop-time-fs 102.3 --snapshot-times 102.3
  --diagnostic-level 1 --diagnostic-interval 1
)
```

$K_{in}=5.25$ MeV，$N_{res}=8$：

```bash
rm -rf ./output/h10_r4_kin_525
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver "${COMMON_H10[@]}" \
  --tail-return-energy-mev 5.25 --tail-return-residence-steps 8 \
  --output-dir ./output/h10_r4_kin_525
```

共享中央档$K_{in}=5.5$ MeV，$N_{res}=8$：

```bash
rm -rf ./output/h10_r4_kin550_res8
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver "${COMMON_H10[@]}" \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 8 \
  --output-dir ./output/h10_r4_kin550_res8
```

$K_{in}=5.75$ MeV，$N_{res}=8$：

```bash
rm -rf ./output/h10_r4_kin_575
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver "${COMMON_H10[@]}" \
  --tail-return-energy-mev 5.75 --tail-return-residence-steps 8 \
  --output-dir ./output/h10_r4_kin_575
```

$K_{in}=5.5$ MeV，$N_{res}=4$：

```bash
rm -rf ./output/h10_r4_res_4
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver "${COMMON_H10[@]}" \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 4 \
  --output-dir ./output/h10_r4_res_4
```

$K_{in}=5.5$ MeV，$N_{res}=16$：

```bash
rm -rf ./output/h10_r4_res_16
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver "${COMMON_H10[@]}" \
  --tail-return-energy-mev 5.5 --tail-return-residence-steps 16 \
  --output-dir ./output/h10_r4_res_16
```

五档执行后，运行：

```bash
for RUN in \
  ./output/h10_r4_kin_525 \
  ./output/h10_r4_kin550_res8 \
  ./output/h10_r4_kin_575 \
  ./output/h10_r4_res_4 \
  ./output/h10_r4_res_16; do
  test "$(wc -l < "$RUN/vpfp_step_diagnostics.dat")" -gt 1 || exit 34
  test ! -s "$RUN/vpfp_failure.dat" || exit 35
  test -n "$(find "$RUN" -maxdepth 2 -name _COMPLETE -print -quit)" || exit 36
done

python3 tools/compare_tail_return_sensitivity.py \
  --threshold-525 ./output/h10_r4_kin_525 \
  --threshold-550 ./output/h10_r4_kin550_res8 \
  --threshold-575 ./output/h10_r4_kin_575 \
  --residence-4 ./output/h10_r4_res_4 \
  --residence-8 ./output/h10_r4_kin550_res8 \
  --residence-16 ./output/h10_r4_res_16 \
  --result ./output/h10_r4_sensitivity.result

cat ./output/h10_r4_sensitivity.result
grep -q '^status=PASS$' ./output/h10_r4_sensitivity.result || exit 37
```

R4比较器采用适用于相邻控制参数档的组合波形门：场相对$L_2\le5\times10^{-3}$、场相关系数
$\ge0.99998$、$\max|\Delta E|/\max|E|\le10^{-2}$；场能代理相对差$\le10^{-3}$、密度相对
$L_2\le10^{-3}$、combined谱相对$L_2\le5\times10^{-3}$。不同参数档不是bitwise replay，禁止
继续使用旧的$10^{-9}$/$10^{-10}$近逐点相等门。所有档还必须具有相同accepted时间窗口、完整
`_COMPLETE`快照、无未解决失败、返回$N/P_x/K$残差$\le10^{-12}$和MPI请求残差$\le10^{-13}$。
比较结果还必须逐档报告末态Tail粒子数、累计返回粒子/物理数、返回包的
$\Pi_{xx}/\Pi_\perp$、平均总wall time和Tail内核wall time。压力矩和耗时用于选择稳定平台，
不要求不同控制参数逐点相等。

最终`h10_r4_sensitivity.result`为`status=PASS`且`parameter_plateau=1`。准生产H10配置登记为
`--tail-return-mode hysteretic --tail-return-energy-mev 5.5 --tail-return-residence-steps 8`；进入新的
长期生产前仍需完成连续边界能量账验证。

最终结果摘要：

| 扫描对 | 场相对$L_2$ | 场相关系数 | 最大局部差/峰值 | 场能相对差 | 密度相对$L_2$ | combined谱相对$L_2$ |
|---|---:|---:|---:|---:|---:|---:|
| $K_{in}:5.25\to5.5$ MeV | $2.339\times10^{-3}$ | 0.99999727 | $3.999\times10^{-3}$ | $2.288\times10^{-4}$ | $2.338\times10^{-4}$ | $5.010\times10^{-4}$ |
| $K_{in}:5.5\to5.75$ MeV | $3.357\times10^{-3}$ | 0.99999437 | $8.078\times10^{-3}$ | $4.568\times10^{-5}$ | $3.809\times10^{-4}$ | $1.102\times10^{-3}$ |
| $N_{res}:4\to8$ | $2.417\times10^{-3}$ | 0.99999708 | $5.449\times10^{-3}$ | $8.701\times10^{-5}$ | $2.271\times10^{-4}$ | $6.835\times10^{-5}$ |
| $N_{res}:8\to16$ | $2.281\times10^{-3}$ | 0.99999740 | $4.283\times10^{-3}$ | $1.566\times10^{-4}$ | $2.398\times10^{-4}$ | $1.154\times10^{-4}$ |

五档均满足：相同accepted窗口、完整`_COMPLETE`快照、无未解决失败、返回不变量最大残差
$4.663\times10^{-15}$、MPI请求残差为0。末态Tail宏粒子数随$K_{in}=5.25/5.5/5.75$ MeV
分别为63,069,952、61,141,128、58,263,961；随$N_{res}=4/8/16$分别为61,012,203、
61,141,128、61,401,317，变化方向符合“阈值越高或驻留越短，返回越多”的控制逻辑。

初次`FAIL`不是生产物理失败，而是比较器误用R3同参数A/B的场相对$L_2\le2\times10^{-3}$。
修复后不只是放宽到$5\times10^{-3}$，还加入相关系数与局部最大差两项形状门；其余物理和守恒
阈值未放宽。因此当前`PASS`同时满足宏观平台与局部波形结构约束。

### 17.14 H10中央档连续准生产验证 [已完成：资源通过/长期物理未通过]

本阶段不重跑R4，也不从原始100 fs `tail_return_mode=none` checkpoint直接启动。必须从R3a已经
预热驻留年龄的100.3 fs hysteretic checkpoint继续，以$K_{in}=5.5$ MeV、$N_{res}=8$推进到
120 fs。该窗口用于验证长期返回、边界能量账和资源增长，不得把100.3 fs以前的历史能量伪装成
本次累计量。

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target fp_solver

mapfile -t CHECKPOINTS_100P3 < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#CHECKPOINTS_100P3[@]}" -eq 1 || {
  echo "expected one hysteretic 100.3 fs checkpoint, found ${#CHECKPOINTS_100P3[@]}";
  exit 38;
}
CHECKPOINT_100P3="${CHECKPOINTS_100P3[0]}"
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 39

rm -rf ./output/h10_quasiprod_100p3_to_120
yhrun -N 5 -n 80 -c 4 --cpu-bind=cores \
  stdbuf -oL -eL ./build/fp_solver \
  --restart-dir "$CHECKPOINT_100P3" \
  --field-boundary dirichlet-phi --phi-left 0 --phi-right 0 \
  --beam-enabled 1 --collision-model moment-closure \
  --bulk-collision-integrator chang-cooper-flux \
  --collision-interface-mode exporting-absorbing \
  --background-tail-mode pic \
  --tail-convert-energy-mev 6.0 --tail-conversion-mode flux-interface \
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
  --checkpoint-times 110,120 \
  --stop-time-fs 120 --snapshot-times 110,120 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir ./output/h10_quasiprod_100p3_to_120
```

运行后先执行最低成本完整性门：

```bash
RUN=./output/h10_quasiprod_100p3_to_120
test -s "$RUN/vpfp_step_diagnostics.dat" || exit 40
test ! -s "$RUN/vpfp_failure.dat" || exit 41
test -n "$(find "$RUN" -maxdepth 2 -name _COMPLETE -print -quit)" || exit 42
tail -n 1 "$RUN/vpfp_step_diagnostics.dat"
```

验收时必须同时报告：

1. 100.3--120 fs全部时间步accepted，NaN/Inf、未解决failure、split和rollback为0；
2. 每步返回$N/P_x/K$残差不超过$10^{-12}$，MPI请求残差不超过$10^{-13}$；
3. 分别累计背景左右边界入/出能量、Beam入/出能量、`collision_reservoir`以及域内
   $U_E+K_e+K_b$变化，明确能量账只覆盖100.3--120 fs；
4. 报告Tail宏粒子总数、单rank最大值、平均/末段wall time和Tail内核时间，并与旧
   controller-off 100--120 fs窗口比较；
5. 比较110/120 fs的combined能谱、密度、电场、平均动能和动量分布，确认迟滞返回未在
   5--6.5 MeV接口附近制造人工谷或宏观波形突变。

实测770步全部接受，split/rollback/failure为0，返回硬守恒通过；末态Tail宏粒子减少44.0%，
平均wall time减少41.6%，资源门通过。但相对controller-off参考，120 fs电场相对$L_2$为
16.29%、combined能谱为5.18%、$u_\parallel$分布为2.44%、bulk $u_\perp$分布为4.38%，
combined背景动能高3.13%。因此五项并未全部通过，禁止直接从0 fs生成H10最终生产结果。

旧诊断重建的100.3--120 fs能量余额约$4.13\times10^7\ \mathrm{J/m^2}$，但缺少固定电势
边界功和accepted-step完整恒等式，当前只能判为`INCOMPLETE`，不能判成数值生能。下一步见
§17.15；不得回退static-cell、开启Population Controller或放宽守恒阈值。

### 17.15 H10能量恒等式与长期表示误差短A/B [短A/B完成，局部修复待集群回归]

#### 17.15.1 已实现诊断

本阶段只增加只读诊断，不修改物理推进、返回投影、碰撞、Poisson、事务接受或checkpoint状态：

1. `OpenElectrostaticSolver::boundary_energy_work()`按接受步前后端面场计算固定电势电极功
   $W_\phi=\phi_L\Delta(\varepsilon_0E_L)-\phi_R\Delta(\varepsilon_0E_R)$。当前
   `DIRICHLET_PHI(0,0)`严格返回0且不增加MPI归约；
2. `VpfpStepLedger`记录步前/步后$U_E+K_{bulk}+K_{tail}+K_{beam}$、背景和Beam开放边界净能流、
   Tail物理出流、碰撞储库、边界功、总源以及单步余额；
3. 能量账采用
   $$
   \Delta U_{domain}=Q_{bkg,bnd}+Q_{beam,bnd}-Q_{tail,out}-Q_{collision,res}+W_\phi+R_E.
   $$
   `energy_balance_residual`即$R_E$，仅诊断，不参与接受；
4. 每次Tail返回输出Eulerian新增减PIC删除的$N/P_x/J_x/K/\Pi_{xx}/\Pi_\perp$有符号差。
   $N/P_x/K$仍是硬事务不变量，其他矩只读累计，禁止据此事后修正状态；
5. `tools/analyze_h10_quasiproduction.py`分别给出数值、守恒、性能、能量和宏观物理Gate，
   防止“性能PASS”被误写成“H10整体PASS”。

#### 17.15.2 构建与低成本前置验收

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver
module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

cmake --build build -j4 --target \
  open_electrostatic_solver_test tail_bulk_return_projection_test fp_solver

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/open_electrostatic_solver_test --case all \
  --field-boundary dirichlet-phi \
  --result ./output/h10_energy_ledger_poisson.result

yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
  ./build/tail_bulk_return_projection_test --case all \
  --result ./output/h10_signed_return_moments.result

grep -q '^status=PASS$' ./output/h10_energy_ledger_poisson.result || exit 43
grep -q '^status=PASS$' ./output/h10_signed_return_moments.result || exit 44
python3 -m py_compile tools/analyze_h10_quasiproduction.py || exit 45
```

#### 17.15.3 同checkpoint短A/B

两档必须从同一个预热100.3 fs hysteretic checkpoint出发并使用相同MPI分解。A档临时关闭
返回需要显式允许只改变return配置；B档沿用checkpoint中的迟滞状态。只推进到105 fs，不创建
新的0--120 fs结果。

以下命令故意不使用`COMMON_H10_105`等Shell数组。每一档都能独立放入新的Slurm作业脚本；
不得只复制末尾的`yhrun`片段而遗漏checkpoint解析和完整物理配置。

**A档：关闭返回。** 输出目录必须不存在；若需要重跑，应改用新的目录名或人工确认后清理旧目录。

```bash
mapfile -t CHECKPOINTS_100P3 < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#CHECKPOINTS_100P3[@]}" -eq 1 || exit 46
CHECKPOINT_100P3="${CHECKPOINTS_100P3[0]}"
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 47
OUTPUT_DIR=./output/h10_energy_ab_none_100p3_to_105
test ! -e "$OUTPUT_DIR" || { echo "output directory already exists: $OUTPUT_DIR" >&2; exit 48; }

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
  --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 \
  --tail-collision-kernel coulomb-nanbu-perez \
  --tail-collision-weight-mode virtual-split \
  --tail-collision-max-substeps 1024 \
  --tail-collision-max-particle-growth 0 \
  --tail-population-control-interval 0 \
  --tail-return-mode none \
  --stop-time-fs 105 --snapshot-times 105 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir "$OUTPUT_DIR"
```

**B档：开启迟滞返回。** 该命令再次独立解析同一checkpoint，且显式传入
`--background-tail-mode pic`，不会再触发“hysteretic requires pic”配置错误。

```bash
mapfile -t CHECKPOINTS_100P3 < <(
  find ./output/h10_r3a_return_100_to_100p3 -maxdepth 1 -type d \
    -name 'checkpoint_target100.3fs*' -print
)
test "${#CHECKPOINTS_100P3[@]}" -eq 1 || exit 46
CHECKPOINT_100P3="${CHECKPOINTS_100P3[0]}"
test -s "$CHECKPOINT_100P3/manifest.dat" || exit 47
OUTPUT_DIR=./output/h10_energy_ab_return_v8_100p3_to_105
test ! -e "$OUTPUT_DIR" || { echo "output directory already exists: $OUTPUT_DIR" >&2; exit 48; }

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
  --tail-flux-quadrature-order 4 \
  --tail-flux-max-supports 7 \
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
  --stop-time-fs 105 --snapshot-times 105 \
  --diagnostic-level 1 --diagnostic-interval 1 \
  --output-dir "$OUTPUT_DIR"
```

#### 17.15.4 分析命令与判定

```bash
python3 tools/analyze_h10_quasiproduction.py \
  --return-run ./output/h10_energy_ab_return_v8_100p3_to_105 \
  --baseline-run ./output/h10_energy_ab_none_100p3_to_105 \
  --result ./output/h10_energy_ab_v8_100p3_to_105.result

cat ./output/h10_energy_ab_v8_100p3_to_105.result
```

必须分别解释五类结果，不能只看末行：

- `numerical_gate_pass=1`：两档完整推进，无failure/split/rollback；
- `conservation_gate_pass=1`：返回$N/P_x/K\le10^{-12}$且MPI请求残差$\le10^{-13}$；
- `energy_ledger_complete=1`：新列真实存在。`absolute_energy_gate_pass`检查完整求解器累计余额是否
  $\le10^{-3}$；`energy_gate_pass`在有同窗关闭档时只检查H10引入的增量余额是否$\le10^{-3}$。
  若A/B均同号同量级，H10可以通过，但`production_ready`必须保持0，直到公共能量账或公共算子修复；
- `return_*_signed_difference_total`：检查$J_x/\Pi_{xx}/\Pi_\perp$偏差是否具有稳定符号并与
  combined动能漂移同步。它们是根因证据，不是新的强制零残差门；
- `physical_gate_pass=1`的短程目标为：场相对$L_2\le2\%$、密度$\le1\%$、combined谱
  $\le2\%$、两类动量边缘分布$\le2\%$。

不足10 fs的窗口中，`performance_gate_evaluated=0`且性能门不参与总判定。此时仍输出粒子数和
wall-time变化供参考，但不能要求已经达到长期20%降幅，也不能据此否定返回算法。能量门不因
窗口较短而放宽。

只有数值、守恒和H10增量能量Gate通过，且宏观差异没有随累计返回单调增长，才允许把同一A/B扩展到
111 fs。若短程已复现单调漂移，下一项不是放宽阈值，而是设计保持受控高阶矩或减少算子切换
突变的局部返回表示；在获得独立收敛测试前不得接入生产。

注意：`status=PASS`表示H10算子级A/B通过，不等于整个求解器生产就绪。只有
`production_ready=1`才表示公共绝对能量门也通过。该拆分禁止用公共离散能量缺陷反复误判H10，
同时也禁止借A/B相消掩盖完整程序的能量问题。

#### 17.15.5 已实施的H10局部表示修复

短A/B证明能量余额是关闭/开启返回共有的问题，但同时暴露了H10自身的压力矩表示偏差。修复
严格限定在Tail到bulk的表示转换，不修改Vlasov、Poisson、Beam、碰撞或时间分裂：

1. 本地返回单元不再直接把观测速度直方图作为最终投影先验。现在先保留完整直方图作为
   reference prior，再调用既有`improve_representation_prior()`，在非负约束下最小化
   $J_x/\Pi_{xx}/\Pi_\perp$表示误差；
2. 最终仍由`project_invariants_near_prior()`或三不变量几何回退严格投影$N/P_x/K$。高阶矩
   优化不得替代三个事务不变量；
3. 对最终候选计算$J_x/\Pi_{xx}/\Pi_\perp$最大相对表示误差。`v7`曾固定要求不超过
   $5\times10^{-3}$；该要求只适用于速度恰好落在网格中心附近的理想组，不能直接用于真实
   连续速度PIC云；
4. 新增`tail_return_projection_representation_incompatible_cells`，把这种物理表示不兼容与
   输入无效、支撑不足、$N/P_x/K$不可行明确分开，并纳入MPI全局统计；
5. 该门控是“是否允许表示转换”的保真门，不是能量补丁。它可能降低粒子回收率，但不会修改
   已接受粒子的速度或对Eulerian分布做事后缩放。

第一版修复的`v6`仍错误使用无量纲常数1作为$J_x/\Pi$矩的尺度。100.3--105 fs回归中
`return_representation_incompatible_cells_total=0`且$\Pi_\perp$累计偏差仍为$-13.6129$，
证明该门在生产量级下没有生效。根因是cell-volume积分矩经常小于1，原
`max(1,|target|,target_l1)`把0.5%相对门退化为绝对误差0.005，并在优化目标中系统性低估
$J_x/\Pi$行。

现已改为按每个矩自身的$\max(|target|,target_{L1})$无量纲化；仅当该尺度为零时使用与候选
列幅值相关的舍入级fallback。全局返回残差也取消无量纲常数1的分母下限。新实现指纹为
`conservative_n_px_k_dimensionless_representation_gate_v7`。集群结果若仍打印`v5`或`v6`
指纹，必须判为旧二进制，不能用于本节验收。

`v7`的100.3--105 fs实算进一步给出：`numerical/conservation/physical`三门均通过，场、密度、
combined谱相对差分别仅为$1.22\%$、$0.123\%$和$0.0365\%$，$\Pi_\perp$累计有符号偏差也从
旧实现约$10.8\%$降到$5.61\times10^{-4}$。但184步内531246个单元贡献被判为表示不兼容，
最后一步1644个返回组只提交6组，Tail粒子仅下降$1.91\%$且wall time反而增加$6.54\%$。
这不是碰撞或返回事务错误，而是固定$0.5\%$门低于当前速度网格的固有离散表示误差。

背景Eulerian状态及生产矩均使用速度单元中心核，而Tail使用连续PIC速度核。历史真实阈值审计
已经测得最近单元表示的$\Pi_\perp$固有相对差约$3.3335\%$；因此固定$0.5\%$会使绝大多数
真实组数学上不可通过。`v8`改为网格感知门：

1. 对每个候选组先用其真实PIC速度构造最近速度单元直方图，并计算该组原生最近单元映射的
   $J_x/\Pi_{xx}/\Pi_\perp$最大相对误差$R_{native}$；
2. 最终投影候选仍严格守恒$N/P_x/K$，其允许误差为
   $\min(5\%,\max(0.5\%,1.05R_{native}))$；
3. 该门不允许候选比本组原生Eulerian离散表示更差，并以$5\%$作为硬上限；若原生误差本身
   超过硬上限，只有优化后降到$5\%$以内才允许返回；
4. 禁止直接把门改成无条件$3.5\%$或$5\%$。局部门必须来自该组实际网格量化误差，否则会
   重新引入旧实现具有稳定符号的压力偏差；
5. 新实现指纹为`conservative_n_px_k_grid_aware_representation_gate_v8`。

`v8`本地六项回归（none、hysteresis、parameter-switch、projection、transaction、MPI单rank
协议路径）均通过。真实MPI和生产分布仍必须在集群重跑以下门：

```bash
cmake --build build -j4 --target \
  tail_bulk_return_none_regression_test \
  tail_bulk_return_hysteresis_test \
  tail_bulk_return_parameter_switch_test \
  tail_bulk_return_projection_test \
  tail_bulk_return_transaction_test \
  tail_bulk_return_mpi_test fp_solver

for TEST in \
  tail_bulk_return_none_regression_test \
  tail_bulk_return_hysteresis_test \
  tail_bulk_return_parameter_switch_test \
  tail_bulk_return_projection_test \
  tail_bulk_return_transaction_test; do
  yhrun -N 1 -n 1 -c 4 --cpu-bind=cores \
    "./build/$TEST" --case all \
    --result "./output/${TEST}_reprfix_v8.result" || exit 51
done

for NP in 1 2 5; do
  yhrun -n "$NP" --cpu-bind=cores \
    ./build/tail_bulk_return_mpi_test --case all \
    --result "./output/tail_bulk_return_mpi_reprfix_v8_n${NP}.result" || exit 52
done
```

随后只需重跑§17.15.3的返回档到105 fs；关闭档未改，可以复用。新返回档必须满足：

- `numerical_gate_pass=1`、`conservation_gate_pass=1`、`physical_gate_pass=1`；
- 新表示不兼容计数显著低于`v7`的531246，最后100步提交组比例不再接近零；
- $N/P_x/K$差仍处于$10^{-12}$门内；
- $\Pi_\perp$有符号累计偏差相对旧值$-13.61$显著下降，并且相对返回量不得超过$5\%$；
- 不以4.67 fs短窗的粒子下降率或wall time作为失败原因。

`v8`回归必须使用新输出目录，不能续写本节已有的`v6/v7`结果。执行§17.15.3的B档完整独立命令，
仅将其中输出目录设置为：

```bash
OUTPUT_DIR=./output/h10_energy_ab_return_v8_100p3_to_105
```

关闭档直接复用，分析命令为：

```bash
python3 tools/analyze_h10_quasiproduction.py \
  --return-run ./output/h10_energy_ab_return_v8_100p3_to_105 \
  --baseline-run ./output/h10_energy_ab_none_100p3_to_105 \
  --result ./output/h10_energy_ab_v8_100p3_to_105.result

cat ./output/h10_energy_ab_v8_100p3_to_105.result
```

启动日志必须包含
`tail_return_projection_schema=conservative_n_px_k_grid_aware_representation_gate_v8`；否则停止验收。

#### 17.15.6 `v8`最新100.3--105 fs A/B验收结论

最新结果文件`output/h10_energy_ab_100p3_to_105.result`显示`status=PASS`。该PASS是H10算子级
A/B通过，不是完整求解器生产门全部通过。定量结论如下：

| 项目 | 最新`v8`结果 | 判定 |
|---|---:|---|
| 接受步数/物理窗口 | 184步，4.674 fs | 短窗有效 |
| `numerical_gate_pass` | 1 | 通过，无failure/split/rollback |
| `conservation_gate_pass` | 1 | 通过 |
| $N/P_x/K$最大相对残差 | $8.41\times10^{-15}$ | 显著优于$10^{-12}$门 |
| 场相对$L_2$ | $1.424\%$ | 通过$2\%$门 |
| 背景密度相对$L_2$ | $0.153\%$ | 通过$1\%$门 |
| combined能谱相对$L_2$ | $0.105\%$ | 通过$2\%$门 |
| $u_\parallel$边缘分布相对$L_2$ | $0.0771\%$ | 通过$2\%$门 |
| bulk $u_\perp$边缘分布相对$L_2$ | $0.0831\%$ | 通过$2\%$门 |
| H10增量能量相对残差 | $7.78\times10^{-4}$ | 通过$10^{-3}$门 |
| 绝对累计能量相对余额 | $4.624\%$ | 未通过，公共问题 |
| `production_ready` | 0 | 尚不可标记为生产就绪 |
| $\Pi_\perp$返回有符号差/真实移除量 | $0.9236\%$ | 通过局部表示硬门 |
| Tail最终粒子数下降 | $4.506\%$ | 仅作短窗信息 |
| 平均wall-time变化 | 慢$3.09\%$ | 仅作短窗信息 |

`v8`相对`v7`的直接改善为：Tail粒子下降率由$1.91\%$提高到$4.51\%$，移除粒子数由
1310724提高到3128649，最后一步提交组由6提高到24；同时宏观量仍全部处于短程门内。
这证明网格感知门解决了“固定0.5%导致几乎完全禁止返回”的主要问题，且没有破坏$N/P_x/K$
事务守恒。

但当前结果也明确表明H10资源目标尚未验收：184步累计尝试277945组，只提交9495组，提交率
约$3.42\%$；仍有490130个单元贡献被判为表示不兼容。该计数较`v7`的531246只下降约
$7.7\%$。因此不得把本次PASS解释为Population Controller已经有效，也不得继续仅靠放宽门限
追求粒子下降。固定$5\%$硬上限保持不变。

以下H10性能判定**暂缓执行**。完整求解器仍有约$4.6\%$的公共绝对能量余额，物理闭合的
优先级高于性能。只有§17.15.7定位并修复公共能量缺口，且短回归确认绝对能量门通过后，才从
同一100.3 fs checkpoint分别运行return-none和`v8` hysteretic到111 fs。输出目录预留为：

```bash
./output/h10_energy_ab_none_100p3_to_111
./output/h10_energy_ab_return_v8_100p3_to_111
```

两档命令沿用§17.15.2--17.15.3的完整参数，只将`--stop-time-fs`和`--snapshot-times`均改为
`111`，并替换输出目录。分析命令为：

```bash
python3 tools/analyze_h10_quasiproduction.py \
  --return-run ./output/h10_energy_ab_return_v8_100p3_to_111 \
  --baseline-run ./output/h10_energy_ab_none_100p3_to_111 \
  --result ./output/h10_energy_ab_v8_100p3_to_111.result

cat ./output/h10_energy_ab_v8_100p3_to_111.result
```

恢复该测试后，只有长窗中`performance_gate_evaluated=1`且`performance_gate_pass=1`，才能确认H10达成资源
目标。若物理与守恒门继续通过但性能门失败，下一步应优化返回分组/投影成本或调整控制器策略，
不能继续扩大表示误差硬上限。

#### 17.15.7 公共能量余额阶段审计：已完成第一轮定位

阶段审计、MPI源项所有权修复和100步none/return A/B已经完成。最新可信结果位于：

```text
output/vpfp_stage_energy_audit_none_100p3_100steps
output/vpfp_stage_energy_audit_return_100p3_100steps
output/vpfp_stage_energy_audit_none_100p3_100steps.result
output/vpfp_stage_energy_audit_return_100p3_100steps.result
```

审计实现曾有两个纯诊断错误，现均已修复：

1. 已在`ConservativePpmRemap`和`CylindricalFpCollision`内部全局归约的背景边界能流与碰撞
   reservoir被阶段审计再次按80 ranks求和；
2. `accepted_n`和`collision_half1`错误读取尚未预测的`beam_work_`，使接受态$K_{beam}$记为0。

修复后，100步A/B均满足：阶段序列完整、接受步数正确、无split/failure、Beam接受态正确，且
$x$输运、两个碰撞半步和H10返回的累计闭合残差均在约$10^{-5}\ \mathrm{J/m^2}$以内。具体为：

| 累计量 | return-none | hysteretic return | 判定 |
|---|---:|---:|---|
| $R_{x_1}$ | $6.26\times10^{-6}$ | $-1.32\times10^{-5}$ | 浮点误差 |
| $R_{x_2}$ | $-1.27\times10^{-5}$ | $-4.85\times10^{-8}$ | 浮点误差 |
| $R_{C_1}$ | $6.01\times10^{-6}$ | $8.36\times10^{-6}$ | 浮点误差 |
| $R_{C_2}$ | $-2.94\times10^{-6}$ | $8.60\times10^{-6}$ | 浮点误差 |
| $R_{H10}$ | $-4.77\times10^{-7}$ | $2.62\times10^{-6}$ | 浮点误差 |
| 完整步累计余额 | $2.9855\times10^6$ | $3.0415\times10^6$ | 仍未闭合 |

两个`.result`当前仍写`status=FAIL`，但最大阶段望远镜绝对误差仅
$1.46\times10^{-6}\ \mathrm{J/m^2}$，与正式完整步余额的最大绝对差仅
$1.91\times10^{-6}\ \mathrm{J/m^2}$。失败来自分析器固定$10^{-12}$相对门在强相消求和下过严，
不是阶段账缺失。分析器必须改成“绝对误差+机器精度缩放”的结构完整性门；该调整不得改变物理
能量门。

100步公共余额已经分解为：

| 组合 | return-none | hysteretic return | 占完整余额 |
|---|---:|---:|---:|
| midpoint/final Poisson净场能项 | $2.0871\times10^6$ | $2.1178\times10^6$ | 约70% |
| bulk/Tail/Beam受力与force后转换组合 | $8.9836\times10^5$ | $9.2369\times10^5$ | 约30% |

其中none档受力组合的分量为：bulk $+7.3891\times10^6$、Tail
$+5.0323\times10^5$、Beam $-8.2907\times10^6$、force后转换
$+1.2967\times10^6\ \mathrm{J/m^2}$。这证明主要问题不是$x$ remap、碰撞或H10，而是
Poisson约束场与三类粒子受力功是否属于同一个离散时间层和离散功系统。

下一步的详细代码改动、只读功审计、$dt/dt/2$判别、验收标准和集群命令已经迁移到独立文档：

```text
docs/VPFP公共能量残差_离散功配对审计与修复实施方案.md
```

在该独立文档完成“审计可信门”和“根因判别门”以前：

- 不得修改Poisson边界条件、H10、碰撞、$x$ remap或转换阈值；
- 不得加入全局能量补丁、粒子速度缩放或场能强制修正；
- 不得仅为`status=PASS`放宽物理能量余额门；
- §17.15.6的111 fs H10性能A/B继续暂停。

---

## 18. 性能和内存要求

### 18.1 设计目标

混合方案的目的之一是避免为稀疏高能尾部在每个 $x$ 单元分配巨大 Eulerian 速度网格。它不自动保证更快，必须控制粒子数量、沉积开销和工作数组。

### 18.2 必须实施的低风险优化

- tail 粒子使用结构化连续数组或 AoSoA，避免每粒子堆分配；
- 所有 density/current/迁移 buffer 预分配并复用；
- 转换掩码、能量和矩权重预计算；
- conversion 扫描与 forbidden-mass/矩统计融合；
- tail push、deposit 和 conversion 的本地循环 OpenMP 并行；
- 每个接受步合并 MPI 全局归约；
- 详细能谱只在快照输出时计算；
- 粒子迁移按邻居打包，不做全局 gather；
- checkpoint 每 rank 独立写，manifest 仅 rank 0 写。
- 碰撞 cell bin 使用 counting-sort 风格的 `cell_count/cell_offset/particle_index` 工作区；只重排索引，不重排粒子主数组；
- counter-based pair permutation 直接作用于 cell index slice，禁止为每个 cell 单独分配 vector；
- Coulomb logarithm 和只依赖 cell combined moments 的碰撞系数按 collision substep/cell 缓存，不对每个 pair 重算。

不要在转换热路径中调用通用 dense NNLS。候选矩阵维度固定后，应使用小规模、预分配、确定性 active-set 求解器；常见成功路径不得进行堆分配。参考 quartet 只在测试和压缩 fallback 中完整物化，正常路径可直接从 cell-mass 累计候选矩阵和权重。

### 18.3 粒子数预算

H4 必须报告：

- 3、12、25、40 fs 的全局 tail 粒子数；
- 每步新增粒子率；
- 每 rank 最大粒子数和负载不平衡；
- tail 模块占 wall time 的比例；
- 按当前增长率外推到 120 fs 的内存。

若外推粒子数不可接受，优先改进守恒聚合和受控 splitting/merging，不得重新扩大整个 Eulerian 速度域掩盖问题。

**120 fs实际记录（2026-08-12）**：40 fs Gate B曾外推全局约
$4.83\times10^7$、最大rank约$5.55\times10^6$个Tail宏粒子和约4.44 h wall time。
实际值分别为$1.4317\times10^8$、$5.30\times10^6$和5.99 h。最大rank估算尚可，
但全局粒子数低估约2.97倍，说明增长不是简单线性且Tail空间负载逐渐铺开。后续资源模型
必须至少使用分段增长率，并显式区分全局粒子数、最大rank粒子数和低能Tail比例。

单步wall time由40--50 fs的2.54 s增至110--120 fs的11.65 s。主要热点是collision、
conversion和Tail push；MPI collective与输出诊断不是主因。H10长期运行已把平均wall time
降低41.6%，证明表示返回能解决主要资源增长；但其长期物理门未通过。禁止因性能收益而忽略
§17.15的能量和高阶矩审计，也禁止仅通过放宽碰撞子步或忽略转换粒子继续提速。

---

## 19. 阶段验收判据

### 19.1 单元测试级

- 浮点有限；
- $N/P_x/K$ 转换误差接近确定性求和精度；
- Poisson 转换前后密度差接近沉积舍入误差；
- MPI rank 数变化不改变全局物理矩；
- pusher 达到预期时间阶。

阈值必须基于量的尺度和求和误差设置，禁止统一写成任意的 `1e-6`。

### 19.2 短生产级

- bulk 最外速度面质量不持续增长；
- combined 密度、动量、动能无转换台阶；
- 阈值附近能谱无明显人工峰谷；
- 转换时场不出现单步尖脉冲；
- 3 至 12 fs 原有合理波形不因 PIC 噪声破坏；
- 阈值和 bin 加密后宏观量收敛；
- tail 粒子数量和 wall time 可外推。
- combined 连续性中的 conversion source 逐 cell 抵消；
- population-controller off/on 不改变阈值附近的电流和压力矩。
- 场边界与域长敏感性不污染预先定义的核心研究区。

### 19.3 碰撞级

- Maxwellian 平衡保持；
- moment-closure 离散平衡温度与网格矩自洽：无 Beam 40fs 长跑能量漂移
  <1e-9（修复前 −39%），碰撞储库无系统符号翻转；
- 漂移/扩散矩与输入系数一致；
- bulk 与 tail 的碰撞系数同源；
- 生产 pair registry 显式开启并分别记账 $C_{bb}$、$C_{bt}$、$C_{tb}$ 和 $C_{tt}$；不允许用“Nanbu 或 SDE 二选一”隐式关闭跨表示 pair；
- combined 能量交换和碰撞储库账本闭合；
- tail 热化不会导致不可控 PIC 核心人口；
- electron-electron bulk-tail 反作用使 combined 动量和能量闭合；
- electron-ion 或外部碰撞交换只进入对应 reservoir，不混入数值误差；
- Coulomb/Landau tail--tail Nanbu--Perez 等权重事件的加权三动量和相对论能量达到求和精度；
- 不等权重 Sentoku--Kemp 模式必须保持宏粒子数和权重不变、单事件加权相对论能量达到求和精度，并以多 cell/多 seed 集合报告动量均值、方差和长期随机游走；禁止把统计动量闭合误写为单事件严格闭合，也禁止出现 `conversion_N=0` 时由碰撞引起的粒子增长；
- 至少使用两组独立 counter-based seed 验证 SDE 弱矩或 Nanbu--Perez 散射统计处于预期采样误差，并进行 $dt/dt/2$ 碰撞收敛测试；
- `TraceStationaryBackground` 结果不得被标记为完整 electron--electron 守恒通过。
- manifest 中的 field boundary、collision backend、pair mask、weight mode、子步和 population controller 必须与实际运行一致；溯源不一致的结果只能作为调试数据。
- H9 Beam 40 fs 必须同时报告 tail 物理数分数、宏粒子数、最大 rank 占比、有效宏粒子数/噪声和 controller off/on A/B；仅仅“没有 OOM”不算资源验收通过。
- 阈值接口验收使用带明确 bin 边界的 $dN/dK$，不直接比较宽度可能不同的 raw count。

**当前H9状态（2026-08-11）**：17A--17F、Beam 12 fs阻断门、Beam 3 fs、no-Beam
40 fs、25.5 fs故障专项和Beam 40 fs controller-off生产门均已通过。最新40 fs运行
完成1564步且0 rejected step，转换、Gauss、tail账本、全局/局部粒子预算与wall-time
validator全部通过。转换失败共识已经接入所有生产conversion路径，25.5 fs与40 fs
均未再出现rank控制流分叉。

H9的**数值推进与既定资源门已通过**，但完整物理验收仍保留两项：使用最新
flux-interface快照重新计算6 MeV附近cell-volume bulk + PIC tail统一能谱；按真实增长率
给出120 fs粒子数、单rank内存与wall-time上界。旧static-cell运行中的人工谱谷不能直接
判到新结果上，也不能在未分析新快照前宣称已经消失。population controller继续关闭。

2026-08-08 最终补充：新12 fs全局源谱回归证明碰撞人口、Gauss和转换总矩稳定；
真实转换事件矩审计随后完成。阈值以下转换质量仅为舍入量级，说明阈值选择和事务
闭合正确；但低 $u_\perp$ 的 $\Pi_\perp$ 中心/体积表示差达到累计3.26%、单cell
11.3%，743个真实cell的展开非负六矩目标全部不可行。因此当前阶段门准确状态为：

```text
collision boundedness      PASS
conversion global moments  PASS
transaction/Poisson        PASS
edge-group conservation    PASS
global source-spectrum     PASS
source number chain        PASS
near-axis center support   FAIL/CONFIRMED
cell-volume spectrum       PASS
subcell solver controls    PASS
center-to-spread six-moment feasibility FAIL AT K
subcell PIC loading        FAIL/DISABLED
real-event moment audit    PASS/RED BRANCH
below-threshold mass       PASS (1.11e-16)
Jx/K/Pixx representation  PASS (about 1e-6)
Piperp representation     FAIL (3.26e-2 cumulative)
volume target feasibility FAIL (0/743)
audit empty-event filter  FAIL/DIAGNOSTIC ONLY
fine-spectrum fidelity     FAIL
H9 overall                 PARTIAL PASS
```

下一次状态升级只能由§7.11.16的诊断修复和固定 $(192,64)$ 局部速度网格离线原型
触发。原subcell六矩候选已经失败，不得继续执行其旧第11、12步；也不得把
$\Pi_\perp$ 误差通过全局能量补丁或增加粒子数掩盖。

### 19.4 最终生产级 [部分完成/待能量边界账与H10 A/B]

比较：

- 电场能量；
- combined 背景动能和平均动能；
- Beam 能量；
- 总能量账；
- combined 能谱和动量分布；
- 宏观波包包络、相位和传播区域；
- 参数收敛性。

EPOCH 的 PIC 噪声和本程序的混合噪声结构不同，不要求逐点波形重合。

2026-08-12 的controller-off 120 fs结果给出以下阶段结论：

- **通过**：120 fs有限推进；Gauss、Tail数目和conversion矩闭合；阈值combined谱没有
  人工硬谷；combined背景平均动能由约770 eV升至2265 eV；
- **待核对**：电场能由约$7.32\times10^7$升至$8.68\times10^7\ \mathrm{J/m^2}$，
  不能在缺少背景开放边界能流时单独判为场残留错误；
- **未完成**：累计碰撞储能和背景边界能流组成的完整总能量账；
- **旧controller-off基线未通过最终资源门**：低于6 MeV的Tail占Tail总数58.80%，导致
  1.4317亿宏粒子和约12 s/step的末期成本；
- **H10长期资源门通过但物理门未通过**：H10把末态Tail粒子数降低44.0%、平均wall time降低
  41.6%，但120 fs场相对$L_2$为16.29%、combined谱为5.18%、combined背景动能高3.13%。
  因此`none`和当前`hysteretic`都不是最终签字配置；
- **待外部对比**：宏观波包包络、相位、能谱和动量分布与EPOCH的统计级比较。不得用
  PIC噪声逐点差异作为失败标准。

最终签字条件更新为：H10返回事务、参数平台和资源收益已经完成；当前仍需§17.15关闭完整连续
能量账，并解释返回六矩表示差与长期宏观漂移。只有短A/B通过后才能扩展到111 fs，再决定是
否重跑新的120 fs最终生产。与EPOCH的统计级宏观物理验收位于上述内部一致性验收之后。

---

## 20. 自动编码模型执行纪律

### 20.1 每次只执行一个阶段

自动编码模型每次工作必须：

1. 阅读本阶段涉及的现有头文件和实现；
2. 列出将保留的接口和将新增的接口；
3. 先实现生产代码；
4. 测试直接调用生产实现，禁止复写离散公式；
5. 编译本阶段 target；
6. 运行本阶段最小验收；
7. 报告数值结果和未覆盖风险；
8. 未通过时只修当前阶段，不提前接入下一阶段。

### 20.2 禁止行为

- 不得为了编译通过创建空实现或始终返回成功；
- 不得放宽保护隐藏 conversion 不守恒；
- 不得修改 Beam 以迁就 tail；
- 不得把 tail 密度漏出 midpoint Poisson；
- 不得在 checkpoint 读取时静默重新采样粒子；
- 不得以全局能量修正弥补局部转换错误；
- 不得同时重写 remap、Poisson、Beam 和混合转换；
- 不得在 H4 前运行长生产并据此调参。

### 20.3 代码审查清单

每阶段结束时确认：

- 数组所有权和生命周期清楚；
- accepted/trial 状态没有别名污染；
- MPI collective 顺序在所有 rank 一致；
- OpenMP 循环无共享写冲突；
- 粒子 ID 和 RNG 可重启；
- 诊断只读取已接受态；
- 所有新增 CLI 参数写入 manifest；
- 所有新文件加入显式 CMake target；
- 没有遗留 TODO 作为生产路径。

---

## 21. 最终架构结论

最终生产程序应是：

```text
开放 reservoir 背景 bulk Eulerian Vlasov/FP
        +
同一背景物种的 1D3V PIC 高能 tail
        +
开放注入/开放流出的 Beam PIC
        +
非周期 Gauss/Poisson 静电场
        +
按物理核选择的 bulk-FP/BGK 后端、tail--tail Nanbu--Perez 与 tail--bulk SDE/反作用
        +
跨 bulk--tail 表示的局部反作用与守恒账本
```

该架构保留了当前已经验证的开放空间边界、非周期场、保守 remap、Beam 推进和事务调度，同时避免用全局高分辨率速度网格追逐稀疏高能尾部。

它能否成为最终最优实现，取决于三个可测条件：

1. bulk-to-tail 转换是否在离散电荷、数目、动量和能量上闭合；
2. tail 粒子数量和噪声是否在 120 fs 资源预算内；
3. 碰撞开启后是否需要双向返回。

在这三个条件通过前，不应宣称重构完成，也不应继续围绕旧阶段 5 的 $u_{\max}$ 扩展方案做局部修补。

EPOCH 4.20.1 对最终实现的定位是：它为 tail 的相对论 1D3V 粒子推进、staggered 场插值、形函数差分电流、开放粒子边界、cell-local 随机配对和 Coulomb Nanbu--Perez 碰撞提供源码级参考；它不替代本项目的 Eulerian bulk、Poisson 场、表示转换、跨表示碰撞反作用和事务接受机制。
