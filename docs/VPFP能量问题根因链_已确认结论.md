# VPFP 能量问题根因链：已确认结论

> 文档性质：历史审计结论的单一入口，不是新的生产修改任务书。
>
> 证据来源：`docs/archive/VPFP公共能量残差_离散功配对审计与修复实施方案.md`、
> `docs/archive/VPFP场粒子离散功同源闭合_根因判别与修复实施方案.md`、
> `docs/archive/VPFP时间中心不一致_JC修复与验收实施方案.md`、
> `docs/archive/VPFP_P3联合x_u_Poisson功通量重构实施方案.md`。
>
> 当前后续实现入口：`docs/VPFP_J2开放背景联合中点离散与验收实施方案.md`。

## 1. 最终结论

旧生产推进的公共能量残差属于**结构性离散不闭合**。它不是某个能量诊断漏记，也不是单纯把时间步缩小即可消除的普通时间截断误差。

根因是旧背景电子推进将下列三个对象分别构造：

1. `x` 输运有限体积通量产生的电荷连续性电流 $J_N$；
2. `u` 方向受力推进产生的动能变化及其功电流 $J_E$；
3. 中点 Poisson 场、粒子受力场和末态 Poisson 场。

它们不属于同一个时间中心的相空间离散算子，之间也没有共同的离散乘积法则、伴随关系或反对称 Poisson bracket。因此，即使每个子算子各自守恒，整个背景 Vlasov--Poisson/VPFP 更新也没有必然满足：

$$
\Delta K_{\rm bulk}+\Delta U_E-W_{\rm electrode}=0
$$

或包含 Beam、Tail、转换和碰撞后的完整离散能量恒等式。

## 2. 观测到的宏观后果

在长时间生产计算中，公共能量余额逐步积累，并在非线性阶段表现为：

- 电场能量没有按照参照计算中的能量转移趋势充分下降；
- 背景电子总能量和平均动能偏低，能量留在场或离散公共余额中；
- 波包包络、相位和后期能量分配逐步偏离；
- 中点迭代与 FCT active set 在后期更难收敛，计算成本显著增加。

这些现象相互关联，但宏观波形偏差不是根因本身；其上游原因是每个时间步重复产生的离散场--粒子功缺口。

## 3. 直接证据链

### 3.1 公共余额集中于场--粒子配对块

阶段功率分解显示，公共余额几乎完全落在：

$$
R_{\rm field\text{-}particle}
=\Delta U_E+W_{\rm bulk}+W_{\rm tail}+W_{\rm beam}-W_{\rm electrode}.
$$

归档 Gate H 的粗、细时间步窗口中，该项解释完整能量余额约 $99\%$。这说明剩余能量不是散落在多个独立小误差中，而是集中在场能变化与粒子受力功的共同配对关系。

### 3.2 组分账本与源项所有权不是主导缺口

以下对象已分别通过或被证明不足以解释公共余额：

- Bulk、Tail、Beam 的局部数目连续性与接受态 source/owner 对应；
- Beam 轨迹沉积、开放出流、尾粒子转换和 return 事务；
- Tail kick、Beam kick、bulk $u_\parallel$ 通量的局部功账；
- MPI shared-face/owner 合并与 checkpoint/restart 事务；
- Poisson/Gauss 空间恒等式及电极功边界项。

这些模块仍必须持续回归，但当前证据不支持把它们视为公共能量残差的主因。

### 3.3 Poisson 空间算子不是修复目标

固定状态的 Poisson work identity 审计通过，表明当前非周期 Poisson 空间离散、Dirichlet 电势边界和电荷--场能关系本身没有构成主导缺口。

因此禁止以如下方式“修复”能量账：

- 改写 Poisson 空间差分以拟合能量；
- 每步全局能量投影、全局电流补偿或人为电场缩放；
- 强制扣除电荷/电流零模；
- 以全局能量补丁覆盖局部离散不一致。

这些操作会改变所求物理系统，且掩盖真正的时间中心/相空间通量缺口。

### 3.4 缩小时间步不能被当作根治方案

历史 $dt$、$dt/2$、$dt/4$ 对照没有给出公共残差按预期阶数稳定下降的证据；随机碰撞和转换事件会使严格 Richardson 阶数比较受污染，但两档中场--粒子配对块仍保持主导。

因此，减小 $dt$ 可以降低部分局部误差或改善稳定性，却不能替代离散结构修复。将生产时间步无限缩小只会以不可接受的性能代价掩盖问题。

## 4. 离散层面的根因链

旧生产步本质上是分裂结构：

$$
T_x(\Delta t/2)
\;\rightarrow\;
T_u(E,\Delta t)
\;\rightarrow\;
T_x(\Delta t/2),
$$

其问题链如下：

1. `x` 有限体积扫掠通量可严格给出密度变化和 $J_N$，但该通量并不自动等于受力功所需的能量共轭电流；
2. `u` 受力通量可给出动能变化和 $J_E$，但其时间层、分布状态及离散速度导数与 `x` 电流构造独立；
3. 中点场用于部分受力/试探计算，而末态电荷再通过 Poisson 得到末态场；
4. 因而 $E^{n+1/2}$、$f^{n+1/2}$、$J_N^{n+1/2}$、$J_E^{n+1/2}$ 不是由同一隐式相空间残差系统共同解出；
5. 最终 $\Delta U_E$ 与粒子受力功不构成一个离散 telescope identity，产生每步公共功缺口；
6. 该缺口长期累计，最终改写电场--背景电子之间的真实能量交换。

FCT、PPM、正性限制、Tail 转换和边界强梯度会改变缺口的局部分布、非线性收敛性和可见宏观后果，但它们不是上述公共余额的第一性来源。

## 5. 为什么 P3 局部修正路线终止

P3 路线曾尝试通过 $u$ 方向能量共轭离散速度表、局部 $x/u/Poisson$ 功通量修正来闭合缺口。其审计结果确认：

- 单速度功电流的符号、单位、cell volume、MPI shared-face 和基本时间层口径没有错误；
- 调整局部速度矩无法使联合 $x$ 输运、$u$ 受力与场能变化成为同一个离散系统；
- 局部 dual/current correction 若不随网格收敛，就会成为非一致物理补丁。

故 P3 不应继续叠加局部修正。它的价值是排除了“只需修一个离散速度表或某个 Poisson 面映射”的假设。

## 6. 已排除与仍需保持的边界

### 已排除为主根因

- 非周期 Poisson 的空间算子或 Dirichlet 电势边界；
- Beam 开放注入/出流的轨迹沉积账本；
- Tail--bulk 转换、Tail return 或碰撞反作用账本；
- 单纯 MPI owner/face 合并错误；
- 单纯诊断公式错误；
- 单纯时间步过大。

### 仍需保留为独立质量门

- 所有物种连续性、Gauss/Poisson 恒等式、局部功账和重启等价性；
- FCT 正性和速度域边界对宏观波形的影响；
- Tail 转换阈值、粒子数资源门及碰撞统计收敛；
- 非线性中点求解的正域收敛与事务接受。

这些质量门可能造成额外误差或性能问题，但不能被用来取代本根因链。

## 7. 当前唯一正确的修复方向

应构建并验证新的联合 $x/u/Poisson$ 时间中心相空间离散：

1. 在同一候选 $f^{n+1/2}$ 与 $E^{n+1/2}$ 上定义 `x` 和 `u` 面通量；
2. 由最终的同一通量同时构造密度连续性、charge current 与 kinetic-work current；
3. 使用与 Poisson/场能离散一致的场时间中心和边界电极功；
4. 将背景分布、场和必要的 Lagrange/Poisson 变量置于同一个非线性残差系统中求解；
5. 只接受有限、正域可接受、Gauss 通过且联合能量残差通过的候选态；
6. Beam、Tail、转换、碰撞和 FCT 先保持为冻结的外部事务，待 bulk 联合核心闭合后再逐项接入。

当前具体执行和验收要求由：

- `docs/VPFP_联合相空间时间中心能量闭合重构实施方案.md`；
- `docs/J1联合中点当前失败分析_2026-08-21.md`

共同规定。

## 8. 当前状态与决策规则

联合中点 J1 原型已证明其非线性/Poisson 残差可以下降，但仍存在两个未解决问题：

1. unit 中联合能量闭合仍触发 code 75，说明当前 center-trace flux 尚未建立真正的公共能量恒等式；
2. smoke 的 line search 曾接受巨大负质量候选，说明非线性求解仍未保持正域。

在定位 J1 的首个不闭合能量项并修复正域求解前：

- 不进入 J2/J3；
- 不将 joint mode 设为生产默认；
- 不继续修改归档的 JC/P3 局部修补路线；
- 不使用能量投影、全局补偿或仅缩小 $dt$ 作为替代修复。

## 9. 2026-08-22 F10 对 J1 根因链的最新约束

最新 F10 单 rank 结果并未推翻本文件关于旧 Strang/PPM 结构性不闭合的结论，但将**当前 J1 center-trace 原型**的首个可测缺口显著收窄。

平衡 `smooth-background` 已通过；只有非平凡 `smooth-perturbed-background` 触发 code 75。其能量分解为：

$$
W_u=28.891863089349869,
\qquad
W_F=28.891863089349886,
$$

$$
W_J=22.224220061661196,
\qquad
R_{uJ}=6.6676430276886727\ \mathrm{J/m^2},
$$

$$
R_{PJ}=-5.2950189122213942\times10^{-6}\ \mathrm{J/m^2}.
$$

其中 $W_F$ 使用同一 accepted midpoint mass、同一 pairing cell field 和生产 `vH` 构造。因此：

$$
W_u-W_F\simeq-1.7\times10^{-14}\ \mathrm{J/m^2},
$$

说明当前 J1 的 u-face 功、u trace、`vH` 和 midpoint mass 已通过其直接共轭检查。主缺口现在确定在：

$$
W_F-W_J,
$$

即 cell force current 到由 x flux 构造的 `charge_current_face` 的映射、face averaging/ownership 或其与 pairing face 的离散内积。

这不是修改 Poisson 空间算子、放宽 Newton 或调节能量 gate 的依据。并且，已记录的 endpoint `W_J` 本身接近零并不能排除周期接缝影响首末 cell force work；在没有 boundary-cell force-work 分解前，不得将问题简单归为“纯内部面错误”。

所以当前 J1 的下一步严格遵循主方案第 11 节情形 B：审计生产 `x_flux_rate`、`charge_current_face`、midpoint face trace、首末 cell/seam 映射和 MPI face ownership。不得进入 F11、J1 MPI、J2 或 J3。

## 10. 2026-08-23 J1 最终闭合结论

本节是当前状态的最终结论，覆盖第 8、9 节中“J1 尚未解决”的时效性描述；第 8、9 节继续保留为历史定位记录。

### 10.1 情形 B：periodic seam 加权伴随缺口

J1 manufactured 测试使用 periodic x current topology，但 Poisson pairing face 保留两个非周期端点和 endpoint half-weight。旧 cell gather：

$$
E_i^{\mathrm{cell}}=\frac{E_{i-1/2}+E_{i+1/2}}{2}
$$

不是 periodic current map 在该 face quadrature 下的加权转置，导致：

$$
W_F-W_J=6.6676430276886904\ \mathrm{J/m^2}.
$$

seam 解析预测在舍入误差内完全解释该缺口。修复为 J1 测试拓扑专用的：

$$
E_{\mathrm{cell}}=G^*E_{\mathrm{pair}}.
$$

修复后：

$$
W_u-W_F=O(\epsilon_{\mathrm{mach}}),
\qquad
W_F-W_J=O(\epsilon_{\mathrm{mach}}).
$$

该 helper 只能用于 J1 periodic manufactured topology，不能直接用于最终 open/reservoir background。

### 10.2 情形 A：近中性电荷装配消减

J1 原先分别构造：

$$
\rho^n=q_e(n_i-n_e^n),
\qquad
\rho^{n+1}=q_e(n_i-n_e^{n+1}),
$$

再计算两者之差。近中性大数相减使电荷连续性诊断出现约：

$$
3.26\times10^{-11}\ \mathrm{C/m^2}
$$

的装配舍入项。J1 固定离子、无 Beam/Tail/source 时，改用数学等价的稳定增量装配：

$$
\rho_i^{n+1}
=\rho_i^n+
\frac{q_e}{\Delta x}
\sum_{j,k}\left(M_{i,j,k}^{n+1}-M_{i,j,k}^{n}\right).
$$

装配 mismatch 降至约：

$$
7.5\times10^{-17}\ \mathrm{C/m^2}.
$$

这不是能量补丁，也不是由 current divergence 覆盖电荷；candidate rho 仍由 candidate mass 决定。

### 10.3 Poisson scalar identity 的稳定求和

旧 `evaluate_work_identity()` 先分别求两个约为：

$$
1.39\times10^5\ \mathrm{J/m^2}
$$

的场能总量，再相减得到约：

$$
17.78\ \mathrm{J/m^2}.
$$

普通 double 顺序求和产生抵消放大。修复采用逐 cell 因式分解的场能差、`long double` Neumaier 累加和 `MPI_LONG_DOUBLE` 归约；未修改 Poisson solve、stencil、边界或场状态。

零端点生产配置的 scalar identity 已通过。非零 Dirichlet endpoint 的完整离散边界功仍是独立已知限制，但不属于当前固定：

```text
phi_left=0
phi_right=0
```

的生产模型。

### 10.4 情形 A-N：非线性 residual 的功投影

phase residual 与 Poisson residual 已经很小，但 potential-weighted continuity residual 仍可超出最终能量门。J1 收敛条件因此由两门扩展为：

```text
phase residual
Poisson residual
Poisson-current pairing residual
```

三者必须同时通过。pairing 门最终保持：

$$
10^{-9},
$$

总能量门保持：

$$
10^{-8}.
$$

没有通过缩放场、电流、分布或能量投影改变候选态。

### 10.5 多步 signed-state 契约

Newton trial 使用 signed residual domain，最终候选由 code-76 正性门验收。第一步接受态可能包含 code-76 容差内的 signed roundoff mass，因此下一步初始 candidate 也必须允许同一 signed residual domain；否则第二步会在初始 residual 阶段错误返回 code 71。

修复只统一 residual 定义域，最终 code-76、`negative_tolerance` 和拒绝逻辑保持不变。

### 10.6 最终验收结果

以下均已通过：

- J0 全部离散恒等式；
- J1 periodic seam weighted-adjoint；
- Poisson scalar identity；
- charge continuity 与 potential-weighted prediction；
- `smooth-background`；
- `smooth-perturbed-background`；
- `dt` 与 `dt/2`；
- 10 步累计残差门；
- 1/2/5 rank MPI ownership；
- 情形 A 全量 A8 回归。

MPI rank 比较必须区分两类尺度：结构恒等式按舍入误差验收；非线性停止量按 solver tolerance 验收。不同 rank 的迭代次数允许不同，只要都进入同一个收敛球且结构恒等式保持闭合。

### 10.7 当前唯一下一步

J1 只是 periodic manufactured 核心测试。最终生产模型必须回到：

```text
reservoir/open background
nonperiodic DIRICHLET_PHI Poisson
phi_left=0
phi_right=0
```

因此下一步是 J2：为开放 x flux、reservoir source、边界 number/current/kinetic ledger 和非周期 Poisson pairing 重新建立同源离散。禁止把 J1 periodic `G^*` 直接接入生产开放边界。
