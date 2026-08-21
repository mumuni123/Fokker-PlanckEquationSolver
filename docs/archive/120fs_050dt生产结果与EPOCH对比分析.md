# 120 fs、0.5 dt 生产结果与 EPOCH 对比分析

## 1. 分析对象

本次分析对比以下两组无碰撞结果：

- 当前程序生产结果：
  `results/production_120fs_050dt_20260727_173708`
- EPOCH 参考结果：
  `EPOCH-Compare/noCollision`
- 当前程序对应的原始结果与诊断：
  `output/production_120fs_050dt_20260727_173708`

两组结果都包含从约 $0$ 到 $120\ \mathrm{fs}$ 的 $201$ 个快照，输出间隔约为
$0.6\ \mathrm{fs}$，编号可以按 `00000-00200` 一一对应。

为排除周期接缝和端点尖峰对核心区统计的主导影响，本文将
$0.2\ \mu\mathrm{m}\le x\le7.8\ \mu\mathrm{m}$ 定义为核心区。归一化参数采用：

- 电场归一化值：$E_0=2.990908\times10^{11}\ \mathrm{V/m}$；
- 背景电子初始密度：$n_{e0}=1.2\times10^{29}\ \mathrm{m^{-3}}$。

## 2. 总体结论

本次 $0.5\,dt$ 长跑具有清晰的时间分段：

1. $0-12\ \mathrm{fs}$：当前程序与 EPOCH 定量接近，主波包位置、幅度和能量演化均较好。
2. $12-24\ \mathrm{fs}$：宏观量仍接近，但局部相位开始出现差异。
3. $24-36\ \mathrm{fs}$：电场总幅度仍接近，但背景电子加热已经偏低。
4. $36-48\ \mathrm{fs}$：进入系统性非线性分叉。
5. $48-72\ \mathrm{fs}$：能量账快速恶化，高能尾部发展不足，电场衰减明显慢于 EPOCH。
6. $72-120\ \mathrm{fs}$：当前程序维持过强、过于相干的电场，并出现较强的 Beam 滞留和再聚束。

因此，$60\ \mathrm{fs}$ 之后的结果不能再视为 EPOCH 的定量复现。该差异也不能全部归因于
PIC 噪声与 Eulerian Vlasov 方法之间的自然差别，因为背景平均动能、高能尾部、核心区场能和总能量账均出现了系统性偏差。

## 3. 电场定量对比

核心区电场均方根定义为：

$$
E_{\mathrm{rms,core}}
=
\sqrt{\frac{1}{N_{\mathrm{core}}}
\sum_{i\in\mathrm{core}}E_{x,i}^2}.
$$

下表中的“核心区 $E_x$ RMS 比值”表示当前程序除以 EPOCH：

| 时间 | 核心区 $E_x$ RMS 比值 | 当前程序背景平均动能 | EPOCH 背景平均动能 | 当前程序累计能量账误差 |
|---:|---:|---:|---:|---:|
| $12\ \mathrm{fs}$ | $1.02$ | $246\ \mathrm{eV}$ | $253\ \mathrm{eV}$ | $-0.67\times10^6\ \mathrm{J}$ |
| $24\ \mathrm{fs}$ | $0.98$ | $1466\ \mathrm{eV}$ | $1542\ \mathrm{eV}$ | $-4.29\times10^6\ \mathrm{J}$ |
| $36\ \mathrm{fs}$ | $1.00$ | $2328\ \mathrm{eV}$ | $2820\ \mathrm{eV}$ | $-9.71\times10^6\ \mathrm{J}$ |
| $48\ \mathrm{fs}$ | $1.38$ | $3206\ \mathrm{eV}$ | $4314\ \mathrm{eV}$ | $-1.73\times10^6\ \mathrm{J}$ |
| $60\ \mathrm{fs}$ | $2.05$ | $3313\ \mathrm{eV}$ | $4623\ \mathrm{eV}$ | $-73.5\times10^6\ \mathrm{J}$ |
| $72\ \mathrm{fs}$ | $2.15$ | $3002\ \mathrm{eV}$ | $4744\ \mathrm{eV}$ | $-134.9\times10^6\ \mathrm{J}$ |
| $96\ \mathrm{fs}$ | $2.84$ | $3374\ \mathrm{eV}$ | $4871\ \mathrm{eV}$ | $-145.9\times10^6\ \mathrm{J}$ |
| $120\ \mathrm{fs}$ | $3.59$ | $3348\ \mathrm{eV}$ | $4909\ \mathrm{eV}$ | $-146.3\times10^6\ \mathrm{J}$ |

主要现象如下：

- $12\ \mathrm{fs}$ 时，核心区电场 RMS 只相差约 $2\%$。
- $24\ \mathrm{fs}$ 时，整体幅值仍接近，但主波包局部相位已经发生移动。
- $36\ \mathrm{fs}$ 时，两者核心区电场能几乎相同，但背景平均动能已相差约 $17\%$。
- $48\ \mathrm{fs}$ 后，当前程序保留明显更强的相干波包。
- $60-120\ \mathrm{fs}$，EPOCH 的电场逐渐衰减并转为低幅、宽带的非规则结构，而当前程序仍保持清晰、规则且幅度较高的波列。

关键图像：

- 当前程序电场：`results/production_120fs_050dt_20260727_173708/fields`
- EPOCH 电场：`EPOCH-Compare/noCollision/Electric_Field_Ex`
- EPOCH 电场时空图：`EPOCH-Compare/noCollision/space_time_maps/Ex_space_time.png`

## 4. 电场能量演化

电场能量按以下公式计算：

$$
U_E
=
\frac{\epsilon_0}{2}
\sum_i E_{x,i}^2\Delta x.
$$

约 $30\ \mathrm{fs}$ 时，两者电场能峰值都接近：

$$
U_E^{\mathrm{peak}}\simeq3.7\times10^8\ \mathrm{J}.
$$

但到 $120\ \mathrm{fs}$：

| 量 | 当前程序 | EPOCH |
|---|---:|---:|
| 电场能量 | $1.001\times10^8\ \mathrm{J}$ | $6.92\times10^6\ \mathrm{J}$ |
| 背景电子平均动能 | $3.348\ \mathrm{keV}$ | $4.909\ \mathrm{keV}$ |
| 域内总能量 | $6.152\times10^8\ \mathrm{J}$ | $7.670\times10^8\ \mathrm{J}$ |

当前程序并不是把过多能量全部转移给了背景电子。实际情况恰好相反：

1. 背景电子获得的动能不足；
2. 大量能量继续滞留在相干电场中；
3. 系统同时出现约 $1.46\times10^8\ \mathrm{J}$ 的累计能量缺口。

当前程序的能量账定义为：

$$
E_{\mathrm{accounted}}
=
E_{\mathrm{total}}
-E_{\mathrm{beam,in}}
+E_{\mathrm{beam,out}}
-E_{\mathrm{collision}}.
$$

累计能量账残差为：

$$
R_E(t)
=
E_{\mathrm{accounted}}(t)
-E_{\mathrm{accounted}}(0).
$$

到 $120\ \mathrm{fs}$：

$$
R_E\simeq-1.4625\times10^8\ \mathrm{J}.
$$

这一缺口主要在 $48-72\ \mathrm{fs}$ 快速形成，与波形开始明显偏离的时间段一致。

相关图像：

- `results/production_120fs_050dt_20260727_173708/electron_kinetics/electric_field_energy_history.png`
- `results/production_120fs_050dt_20260727_173708/electron_kinetics/domain_energy_history.png`
- `results/production_120fs_050dt_20260727_173708/electron_kinetics/energy_balance_history.png`
- `EPOCH-Compare/noCollision/energy_history/field_energy_history.png`
- `EPOCH-Compare/noCollision/energy_history/domain_energy_history.png`

## 5. 边界尖峰不是核心偏差的主要来源

虽然当前程序在后期的左右端点出现了比 EPOCH 更强的尖峰，但将
$x<0.2\ \mu\mathrm{m}$ 和 $x>7.8\ \mu\mathrm{m}$ 排除后，核心区差异仍然很大。

$120\ \mathrm{fs}$ 时：

- 当前程序核心区场能：$8.52\times10^7\ \mathrm{J}$；
- EPOCH 核心区场能：$6.62\times10^6\ \mathrm{J}$；
- 核心区场能比值约为 $12.9$；
- 当前程序边界区域约占自身总场能的 $14.9\%$；
- EPOCH 边界区域约占自身总场能的 $4.2\%$。

当前程序相对于 EPOCH 的场能超额约为：

$$
\Delta U_E
\simeq
9.32\times10^7\ \mathrm{J}.
$$

其中约 $84\%$ 位于核心区。因此：

> 边界尖峰会增加总场能并污染近边界区域，但不能解释核心区电场在后期高出 EPOCH 数倍这一主要差异。

## 6. 背景密度对比

核心区背景密度扰动 RMS 定义为：

$$
\delta n_{\mathrm{rms}}
=
\sqrt{
\frac{1}{N_{\mathrm{core}}}
\sum_{i\in\mathrm{core}}
\left(
\frac{n_{e,i}}{n_{e0}}-1
\right)^2
}.
$$

对比结果如下：

| 时间 | 当前程序 $\delta n_{\mathrm{rms}}$ | EPOCH $\delta n_{\mathrm{rms}}$ |
|---:|---:|---:|
| $12\ \mathrm{fs}$ | $0.0188$ | $0.0210$ |
| $36\ \mathrm{fs}$ | $0.0935$ | $0.1216$ |
| $60\ \mathrm{fs}$ | $0.0661$ | $0.0898$ |
| $72\ \mathrm{fs}$ | $0.0520$ | $0.1045$ |
| $96\ \mathrm{fs}$ | $0.0380$ | $0.1139$ |
| $120\ \mathrm{fs}$ | $0.0339$ | $0.1110$ |

EPOCH 是 PIC 程序，其后期密度 RMS 中包含粒子噪声，因此不能要求两者逐点相同。但趋势仍然明确：

- 当前程序的核心电场更强；
- 当前程序的密度细结构却更弱、更规则；
- 当前程序没有像 EPOCH 那样形成足够宽的相空间相混和高能响应。

这说明后期差异已经进入背景电子动力学，而不是单纯的绘图噪声或端点电场问题。

## 7. Beam 密度对比

Beam 的演化可以分为三个阶段：

### 7.1 $12-48\ \mathrm{fs}$

- 注入前沿位置基本相同；
- Beam 聚束的主要空间范围接近；
- 当前程序峰值略高，但仍处于相同数量级。

### 7.2 $60-72\ \mathrm{fs}$

- EPOCH 中 Beam 主要继续向右端稀释和出流；
- 当前程序在 $x\simeq5-7\ \mu\mathrm{m}$ 形成更明显的滞留和再聚束；
- $72\ \mathrm{fs}$ 时，当前程序 Beam 峰值约为 EPOCH 的两倍。

### 7.3 $96-120\ \mathrm{fs}$

- 当前程序在整个计算域中仍保留明显的 Beam 密度；
- $96\ \mathrm{fs}$ 时，局部峰值约为 EPOCH 的 $3-4$ 倍；
- $120\ \mathrm{fs}$ 时，当前程序 Beam 密度的空间形状仍与 EPOCH 明显不同。

这形成如下反馈：

$$
E_x\text{ 衰减不足}
\rightarrow
\text{Beam 减速、俘获或再聚束}
\rightarrow
J_{\mathrm{beam}}\text{ 持续存在}
\rightarrow
E_x\text{ 继续维持}.
$$

因此，Beam 的后期差异更可能是过强相干场所形成的反馈结果，而不是最初的单一误差源。

## 8. 背景电子能谱

能谱统计表明，当前程序的主要问题不是低能电子数量严重错误，而是高能尾部发展不足。

### 8.1 高能电子比例

| 时间 | 程序中 $E>100\ \mathrm{keV}$ | EPOCH 中 $E>100\ \mathrm{keV}$ | 程序中 $E>1\ \mathrm{MeV}$ | EPOCH 中 $E>1\ \mathrm{MeV}$ |
|---:|---:|---:|---:|---:|
| $49.8\ \mathrm{fs}$ | $1.35\times10^{-3}$ | $3.84\times10^{-3}$ | $1.53\times10^{-4}$ | $4.15\times10^{-4}$ |
| $75\ \mathrm{fs}$ | $1.81\times10^{-3}$ | $4.43\times10^{-3}$ | $2.00\times10^{-4}$ | $5.90\times10^{-4}$ |
| $100.2\ \mathrm{fs}$ | $1.48\times10^{-3}$ | $4.60\times10^{-3}$ | $1.23\times10^{-4}$ | $5.89\times10^{-4}$ |
| $120\ \mathrm{fs}$ | $9.04\times10^{-4}$ | $4.57\times10^{-3}$ | $3.80\times10^{-5}$ | $5.81\times10^{-4}$ |

$120\ \mathrm{fs}$ 时：

- 当前程序中 $E>100\ \mathrm{keV}$ 的比例约为 EPOCH 的 $1/5$；
- 当前程序中 $E>1\ \mathrm{MeV}$ 的比例约为 EPOCH 的 $1/15$。

### 8.2 高能窗口平均能量

对 $E>20\ \mathrm{keV}$ 的电子：

| 时间 | 当前程序平均能量 | EPOCH 平均能量 |
|---:|---:|---:|
| $49.8\ \mathrm{fs}$ | $58.4\ \mathrm{keV}$ | $120.9\ \mathrm{keV}$ |
| $75\ \mathrm{fs}$ | $59.9\ \mathrm{keV}$ | $119.6\ \mathrm{keV}$ |
| $100.2\ \mathrm{fs}$ | $49.2\ \mathrm{keV}$ | $112.8\ \mathrm{keV}$ |
| $120\ \mathrm{fs}$ | $35.7\ \mathrm{keV}$ | $109.7\ \mathrm{keV}$ |

当前程序中 $E>20\ \mathrm{keV}$ 的电子数量仍在增加，但其平均能量在
$75\ \mathrm{fs}$ 后反而下降。这表明电子被堆积在较低的高能窗口中，却不能持续进入
$100\ \mathrm{keV}$ 以上和 MeV 能段。

相关图像：

- `results/production_120fs_050dt_20260727_173708/electron_kinetics/electron_energy_spectrum.png`
- `results/production_120fs_050dt_20260727_173708/electron_kinetics/energy_spectrum_low`
- `results/production_120fs_050dt_20260727_173708/electron_kinetics/energy_spectrum_high`
- `EPOCH-Compare/noCollision/energy_spectrum_low`
- `EPOCH-Compare/noCollision/energy_spectrum_high`

## 9. 平行与垂直动量分布

$120\ \mathrm{fs}$ 时平行动量均方根为：

$$
u_{\parallel,\mathrm{rms}}^{\mathrm{program}}
\simeq0.116,
$$

$$
u_{\parallel,\mathrm{rms}}^{\mathrm{EPOCH}}
\simeq0.186.
$$

对于 $|u_\parallel|>8$ 的粒子：

$$
f_{|u_\parallel|>8}^{\mathrm{program}}
\simeq4.2\times10^{-8},
$$

$$
f_{|u_\parallel|>8}^{\mathrm{EPOCH}}
\simeq1.13\times10^{-4}.
$$

当前程序的极高平行动量尾部比 EPOCH 低约三个数量级。

与此同时，当前程序的垂直动量分布在 $0-120\ \mathrm{fs}$ 基本重合。这与当前无碰撞、
纯纵向 $E_x$ 作用的模型一致，因为主要能量交换应发生在平行动量方向。

因此，主要差异可以进一步定位为：

> 当前程序的平行速度空间加速、相混或高能尾部输运受到抑制，而不是垂直方向出现了非物理加热。

## 10. FCT 介入情况

区域诊断表明，全局 limiter 覆盖率较高，但绝大部分发生在速度尾部。

### 10.1 $u$ 方向

| 时间 | 速度核心 active fraction | 速度尾部 active fraction |
|---:|---:|---:|
| $36\ \mathrm{fs}$ | $3.29\%$ | $89.1\%$ |
| $48\ \mathrm{fs}$ | $2.06\%$ | $90.1\%$ |
| $72\ \mathrm{fs}$ | $0.63\%$ | $90.9\%$ |
| $120\ \mathrm{fs}$ | $0.38\%$ | $90.2\%$ |

### 10.2 $x$ 方向

| 时间 | 速度核心 active fraction | 速度尾部 active fraction |
|---:|---:|---:|
| $36\ \mathrm{fs}$ | $1.58\%$ | $39.2\%$ |
| $48\ \mathrm{fs}$ | $0.90\%$ | $38.8\%$ |
| $72\ \mathrm{fs}$ | $0.11\%$ | $35.5\%$ |
| $120\ \mathrm{fs}$ | $0.054\%$ | $34.1\%$ |

速度尾部长期存在 $\alpha=0$ 的面。这意味着部分高阶通量完全退化到低阶通量。

虽然 limiter 在速度核心中的覆盖率并不高，但不能因此认为它对宏观物理无影响。原因是：

- 高能尾部粒子数较少；
- 但其单粒子动能很高；
- 场能向粒子能的转移和 Landau 型相混对速度尾部十分敏感；
- 尾部的低阶退化可以明显改变能量矩，而不显著改变总粒子数。

因此，当前结果支持以下判断：

> FCT 没有大范围直接压制低能核心分布，但它对高能尾部的长期强介入很可能抑制了背景电子吸收场能和形成高能尾部的过程。

## 11. 中点迭代与软接受

首次软接受出现在：

$$
t\simeq37.69\ \mathrm{fs}.
$$

各时间窗内的软接受比例和平均耦合迭代次数为：

| 时间窗 | 软接受比例 | 平均迭代次数 |
|---:|---:|---:|
| $0-12\ \mathrm{fs}$ | $0$ | $8.35$ |
| $12-24\ \mathrm{fs}$ | $0$ | $9.32$ |
| $24-36\ \mathrm{fs}$ | $0$ | $9.27$ |
| $36-48\ \mathrm{fs}$ | $51.8\%$ | $25.3$ |
| $48-60\ \mathrm{fs}$ | $30.4\%$ | $17.5$ |
| $60-72\ \mathrm{fs}$ | $7.1\%$ | $13.4$ |
| $72-96\ \mathrm{fs}$ | $51.8\%$ | $31.5$ |
| $96-120\ \mathrm{fs}$ | $80.4\%$ | $37.2$ |

$120\ \mathrm{fs}$ 时：

- 场残差约为 $7.3\times10^{-8}$；
- 背景电流残差约为 $7.7\times10^{-5}$；
- 背景电流严格阈值为约 $10^{-5}$；
- 耦合迭代达到 $40$ 次上限并软接受。

这说明后期收敛瓶颈主要来自 $J_{\mathrm{bkg}}$，而不是场变量本身。

软接受会使最终的 $f^{n+1}$、$J_{\mathrm{bkg}}^{n+1/2}$ 和
$E_x^{n+1}$ 不再属于同一个严格收敛的中点解，从而累计相位和能量交换误差。

但软接受不是全部能量缺口的唯一来源。$60-72\ \mathrm{fs}$ 的软接受比例只有
$7.1\%$，能量账却仍然快速恶化。因此还必须检查 FCT 最终通量及速度尾部边界的能量闭合。

## 12. 有限速度域与高能边界

当前程序的平行动量范围约为：

$$
u_\parallel\in[-10,10],
$$

而 EPOCH 结果中的尾部可以延伸到约：

$$
u_x\in[-12,12].
$$

当前程序高能谱的最高网格能量约为 $6.36\ \mathrm{MeV}$，EPOCH 高能诊断延伸到约
$8.0\ \mathrm{MeV}$。当前程序的高能谱还存在较明显的离散梳齿，说明有限的二维速度网格投影已经影响高能端分辨率。

背景电子总数到 $120\ \mathrm{fs}$ 下降约：

$$
\Delta N_{\mathrm{bkg}}
\simeq-2.50\times10^{20},
$$

相对变化约为：

$$
\frac{\Delta N_{\mathrm{bkg}}}{N_{\mathrm{bkg},0}}
\simeq-2.6\times10^{-4}.
$$

粒子数变化很小，不能单独解释宏观密度。但如果将最终能量缺口除以该粒子数变化，可得到：

$$
\frac{1.46\times10^8\ \mathrm{J}}
{2.50\times10^{20}}
\simeq3.65\ \mathrm{MeV/electron}.
$$

该能量尺度接近当前速度域的高能边界。因此，速度边界或高能 remap/FCT 对少量高能粒子的处理是一个强嫌疑点。

不过，现有 `step_diagnostics.dat` 中的速度边界 loss 列没有给出相应的非零累计量。故目前只能得出：

> 有限速度域及高能边界处理很可能参与了能量缺口，但尚未由现有诊断单独证实，不能直接断言所有缺失能量都来自速度边界删除。

## 13. 根因链

目前证据最支持以下分层根因链：

$$
\text{进入非线性阶段}
\rightarrow
\text{速度尾部高阶通量产生强梯度和细结构}
\rightarrow
\text{尾部 FCT 长期大范围介入并出现 }\alpha=0
$$

$$
\rightarrow
\text{高能尾部输运、加速和相混不足}
\rightarrow
\text{背景电子吸收的场功不足}
\rightarrow
\text{高能电子比例和平均动能偏低}
$$

$$
\rightarrow
\text{大量能量滞留于相干 }E_x
\rightarrow
\text{Beam 被进一步减速、俘获或再聚束}
\rightarrow
J_{\mathrm{beam}}\text{ 持续维持电场}
$$

$$
\rightarrow
\text{中点映射更不光滑，背景电流难以严格收敛}
\rightarrow
\text{软接受和相位误差继续累计}.
$$

与该动力学链同时发生的数值问题是：

$$
\text{FCT/最终接受通量的离散能量残差}
\rightarrow
R_E\text{ 在 }48-72\ \mathrm{fs}\text{ 快速下降}
\rightarrow
\text{累计能量缺口约 }1.46\times10^8\ \mathrm{J}.
$$

## 14. 各因素的重要性排序

### 14.1 主要问题

1. 高能和平行动量尾部发展严重不足；
2. 场能向背景电子动能的转换不足；
3. $48-72\ \mathrm{fs}$ 形成显著累计能量缺口；
4. 后期中点耦合长期软接受；
5. 过强相干场造成 Beam 滞留和再聚束。

### 14.2 次要但需要保留关注的问题

1. 周期接缝及端点场尖峰；
2. 边界背景密度异常；
3. 有限速度域和高能边界；
4. PIC 噪声与 Eulerian Vlasov 平滑结果之间的自然差异。

### 14.3 可以排除的单一解释

以下解释均不足以单独说明当前结果：

- 不能只用 EPOCH 的 PIC 噪声解释；
- 不能只用边界尖峰解释；
- 不能只用背景总粒子数漂移解释；
- 不能只用软接受解释；
- 不能只用 Beam 注入或出流误差解释。

## 15. 可信范围

### 15.1 $0-12\ \mathrm{fs}$

主波包幅度、位置、场能和背景加热均较接近 EPOCH，可认为定量可信。

### 15.2 $12-24\ \mathrm{fs}$

宏观幅值和能量仍具有较高参考价值，但不应要求局部相位逐点一致。

### 15.3 $24-36\ \mathrm{fs}$

电场幅度仍接近，但背景加热已偏低。该阶段可用于比较波包结构，但能量分配需谨慎。

### 15.4 $36-48\ \mathrm{fs}$

已经进入非线性分叉阶段，只能作定性比较。

### 15.5 $60-120\ \mathrm{fs}$

当前程序与 EPOCH 在场能、平均动能、高能尾部、Beam 残留和空间结构上均存在系统差异，不能视为 EPOCH 的定量等价结果。

## 16. 最终判断

本次生产长跑证明了当前程序可以稳定推进到 $120\ \mathrm{fs}$，并且早期结果明显优于过去多个版本。但长时间稳定推进并不等于后期物理结果已经正确。

最重要的结论是：

> 当前程序在约 $48\ \mathrm{fs}$ 后出现了可量化的高能尾部不足、背景加热不足、场能衰减不足和累计能量缺失。这些差异中包含明确的数值离散影响，不能全部归因于 EPOCH PIC 与当前 Vlasov 求解方法不同。

下一步问题定位应优先围绕以下三项展开：

1. 对 $36-72\ \mathrm{fs}$ 的最终已接受态逐步累计 FCT 能量残差；
2. 单独统计高能速度边界的粒子、动量和能量通量；
3. 区分严格收敛步与软接受步对 $J_{\mathrm{bkg}}$、场功和高能尾部的影响。

在完成这些审计前，不宜继续通过修改边界尖峰、强行匹配 EPOCH 波形或增加全局能量补丁来掩盖后期偏差。
