# 57 fs 波形偏差、时间步 A/B、FCT 与 face 配对回归综合分析

## 1. 文档状态与目标

本文档汇总截至当前已经完成的证据链，并以最新代码已经回退到稳定cell-capacity基线为前提：

1. `0--57 fs` 的 `dual-u + 0.6dt` 生产长跑及其与 EPOCH 的对比；
2. 45 fs checkpoint 的 `0.5dt/0.6dt` A/B；
3. 36 fs checkpoint 的 `40/80` 次中点迭代 A/B；
4. 24 fs、30 fs 和 36 fs checkpoint 的 `x-FCT/u-FCT` 分方向宏观加权审计；
5. `diagnostic-level=0/2` 的一步状态等价性检查。
6. 最终Yee面`J_N-G*J_E`直接配对审计；
7. “周期面最小范数伪逆 + 逐cell容量截断”的24/30/36 fs失败回归及其回退。

当前目标按优先级排序为：

1. 解释约36 fs后逐渐出现的核心区波形幅度、相位和包络偏差，并区分“已证实的离散缺陷”与“尚未证实的波形因果关系”；
2. 保持 3--27 fs 已经较好的物理结果；
3. 保持当前稳定的cell-capacity配对生产基线，不再用强制能量闭合直接改写物理状态；
4. 仅通过离线、受限、可回退的方式研究真实Yee面配对残差；
5. 在因果证据充分后再减少u-FCT的核心数值耗散；
6. 最后处理后期中点迭代变慢和软接受问题。

边界尖峰本身不是本文档的主要研究对象。只在边界误差进入核心区、能量账或全局收敛范数时处理。

本文档采用以下证据等级：

- **已证实事实**：由守恒恒等式、逐面账本或可复现A/B直接支持；
- **强相关假设**：空间、时间和宏观矩演化一致，但尚未通过只改变单一算子的A/B建立因果；
- **已否定路线**：实现后在24/30/36 fs短回归中明确恶化残差、动能、收敛或波形。

---

## 2. 数据来源

### 2.1 生产长跑

```text
output/dual_u_prod_dt060
results/优化后udal-u
```

该计算使用：

```text
background-coupling-mode = dual-u
dt = 0.6 * 原始时间步
midpoint-max-iters = 40
```

计算推进至约 57 fs。

### 2.2 EPOCH 参照

```text
noCollision_fields
noCollision_Beam_denstiy
noCollision_Bkgelectron_denstiy
noCollision_full_Bkgelectron_denstiy
noCollision_current_Jx
```

主要时间映射：

| VF/EPOCH编号 | 物理时间 |
|---|---:|
| `00005 / 110005` | 3 fs |
| `00015 / 110015` | 9 fs |
| `00030 / 110030` | 18 fs |
| `00045 / 110045` | 27 fs |
| `00060 / 110060` | 36 fs |
| `00065 / 110065` | 39 fs |
| `00070 / 110070` | 42 fs |
| `00075 / 110075` | 45 fs |
| `00090 / 110090` | 54 fs |
| `00095 / 110095` | 57 fs |

### 2.3 时间步与中点迭代 A/B

```text
output/checkpoint45_ab
output/wave_bias_36fs/A_dt050_iter40
output/wave_bias_36fs/B_dt050_iter80
```

### 2.4 FCT 宏观拆分审计

```text
output/wave_bias_36fs/C_fct_audit/24_C
output/wave_bias_36fs/C_fct_audit/30_C
output/wave_bias_36fs/C_fct_audit/36_C
```

审计文件：

```text
fct_macro_budget_by_region.dat
accepted_coupling_residual_by_region.dat
step_diagnostics.dat
scalars.dat
performance_summary.result
```

FCT 审计按以下维度拆分：

```text
flux_direction = x / u
x_region       = B_left / core / B_right
velocity_region= velocity_core / velocity_tail
```

### 2.5 最终面配对与失败修正回归

主要数据仍位于：

```text
output/wave_bias_36fs/C_fct_audit/24_C
output/wave_bias_36fs/C_fct_audit/30_C
output/wave_bias_36fs/C_fct_audit/36_C
```

重点文件：

```text
fixed_midpoint_face_pairing.dat
final_dual_u_pairing.dat
step_diagnostics.dat
scalars.dat
performance_summary.result
```

最新一轮曾将最终面残差通过周期`G*`伪逆投影为cell功电流修正，再按dual-u容量逐cell截断。该实现的24/30/36 fs短回归均失败，随后已经从生产代码回退。

因此本文档后文必须区分：

1. 回退前稳定cell-capacity基线；
2. 已失败且禁止恢复的周期面伪逆方案；
3. 尚未接入生产的带正则、容量和信赖域联合约束原型。

---

## 3. 与 EPOCH 的宏观对比

### 3.1 3--27 fs

这一阶段当前程序与 EPOCH 的以下量总体接近：

- 主波包位置；
- 主峰幅度；
- 波长和主相位；
- Beam 密度前沿和整体包络；
- 背景密度的主要相干响应。

因此，早期物理响应可以作为后续修改必须保留的基准。

任何修复如果明显改变 3--27 fs 的主峰、相位或包络，均不能仅凭能量账改善而接受。

### 3.2 27--36 fs

差异开始形成：

- VF 的波形比 EPOCH 更规则；
- `2.5--5.5 um` 的非规则调制开始减少；
- 背景密度同样变得更平滑；
- 主峰幅度和局部相位开始出现系统偏差。

36 fs 附近仍可认为定性接近，但已不能认为逐点一致。

### 3.3 39--42 fs

36 fs checkpoint 的 `0.5dt` 结果显示：

| 时间 | VF核心最大场 | EPOCH图像约值 |
|---|---:|---:|
| 36 fs | \(41.7E_0\) | 约 \(40--47E_0\) |
| 39 fs | \(48.9E_0\) | 约 \(45E_0\) |
| 42 fs | \(50.3E_0\) | 约 \(45E_0\) |

更重要的差异不是单个峰值，而是：

- VF 后续波列更加连续、平滑和规则；
- EPOCH 中存在更明显的非规则调制；
- VF 的背景密度也出现相同方向的平滑化。

EPOCH 是 PIC，部分细波来自粒子噪声，不能要求逐点重合。但主峰、包络宽度、相位和能量分配的系统偏差仍然具有物理意义。

### 3.4 45--57 fs

后期偏差已经宏观可见：

- 主波包峰值和相位持续漂移；
- 后续包络分布与 EPOCH 不同；
- 57 fs 时 VF 核心峰值约为 \(53.5E_0\)，而 EPOCH 图像约为 \(20--25E_0\)；
- 排除 `x<0.2 um` 和 `x>7.8 um` 后，偏差仍然存在。

因此，后期问题不是单纯的边界尖峰绘图污染。

上述对比能够证明宏观偏差存在，但不能单凭图像把偏差唯一归因于x-FCT、u-FCT或时间步。后文的算子审计用于定位离散缺陷，只有受控A/B才能进一步建立“某一缺陷导致波形偏差”的因果关系。

---

## 4. 已基本排除的来源

### 4.1 Beam 不是主要直接来源

Beam 密度和前沿在较长时间内仍与 EPOCH 接近，且 Beam 连续性、注入和出流账保持在较小误差范围。

因此不能把核心场偏差首先归因于：

- Beam 注入率错误；
- Beam 轨迹沉积再次失效；
- 右端开放出流失效。

Beam 与背景场的反馈仍可能放大误差，但不是当前证据指向的首要离散缺陷。

### 4.2 负分布爆炸不是本次后期偏差来源

接受态的 `f_negativity_monitor.dat` 显示：

```text
min_f = 0
neg_mass_total = 0
neg_cell_count = 0
```

因此当前偏差不是由接受态负分布突然爆炸引起。

### 4.3 40次中点上限不是36--42 fs偏差来源

36 fs checkpoint 的 A/B：

| 指标 | 40次上限 | 80次上限 |
|---|---:|---:|
| 接受步 | 8994 | 8994 |
| 严格接受 | 8994 | 8994 |
| 软接受 | 0 | 0 |
| 平均迭代数 | 13.0586 | 13.0586 |
| 采样最大迭代数 | 14 | 14 |

两组的：

- `fields_*.dat`；
- `density_*.dat`；
- `current_*.dat`；
- `scalars.dat`；
- `step_diagnostics.dat`；
- 背景最终状态哈希；
- 场最终状态哈希；

均一致。

核心区 `Ex`、`n_bkg` 和 `J_bkg_ampere` 的 A/B 相对 L2 与 Linf 差均为零。

结论：

> 36--42 fs 阶段的40次上限没有被触及。提高到80次不会修复波形，也没有生产价值。

### 4.4 直接周期面伪逆不是可行修复

24/30/36 fs回归已经排除以下路线：

```text
真实面残差
-> 去掉严格交替模
-> 对G*做最小范数伪逆
-> 得到全局cell修正
-> 按cell正性容量事后截断
```

该方案使真实面残差放大约6--7倍，并在0.2 fs内产生`1e9 J/m2`量级的异常背景动能变化。36 fs核心区`Ex/n_bkg/J_bkg`也发生数十个百分点的变化。

所以失败原因不是保护阈值太严，也不是中点迭代次数不足，而是算法结构本身错误：近Nyquist弱可表示模被伪逆放大，逐cell截断又破坏全局解的空间相关性。该方案已经回退，后续不得以“换一个伪逆容差”的方式恢复。

---

## 5. 时间步 A/B 的结论

### 5.1 45--51 fs 的 `0.5dt/0.6dt`

| 指标 | 原始0.5dt | 原始0.6dt |
|---|---:|---:|
| 平均迭代数 | 14.12 | 19.91 |
| 软接受比例 | 5.97% | 20.20% |
| 单位物理时间耗时 | 较低 | 约为0.5dt的1.51倍 |

尽管 `0.5dt` 步数更多，它仍然更快，因为迭代数和软接受显著减少。

51 fs 两组核心状态差：

```text
Ex relative L2    ≈ 0.72%
Ex relative Linf  ≈ 1.80%
n_bkg relative L2 ≈ 0.066%
J_bkg relative L2 ≈ 0.56%
```

两组新增能量误差只相差约2.3%。

### 5.2 时间步结论

`0.5dt` 比 `0.6dt`：

- 更容易严格收敛；
- 单位物理时间更快；
- 后期结果略稳定。

但45 fs checkpoint已经继承波形偏差，因此从45 fs缩小时间步无法恢复到 EPOCH。

当前生产如果必须继续，应优先使用原始 `0.5dt`。但时间步不是36 fs偏差的根本修复。

---

## 6. FCT覆盖率的正确解释

生产长跑中的原始总体覆盖率大致为：

| 时间段 | FCT active fraction |
|---|---:|
| 0--12 fs | 0.122 |
| 12--24 fs | 0.339 |
| 24--36 fs | 0.429 |
| 36--48 fs | 0.418 |
| 48--57 fs | 0.395 |

30 fs和36 fs的拆分短审计则给出：

| 指标 | 30 fs | 36 fs |
|---|---:|---:|
| 总 limiter 覆盖率 | 42.88% | 41.68% |
| x-limiter 覆盖率 | 16.31% | 15.23% |

因此：

1. FCT从早期到24--30 fs明显扩大；
2. 30--36 fs并没有继续扩大，反而略有下降；
3. “36 fs附近FCT覆盖率突然增加导致波形偏差”这一说法不成立；
4. FCT仍可能通过长期累计和宏观矩重分配造成偏差，但机制不是覆盖率在36 fs突增。

`performance_summary.result`中旧的`x_limiter_active_fraction_mean=0`与逐步诊断矛盾。新增审计已经给出可信的实际x覆盖率，后续不再使用旧的零值判断。

同时必须避免另一种过度解释：面数覆盖率没有质量、动量或能量权重。尾部大量低密度面被限制，并不等价于相同比例的物理核心被低阶化。判断u-FCT是否破坏波形，应优先使用质量/能量加权覆盖率和`delta_N/P/K`，不能只看active fraction。

---

## 7. FCT诊断验收

### 7.1 输出完整性

30 fs审计：

```text
299 accepted steps
0 soft accepted
30.0015--30.2003 fs
```

36 fs审计：

```text
300 accepted steps
0 soft accepted
36.0008--36.2003 fs
```

每个接受步均有：

```text
2 flux directions
× 3 x regions
× 2 velocity regions
= 12 rows
```

没有混入失败候选态。

### 7.2 守恒与独立闭合

`x-FCT`全局满足：

$
\sum\Delta N_x\simeq0,\qquad
\sum\Delta J_x\simeq0,\qquad
\sum\Delta K_x\simeq0.
$

`u-FCT`全局满足：

$
\sum\Delta N_u=0,
\qquad
\Delta K_u-E\cdot J_u\simeq0.
$

各分区`delta_K_fct`之和与独立的
`stage5_R_couple_fct_global`闭合：

```text
30 fs 最大误差 ≈ 9.88e-7 J/m2
36 fs 最大误差 ≈ 7.44e-7 J/m2
```

因此，方向拆分和能量分区统计可以用于根因判断。

### 7.3 diagnostic-level哈希验收

一步测试结果：

```text
background hash : identical
field-face hash : identical
Beam hash       : different
```

同时以下宏观文件逐字节一致：

```text
fields_00001.dat
fields_face_00001.dat
density_00001.dat
current_00001.dat
fv_bkg_e_00001.dat
```

重复运行中Beam hash本身也会变化，而Beam宏观量一致。这说明当前Beam hash不是稳定的、粒子排列无关的哈希。

结论：

- 详细FCT诊断对背景和场无副作用，验收通过；
- Beam哈希的字面验收未通过，但证据指向哈希方法不稳定，而非FCT诊断改变Beam物理状态；
- 后续应将Beam哈希改为按稳定粒子ID排序，或使用排列无关哈希；
- 在哈希修正前，Beam需同时比较总权重、动能、位置/动量矩、沉积密度和沉积电流。

---

## 8. x-FCT拆分结果

### 8.1 全局性质

`x-FCT`对全局质量、电流矩和动能几乎严格守恒，但产生非零局域场功：

$\Delta K_x^{\rm FCT}\simeq0,\qquad\Delta t\sum E_{\rm mid}J_x^{\rm FCT}\,dx\ne0.$

因此：
$R_{x,\rm FCT}=\Delta K_x^{\rm FCT}-\Delta t\sum E_{\rm mid}J_x^{\rm FCT}\,dx$
为持续非零。

| 审计窗口 | \(E\cdot J_x^{FCT}\) | \(R_{x,FCT}\) |
|---|---:|---:|
| 30.0--30.2 fs | \(-7.526\times10^5\) J/m² | \(+7.526\times10^5\) J/m² |
| 36.0--36.2 fs | \(-3.580\times10^5\) J/m² | \(+3.580\times10^5\) J/m² |

两组中`R_x`几乎每一步均为正。

### 8.2 空间分布

| 窗口 | 核心区 \(R_x\) | 两侧边界 \(R_x\) |
|---|---:|---:|
| 30 fs | \(+4.376\times10^5\) | \(+3.150\times10^5\) J/m² |
| 36 fs | \(+1.458\times10^5\) | \(+2.122\times10^5\) J/m² |

因此：

1. 30 fs时x-FCT残差的核心贡献大于边界；
2. 36 fs时边界净贡献更大；
3. x-FCT误差不能简单归结为边界；
4. x-FCT残差从30到36 fs下降约52%，不是36 fs附近新产生的误差。

### 8.3 物理解释

背景密度连续性要求Ampere使用来自x输运的电荷守恒面电流。

但是，x-FCT的数值修正电流会进入Ampere并产生场功，而周期空间输运本身不改变全局背景动能。如果速度受力功电流中没有同源的对应修正，就出现：

$J_N^{\rm FCT} \ne G^*J_E^{\rm FCT}.$

这不是“FCT是否保质量”的问题，而是：

> x-FCT连续性电流与u方向受力功电流没有属于同一个受限离散系统。

x-FCT是当前离散能量不配对的明确来源，这是已证实事实。但现有结果尚不能证明“将该残差压到零就会恢复EPOCH波形”。最新失败的face修正恰好说明：如果修正没有同时尊重`G*`频谱、dual-u容量和完整候选态，能量账可以被强行改变，而动力学会更差。

---

## 9. u-FCT拆分结果

### 9.1 能量配对

u-FCT满足：
$\Delta K_u^{\rm FCT}=\Delta t\sum E_{\rm mid}J_u^{\rm FCT}\,dx$


到约 $10^{-8}\,\mathrm{J/m^2}$ 的舍入误差。

因此：

> u-FCT不是当前全局能量残差的直接来源。

不能为了修能量账而破坏现有u-FCT功配对。

### 9.2 速度核心到尾部的重分配

核心空间区的累计修正：

| 0.2 fs窗口 | 速度核心 \(\Delta K_u\) | 速度尾部 \(\Delta K_u\) | 核心净值 |
|---|---:|---:|---:|
| 30 fs | \(-3.865\times10^6\) | \(+3.682\times10^6\) | \(-1.827\times10^5\) J/m² |
| 36 fs | \(-2.899\times10^6\) | \(+2.495\times10^6\) | \(-4.046\times10^5\) J/m² |

对应的粒子数修正表现为：

```text
velocity_core : 大量负修正
velocity_tail : 近等量正修正
```

这说明u-FCT持续把分布从速度核心向尾部重分配。

### 9.3 覆盖率

典型激活比例：

| 方向 | 区域 | 30 fs速度核心 | 30 fs尾部 | 36 fs速度核心 | 36 fs尾部 |
|---|---|---:|---:|---:|---:|
| x | core | 1.73% | 39.56% | 1.24% | 38.64% |
| u | core | 4.57% | 88.63% | 3.14% | 89.12% |

多数受限面位于速度尾部，但少量速度核心面具有较大的宏观权重。

### 9.4 物理影响

u-FCT即使能量闭合，也可能通过数值扩散：

- 降低速度核心中的相干响应；
- 增加尾部人口；
- 平滑背景密度；
- 降低或改变相干背景电流；
- 改变波包包络和相位。

36 fs时核心区u-FCT净能量损失比30 fs增大约2.2倍，说明其对核心动力学的偏置没有随覆盖率下降而消失。

这是一条强相关证据，不是已经完成的因果验证。因为尚未进行“只改变u-FCT预算、其他算子完全相同”的30--42 fs受控A/B，当前不能断言u-FCT就是36 fs后波形偏差的唯一主因。特别是未经加权的尾部高覆盖率不能直接代表核心物理被同等强度压制。

---

## 10. 能量账的更新判断

生产长跑累计能量误差约为：

| 时间 | 能量账误差占总能量 |
|---|---:|
| 12 fs | \(-0.008\%\) |
| 24 fs | \(+0.16\%\) |
| 36 fs | \(+0.88\%\) |
| 42 fs | 约 \(1\%\) |
| 54 fs | \(-2.56\%\) |
| 57 fs | \(-4.50\%\) |

短窗口结果：

| 窗口 | 段内能量账误差 | x-FCT残差 |
|---|---:|---:|
| 30.0--30.2 fs | \(+5.420\times10^5\) | \(+7.526\times10^5\) J/m² |
| 36.0--36.2 fs | \(-3.720\times10^5\) | \(+3.580\times10^5\) J/m² |

30 fs时x-FCT残差与段内能量误差同量级，说明它是重要来源。

36 fs时总能量误差变为负，而x-FCT残差仍为正，说明该阶段还有其他反号贡献。不能用单一x-FCT项解释总账符号，但x-FCT的不配对事实已经成立。

更新后的判断是：

1. x-FCT产生持续同号的局域场功不配对；
2. u-FCT能量闭合，但会重塑速度分布；
3. 完整生产能量误差还包含中点截断误差、Beam出流后的累计账和其他离散项；
4. 不能通过加入一个全局能量补偿项掩盖x/u不配对；
5. 能量残差是必要诊断量，但不是可以脱离状态变化单独优化的目标函数。

失败的周期面伪逆回归提供了决定性的反例：即使目标是减小`J_N-G*J_E`，不受控的功电流修正仍会在0.2 fs内造成约`1.4e9--2.4e9 J/m2`的异常背景动能增长，并显著改变`Ex/n_bkg/J_bkg`。因此后续验收必须同时检查真实面残差、完整背景动能、场功、中点收敛和宏观波形，不能只检查某一列能量账是否更接近零。

---

## 11. 中点收敛问题的正确位置

36--42 fs：

- 所有步骤严格接受；
- 最大实际迭代数约14；
- 40/80次A/B完全一致。

因此，这一阶段波形偏差不是软接受造成。

48 fs以后：

- 平均迭代数明显增加；
- `0.6dt`软接受比例上升；
- 每步可能执行满40次昂贵算子；
- 性能和离散误差进一步恶化。

所以中点问题是：

1. 后期性能恶化的直接原因；
2. 48 fs以后误差的放大器；
3. 不是30--42 fs波形偏差的起点。

修复顺序不能颠倒，但也不能再直接修改生产算子。应先用离线固定状态原型判断x/u面残差在容量和频谱约束下究竟有多少可安全降低；通过后才做生产A/B。u-FCT核心耗散必须单独测试。最后再优化后期非线性求解器。

---

## 12. 当前根因链

现有证据支持两条并行假设链和一条已证实的后期放大链，但证据强度不同。

### 12.1 已证实的离散能量不配对链

$
\text{x高阶通量被FCT修正}
\rightarrow
J_N^{\rm FCT}\text{进入Ampere}
\rightarrow
E\cdot J_N^{\rm FCT}\ne0
\rightarrow
\text{u受力功中缺少同源对应}
\rightarrow
\text{能量误差长期累计}.
$

其中“x-FCT产生非零配对残差”已经由逐面账本证明；“该累计是36 fs后波形偏差的主要原因”仍是强相关假设，尚需受限A/B验证。

### 12.2 强相关但尚未完成因果验证的动力学平滑链

$
\text{u-FCT长期高尾部覆盖}
\rightarrow
\text{速度核心向尾部重分配}
\rightarrow
\text{相干背景响应被平滑}
\rightarrow
n_{\rm bkg},J_{\rm bkg},E_x
\text{包络和相位逐渐偏离}.
$

目前已证明u-FCT在速度核心和尾部间进行宏观可见的质量/能量重分配，也证明背景密度和场相对EPOCH更平滑；尚未证明只改变u-FCT预算就能恢复波形。

### 12.3 后期放大链

$
\text{非线性结构增强}
\rightarrow
\text{FCT active set和固定点映射不光滑}
\rightarrow
\text{迭代数增加}
\rightarrow
\text{软接受增多}
\rightarrow
\text{48 fs后性能和波形进一步恶化}.
$

该链在48 fs后由迭代数、软接受比例和耗时共同支持，但不能反推30--42 fs的偏差也由软接受产生。

### 12.4 已否定的“强制face闭合”链

$
\text{面残差包含近Nyquist弱可表示模}
\rightarrow
\text{无约束伪逆产生大幅振荡cell修正}
\rightarrow
\text{逐cell容量截断破坏空间相关性}
\rightarrow
\text{真实面残差反而放大}
\rightarrow
\Delta K_{\rm bkg},E_x,n_{\rm bkg},J_{\rm bkg}
\text{迅速失真}.
$

这条链已经由24/30/36 fs短回归直接证实。因此下一步不能以“更严格闭合”为目标，而应以“在容量、频谱和信赖域内，真实残差单调下降且物理状态不恶化”为目标。

需要明确：

- FCT不是在36 fs突然开始工作；
- x-FCT残差30 fs已经很强，36 fs反而下降；
- 后期偏差可能是早期累计误差在非线性阶段的延迟显现；
- 不能仅凭覆盖率判断宏观影响；
- 不能仅凭能量残差变小判断修复正确；
- 当前证据不支持把某一个算子缺陷宣布为全部波形偏差的唯一根因。

---

## 13. 下一步工作

### 13.0 当前状态与唯一基线

当前代码已经回退到face伪逆修改之前的稳定cell-capacity实现。后续所有A/B的基线必须是该版本，而不是失败实验留下的24/30/36 fs输出。

状态划分如下：

| 项目 | 状态 | 是否可作为生产基线 |
|---|---|---|
| 24/30/36 fs FCT宏观拆分和逐面账本 | 已完成、可信 | 是，作为诊断基线 |
| cell平均目标 + `apply_final_limited_capacity_pairing()` | 已恢复、稳定 | 是 |
| 周期面无约束伪逆 + 逐cell截断 | 已失败并回退 | 否 |
| 带正则、容量和信赖域联合约束 | 尚未实现 | 否，只能先做离线原型 |
| u-FCT核心预算优化 | 尚未做因果A/B | 否，排在13.3之后 |

如果现有`output/wave_bias_36fs/C_fct_audit/24_C`、`30_C`、`36_C`已经被失败实验覆盖，必须保留失败输出用于反例分析，同时从相同checkpoint用当前回退代码重新生成一份明确命名的`stable_cell_baseline`。两类输出不能混写或覆盖。

### 13.1 第一步：稳定基线的24 fs短审计（历史结果已完成；必要时重建输出）

目的：

- 确定x-FCT不配对在波形仍较好时是否已经存在；
- 判断30 fs的大残差是长期存在还是24--30 fs开始增强；
- 不需要从24 fs连续跑到30 fs。

只从24 fs checkpoint推进0.2 fs。

提交脚本：

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J fct24audit
#SBATCH -o fct24audit_%j.out
#SBATCH -e fct24audit_%j.err

set -euo pipefail

cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

CHECKPOINT_24=/完整路径/checkpoints/t_24fs

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4

yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --midpoint-trace-window-fs 24.0,24.2 \
  --stop-time-fs 24.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/wave_bias_36fs/stable_cell_baseline/24_C \
  --overwrite-output
```

验收指标与30/36 fs完全相同：

```text
x/u分方向的 delta_N, delta_J, delta_K, E_dot_J, R_fct_E
空间边界/核心
速度核心/尾部
limiter覆盖率
段内dKE_bkg和能量账
```

#### 13.1.1 运行有效性

以下数值属于face伪逆实验之前的稳定cell-capacity基线。24 fs checkpoint实际推进至约24.2 fs，共接受299步：

- `strict_accepted_steps=299`；
- `soft_accepted_steps=0`；
- 平均中点迭代数为13；
- 最终场残差为 $6.24\times10^{-7}$；
- 最终背景电流残差为 $5.19\times10^{-6}$。

因此本段结果不是由软接受或未收敛中点状态造成，能够用于与30/36 fs审计进行同口径比较。

生产性能口径下：

- 总limiter平均覆盖率为43.74%；
- x-limiter平均覆盖率为17.31%；
- 两者的最小alpha均达到0。

#### 13.1.2 x-FCT结果

24--24.2 fs累计：

| 区域 | $R_{x,\mathrm{FCT}}$ / J/m2 |
|---|---:|
| 左边界 | $+1.226\times10^5$ |
| 核心区 | $+3.627\times10^5$ |
| 右边界 | $+1.175\times10^5$ |
| 全域 | $+6.029\times10^5$ |

其中核心区占全域x-FCT残差约60.2%，不能把该项简单归结为边界异常。

全域宏观账为：

$
\Delta K_{x,\mathrm{FCT}}\simeq-3.97\times10^{-7}\ {\rm J/m^2},
$

$
W_{x,\mathrm{FCT}}=-6.029\times10^5\ {\rm J/m^2},
$

$
R_{x,\mathrm{FCT}}
=\Delta K_{x,\mathrm{FCT}}-W_{x,\mathrm{FCT}}
=+6.029\times10^5\ {\rm J/m^2}.
$

x-FCT对全局粒子数、动量矩和动能几乎守恒，但其连续性电流进入Ampere后产生了没有同源背景动能变化对应的场功。该不配对在24 fs、即波形仍相对合理时已经存在。

#### 13.1.3 u-FCT结果

24--24.2 fs累计：

$
\Delta K_{u,\mathrm{FCT}}
=-4.01126566\times10^5\ {\rm J/m^2},
$

$
W_{u,\mathrm{FCT}}
=-4.01126566\times10^5\ {\rm J/m^2},
$

$
R_{u,\mathrm{FCT}}\simeq1.0\times10^{-8}\ {\rm J/m^2}.
$

因此u-FCT的全局功配对保持到舍入误差。但速度分区显示：

| 速度区域 | $\Delta K_{u,\mathrm{FCT}}$ / J/m2 |
|---|---:|
| 速度核心 | $-4.077\times10^6$ |
| 速度尾部 | $+3.676\times10^6$ |
| 全域净值 | $-4.011\times10^5$ |

这再次证明“全局能量闭合”不等于“没有动力学平滑”：u-FCT仍在大量重分配速度核心和尾部能量。

#### 13.1.4 24/30/36 fs时间演化

| 起始时间 | x-FCT残差 / J/m2 | 直接x/u-FCT配对残差 / J/m2 |
|---|---:|---:|
| 24 fs | $+6.029\times10^5$ | $-2.017\times10^5$ |
| 30 fs | $+7.526\times10^5$ | $-6.513\times10^5$ |
| 36 fs | $+3.580\times10^5$ | $-1.985\times10^5$ |

据此可得：

1. x/u不配对不是30 fs才出现，而是在24 fs已经存在；
2. 该误差在30 fs附近明显增强；
3. 36 fs时瞬时残差反而下降，不能把36 fs后波形偏差解释为该时刻FCT突然增强；
4. 更合理的解释是历史累计、30 fs附近的增强、u-FCT速度空间重分配和后续非线性相位演化共同作用。

### 13.2 第二步：直接FCT配对残差（已完成）

`fixed_midpoint_face_pairing.dat`已经直接输出：

$J_{N,\rm FCT},\qquad G^*J_{E,\rm FCT},\qquad J_{N,\rm FCT}-G^*J_{E,\rm FCT}.$

该账本满足原定审计要求：

1. 使用最终接受态实际采用的x/u FCT面通量；
2. 不能重放算子或重新近似构造；
3. 按`B_left/core/B_right`输出；
4. 同时逐面输出并汇总：

$\Delta t\sum E_{\rm mid}\left(J_{N,\rm FCT}-G^*J_{E,\rm FCT}\right)dx.$

#### 13.2.1 24 fs直接配对结果

24--24.2 fs累计：

$
\Delta R_{x,\mathrm{FCT}}
=-6.02872\times10^5\ {\rm J/m^2},
$

$
\Delta R_{u,\mathrm{FCT}}
=+4.01127\times10^5\ {\rm J/m^2},
$

$
\Delta R_{\mathrm{FCT,pair}}
=\Delta R_{x,\mathrm{FCT}}+\Delta R_{u,\mathrm{FCT}}
=-2.01745\times10^5\ {\rm J/m^2}.
$

这里x项与13.1宏观预算中的 $R_{x,\mathrm{FCT}}$ 大小一致、符号相反，来自两份诊断对残差方向的定义不同；数值大小的一致完成了直接逐面账本与宏观FCT预算的交叉验证。

直接配对残差的空间分布为：

| 区域 | $\Delta R_{\mathrm{FCT,pair}}$ / J/m2 |
|---|---:|
| 左边界 | $-1.457\times10^5$ |
| 核心区 | $+8.691\times10^4$ |
| 右边界 | $-1.429\times10^5$ |
| 全域 | $-2.017\times10^5$ |

核心与两侧边界存在明显反号抵消。因此只看全局净残差会低估局部配对误差，也不能把全部误差仅归因于边界。

加入低阶、高阶中心/迎风和FCT全部分量后，最终完整配对残差累计为：

$
R_{\mathrm{pair,final}}=-8.960\times10^4\ {\rm J/m^2}.
$

其空间分布为：

| 区域 | $R_{\mathrm{pair,final}}$ / J/m2 |
|---|---:|
| 左边界 | $-7.256\times10^4$ |
| 核心区 | $+1.052\times10^5$ |
| 右边界 | $-1.223\times10^5$ |

这说明FCT配对残差还会被其他离散分量部分抵消，但最终受限系统仍没有严格闭合。

#### 13.2.2 账本可信度

- 299个接受步全部覆盖4000个唯一x面；
- `R_pair_reconstructed_final-R_pair_final`累计仅约 $6.4\times10^{-10}\ {\rm J/m^2}$；
- 最大逐面重构误差约 $4.3\times10^{-11}\ {\rm J/m^2}$；
- 其中293步`ledger_match=1`；
- 其余6步的独立求和差仅为约 $1.9\times10^{-11}$ 至 $1.06\times10^{-10}\ {\rm J/m^2}$，只是超过了当前按机器精度设置的极严诊断阈值，不是物理账本失配。

因此13.2已经完成，不需要继续扩大通用审计范围。后续只需将这6次诊断假阴性改为稳定求和或与实际累计规模匹配的舍入容差，不能据此修改生产物理算子。

### 13.3 第三步：重做x/u受限离散配对，但禁止直接周期面伪逆

#### 13.3.1 最新回归对原方案的否定

已经实际尝试过以下路线：

1. 在周期唯一面上构造
   \[
   r_f=J_{N,f}-G^*J_{E,f};
   \]
2. 去掉严格不可表示的周期交替模；
3. 对剩余面残差做最小范数伪逆，得到cell功电流修正；
4. 再按每个cell可用的dual-u正性容量独立截断。

该实现已经回退，不能继续通过调小阈值或增加迭代复用。24/30/36 fs短回归给出的结果是：

- 修正后/修正前面残差中位数分别约为`6.93`、`7.28`和`6.17`，即残差不是下降，而是被放大约6--7倍；
- 24、30、36 fs中真正改善的步数分别仅为`0/299`、`5/299`和`0/300`；
- 可表示残差仍约为`1e19`量级，严格交替模仅约为`1e15`量级，说明失败不由单个严格零空间模主导；
- 每步受容量限制的cell平均约为813、670和566个，最小缩放中位数仅约`0.0047`、`0.0212`和`0.0208`；
- 仅推进0.2 fs，背景动能就异常增加约`2.4e9`、`1.6e9`和`1.4e9 J/m2`；
- 36 fs相对回退前稳定基线，核心区`Ex`、`n_bkg`和`J_bkg`的相对L2差已分别约为46%、34%和79%。

根因是：周期`G*`除严格交替零模外还存在一组近Nyquist弱可表示模。直接伪逆会用很大的cell修正追逐这些面模态；随后对各cell独立做容量裁剪，又破坏了伪逆解赖以成立的空间相关性。因此“无约束全局伪逆 + 逐cell事后截断”在数学上和物理上都不可接受。

#### 13.3.2 当前生产基线

在新方案通过全部离线验收前，生产程序必须保留已恢复的稳定实现：

```cpp
dual_target_jn_cell[ix] =
    0.5 * (dual_jn_high_face[ix] + dual_jn_high_face[ix + 1]);
```

并继续调用现有`apply_final_limited_capacity_pairing()`完成cell级容量受限配对。

该基线的局限必须明确：

- 它能稳定地完成cell标量功矩配对；
- 它保持现有x连续性电流、Yee面Ampere、u边界零通量和正性容量；
- 它没有证明逐个Yee面满足`J_N=G*J_E`；
- `fixed_midpoint_face_pairing.dat`中的真实面残差继续作为审计量，不能宣称已经实现严格face闭合。

这是一条稳定基线，不是最终数学闭合方案。当前不得为了减小能量账再次直接修改生产场或背景分布。

#### 13.3.3 最佳修复形式：带容量和频谱正则的联合约束问题

下一版不能先求无约束解再裁剪，而应把面残差、平滑正则、修正幅度和dual-u容量放进同一个问题：

$
\min_{\delta J_E}
\frac12\left\|W_f^{1/2}
\left(G^*\delta J_E-r_f^{\rm rep}\right)\right\|_2^2
+
\frac{\lambda}{2}\left\|D_x\delta J_E\right\|_2^2
+
\frac{\eta}{2}\left\|\delta J_E\right\|_2^2 ,
$

约束为：

$
\delta J_{E,i}\in
[\delta J_{E,i}^{\min},\delta J_{E,i}^{\max}],
$

其中上下界必须直接来自该cell的dual-u可用正、负通量容量，并同时满足：

1. 每个x cell内u修正质量矩严格为零；
2. u速度边界通量保持为零；
3. 共享u面只使用一个修正值；
4. 修正后的候选分布满足现有正性预算；
5. 不修改`J_N`，Ampere仍使用真实x输运连续性电流；
6. 不修改场零模，不做Poisson或全局能量回填。

这里的`r_f_rep`不能只减去一个严格交替模。应对`G*`做一次离线谱审计，将奇异值低于

$
\sigma_k < \tau_\sigma \sigma_{\max}
$

的严格零模和近零模统一列入`r_f_unresolved`。第一版建议从`tau_sigma=1e-6,1e-8,1e-10`三档做敏感性测试，而不是根据EPOCH波形拟合。未解析残差必须原样记账：

$
r_f=r_f^{\rm rep}+r_f^{\rm unresolved}.
$

求解器应使用带边界约束的投影CG、LSQR或active-set方法；不允许显式形成稠密伪逆。`lambda`和`eta`的作用分别是抑制近Nyquist振荡修正和限制修正幅度，不能把目标设为机器精度意义上的残差归零。

#### 13.3.4 信赖域与单调接受

即使联合约束问题收敛，也只能在同时满足以下条件时接受修正：

$
\|r_f^{\rm after}\|_{W_f}
<
\|r_f^{\rm before}\|_{W_f},
$

并且：

- 核心区和全域面残差都下降，不能依靠边界/核心反号抵消；
- `max(|delta_JE|/J_scale)`不超过预先给定的信赖域；
- 该步`|delta_KE_bkg|`和`|delta_W_bkg|`不超过未修正基线对应量的预设小比例；
- `candidate_min_f`、背景质量误差、中点残差和迭代次数不恶化；
- 重新由修正后的最终u通量计算`J_E`，不得直接覆盖诊断电流数组。

第一版信赖域应保守，例如只允许修正当前可用dual-u容量的10%，再依次测试25%和50%。若任何一项失败，该步必须无条件回退到稳定cell-capacity基线，并记录：

```text
face_pairing_attempted
face_pairing_accepted
face_residual_before
face_residual_after
unresolved_mode_l2
capacity_active_cells
trust_region_active_cells
delta_ke
delta_work
fallback_to_cell_baseline
```

回退必须恢复完整候选态，不能只恢复`J_E`而保留已经修改的u通量。

#### 13.3.5 实施与验收顺序

该修复必须分阶段进行：

1. **独立算子测试**：在周期均匀网格上构造常数、低频、接近Nyquist和严格交替面模态，验证谱分解、伴随关系、容量约束和回退逻辑。
2. **固定checkpoint离线原型**：从24/30/36 fs最终中点状态读取真实`r_f`和容量，不推进时间，只求一次联合约束修正。
3. **单步A/B**：基线与新方案各推进一步，要求真实Yee面残差下降且完整`f/E/J`变化受信赖域约束。
4. **0.2 fs短回归**：分别从24/30/36 fs运行。任何一个时间点出现面残差放大、动能突增、软接受增加或波形L2大幅变化，立即判废。
5. **30--42 fs物理回归**：只有前三个checkpoint均通过，才比较EPOCH的主峰、相位和包络。
6. **生产接入**：默认仍关闭新face修正；完成全部回归后才允许通过显式参数启用，最后再考虑设为默认。

通过目标不是把`J_N-G*J_E`压到舍入误差，而是在不改变已有合理波形的前提下稳定降低其可表示部分。若带容量约束后只能获得很小改善，应接受并记录剩余残差，不能继续提高修正强度。

#### 13.3.6 明确禁止

- 恢复已失败的周期面最小范数伪逆；
- 先无约束求解、后逐cell独立裁剪；
- 直接扣除电流零模或交替模；
- 每步Poisson投影；
- 强制全局能量回填；
- 从Ampere删除x-FCT电流；
- 为了减小能量残差直接覆盖`J_bkg`；
- 用EPOCH波形反向拟合`lambda`、`eta`或信赖域。

#### 13.3.7 各阶段编译与运行命令

以下命令分为两类：

- 标记为“当前可运行”的命令只使用当前已有参数；
- 标记为“实现后运行”的命令依赖13.3新增接口。接口尚未实现前，不能把示例参数直接传给`fp_solver`。

所有`yhrun`命令均应放在以下5节点提交脚本环境中，不能直接在登录节点运行生产求解器：

```bash
#!/bin/bash
#SBATCH -p xyfree
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -c 4
#SBATCH -J pairing_ab
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
```

##### A. 必须新增的13.3接口

建议新增两个独立测试目标：

```text
face_pairing_regularized_unit_test
face_pairing_checkpoint_audit
```

并给`fp_solver`新增：

```text
--face-pairing-mode cell-baseline|regularized
--face-pairing-sigma-cutoff <value>
--face-pairing-lambda <value>
--face-pairing-eta <value>
--face-pairing-trust-fraction <value>
```

要求：

- 默认必须是`cell-baseline`；
- `regularized`只有显式指定时启用；
- 所有失败都完整回退到`cell-baseline`；
- 输出目录中必须记录全部face参数和回退次数。

`face_pairing_checkpoint_audit`至少支持：

```text
--operator-audit <final_midpoint目录>
--sigma-cutoff <value>
--lambda <value>
--eta <value>
--trust-fraction <value>
--output-dir <目录>
```

##### B. 13.3统一编译（实现后运行）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

rm -rf build
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DFP_BUILD_BACKGROUND_COUPLING_TESTS=ON \
  -DFP_BUILD_PERIODIC_STAGGERED_TESTS=ON

cmake --build build -j4 --target \
  fp_solver \
  periodic_staggered_adjoint_test \
  final_dual_u_limited_pairing_test \
  face_pairing_regularized_unit_test \
  face_pairing_checkpoint_audit
```

如果`face_pairing_regularized_unit_test`或`face_pairing_checkpoint_audit`目标不存在，说明13.3独立原型尚未实现，不能进入后续生产A/B。

##### C. 门A：独立算子测试

当前已有测试：

```bash
yhrun -N 1 -n 16 --cpu-bind=cores ./build/periodic_staggered_adjoint_test \
  > output/face_pairing_gateA_periodic_adjoint.result 2>&1

./build/final_dual_u_limited_pairing_test \
  > output/face_pairing_gateA_capacity.result 2>&1
```

实现后新增测试：

```bash
./build/face_pairing_regularized_unit_test \
  > output/face_pairing_gateA_regularized.result 2>&1
```

三个文件均显示`PASS`后才能进入门B。新增测试必须分别报告常数、低频、近Nyquist、严格交替模、容量饱和和回退用例，不能只给一个总`PASS`。

##### D. 生成24/30/36 fs固定中点状态（当前可运行）

先在提交脚本中设置真实checkpoint路径：

```bash
CHECKPOINT_24=/完整路径/checkpoints/t_24fs
CHECKPOINT_30=/完整路径/checkpoints/t_30fs
CHECKPOINT_36=/完整路径/checkpoints/t_36fs
```

以下命令使用`--dt-scale 0.8333333333333333`的前提是：24/30/36 fs checkpoint均由`0.6dt`长跑生成，需要转换到原始`0.5dt`。若checkpoint保存时已经使用目标时间步，则改为：

```text
--dt-scale 1.0
```

三组A/B必须使用相同checkpoint和相同`dt-scale`。

24 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --dump-final-midpoint \
  --output-dir ./output/face_pairing_gateB/24_midpoint \
  --overwrite-output
```

30 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --dump-final-midpoint \
  --output-dir ./output/face_pairing_gateB/30_midpoint \
  --overwrite-output
```

36 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_36" \
  --background-coupling-mode dual-u \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --dump-final-midpoint \
  --output-dir ./output/face_pairing_gateB/36_midpoint \
  --overwrite-output
```

必须确认以下目录分别存在且manifest中的MPI规模、`nx/Nu/Nuperp`一致：

```text
output/face_pairing_gateB/24_midpoint/final_midpoint
output/face_pairing_gateB/30_midpoint/final_midpoint
output/face_pairing_gateB/36_midpoint/final_midpoint
```

##### E. 门B：固定checkpoint离线联合约束审计（实现后运行）

先只测试最保守参数：

```text
sigma-cutoff = 1e-6
trust-fraction = 0.10
```

24 fs：

```bash
yhrun --cpu-bind=cores ./build/face_pairing_checkpoint_audit \
  --operator-audit ./output/face_pairing_gateB/24_midpoint/final_midpoint \
  --sigma-cutoff 1e-6 \
  --lambda 1e-2 \
  --eta 1e-4 \
  --trust-fraction 0.10 \
  --output-dir ./output/face_pairing_gateB/24_regularized
```

30 fs：

```bash
yhrun --cpu-bind=cores ./build/face_pairing_checkpoint_audit \
  --operator-audit ./output/face_pairing_gateB/30_midpoint/final_midpoint \
  --sigma-cutoff 1e-6 \
  --lambda 1e-2 \
  --eta 1e-4 \
  --trust-fraction 0.10 \
  --output-dir ./output/face_pairing_gateB/30_regularized
```

36 fs：

```bash
yhrun --cpu-bind=cores ./build/face_pairing_checkpoint_audit \
  --operator-audit ./output/face_pairing_gateB/36_midpoint/final_midpoint \
  --sigma-cutoff 1e-6 \
  --lambda 1e-2 \
  --eta 1e-4 \
  --trust-fraction 0.10 \
  --output-dir ./output/face_pairing_gateB/36_regularized
```

只有三组真实面残差均下降且没有动能异常，才依次测试：

```text
sigma-cutoff = 1e-8, 1e-10
trust-fraction = 0.25, 0.50
```

不能一次批量提交全部参数。每扩大一档前先分析上一档；一旦残差、动能或容量占用恶化，停止扩大。

##### F. 门C1：24/30/36 fs单步baseline与regularized A/B（实现后运行）

以下以24 fs为例，30/36 fs只替换checkpoint和输出目录。

baseline：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --dump-final-midpoint \
  --output-dir ./output/face_pairing_gateC/24_baseline_1step \
  --overwrite-output
```

regularized：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --face-pairing-mode regularized \
  --face-pairing-sigma-cutoff 1e-6 \
  --face-pairing-lambda 1e-2 \
  --face-pairing-eta 1e-4 \
  --face-pairing-trust-fraction 0.10 \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-after-steps 1 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --dump-final-midpoint \
  --output-dir ./output/face_pairing_gateC/24_regularized_1step \
  --overwrite-output
```

30 fs和36 fs必须分别独立重复，不能用24 fs通过替代。

##### G. 门C2：24/30/36 fs各推进0.2 fs（实现后运行）

24 fs baseline：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 24.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/face_pairing_gateC/24_baseline_02fs \
  --overwrite-output
```

24 fs regularized：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --face-pairing-mode regularized \
  --face-pairing-sigma-cutoff 1e-6 \
  --face-pairing-lambda 1e-2 \
  --face-pairing-eta 1e-4 \
  --face-pairing-trust-fraction 0.10 \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 24.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/face_pairing_gateC/24_regularized_02fs \
  --overwrite-output
```

30 fs对应使用：

```text
CHECKPOINT_30
--stop-time-fs 30.2
output/face_pairing_gateC/30_baseline_02fs
output/face_pairing_gateC/30_regularized_02fs
```

36 fs对应使用：

```text
CHECKPOINT_36
--stop-time-fs 36.2
output/face_pairing_gateC/36_baseline_02fs
output/face_pairing_gateC/36_regularized_02fs
```

门C必须比较：

```text
fixed_midpoint_face_pairing.dat
final_dual_u_pairing.dat
step_diagnostics.dat
scalars.dat
fields_*.dat
density_*.dat
current_*.dat
performance_summary.result
```

##### H. 门D：30--42 fs物理回归（仅门A--C全部通过后）

baseline：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --checkpoint-times 36,42 \
  --stop-time-fs 42 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir ./output/face_pairing_gateD/30_to_42_baseline \
  --overwrite-output
```

regularized：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode regularized \
  --face-pairing-sigma-cutoff 1e-6 \
  --face-pairing-lambda 1e-2 \
  --face-pairing-eta 1e-4 \
  --face-pairing-trust-fraction 0.10 \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --checkpoint-times 36,42 \
  --stop-time-fs 42 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir ./output/face_pairing_gateD/30_to_42_regularized \
  --overwrite-output
```

门D不通过时必须保留`cell-baseline`为生产默认，并停止13.3；不能继续增大trust fraction追求残差。

### 13.4 第四步：在13.3稳定后再降低u-FCT核心耗散

u-FCT的全局功配对目前已达到舍入级，因此13.4不是优先于13.3的紧急修复。它的目标也不是减少`limiter_active_fraction`这个未加权数字，而是减少物理速度核心中不必要的低阶退化，同时保留防止负分布失控的能力。

#### 13.4.1 先补齐不改变生产推进的诊断

将u-limiter统计拆成速度核心、过渡区和尾部，并同时输出：

- 面数加权`alpha`直方图；
- 正质量加权和动能加权的`alpha`直方图；
- `alpha=0`面对应的正质量、动能和电流权重；
- 各区间`delta_N_u`、`delta_P_u`、`delta_K_u`；
- 低阶通量、高阶反扩散通量和最终通量各自的能量矩；
- 触发限制的原因：donor正性预算、receiver上界预算、dual-u容量或边界通量。

只有当质量/能量加权结果确认核心区存在宏观显著退化时，才进入生产修改。大量低密度尾部面被限制，不能单独作为改算法的依据。

#### 13.4.2 只修改高阶反扩散预算

保持以下内容不变：

1. 当前保正低阶u通量；
2. 当前高阶候选通量公式；
3. u边界零通量；
4. 每个共享u面唯一的`alpha_face`；
5. 最终u通量用于动能矩和`J_E`的同源关系。

需要审计并修改的是`F_high-F_low`的反扩散预算：

- donor和receiver预算必须来自同一个最终中点候选态；
- 预算只在共享面的两个相邻cell间分配，避免过宽邻域中的极小尾部值把速度核心一起压到低阶；
- 对同一cell的所有入流/出流反扩散通量联合分配预算，不能逐面独立耗尽同一份正质量；
- `alpha_face`取面两侧约束的最小值，更新两侧时使用同一个值；
- dual-u功矩修正必须使用FCT结束后的剩余容量，不能与FCT分别消费同一份容量；
- 不允许通过事后`max(f,0)`、负值裁剪或质量回填补救。

速度核心和尾部不能使用不同物理通量公式。允许不同的只是诊断分区和容量分配优先级：在总正性约束不变时，优先保留质量/能量权重显著的核心反扩散通量，尾部在必要时更早退化。

#### 13.4.3 修改必须是单调、小步和可回退的

每次只改变一个因素，建议顺序为：

1. 修正过宽邻域预算；
2. 修正同一cell预算被重复消费；
3. 将dual-u和FCT容量统一结算；
4. 最后才测试核心优先的容量分配。

每一步都保留旧实现作为A/B模式。若出现以下任一情况，立即回退：

- u-FCT功配对残差不再是舍入级；
- 背景质量误差增大；
- `min_f`或负质量出现持续增长；
- x方向配对残差增大；
- 中点平均迭代数或软接受比例上升；
- 3--27 fs波形幅度被压低，或24/30/36 fs短回归与稳定基线出现显著宏观差异。

#### 13.4.4 验收顺序

1. 单cell和两cell共享u面的质量、正性、能量矩单元测试；
2. 非均匀u网格与正负电场测试；
3. 固定24/30/36 fs中点状态的limiter重放；
4. 三个checkpoint分别推进0.2 fs；
5. 比较质量/能量加权limiter覆盖率，而不是只看原始面数覆盖率；
6. 30--42 fs包络、相位和能量账回归；
7. 只有核心耗散下降、面配对不恶化且波形改善，才进入0--57 fs长回归。

13.3与13.4不能同时大改。应先完成13.3的受限原型判定，再以稳定版本为唯一基线执行13.4，否则无法区分面配对修正与FCT预算变化各自造成的影响。

#### 13.4.5 各阶段编译与运行命令

##### A. 必须新增的13.4接口

建议新增：

```text
u_fct_budget_unit_test
u_fct_checkpoint_replay
```

并给`fp_solver`新增：

```text
--u-fct-budget-mode baseline|local-joint|unified-capacity|core-priority
--u-fct-weighted-audit
```

约束：

- 默认必须为`baseline`；
- `--u-fct-weighted-audit`只增加已接受态诊断，不改变推进；
- `local-joint`只修正过宽邻域和同一cell预算重复消费；
- `unified-capacity`在`local-joint`基础上统一FCT与dual-u剩余容量；
- `core-priority`最后测试，不能先于前两者启用；
- 每个模式必须写入结果文件，防止输出来源混淆。

`u_fct_checkpoint_replay`至少支持：

```text
--operator-audit <final_midpoint目录>
--u-fct-budget-mode <mode>
--output-dir <目录>
```

##### B. 13.4编译（实现后运行）

```bash
cd /HOME/sztu_lbju/sztu_lbju_1/HDD_POOL/Chendong/FokkerPlanckEquationSolver

module purge
module load mpi/mpich/4.1.2-gcc-11.4.0-ch4
module load cmake/3.30.1-gcc-11.4.0

rm -rf build
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DFP_BUILD_BACKGROUND_COUPLING_TESTS=ON \
  -DFP_BUILD_PERIODIC_STAGGERED_TESTS=ON

cmake --build build -j4 --target \
  fp_solver \
  background_coupling_dual_u_muscl_fct_test \
  final_dual_u_limited_pairing_test \
  u_fct_budget_unit_test \
  u_fct_checkpoint_replay
```

##### C. 13.4.1加权诊断基线

13.4必须冻结13.3的最终结论。若13.3未通过，以下命令固定使用：

```text
--face-pairing-mode cell-baseline
```

如果13.3新增CLI尚未接入，而代码已经回退到当前稳定基线，则运行时删除这一行；当前程序本身即等价于`cell-baseline`。不得用已经失败的face伪逆版本执行13.4。

若13.3通过，则把所有13.4 A/B统一替换为同一组已经验收的regularized参数，不能只在某一组启用。

24 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_24" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode baseline \
  --u-fct-weighted-audit \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 24.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateA/24_weighted_baseline \
  --overwrite-output
```

30 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode baseline \
  --u-fct-weighted-audit \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 30.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateA/30_weighted_baseline \
  --overwrite-output
```

36 fs：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_36" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode baseline \
  --u-fct-weighted-audit \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 36.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateA/36_weighted_baseline \
  --overwrite-output
```

若加权统计显示核心区贡献很小，13.4应停止，不得仅凭尾部面数覆盖率继续改FCT。

##### D. 13.4.2单元测试

```bash
yhrun -N 1 -n 16 --cpu-bind=cores ./build/background_coupling_dual_u_muscl_fct_test \
  > output/u_fct_gateB_existing_fct.result 2>&1

./build/final_dual_u_limited_pairing_test \
  > output/u_fct_gateB_existing_capacity.result 2>&1

./build/u_fct_budget_unit_test \
  > output/u_fct_gateB_new_budget.result 2>&1
```

新增测试必须覆盖：

```text
单cell多面共同消费预算
两cell共享u面唯一alpha
正负电场
非均匀u网格
u端点零通量
FCT与dual-u统一容量
candidate_min_f
delta_N/delta_P/delta_K
```

##### E. 固定中点重放，按模式逐个测试

以30 fs为例，必须按以下顺序运行并分析，不能并行一次性提交全部模式。

baseline：

```bash
yhrun --cpu-bind=cores ./build/u_fct_checkpoint_replay \
  --operator-audit ./output/face_pairing_gateB/30_midpoint/final_midpoint \
  --u-fct-budget-mode baseline \
  --output-dir ./output/u_fct_gateC/30_baseline
```

第一项修改`local-joint`：

```bash
yhrun --cpu-bind=cores ./build/u_fct_checkpoint_replay \
  --operator-audit ./output/face_pairing_gateB/30_midpoint/final_midpoint \
  --u-fct-budget-mode local-joint \
  --output-dir ./output/u_fct_gateC/30_local_joint
```

只有`local-joint`通过后才运行`unified-capacity`：

```bash
yhrun --cpu-bind=cores ./build/u_fct_checkpoint_replay \
  --operator-audit ./output/face_pairing_gateB/30_midpoint/final_midpoint \
  --u-fct-budget-mode unified-capacity \
  --output-dir ./output/u_fct_gateC/30_unified_capacity
```

只有前三种模式都通过，且核心加权耗散仍明显，才运行`core-priority`：

```bash
yhrun --cpu-bind=cores ./build/u_fct_checkpoint_replay \
  --operator-audit ./output/face_pairing_gateB/30_midpoint/final_midpoint \
  --u-fct-budget-mode core-priority \
  --output-dir ./output/u_fct_gateC/30_core_priority
```

30 fs选出候选模式后，还必须对24 fs和36 fs分别重放同一模式。

##### F. 24/30/36 fs的0.2 fs生产A/B

以下以候选`local-joint`和30 fs为例。

baseline：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode baseline \
  --u-fct-weighted-audit \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 30.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateD/30_baseline_02fs \
  --overwrite-output
```

candidate：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode local-joint \
  --u-fct-weighted-audit \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --stop-time-fs 30.2 \
  --diagnostic-level 2 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateD/30_local_joint_02fs \
  --overwrite-output
```

随后用同一候选模式分别执行：

```text
24 fs: --restart "$CHECKPOINT_24" --stop-time-fs 24.2
36 fs: --restart "$CHECKPOINT_36" --stop-time-fs 36.2
```

每个时间点都必须有独立baseline，不能拿历史不同代码版本输出直接比较。

##### G. 30--42 fs物理回归

baseline：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode baseline \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --checkpoint-times 36,42 \
  --stop-time-fs 42 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateE/30_to_42_baseline \
  --overwrite-output
```

candidate：

```bash
yhrun --cpu-bind=cores stdbuf -oL -eL ./build/fp_solver \
  --restart "$CHECKPOINT_30" \
  --background-coupling-mode dual-u \
  --face-pairing-mode cell-baseline \
  --u-fct-budget-mode local-joint \
  --dt-scale 0.8333333333333333 \
  --midpoint-acceleration none \
  --midpoint-max-iters 40 \
  --checkpoint-times 36,42 \
  --stop-time-fs 42 \
  --diagnostic-level 1 \
  --beam-ledger-mode summary \
  --output-dir ./output/u_fct_gateE/30_to_42_local_joint \
  --overwrite-output
```

若最终候选不是`local-joint`，必须将命令中的模式和输出目录改成实际候选，不能保留误导性名称。

### 13.5 非阻塞维护：修复Beam哈希测试

该项只修改诊断，不改Beam推进，也不是13.3/13.4的前置条件。Beam宏观沉积、连续性和账本正常时，不能因为排列相关哈希不一致阻塞物理A/B。

推荐：

1. 若粒子具有稳定ID，按ID排序后哈希；
2. 否则使用排列无关的多重哈希，例如同时归约：

```text
sum(hash_particle)
xor(hash_particle)
sum(hash_particle^2)
particle_count
```

同时记录：

```text
total_weight
sum(w*x)
sum(w*p_x)
sum(w*gamma)
deposited_density_hash
deposited_current_hash
```

### 13.6 第六步：修改后的回归顺序

执行顺序分为四道门，不能跳级：

1. **门A：不推进时间的独立测试**
   - `G/G*`常数、低频、近Nyquist和严格交替模谱审计；
   - dual-u正负容量和u边界零通量测试；
   - 联合约束求解器的KKT残差、边界约束和单调回退测试。
2. **门B：固定checkpoint离线测试**
   - 对24/30/36 fs最终中点状态只求一次候选修正；
   - 比较真实面残差、不可解析模、修正范数、容量占用和候选动能；
   - 任一checkpoint不满足单调下降即停止，不编译进生产默认路径。
3. **门C：生产内核短A/B**
   - `diagnostic-level=0/2`一步背景和场哈希；
   - 24、30、36 fs各单步，然后各推进0.2 fs；
   - 比较完整`f/E/J/n_bkg`、中点迭代、软接受、x/u FCT宏观矩和面配对残差。
4. **门D：物理与性能回归**
   - 从30 fs推进到42 fs，比较EPOCH主峰、包络和相位；
   - 只有30--42 fs改善且3--27 fs基准未破坏，才进行0--57 fs；
   - 最后才做性能优化和120 fs生产。

13.3与13.4必须分别走完门A--D，不能在同一版本同时修改。失败版本的输出必须单独保存，不能用新运行覆盖稳定基线。

---

## 14. 判定门槛

### 14.1 x/u配对修复通过条件

- 固定状态离线求解中，`J_N_FCT-G*J_E_FCT`的可表示部分在24/30/36 fs均严格下降；
- 生产短回归中，绝大多数接受步的真实面残差下降，且任何单步都不得出现失败伪逆版本那样的倍数级放大；
- 核心区和全域残差同时下降，不能依靠空间反号抵消；
- 不可表示及近零模残差被显式记录，不能通过病态伪逆强制消除；
- 修正幅度、动能变化和场功变化均处于预设信赖域内；
- 修正后的`J_E`必须由最终u通量重新取矩得到，诊断重构误差保持舍入级；
- 背景质量误差不恶化；
- 接受态保持有限且无宏观负分布；
- 平均中点迭代数和软接受比例不恶化；
- 24/30/36 fs的`Ex/n_bkg/J_bkg`相对稳定基线只能出现与残差改善相称的小变化，不能再次出现数十个百分点漂移；
- 3--27 fs主波包不被压制；
- 30--42 fs包络和相位更接近EPOCH。

前四项属于离散验收，后六项属于物理和求解器验收。任一类失败都必须回退，不能用另一类指标改善抵消。

### 14.2 u-FCT耗散改进通过条件

- 速度核心的`|delta_N_u|`和`|delta_K_u|`显著降低；
- 质量和能量加权的核心limiter介入下降，而不是只降低未加权面数覆盖率；
- u-FCT的`R_u`仍保持舍入级；
- 最终u通量继续满足质量守恒、共享面唯一alpha和u边界零通量；
- limiter不转移成x方向更大的配对残差；
- 不增加软接受和子循环数量；
- 物理核心的电场与背景密度不被进一步平滑。

### 14.3 不能只看

以下任一项单独改善都不足以验收：

- 总能量误差变漂亮；
- limiter覆盖率下降；
- `min_f`保持非负；
- 与EPOCH某一个时刻峰值更接近；
- 单步残差下降；
- cell标量配对残差变为零；
- 只去掉严格周期交替模后的伪逆残差；
- 程序运行更快。

必须同时检查连续性、真实Yee面配对、完整背景动能、宏观波形、中点收敛和接受态稳定性。

---

## 15. 当前生产建议

在上述修复完成前：

1. 先确认集群重新编译的是已经回退face伪逆后的稳定cell-capacity代码；
2. 先从24/30/36 fs checkpoint各做一步或极短回归，确认结果恢复到稳定基线；
3. 不建议继续`0.6dt`的0--120 fs正式生产；
4. 如必须推进，使用稳定基线、原始`0.5dt`和40次上限；
5. 3--27 fs可作为较可信区间；
6. 27--36 fs只适合定性参考；
7. 36 fs以后不能用于高精度定量结论；
8. 48 fs以后同时受到波形偏差、能量累计误差和软接受影响。

失败face伪逆版本产生的24/30/36 fs结果只能作为反例，不能继续接力生产，也不能与稳定基线拼接。

不应继续：

- 提高到80次中点上限；
- 恢复无约束周期面伪逆；
- 直接整体放松或删除FCT；
- 用全局能量回填修饰能量账；
- 为了性能放宽软接受；
- 同时修改face配对和u-FCT预算；
- 在根因未修复前进行120 fs长跑。

---

## 16. 最终结论

当前可以明确区分三类结论。

**已证实的离散事实：**

> x-FCT连续性电流产生了没有同源背景动能变化对应的局域场功；u-FCT全局功配对保持到舍入级，但持续在速度核心与尾部之间重分配分布。

24/30/36 fs短审计表明：x-FCT残差分别约为 $6.03\times10^5$、$7.53\times10^5$ 和 $3.58\times10^5\ {\rm J/m^2}$，直接x/u-FCT配对残差绝对值分别约为 $2.02\times10^5$、$6.51\times10^5$ 和 $1.98\times10^5\ {\rm J/m^2}$。问题在24 fs已经存在、30 fs附近增强、36 fs又下降，因此36 fs偏差不是limiter在该时刻突然增强。更可能的过程是：

$
\text{24--30 fs前后已存在的x/u离散不配对和速度扩散}
\rightarrow
\text{误差长期累计}
\rightarrow
\text{在非线性波包阶段延迟显现}
\rightarrow
\text{36 fs后幅相和包络逐渐偏离}
\rightarrow
\text{48 fs后中点软接受进一步放大问题}.
$

这条从离散误差到波形偏差的传播链目前属于**强相关假设**，不是唯一根因已经被证明。尤其是u-FCT对波形的净影响仍缺少单因素A/B。

**已经否定的修复：**

> 无约束周期面伪逆即使去掉严格交替模，仍会放大近Nyquist弱可表示模；后续逐cell容量截断进一步破坏空间相关性。

该版本使真实面残差放大约6--7倍，在0.2 fs内产生`1e9 J/m2`量级异常背景动能变化，并显著改变核心`Ex/n_bkg/J_bkg`，因此已经回退。这个结果证明“能量账更闭合”不能凌驾于真实动力学之上。

**当前最佳下一步：**

1. 保持回退后的稳定cell-capacity生产基线；
2. 按13.3先构造带谱正则、dual-u容量约束、信赖域和完整回退的离线联合优化原型；
3. 只有真实Yee面残差在24/30/36 fs均单调下降，且动能、中点收敛和宏观状态不恶化，才进入生产短A/B；
4. 若可安全降低的残差很小，则保留并记录剩余不可解析残差，不继续增强修正；
5. 13.4的u-FCT预算优化在13.3完成判定后单独执行。

因此，当前工作的中心不再是“追求严格完全能量闭合”，而是：

> 在连续性、正性容量和已有合理波形不被破坏的前提下，只减少能够被稳定、局部且受限地消除的配对误差；其余误差诚实记账。
