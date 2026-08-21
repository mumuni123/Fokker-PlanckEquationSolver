# 相对论电子束–等离子体输运方程组及无碰撞/碰撞数值方法

## 1. 文档目的

本文针对一维空间、轴对称二维动量空间中的相对论电子束–等离子体输运问题，整理两套数值求解方案：

1. **无碰撞 Vlasov–Ampère–Gauss + 束流 PIC 方法**；
2. **包含碰撞项的 Vlasov–Fokker–Planck–Ampère–Gauss 方法**。

本文同时列出两种模型所求解的连续方程组、推荐的离散方法、时间推进方法、边界条件、守恒要求及实现路线。

当前程序的实际生产路径采用柱坐标动量网格

$$
(x,u_\parallel,u_\perp),
$$

背景电子和电场采用周期空间边界，束流电子从左边界注入并在全局边界处开放流出。当前默认关闭碰撞项，因此实际求解的是无碰撞 Vlasov–Ampère–Gauss 与 PIC 束流耦合系统。

---

## 2. 公共变量与相对论关系

### 2.1 归一化动量

定义电子归一化动量

$$
u_\parallel=\frac{p_x}{m_ec},
\qquad
u_\perp=\frac{p_\perp}{m_ec}.
$$

总归一化动量大小为

$$
u=\sqrt{u_\parallel^2+u_\perp^2}.
$$

相对论因子为

$$
\gamma(u_\parallel,u_\perp)
=
\sqrt{1+u_\parallel^2+u_\perp^2}.
$$

纵向速度为

$$
v_x
=
c\frac{u_\parallel}{\gamma}.
$$

### 2.2 柱坐标动量体积元

在轴对称柱坐标动量空间中，

$$
\mathrm d^3u
=
2\pi u_\perp\,
\mathrm du_\perp\mathrm du_\parallel.
$$

背景电子分布函数记为

$$
f_e=f_e(x,u_\parallel,u_\perp,t).
$$

电子数密度为

$$
n_e(x,t)
=
2\pi
\int_0^\infty
\int_{-\infty}^{\infty}
f_e(x,u_\parallel,u_\perp,t)
u_\perp\,
\mathrm du_\parallel\mathrm du_\perp.
$$

背景电子电流密度为

$$
J_{\mathrm{bkg}}(x,t)
=
q_e\,2\pi
\int_0^\infty
\int_{-\infty}^{\infty}
v_x f_e
u_\perp\,
\mathrm du_\parallel\mathrm du_\perp,
$$

其中

$$
q_e=-e.
$$

---

# 第一部分：无碰撞方法

## 3. 无碰撞连续方程组

### 3.1 背景电子 Vlasov 方程

背景电子满足相对论 Vlasov 方程

$$
\frac{\partial f_e}{\partial t}
+
v_x\frac{\partial f_e}{\partial x}
-
\frac{eE_x(x,t)}{m_ec}
\frac{\partial f_e}{\partial u_\parallel}
=0.
$$

定义归一化动量加速度

$$
a_e(x,t)
=
\frac{q_eE_x(x,t)}{m_ec}
=
-\frac{eE_x(x,t)}{m_ec},
$$

则守恒形式为

$$
\boxed{
\frac{\partial f_e}{\partial t}
+
\frac{\partial}{\partial x}
\left(v_xf_e\right)
+
\frac{\partial}{\partial u_\parallel}
\left(a_ef_e\right)
=0.
}
$$

因为纯纵向静电场不改变 $u_\perp$，无碰撞方程中不存在 $u_\perp$ 方向输运项。

---

### 3.2 束流 PIC 方程

每个束流宏粒子满足

$$
\frac{\mathrm dx_p}{\mathrm dt}
=
\frac{p_{x,p}}{\gamma_pm_e}
=
c\frac{u_p}{\sqrt{1+u_p^2}},
$$

$$
\frac{\mathrm dp_{x,p}}{\mathrm dt}
=
-eE_x(x_p,t),
$$

其中

$$
u_p=\frac{p_{x,p}}{m_ec}.
$$

束流数密度由粒子形函数沉积：

$$
n_b(x,t)
=
\sum_pW_pS(x-x_p).
$$

若使用电荷守恒轨迹沉积，束流电荷与电流满足离散连续性关系

$$
\frac{\partial\rho_b}{\partial t}
+
\frac{\partial J_b}{\partial x}
=
S_{\rho,b},
$$

其中 $S_{\rho,b}$ 表示开放边界注入和流出引起的电荷源。

---

### 3.3 总电荷密度

固定离子背景满足

$$
Zn_i(x)=n_0.
$$

总电荷密度为

$$
\boxed{
\rho(x,t)
=
e\left[
Zn_i(x)-n_e(x,t)-n_b(x,t)
\right].
}
$$

---

### 3.4 Gauss 方程

对于周期电场，非零平均电荷不能由周期电场拓扑支持，因此当前模型使用

$$
\boxed{
\frac{\partial E_x}{\partial x}
=
\frac{
\rho(x,t)-\langle\rho\rangle_x
}{\varepsilon_0},
}
$$

其中

$$
\langle\rho\rangle_x
=
\frac{1}{L_x}
\int_0^{L_x}\rho(x,t)\,\mathrm dx.
$$

初始时刻通常由该方程构造与初始电荷分布相容的电场。

---

### 3.5 Ampère 方程

纵向电场由总电流推进：

$$
\boxed{
\frac{\partial E_x}{\partial t}
=
-\frac{
J_{\mathrm{bkg}}+J_b
}{\varepsilon_0}.
}
$$

若不主动去除平均电流，则平均电场满足

$$
\frac{\mathrm d\langle E_x\rangle_x}{\mathrm dt}
=
-
\frac{\langle J_{\mathrm{bkg}}+J_b\rangle_x}
{\varepsilon_0}.
$$

---

## 4. 无碰撞推荐数值方法

## 4.1 背景电子：守恒有限体积法

推荐直接存储相空间单元积分质量

$$
M_{i,j,k}
=
\int_{\Delta x_i}
\int_{\Delta u_{\parallel,j}}
\int_{\Delta u_{\perp,k}}
f_e\,
2\pi u_\perp\,
\mathrm du_\perp
\mathrm du_\parallel
\mathrm dx.
$$

相空间单元体积为

$$
\Delta V_{u,jk}
=
\Delta u_{\parallel,j}
\pi
\left(
u_{\perp,k+1/2}^2-u_{\perp,k-1/2}^2
\right).
$$

因此

$$
M_{i,j,k}
\approx
f_{e,i,j,k}
\Delta x_i
\Delta V_{u,jk}.
$$

对一个时间子步 $h$，二维非分裂有限体积更新为

$$
M_{i,j,k}^{m+1}
=
M_{i,j,k}^{m}
-h
\left(
\Phi^x_{i+1/2,j,k}
-
\Phi^x_{i-1/2,j,k}
\right)
-h
\left(
\Phi^u_{i,j+1/2,k}
-
\Phi^u_{i,j-1/2,k}
\right).
$$

其中：

- $\Phi^x$ 为实空间质量通量；
- $\Phi^u$ 为 $u_\parallel$ 方向质量通量。

---

## 4.2 低阶正性通量

低阶基线建议使用 donor-cell 迎风格式：

$$
\Phi^{x,L}
=
v_x
\frac{M_{\mathrm{upwind}}}{\Delta x},
$$

$$
\Phi^{u,L}
=
a_e
\frac{M_{\mathrm{upwind}}}
{\Delta u_{\parallel,\mathrm{donor}}}.
$$

低阶状态应构造成非负质量项之和，从而作为后续高阶修正的正性基线。

---

## 4.3 高阶重构

可使用以下方法：

1. MUSCL；
2. PPM；
3. WENO-Z；
4. 高阶 DG。

实际工程中推荐：

$$
\boxed{
\text{有限体积}
+
\text{PPM 或 WENO-Z}
+
\text{正性限制器}.
}
$$

MUSCL 实现简单，但数值耗散较大。PPM 是精度与成本之间较好的折中。WENO-Z 更适合相空间中出现强梯度或细结构的情况。

---

## 4.4 FCT 正性限制

将高阶通量写成

$$
\Phi^H=\Phi^L+A,
$$

其中 $A$ 为反扩散通量。最终通量取

$$
\Phi^{\mathrm{final}}
=
\Phi^L+\alpha A,
\qquad
0\leq\alpha\leq1.
$$

限制系数 $\alpha$ 根据供体单元可用的低阶非负质量确定，以同时保证：

- 单元质量非负；
- 尽可能保留高阶精度；
- 共享面通量严格守恒。

---

## 4.5 CFL 与子循环

组合 CFL 数为

$$
C_{\mathrm{comb}}
=
\max_{i,j,k}
\left[
\frac{\Delta t|v_{x,jk}|}{\Delta x_i}
+
\frac{\Delta t|a_{e,i}|}
{\Delta u_{\parallel,j}}
\right].
$$

内部子步数可取

$$
N_{\mathrm{sub}}
=
\max
\left[
1,
\left\lceil
\frac{C_{\mathrm{comb}}}{C_{\mathrm{target}}}
\right\rceil
\right],
$$

其中通常取

$$
C_{\mathrm{target}}\approx0.8.
$$

子步长度为

$$
h=\frac{\Delta t}{N_{\mathrm{sub}}}.
$$

---

## 4.6 Vlasov–Ampère 中点耦合

推荐使用时间中心耦合：

$$
E^{n+1/2,(r)}
=
\frac{
E^n+E^{n+1,(r)}
}{2},
$$

$$
f^{n+1/2,(r)}
=
\frac{
f^n+f^{n+1,(r)}
}{2}.
$$

背景电流应直接由同一次有限体积空间通量计算：

$$
J_{\mathrm{bkg},i+1/2}^{n+1/2}
=
q_e
\sum_{j,k}
\overline{\Phi^x_{i+1/2,j,k}},
$$

其中横线表示整个时间步内的时间平均。

电场更新为

$$
E_{i+1/2}^{n+1}
=
E_{i+1/2}^{n}
-
\frac{\Delta t}{\varepsilon_0}
\left(
J_{\mathrm{bkg},i+1/2}^{n+1/2}
+
J_{b,i+1/2}^{n+1/2}
\right).
$$

这样得到的电流与离散连续性方程严格一致。

---

## 4.7 束流推进的两种选择

### 方法 A：保留 PIC 束流

推荐采用：

$$
\boxed{
\text{守恒 Vlasov 背景}
+
\text{电荷守恒 PIC 束流}
+
\text{全耦合中点迭代}.
}
$$

在每一次场的非线性迭代中，都应使用当前中点电场重新推进束流并重新沉积守恒电流，避免束流电流对场迭代显式滞后。

### 方法 B：束流也改为 Eulerian Vlasov

若主要研究微弱电场结构和增长率，推荐将束流写成

$$
f_b=f_b(x,u_b,t),
$$

并求解

$$
\frac{\partial f_b}{\partial t}
+
\frac{\partial}{\partial x}
\left[
c\frac{u_b}{\sqrt{1+u_b^2}}f_b
\right]
+
\frac{\partial}{\partial u_b}
\left[
-\frac{eE_x}{m_ec}f_b
\right]
=0.
$$

束流密度和电流为

$$
n_b(x,t)
=
\int f_b(x,u_b,t)\,\mathrm du_b,
$$

$$
J_b(x,t)
=
-e c
\int
\frac{u_b}{\sqrt{1+u_b^2}}
f_b(x,u_b,t)\,
\mathrm du_b.
$$

这种方法不存在有限宏粒子数导致的统计噪声，更适合观察精细电场演化。

---

## 4.8 非线性求解器

推荐顺序为

$$
\boxed{
\text{Picard}
\rightarrow
\text{Anderson-Picard}
\rightarrow
\text{JFNK}.
}
$$

第一版可使用 Picard 迭代与松弛；随后建议增加 Anderson acceleration。对于强耦合和大时间步，可采用 Jacobian-Free Newton–Krylov 方法。

若非线性迭代不收敛，不应直接接受最后一次迭代状态，而应：

$$
\Delta t
\leftarrow
\frac{\Delta t}{2},
$$

然后重新计算该时间步。

---

## 5. 无碰撞方法的主要优缺点

### 优点

- 背景电子粒子数可严格守恒；
- 可以保持分布函数正性；
- 电流可与离散连续性方程一致；
- 适合研究无碰撞相空间输运、双流不稳定性和电场增长；
- 相比完整碰撞模型，计算成本和实现难度较低。

### 局限

- 无法描述电子–电子和电子–离子碰撞；
- 无法描述束流慢化、俯仰角散射和碰撞热化；
- 高密度、低温等离子体中可能遗漏主要输运机制；
- 周期电场与开放束流边界之间需要额外的全局电荷、电流闭合模型。

---

# 第二部分：包含碰撞项的完整 Fokker–Planck 方法

## 6. 速度空间 Fokker–Planck 方程

用户给出的速度空间形式为

$$
\frac{\partial f_e}{\partial t}
+
v_x
\frac{\partial f_e}{\partial x}
-
\frac{eE(x,t)}{m_e}
\frac{\partial f_e}{\partial v_x}
=
-\frac{\partial}{\partial\boldsymbol v}
\left(
f_e
\langle\Delta\boldsymbol v\rangle
\right)
+
\frac12
\frac{\partial^2}
{\partial\boldsymbol v\partial\boldsymbol v}
:
\left(
f_e
\langle
\Delta\boldsymbol v\Delta\boldsymbol v
\rangle
\right).
$$

漂移矩为

$$
\left\langle\Delta\boldsymbol v\right\rangle
=
\frac{1}{\Delta t}
\int
\omega
\left(
\boldsymbol v,\Delta\boldsymbol v
\right)
\Delta\boldsymbol v\,
\mathrm d(\Delta\boldsymbol v),
$$

扩散矩为

$$
\left\langle
\Delta\boldsymbol v\Delta\boldsymbol v
\right\rangle
=
\frac{1}{\Delta t}
\int
\omega
\left(
\boldsymbol v,\Delta\boldsymbol v
\right)
\Delta\boldsymbol v
\Delta\boldsymbol v\,
\mathrm d(\Delta\boldsymbol v).
$$

该表达式是 Kramers–Moyal 展开保留到二阶后的漂移–扩散形式。

---

## 7. 相对论动量空间形式

由于当前模型采用相对论动量变量，更一致的写法是

$$
\boxed{
\frac{\partial f_e}{\partial t}
+
c\frac{u_\parallel}{\gamma}
\frac{\partial f_e}{\partial x}
-
\frac{eE_x}{m_ec}
\frac{\partial f_e}{\partial u_\parallel}
=
C_{\mathrm{rFP}}[f_e].
}
$$

相对论 Fokker–Planck 算子写为

$$
C_{\mathrm{rFP}}[f_e]
=
-
\frac{\partial}{\partial u_i}
\left(
A_i f_e
\right)
+
\frac12
\frac{\partial^2}{\partial u_i\partial u_j}
\left(
D_{ij}f_e
\right),
$$

其中

$$
A_i
=
\frac{
\langle\Delta u_i\rangle
}{\Delta t},
$$

$$
D_{ij}
=
\frac{
\langle\Delta u_i\Delta u_j\rangle
}{\Delta t}.
$$

对于完整相对论库仑碰撞，应使用相对论 Landau 或 Beliaev–Budker 型碰撞核，而不能简单地把非相对论速度空间系数直接替换成动量空间系数。

---

## 8. 完整碰撞方程组

包含碰撞后，可写成

$$
\boxed{
\frac{\partial f_e}{\partial t}
+
\frac{\partial}{\partial x}
\left(
v_xf_e
\right)
+
\frac{\partial}{\partial u_\parallel}
\left(
a_ef_e
\right)
=
C_{ee}[f_e,f_e]
+
C_{ei}[f_e,f_i].
}
$$

其中

$$
a_e
=
-\frac{eE_x}{m_ec}.
$$

电子密度为

$$
n_e(x,t)
=
2\pi
\int_0^\infty
\int_{-\infty}^{\infty}
f_e
u_\perp\,
\mathrm du_\parallel
\mathrm du_\perp.
$$

背景电流为

$$
J_{\mathrm{bkg}}(x,t)
=
-e\,2\pi
\int_0^\infty
\int_{-\infty}^{\infty}
v_xf_e
u_\perp\,
\mathrm du_\parallel
\mathrm du_\perp.
$$

总电荷密度仍为

$$
\rho
=
e
\left[
Zn_i-n_e-n_b
\right].
$$

Gauss 方程为

$$
\frac{\partial E_x}{\partial x}
=
\frac{
\rho-\langle\rho\rangle_x
}{\varepsilon_0},
$$

Ampère 方程为

$$
\frac{\partial E_x}{\partial t}
=
-
\frac{
J_{\mathrm{bkg}}+J_b
}{\varepsilon_0}.
$$

若束流仍采用 PIC，则粒子方程保持

$$
\frac{\mathrm dx_p}{\mathrm dt}
=
c
\frac{u_p}{\sqrt{1+u_p^2}},
$$

$$
\frac{\mathrm du_p}{\mathrm dt}
=
-\frac{eE_x(x_p,t)}{m_ec},
$$

但还需要额外加入束流与背景之间的碰撞模型。

---

## 9. 线性与非线性碰撞算子

## 9.1 固定背景测试粒子算子

若漂移和扩散系数由固定 Maxwell 或 Maxwell–Jüttner 背景决定，则

$$
A_i=A_i(\boldsymbol u),
\qquad
D_{ij}=D_{ij}(\boldsymbol u),
$$

碰撞算子对 $f_e$ 为线性算子。

这种方法适合描述：

- 束流减速；
- 能量扩散；
- 俯仰角散射；
- 向固定热背景弛豫。

但是它通常不能同时保证被模拟电子系统自身的总动量和总能量守恒。

---

## 9.2 完整非线性电子–电子算子

对于自洽电子–电子碰撞，

$$
A_i=A_i[f_e],
\qquad
D_{ij}=D_{ij}[f_e].
$$

因此

$$
C_{ee}[f_e,f_e]
$$

是非线性算子。

它必须满足

$$
\int
C_{ee}
\,\mathrm d^3u
=0,
$$

$$
\int
\boldsymbol p
C_{ee}
\,\mathrm d^3u
=0,
$$

$$
\int
(\gamma-1)m_ec^2
C_{ee}
\,\mathrm d^3u
=0.
$$

即电子–电子碰撞应守恒：

1. 电子数；
2. 总动量；
3. 总能量。

---

## 9.3 电子–离子碰撞

完整碰撞项写为

$$
C_e
=
C_{ee}
+
C_{ei}.
$$

在无限重离子近似下，电子–离子碰撞主要表现为俯仰角散射。若使用球坐标动量变量，其典型形式为

$$
C_{ei}[f_e]
=
\frac{\nu_D(u)}{2}
\frac{\partial}{\partial\mu}
\left[
(1-\mu^2)
\frac{\partial f_e}{\partial\mu}
\right].
$$

该算子守恒电子数，主要改变电子运动方向，并把平行动量传递给离子背景。

---

## 10. 推荐动量坐标：$(u,\mu)$

对完整碰撞算子，推荐定义

$$
u
=
\sqrt{
u_\parallel^2+u_\perp^2
},
$$

$$
\mu
=
\frac{u_\parallel}{u},
\qquad
-1\leq\mu\leq1.
$$

于是

$$
u_\parallel=u\mu,
$$

$$
u_\perp
=
u\sqrt{1-\mu^2}.
$$

动量体积元变成

$$
\mathrm d^3u
=
2\pi u^2
\mathrm du
\mathrm d\mu.
$$

相对论因子为

$$
\gamma=\sqrt{1+u^2},
$$

纵向速度为

$$
v_x
=
c
\frac{u\mu}{\gamma}.
$$

---

## 11. $(x,u,\mu)$ 坐标中的 Vlasov–Fokker–Planck 方程

令

$$
a(x,t)
=
-\frac{eE_x(x,t)}{m_ec}.
$$

纵向电场在球坐标动量空间中产生

$$
\dot u=a\mu,
$$

$$
\dot\mu
=
a
\frac{1-\mu^2}{u}.
$$

因此非守恒形式为

$$
\boxed{
\frac{\partial f_e}{\partial t}
+
c\frac{u\mu}{\gamma}
\frac{\partial f_e}{\partial x}
+
a\mu
\frac{\partial f_e}{\partial u}
+
a
\frac{1-\mu^2}{u}
\frac{\partial f_e}{\partial\mu}
=
C_e[f_e].
}
$$

带动量 Jacobian

$$
J_u=2\pi u^2
$$

的守恒形式为

$$
\boxed{
\begin{aligned}
\frac{\partial}{\partial t}
\left(
J_uf_e
\right)
&+
\frac{\partial}{\partial x}
\left(
J_u
c\frac{u\mu}{\gamma}
f_e
\right)
\\
&+
\frac{\partial}{\partial u}
\left(
J_ua\mu f_e
\right)
\\
&+
\frac{\partial}{\partial\mu}
\left[
J_u
a\frac{1-\mu^2}{u}
f_e
\right]
=
J_uC_e[f_e].
\end{aligned}
}
$$

相空间单元质量可定义为

$$
M_{i,j,k}
=
\int_{\Delta x_i}
\int_{\Delta u_j}
\int_{\Delta\mu_k}
f_e
2\pi u^2
\,\mathrm d\mu
\mathrm du
\mathrm dx.
$$

---

## 12. 碰撞项的守恒通量形式

碰撞算子应写成动量空间通量散度：

$$
\boxed{
C_e[f_e]
=
-
\frac{1}{J_u}
\left[
\frac{\partial}{\partial u}
\left(
J_u\Gamma_u^c
\right)
+
\frac{\partial}{\partial\mu}
\left(
J_u\Gamma_\mu^c
\right)
\right].
}
$$

碰撞通量一般写为

$$
\Gamma_\alpha^c
=
A_\alpha f_e
-
D_{\alpha u}
\frac{\partial f_e}{\partial u}
-
D_{\alpha\mu}
\frac{\partial f_e}{\partial\mu}.
$$

二维动量有限体积离散为

$$
M_{j,k}^{n+1}
=
M_{j,k}^{n}
-
\Delta t
\left(
F^u_{j+1/2,k}
-
F^u_{j-1/2,k}
\right)
-
\Delta t
\left(
F^\mu_{j,k+1/2}
-
F^\mu_{j,k-1/2}
\right).
$$

共享面只使用一个碰撞通量，可自然保证粒子数守恒。

---

## 13. 碰撞项的时间推进

碰撞扩散项若显式推进，会受到

$$
\Delta t
\lesssim
\frac{\Delta u^2}{D_{uu}},
$$

以及

$$
\Delta t
\lesssim
\frac{\Delta\mu^2}{D_{\mu\mu}}
$$

的稳定性限制。

因此碰撞项应采用隐式方法。

推荐顺序为

$$
\boxed{
\text{SDIRK2 或 TR-BDF2}
>
\text{BDF2}
>
\text{后向 Euler}
>
\text{Crank–Nicolson}.
}
$$

### 13.1 后向 Euler

$$
\frac{
f^{n+1}-f^n
}{\Delta t}
=
C[f^{n+1}].
$$

优点是稳定、实现简单、适合第一版程序；缺点是一阶时间精度。

### 13.2 BDF2

$$
\frac{
3f^{n+1}
-
4f^n
+
f^{n-1}
}{2\Delta t}
=
C[f^{n+1}].
$$

BDF2 为二阶方法，适合长期生产计算。

### 13.3 SDIRK2

SDIRK2 兼具二阶精度和良好的刚性衰减能力，适合强碰撞扩散问题。

---

## 14. Vlasov 与碰撞的 Strang 分裂

将方程写成

$$
\frac{\partial f_e}{\partial t}
=
\mathcal V(f_e,E_x)
+
\mathcal C(f_e),
$$

其中：

- $\mathcal V$ 表示实空间输运、电场加速和 Ampère 耦合；
- $\mathcal C$ 表示碰撞算子。

推荐采用二阶 Strang 分裂：

$$
\boxed{
C_{\Delta t/2}
\rightarrow
(Vlasov\text{--}Ampère)_{\Delta t}
\rightarrow
C_{\Delta t/2}.
}
$$

### 第一步：碰撞半步

$$
\frac{
f^{(1)}-f^n
}{\Delta t/2}
=
C[f^{(1)}].
$$

### 第二步：Vlasov–Ampère 整步

$$
f^{(1)}
\longrightarrow
f^{(2)},
$$

同时使用通量一致电流推进

$$
E^{n+1}
=
E^n
-
\frac{\Delta t}{\varepsilon_0}
J^{n+1/2}.
$$

### 第三步：碰撞半步

$$
\frac{
f^{n+1}-f^{(2)}
}{\Delta t/2}
=
C[f^{n+1}].
$$

碰撞步骤在每个空间单元中独立求解二维动量问题，因此非常适合：

- MPI 沿 $x$ 分解；
- OpenMP 并行；
- GPU 批量求解；
- 每个空间单元独立构造隐式矩阵或预条件器。

---

## 15. 碰撞隐式非线性求解器

对非线性电子–电子碰撞，隐式方程为

$$
R(f^{n+1})
=
f^{n+1}
-
f^\ast
-
\Delta t
C_{ee}[f^{n+1}]
=0.
$$

推荐采用

$$
\boxed{
\text{Picard}
+
\text{Anderson acceleration}
}
$$

作为第一版非线性求解器。

一次 Picard 迭代可写成：

1. 用 $f^{(r)}$ 计算碰撞势或碰撞系数；
2. 构造 $A^{(r)}$ 和 $D^{(r)}$；
3. 固定系数求解线性隐式漂移–扩散方程；
4. 得到 $f^{(r+1)}$；
5. 进行守恒修正；
6. 检查非线性残差和守恒矩残差。

更高级的方法为

$$
\boxed{
\text{Jacobian-Free Newton–Krylov}
+
\text{物理预条件器}.
}
$$

---

## 16. 碰撞离散必须满足的性质

### 16.1 粒子数守恒

$$
\sum_{j,k}
C_{j,k}
\Delta V_{u,j,k}
=0.
$$

### 16.2 动量守恒

对于电子–电子碰撞，

$$
\sum_{j,k}
p_{\parallel,j,k}
C_{j,k}
\Delta V_{u,j,k}
=0.
$$

### 16.3 能量守恒

$$
\sum_{j,k}
(\gamma_{j,k}-1)
m_ec^2
C_{j,k}
\Delta V_{u,j,k}
=0.
$$

### 16.4 正性

$$
f_{i,j,k}^{n+1}\geq0.
$$

### 16.5 H 定理

对于封闭碰撞系统，

$$
\frac{\mathrm d}{\mathrm dt}
\int
f_e\ln f_e
\,\mathrm d^3u
\leq0.
$$

最终平衡态应趋向 Maxwell–Jüttner 分布，而不是数值离散产生的伪平衡态。

---

## 17. 正性与平衡保持离散

对于一维漂移–扩散方程

$$
\frac{\partial f}{\partial t}
=
\frac{\partial}{\partial u}
\left[
D(u)
\frac{\partial f}{\partial u}
-
A(u)f
\right],
$$

可使用 Chang–Cooper 或 Scharfetter–Gummel 型指数拟合通量。

这些方法有利于：

- 保持 $f\geq0$；
- 准确保持给定平衡态；
- 处理摩擦主导区域；
- 处理扩散主导区域。

对于二维碰撞算子，还需正确处理交叉扩散

$$
D_{u\mu}\neq0.
$$

因此不能简单地分别对 $u$ 和 $\mu$ 使用两个独立的一维 Chang–Cooper 更新。推荐：

1. $u$ 方向使用指数拟合漂移–扩散通量；
2. $\mu$ 方向使用对称有限体积角扩散通量；
3. 交叉扩散采用共享角点梯度或混合通量；
4. 施加离散动量和能量守恒修正；
5. 每一步检查最小分布函数值和碰撞矩误差。

---

## 18. 背景电子和束流电子的碰撞耦合

若仍将背景电子和束流电子分开，则完整碰撞项应包括

$$
C_{\mathrm{bkg,bkg}},
$$

$$
C_{\mathrm{bkg,b}},
$$

$$
C_{b,\mathrm{bkg}},
$$

$$
C_{b,b}.
$$

交叉碰撞必须满足

$$
\Delta P_{\mathrm{bkg,b}}
+
\Delta P_{b,\mathrm{bkg}}
=0,
$$

$$
\Delta W_{\mathrm{bkg,b}}
+
\Delta W_{b,\mathrm{bkg}}
=0.
$$

如果背景采用 Eulerian Vlasov，而束流采用 PIC，则需要：

- 根据 Eulerian 背景计算粒子碰撞；
- 把粒子损失的动量和能量反馈给网格；
- 保证网格–粒子之间的动量和能量闭合；
- 控制 Monte Carlo 碰撞噪声。

该方案实现难度较高。

更干净的方法是将背景和束流写成一个总电子分布：

$$
f_e
=
f_{\mathrm{thermal}}
+
f_{\mathrm{beam}},
$$

束流由左边界流入条件产生。这样电子–电子碰撞直接作用于总分布：

$$
C_{ee}[f_e,f_e].
$$

但这要求重新设计当前背景电子周期边界与束流开放边界不一致的问题。

---

## 19. 含碰撞方法的推荐开发路线

### 阶段 1：守恒模型碰撞算子

先实现简化的 Dougherty 或 Lenard–Bernstein 型模型：

$$
C_M[f]
=
\nu
\nabla_u
\cdot
\left[
(\boldsymbol u-\boldsymbol U)f
+
\Theta\nabla_u f
\right],
$$

其中 $\boldsymbol U$ 和 $\Theta$ 由当前分布矩计算。

主要验证：

- 粒子数守恒；
- 总动量守恒；
- 总能量守恒；
- 正性；
- 熵变化；
- 向平衡分布弛豫；
- Strang 分裂二阶精度。

### 阶段 2：电子–离子 Lorentz 散射

实现

$$
C_{ei}
=
\frac{\nu_D(u)}{2}
\frac{\partial}{\partial\mu}
\left[
(1-\mu^2)
\frac{\partial f}{\partial\mu}
\right].
$$

### 阶段 3：线性化电子–电子算子

以局域 Maxwell–Jüttner 分布为背景：

$$
f_e=f_M+\delta f,
$$

并实现

$$
C_{ee}^{\mathrm{lin}}[\delta f]
=
C_{\mathrm{test}}[\delta f]
+
C_{\mathrm{field}}[\delta f].
$$

### 阶段 4：完整非线性相对论 Landau 算子

实现由 $f_e$ 自洽决定的摩擦和扩散系数，使用全隐式非线性求解，并在离散层面强制粒子数、动量和能量守恒。

### 阶段 5：必要时加入大角碰撞

Fokker–Planck 算子主要描述大量小角库仑碰撞。对于高能电子的大角散射或 knock-on 过程，必要时应加入额外 Boltzmann 碰撞项或源项。

---

# 第三部分：两种方法的对比与最终建议

## 20. 方法对比

| 项目 | 无碰撞方法 | 完整 Fokker–Planck 方法 |
|---|---|---|
| 主方程 | Vlasov–Ampère–Gauss | Vlasov–Fokker–Planck–Ampère–Gauss |
| 背景电子 | Eulerian Vlasov | Eulerian Vlasov–Fokker–Planck |
| 束流 | PIC 或 Eulerian Vlasov | 最好与背景统一为总电子分布 |
| 动量坐标 | $(u_\parallel,u_\perp)$ 较方便 | $(u,\mu)$ 更方便 |
| 空间离散 | 守恒有限体积 | 守恒有限体积 |
| Vlasov 时间推进 | 中点迭代 | 中点迭代 |
| 碰撞时间推进 | 无 | 隐式 SDIRK2/BDF2 |
| 总体耦合 | Vlasov–Ampère 全耦合 | Strang 分裂或全耦合 IMEX |
| 主要限制 | PIC 噪声、边界闭合 | 刚性、非线性、守恒和正性 |
| 计算成本 | 较低 | 较高 |
| 适用问题 | 无碰撞不稳定性、场演化 | 慢化、散射、热化、电阻和碰撞阻尼 |

---

## 21. 无碰撞方案的最终推荐

$$
\boxed{
\begin{aligned}
&\text{背景电子：}
(x,u_\parallel,u_\perp)
\text{ 守恒有限体积};\\
&\text{高阶通量：PPM/WENO-Z + FCT 正性限制};\\
&\text{束流：守恒 PIC，或改为 }(x,u_b)\text{ Vlasov};\\
&\text{电场：通量一致 Ampère 更新};\\
&\text{时间：全耦合隐式中点};\\
&\text{非线性求解：Anderson-Picard 或 JFNK};\\
&\text{失败时间步：拒绝并减小 }\Delta t.
\end{aligned}
}
$$

如果主要研究电场精细结构，推荐束流也改成 Eulerian Vlasov，以消除 PIC 统计噪声。

---

## 22. 含碰撞方案的最终推荐

$$
\boxed{
\begin{aligned}
&\text{分布函数：}
f_e(x,u,\mu,t);\\
&\text{Vlasov 离散：守恒有限体积 + 正性限制};\\
&\text{碰撞算子：相对论 Landau/Fokker--Planck};\\
&\text{电子--离子：Lorentz 角散射或有限质量算子};\\
&\text{碰撞时间推进：SDIRK2、TR-BDF2 或 BDF2};\\
&\text{总体耦合：二阶 Strang 分裂};\\
&\text{非线性求解：Anderson-Picard 或 JFNK};\\
&\text{离散约束：粒子数、动量、能量、正性和 H 定理};\\
&\text{失败时间步：拒绝并减小 }\Delta t.
\end{aligned}
}
$$

---

## 23. 对当前程序的具体修改建议

1. **无碰撞版本可以继续保留当前柱坐标单元积分质量结构。**
2. **束流 predictor 不应在背景电子–电场 Picard 迭代中冻结。**
3. **中点迭代失败后应拒绝当前时间步，而不是软接受。**
4. **应增加离散总能量和边界能流残差诊断。**
5. **当前旧碰撞模块使用旧球坐标数据布局，不能直接接入现有柱坐标质量数组。**
6. **若正式加入完整碰撞，建议统一重构为 $(x,u,\mu)$ 数据结构。**
7. **第一版先实现 Lorentz 角散射和守恒模型碰撞算子。**
8. **验证完成后再实现线性化和完整非线性相对论 Landau 算子。**
9. **必须重新审视周期背景电子/电场与开放束流之间的电荷和电流闭合。**
10. **若束流与背景统一为总电子分布，则空间边界条件也必须统一重构。**

---

## 24. 总结

不含碰撞时，最合适的方法是

$$
\boxed{
\text{守恒正性 Vlasov 有限体积}
+
\text{守恒 PIC 或 Eulerian 束流}
+
\text{通量一致 Ampère}
+
\text{全耦合中点迭代}.
}
$$

包含完整碰撞时，最合适的方法是

$$
\boxed{
\text{守恒 Vlasov 有限体积}
+
\text{隐式相对论 Fokker–Planck}
+
\text{Strang 分裂}
+
\text{严格粒子数、动量、能量和正性控制}.
}
$$

加入碰撞后，数值问题的核心从单纯的双曲输运转变为：

$$
\boxed{
\text{双曲输运}
+
\text{刚性非线性动量扩散}
+
\text{电场自洽耦合}.
}
$$

因此，完整程序的关键不是简单地在右端增加一个碰撞项，而是重新设计动量坐标、碰撞通量、隐式求解器、守恒修正和整体时间耦合结构。
