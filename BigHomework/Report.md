# Report

> PB20020480  王润泽
>

## 1 Question

​	研究2维XY模型+一维Ising相互作用体系的相变

## 2 Background

### 2.1 Ising模型

​	Ising模型是统计物理中一个重要的模型，用于描述在不同温度下自发磁化的现象。它由二维晶格上的一组自旋组成，每个自旋取值为+1或-1。自旋与周围自旋的相互作用通过两个参数描述：J代表正的相互作用，也就是说，两个自旋越相似（两个+1或两个-1），其相互作用越大；其哈密顿量为:
$$
H=−J\sum_{\left<i,j\right>}{\sigma_i\sigma_j}
$$
​	对于一维Ising模型没有相变解，而二维时存在相变。

### 2.2 二维XY模型

​	2维XY模型是一个简单的经典统计物理系统，其中二维正方格子上的自旋在连续空间范围内取值。在此模型中，自旋之间存在相邻交互作用，其哈密顿量可以写成以下形式：
$$
H=−J\sum_{\left<i,j\right>}{\cos(\theta_i−\theta_j)}
$$
​	其中，$J$ 表示耦合常数，$\theta_i$ 是位于格点 $(i,j)$ 处的自旋角度。

​	这个模型存在一个连续的相变，当温度低于某个临界温度时，系统的自旋将具有长程有序性。尽管系统的平均磁化强度为零，但系统可以存在一种亚稳态，自旋的排列形成涡旋。在相变温度以上，涡旋是自由的。在相变温度之下，自旋涡旋是成对出现 的，并且对于$T < T_{KT}$ 的所有温度系统都和 $T=T_{KT}$ 时一样，因此临界点实际上是临界线。XY模型模拟结果如下：

<center class="half">
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\XY_model\T_0.1.png" style="zoom:30%;" />
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\XY_model\xy_Cv.png" style="zoom:40%;" />
<center><p>图1 左：稳定下XY模型旋度；右：热容随温度变化曲线

<div STYLE="page-break-after: always;"></div>

## 3 Experiment

### 3.0 实验方法---Monte Carlo 模拟

​	Monte Carlo模拟的主要过程是Metropolis 演化算法，一共有N个粒子，每个粒子只与相邻的粒子产生相互作用。则以下为算法的基本形式

1. 随机选出一个自旋处于 $\mu$ 状态，计算所有和其自旋相关的能量贡献，为 $H_{\mu}$ 。(也就是所有和它相邻自旋交互作用的能量贡献)
2. 随机的改变其自旋状态，变为状态 $\nu$ ，再次计算一次所有和它相关的能量贡献，为 $H_\nu$
3. 如果其能量贡献下降，状态改变为 $\nu$。
4. 如果贡献能量上升，则令自旋有 $\exp{-\beta(H_\nu-H_\mu)}$的概率接受改变。如果接受，状态变为 $\nu$，否则保持原状为 $\mu$
5. 重复步骤一，直到能量趋于稳定

假设迭代经过n步骤，去掉前面热化的m步，则剩下$N_E=n-m$的状态作为统计量统计，能量E，磁化强度M，热容 $Cv$，寻找相变点。
$$
\left<E\right>=\frac1{N_E}\sum_{k=1}^{N_E}E_k\\
\left<M\right>=\frac1{N_E}\sum_{k=1}^{N_E}(\mu_B\sum_{i=1}^N\sigma_i)\\
C= \frac{\Delta E^2}{k_BT^2}=\frac{\left<E^2\right>-\left<E\right>^2}{k_BT^2}
$$


### 3.1 模型一

​	模型一，主要构建的是如下类型的相互作用体系

![](F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model1\model.png)

<center>图2 模型一	

​	即一维Ising模型在上方(z=1)，XY模型在 z=0 的平面上，整个体系产生相互作用。能量表达式为
$$
H=−J_1\sum_{\left<i,j\right>}{\cos(\theta_i−\theta_j)}-J_2\sum_{\left<k,l\right>}\sigma_k\sigma_l-J_3\sum_{\left<u,v\right>}\cos\theta_u \sigma_v
$$
​	其中 $J_1,J_2,J_3$分别代表XY模型耦合常数，一维Ising模型耦合常数，XY模型与Ising模型相互作用常数；$i,j,u$粒子属于XY模型，$k,l,v$粒子属于Ising模型；XY模型只与周围4个点有相互作用，Ising模型只与其周围两个点有相互作用；两个模型之间，Ising模型的粒子只与其下方投影点的XY模型粒子起相互作用。

​	根据上述模型，设 $\mu_B=1，J_1=J_2=J_3=1,k_B=1$ 得到以下结果

​	![](F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model1\output.png)

<center>图3：模型一的能量，磁场，热容随温度变化图像
​	其中比热的临界温度大致在 $(1.1,1.2)$范围内，这与XY模型的临界温度接近，可以认为这类模型中XY模型由于粒子数众多，占据了主导作用，导致最终相变点任然由XY模型决定。



### 3.2 模型二

​	模型二，主要构建的是如下类型的相互作用体系

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model2\model2.png" style="zoom:67%;" />

<center>图4:模型二	


​	即一维Ising模型在 XY模型中轴上，整个体系产生相互作用，该模型主要模拟的是XY模型中参入一维Ising模型杂质后的表现。能量表达式为
$$
H=−J_1\sum_{\left<i,j\right>}{\cos(\theta_i−\theta_j)}
$$
​	其中，$i=width/2$处的粒子 $\theta\in\{0,\pi\}$，其他的地方的粒子 $\theta\in[0,2\pi)$

​	根据上述模型得到以下结果

![](F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model2\output.png)

<center>图5：模型二能量,热容随温度变化图像

​	热容的相变点依然和XY模型一致.但涡旋图像发生改变，如下图:

<center class="half">
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model2\T_0.01.png" style="zoom:67%;" />
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\XY_model\T_0.01_XY.png" style="zoom: 67%;" />
<center> 图6：左：模型二极速冷却后旋度图；右：XY模型极速冷却后旋度图
​	由于XY模型与一维Ising模型的相互作用，在低温下，Ising模型不再只朝一个方向。同时与XY模型在相同温度下对比可以看到，模型二有更多的的涡旋数目。

<div STYLE="page-break-after: always;"></div>

### 3.3 模型三

​	模型三，主要构建的是如下类型的相互作用体系

<center class="half">
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\model3_cube.png"" style="zoom:70%;" />
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\model3_plane.png"" style="zoom:70%;" />
<center><p>图7 左：模型三的3维图像；右：模型三的2维平面图

​	该模型是在XY模型的基础上，为平面上的粒子增加了Z轴的维度，使得其不仅可以在平面转动，还可以上下翻转，这更加接近电子自旋的实际状态。在该模型中，粒子自旋方向实际上是与Z轴呈 $\pi/4$ 或 $3\pi/4$的旋转的，二维平面图更加直观描述了该模型的样子。该模型实际上属于是二维Ising模型与XY模型的结合。

​	体系能量表达式为：
$$
H=−\frac J{\sqrt2}\sum_{\left<i,j\right>}[{\cos(\theta_i−\theta_j)}+\sigma_i\sigma_j]
$$
​	根据上述模型得到以下结果：

![](F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\output.png)

<center>图9：模型三的能量，磁场，热容随温度变化图像
​	从该热容随温度变化的图像可以看到，存在着两个尖峰，分别是XY模型相变的临界点 $T_1$ ，和二维Ising模型的临界点 $T_2$ 。从磁场变化看出，在 $T_2$临界点处，磁场变化最明显，可以认为是相变的温度。

<div STYLE="page-break-after: always;"></div>

​	从旋度图可以对比一下温度变化的过程：

<center class="half">
    <img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\T_0.010_model3.png"width=250/>
    <img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\T_1.100_model3.png"width=250/>
</center>
<center class="half">
    <img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\T_2.300_model3.png" width=250/>
    <img src="F:\MyDocuments\Physics\Computational Physics\Homework\BigHomework\model3\T_4.000_model3.png" width=250/>
</center>
<center><p>图10 模型三中不同温度下旋度的对比

​	当温度 $T< T_1$时粒子有旋度，局部也会有Z轴翻转；当温度 $T\approx T_1$时旋度较多，开始出现混乱，但Z轴方向依然较稳定；当温度 $T\approx T_2$时，粒子在XY平面方向上几乎混乱，而Z轴出现相变，开始有大面积翻转现象; 当 $T\gg T_2$后粒子在Z轴上也表现出混乱状态，整体为随机状态。

## 4 Summary

​	本次实验主要是探究XY模型+Ising模型的相互作用体系的相变，由于没有规定模型具体耦合形式，故探究了三种类型的模型。

1. 第一种模型主要是探究二者相互作用，即一个模型激发出磁场对另一个模型产生影响，彼此相互作用最终稳定。从结果可以看到，体系的相变主要还是取决于XY模型的相变，因为其维数更大，占据了主导作用。
2. 第二种模型主要是探究XY模型中参入一维Ising模型杂质后的表现。从结果可以看出参入Ising模型后，在低温下XY模型的涡旋数目明显增多，但比热容达到峰值的位置保持不变。
3. 第三种模型主要是探究XY模型+2维Ising模型相互作用体系，更加真实的模拟电子自旋的物理图像。实验结果展示出，在低温下，由于XY模型的作用，系统局部会发生Z轴翻转的现象；而从磁场变化图像中看到，2维Ising模型对磁场变化起到了相变主导的作用；比热容随温度变化图像有两个尖峰，分别对于着XY模型和2维Ising模型的临界温度。
