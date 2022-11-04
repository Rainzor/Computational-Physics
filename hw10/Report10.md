# Report10

> PB20020480 王润泽

## 1. Question

​	Monte Carlo 方法研究二维平面上荷电粒子在正弦外电场( ~ sin𝜔t )中的随机行走。 推导速度自相关函数的表达式。它随时间的变化是怎样的行为？能否模拟得到该自相 关函数的曲线?是的话与理论曲线进行比较，否的话讨论理由

## 2. Algorithm

#### 2.1 速度自相关函数

​	由于题目未给定电场方向，不妨假设电场沿着x方向，那么可以类比**Langevin方程**得到
$$
&\frac{dv_x}{dt}=-\frac1{\tau}v_x+A_x(t)+B\sin(\omega t) \tag{1}\\
$$

$$
\frac{dv_y}{dy}=-\frac1{\tau}v_y+A_y(t) \tag{2}
$$

​	与讲义一致$\tau = \frac{m}{6\pi\eta a} $  ，$\bold A$代表了一种涨落力，满足
$$
\left<\bold A(t)\right>=0\\
\left<\bold A(t)\cdot\bold A(0) \right>=D\delta(t)
$$
​	由（1）（2）式，解出微分方程
$$
v_x(t)=v_x(0)e^{-t/\tau}+e^{-t/\tau}\int^t_0e^{t'/\tau}\left[ A_x(t')+B\sin(\omega t')\right]dt'\\
v_y(t)=v_y(0)e^{-t/\tau}+e^{-t/\tau}\int^t_0 e^{t'/\tau} A_y(t')dt'
$$
**以下公式默认 ** $\left<\bold r\right>=\left<\bold v\right>=0$	

由于是二维平面，所以距离平方的期望是
$$
\left<\bold r^2(t)\right>=4Dt
$$
这样根据定义，二维速度自相关函数为
$$
C(t)=\frac12\left<\bold v(t)\cdot \bold v(0)\right>=\frac12\left<v_x(t)v_x(0)\right>+\frac12\left<v_y(t)v_y(0)\right>\\
\left<v_x(t)v_x(0)\right>=\left< v_x^2(0)\right>e^{-t/\tau}+e^{-t/\tau}\int^t_0e^{t'/\tau}\left[ \left<v_x(0)A_x(t')\right>+\left<v_x(0)B\sin(\omega t')\right>\right]dt'\\
\left<v_y(t)v_y(0)\right>=\left< v_y^2(0)\right>e^{-t/\tau}+e^{-t/\tau}\int^t_0\left[e^{t'/\tau} \left<v_y(0)A_y(t')\right>\right]dt'
$$
​	很显然，$\bold A$作为随机力，与速度 $\bold v$无关，同时 $\left<A\right>=0$

​	所以可以推导出二维速度自相关函数的表达式为：
$$
\begin{aligned}
C(t)&=\frac12\left<v_x(t)v_x(0)\right>+\frac12\left<v_y(t)v_y(0)\right>\\
&=\frac12(\left< v_x^2(0)\right>+\left< v_y^2(0)\right>)e^{-t/\tau}\\
&=\frac{\left< \bold v^2(0)\right>}2e^{-t/\tau}
\end{aligned}\tag{3}
$$
​	注：$v_x(0)$与 $v_y(0)$在初态应该属于相同的均匀分布

​	**速度自相关系数的变化符合e指数衰减，与外界作用力本身无关。**

#### 2.2 模拟验证的算法

##### 2.2.1数值公式

​	为了验证理论的可靠性，得用数值计算模拟出粒子速度随时间的方程 $\bold v(t)$

​	仍然根据（1）（2）所得的**Langevin方程**
$$
\frac{dv_x}{dt}=-\frac1{\tau}v_x+A_x(t)+B\sin(\omega t)
$$

$$
\frac{dv_y}{dy}=-\frac1{\tau}v_y+A_y(t) 
$$

​	为了计算微分方程，采取**微分近似为差分**的思想，采用 **Euler-梯形公式** 得，其中时间间隔 $h=t_{n+1}-t_n$
$$
v_{x,n+1}=v_{x,n}+{h}\left[-\frac1{\tau}v_{x,n}+A_x(t_n)+B\sin(\omega t_n)\right]\tag{4}\\
$$

$$
v_{y,n+1}=v_{y,n}+{h}\left[-\frac1{\tau}v_{y,n}+A_x(t_n)\right]\tag{5}
$$

​	得到 数值方程解：
$$
v_{x,n+1}=(-\frac{h}{\tau}+1)v_{x,n}+{h}\left[A_x(t_n)+B\sin(\omega t_n)\right]\tag{6}\\
$$

$$
v_{y,n+1}=(-\frac{h}{\tau}+1)v_{y,n}+{h}\left[A_y(t_n)\right]\tag{7}
$$

##### 2.2.2 参数选择

1. 对于 $\tau= \frac{m}{6\pi\eta a}=\frac{mD}{kT}$，由于不知道粘滞阻力等具体数值，编程中设置为5，保证粒子的一定随机性，而不至于导致速度快速减小到0附近，只有电场的周期震荡。	

$$
\tau=5
$$
2. 对于指数衰减来说，总时间取 $T=10\tau$较为合适，此时指数已衰减到很小

3. 对于布朗运动来说，无论时间间隔 $h$取的多么小，都不会影响整体运动轨迹，不妨为了方便采样取

$$
h=\tau/1000
$$
4. 对于电磁场周期频率，为了让实验符合正常观察情况，至少应该确保在**总时间T**内有多个个周期的电磁场作用力被观测到，取

$$
\omega = 2\pi
$$
5. 为了保证结果的准确性，在数值计算平均值式，应当取以在相同时刻t取100次，以平均作为检验统计量来表示 $\left<\bold v(t)\cdot \bold v(0)\right>$

6. 在实际情况中，一般**随机力应该小于电磁作用力**，故我选择以下表达式 $B=20A$，然而由于统计检验精度有限，对于 $\bold v(0)\in[-1,1]$内，$<\bold v(0)>\approx0.01$，故选择作用力时大小范围应当合适，所以有

$$
A_{max}=0.05\\
B\sin(wt)=\sin(2\pi t)
$$

##### 2.2.3 流程图

```flow
st=>start: 初始化速度与位置,随机力
 
op=>operation: 更新速度与位置
 
cond=>condition: 迭代是否结束(是或否?)
 
sub1=>subroutine: 产生随机力

io=>inputoutput: 计算速度关联函数

plt=>operation: 画图

e=>end: 结束
 
st->op->cond
 
cond(yes)->io->plt->e
 
cond(no)->sub1(right)->op
```

## 3. Experiment

#### 3.1粒子随机游走图像

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw10\random_walk.png" style="zoom:50%;" />

<center><p>图1：在强电场中粒子随机游走
得到如上图像
**解释**：一开始速度比较大，而随机力比较小，所以随机游走不明显，呈直线行走状态；但是随着时间增长，由于阻力的作用的效果，使得速度逐渐降下来，导致粒子做衰减的简谐振动，同时Y方向仍然有一定的随机性。

#### 3.2 速度相关函数

按照之前的推导得到的(3)式，为理论的速度相关系数
$$
C_1(t)=\frac{\left< \bold v^2(0)\right>}2e^{-t/\tau}\tag{8}
$$
实验模拟所得速度相关函数，采用速度相关系数的定义式
$$
C_2(t)=\frac12\left<\bold v(t)\cdot \bold v(0)\right>=\frac12\left<v_x(t)v_x(0)\right>+\frac12\left<v_y(t)v_y(0)\right>\tag{9}
$$
实验与理论比较得下图

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw10\covariance_v.png" style="zoom: 50%;" />

<center><p>图2：在强电场中，速度相关系数理论与模拟图像

​	图中尾处有一定的波动情况，这是由于初始速度平均值 $<\bold v(0)>$ 并不精准为0导致的。但在误差涨落允许范围内可以看到，图像的整体趋势实验模拟与理论吻合较好，即速度自相关函数确实为

$$
C(t)=\frac{\left< \bold v^2(0)\right>}2e^{-t/\tau}
$$


#### 3.3 调整随机力大小和粘滞阻力

加大粘滞阻力与随机力
$$
\tau=1\\
A_{max}=10\\
B\sin(wt)=\sin(2\pi t)
$$

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw10\random_walk2.png" style="zoom:50%;" />

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw10\covariance_v2.png" style="zoom:50%;" />





​	由图像可见，当加大随机力后，图像更加接近随机游走，局部有微小的振动；而由于粘滞阻力的作用，粒子相关系数仍然会逐渐收敛到0，但由于实验模拟精度有限，右图像末端会有一些涨落，但整体趋势仍然符合理论。



### 4. Summary

​	实验模拟外加周期力的情况下粒子的随机行走现象。

​	实验可见当外加电场相对大小较大的时候，并且粒子有向着初始方向漂移的趋势。直到因为粘滞阻力效果，粒子最终在某处沿x方向来回周期运动，但在y方向，仍然保持一定的随机性。

​	在外加电场相对大小较小的情况下粒子主要体现随机运动，在局部有微小的振动。

​	无论是较小的电场还是较大的电场，随机运动模型下粒子的自相关函数和粒子所受的外力无关，粒子的自相关函数只取决于粒子所受的衰减力的大小。