# Report15

## 1 Question

​	设体系的能量为 $H(x,y) = -2(x^2+y^2)+\frac12(x^4+y^4)+\frac12(x-y)^4$，取$\beta=0.2，1， 5$，采用 Metropolis 抽样法计算$\left< x^2 \right> ,\left< y^2 \right> , \left< x^2 +y^2\right>$。抽样时在 2 维平面 上依次标出 Markov 链点分布，从而形象地理解 Markov 链

## 2 Analysis

### 2.1 概率密度分布

​	根据玻尔兹曼分布可知，体系的概率分布函数正比于
$$
F(x,y) = \exp(-\frac H{kT})=\exp\{-\beta\times[-2(x^2+y^2)+\frac12(x^4+y^4)+\frac12(x-y)^4]\}\tag1
$$
其中 $\beta=\frac{1}{KT}$。

​	对 $F(x,y)$求偏导，求出其极值点，得到以下几组解
$$
\{x\to 0,y\to 0\},\left\{x\to -\sqrt{2},y\to -\sqrt{2}\right\},\left\{x\to \frac{\sqrt{2}}{3},y\to -\frac{\sqrt{2}}{3}\right\},\left\{x\to -\frac{\sqrt{2}}{3},y\to \frac{\sqrt{2}}{3}\right\},\left\{x\to \sqrt{2},y\to \sqrt{2}\right\}
$$
​	从解的结果来看，与 $\beta$ 的取值无关，在固定的点处取到极值，为了使得结果更加形象，不妨取 $\beta = 1$，得到

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\Fxy_beta1.png" style="zoom: 33%;" />

<center> 图1：F(x,y)图像


### 2.2 建议分布 $T(x',y')$

​	如图1所示的图像，可以明显的看到 $(-\sqrt2,-\sqrt2)(\sqrt2,\sqrt2)$处有两组极大峰值，故选择下列正态分布作为建议分布来抽样，其中
$$
G_1(x,y)=\frac{1}{2\pi\sigma^2}\exp(-\frac{(x-\mu)^2+(y-\mu)^2}{2\sigma^2})\\
G_2(x,y)=\frac{1}{2\pi\sigma^2}\exp(-\frac{(x+\mu)^2+(y+\mu)^2}{2\sigma^2})
$$
​	上式中取 $\mu=\sqrt2,\sigma=\frac{1}{2\sqrt\beta}$

设建议分布 $T$为
$$
T(x\rarr x',y\rarr y')=T(x',y')=\frac12(G_1(x',y')+G_2(x',y'))
$$
则建议分布的图像如下所示

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\Gauss_beta1.png" style="zoom:33%;" />

<center> 图2：T(x,y)图像

### 2.3 算法

那么根据 Metropolis-Hasting 抽样方法可以得到以下算法：

设初始点为 $(x_0,y_0)=(0,0)$

1. 设$M$为Markov链集合 $\{(x_k,y_k)\}$，取最后一个点为 $(x_{n},y_n)$

2. 在二维正态分布 $T(x',y')$ 下抽样得到试探点 $(x',y')$

3. 根据Metropolis-Hasting方法，令
   $$
   r_1=\frac{p_jT_{ji}}{p_i T_{ij}}=\frac{F(x',y')T(x_n,y_n)}{F(x_n,y_n) T(x',y')}
   $$

4. 在[0,1]均匀分布中抽样出随机点$\xi$，若 $\xi<\min(1,r)$，取 $(x_{n+1},y_{n+1})=(x',y')$；若**$\xi>\min(1,r)$**，取 $(x_{n+1},y_{n+1})=(x_{n},y_{n})$

5. 将 $x_{n+1}$加入到 $M$集合，循环到步骤1，直到集合总数大于 N

​	当得到有N个点的Markov链集合后，将序列中热化过程的前M个构型舍去，那么对应的平均值如下所示
$$
\left< x^2 \right> =\frac1{N-M}\sum_{i=M+1}^N x_i^2\\
\left< y^2 \right> =\frac1{N-M}\sum_{i=M+1}^N y_i^2\\
\left< x^2 +y^2\right>=\frac1{N-M}\sum_{i=M+1}^N x_i^2+y_i^2
$$

### 2.4 Markov Chain 形象化展示算法

​	为了绘制出 Markov Chain 的轨迹，上述建议分布 $T(x,y)$ 会在两个高斯分布间来回跳转，规避了粒子落入势阱中，无法逃离的情况，但是，无法展示出 Markov Chain 形象的过程，故在热化过程中，采取了另一种采样的方式，算法如下:

设初始点为 $(x_0,y_0)=(5,5)$

1. 设$B$为Markov链热化点集合 $\{(x_k,y_k)\}$，取最后一个点为 $(x_{n},y_n)$

2. 随机向前试探一步：
   $$
   (x_t,y_t)=(x_n+\xi_x*s,y_n+\xi_y*s)
   $$
   其中 $\xi_x,\xi_y$ 是[-1,1]上均匀分布分布的随机数，s=0.1为固定步长

3. 根据Metropolis抽样规则，设 能量改变为 $\Delta H$，选择概率为 $r_1$
   $$
   \Delta H = H(x_t,y_t)-H(x_n,y_n)\\
   r_1=\exp(-\beta\Delta H)\\
   $$

4. 在[0,1]均匀分布中抽样出随机点$\xi$，若 $\xi<\min(1,r)$，取 $(x_{n+1},y_{n+1})=(x_t,y_t)$；若**$\xi>\min(1,r)$**，取 $(x_{n+1},y_{n+1})=(x_{n},y_{n})$

5. 将 $x_{n+1}$加入到 $B$集合，循环到步骤1，直到集合总数大于 M



## 3 Experiment

### 3.1 $\beta=0.2$

​	数值计算可得以下结果
$$
\left<x^2\right> = 1.556912649927367,\\ \left<y^2\right> = 1.559989388948787, \\\left<x^2+y^2\right> = 3.116902038876154
$$
​	理论值计算是：$< x^2 > = 1.56189， < y^2 >= 1.56189$，结果比较准确

​	此时对应的是高温下的情况。由图 3 可见，分布的极大值处于原点附近。，但是原点和离原点较远的地方分布都比较小，抽样结果来看，基本与原分布吻合

<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\Boltzmann_0.2.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\sample_0.2.png" style="zoom:33%;" />
<center> 图3：分布与抽样图像比较

​	根据 $2.4$ 其热化得到的 Markov Chain如下所示

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\burn_points_0.2.png" style="zoom:50%;" />

<center> 图4：热化时的Markov Chain 

## 3.1 $\beta=1$

​	数值计算可得以下结果
$$
\left<x^2\right> = 1.6822846190170353,\\ \left<y^2\right> = 1.679052276012874, \\\left<x^2+y^2\right> = 3.36133689502991
$$
​	理论计算值是：$< x^2 > = 1.68247， < y^2 >= 1.68247$，结果比较准确

​	对应的图像对比与Markov Chain如图5，6所示

<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\Boltzmann_1.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\sample_1.png" style="zoom:33%;" />
<center> 图5：分布与抽样图像比较


<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15\burn_points_1.png" style="zoom:50%;" />

<center> 图6：热化时的Markov Chain 

