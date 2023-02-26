# Report15

> Rainzor

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

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\Fxy_beta1.png" style="zoom: 33%;" />

<center> 图1：F(x,y)图像


### 2.2 建议分布 $T(x',y')$

​	如**图1**所示的图像，可以明显的看到 $(-\sqrt2,-\sqrt2)(\sqrt2,\sqrt2)$处有两组极大峰值，故选择下列**正态分布作为建议分布**来抽样，其中
$$
G_1(x,y)=\frac{1}{2\pi\sigma^2}\exp(-\frac{(x-\mu)^2+(y-\mu)^2}{2\sigma^2})\\
G_2(x,y)=\frac{1}{2\pi\sigma^2}\exp(-\frac{(x+\mu)^2+(y+\mu)^2}{2\sigma^2})
$$
​	上式中取 $\mu=\sqrt2,\sigma=\frac{1}{2\sqrt\beta}$

​	$\mu$是固定在极大值处的，而$\sigma$ 与 $\beta$相关的。这是由于随着 $\beta$增大，极值点不改变，而待抽样的分布会逐渐收敛到极值点附近，相应**图1**中极值点附近的弥散度下降，峰变得尖锐，这可以在后续报告中的 **图3，图5，图7**中可以看出，故**随着 $\beta$增大，标准差$\sigma$要逐渐减小**

设建议分布 $T$为
$$
T(x\rarr x',y\rarr y')=T(x',y')=\frac12(G_1(x',y')+G_2(x',y'))
$$
则建议分布的图像如下所示<a id="ch2.3"> </a>

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\Txy_beta1.png" style="zoom:33%;" />

<center> 图2：T(x,y)图像

### 2.3 抽样算法

那么根据 **Metropolis-Hasting 抽样方法**可以得到以下算法：

设初始点为 $(x_0,y_0)=(0,0)$

1. 设$M$为Markov链集合 $\{(x_k,y_k)\}$，取最后一个点为 $(x_{n},y_n)$

2. 在二维正态分布 $T(x',y')$ 下抽样得到试探点 $(x',y')$

3. 根据Metropolis-Hasting方法，令
   $$
   r_1=\frac{p_jT_{ji}}{p_i T_{ij}}=\frac{F(x',y')T(x_n,y_n)}{F(x_n,y_n) T(x',y')}
   $$

4. 在[0,1]均匀分布中抽样出随机点$\xi$，若 $\xi<\min(1,r)$，取 $(x_{n+1},y_{n+1})=(x',y')$；若**$\xi>\min(1,r)$**，取 $(x_{n+1},y_{n+1})=(x_{n},y_{n})$

5. 将 $x_{n+1}$加入到 $M$集合，循环到步骤1，直到集合总数大于 N 

​	当得到有N个点的Markov链集合后，将序列中热化过程的前M个构型舍去，那么对应的平均值如下所示<a id="ch2.4"> </a>
$$
\left< x^2 \right> =\frac1{N-M}\sum_{i=M+1}^N x_i^2\\
\left< y^2 \right> =\frac1{N-M}\sum_{i=M+1}^N y_i^2\\
\left< x^2 +y^2\right>=\frac1{N-M}\sum_{i=M+1}^N x_i^2+y_i^2
$$

### 2.4 Markov Chain 形象化展示算法

​	为了绘制出 Markov Chain 的轨迹，上述建议分布 $T(x,y)$ 会在两个高斯分布间来回跳转，规避了粒子落入势阱中，无法逃离的情况，但是，无法展示出 Markov Chain 形象的过程，故**在热化过程中，采取了另一种采样的方式**，算法如下:

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

​	此时对应的是高温下的情况。

​	数值计算可得以下结果
$$
\left<x^2\right> = 1.556912649927367,\quad \left<y^2\right> = 1.559989388948787, \\\left<x^2+y^2\right> = 3.116902038876154
$$
​	理论值计算是：$< x^2 > = 1.56189， < y^2 >= 1.56189$，结果比较准确

​	由**图 3** 可见，分布的极大值处于原点附近，但是原点和离原点较远的地方分布都比较小，抽样结果来看，基本与原分布吻合

<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\Boltzmann_0.2.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\sample_0.2.png" style="zoom:30%;" />
<center> 图3：高温分布与抽样图像比较


​	根据<a href="#ch2.4">2.4</a> 其热化得到的 Markov Chain如**图4**所示，粒子较为均匀在两个势阱间穿梭<a id="graph4"> </a>

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\burn_points_0.2.png" style="zoom: 33%;" />

<center> 图4：高温热化时的Markov Chain 

### 3.2 $\beta=1$

​	此时对应的是中等温度下的情况。当β较大时马尔科夫链开始向分布极值靠近。

​	数值计算可得以下结果
$$
\left<x^2\right> = 1.6822846190170353,\quad \left<y^2\right> = 1.679052276012874, \\\left<x^2+y^2\right> = 3.36133689502991
$$
​	理论计算值是：$< x^2 > = 1.68247， < y^2 >= 1.68247$，结果比较准确

​	对应的图像对比与Markov Chain如**图5，6**所示

<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\Boltzmann_1.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\sample_1.png" style="zoom:27%;" />
<center> 图5：中等温度分布与抽样图像比较


<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\burn_points_1.png" style="zoom: 33%;" />

<center> 图6：中等温度热化时的Markov Chain 

### 3.3 $\beta=5$

​	此时对应的是低温的情况。

​	数值计算可得以下结果
$$
\left<x^2\right> = 1.9464660730764258,\quad \left<y^2\right> = 1.950173563246566, \\\left<x^2+y^2\right> = 3.896639636322991
$$
​	理论计算值是：$< x^2 > = 1.94783， < y^2 >= 1.94783$，结果比较准确

​	对应的分布与抽样图像对比如**图7**所示
<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\Boltzmann_5.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\sample_5.png" style="zoom:27%;" />
<center> 图7：低温分布与抽样图像比较

​	为比较<a href="#ch2.3">2.3</a>与<a href="#ch2.4">2.4</a>提出的算法，取在 $\beta=5$时比较更为形象，因为此时在极值点处是势阱很深，如**图7**所示，粒子基本只在极值点附近，而极少的有“隧穿”过程，比较结果如下所示

<center><half>
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\burn_points_5.png" style="zoom:33%;" />
	<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw15_Markov链\mcmc_5.png" style="zoom:33%;" />
<center> 图8：低温热化与抽样时的Markov Chain 对比

​	当采取<a href="#ch2.4">2.4</a>中的算法，如**图8左**所示，在从（5，5）出发的点，当靠近极值点后，基本活动范围只在第一象限内的势阱中，几乎不可能到另一边的势阱，这与<a href="#graph4">图4</a>展现的 $\beta$ 较小（高温）时热化的Markov Chain很不一样。

​	而抽样采取 <a href="#ch2.3">2.3</a> 中提出的两个高斯分布作为建议抽样分布，可以规避这类情况的发生，如**图8右**所示，粒子在两个势阱间来回穿梭，不会只局限在一个势阱中，可以遍历全局，抽样出完整的分布。

## 4 Summary

- 利用 Metropolis 方法对玻尔兹曼分布积分。当温度较高的时候马克科夫链会同时 往远离原点和靠近原点的方向扩散，但是由于离原点较远的点对$< x^2 >$贡献远大于靠近 原点的点，所以积分值会偏高，随着温度降低，马尔科夫链将收拢与分布极值点，使得积分收 敛于一个稳定值。
- 比较了两种 Metropolis 抽样方法，对于普通的Metropolis 方法，粒子在遇到较大势阱时，很容易在一个势阱中无法逃脱。根据Metropolis-Hasting 方法，应当采取与待抽样分布近似的建议分布，使得让粒子有机会挣脱势阱束缚，在全局进行遍历，这样才能采样出完整的分布。
- 在待抽样的Boltzmann分布的参数 $\beta$ 在变化时，相应的建议分布 $T(x,y)$也要随之改变，在实验中发现，在调整标准差 $\sigma=\frac{1}{2\sqrt\beta}$ 时，建议分布于待抽样分布比较接近，抽样的结果较好。
