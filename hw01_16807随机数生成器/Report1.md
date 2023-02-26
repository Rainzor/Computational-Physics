# Report1

> Rainzor
>

### 1.Question

​	用Schrage方法编写随机数子程序，用指定间隔（非连续 l >1） 两个随机数作为点的坐标值绘出若干点的平面分布图。再用 $x^k$ 测试均匀性（取不同量级的N值，讨论偏差与N的关系）、 $C(l)$ 测试其2维独立性（总点数$N > 10^7$）

### 2.Method

#### 2.1 随机数的生成

​	Schrage 方法主要解决的是线性同余法产生随机数时可能出现的越界问题。其递推关系为:
$$
I_{n+1}=\left\{
\begin{aligned}
&a(I_n~mod~q)-r[I_n/q],&if~I_{n+1}>0\\
&a(I_n~mod~q)-r[I_n/q]+m,&otherwise
\end{aligned}
\right.
$$
​	实际计算中取q=127773， r=2836，a=16807，m=2147483647。当提供了第一个小于 m 的初始值 $I_0$ 之后，就可以 根据次递推公式得到后续的伪随机数。

​	代码中用 `Class Schrage16807`封装实现

#### 2.2 种子生成

​	采用计算机中的时间，来生成初始值
$$
I_0=(i_y-2000)+70(i_m+12\{i_d+31[i_h+23(i_n+59i_s)]\})
$$
其值所在区间为 $[0,2^{31}-1]$

​	代码中用 `seed_time()`函数封装实现



### 3.Experiment

#### 3.1 随机数平面分布图

​	采用相邻间隔随机数 $(l_n,l_{n+l})$ 生成二维平面,取5000个数,画出散点图.

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw01_16807随机数生成器\随机数分布图.png"
	style="zoom:67%;" />

<center><p>图1：随机数分布</p></center>

​	图中取间隔为2,由图可见未有散点存在明显的分层分区.可见在此检验下,Schrage随机数生成器不具有明显相关性,性能较好.

#### 3.2 检验随机数均匀性

​	本实验中,我采取了两种方式检验随机性:关于$<x^k>$k阶矩与期望值的比较以及卡方分布检验

##### 3.2.1 k阶矩检验

​	随机数均匀性测试结果如图 2 所示。N表示随机数个数，k表示k阶矩，下图展现一部分结果

​	其期望值满足
$$
E(x^k)=\int_0^1x^kdx=\frac{1}{1+k}
$$
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw01_16807随机数生成器\均匀性检验程序输出.png" style="zoom:50%;" />

<center><p>图2：不同N，k下，实验值与理论值结果


​	总体上可见随着 N 阶 数的增加，各 k 阶矩的实验值越来越接近理论值。可见 Schrage 方法得到的随机数均匀性不错。

​	以k=5为例，如图3所示，样本矩和理想值的偏差大致满足满足 $O（1/\sqrt{N}$）的关系

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw01_16807随机数生成器\均匀性.png" style="zoom:50%;" />

<center><p>图3：5阶矩实验偏差与样本数N的关系</p></center>

##### 3.2.2 卡方检验

​	将区间[0,1]  分为 K 个子区间，统计随机数落在第k 个子区间的实际频数 $n_k$ ， 它应当趋近于理论频数$m_k=N/K,(k=1,....,K)$,令统计量为
$$
\chi^2=\sum_{k=1}^{K}\frac{(n_k-m_k)^2}{m_k}
$$
​	$\chi^2$服从自由度为K-1的卡方分布，在显著水平 $\alpha=0.05$的条件下，对不同k，与N进行检验，得到图4结果（下图仅给出 $N\ge10^5$的结果 ），N是随机数个数，k为子区间个数，statistics为检验统计量，Percent_point为分布点

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw01_16807随机数生成器\卡方分布检验图.png" style="zoom:50%;" />

<center><p>图4：卡方分布检验</p></center>

​	由图可见，在置信度 $1-\alpha=0.95$的条件下，满足随机数满足均匀性，性能很好

#### 3.3 二维独立性检验

​	根据公式，得到检验统计量：
$$
C(l)=\frac{<x_nx_{n+l}>-<x_n>^2}{<x^2>-<x_n>^2}
$$
​	代码实现在 `covariance2`函数中。

​	实验当中 N 取 $2\times10^7$，在 l 处于 1 到 9 时，可见 C(l)接近于 0。即数据的二维独立性非常高，可见 Schrage 方法得到的 数据独立性不错。

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw01_16807随机数生成器\关联性检验程序输出.png" style="zoom:80%;" />

### 4 Summary

​	本实验掌握了随机数的产生方式和判断随机数产生器优劣的标准。

​	了解用 schrage 方法编写的 16807 产生器产生的随机数具有很弱的相关性以及很好的均匀性。
