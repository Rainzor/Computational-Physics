# CH1 Monte Carlo方法基础

**蒙特卡罗方法**（英语：Monte Carlo method），也称**统计模拟方法**，是1940年代中期由于科学技术的发展和电子计算机的发明，而提出的一种**以概率统计理论为指导的数值计算方法**。是指使用随机数（或更常见的伪随机数）来解决很多计算问题的方法

## 1 随机数产生器

### 1.1 均匀随机数产生器

#### a.随机数要求

- 伪随机数：效率高，有周期，适用于实验数值模拟

- 真随机数：效率低，无周期，应用于游戏等领域

**优秀的伪随机数产生器**

- 随机数序列周期长

- 随机数之间没有很大的关联性

- 随机数分布均匀

- 随机数产生快


随机数分布函数
$$
p(x)=\left \{
\begin{aligned}
&1&x\in[0,1]\\
&0&x\notin[0,1]
\end{aligned}
\right.
$$


#### b. 随机数生成方法

##### Lehmer 线性同余

$$
I_{n+1} = (aI_n+b) \text{mod}m\\
x_n= I_n/m
$$

##### 16807产生器

线性同余法的特例

$$
a =7^5=16807,b=0,m=2^31-1=2147483647
$$

理论与实际效果来看，$b=0$和$b\neq0$效果差不多

##### **Schrage方法**：

由于32位计算机对于超过$2^{31}-1$的数无法计算，而$aI_n$是有可能超过的，所以选择合适的计算方法很重要
$$
m=2147483647=16807\times q+r\\
q=127773\\
r=2836\\
aI_n\,\text{mod}\,m=\left\{
\begin{aligned}
&a(z\,\text{mod}\,q)-r[z/q],&if\geq0 \\
&a(z\,\text{mod}\,q)-r[z/q]+m,&otherwise\\
\end{aligned}
\right.
$$

##### 种子值

根据电脑内年月日时分秒时间来设置种子值
$$
I_n=i_y+70\{i_m+12[i_d+31(i_h+23(i_n+59i_s))]\}
$$
范围是：$[0,2^{31}-1]$



##### 位移计数器法

首先用其它方法产生一个随机的整数序列，然后对两个整数进行 XOR（与或）操作以 产生一个新的随机整数：
$$
I_n=I_{n-p}\oplus I_q
$$
其中的[ ] p q, 是一对整数，最佳的选择是满足条件
$$
p^2+q^2+1=prime
$$
如：[31，3]、[98，27]、[250，103]、[1279，(216,418)]、[9689， (84,471,1836,2444,4187)]等



其中的 R250 ( p = 250 , q = 103) 是最为常用的产生 器。一般来说[p, q] , 值越大，产生的随机数质量越好，而且起始的随机整数表的质量也很重要，因此常用同余法产生有 p 个数的初始值表

##### Fibonacci 延迟产生器

位移计数器法是更一般的 Fibonacci 延迟产生器的一个特例

Fibonacci延迟产生器思想是：
$$
I_n=I_{n-p}\otimes I_{n-q}\mod{m}
$$
其操作符 $\otimes$可以是：加、减、乘、XOR，整数对[ p q, ] 表示延迟。



Fibonacci 延迟产生器较线性同余法的优势在于它的周期非常长，32 位机上的最大周期为$(2^p-1)2^{32}$



#### c.  伪随机数的统计检验

##### 独立性检验

##### 均匀性检验



## 2 已知分布的随机抽样

从均匀分布到特定概率分布

### 2.1 直接抽样法

#### a.离散型变量分布

设变量 x 是离散型的，取值为 $x_1,x_2,...$，相应值出现的几率为 $P_1,P_2,...$

如果从[0,1]区间中均匀抽样得到的随机数$\xi$满足下式时
$$
\sum_{i=1}^{n-1}P_i<\xi\leq\sum_{i=1}^{n}P_i
$$
则物理量 x 取值为$x_n$

#### b.连续型变量分布

思想：连续问题离散化
$$
p(x)=\sum_i p_i\delta(x-x_i)\\
P(x,x+dx)=p(x)dx\\
\int p(x)\text{d}x=1
$$
==注意！！由上式可见$p(x)$量纲是$[x]^{-1}$！！！！==



对连续随机变量$x\in[a,b]$的累计函数：
$$
\xi(x)=\int_a^xp(t)\text{d}t
$$
性质：$\xi(a)=0,\xi(b)=1,\xi(x)$单调递增，无量纲

所以每一个$x对应唯一的\xi(x)$,因此可求其反函数$x(\xi)$,通过插值，对[0，1]区间内随机数$\xi $，可以唯一确定x的值，x的概率分布满足$p(x)$



- **例：**粒子随机运动的自由程分布为指数分布

$$
p(x)=\left\{
\begin{aligned}

&\lambda^{-1}e^{-x/\lambda},&x>0\\
&0,&otherwise
\end{aligned}
\right.
$$
其中的$\lambda$为运动平均自由程
$$
\xi =\int_0^x\lambda{-1}e^{-t/\lambda}=1-e^{-x/\lambda}
$$
那么可以求反函数得
$$
x=-\lambda\ln{(1-\xi)}
$$
注意到$\xi$与$1-\xi$都是均匀随机变量，那么对采样结果不产生影响，所以
$$
x=-\lambda\ln{\xi}
$$
类似的
$$
p(x) = \pi^{-1}(1-x^2)^{-1/2},(-1\leq x\leq1)\\
x=\sin{(\pi \xi -\pi/2)}\iff \sin{2\pi \xi}\iff \cos{2\pi \xi}
$$
**直接抽样法缺陷：**对于复杂的概率分布，积分与反函数很难用数学形式表达，比如经典的Gauss函数



### 2.2 变换抽样法

#### a.一般方法

- **思想**：将一个复杂分布$p(x)$的抽样转换为已知简单分布$g(y)$，比如[0, 1]区间中的均匀分布。那么只要采样简单分布$g(y)$，通过下式，即可立马得到$p(x)$

$$
p(x)dx=g(y)dy \Rightarrow p(x)=\left |\frac{dy}{dx}\right|g(y)
$$
当g(y)取均匀随机分布时，问题转换为：

寻找$y(x)$，使其导数为$p(x)$，然后在[0, 1]  区间中对变量 y 抽样得到均匀分布的随机数，再由 $x ( y) $关系得到对应几率密度函数 $p(x)$的随机抽样 x 。

- **一维变换法：**a.累积函数本身就是一维变换抽样法的一类特殊情况，只要设$\xi(x)=y(x)$，那么采样$\xi$并通过变换关系，就可得到了x的分布。

这样根据 $p(x)= \left|\frac{d\xi}{dx}\right|$，可以反解出 $x(\xi)$关系式，即可使用直接抽样法得到抽样分布。 **然而这对于复杂分布仍然不可解**

- **多维变换法：**有两个变量 x 和 y 的联合分布密度函数$p(x)$，欲变换至变量u 和v ，它们的联合分布密度函数为$g(u,v)$，则推广方程：

$$
p(x,y)dxdy=g(u,v)dudv=g(u(x,y),v(x,y))\left| \frac{\part(u,v)}{\part(x,y)}\right|dxdy
$$
取联合分布密度函数$g(u,v)$为均匀分布，则只要满足$p(x,y)=\left| \frac{\part(u,v)}{\part(x,y)}\right|$

问题转换为对均匀随机变量(u, v)进行抽样，代入$x=x(u,v),y=y(u,v)$,得到x，y的抽样

#### b. Box-Muller法

对Guass正态分布，考虑标准型
$$
p(x)=\frac{1}{\sqrt{2\pi}}e^{-x^2/2}
$$
假设x,y独立，都是正态分布
$$
p(x,y) = p(x)p(y)=\frac{1}{{2\pi}}e^{-(x^2+y^2)/2}
$$
==那么只要找到均匀随机变量(u, v)的变换关系，那么抽样得到的(x，y)都满足正态分布==
$$
\begin{aligned}
p(x,y)dxdy&=p(r,\varphi)rdrd\varphi&(极坐标)\\
&=p(t,\varphi)d(t/2)d\varphi &(t=r^2)\\
\end{aligned}
$$
所以$p(t,\varphi)=\frac{1}{2}exp(-\frac{t}{2})\times(\frac{1}{2\pi})$，$t\ge0,\varphi \in[0,2\pi]$

则对上式抽样是简单已知的,采取直接抽样法，(u, v)是[0, 1]均匀分布

t是指数分布（$\lambda =2$）: $t = r^2 = x^2+y^2=-2\ln u$

$\varphi$是$[0,2\pi]$的均匀分布 : $\varphi =\arctan{y/x} =2\pi v$

以上，得到$(u, v)\iff(x,y)$的关系
$$
x=\sqrt{-2\ln u}\cos{2\pi v}, y =\sqrt{-2\ln{u}}\sin{2\pi v}\\
\left| \frac{\part(u,v)}{\part(x,y)}\right|=\frac{1}{{2\pi}}e^{-(x^2+y^2)/2}=p(x,y)
$$

#### c. Marsaglia 方法

- **背景**：应用中经常遇到求在圆上（二维）或球面（高维）上均匀分布的抽样问题，当然最简单的是用极坐标或球坐标对角度进行抽样，然后再用坐标变换变到直角坐标下。例如，二维时首先取极角 $\varphi \in(0, 2π )$ 的均匀分布抽样，再计算 $x =r\cosφ 和 y= r\sinφ$。但是，三角函数的计算耗时较多，一般不希望采用这样的抽样方式。

- **思想**：对在圆外不符合要求的点抛弃，重新取样，代替三角函数。

- **方法：**

（1）随机抽样一对均匀分布的随机数 $(u,v)\in[-1,1]$

（2）计算 $r^2=u^2+v^2$ ， 如果 $r^2>1$ 则重新抽样直至 $r^2\le1$；

（3）则 $\cos\varphi=u/r,\sin\varphi=v/r$，且$r^2$也符合均匀分布特点。

- **评价：** 按此方式需要 2 个均匀随机数抽样，而且在第二步中可能舍去不符要求的抽样，抽样效率为π /4。即使这样，计算的效率仍较计算三角函数的为高。

- **应用**：Gauss分布抽样

  由于 $u^2+v^2$服从[0，1]的均匀分布，且$(u/r,v/r)$是圆环上的均匀分布点，所以得
  $$
  x=\frac{u}{r}\sqrt{-2\ln r^2},y=\frac{v}{r}\sqrt{-2\ln r^2}
  $$
  **则x,y均符合高斯分布**

- **推广：**三维球面分布的Marsaglia 方法

（1）随机抽样一对均匀分布的随机数，$(u,v)\in[-1,1]$ ；

（2）计算 $r^2=u^2+v^2$ ， 如果 $r^2>1$ 则重新抽样直至 $r^2\le1$；

（3）得$x_1=2u\sqrt{1-r^2},y=2v\sqrt{1-r^2},z=1-2r^2$

- **推广：**4 维超球面上分布的 Marsaglia 方法

（1）随机抽样一对均匀分布的随机数，$(y_1,y_2)\in[-1,1]$，直至满足$r_1^2=y_1^2+y^2_2\le1$；

（2）随机抽样一对均匀分布的随机数，$(y_3,y_4)\in[-1,1]$，直至满足$r_2^2=y_3^2+y^2_4\le1$；

（3）得$x_1=y_1,x_2=y_2,x_3=(y_3/r_2)\sqrt{1-r_1^2},x_4=(y_4/r_2)\sqrt{1-r_1^2}$

### 2.3 舍选抽样法

#### a. 简单分布

- **定义**：密度分布函数 p(x)定义在有限区域[a, b]内且是有界的，M为上界，则
- **舍选方式：**设g(x，y)在$[a, b]\times[0,M]$内为均匀分布，且

（1）产生一对[0,1]区间中均匀分布的随机抽样值 $(\xi_1,\xi_2)$，由 g(x,y)得抽样表示式$\xi_1=(x-a)/(b-a),\xi_2=y/M$

（2）判断条件 $M\xi_2\le p(a+(b-a)\xi_1)$ 是否成立，否，则舍；

（3）是，则取 $x=a+(b-a)\xi_1$

- **解释：**x 的取值落在区间 ( x, x+dx ) 内的概率等于面积比

  <img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20220919204836803.png" alt="image-20220919204836803" style="zoom:50%;" />

- 缺陷：当曲线p(x)呈尖峰形状时， 抽样效率很低。这时需要把变换抽样与舍选法结合起来。

- 改进：将上面的 y = M直线改为一个形状已知且是可积分的函数F(x)， 曲线形状与 p(x) 类似但处处比 p(x)大： $F(x)>p(x)$ ， F(x)称为比较函数。在比较函数 的面积区内产生随机点(x,y) ，由反函数法推出$x=x(\xi_1)$。
  $$
  \xi_1=\frac{\int^{x}_aF(x)dx}{\int^b_aF(x)dx},y=\xi_2F(x)
  $$
  

  对于较大的 $F(x)$ ，由此得到的 y 轴上的各个 y 点分布较稀疏，但 x 轴方向上 x处附近点的分布较密，因此==单位面积内(x ,y ) 点的分布仍是均匀的==，且全部落 在 F (x ) 曲线下的面积内。

  如果该点在 p(x)的面积区内，即$y=\xi_2F[x(\xi_1)]< p[x(x_1)]$ ,则取$x$，否则舍去该点。

  比较函数的具体形式不影响抽样的准确性，但抽样效率即有效选 取的点数为 p 与 F 的面积比，与F(x)的形状有关。实际往往去F(x)为阶梯函数形式

  <img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20220919205237523.png" alt="image-20220919205237523" style="zoom:50%;" />

#### b. 一般形式 

- **背景**：采用直接抽样法和变换抽样法经常遇到很大困难，主要是各种解析表达式不 易给出，甚至密度分布函数本身就是以数值表的形式给出的。

- **思想方法：** $g(x,y)$为任意 2 维联合分布几率密度函数，$h(x)$ 是任意函数，则对于可以表示成如下积分形式的分布$p(x)$
  $$
  p(x)=\frac{\int^{h(x)}_{-\infin}g(x,y)dy}{\int^{+\infin}_{-\infin}dx\int^{h(x)}_{-\infin}g(x,y)dy}
  $$
  （1）由g(x,y)产生一对随机抽样值(x, y)，由直接抽样法规则，对均匀分布的 $(\xi_1,\xi_2)$,求出下式的反函数
  $$
  \xi_1(x)=\int^x_{-\infin}dx\int^{+\infin}_{-\infin}dyg(x,y),\xi_2(y)=\int^{+\infin}_{-\infin}dx\int^y_{-\infin}dxg(x,y)
  $$
  （2）判断条件 $y(\xi_2)\le h(x(\xi_1))$是否成立

  （3）否，则舍去；是，则取x为p(x)的随机抽样

- 实例：Marsaglia 方法

#### c. 乘分布

乘分布的一般形式是 $p(x)=h(x)q(x)$，其中：$\int q(x)dx=1$，$h(x)$的上界是M。有g(x)（即一般形式下的表达式g）

$$
g(x,y)=\left\{
\begin{aligned}
&q(x)/M,&if~~0\le y\le M\\
&0,& otherwise
\end{aligned}
\right.
$$
按照舍选方法的一般形式有：

1. 产生服从分布 $q(x)$的随机抽样 $x$，$\xi_1$是（0，1）的均匀分布，x满足

$$
\int^x_{-\infin}dx\int^{+\infin}_{-\infin}dy~g(x,y)=\int^x_{-\infin}q(x)dx=\xi_1
$$

2. 另外再产生一个[0, 1] 区间中均匀分布的随机抽样值 $\xi_2$ ，
   $$
   \xi_2(y)=\int^y_{-\infin}dy\int^{+\infin}_{-\infin}dx~~g(x,y)=y/M
   $$
   
2. 判断条件  $M\xi_2\le h(x)$ 是否成立。

2. 是，则取x；否，则舍去

**优点：**可以根据一定的概率密度选择抽样，这样提高抽样效率

## 3 定积分的计算

### 3.1 简单抽样

#### a.掷石法

​	用 Monte Carlo 方法计算定积分的原理，f(x)的定积分值即为曲线下的面积值，其中 $m\le f(x)\le M$
$$
S = \int^b_af(x)dx\approx S_0\times(n/N)
$$
​	N是随机数总数，均匀分布在 $[a,b]\times[m,M]$，n为在f(x)内的点数

#### b.平均值法

​	各种定积分的数值计算方法中，都要在按某种方式固定划分的网格上算出 $f(x_i)$ 值，而采用 Monte Carlo 方法计算时，$x_i$ 可以是在[a , b ]区间中均匀随机选取的。根据积分的平均值定理
$$
\int^b_af(x)dx=(b-a)\left<f\right>
$$
​	平均值又可从下式得到
$$
\left<f\right>\approx\frac{1}{N}\sum_i^Nf(x_i)
$$
​	联立即可得f的积分表达式

**缺陷：**在简单抽样中，我们由均匀分布中选取随机数，并不考虑被积函数的具体情况。因此，被积函数的 极大处和极小处有相同的抽样权重，而对积分贡献较大的更多在函数值较大处， 故直观上就可以看出，**非权重的简单抽样效率较低。**即为了获得较高计算结果的精度，需要大量的抽样。

#### c.大数定理与中心极限定理

大数定理指出：
$$
\lim_{N\rarr\infin}\frac{1}{N}\sum f_i\rarr\mu
$$
中心极限定理指出

​	当抽样N次后，对于x服从分布 $f(x)$，平均值为 $\mu$，方差为 $\sigma_f$，那么
$$
\frac{x-\mu}{\sigma_f/\sqrt{N}}\rarr N(0,1)
$$
​	这意味着，随机抽样得到的x的平均值$\left<f\right>$与真实值 $\mu$的偏差 $\sigma_s$随抽样点的变化为
$$
\sigma_s=\left|\left<f\right>-\mu\right|\propto\frac{\sigma_f}{\sqrt N}=\sqrt{\frac{{\left<f^2\right>-\left<f\right>^2}}{N}}
$$
​	上述体现了Monte Carlo法的**重要方面**：

- $\sigma_s$ 随 $1/\sqrt N$ 变化， 即当样本点增加100倍时误差缩小10倍，反过来说，要达到一定的计算精度，必须以平方的方式增加总样本点数，这是 Monte Carlo 计算的一个固有弱点，即和其它积分数值计算方法相比收敛速度慢。

  但这只限于低维情形，**在高维积分和奇异积分下，Monte Carlo 方法有巨大的优势**。

-  当 $\sigma_f$ 小即函数平坦时，计算精度可以提高。极限情况下，对于常数 f (x) ，只取一点就可得到积分值。另一极 限情况下，对于 $\delta(x)$ 函数，只有很少的样本点能被选中，误差将非常大。

#### d. 多重积分

​	用简单抽样的 Monte  Carlo 方法计算多重积分
$$
\int_{a_1}^{b_1}dx_1...\int_{a_n}^{b_n}dx_n~f(x_1,...x_n)\approx\frac{1}{N}\left[\prod_{j=1}^n(b_j-a_j)\right]\sum_{i=1}^Nf(x_{1_i},....x_{n_i})=V\left<f\right>
$$
其中对每个坐标的抽样是在响应区间内均匀抽取的。

​	**优势：**对于固定的样本数 N 值，Monte Carlo 方法给出的误差~ $1/\sqrt N$ ，而对于固定网格点法，由于每一维上的平均点数为 $N^{1/d}$ ，性能较好的方法如复化Simpson积分公式，给出的误差为 $1/N^{2/d}$ ，没有Monte Carlo方法号。

​	当多重积分的维数d > 4时，几乎没有其他数值计算方法可以超过 Monte Carlo 方法。

​	**劣势：**如果积分中的被积函数不光 滑时，如果函数式中的变化急剧时，则用目前的非权重简单抽样 Monte Carlo 方法 得到的精度可能很差。简单抽样中如果要提高精度，必须增大抽样 量。这种办法是不可取的，因为这样的计算效率很低。引入带权重的重要抽样法， 将既可以保证计算精度又有很高的效率

#### e. 提取法

​	由于$\sigma_s$ 正比于 $\sigma_f$ ，因此将被积函数变为平坦还可以有效地提高 Monte Carlo 方法的精度，因此发展了几种抽样技巧以减小方差，包括**提取法**和**重要抽样法。** 

​	对于一维积分，设我们能够构造一个与被积函数 f x( ) 形状相似的函数 g x( ) ，且它的积分值已知
$$
\left|f(x)-g(x)\right|<\epsilon
\int_a^bg(x)dx=J
$$
​	对于提取法，可以把 $f(x)$转换乘较为平缓的 $h(x)=f(x)-g(x)$,然后再进行简单抽样。
$$
\int_a^b f(x)dx=\int_a^bh(x)dx+J
$$

#### f. 奇异积分

​	奇异积分是指在积分限内被积函数有发散点的情形，尽管在奇点上被积函数是无穷大，但积分值可以是有限的，如 $\int^1_0dx/\sqrt{1-x^2}=\pi/2$ 中的被积函数在上限处发散的。

​	在数值计算方法中，奇点必须避开，但因为奇点附近的贡献是相当大 的，不可忽略，因此要尽可能的逼近奇点，这在固定网格划分形式下的其它数值 方法中是很难做到的。

​	但用 Monte Carlo 方法就比较容易，这是因为，随机选取 的 x 恰为奇点的可能性不大，还可以在包含奇点的附近选择一非常小的排斥区域以避开奇点。

​	不过，简单抽样的效率很低。

### 3.2 重要抽样

#### a. 重要抽样方法

​	**思想：**由于$\sigma_s$ 正比于 $\sigma_f$ ，因此将被积函数变为平坦还可以有效地提高 Monte Carlo 方法的精度，对于重要抽样方法，是利用除法来改造被积函数。

​	设有一个几率分布 p(x) 与 f(x) 形状相似，则有
$$
\frac{f(x)}{p(x)}=\Theta(1)\\
\int_a^bp(x)dx=1\\
\int_a^bf(x)dx=\int^b_a\frac{f(x)}{p(x)}p(x)dx
$$
<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221002195255721.png" alt="image-20221002195255721" style="zoom:50%;" />

​	如果 x 在[a , b] 区间内的随机抽样不是均匀选取的，而是按照几率分布 g(x) 选取的则**重要抽样法下的 Monte Carlo 积分**为
$$
\int f(x)dx =\left<f/p\right>\approx \frac{1}{N}\sum_i^N\frac{f(x_i)}{p(x_i)}
$$
​	对于[a , b]区间内的均匀分布，显然应该有 $p(x)=1/(b-a)$，代入上式中，则回到简单抽样公式

#### b. 带权重的Monte Carlo 积分

​	如果想要变换积分区间，从[a, b]到[0, 1]，则需要通过适当的变量代换。
$$
y = \int_a^xp(t)dt=F(x)\\
dy = p(x)dx
$$
显然变量代换函数为累积函数。于是积分公式变为：
$$
\int_a^bf(x)dx = \int_a^b\frac{f(x)}{p(x)}p(x)dx=\int_0^1\frac{f(y)}{p(y)}dy
$$
​	上式对 y 的简单抽样即为原重要抽样公式。

​	可见，**重要抽样法的精神是**：将随机变量 x 变换为另一随机变量 y ，而变换关系 y(x) 是表征被积函数 f(x) 曲线变化特征的 p(x) 的累积函数，对 y 的简单抽样就是对 x 带权重 p(x) 的抽样。

​	更一般地，如选择的  p(x) 本身不是归一化的 $c=\int_a^bp(x)dx$，则归一化步 骤相当于将式的上下限作一简单的替换：
$$
I= \int_0^c\frac{f(y)}{p(y)}dy
$$
​	其中对 y 的均匀抽样是在[0, c]区间。

​	**在实际中常常使用泰勒展开，或者常见的分布函数作为p(x)。**



## 4 布朗运动-----随机行走与生长问题

扩散过程是粒子在介质中随机行走，无规性来自运动的随机性

### 4.1 随机行走

#### a. Brown运动

系综：相同条件下的粒子系统的copy，用来模拟大量相同状态的粒子的平均

- 一维坐标系，从原点出发，以p概率向左走1步，1-p概率向右走1步，大量粒子平均后，有

==$\left<\Delta x\right>=0$==

==$\left<\Delta x^2\right>=2D\tau=N$==

即$\tau$时间内，**粒子偏离原点的距离平均** $\sqrt{\left<\Delta x^2\right>}$ 与 **$\sqrt \tau$**成正比，

- 或理解为与行走的步数 $N=2D\tau$的平方根 $\sqrt N$ 成正比（假设匀速运动）。

假设每一步走 $l$步长，向左向右走概率相等，每一步为 $X_i$，且**相互独立**，那么
$$
\left<X_i\right>=\frac12l+\frac{1}2(-l)=0\\
\left<X\right>=\sum \left<X_i\right>=0\\
\left<X_i^2\right>=\frac12l^2+\frac{1}2(-l)^2=l^2\\
\left<X^2\right>=\sum \left<X_i^2\right>=Nl^2\\
$$

- 或理解为与 二项分布（撒豆子模型），最后豆子的落点位置:$\sigma^2=Np(1-p)$

- 扩散方程

对于连续变量，粒子在时间t，出现在x处的概率有
$$
p(x,t)=\frac{1}{\sqrt{4\pi Dt}}exp{\frac{x^2}{4Dt}}
$$
这等价于扩散方程的解（初态为 $\delta(x)$ ）
$$
\frac{\part p}{\part t}=D\frac{\part^2p}{\part x^2}
$$
有限边界热传导方程：
$$
\left\{
\begin{aligned}
u_t(x,t)=Du_{xx}(x,t)\\
u(x,t=0)=f(x)
\end{aligned}
\right.
$$
则函数的解是
$$
u=f*p
$$

#### b. Einstein理论

扩散（涨落）与粘滞（耗散）之间存在联系
$$
D=\frac{kT}{6\pi\eta a}=\frac R{N_A}\frac T{6\pi\eta a }
$$

#### c 三维

$$
\left<r(t)\right>=0\\
\left<r^2(t)\right>=\int \bold{r}^2p(r,t)~d\bold{r}=\int r^4p(r,t)dr=\left<x^2\right>+\left<y^2\right>+\left<z^2\right>=6Dt\\
$$



#### d. 涨落与自相关系数

**方差 variance**：$var(x)=\left<(x-\left<x\right>)^2\right>=\left<x^2\right>-\left<x\right>^2$

**协方差 covariance**：$cov(x,y)=\left<(x-\left<x\right>)\right>\left<(y-\left<y\right>)\right>=\left<xy\right>-\left<x\right>\left<y\right>$

**相关系数**：$cor(x,y)=\frac{cov(x,y)}{var(x)var(y)}$

##### **涨落**

$\delta A(t)=A(t)-\left<A\right>$，where $\left <A\right>=\left<A(t)\right>=\left<A(0)\right>$

##### **自相关系数**

$$
\begin{aligned}
C(t)&=cov(A(t),A(0))=\left<\delta A(t)\delta A(0)\right>\\
&=\left<A(t)A(0)\right>-\left<A(t)\right>\left<A(0)\right>\\
\end{aligned}
$$

一般会默认为$\left<A\right>=0$，所以自相关系数有时也直接写为
$$
C(t)=\left<A(t)A(0)\right>
$$
在初态和时间无穷的稳态
$$
C(0)=\left<(\delta A)^2\right>=\left<A^2\right>-\left<A\right>^2
$$

$$
\lim_{t\rarr\infin}C(t)=\lim_{t\rarr\infin}\left<\delta A(t)\delta A(0)\right>=\left<\delta A(t)\right>\left<\delta A(0)\right>=0
$$

即在无穷处，可以认为二者几乎没有相关性

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221010213823685.png" alt="image-20221010213823685" style="zoom:67%;" />

##### **第一涨落耗散定理**

由 $dX=Vdt$，$\frac{\left<X\right>}{dt}=\left<V\right>$，$\frac{\left<X^2\right>}{dt}=2\left<XV\right>$

可推导出结论：$\frac{d}{dt}var(X)=\frac d{dt}(\left<X^2\right>-\left<X\right>^2)=2\left<XV\right>-2\left<X\right>\left<V\right>=2cov(X,V)$

那么，再根据 $r(t)-r(0)=\int^t_0v(t')dt'$，$\left<\bold r^2\right>=var(\bold r(t)-\bold r(0))=6Dt$

可以得到 （设$\bold r(0)=0$，且 $\left<r(t)-r(0)\right>=0$)
$$
\begin{aligned}
\frac d{dt}var(r(t)-r(0))&=2cov(\bold r(t),\bold v(t))\\
&=2\left<\bold r(t)·\bold v(t)\right>-2\left<\bold r(t)\right>\left<\bold v(t)\right>\\
&=2\int^t_0\left<\bold v(t')·\bold v(t)dt'\right>\\
&=2\int^t_0\left<\bold v(0)·\bold v(t-t')dt'\right>=2\int^t_0\left<\bold v(0)·\bold v(t^{"})dt^{"}\right>(积分自变量变换t''=t-t')
\end{aligned}
$$
而当时间趋于无穷时，$\frac d{dt}var(\bold r(t)-\bold r(0))=\frac d{dt}\left<\bold r^2(t)\right>=6D$



于是得到**三维下第一涨落耗散定理**
$$
D=\frac13\int^{\infin}_0\left<\bold v(t)·\bold v(0)\right>dt=\int^{\infin}_0C(t)dt\\
定义自相关函数：C(t)=\frac13 \left<\bold v(t)\cdot\bold v(0)\right>
$$
**物理含义**：宏观扩散系数是微观速度自相关函数的积分



#### e. Langevin理论与Brownian动力学

即粒子受到随机涨落力作用下的动力学方程

##### **涨落力**

$$
\left<\bold F(t)\right>=0\\
\left<\bold F(t)\cdot\bold F(0) \right>=D\delta(t)
$$

##### **Brownian动力学**

假设Brown 粒子处于流体环境中，受到

1. 粘滞阻力，$-v/B$（$B=\frac 1{6\pi\eta a}$  是迁移率）
2. 快速涨落力 F(t) （系统特征时间 $\tau$ >>碰撞时间）

则得到著名的Langevin方程

##### **Langevin方程（随机微分方程）**

$$
m\frac{d\bold v}{dt}=-\frac 1B\bold v+\bold F(t)\\
\frac{d\bold v}{dt}=-\frac 1\tau\bold v+\bold A(t)(等价表达)
$$

其中 $\tau = Bm$

那么自然的有
$$
\bold v(t)=\bold v(0)e^{-t/\tau}+e^{-t/\tau}\int^t_0e^{t'/\tau}\bold A(t')dt'
$$
则速度自相关系数为
$$
\begin{aligned}
C(t)&=\frac13\left<\bold v(t)\cdot \bold v(0)\right>\\
&=\frac13\left< \bold v(0)^2\right>e^{-t/\tau}+e^{-t/\tau}\int^t_0e^{t'/\tau}\left< \bold v(0)\cdot \bold A(t')\right>dt'\\
&=\frac13\left< \bold v(0)^2\right>e^{-t/\tau}（v(0)与A(t)无关）\\
&=C(0)e^{-t/\tau}

\end{aligned}
$$
那么可以得到D的准确值
$$
D=\int^{\infin}_0C(t)dt=\frac13\int^{\infin}_0\left< \bold v(0)^2\right>e^{-t/\tau}dt=\frac13\left<\bold v^2\right>{\tau}=\frac 13 \frac{3kT}{m}\tau=kTB=\frac{kT}{6\pi\eta a}
$$
即 **Einstein关系**。

PS：粒子的速度分布符合Maxwell 速度分布 $p(v)=\sqrt{\frac{m}{2\pi kT}}\exp{-\frac{mv^2}{2kT}}$



#### f. **第二涨落定理**：

阻力摩擦力等与涨落有关。摩擦越小，涨落力越大

对Langevin方程作点积
$$
\begin{aligned}
\bold r\cdot\frac{d\bold v}{dt}&=\frac12\frac{d^2}{dt^2}r^2-v^2\\
&=-\frac 1\tau\bold r\cdot\frac{d\bold r}{dt}+\bold r \cdot\bold A(t)
\end{aligned}
$$
取平均，易知 $\left<\bold r\cdot \bold A\right>=0$，所以
$$
\frac{d^2}{dt^2}\left<r^2\right>+\frac1\tau\frac d{dt}\left<r^2\right>=2\left<v^2\right>=\frac{6kT}m
$$
解得微分方程
$$
\left<r^2\right>&=\frac{6kT}{m}\tau^2[\frac t\tau-1+\exp(-t/\tau)]\\
&\approx\left\{
\begin{aligned}
&\frac{3kT}mt^2=\left<v^2\right>t^2\quad(t\ll\tau)\\
&\frac{6Kt}m\tau t=6Dt\quad(t\gg \tau)

\end{aligned}
\right.
$$
以上方程说明了

在特征时间内，粒子是直线行走的；较长尺度下是随机行走的。



### 4.2 自规避随机行走模型

模拟高分子物理实体，即轨迹不能相交（贪吃蛇）

#### a. SAW模型

​	Brown 运动的粒子在作随机行走时，粒子先前在空间中漫游时任何位置的记忆对当前的运动没有任何限制，可以与自己的历史路径相交；而一条柔性高分子链是一个占据有限空间的物理实体，它不可能自相交。

​	在粒子的自 规避行走 SAW（self-avoiding walks）模型中，粒子要记住以往走过的格点位置， 在禁止返回上一位置的同时，当它与历史位置相交时粒子将死亡从而终止行走。

​	

#### b.指标计算

​	**指数律式**
$$
\left<r^2(N)\right>=aN^{2\nu}(1+bN^{-\Delta})
$$
​	由于自规避的排斥效应，尽管SAW行走结果造成的路径范围要较理想的RW 行走要大，但是大 N 的路径数目较少，因此在由**指数律式**求指数 $\nu$ 时， 要通过双对数曲线 $\log{\left<r^2(N)\right>-\log N}$ 求斜率。通常是由比值法求斜率，忽略小量，那么 $\nu$的计算方法为：
$$
\nu(N)=\frac12 \frac{\ln{[<r^2(N+1)>/<r^2(N-i)>]}}{\ln{[(N+i)/(N-i)]}}
$$
式子中 $i\ll N$，选取i 时要使它大到可以忽略附近的涨落，但又要比 N 小很多以 使修正项对指数的计算影响很小

#### c. 径向分布函数

​	径向分布函数的物理直观意义是：以一个原子为原点，它的径向上其它周围原子的分布概率

​	对于随机行走问题，我们最关心的是行走前后的净距离或首末端距长度r ， 将它的方均根值记为 $<r^2>=r_{rms}^2$。已知对于 **RW 模型**，标度指数律 $r_{rms}=lN^{1/2}$ 对于空间维数 d =1, 2,3均成立。在无偏压情况下，我们可将矢量r 认为是各向同性的。（1.4.1.4-6）式可改写成三维空间中的 Gauss 型概率密度分布函数

​	d维空间中随机行走 $RW$的径向分布函数是
$$
RDF_d(r)=C_d·r^{d-1}\exp[-r^2d/(2l^2N)]
$$


#### d. Flory-Fisher平均场理论

​	
$$
Z(R)=RDF_d(R)<\exp(-\beta U)>
$$
后面一项代表了玻尔兹曼分布分配的比例

### 4.3 非平衡生长

#### a. 生长模型的概念

​	生长是指一个对象随时间长大的过程

​	由于生长过程非常复杂，任何模型都不可能全面 地概括所有实际生长过程中的因素，一个模型中只要抓住最基本的因素即可，这正是我们在上面对高分子构型所做的事情。.

#### b.**Laplace生长**

$$
\frac{\part f }{\part t}=v\frac{\part f}{\part x}\\
\frac{1}{l_D}\frac{\part f }{\part x}=\frac{\part^2 f}{\part x^2}\\
扩散长度：l_D=D/v
$$

#### c. 扩散受限聚集模型（DLA）

​	它是在1981年由 T.A. Witten 和 L.M. Sander 提出的，最初是用来研究悬浮于大气 中的灰尘或粉末的扩散聚集过程。DLA 模型是用一个简单算法产生复杂无序图 形的最典型例子

​	格点 DLA 的模拟规则是，取一个2维的方形点阵，在点阵中 央原点处放置一个粒子作为生长的种子，然后从距原点足够远的圆周界处释放一 个粒子，让它作 Brown 运动或随机行走，其结果 是：该粒子走到种子的最近邻位置与种子相碰， 这时让粒子粘结到种子上不再运动；或者粒子走 到大于起始圆的更远处（如2-3倍的半径处）或干 脆走到点阵边界，这时认为粒子走了一条无用的 轨迹，取消该粒子，把它重新放回原点。因此， 那些有用的粒子与种子相粘结后形成不断生长的 聚集集团。

<center>	<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221017105502121.png" alt="image-20221017105502121" style="zoom:30%;" />
**算法：**

1. 随机生成生成一个粒子 $x$

   - 若 $R<N/2$，则在半径为 $R$的圆上，随机生成一个粒子

   - 否则，在方形点阵的边界上随机生成一个粒子
2. 让粒子 $x$随机行走，判断周围8个点是否有粒子集团
   - 若碰到点阵边界，则返回1，重新生成一个粒子

   - 若附近有集团粒子，依据x黏度 **`stickness`**判断是加入集团还是继续随机行走
3. 当粒子加入集团，更新集团半径 $r$。
4. 返回1；或当集团足够大时，程序结束



#### d. Levy随机行走或飞行

**Brown运动步长**限制，其概率由下式得到，$p(x)$未归一
$$
p(x)=\left\{
\begin{aligned}
s^{-\alpha},&\quad if~s\ge1\\
1,&\quad if ~s<1
\end{aligned}
\right.
$$

### 4.4 粒子输运

​	Monte Carlo 方法应用的一个非常重要的研究领域是物质中粒子输运过程以及粒子碰撞过程的模拟。

​	以材料中产生的中子为例，它从诞生至死亡的整个过程可以用一个随机过程描述。起始时，粒子以一定几率随机产生，其初始位置（源有一定大小）、能量 （非单色源情形下）、发射方向**均是随机变量**。在其后的历史运动过程中，其随机性一直保持，如行走一步的距离、紧跟的相互作用类型、散射后的方向和能量等等。根据散射机制，还可能有粒子的湮没和二次粒子的产生。物理上，我们只能够确定这些变化的几率或散射/反应截面，而不能对每个粒子确定它在何时何地确切地将发生什么事情，这就是我们为什么要用 Monte Carlo 方法研究粒子输 运过程的根本原因。

​	Monte Carlo数值模拟可以处理采用Boltzmann方程等解析方法难以计算的问题，可以把各种详尽的相互作用机制考虑在内，还可以方便处理材料的结合边界

**步骤：**

- 粒子的相互作用作用机制
- 数值模拟
- 要注意物理量量纲

#### a. 电子

**弹性散射**

**相对论修正的弹性散射**

**非弹性散射的离散时间方法**

**Monte Carlo 模型和步骤**

​	主要思想是：找到粒子随机行走的步长和散射角度，步长一般服从指数分布或是平均自由程，散射角度由微分散射界面来确定

## 5 逾渗问题

逾渗描述流体流体在无序介质中作随机的扩展和流动，

### 5.1 逾渗模型

#### a. 概念

​	逾渗（percolation）这个词汇是由数学家 J.M. Hammersley 在1957年创造的， 其目的是为了描述流体在无序介质中作随机的扩展和流动，无规性源自介质本身的无序结构。逾渗过程是一种极好的教学模型，可借以阐明相变和临界现象的一些最重要的物理概念

**逾渗分类**

- 座逾渗（site）：每个座是无规连接性的：每一个座或 者是被占据的（连接的、畅通的），或 者是不被占据的（不连接的、堵塞的）， 相应的概率分别为 p 和1− p
- 键逾渗（bond）：每条键或者是连接的，或者是不连接的。连接的概率为 p ，不连接的概率为 1− p 

**连接**：左为连接，右为非连接

​	<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024100640097.png" alt="image-20221024100640097" style="zoom:50%;" />

**集团**：所谓集团是指座与座（或键与键）之间互相连接而没有间断

​	对座逾渗，若两个 已占座可以通过由一系列最近邻的已占座连成的路径连接起来，则这两个 座属于同一集团；

​	对键逾渗，若两条键可以通过至少一条由联键连成的路径连接起来，则这两条键属于同一集团

#### b. 逾渗阈值

​	假定有一个大到可以忽略其边界效应的二维方形点阵，点阵上的座以概率 p 被随机地占有， 若相邻的座都被占据时，这些占座就可以连接成为一个集团。

​	当 p 增加时属于s ≥ 2 的集团的占座的比例也将增加，集团的大小也相应地增大，但仍然是有限的。当 p 增加到0.5时，集团进一步长大，有些集团连在一起形成相当大的集团。尽管与小浓度 时的情况相差甚大，但此时的关键性质并未改变，即所有集团的大小都是有限的。

​	再进 一步增加占据的浓度（下图），我 们可以看到一个很大的集团，它扩张到整个系统，从顶到底，从左到右。对于有限大小 的系统，这个扩张的集团称为**跨越集团**，跨越集团随点阵系统尺度的增长而增长，直到 无穷大，这个无限扩张的或无界的集团称为 **逾渗集团或逾渗通路**。注意，虽然逾渗集团 是无限大的，即s → ∞，但它并非占据全部 点阵（只有当 p =1的高密度极限时棋盘的所 有座才被全部占据），它与一些有限大小的 集团以及空座所形成的岛同时并存。

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024102240581.png" alt="image-20221024102240581" style="zoom:67%;" />

​	计算结果显示，在 p = 0.59时逾渗 通路开始出现，0.59就是在正方形点阵上座 逾渗的临界浓度 $p_c$ 。在 $p_c$ 以上，存在逾渗 通路，在此以下就不存在逾渗通路。从 $p_c$  到 p =1,逾渗通路不断丰满，最后占满整个点 阵。

​	**逾渗阈值 **$p_c$ ： 关联长度达到无穷大，几何上理解为上下导通了

​	对于相同的点阵（图1.6.1.3-3），座逾渗的 $p_c$ 值都要比键逾渗的值大，这表明键 逾渗要比座逾渗容易发生，实际上，逾渗阈值由配位数决定，配位数增加，$p_c$减小

<center>	
<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024102400992.png" alt="image-20221024102400992" style="zoom:60%;" />

### 5.2 逾渗与相变

#### a. 物理量描述

s：大小为s的集团

$N_s$：大小为s的集团数目

N：点阵格点总数

平均集团大小的分布：$n_s(p)=N_s/N$

某个点属于s集团的概率: $P_s(p)=s*ns_(p)$

点阵中所有被占据格点数：$N\sum_ssn_s$

随机选择的某一占据格点属于大小为s 的集团（排除s = ∞ 的无限大集团的概率 $w_s=sn_s/\sum_ssn_s$

集团平均大小：$S=\sum_ssw_s=\sum_ss^2n_s/\sum_ssn_s$

平均跨越长度 $\xi$

s 个占据格点的单一集团的回转半径 $R_s$

#### b. 集团的标识——判断导通方式

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024103141426.png" alt="image-20221024103141426" style="zoom:80%;" />

**导通**：判断上下是否由相同标号的集团

#### c. **一维逾渗**



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024105400602.png" alt="image-20221024105400602" style="zoom:67%;" />



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024105436952.png" alt="image-20221024105436952" style="zoom:67%;" />



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024110154505.png" alt="image-20221024110154505" style="zoom:67%;" />

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024110345429.png" alt="image-20221024110345429" style="zoom: 33%;" />





<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024105041570.png" alt="image-20221024105041570" style="zoom:67%;" />



<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024110925317.png" alt="image-20221024110925317" style="zoom:67%;" />



#### d. 有限尺度标度法

#### e. 标度律

#### f. 边缘维数

### 5.3  数值重整化

#### a. 重整化群概念

​	Kadanoff 在1966年按照粗粒平均的思想引进了重整化群理论中的基本物理概 念使之运用于提取临界区域的临界指数，而无需研究配分函数。

​	重整化群方法的一般思路是将描写一个物理问题的参数用另外一组更 为简单的参数表示出来，但保持问题中感兴趣的物理性质不变。

​	重整化群方法有两类，一种是与量子场论中的方法类似的k 空间方法，另一 种是简单明了的实空间重整化群方法。

​	为了引入重整化群的概念，我们先考虑不同 p 值下生成的逾渗集团图形在尺度变换时的行为。

​	对于 $p<p_c$ 形成的逾渗图案，我们认为占据格子是黑格子，未占据格子是白格子。则当我们的肉眼逐步远离该图像时，由单独的 黑格子形成的小集团变得不可见了，只有大集团才可以看得见，并且这些大集团 中间连接各小集团的细窄跨桥或集团边缘上的细小分枝也是不可见的，即大集团 在这个尺度下变成分立的小集团。再进一步远离图像，则这些小集团也会逐步不 可识别。这个图像在尺度变换下最终趋于不动点 p = 0的图像。

​	对于$p>p_c$  的情 形结果应该是怎样呢？这时，图案中黑格子较多，居于统治地位，按类似相同的 理由，当逐步远离图像时，由白格子形成的小集团变得不可见了，也就是说图像 逐步趋于黑格子形成的大集团，最终趋于不动点 p =1的图像。

​	那么，在黑白格 子比重均衡时的逾渗阈值$p=p_c$ 下，尺度变换的结果又是什么？这时，图像的性 质不随尺度的变化而变化，与距离无关，可以近似地认为是缩小的原图像，也即 趋于另外一个不稳定的不动点。

​	![image-20221024112633017](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024112633017.png)

​	它的基本思想就是对体系的长度尺 度连续不断地做变换，将体系元胞尺度由 a 变换成$b*a$ （$b*a$ 应小于体系的相关长度 ξ ），相继标度变换的结果产生出一个流向图，空间流向场趋向于若干特殊的不 动点，这些点在标度变换下保持不变



#### b. 标度变换法

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024113341034.png" alt="image-20221024113341034" style="zoom:80%;" />

#### c. Monte Carlo重整化群方法

![image-20221024113405605](C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221024113405605.png)

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221031095835085.png" alt="image-20221031095835085" style="zoom:50%;" />
