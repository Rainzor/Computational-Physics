# Report8

> PB20020480 王润泽

### 1 Question

用 Monte Carlo 方法计算如下定积分，并讨论有效数字位数。
$$
\begin{aligned}
I_1 =\int f_1(x)dx =\int_0^5dx\sqrt{x^2+2\sqrt x}\\
I_2=\int f_2(x,y,z,u,v)dxdydzdudv=\int_0^{7/10}dx\int_0^{4/7}dy\int_0^{9/10}dz\int_0^2du\int_0^{13/11}dv(5+x^2-y^2+3xy-z^2+u^3-v^3)
\end{aligned}
$$



### 2 Algorithm

#### 2.1 积分 $I_1$——重要抽样法

​	采取 **重要抽样方法**，选择的 $p(x)=(x-2.5)+f_1(2.5)=0.567944 + x$

​	被积函数 $f_1(x)与p(x)$ 如下图所示

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\f1(x).png" style="zoom:67%;" />

<center><p>图1：f(x)与g(x)图像

​	这样将积分转化为： $\int_0^5f(x)dx=\int^5_0\frac{f(x)}{g(x)}g(x)dx$

​	为了方便积分，对g(x)进行归一化处理
$$
p(x)=\frac{g(x)}{\int_0^5g(x)dx}=\frac{g(x)}{15.3397}
$$
那么算法如下:

1. 随机生成N个以p(x)分布的采样点(舍选法)
2. 计算对应的 $\{y_i\}$，$y_i=f(x_i)/p(x_i)$
3. 求得数值积分表达式 $I(N)=\frac{\sum y_i}{N}=\left<y\right>$
4. 对 $N\in\{2^k\}_{k=5}^{20}$的不同取值，计算积分误差

$$
e_k= \left|I(N_k)-I_2\right|\approx|I(N_k)-15.4390107|
$$


$$
e_k= \left|I(N_k)-I_1\right|\approx5\frac{\sqrt{\left<y^2\right>-\left<y\right>^2}}{\sqrt {N_k}}
$$

**注：**大数定理适用于N较大情况，所以对于N较小的情况(N<50)没有考虑

#### 2.2 积分 $I_2$——简单抽样法

​	**采取平均值法积分**，算法如下：

1. 分别在 [0,7/10]上的均匀分布的抽样值$\{x_i \}$，[0,4/7]上均匀分布的抽样值 $\{y_i\}$ ， [0,9/10]上均匀分布的抽样值 $\{z_i\}$ ，[0,2]上均匀分布的抽样值 $\{u_i\}$ ，[0,13/11]上均匀分布的抽样值 $\{v_i\}$生成 $N$个随机值
2. 计算 $f_i =5+x^2-y^2+3xy-z^2+u^3-v^3=5+(x+\frac{3}{2}y)^2-\frac{13}{4}y^2-z^2+u^3-v^3$
3. 求数值积分 $I(N) = \frac{7}{10}*\frac{4}{7}*\frac{9}{10}*2*\frac{13}{11}\frac{\sum f_i}{N}=\frac{234}{275}\left<f\right> $
4. 对 $N\in\{2^k\}_{k=5}^{20}$的不同取值，计算积分误差

$$
e_k= \left|I(N_k)-I_2\right|\approx|I(N_k)-5.67712092|
$$


$$
e_k= \left|I(N_k)-I_2\right|\approx\frac{234}{275}\frac{\sqrt{\left<f^2\right>-\left<f\right>^2}}{\sqrt {N_k}}
$$

### 3. Experiment

#### 3.1 积分 $I_1$结果

​	对于一维积分，实验中采取了带权抽样方法，使得图像更为平缓

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221003133855037.png" alt="image-20221003133855037" style="zoom:67%;" />

<center><p>图2：f1(x)/g(x)图像

​	实验结果如下图所示，可以看到**在N>5000后，有效位有4位**

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\result1.png" style="zoom: 80%;" />

<center><p>图3

​	为了验证积分结果的有效位数精度的变化趋势，有下图

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\error1.png" style="zoom: 80%;" />

<center><p>图4：Monte Carlo方法的误差趋势

可见确实符合 $O(1/\sqrt N)$的趋势，从另一方面佐证了实验的可靠性。

​	但同样在实验中会发现，在N较小时，误差会有较大波动性。

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\BadError1.png" style="zoom:80%;" />

<center><p>图5：不稳定的结果

​	可能是由于一开始误差就降到0.05以内，导致后续存在了一些波动。在N足够大（N>10000）时，仍然会趋于稳定。

​	这也间接的说明了Monte Carlo方法可以用较少的点，就可以快速逼近精确值。

### 3.2 积分 $I_2$结果

​	实验结果如下，在**N>20000后，至少有2位的有效位数**

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\result2.png" style="zoom:67%;" />

<center><p>图6

​	误差变化趋势如下

<img src="F:\MyDocuments\Physics\Compututation Physics\Homework\hw08\error2.png" style="zoom:80%;" />

<center><p>图7

​	积分误差符合 $O(1/\sqrt N)$的趋势，在数值N<1000时，会存在一些波动，但当N较大时，仍然会趋于稳定。

### 4.Summary

​	本次实验中利用 Monte Carlo 方法获得了两个定积分的近似解，估计了计算结果的有效位数。

​	最后验证了积分误差和抽样点的个数满足正比于 $1/\sqrt N$ 的 关系，并验证了在N较大时，Monte Carlo方法具有精确性、稳定性。
