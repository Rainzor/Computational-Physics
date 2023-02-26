# Report9

> Rainzor

## 1.Question

​	考虑泊松分布、指数分布，并再自设若干个随机分布（它们有相同或不同的 μ 和 𝜎 2 )，通过 Monte Carlo 模拟，验证中心极限定理成立（N=2、5、10)。

## 2. Algorithm

​	由中心极限定理：
$$
X_N=\frac{\sum^N_i\frac{x_{i}}N-\mu}{\sigma/\sqrt N}\sim N(0,1)
$$
​	其中 n是统计量 $x_N$的个数，$\mu$是分布的期望，$\sigma$是分布的标准差，N(0, 1)是标准正态分布

​	要验证中心极限定理，就要得到各个分布的随机数 $r$，然后以 $x_N$作为统计量,进行标准化处理，验证中心极限定理。
$$
x_N= \frac{\sum_i^Nr_i}N
$$
​	**算法如下**：

1. 生成满足概率分布为 $f$的抽样样本 $r$，得到分布对应的期望与方差
2. 对N（2，5，10，50）个样本计算其统计量 $x_N$，标准化处理后，将结果 $X_n$ 存储到文件中。
3. 重复统计M次，绘制相关图像，与标准正态分布比较



​	实验用到的**分布**有

- 指数分布：$f(x)=e^{-x}$

- 泊松分布：$P_n(\lambda=5 )=e^{-\lambda}\frac{\lambda^n}{n!} $

- 二项分布：$B(k,n=10,p=\frac12)=C_n^kp^k(1-p)^{n-k}$

- 均匀分布：$U[0,1]$

- 余弦分布：$Cos(x),x\in[0,\pi/2]$

  

### 3 Experiment

#### 3.1 指数分布

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw09_大数定理验证\Exponential.png)
<center><p> 图1：指数分布验证中心极限定理

#### 3.2 泊松分布

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw09_大数定理验证\Poisson.png)

<center><p> 图2：泊松分布验证中心极限定理

#### 3.3 二项分布

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw09_大数定理验证\Binomial.png)

<center><p> 图3：二项分布验证中心极限定理

#### 3.4 均匀分布

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw09_大数定理验证\Uniform.png)

<center><p> 图4：均匀分布验证中心极限定理

#### 3.5 余弦分布

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw09_大数定理验证\Cos.png)

<center><p> 图5：余弦分布验证中心极限定理

### 4.Summary

​	实验中可以看到,5类分布最终当N取得越大时,越接近正态分布,从而论证了中心极限定理的正确性.
