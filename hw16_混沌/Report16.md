# Report16

> 王润泽 PB20020480

## 1 Question

​	以 $x_{n+1}=\lambda\sin(\pi x_n)$ 为迭代方程进行迭代

(1). 画出系统状态随参数入的变化图，要求在图中体现出定值状态、倍周期分叉和混沌状态;

(2). 列出各个倍周期分叉处的 $\lambda$ 值，求相应的 Feigenbaum 常数。

## 2 Analysis

​	主要使用迭代法进行计算，对于 $\lambda=0.8$的系统，大致呈现如图1所示的类型，与一维 $Logistic$方法类似

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw16\compare.png)
<center> 图1：迭代方程形式
对于有固定周期的态，结果与初始值无关，且在有限步迭代后必然收敛。然而对于混沌的体系，周期无穷，那么最终绘图的结果与初值也无关，不妨取初值为 $x_0=0.5$

​	设$\lambda\in(0,1)$，步长为 $\delta\lambda=10^{-3}$，不断增大 $\lambda$ 的值进行如下迭代：

初始值 $x_0=0.5$,$k=1$，$\lambda_0=0$，$M=5000$，$N=1000$

1. 第k个$\lambda$取值为： $\lambda_k=\delta\lambda+\lambda_{k-1}$

2. 根据公式 $x_{n+1}=\lambda\sin(\pi x_n)$对 $x_n$进行M步迭代，再舍去，热化掉不稳地的点

3. 创建变量 `data=[]` ，再进行N步迭代，将新产生的点 $x$ 加入到 `data`中。

   但考虑到由于计算机有精度限制，所以在迭代时会有 $|x_i-x_j|$极小的情况，故为了避免这种情况要对 $|x_i-x_j|<10^{-5}$的点舍去 

4. 将当前的 $\lambda_k$与系统状态存储在 `lambs`和`data_list`中

## 3 Experiment

### 3.1 状态变换图

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw16\Bifurcation_diagram.png)
<center> 图2：倍周期分岔到混沌示意图

​	图 2 是将 $\lambda$ 从 0 变化到 1 所生成出来的一个图像。当 $\lambda$ <0.3 时，迭代序列稳定地收敛到 0。当 $\lambda$ 在 0.3 附近时，短暂地出现了一次倍周期分岔。当0.33 < $\lambda$ < 0.7时，迭代序列收敛到单一值。当 $\lambda$ > 0.7时，迭代序列分别展现出倍周期分岔和混沌的形态。

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw16\Bifurcation_diagram_expend.png)
<center> 图3：局部放大图

​	将图 2 在0.7 <  $\lambda$< 1的位置放大，得到图 3，从图中可以清晰地看到有 T=2、 T=4、T=8、T=16 的分岔点。

### 3.2 分叉处的 $\lambda$ 值与 Feigenbaum 常数

#### 3.2.1 横轴方向倍周期分岔中的标度行为

​	$\lambda_m$ 按以下的几何级数（幂函数）收敛到 $\lambda_{\infin}$
$$
\lambda_{\infin}-\lambda_m=A\delta^{-m}(m\gg1)
$$
​	A 是依赖于迭代函数的常数，而 $\delta$ 是不依赖于 迭代函数的普适常数	

​	为了找到分叉点，可以调整 $\lambda$的区间和精度，放大局部图像，使得可求出更加准确的分叉处$\lambda$值。

​	此时取 $\lambda\in[7.15,7.25]\cup[0.83,0.835]\cup[0.85,0.87],\delta\lambda=10^{-6}$，程序可输出对应周期数与 $\lambda$值 ，数据存储在 `Feigenbaum_size.csv`文件中，在分叉处，或许是因为计算精度有限，周期会出现波动，故要在周期稳定时取值较为合适，观察表格中的数据与公式得到下表 
$$
\text{Feigenbaum }\delta= \frac{\lambda_m-\lambda_{m-1}}{\lambda_{m+1}-\lambda_m}=\frac{\Delta\lambda_{m}}{\Delta\lambda_{m+1}}
$$
| m    | T           | $\lambda_m$ | $\Delta\lambda_m$ | $\delta$ |
| ---- | ----------- | ----------- | ----------------- | -------- |
| 1    | $1\rarr2$   | 0.720043    |                   |          |
| 2    | $2\rarr4$   | 0.833299    | 0.113256          |          |
| 3    | $4\rarr8$   | 0.858621    | 0.025322          | 4.472632 |
| 4    | $8\rarr16$  | 0.864089    | 0.005468          | 4.630944 |
| 5    | $16\rarr32$ | 0.865261    | 0.001172          | 4.665529 |
| 6    | $32\rarr64$ | 0.865514    | 0.000253          | 4.632411 |
<center>  表 1：分岔横向距离和 Feigenbaum 常数计算

​	即前后分岔间距的比值趋向于一个常数：$4.669 201 609$

#### 3.2.2 纵轴方向倍周期分岔的标度行为

如下图所示，取 $x=0.5$,其纵向标度的比值趋于一个常数，即
$$
\frac{d_m}{d_{m+1}}\rarr\alpha=2.502 907 875 095 892 822 283 902 87\quad(m\rarr\infin)
$$
<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw16\mid.png" style="zoom: 25%;" />
<center> 图4：纵向距离交点示意图

​	对于此题如果取在 $x=0.5$处的值，计算分岔的纵向间距比值，可以看到Feigenbaum $\alpha$是趋于理论值的。

​	而计算机精度有限，故实验中选择了 $0.4995<x<0.5005$之间的值输出。同时要注意相同周期数，不同分叉之间纵向距离不一定相同，所以要取在 $x=0.5$分叉处的纵向距离才能得到正确的结果。

​	数据保存在 `Feigenbaum_alpha.csv`中，如下表所示




| m    | T    | $d_m$    | $\alpha$ |
| ---- | ---- | -------- | -------- |
| 1    | 2    | 0.277659 |          |
| 2    | 4    | 0.047581 | 5.83550  |
| 3    | 8    | 0.019911 | 2.52277  |
| 4    | 16   | 0.007477 | 2.66296  |
| 5    | 32   | 0.003082 | 2.42602  |
<center>  表 2：分岔纵向距离和 Feigenbaum 常数计算

​	结果已经比较接近理论值：$2.502 907 875$

### 3.3 更大范围内变化

如果将区间范围调整到 $\lambda\in[-2,2]$范围，可以看到迭代序列出现了负值，这和迭代函数的函数取值是吻合的，同时呈现处一定的对称性

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw16\Bifurcation_diagram_big_scale.png" style="zoom: 33%;" />

<center> 图5：更大范围内状态的变化

## 4 Summary

​	利用迭代序列 $x_{n+1}=\lambda\sin(\pi x_n)$ 进行迭代，观察到了迭代序列随 $\lambda$ 会从定值变化到周期再变化到混沌。而且验证了 Feigenbaum 常数是一个和迭代序列无关的普适常数。
