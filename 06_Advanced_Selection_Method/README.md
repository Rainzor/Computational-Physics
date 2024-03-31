# README

## Question

对两个函数线型（Gauss 分布和 类Lorentz 型分布），设其一为 p(x)，另一为F(x)，其中常数 $a\neq b\neq 1$  , 用舍选法对 p(x) 抽样。将计算得到的归一化频数分布直方图与理论曲线 p(x) 进行比较，讨论差异，讨论抽样效率。

$$
\begin{aligned}
    &\text{Gaussian}\sim \exp(-ax^2)\\
    &\text{Lorentzian}\sim \frac{1}{1+bx^4}
\end{aligned}
$$

## Code

本实验对16807Schrage产生器进行单独封装，其代码在“code/Schrage_16807.py”中。

主程序在“code/hw6.py”中。

## Data

提供了测试集数据“data_06.csv”，存储着标准正态分布的采样点。
代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数。

## Import

主要用到numpy、time、pandas、matplotlib包