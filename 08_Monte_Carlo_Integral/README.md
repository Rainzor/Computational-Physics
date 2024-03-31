# README

## Question

用 Monte Carlo 方法计算如下定积分，并讨论有效数字位数。

$$
\begin{aligned}
I_1 &=\int f_1(x)dx =\int_0^5dx\sqrt{x^2+2\sqrt x}\\
    I_2&=\int f_2(x,y,z,u,v)dxdydzdudv\\
    &=\int_0^{7/10}dx\int_0^{4/7}dy\int_0^{9/10}dz\int_0^2du\int_0^{13/11}dv(5+x^2-y^2+3xy-z^2+u^3-v^3)
\end{aligned}
$$

## Code

本实验对16807Schrage产生器进行单独封装，其代码在"code/Schrage_16807.py" 中。

主程序在"code/hw8_1.py"和"code/hw8_2.py"中，进行两种类型的积分计算与结果。

程序输出积分结果与误差，并绘制误差趋势图像。

## Data

压缩包中提供了测试集数据"data_08_1.csv"，"data_08_2.csv"，存储：样本个数N,积分值Integral，误差Error

代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数。

## Explanation

代码中定义了 `simple_monte_carlo_int`、`weight_monte_carlo_int`和`monte_carlo_iint`函数，用来计算数值积分

主要用到numpy、time、pandas、matplotlib、scipy包