# README

## Question

对于球面上均匀分布的随机坐标点，给出它们在（x, y）平面上投影的几率分布函数。

并由此验证Marsaglia抽样方法:

$$
\begin{aligned}
x=2u\sqrt{1-r^2}\\
y=2v\sqrt{1-r^2}\\
z=1-2r^2
\end{aligned}
$$

确为球面上均匀分布的随机抽样

## Code

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

主程序在“hw5.py”中。

## Data

提供了测试集数据“data_05.csv”，存储着球面坐标

代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数。

## Import

主要用到numpy、time、pandas、matplotlib包