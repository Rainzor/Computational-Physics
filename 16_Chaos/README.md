# README

## Question

以 $x_{n+1}=\lambda\sin(\pi x_n)$ 为迭代方程进行迭代

(1). 画出系统状态随参数入的变化图，要求在图中体现出定值状态、倍周期分叉和混沌状态;

(2). 列出各个倍周期分叉处的 $\lambda$ 值，求相应的 Feigenbaum 常数。

## Code

主程序在“hw16.py”中，主要用到'numpy'、'pandas'、'matplotlib'包

主程序会输出一张lambda在[0,1]之间变化的图像
"hw16.py" 中有保存数据的操作，已注释掉了

## Data

data中提供了实验数据

data.csv 是最原始的lambda与相应状态x的取值，方便画散点图

Feigenbaum_size.csv 保存的是不同lambda与对应的周期数

Feigenbaum_alpha.csv 保存了在x=0.5附近时，lambda取值

