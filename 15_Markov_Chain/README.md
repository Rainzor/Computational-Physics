# README

## Question

设体系的能量为 

$$
H(x,y) = -2(x^2+y^2)+\frac12(x^4+y^4)+\frac12(x-y)^4
$$

取 $\beta=1/kT=0.2, 1, 5$

采用 Metropolis 抽样法计算 $\left< x^2 \right> ,\left< y^2 \right> , \left< x^2 +y^2\right>$ 。抽样时在 2 维平面上依次标出 Markov 链点分布，从而形象地理解 Markov 链。

## Code

主程序在“hw15.py”中
程序中含有”MetroSample“类，使用Metroplis-Hastings 算法抽样特定的分布。

类中包含有各种方法，详细见注释

"hw15.py" 中写好了 `plot_mcmc()`,`plot_born()`函数，用来展现 Markov链的过程，您可以取消注释查看图像

code中还提供了matlab文件，用来检验数据和生成图像

## Data

data中提供了实验数据

data_beta_0.2.csv,data_beta_1.csv,data_beta_1.csv 分别是beta取不同值时，抽样得到10000个点分布，存储其（x,y）坐标值

burn_beta_1.csv 是在 beta = 1下利用简单Metropolis方法抽样得到的Markov Chain轨迹坐标值

## Explanation

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

内部包含“seed_time”随机数产生种子函数，“Schage16807”类，故代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数，

主要用到numpy、time、pandas、matplotlib包



