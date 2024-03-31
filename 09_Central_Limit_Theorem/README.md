# README

## Question

考虑泊松分布、指数分布，并再自设若干个随机分布（它们有相同或不同的 $\mu$ 和 $\sigma^2$ ，通过 Monte Carlo 模拟，验证中心极限定理成立 $(N=2、5、10)$ 

## code

主程序在 "hw9.py" 中

## data

data中提供了测试集数据“data_09_2.csv”,“data_09_5.csv”,“data_09_10.csv”,“data_09_10.csv”存储不同N得到的标准化统计量
每个csv中含有不同分布所得的统计量.

代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数，也可以直接调用程序即可画图

## Explanation

主要用到numpy、time、pandas、matplotlib、abc、math包


本实验对16807Schrage产生器进行单独封装，其代码在"Schrage_16807.py"中。内部包含"seed_time"随机数产生种子函数，"Schage16807"类

同时对各类分布随机数产生以及相应的概率密度函数进行的封装，其代码在"myStats.py"中

内部包含"norm"正态分布，"expon"指数分布，"poisson"泊松分布，"uniform"均匀分布，"bernulli"伯努利分布，"binom"二项分布，"Cos"余弦分布

建构逻辑仿照了标准Scipy库中"`stats`"的方式:

- `rvs` 是产生随机数的方法

- `ppd/pmf` 为概率密度函数

- `cdf` 为累积分布函数

- `pmf` 为分位点函数

- `mean` 为分布平均值，var为方差，std为标准差，mid为中位数