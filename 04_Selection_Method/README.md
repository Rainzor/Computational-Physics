# README

## Question

设pdf函数满足关系式：

$$
p'(x)=a\delta(x)+b\exp(-cx),\quad x\in[-1,1],a\neq0
$$

讨论该函数性质并给出抽样方法。

## Code

本实验对16807Schrage产生器进行单独封装，为的是以后实验方便继续调用，其代码在“Schrage_16807.py”中。

主程序在“hw4.py”中。

## Data

提供了测试集数据data_04.csv，舍选法产生的点

代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数。

## Import

主要用到numpy、time、pandas、matplotlib包