# README

## Question

Monte Carlo 方法研究二维平面上荷电粒子在正弦外电场 $( sin\omega t )$ 中的随机行走。 推导速度自相关函数的表达式。它随时间的变化是怎样的行为？能否模拟得到该自相关函数的曲线?是的话与理论曲线进行比较，否的话讨论理由。

## Code

主程序在“hw10.py”中

主要用到numpy、time、pandas、matplotlib包

## Data

data中提供了测试集数据“Brown1.csv”,“Brown2.csv”,“cov_velocity1.csv”,“cov_velocity2.csv”

其中含有粒子随机游走的位置图像数据、速度相关函数理论与模拟数据

## Explanation

代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数，可以直接调用程序即可画图

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

内部包含“seed_time”随机数产生种子函数，“Schage16807”类

