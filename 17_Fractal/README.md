# README

## Question

进行单中心DLA模型的模拟(可以用圆形边界，也可 以用正方形边界)，并用两种方法计算模拟得到的DLA图形的分形维数，求分形维数时需要作出双对数图

## Code

主程序在“DLA.py”中

DLA.py: 含有 `Growth_Model` 抽象类，DLA类继承了Growth_Model类

主要用到numpy、time、pandas、matplotlib包

## Data

**code**中提供了测试集数据“DLA.csv”。其中含有生长点的坐标

## Explanation

代码中直接调用数据，没有现场形成输出，因为产生数据时间较长，如果要复现，可以取消代码注释得到数据

可以直接调用程序即可得到分形维数双对数图像和输出数据

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

内部包含“seed_time”随机数产生种子函数，“Schage16807”类

