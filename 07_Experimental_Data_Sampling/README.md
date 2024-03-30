# README

## Code

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

主程序在“hw7.py”中

## Data

“data.TXT”为实验原始数据

提供了测试集数据“data_07.csv”，存储直接法与舍选法得到的采样数据。

代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数。

## Explanation

代码对实验所得离散分布所形成的概率密度分布函数进行的类的封装“MyDistribution”，其中含有自己构建的常用关于分布的方法，与标准`scipy.stats`一致

主要用到numpy、time、pandas、matplotlib包