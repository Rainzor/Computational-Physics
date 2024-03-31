# README

## Question

用16807产生器测试随机数序列中满足关系 $X_{n-1} > X_{n+1} > X_{n}$。

讨论Fibonacci延迟产生器中出现这种关系的比重。


## Code

本实验对16807Schrage产生器进行单独封装，为的是以后实验方便继续调用。其代码在“Schrage_16807.py”中

主程序在“hw2.py”中，其中对Fibonacci产生器进行的类的封装

## Data

提供了测试集数据“test_dataset.csv”，但实验中没有调用，实验采取“时间-种子”方式，现场形成随机数。

为了保证实验结果的可靠性，进行了10次迭代，每次迭代程序会暂停4秒钟，最后打印出结果。

## Import

主要用到queue、numpy、time、pandas包
