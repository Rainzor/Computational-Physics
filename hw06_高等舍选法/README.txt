本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

主程序在“hw6.py”中。

压缩包中提供了测试集数据“data_06.csv”，存储着标准正态分布的采样点。
代码中没有直接调用，而是采取“时间-种子”方式，现场形成随机数。

主要用到numpy、time、pandas、matplotlib包