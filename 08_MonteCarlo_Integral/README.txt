本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。

主程序在“hw8_1.py”和”hw8_2“中，进行两种类型的积分计算与结果。
程序输出积分结果与误差，并绘制误差趋势图像。

压缩包中提供了测试集数据“data_08_1.csv”，“data_08_2.csv”，存储：样本个数N,积分值Integral，误差Error
代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数。

代码中定义了”simple_monte_carlo_int“、”weight_monte_carlo_int“和”monte_carlo_iint“函数，用来计算数值积分

主要用到numpy、time、pandas、matplotlib、scipy包