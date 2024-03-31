# README

## Question

​用 Metropolis-Hasting 抽样方法计算积分：

$$
I=\int_0^{\infin}(x-\alpha \beta)^2f(x)dx=\alpha \beta^2
$$

​其中 $f(x)=\frac1{\beta\Gamma(\alpha)}(\frac x \beta)^{\alpha-1}\exp(-x/\beta)$  。

设权重函数为： $p(x)=f(x)$和 $p(x)=(x-\alpha \beta)^2f(x)$ 。给定参数 $\alpha,\beta$ 用不同的 $\gamma$ 值计算积分，讨论计算精度和效率。


## Code

主程序在“code/Metroplis.py”中

程序中含有”MetroSample“类，使用Metroplis-Hastings 抽样算法特定的分布

类中包含有各种方法，详细见注释

主要用到numpy、time、pandas、matplotlib包

”Metroplis.py“ 中写好了 `plot_data(filename)`函数，您可以对给定的`data_error_rate_gamma_1.csv`,`data_error_rate_gamma_2.csv`文件，调用函数直接绘制图像

## Data

data中提供了实验数据
data_error_rate_gamma_1.csv,data_error_rate_gamma_2.csv为随着 gamma改变，误差与效率的变化

data_best_gamma_1.xlsx, data_best_gamma_2.xlsx则是不断改变 alpha beta值，找各个最优gamma的结果

## Explanation

主程序中只保留了两种情况积分的输出结果。

注释了很多其他的代码，是报告中实验实现的具体细节，您可以依次取消注释，自行复现实验结果，调用程序即可画图


本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。
内部包含“seed_time”随机数产生种子函数，“Schage16807”类，故代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数，

