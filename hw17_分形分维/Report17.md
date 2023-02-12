# Report17

> 姓名：王润泽
>
> 学号：PB20020480

## 1 Question

​	进行单中心DLA模型的模拟(可以用圆形边界，也可 以用正方形边界)，并用两种方法计算模拟得到的DLA图形的分形维数，求分形维数时需要作出双对数图

## 2 Algorithm

### 2.0 DLA

​	DLA模型的生成方式在作业11中已经实现过了，大致步骤是：

- 取一个2维的方形点阵，在点阵中央原点处放置一个粒子作为生长的种子。
- 然后从距原点足够远的圆周界处释放一 个粒子，让它作随机行走
- 该粒子走到种子的最近邻位置与种子相碰， 这时让粒子粘结到种子上不再运动；
- 当粒子走到点阵边界，这时认为粒子走了一条无用的轨迹，取消该粒子，重新生成新的粒子。
- 因此， 那些有用的粒子与种子相粘结后形成不断生长的聚集集团。

### 2.1 Sandbox Method

​	Sandbox来计算分形位数的方法，公式为
$$
N(r)\sim r^D\\
D = \frac{\ln N}{\ln r}+C
$$
​	N为方盒子中的像素点数，r为盒子的边长，C为其他常数，统计过程如下图所示

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221207221203711.png" alt="image-20221207221203711" style="zoom:50%;" />

​	按照 $\sqrt 2$幂次增长盒子的边长，统计内部点数，得到对数坐标图形，计算斜率即为分形位数D

### 2.3 密度-密度相关函数法

​	平面上分形图形的密度-密度 相关函数定义:
$$
C(r)=\left<\sum\frac{\rho(r')\rho(r'+r)}{N}\right>\sim r^{-\alpha}
$$
​	$\rho(r’)$ 是图形的密度函数，有图形象素处为 1，无 图形处为 0；N 是总象素数。C(r) 的几何意义是：原始图形和平移 r 后的图形重叠部分的象素数和全部象素数的比值，即在相距 r  处发现另一象素的概率。在实际计算中除了对 N 个不同的 r’ 求平均之外，还对不同方向、 长度相同的 r 求平均

​	特例: 固定 r’， 并把它取为图形的中心，即r’=0，此时
$$
C(r)=\left<\sum\rho(0)\rho(r)\right>\sim r^{-\alpha}
$$
则要做的就是按照 $\sqrt 2$幂次增长圆的半径，统计内部重叠点数，绘制出对数坐标，得到斜率 $\alpha$

​	在C(r) 在回转半径R内积分，R 足够大时，积分值很接近于和图形总象素数 N 成正比
$$
\int_0^RC(r)\text{d}^{dim}r\sim N\\
N\sim R^{dim-\alpha}\\
D=dim-\alpha
$$

 dim是欧氏空间维数，此题中为2；D是分形维数。

## 3 Experiment

### 3.0 DLA

​	实验画出的图像如下

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw17\DLA.png)

### 3.1 Sandbox Method

​	实验结果如下图所示

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw17\Sandbox.png)	

​	实验计算得到图像的斜率为
$$
D_1=1.6666788360254807
$$
该斜率即为分形维数。

​	分形维数结果在 $1.6\sim1.7$之间，与理论比较符合。

<div STYLE="page-break-after: always;"></div>

### 2.3 密度-密度相关函数法

​	实验结果如下图所示

![](F:\MyDocuments\Physics\Computational Physics\Homework\hw17\Density.png)

​	图像斜率与分形维数是：
$$
k=-\alpha  = -0.35808654164728
\\D_2=d-\alpha=1.64191345835272
$$
​	分形维数结果在 $1.6\sim1.7$之间，与理论比较符合。

## 4 Summary

​	本次实验利用 DLA 的生长规则模拟出了 DLA 的轨迹图形。同时用 sandbox 法和密度-密度相关函数法求出了 DLA 的维数，两个方法所求出来的维数都在1.6~1.7内和理论值较为贴近。