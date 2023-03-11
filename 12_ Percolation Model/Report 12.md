# Report 12

> Rainzor

## 1. Question

推导正方格子点阵上键逾渗的重整化群变换表达式 $p’ = R(p)$，求临界点 $p_c$ 与临界指数 $\nu$ ，与正确值（下表）相比较。

<img src=".\img\table.png" alt="image-20221025001847845" style="zoom:50%;" />

<center><p>表1

## 2. Answer

### 2.1 临界点 $p_c$

​	本题要讨论的是正方形点阵上的键逾渗，而键逾渗有6条相邻的键，如下图所示

<img src=".\img\键逾渗.png" alt="键逾渗" style="zoom: 25%;" />

<center><p>图1：方阵键逾渗

​	且与两边有两个点联通，与上下有四个点联通，但实际上这个结构与三角形的座逾渗结构十分相似，如下图所示：

<img src=".\img\三角点阵.png" alt="三角点阵" style="zoom: 25%;" />

<center><p>图2：三角阵座逾渗

​	与图1的结构一样，中心的点有6个相邻的座，二者是同构。

​	那么为了方便问题的讨论，我们就可以把问题转化成：讨论三角格子点阵上座逾渗的临界点 $p_c$ 与临界指数 $\nu$。

​	我们重整化方式如下，取尺度放大因子为2的元胞。

<img src=".\img\4.png" alt="4" style="zoom:25%;" />
<center><p>图3：b=2元胞


​	一共有三种情况可以联通

1. 四个点都占据

   <img src=".\img\4.png" alt="4" style="zoom:25%;" />
<center><p>图4：四个点占据

2. 三个点占据

<center><half>
   <img src=".\img\3_1.png" alt="3_1" style="zoom:10%;" />
	<img src=".\img\3_2.png" alt="3_2" style="zoom:10%;" />
    <img src=".\img\3_3.png" alt="3_3" style="zoom:10%;" />
    <img src=".\img\3_4.png" alt="3_3" style="zoom:10%;" />


<center><p>图5：三个点占据
3. 两个点占据
<center><half>
​	<img src=".\img\2_1.png" alt="2_1" style="zoom:15%;" />
    <img src=".\img\2_2.png" alt="2_2" style="zoom:15%;" />
     <img src=".\img\2_3.png" alt="2_3" style="zoom:15%;" />
<center><p>图6：两个点占据


​	假设格点的占据几率为 $p$，对于 **b = 2** ，上下端连接的图形有7个（图1.6.3.3-1），其变换表达式为
$$
p' = R(p|b=2)=p^4+4p^3(1-p)+3p^2(1-p)^2 \tag1
$$
​	一般来说，重整化后的格子点阵占据几率 p ' 相异于原格子点阵的占据几率 p，但对于零界点 $p_c$，它满足关系式：
$$
p_c= R(p_c)\tag2
$$
​	由（1）（2）解出非平凡的值为 $\bold{p_c=0.5}$，与表1对比可以看出，与正确值一样

### 2.2 临界指数 $\nu$

​	为了计算临界指数 $\nu$，对（1）式求导
$$
R'(p)=4p^3+12p^2(1-p)-4p^3+6p(1-p)^2-6p^2(1-p)=6(1-p)p\tag3
$$
​	那么在 $p_c$处的值为 $R’(p_c) = 1.5$，那么可以解得临界指数为
$$
\nu=\frac{\ln b}{\ln k}=\frac{\ln2}{\ln 1.5}\approx1.7095
$$
​	与正确的值对比 $p^*=0.5,\nu^*=4/3$,可见 $p_c$求得了精确值，而 $\nu$的值有所差距，但对于b=2的简单情况来说，近似结果已经不错了。可能对于较大的b，结果可能会得到更好的改善。



​	
