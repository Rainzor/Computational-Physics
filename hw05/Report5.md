# Report5

> Rainzor

### 1. Question

​	对于球面上均匀分布的随机坐标点，给出它们在(x, y)平面上投影的几率分布函 数。并由此验证 Marsaglia 抽样方法$x = 2u\sqrt{1-r^2},y = 2v\sqrt{1-r^2}, z=1-r^2$ 确为球面上均匀分布的随机抽样

### 2. Algorithm

#### 2.1 Marsaglia抽样方法

​	三维球面上分布的 Marsaglia 方法为：

1. 随机抽样一对均匀分布的随机数，$(u,v)\in[-1,1]$
2. 计算 $r^2=u^2+v^2$ ，如 果 $r^2>1$ 则重新抽样直至 $r^2\le1$ ；
3. 得三维坐标抽样为：

$$
x=2u\sqrt{1-r^2},y=2v\sqrt{1-r^2},z=\sqrt{1-x^2-y^2}=(1-2r^2),r^2=u^2+v^2
$$

####   2.2 证明

​	为证明抽样结果确为球面均匀分布。

​	首先，Marsaglia抽样方法所得关于 $(u,v)$的联合分布满足为:在单位圆内均匀分布
$$
\int p(u,v)dudv=\int_{u^2+v^2\le1}Cdudv=1
$$
​	$C=1/\pi$， 所以 $(u,v)$满足联合密度分布为：

$$
p(u,v)=1/\pi
$$
​	可以得到（u，v）坐标系到（x，y）坐标系的变换关系
$$
\begin{aligned}
\frac{\part(u,v)}{\part(x,y)}&=(\frac{\part(x,y)}{\part(u,v)})^{-1}\
=\left|\begin{array}{cc} 
    \part x/\part u &\part x/\part v     \\ 
    \part y/\part v&\part y/\part v    
\end{array}\right|^{-1}\\
&=\left|\begin{array}{cc} 
    2\sqrt{-r^2+1}-2\frac{2u^2}{\sqrt{1-r^2}} &-\frac{2 u v}{\sqrt{-r^2+1}}     \\ 
    -\frac{2 u v}{\sqrt{-r^2+1}}&   2 \sqrt{-r^2+1}-\frac{2 v^2}{\sqrt{-r^2+1}}
\end{array}\right|^{-1}\\
&=\frac{1}{4(1-2r^2)}
=\frac{1}{2\sqrt{1-x^2-y^2}}
\end{aligned}
$$
​	所以可以得到$(x,y)$概率密度函数
$$
p(x,y)=\frac{\part(u,v)}{\part(x,y)}p(u,v)=\frac{1}{2\pi\sqrt{1-x^2-y^2}}
$$

​	直接看球面上均匀分布，在球坐标系中的联合分布函数：
$$
\int p(\theta,\varphi)d\theta d\varphi=\int^{2\pi}_0\int^{\pi}_0 C\sin{\theta}d\theta d\varphi=1
$$
​	$C=1/4\pi$，所以 $(\theta,\varphi)$满足联合密度分布为
$$
p(\theta,\varphi)=\frac{\sin \theta}{4\pi}
$$
​	根据关系
$$
x=\sin\theta\cos\varphi ,y=\sin\theta\sin\varphi
$$
​	同理得到(x,y)概率密度函数
$$
p(x,y)=\frac{\part(\theta,\varphi)}{\part(x,y)}p(\theta,\varphi)=\frac{1}{2\pi\sqrt{1-x^2-y^2}}
$$
​	从不同角度出发，得到相同的概率的密度，故 Marsaglia 抽样方法得到的确实是球面上的均匀分布。

### 3 Experiment

​	在实验中，我们利用Marsaglia方法，绘制出如下图像

<img src="F:\MyDocuments\Physics\Computational Physics\Homework\hw05\三维球面.png" style="zoom: 50%;" />

​	 左图可见得到的确实是均匀分布的球面。且可以看到，在XY平面内投影并不均匀，从而佐证了原图像是球面均匀分布。

### 4 Summary

​	本实验了解另一种在球面上采样的方法，Marsaglia抽样方法。且其计算效率高于直接的三角函数法，为后续实验打下基础
