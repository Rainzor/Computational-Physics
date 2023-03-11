# Note3 

如果一个物理量的表达式含有物理量其本身，即
$$
x=f(x)
$$
则求解这个物理量的通常采用数值计算上的迭代法

## 3.1计算物理量的迭代方法

#### 1. 直接迭代法

迭代法：
$$
x_{n+1}=f(x_n)
$$
直到，满足相对误差精度 $\epsilon$
$$
\left|\frac{x_{n+1}-x_n}{x_n}\right|<\epsilon
$$

#### 2. 牛顿迭代

利用一阶导数提高精度
$$
g(x)=f(x)-x\\
x_{n+1}=x_n-\frac{g(x_n)}{g'(x_n)}
$$
除了要满足相对误差，还要有绝对误差精度 $\delta$
$$
|g(x_n)|<\delta
$$

#### 3 混沌输入迭代

即输入变量步数前一步，而是前面很多步，比如
$$
x'_{n}=ax_n+(1-a)x_{n-1}\\
f(x_n')=x_{n+1}
$$

#### 4 多元变量

多变量问题
$$
{\phi_i(x)}=f({\phi_i(\bold x)})
$$

## 3.2 混沌

混沌是一种在确定性系统（决定论）中出现的类似随机的过程

系统的本质是确定的，对初值是敏感的

#### 1 一维迭代Logistic方程

$$
f(x)=\lambda x(1-x),\quad (0\le x\le1,0\le\lambda\le4)
$$

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221121102253264.png" alt="image-20221121102253264" style="zoom:67%;" />

#### 2 初值敏感性

初值敏感性不是由任何外界随机因素引起的，而是由确定性系统内部产生的。

初值敏感性是混沌中类似随机现象的根源，因为很难对初值确定得非常精确，于是近似 相同的初值给出了很不相同的貌似随机的结果。

<img src="C:\Users\Lenovo\AppData\Roaming\Typora\typora-user-images\image-20221121104125103.png" alt="image-20221121104125103" style="zoom:67%;" />

#### 3 混沌的性质

- 内随机性
- 分维性
- 普适性和Feigenbaum 常数

#### 4 二维迭代 $Henon$方程

它是 Logistic 方程的推广
$$
x_{n+1}=1-ax^2_n+y_n\\
y_{n+1}=bx_n
$$

#### 5 吸引子和奇异吸引子

#### 6 $Ro$