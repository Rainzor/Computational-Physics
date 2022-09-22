# Hw4

对一维随机变量，设 $\xi$是[0, 1]上的均匀分布，对x的分布函数 $p(x)$满足
$$
p'(x)=a\delta(x)+b\exp(-cx),x\in[-1,1],a\neq0\\
p(x) = \int_{-\infin}^x p'(t)dt
$$
​                                                                                            

由于 $\delta(x)$积分后得到的 阶梯函数$H(x)$在x=0处不连续，所以对p(x)分段处理
$$
p(x)=\left\{
\begin{aligned}
&\frac{b}{c}(\exp{c}-\exp{(-cx)})&-1\le x<0\\
&a+\frac{b}{c}(\exp{c}-\exp{(-cx)})&0\le x\le1
\end{aligned}
\right.
$$
   概率密度要求非负且积分归一，且因为p(x)分段单调，那么分别有

1. $$
   p(-1)=0\ge0
   $$

2. $$
   p(0-)=b/c*(\exp c-1)\ge0
   $$

3. $$
   p(0+)=a+b/c*(\exp c-1)\ge0
   $$

4. $$
   p(1)=a+b/c*(\exp c-\exp{(-c)})\ge0
   $$

5. $$
   P(x)=\int vp(x)dx=(\int^{0}_{-1}+\int^1_0) p(x)dx=a+\frac{2b(c\exp c-\sinh{c})}{c^2}=1
   $$

由上面的条件可得，一个关于a,b,c的**充分条件**
$$
a = 1- \frac{2b}{c}e^c+\frac{2b}{c^2}\sinh c\\
b>0\\
c>0\\
1-\frac{b}{c}(1+\exp{c})+\frac{2b}{c^2}\sinh c\ge0\\
1-\frac{2b}{c}(\cosh c)+\frac{2b}{c^2}\sinh c\ge0
$$
可以取 
$$
b=0.5,c=0.5,a=1-\frac{2b\exp{(-c)}}{c}\approx-0.213061
$$


考虑到p(x)的不连续性，所以无法写出累积函数的解析表达式，我们采取舍选法：

1. 设p(x)在$[-1，0)$上界$p(0-)<M_1=0.65$；$[0，1]$上界$p(1)<M_2=0.83$

2. 随机选择两个在[0，1]之间均匀分布的随机抽样 $(\xi_1,\xi_2)$

3. 判断条件：
   $$
   \begin{aligned}
   &M_1\xi_2\le p(-1+2\xi_1)&-1+2\xi_1\in[-1,0)\\
   &M_2\xi_2\le p(-1+2\xi_1)&-1+2\xi_1\in[0,1]
   \end{aligned}
   $$
   是否成立

4. 否，则舍去；是，则取 $x=-1+2\xi_1$

