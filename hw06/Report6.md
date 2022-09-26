# Report6
 
1. 由于可以找到F(x)>p(x)始终成立，（切线法）

$$
F(x)=\left\{
\begin{aligned}
&F_1(x)=k_1(x+1),&x\in[-1,0)\\
&F_2(x)=k_2x+b, &x\in[0,1]
\end{aligned}
\right.
$$

其中$k_1=0.824361, k_2=0.5, b = 0.43566$


2. 通过分段分别处理，随机选择两个在[0，1]之间均匀分布的随机抽样 $(\xi_1,\xi_2)$

$x\in[-1,0)$时的分布应当满足：
$$
\xi_1 = \frac{\int^x_{-1}F_1(t)dt}{\int^0_{-1}F_1(t)dt}=(1+x)^2,y=\xi_2F_1(x)\\
x=1-\sqrt{\xi_1}
$$
$x\in[0,1]$时的分布满足：
$$
\xi_1 =\frac{\int^x_{0}F_2(t)dt}{\int^1_0F_2(t)dt},y=\xi_2F_2(x)\\
x=\sqrt{\xi_1(1+\frac{2b}{k_2})+\frac{b^2}{k^2_2}}-\frac{b}{k_2}
$$

3.   判断条件：

$$
y<p(x)
$$

4.  否，则舍去；是，则取 $x$