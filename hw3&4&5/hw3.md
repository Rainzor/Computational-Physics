# hw3

设 $p(\theta,\varphi)$是球面上均匀的概率密度分布函数，**以球坐标为系**，半径r=1，那么有
$$
\int pds=\int^{2\pi}_{0}\int^{\pi}_0 p(\theta,\varphi)\sin\theta d\theta d\varphi=1
$$
解出 $p(\theta,\varphi )=1/4\pi$

那么对于 $(\theta,\varphi)$来说，其符合的概率密度函数为
$$
g(\theta,\varphi)=\frac{sin\theta}{4\pi}=\frac{\sin\theta}{2}\times\frac{1}{2\pi}=f(\theta)\times h(\varphi)
$$
对 $(\theta,\varphi)$进行直接抽样法抽样，易知，$\varphi$是 $[0,2\pi]$上均匀分布

而 $\theta$的累积函数 $F(\theta)$为
$$
F(\theta)=\int^{\theta}_0\frac{\sin t}{2}dt=\frac{1-\cos\theta}{2}
$$
故为了得到 $\theta$的抽样值:

需要先得到(0，1）上均匀分布的 $\xi_1$代入 $\theta = \arccos(1-2\xi_1)$，考虑到随机数性质，这等价于 $\theta = \arccos(2\xi_1-1)$

