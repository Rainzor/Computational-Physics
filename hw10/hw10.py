from cProfile import label
from statistics import mean
from matplotlib.lines import lineStyles
from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt


def plot_brown(x):
    #绘制位移图像———加电场的布朗运动
    plt.style.use('classic')
    fig, ax = plt.subplots()
    #颜色的渐变
    num_points_range=range(len(x))
    #绘制点———云彩
    ax.scatter(x[:,0], x[:,1],
              c=num_points_range,
              cmap='viridis',
              edgecolors='none',s=3)
    
    ax.plot(x[:,0],x[:,1],linewidth=0.3)
    #起点与终点
    ax.scatter(x[0,0], x[0,1],c='red',edgecolors='none',s=100,label='start')
    ax.scatter(x[-1,0], x[-1,1],c='yellow',edgecolors='none',s=100,label='end')
    ax.legend()
    #隐藏坐标轴
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title("random walk")

#计算Bsin
def Bsin(x,T=10):
    return 1*np.sin(x*2*np.pi/T)

#数值计算相关函数
def cov_velocity(r,h,tau,T,a=0.01):
    t = np.arange(start=0,stop=T,step=h,dtype=float)
    N = len(t)
    #初始化速度和位移
    v = np.zeros((N,2),dtype = float)
    x = np.zeros((N,2),dtype = float)
    v[0] = r.rand(2)*2-1

    #初始化随机力
    A = r.rand(N,2)
    A = (A*2-1)*a
    F = np.zeros((N,2),float)
    F[:,0] = A[:,0]+Bsin(t,T = 1)
    F[:,1] = A[:,1]
    #print(F)
    #迭代计算速度和位移
    
    k = -h/tau+1
    for i in range(0,N-1):
        v[i+1] = k*v[i] + h*F[i]
        x[i+1] = x[i] + v[i]*h        
        #print(v0)
    v0_t = v[0]
    #计算速度相关函数
    cov_v = (v @ v0_t)/2

    return x,v,cov_v

#理论相关函数
def cov_theory(t,c):
    return c*np.exp(-t/tau)

def plot_cov(t,cov_v_theory,cov_v_sim):
    #绘图
    plt.plot(t, cov_v_theory, label='theory')
    plt.plot(t, cov_v_sim, 'r--', label='simulation')
    plt.title('covariance of velocity')
    plt.xlabel('t(s)')
    plt.ylabel('covariance')
    plt.legend()

if __name__=="__main__":
    #初始化
    r = Schrage16807(seed_time())
    n = 50

    #较大的随机力
    # a = 10

    #较小的随机力
    a = 0.01

    #较大的粘滞阻力
    #tau = 1

    #较小的粘滞阻力
    tau = 5
    T = 10*tau
    h = tau/1000

    #绘制粒子运动图像
    x,v,cov_v = cov_velocity(r,h,tau,10*T,a) 
    plot_brown(x)
    
    t = np.arange(start=0, stop=T, step=h, dtype=float)#时间

    # #保存数据
    # df = pd.DataFrame({'x':x[:,0],'y':x[:,1]})
    # df.to_csv('Brown.csv',index=False)

    #计算相关函数
    #迭代n次，计算平均值
    cov_v = 0
    for i in range(n):
        x_t,v_t,cov_v_t = cov_velocity(r,h,tau,T,a)
        cov_v = cov_v + cov_v_t
    cov_v_sim = cov_v/n#模拟值

    #统计数据
    c = cov_v_sim[0]
    cov_v_theory = cov_theory(t,c)#理论值

    # #保存数据
    # df = pd.DataFrame({'t':t,'cov_v_theory':cov_v_theory,'cov_v_sim':cov_v_sim})
    # df.to_csv('cov_velocity.csv',index=False)

    plt.figure()

    plot_cov(t,cov_v_theory,cov_v_sim)

    plt.show()