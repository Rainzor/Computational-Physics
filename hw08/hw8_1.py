from scipy.optimize import curve_fit
from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#返回N个分布为p的随机数
def rvs(start,stop,N,ppf,M):
    """
    random variate 

    Parameters:
        start(float): the start of distribution
        end(float): the end of distribution
        N(int): the number of sampling
        p(function): distribution of probility function
        M(float): ceiling bound 

    Returns:
        random variate
    """
    a = start
    b = stop
    seed = seed_time()
    r = Schrage16807(seed)
    x = np.zeros(N)
    i = 0
    while i<N:
        xi_x = r.rand()*(b-a)+a
        xi_y = r.rand()
        if M*xi_y<=ppf(xi_x):
            x[i] = xi_x
            i = i+1
    return x
#重要抽样法
def weight_monte_carlo_int(f,start,stop,N,ppf,M):
    """Using weight monte carlo method to calculate the numerical integral

        Parameters:
            f(function): the function will be integrated
            start(float): the start of integral
            end(float): the end of integral
            N(int): the number of sampling
            p(function): weight function
            M(float): the ceil limit for weight function
        Returns:
            intergral: the result of integral
    """
    x = rvs(start,stop,N,ppf,M)
    y = f(x)/ppf(x)
    integeral = y.mean()
    return integeral

#简单抽样法（只实现，未调用）
def simple_monte_carlo_int(f,start,stop,N=1):
    """Using simple monte carlo method to calculate the numerical integral

        Parameters:
            f(function): the function will be integrated
            start(float): the start of integral
            end(float): the end of integral
            N(int): the number of sampling
        
        Returns:
            intergrate: the result of intergration
    """
    seed = seed_time()
    r = Schrage16807(seed)
    x = r.rand(N)
    x = x*(stop-start)+start
    y = f(x)
    y_mean = y.mean()
    integral = (stop-start)*y_mean
    return integral

if __name__=="__main__":
    a = 0
    b = 5
    M = 0.363
    #被积函数
    func = lambda x:np.sqrt(x**2+2*np.sqrt(x))
    #概率密度分布函数
    ppf = lambda x:(0.567944 + x)/15.3397
    
    N = [2**k for k in range(5,21)]#以2的幂次生成
    mc_int = np.array([weight_monte_carlo_int(func,a,b,n,ppf,M) for n in N])#计算数值积分
    mc_error = abs(mc_int-15.4390107355675)#精确值与Monte Carlo方法误差

    df = pd.DataFrame({
        "N":N,
        "Integral":mc_int,
        "Error":mc_error
    })
    print(df)
    # #导入数据
    df.to_csv("data_08_1.csv")

    #曲线拟合
    trend_f = lambda x,c:c/np.sqrt(x)
    popt, pcov = curve_fit(trend_f, N, mc_error)
    c = popt[0]

    #画图
    k = np.linspace(5,20,1000)
    e = c/(2**(k/2))
    plt.plot(range(5,21),mc_error,label="Monte Carlo int error")
    plt.plot(k,e,label=r"$O(\frac{1}{\sqrt{N}})$")
    plt.xlabel(r"$\log_2{N}$")
    plt.ylabel("Error")
    plt.title("The error trend in Monte Carlo Method")
    plt.legend()
    plt.show()