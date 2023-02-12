from scipy.optimize import curve_fit
from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#多重积分的Mente Carlo方法
def monte_carlo_iint(f,domain_func,N=1):
    """Using simple monte carlo method to calculate the numerical integral

        Parameters:
            f(function): the function will be integrated
            domain_func(list of trulpe): the domain of the function e.g. [(a,b)]; [(a,b),(c,d)].
            N(int): the number of sampling
        
        Returns:
            intergral: the result of intergral
    """
    seed = seed_time()
    r = Schrage16807(seed)
    d = len(domain_func)
    x = r.rand(d,N)
    i = 0
    v = 1
    for start,stop in domain_func:
        x[i] = x[i]*(stop-start)+start
        v = v*(stop-start)
        i = i+1
    y = f(x)
    y_mean = y.mean()
    integral = v*y_mean
    return integral

#被积函数
def func(variable):
    """the multivariate function

    Parameters:
        variable(float):the function has 5 varitates, so x is either a vector(5*1) or a matrix(5*N)
    
    Returns:
        f:a number or a vertor(1*N)
    
    """
    x = variable[0]
    y = variable[1]
    z = variable[2]
    u = variable[3]
    v = variable[4]
    return 5+x**2-y**2+3*x*y-z**2+u**3-v**3


if __name__ == "__main__":
    #定义域
    domain_func = [(0,7/10),(0,4/7),(0,9/10),(0,2),(0,13/11)]

    #数点
    N = [2**k for k in range(5, 21)]
    #数值积分
    mc_int = np.array([monte_carlo_iint(func,domain_func,n) for n in N])
    #误差
    mc_error = abs(mc_int-5.67712092042336)  # 精确值

    df = pd.DataFrame({
        "N":N,
        "Integral":mc_int,
        "Error":mc_error
    })
    print(df)
    #df.to_csv("data_08_2.csv")

    #曲线拟合
    trend_f = lambda x,c:c/np.sqrt(x)
    popt, pcov = curve_fit(trend_f, N, mc_error)
    c = popt[0]

    #绘图
    k = np.linspace(5,20,1000)
    e = c/(2**(k/2))
    plt.plot(range(5,21),mc_error,label="Monte Carlo int error")
    plt.plot(k,e,label=r"$O(\frac{1}{\sqrt{N}})$")
    plt.xlabel(r"$\log_2{N}$")
    plt.ylabel("Error")
    plt.legend()
    plt.title("The error trend in Monte Carlo Method")
    plt.show()
