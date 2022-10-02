#Using UTF-8 Code Python
#Runze Wang
#PB20020480

from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#概率分布函数 Probability density function.
def pdf(x):
    b = 0.5
    c = 0.5
    a = 1 - np.exp(-c)*2*b/c
    if -1<=x and x<0:
        return b/c*(np.exp(c)-np.exp(-c*x))
    elif 0<=x and x<=1:
        return a + b/c*(np.exp(c)-np.exp(-c*x))
    else:
        return 0

#舍选法，选择合适的数据
def rejection_method(d):
    M1 = 0.65
    M2 = 0.83
    x=d[0]
    y=d[1]
    if -1<=x and x<0:
        return M1*y<=pdf(x)
    else:
        return M2*y<=pdf(x)


if __name__=='__main__':
    M1 = 0.65
    M2 = 0.83
    N = 100000
    P = M1/(M1+M2)
    seed = seed_time()
    r = Schrage16807(seed)
    xi = np.zeros((N,2))
    #生成二维均匀随机变量

    #按比例抽样
    for i in range(N):
        t = r.rand()
        if 0<t and t<P:
            xi[i,0] = r.rand()-1
        else:
            xi[i,0] = r.rand()

    xi[:,1] = r.rand(N,1)
    #print(xi[:30,0])d

 

    #舍选法
    result = [d[0] for d in xi if rejection_method(d)]
    result = np.array(result)

    df = pd.DataFrame(data=result)
    df.to_csv("data_04.csv")

    #绘制直方图与概率密度分布图像
    fig, ax = plt.subplots()
    num_bins = 100
    #the histogram of the data
    n, bins, patches = ax.hist(result, num_bins, density=True,label='histogram of data')

    # add a 'best fit' line
    y = np.array([pdf(x_i) for x_i in bins])
    ax.plot(bins, y, '--',label="best fit")
    ax.set_xlabel('X')
    ax.set_ylabel('Probability density')
    # Tweak spacing to prevent clipping of ylabel
    ax.legend()
    fig.tight_layout()
    plt.show()



    

