from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


if __name__=='__main__':
    #define some functions
    PI = np.pi  
    lorentz = lambda x:1.01/(np.sqrt(2*PI)*(1+0.25*x**4))
    gauss = lambda x:1/np.sqrt(2*PI)*np.exp(-(x**2)/2)
    qdf = lambda x: 1/(PI*(1+x**2))
    h_x = lambda x: lorentz(x)/qdf(x)

    N = 200000
    M = 2.05
    seed = seed_time()
    r = Schrage16807(seed)
    xi_1 = r.rand(N)
    xi_2 = r.rand(N)
    xi_1 = PI*xi_1 - 1/2
    # from -infinity to +infinity,实际上infinity=1.63e+16
    xi_1 = np.tan(xi_1)

    #为了防止范围过大，且|x|>5之后概率密度非常小，故将超过这部分的值进行替代
    xi_1 = np.clip(xi_1,-5,5)

    #对x1 lorentz分布进行舍选
    x_1 = [x for x,y in zip(xi_1,xi_2) if  M*y < h_x(x)]
    sample_rate1 = len(x_1)/N;
    x_1 = np.array(x_1)

    xi_3 = r.rand(len(x_1))
    y_1 = xi_3*lorentz(x_1)
    #对x_result gauss分布进行舍选
    x_result = [x for x,y in zip(x_1,y_1) if y<gauss(x)]
    sample_rate2 = len(x_result)/len(x_1);

    #输出采样效率
    print("First rate of sampling：",sample_rate1)
    print("Second rate of sampling：",sample_rate2)

    #保存数据
    df = pd.DataFrame(data=x_result)
    df.to_csv("data_06.csv")

    #绘制直方图与概率密度分布图像
    fig, ax = plt.subplots()
    num_bins = 100
    #the histogram of the data
    n, bins, patches = ax.hist(x_result, num_bins, density=True,label='histogram of data')

    # add a 'best fit' line
    y = np.array([gauss(x_i) for x_i in bins])
    ax.plot(bins, y, '--',label="Gauss")
    ax.set_xlabel('X')
    ax.set_ylabel('Probability density')
    # Tweak spacing to prevent clipping of ylabel
    ax.legend()
    fig.tight_layout()
    plt.show()


