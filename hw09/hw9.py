from myStats import norm,binom,expon,poisson,uniform,Cos
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Schrage_16807 import*


def stats_x(r, size , n):
    """statistic x of n nubers using a distribution

    Parameters
    ----------
    n : int
        numbers of statistic 
    size:
        REQURIE: size%n=0

    r : Array
        random array
    """
    m = size//n
    if size%n!=0:
        raise ValueError("size%n not eq to 0")
    else:
        t = r.reshape(m,n)
        #print(t)
        return t.mean(axis=1)


def plot_big_data(SIZE, d, **kwargs):
    d_name = kwargs["name"]
    d_type = kwargs["type"]
    r = d.rvs(SIZE)
    d_mean, d_var = d.stats

    #标准正态分布
    gaussian = norm()
    x = np.linspace(-10, 10, 200, dtype=float)
    y_gaussian = gaussian.pdf(x)
    # df = pd.DataFrame({
    #     "x":x,
    #     "y":y_gaussian
    # })
    # df.to_csv("data_gaussian.csv")


    #计算统计量，并画图
    N = (2, 5, 10, 50)
    i = 0
    fig, ax = plt.subplots(2, 2)
    ax = ax.reshape(4)
    for n in N:
        s_x = stats_x(r, SIZE, n)
        # print(s_x)
        s_x = np.sqrt(n)*(s_x-d_mean)/np.sqrt(d_var)


        # #得到数据
        # try:
        #     df = pd.read_csv("data_09_{}.csv".format(n))
        # except FileNotFoundError:
        #     df = pd.DataFrame({d_name:r})
        #     #df.to_csv("data_09.csv")
        # else:
        #     df[d_name]=r
        # df.to_csv("data_09_{}.csv".format(n), index=False)



        if d_type== 'continuous':
            num_bins = 100
        elif d_type == "discrete":
            num_bins = 40
        else:
            raise ValueError

        #绘制直接抽样法得到的图像
        #the histogram of the data
        nt,bins,patches= ax[i].hist(s_x, bins=num_bins, density=True,
                   alpha=0.5, range=(-10, 10), label=d_name+" Histogram of Data")
        # add a 'best fit' line
        ax[i].plot(x, y_gaussian, '--', label="Standard Gaussian")
        # add a 'best fit' line
        ax[i].set_xlabel('x')
        ax[i].set_ylabel('Probability density')
        ax[i].set_title("N:{}".format(n))
        # Tweak spacing to prevent clipping of ylabel
        ax[i].legend(loc='best', frameon=False)
        print("N:{}".format(n))
        print("nt={}, nt num = {}".format(nt,len(nt)))
        print("bins={}, bins num = {}".format(bins,len(bins)))
        print("patches={}".format(patches))
        i = i+1

    # manager = plt.get_current_fig_manager()
    # manager.window.showMaximized()




if __name__=="__main__":
    SIZE = 100000
    # plot_big_data(SIZE,expon(), type = "continuous",name = "Exponential")
    # plot_big_data(SIZE, poisson(5), type="discrete", name='Poisson')
    # plot_big_data(SIZE,binom(10),type= "discrete",name = "Binomial" )
    # plot_big_data(SIZE,Cos(),type="continuous",name = "Cos" )
    plot_big_data(SIZE,uniform(),type="continuous",name = "Uniform" )
    
    plt.show()

    




