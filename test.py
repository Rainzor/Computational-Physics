# from myStats import norm,binom,bernulli,expon,poisson,Cos
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd

# df = pd.read_csv("data_09.csv")

# if __name__=="__main__":
#     N =10000
#     n = Cos()
#     r = n.rvs(N)
#  #绘制直方图与概率密度分布图像
#     num_bins =50


#     #绘制直接抽样法得到的图像
#     fig, ax = plt.subplots()

#     #the histogram of the data
#     ax.hist(r, bins =num_bins, density=True, histtype='stepfilled', alpha=0.5,range=(-2,2))
    
#     x = np.linspace(n.ppf(0.001), n.ppf(0.999), 100)
#     ax.plot(x, n.pdf(x),'r-', lw=5, alpha=0.6, label="Best fit")
#     # add a 'best fit' line
#     ax.set_xlabel('x')
#     ax.set_ylabel('Probability density')

#     # Tweak spacing to prevent clipping of ylabel
#     ax.legend(loc='best', frameon=False)
#     fig.tight_layout()
#     plt.show()



