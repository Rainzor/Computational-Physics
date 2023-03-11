from Schrage_16807 import*
from math import gamma as gammaFunc
import matplotlib.pyplot as plt
import pandas as pd

# Metroplis-Hastings 抽样算法
class MetroSample(object):
    """The matropolis sampling method

    Parameters
    ----------
    pdf : function
        The probility distribution to be sampled
    **kwargs : dict
        The parameters of pdf and T :alpha beta, gamma

    Methods
    -------
    Tx : function
        The proposal distribution

    Ax :
        The acceptance probability

    ppf_Tx : function
        The percent point function of the proposal distribution
    
    _generate_data : (private)
        Generate N data points from the distribution
    
    sample : (public)
        Sample N points from the distribution

    plot_hist :
        Plot the histogram of the samples
    
    integrate :
        Integrate the hx function from 0 to infinity
    
    """
    def __init__(self,pdf, **kwargs):
        if len(kwargs)==0:
            raise ValueError("No parameters")
        self.pdf = pdf
        self.alpha = kwargs['alpha']
        self.beta = kwargs['beta']
        self.gamma = kwargs["gamma"]
        self.r = Schrage16807(seed=seed_time())
        self.sample_rate = 1
        self.size = 0
        self.markov_chain = None
        self.C = 1 #归一化常数
    # 建议分布
    def Tx(self,x):
        return np.exp(-x/self.gamma)/self.gamma;
    #建议分布累积分布函数的反函数，用来作简单抽样调用
    def ppf_Tx(self,x):
        return -self.gamma*np.log(x)
    # 接受分布
    def Ax(self,x,xt):
        return self.px(xt)*self.Tx(x)/(self.Tx(xt)*self.px(x))
    #待抽样的概率分布（可能未归一化）
    def px(self,x):
        return self.pdf(x,self.alpha,self.beta)

    #核心：抽样方法
    def _generate_data(self, N,x0=1):
        """Generate N data points from the distribution

        Parameters
        ----------
        N : int
            The number of data points to be generated

        Returns
        -------
        data : array
            The data points from the distribution
        """
        if(x0==1):
            m = 2000  # the number of points to be discarded
        else:
            m = 1
        x = np.zeros(N+m)
        x[0] = x0

        n = 1
        s = 0 #the number of points be selected
        while n<N+m:
            R = self.r.rand()
            x_temp = self.ppf_Tx(R)
            r1  = self.Ax(x[n-1],x_temp)
            xi = self.r.rand()
            # 舍选过程
            if xi<min(1,r1):
                x[n] = x_temp
                s += 1
            else:
                x[n] = x[n-1]
            n += 1
        # self.sample_rate = self.sample_rate*(self.size/(self.size+N)) + s/n*(N/(self.size+N))
        self.sample_rate = s/n
        return x[m:]

    #公有接口，对外使用
    def sample(self, N):
        """Sample N points from the distribution

        Parameters
        ----------
        N : int
            The number of points to be sampled

        Returns
        -------
        x : array
            The samples from the distribution
        """
        N = int(N)
        if(self.markov_chain is None):
            self.markov_chain = np.zeros(N)
            self.markov_chain = self._generate_data(N)
            self.size = N
            return self.markov_chain
        elif(self.size<N):
            self.markov_chain = np.concatenate((self.markov_chain,self._generate_data(N-self.size,x0=self.markov_chain[-1])))
            self.size = N
            return self.markov_chain
        else:
            return self.markov_chain[:N]

    
    # 绘制N个采样点的 采样直方图 和 待抽样函数
    def plot_hist(self, N):
        """Plot the histogram of the samples"""
        N = int(N)
        x = self.sample(N)        
        tlim = 8*self.alpha*self.beta
        t = np.arange(0,tlim,0.2)
        bins = len(t)
        ax = plt.subplot(111)
        #绘制直方图
        ax.hist(x, bins=t, density=True,range=(0,tlim),label = "sample")
        #绘制带抽样函数
        ax.plot(t,self.px(t),label = "p(x)")
        ax.set_xlabel("x")
        ax.legend(loc = 0,fontsize=20)
        plt.show()

    # 数值积分计算 hx是被积函数
    def integrate(self,hx,N):
        """Integrate the function hx"""
        N = int(N)
        x = self.sample(N)
        #重要抽样方法（实验1）
        if(hx != 1):
            I1 = np.mean(hx(x, self.alpha, self.beta))
            return I1

        # 比值法（实验2）
        else:
            epsilon = 1e-3
            step = 0.2
            bins = np.arange(0, 8*self.alpha*self.beta, step)
            se = pd.Series(x)
            counts = se.value_counts(bins=bins, normalize=True)
            nts = counts.values/step
            # print(np.sum(nts))
            bins = counts.index
            x = np.array([b.mid for b in bins])
            # plt.bar(x,nts,width=step)
            # plt.show()
            xj = np.array([b.mid for b in bins[nts > epsilon]])
            gj = np.array(nts[nts > epsilon])

            pj = self.px(xj)
            I2 = np.mean(pj/gj)
            self.C  = I2
            return I2


        
# 在一定范围内, 改变gamma, 计算相对误差，采样率，最优gamma
def find_best_gamma(px,hx,N,alpha,beta,proposal_gamma,method = 1):
    """find the best gamma as alpha and beta were set
    
    Parameters
    ----------
    px : function
        The probility distribution to be sampled
    alpha : float
        The parameter of the gamma distribution
    beta : float
        The parameter of the gamma distribution

    Suggestions:
    ------------
    1< alpha < 10

    1< beta < 10

    Returns
    -------
    
    """
    I0 = alpha*beta**2#只适用于题目中的积分结果

    # 合理的构造 gamma 的范围与步长
    if(method == 1):
        gamma_max = proposal_gamma*4
        if proposal_gamma> 5:
            if(gamma_max<80):
                l1 = np.arange(0.2,5,0.2)
                l2 = np.arange(5,gamma_max,0.5)
                gamma = np.concatenate((l1,l2),axis=0)
            else:    
                l2 = np.arange(5, proposal_gamma*2, 1)
                gamma = l2

        else:
            gamma = np.arange(0.2,gamma_max,0.2)
    else:
        gamma_max = proposal_gamma+4.0
        gamma_min = max(1,proposal_gamma-4.0)
        gamma = np.arange(gamma_min,gamma_max,0.1,float)

    l = len(gamma)
    sample_rate  = np.zeros(l)
    relative_error = np.zeros(l)
    max_sample_rate = 0
    best_gamma = 0
    # print(I0)
    for i in range(l):
        m = MetroSample(pdf = px,alpha=alpha,beta=beta,gamma=gamma[i])
        I= m.integrate(hx,N)
        # print(I)
        sample_rate[i] = m.sample_rate
        if(sample_rate[i]>max_sample_rate):
            max_sample_rate = sample_rate[i]
            best_gamma = gamma[i]
        relative_error[i] = np.abs(I-I0)/I0

    df = pd.DataFrame({"gamma": gamma,"sample_rate":sample_rate,"relative_error":relative_error})
    # 返回gamma取值、采样率、相对误差和最佳的 gamma
    return df,best_gamma

#改变alpha和beta的值，可以得到不同最佳gamma值，输出excel文件“data_best_gamma_0.xlsx”
def find_best_gamma_in_range(px, hx, N, a_list, b_list, proposal_gamma, method=1):
    """find the best gamma in a range of alpha and beta

    Parameters
    ----------
    method : int
        1: use the method of in large range of gamma

        2: use the method of in small range of gamma

    """
    g_list = []
    for a in a_list:
        for b in b_list:
            df, best_gamma = find_best_gamma(px, hx, N, a, b, proposal_gamma ,method=method)
            print("a = {}, b = {}, best_gamma = {}".format(a, b, best_gamma))
            g_list.append(best_gamma)
    g_list = np.array(g_list)
    g_list = g_list.reshape((len(a_list), len(b_list)))
    df2 = pd.DataFrame(g_list, index=a_list, columns=b_list)
    df2.to_excel("data_best_gamma_0.xlsx")

# 绘制采样率和相对误差图
def plot_sample_gamma(df,best_gamma,alpha,beta):
    """Plot the sample from the gamma parameter

    """
    gamma = df["gamma"]
    sample_rate = df["sample_rate"]
    relative_error = df["relative_error"]
    index = df[df["gamma"]==best_gamma].index
    max_sample_rate = df.loc[index,"sample_rate"]

    ax = plt.subplot(111)
    l1 =  ax.plot(gamma,sample_rate,'g',label = "sample rate")
    ax.plot(best_gamma,max_sample_rate,'ro')
    ax.annotate("(%d,%f)" % (best_gamma, max_sample_rate), xy=(
        best_gamma, max_sample_rate), xytext=(best_gamma+0.1, max_sample_rate),fontsize=20)
    ax.axvline(best_gamma,ymin=0,ls = '--',color = 'r')
    ax.set_ybound(lower=0)
    ax.set_xlabel(r"$\gamma$",fontsize=20)
    ax.set_ylabel("sample rate",fontsize=20)
    ax.set_title(r"$\alpha$ = {}, $\beta$ = {}".format(
        alpha, beta), fontsize=20)
    ax.set_xbound(lower=0)

    ax2 = ax.twinx()
    l2 = ax2.plot(gamma, relative_error, 'b',label="relative error")
    ax2.set_ybound(lower=0)
    ax2.set_ylabel("relative error",fontsize=20)
    ax2.set_xbound(lower=0)
    l = l1+l2
    labs = [l.get_label() for l in l]
    ax.legend(l,labs,loc = 0,fontsize=20)
    plt.show()
    return best_gamma


#根据提供的数据直接画出图像(程序中未使用，仅供助教测试数据使用)
def plot_data(fime_name):
    df = pd.read_csv(fime_name)
    gamma = df["gamma"]
    sample_rate = df["sample_rate"]
    relative_error = df["relative_error"]


    ax = plt.subplot(111)
    l1 =  ax.plot(gamma,sample_rate,'g',label = "sample rate")
    ax.set_ybound(lower=0)
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel("sample rate")

    ax2 = ax.twinx()
    l2 =ax2.plot(gamma,relative_error,'b',label = "relative error")
    ax2.set_ybound(lower=0)
    ax2.set_ylabel("relative error")
    l = l1+l2
    labs = [l.get_label() for l in l]
    ax.legend(l,labs,loc = 0)
    plt.show()



#主函数
if __name__ == '__main__':
    #定义一些函数
    def fx(x,alpha,beta):
        return (x/beta)**(alpha-1) * np.exp(-x/beta)/(beta*gammaFunc(alpha))
    def hx(x,alpha,beta):
        return (x-alpha*beta)**2
    def gx(x,alpha,beta):# hx*fx
        return ((x-alpha*beta)**2)*fx(x, alpha, beta)
    def Tx(x,gamma):
        return np.exp(-x/gamma)/gamma
    
    #给定alpha,beta,gamma 和 N
    a = 2
    b = 3
    g = a*b
    N = int(5e5)
    I0 = a*(b**2)



#####################################################################################
    # 实验1：p(x)=f(x)的情况

    # # 程序输出
    px = fx
    m = MetroSample(px,alpha = a,beta = b,gamma = a*b)

    # # 抽样 并绘图
    # m.plot_hist(N)
    # plt.show()

    # 计算积分
    I1 = m.integrate(hx=hx,N=N)
    print("I={} ,I1 = {}".format(I0,I1))
    print("relative error = {}".format(np.abs(I1-I0)/I0))
    print("sample rate = {}".format(m.sample_rate))

    # 根据提供数据绘图---误差、效率、gamma取值的关系图
    # plot_data("data_error_rate_gamma_1.csv")
    # plt.show()

    # #绘制图像：误差和效率随gamma的变化
    # px=fx
    # df,best_gamma = find_best_gamma(px, hx, 2e4, a, b, proposal_gamma=a*b)
    # plot_sample_gamma(df,best_gamma,a,b)
    
    # 保存数据
    # df.to_csv("data_error_rate_gamma.csv")

    # #寻找最优gamma表达式
    # px=fx
    # size = 5
    # a_list = np.arange(3, 3+size, 1, int)
    # b_list = np.arange(3, 3+size, 1, int)
    # find_best_gamma_in_range(px,hx,int(1e4), a_list, b_list,proposal_gamma=a*b,method = 2)

############################################################################################

    # # 实验2: p(x)=hx*fx的情况

    # 程序输出
    px = gx         #gx = hx*fx
    hx = 1
    g = (a+2)*b
    m = MetroSample(px,alpha = a,beta = b,gamma = g)

    # #抽样并绘图
    # m.plot_hist(N)

    # 计算积分
    I2 = m.integrate(hx=1,N=N)

    print("\nI={} ,I2 = {}".format(I0,I2))
    print("relative error = {}".format(np.abs(I2-I0)/I0))
    print("sample rate = {}".format(m.sample_rate))


    # # 根据提供数据绘图---误差、效率、gamma取值的关系图
    # plt.figure()
    # plot_data("data_error_rate_gamma.csv")
    # plt.show()

    # #绘制图像：误差和效率随gamma的变化
    # px = gx
    # hx = 1
    # df,best_gamma = find_best_gamma(px,hx=1,N=2e4,alpha = a,beta = b,proposal_gamma=(a+2)*b)
    # plot_sample_gamma(df,best_gamma,a,b)

    # #保存数据
    # df.to_csv("data_error_rate_gamma_2.csv")

    # #寻找最优gamma表达式
    # px=gx
    # hx = 1
    # size = 5
    # a_list = np.arange(3, 3+size, 1, int)
    # b_list = np.arange(3, 3+size, 1, int)
    # find_best_gamma_in_range(px, hx, 1e4, a_list, b_list,proposal_gamma=(a+2)*b, method = 1)


