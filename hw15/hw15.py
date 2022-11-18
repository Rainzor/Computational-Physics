from Schrage_16807 import*
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
        The parameters of pdf and T :beta, hamiltonian

    Methods
    -------
    Tx : function
        The proposal distribution

    Ax :
        The acceptance probability

    
    _generate_data : (private)
        Generate N data points from the distribution
    
    sample : (public)
        Sample N points from the distribution

    plot_hist :
        Plot the histogram of the samples
    
    integrate :
        Integrate the hx function from 0 to infinity
    
    """

    def __init__(self, pdf, **kwargs):
        if len(kwargs) == 0:
            raise ValueError("No parameters")
        self.pdf = pdf
        self.hamiltonian = kwargs.get('hamiltonian', None)
        self.beta = kwargs.get('beta', 1)
        self.sigma = 1/np.sqrt(4*self.beta)
        self.r = Schrage16807(seed=seed_time())
        self.burn = None
        self.sample_rate = 1
        self.size = 0
        self.markov_chain = None
        self.C = 1  # 归一化常数

    #待抽样的概率分布（可能未归一化）
    def px(self, x):
        return self.pdf(x[0], x[1], self.beta)

    def _inv_px(self, x):
        return 


    #高斯分布
    def gaussian(self, x,y, mu, sigma):
        return np.exp(-((x-mu)**2+(y-mu)**2)/(2*sigma**2))/(2*np.pi*sigma**2)

    # 建议分布
    def Tx(self, x):
        return 1/2*(self.gaussian(x[0],x[1], np.sqrt(2), self.sigma)+self.gaussian(x[0], x[1], -np.sqrt(2), self.sigma))

    # 抽样建议分布Tx
    def _sample_Tx(self):
        """Sample N points(x,y) from the proposal distribution"""
        if self.r.rand() < 0.5:
            return self.r.normal(size=2,loc=np.sqrt(2),scale= self.sigma)
        else:
            return self.r.normal(loc=-np.sqrt(2),scale= self.sigma, size=2)

    # 接受分布
    def Ax(self, x, xt):
        return self.px(xt)*self.Tx(x)/(self.Tx(xt)*self.px(x))


    #核心：Metropolis抽样方法,内部接口
    def _generate_data(self, N, x0=None):
        """Generate N data points from the distribution

        Parameters
        ----------
        N : int
            The number of data points to be generated

        Returns
        -------
        data : array (N,2)
            The data points from the distribution
        """
        if x0 is None:
            if self.hamiltonian is None:
                x = np.zeros((N+1, 2))
                x0 = self.r.normal(size=2)
                x[0] = x0
            else:
                m = 10000  # the number of points to be discarded
                x = np.zeros((N+1,2))
                b = np.zeros((m,2))
                b[0] = np.array([5,5])
                n = 1
                # 丢弃前m个点,热化过程,按照走一步试一步探路的方法热化（为了体现marcov链过程）
                step = 0.1
                
                while n < m:
                    xt = self.r.rand(2)*2-1
                    xt = xt*step+b[n-1]

                    delta_E = self.hamiltonian(xt[0],xt[1])-self.hamiltonian(b[n-1][0],b[n-1][1])
                    r1 = np.exp(-delta_E*self.beta)
                    #r1 = min(1, r1) 

                    xi = self.r.rand()
                    b[n] = xt if xi < r1 else b[n-1] 
                    n += 1

                x[0] = b[-1]
                self.burn = b


        else:
            x = np.zeros((N+1,2))
            x[0] = x0
        
  
        n = 1
        s = 0  # the number of points be selected
        while n < N:
            x_temp = self._sample_Tx()
            r1 = self.Ax(x[n-1], x_temp)
            xi = self.r.rand()
            # 舍选过程
            if xi < min(1, r1):
                x[n] = x_temp
                s += 1
            else:
                x[n] = x[n-1]
            n += 1
        # self.sample_rate = self.sample_rate*(self.size/(self.size+N)) + s/n*(N/(self.size+N))
        self.sample_rate = s/n

        return x[1:]

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
            self.markov_chain = self._generate_data(N)
            self.size = N
            return self.markov_chain
        elif(self.size < N):
            self.markov_chain = np.concatenate(
                (self.markov_chain, self._generate_data(N-self.size, x0=self.markov_chain[-1])),axis=0)
            self.size = N
            return self.markov_chain
        else:
            return self.markov_chain[:N]

    # 数值积分计算 hx是被积函数
    def integrate(self, hx, N):
        """Integrate the function hx"""
        N = int(N)
        x = self.sample(N)
        #重要抽样方法（实验1）
        if(hx != 1):
            I1 = np.mean(hx(x[:,0],x[:,1]))
            return I1

        # 比值法（实验2）
        pass
    
    # 绘制marrkov链
    def plot_mcmc(self,N):
        """Plot the markov chain"""
        N = int(N)
        if(N<500):
            a = 0.3
        else:
            a = 0.1
        x = self.sample(N)
        plt.figure();
        ax = plt.subplot(111)

        ax.plot(x[:,0], x[:,1], 'o', markersize=1,color='red')
        ax.plot(x[:,0], x[:,1], '--', linewidth=0.5, color='blue', alpha=a)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(r'Markov Chain ($\beta$={})'.format(self.beta))

    def plot_burn(self):
        """Plot the burn-in points"""
        plt.figure()
        x = self.burn
        ax = plt.subplot(111)
        # ax.plot(x[:, 0], x[:, 1], 'o', markersize=1, color='red')
        ax.plot(x[:, 0], x[:, 1], '-', linewidth=0.5, color='blue', alpha=0.8)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(r'Markov Chain Burn-in Points ($\beta$={})'.format(self.beta))

if __name__ == '__main__':
    H = lambda x,y: -2*(x**2+y**2) + (x**4+y**4)/2 + ((x-y)**4)/2
    F =lambda x,y,beta: np.exp(-beta*H(x,y))
    h1 = lambda x,y: x**2
    h2 = lambda x,y: y**2
    h3 = lambda x,y: x**2+y**2
    
    N = 100000
    
    # #绘制图像Markov链
    # mcmc = MetroSample(F, beta=5, hamiltonian=H)
    # mcmc.plot_mcmc(200)
    # mcmc.plot_burn()
    # plt.show()

    for b in [0.2,1,5]:
        mcmc = MetroSample(F, beta=b)
        mcmc.sample(N)
        x2 = mcmc.integrate(h1,N)
        y2 = mcmc.integrate(h2,N)
        x2_y2 = mcmc.integrate(h3,N)
        print('beta = {}\t <x^2> = {}, <y^2> = {}, <x^2+y^2> = {}\n'.format(b,x2,y2,x2_y2))

        # ##输出数据
        # n = N if N < 10000 else 10000
        # df = pd.DataFrame({"x":mcmc.markov_chain[:n,0],"y":mcmc.markov_chain[:n,1]})
        # df.to_csv("data_beta_{}.csv".format(b),index=False)
        


