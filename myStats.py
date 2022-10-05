from abc import ABC, abstractmethod

from Schrage_16807 import*
import scipy.special as sc
from math import factorial

class rv_generi(ABC):
    """Class which encapsulates common functionality between rv_discrete
    and rv_continuous.

    Parameters
    ----------  
    seed :
        random object using 16807 Method
    """

    def __init__(self) -> None:
        super().__init__()
        self.seed = Schrage16807(seed = seed_time())
        self._stats = {"mid":0,"mean":0,"var":0}
        
    
    def rvs(self, size=1):
        """Random variates using direct or selection sampling method 

            Parameters
            -----------
            size(int): the size of random variates as np.array format 

            Return
            ------
            rvs : ndarray or scalar
            Random variates of given `size`.

            得到size大小的随机数
        """
        if size ==1:
            return self.ppf(self.seed.rand())
        else:
            xi = np.array([self.ppf(self.seed.rand()) for i in range(size)])
            return xi
    
    @abstractmethod
    def cdf(self, x):
        """Cumulative distribution function of the given RV.
        
        累积分布函数

        Parameters
        ----------
        x : array_like
            quantiles
			
		
        Returns
        -------
        cdf : ndarray
            Cumulative distribution function evaluated at `x`
        """

    @abstractmethod
    def ppf(self, p):
        """Percent point function (inverse of `cdf`) at q of the given RV.

        分位点 (累积分布函数的反函数).

        Parameters
        ----------
        p : array_like
            lower tail probability
        

        Returns
        -------
        x : array_like
            quantile corresponding to the lower tail probability q.
        """

    @property
    def median(self):
        """median of the distribution.
        """
        return self._stats["mid"]


    @property
    def mean(self):
        """Mean of the distribution.
        """
        return self._stats["mean"]

    @property
    def var(self):
        """Variance of the distribution.
        """
        return self._stats["var"]

    @property
    def std(self):
        """Standard deviation of the distribution.
        """
        return np.sqrt(self._stats["var"])
    
    @property
    def stats(self):
        """
        Returns
        -------
        mean,var

        """
        return self.mean,self.var

    #二分查找法
    def _binary_search(self, x, low, high, data):
        size = high-low+1
        while low <= high:
            mid = (low+high)//2
            if data[mid] == x:
                return mid
            elif data[mid] > x:
                high = mid-1
            else:
                low = mid+1
        if low < size:
            return low
        else:
            return 0
        

class rv_discrete(rv_generi):
    """ A generic discrete random variable class meant for subclassing.

    `rv_discrete` 是构造特定分布类的基类 以及离散随机变量的实例

    """

    def __init__(self, a=0, b=np.inf):
        super().__init__()
        self.a = a
        self.b = b  

    def pmf(self, k):
        """Probability mass function at k of the given RV.

        散点处概率分布函数

        Parameters
        ----------
        k : array_like
            Quantiles.
            
        Returns
        -------
        pmf : array_like
            Probability mass function evaluated at k

        """
    def pdf(self, k):
        return self.pmf(k)



class rv_continuous(rv_generi):
    """A generic continuous random variable class meant for subclassing.

    `rv_continuous` 是构造特定分布类的基类以及连续随机变量的实例

    Parameters
    ----------
    a : float, optional
        Lower bound of the support of the distribution, default is minus
        infinity.
    b : float, optional
        Upper bound of the support of the distribution, default is plus
        infinity.
    seed :
        random object using 16807 Method

    """
    def __init__(self,a=-np.inf, b = np.inf) -> None:
        super().__init__()
        self.a = a
        self.b = b

    @abstractmethod
    def pdf(self, x):
       """Probability density function at x of the given RV.

        概率密度函数 

        Parameters
        ----------
        x : array_like
            quantiles
        
        Returns
        -------
        pdf : ndarray
            Probability density function evaluated at x
        """



class norm(rv_continuous):
    r"""A normal continuous random variable.

        正态分布

        The location (``loc``) keyword specifies the mean.
        The scale (``scale``) keyword specifies the standard deviation. 
       
    """
    def __init__(self, loc=0, scale = 1 ) -> None:
        super().__init__()
        self.seed = Schrage16807(seed_time())
        self._stats={}
        self._stats["mid"] = loc
        self._stats["mean"] = loc
        self._stats["var"] = scale**2 


    def rvs(self, size=1):
        #Box-Muller法抽样
        i = 0
        output = np.zeros(size)
        while i < size:
            u = self.seed.rand()*2-1
            v = self.seed.rand()*2-1
            r = np.sqrt(u**2+v**2)
            if r>1:
                continue
            else:
                output[i] = (u/r)*np.sqrt(-4*np.log(r))
                i = i+1
        if size == 1:
            return output[0]
        else:
            return output

    def pdf(self,x):
        return np.exp(-x**2/(2*self.var**2))/np.sqrt(2*np.pi*self.var)
    
    def ndtr(self,x):
        #误差公式
        return sc.ndtr(x)

    def cdf(self, x):  
        return self.ndtr((x-self.mean)/self.std)

    def ndtri(self,p):
        return sc.ndtri(p)

    def ppf(self, p):
        #直接使用了库中的反函数erfinv
        if 0 < p and p < 1:
            return self.ndtri(p)*self.std+self.mean
        elif p == 0:
            return -np.inf
        elif p == 1:
            return np.inf
        else:
            raise ValueError("p is out of range")
    
class expon(rv_continuous):
    """An exponential continuous random variable.

        指数分布

    """

    def __init__(self, lambdaa = 1) -> None:
        super().__init__(a=0)
        self._lambda = lambdaa   
        self._stats={}
        self._stats["mid"] = np.log(2)/lambdaa
        self._stats["mean"] = 1/lambdaa
        self._stats["var"] = 1/lambdaa**2   

    @property
    def lambdaa(self):
        return self.lambdaa

    def pdf(self, x):
        return np.where(x>=0,self._lambda*np.exp(-self._lambda*x),0 )


    def cdf(self,x):
        if x >= 0:
            return 1-np.exp(-self._lambda*x)
        else: return 0
    
    def ppf(self,p):
        if 0 < p and p < 1:
            return -np.log(1-p)/self._lambda
        elif p == 0:
            return 0
        elif p == 1:
            return np.inf
        else:
            print("ERROR!")
    
   


class poisson(rv_discrete):
    r"""A Poisson discrete random variable.
            泊松分布

    Parameters
        ----------
        lambdaa : float>0
            _description_, by default 1
    """

    def __init__(self, lambdaa=1):
        if lambdaa>0:
            super().__init__()
            self._lambda = lambdaa
            self.seed = Schrage16807(seed_time())
            self._stats={}
            self._stats["mid"] = lambdaa+1/3-0.02/lambdaa
            self._stats["mean"] = lambdaa
            self._stats["var"] = lambdaa
        else:
            raise ValueError("lambda need to be larger than 0")
    
    def rvs(self, size=1):
        i = 0
        output = np.zeros(size)
        if self._lambda > 0.01:
            while i<size:
                L = np.exp(-self._lambda)
                k = 0
                p = 1
                while True:
                    k=k+1
                    u = self.seed.rand()
                    p = p*u
                    if(p<=L):
                        break                    
                output[i] = k-1
                i = i+1
        else:
            while i<size:
                k = 0
                p = np.exp(-self._lambda)
                s = p
                u = self.seed.rand()
                while True:
                    k = k+1
                    p = p*self._lambda/k
                    s = s+p
                    if u<=s:
                        break
                output[i] = k
                i = i+1

        if size ==1:
            return output[0]
        else:
            return output


    def pmf(self, k):
        if isinstance(k,int):
            if k>=0:
                return self._lambda**k*np.exp(-self._lambda)/factorial(k)
            else:
                return 0
        else:
            k = np.floor(k)
            output = []
            for kk in k:
                if kk>=0:
                    t = self._lambda**kk*np.exp(-self._lambda)/factorial(kk)
                    output.append(t)
                else:
                    output.append(0)
                
            return np.array(output)
        

    def cdf(self, k):
        if (k>=0):
            l = [(self._lambda**i)/np.math.factorial(i) for i in range(np.ceil(k))]
            return np.exp(-self._lambda)*np.sum(l)
        else:
            raise ValueError("k need to be at least 0")

    def ppf(self, p):
        c = 0
        if self._lambda <15:
            for k in range(30):
                c = c + self.pmf(k)
                if c>p:
                    return k
            return np.ceil(2*self._lambda)
        else:#当lambda足够大时，泊松分布趋近于正态分布
            return np.floor(norm.ndtri(p)*self._lambda+self._lambda) #反函数误差函数erfinv
        
    @property
    def lambdaa(self):
        return self._lambda


class uniform(rv_continuous):
    r"""A uniform continuous random variable.

    In the standard form, the distribution is uniform on ``[0, 1]``. Using
    the parameters ``loc`` and ``scale``, one obtains the uniform distribution
    on ``[loc, loc + scale]``.
    """
    def __init__(self, loc = 0, scale=1) -> None:
        super().__init__(loc, loc+scale)
        self.a = loc
        self.b = loc+scale
        self._stats = {}
        self._stats["mid"] = (self.a + self.b)/2
        self._stats["mean"] =(self.a + self.b)/2
        self._stats["var"] = scale**2/12

    def pdf(self, x):
        if x>super().b or x<super().a:
            return 0
        else:
            return super().pdf(x)

    def cdf(self, x):
        if x<self.a:
            return 0
        elif x>self.b:
            return 1
        else:
            return (x-self.a)/(self._b-self.a)

    def ppf(self, p):
        if 0 < p and p < 1:
            return (self.b-self.a)*p+self.a
        elif p == 0:
            return self.a
        elif p == 1:
            return self.b
        else:
            print("ERROR!")

    @property
    def scale(self):
        return self.a,self.b


class bernulli(rv_discrete):
    """A Bernoulli discrete random variable.

    `bernoulli` takes :math:`p` as shape parameter,
    where :math:`p` is the probability of a single success
    and :math:`1-p` is the probability of a single failure.
    """
    def __init__(self, p=1/2):
        super().__init__(a=0,b=1)
        self._p = p
        self._stats["mid"] = None
        self._stats["mean"] = p
        self._stats["var"] = p*(1-p)

    def pmf(self, k):
        if k==0:
            return 1-self._p
        elif k==1:
            return self._p
        else:
            return 0

    def cdf(self, k):
        if 0<=k and k<1:
            return 1-self._p
        elif k>=1:
            return 1
        else:
            return 0

    def ppf(self,p):
        if p<1-self._p:
            return 0
        else: return 1
    

class binom(rv_discrete):
    """A binomial discrete random variable.

    `binom` takes :math:`n`>0 and :math:`p` as shape parameters,
    where :math:`p` is the probability of a single success
    and :math:`1-p` is the probability of a single failure.
    """
    def __init__(self,n=1,p=1/2):
        super().__init__(a=0,b=n)
        self._n = n
        self._p = p
        self.seed = Schrage16807(seed_time())
        self._stats["mid"] = np.floor(n*p) if abs(np.floor(n*p)-n*p)<min(np.log(2),max(p,1-p)) else np.ceil(n*p)
        self._stats["mean"] =n*p
        self._stats["var"] = n*p*(1-p)
        self.__data = []
        if n<=25:
            i=1
            self.__data.append(self.pmf(0))
            while i<=n:
                t = self.__data[-1]+self.pmf(i)
                self.__data.append(t)
                i = i+1

    

            

    def rvs(self, size=1):
        if size ==1:
            xi = self.seed.rand(size)
            x = np.where(xi<1-self._p,0,1)
            return sum(x)
        else:
            xi = self.seed.rand(d0 = self._n ,d1=  size)
            # print("xi",xi)
            x = np.where(xi < 1-self._p,0,1)
            # print("x:",x)
            return np.sum(x ,axis=0)

    def comb(self,n,m):
        return np.math.factorial(n)//(np.math.factorial(n-m)*np.math.factorial(m))
    def pmf(self, k):
        return self.comb(self._n,k)*(self._p**k)*((1-self._p)**(self._n-k))
    
    def cdf(self, k):
        if k>=self._n:
            return 1
        elif k<0:
            return 0
        else:
            x = int(k)
            c = [self.pmf(i) for i in range(x+1)]
            return sum(c)

    def ppf(self, p):
        if len(self.__data) > 0 :
            return super()._binary_search(p, 0, self._n+1, self.__data)
        else:
            return norm.ndtri(p)*self.std+self.mean

class Cos(rv_continuous):
    """Cos function from -pi/2 to pi/2

    its cumulative function is sin(x)

    its percent point function is arcsin(x)

    """
    def __init__(self) -> None:
        super().__init__(0,np.pi/2)
        self.seed = Schrage16807(seed_time())
        self._stats["mid"] = np.pi/6
        self._stats["mean"] = np.pi/2-1
        self._stats["var"] =-3+np.pi
            

    def pdf(self, x):
        return np.cos(x)
    
    def cdf(self, x):
        return np.sin(x)
    
    def ppf(self, p):
        return np.arcsin(p)
    
