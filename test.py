from math import erf
import numpy as np
import matplotlib.pyplot as plt
from sympy import erfinv

from Schrage_16807 import Schrage16807
from Schrage_16807 import seed_time
#make random data of size, using 16807 method to generate random numbers and seed need to using time
def make_data(size):
    schrage = Schrage16807(seed_time())
    data = np.zeros(size)
    for i in range(size):
        data[i] = schrage.rand()
    return data


#a class make the gaussian distribution of the random data
class norm(object):
    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma
    def __call__(self, x):
        return np.exp(-(x - self.mu)**2/(2*self.sigma**2))
    def rvs(self, size):
        return np.random.randn(size)*self.sigma + self.mu
    def pdf(self, x):
        return np.exp(-(x - self.mu)**2/(2*self.sigma**2))
    def ppf(self,x):
        return self.mu + self.sigma*np.sqrt(2)*erfinv(2*x - 1)
    def cdf(self,x):
        return (1 + np.sign(x - self.mu)*erf(np.abs(x - self.mu)/np.sqrt(2)/self.sigma))/2