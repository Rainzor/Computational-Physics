from Schrage_16807 import*
import numpy as np
import matplotlib.pyplot as plt

#概率分布函数
def porb_x(x):
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
        return M1*y<=porb_x(x)
    else:
        return M2*y<=porb_x(x)

if __name__=='__main__':
    seed = seed_time()
    r = Schrage_16807(seed)
    N = 100000
    #生成二维均匀随机变量
    a = [r.rand(dim =2) for item in range(N)]
    #a = np.random.rand(100000,2)  #内置随机数器，可以与自己的进行比较
    xi = np.asarray(a)
    xi[:,0] = xi[:,0]*2-1

    #python 内置迭代器，舍选法选择合适的数
    result =[x[0] for x in xi if rejection_method(x)]
    result = np.array(result)

    #统计区间内，单位长度（1/num）点个数,方便绘图查看
    num = 100
    intervel = list(np.linspace(start=-1,stop=1,num=num+1))
    count = []
    for i in range(num):
        index = (result>=intervel[i]) & (result<(intervel[i+1]))
        count.append(result[index].size)
    #z柱状图
    plt.bar(range(num),count)
    x_label = ['-1','-0.5','0','0.5','1']
    plt.xticks(np.linspace(0,num,num=5),x_label)
    plt.show()


    

