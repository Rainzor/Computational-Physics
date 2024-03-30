import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def twopart(n):
    #判断n是否为2的幂
    return n&(n-1) == 0


if __name__ == '__main__':
    func = lambda x, lamb: lamb*np.sin(np.pi*x)
    x0 = 0.5
    lower = -2
    upper = 2

    lambs = np.arange(lower,upper,1e-3,dtype=np.float64)   
    M = 5000
    N = 1000 
    
    #   #精确计算分叉点
    # M = 8000
    # N = 3000
    # lamb1 = np.arange(0.715,0.725,1e-6,dtype=np.float64)
    # lamb2 = np.arange(0.83,0.835,1e-6,dtype=np.float64)
    # lamb3 = np.arange(0.85,0.87,1e-6,dtype=np.float64)
    # lambs = np.concatenate((lamb1,lamb2,lamb3))

    size_list = []
    data_list = []
    d_alpha = []
    T_alpha = []    
    for lamb in lambs:
        flag = 0
        x = x0
        data = []
        for i in range(M):
            x = func(x,lamb)
        for j in range(N): 
            x = func(x,lamb)
            tmp = np.array(data)
            if (len(data)>0) and (np.any(np.abs(tmp-x)<1e-5))  :          
                continue
            else:
                data.append(x)  
                if (x<0.5005) & (x>0.4995) & (flag==0):
                    flag = 1
        #寻找x=0.5的纵向距离
        if (flag == 1) & (len(data)>1) & (len(data)<257) & twopart(len(data)):
            index = int(np.log2(len(data)))*2-1
            T_alpha.append(len(data))
            sort_s = sorted(data)
            d_alpha.append(sort_s[index]-sort_s[index-1])
            flag=0
            
        size_list.append(len(data))#记录每个lambda对应的数据点个数
        data_list.append(data)#记录每个lambda对应的数据点



    #拼接数据
    lamb_list = np.repeat(lambs, size_list)
    x_list = []
    for data in data_list:
        x_list += list(data)
    x_list = np.array(x_list)

    #保存数据
    # df = pd.DataFrame({'lambda':lamb_list, 'x':x_list})
    # df.to_csv('../data/data.csv', index=False)

    # 保存lambda与周期的数据
    # df = pd.DataFrame({'lambda':lambs,'size':size_list})
    # df.to_csv('../data/Feigenbaum_size.csv',index=False,float_format = '%.6f')

    # # 保存纵向距离数据
    # df2 = pd.DataFrame({"T":T_alpha,"d":d_alpha})
    # df2.to_csv('../data/Feigenbaum_alpha.csv',index=False,float_format = '%.6f')
    
    # ##绘图
    ax = plt.subplot(111)
    ax.scatter(lamb_list, x_list, s=0.005, c='k')
    # ax.scatter(lamb_list, x_list, s=0.005,c=lamb_list, cmap='RdYlBu')
    ax.set_xlabel(r'$\lambda$', fontsize=20)
    ax.set_ylabel(r'$x$', fontsize=20)
    ax.set_title('Bifurcation diagram', fontsize=20)
    ax.set_xlim(lower,upper)
    #ax.axhline(y=0.5,color='r',linestyle='--')
    plt.show()




