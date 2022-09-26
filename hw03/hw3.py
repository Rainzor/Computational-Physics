#Using UTF-8 Code Python
#Runze Wang
#PB20020480

from Schrage_16807 import*
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    seed = seed_time();
    r = Schrage16807(seed= seed)
    N = 5000
    #得到二维随机数，Schrage16807封装在Schrage_16807.py中
    a = r.rand(N,2)
    #print(a[:30])
    df = pd.DataFrame(columns=['theta','phi'])
    df['theta'] = np.arccos(1-a[:,0])
    df['phi'] = a[:,1]*2*np.pi
    x = np.sin(df['theta'])*np.cos(df['phi'])
    y = np.sin(df['theta'])*np.sin(df['phi'])
    z = np.cos(df['theta'])
    
    # #存储数据
    # df1 = pd.DataFrame({
    #     "x":x,
    #     "y":y,
    #     "z":z,
    # })
    # df1.to_csv("data.csv")

    fig = plt.figure() # 创建一个画布figure，然后在这个画布上加各种元素。
    ax = plt.axes(projection ="3d")#绘制三维图像
    ax.scatter3D(x, y, z,s = 0.5)
    plt.title("Evenly distributed on the upper hemisphere")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")   
    plt.show()
