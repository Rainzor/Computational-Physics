from Schrage_16807 import*
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


if __name__ == "__main__":
    seed = seed();
    r = Schrage_16807(seed= seed)
    a = [r.rand(dim =2) for item in range(5000)]
    a = np.asarray(a)
    # a = np.random.randn(10000,2)
    df = pd.DataFrame(columns=['theta','phi'])
    df['theta'] = np.arccos(a[:,0]*2-1)
    df['phi'] = a[:,1]*2*np.pi
    print(df.head())
    x = np.sin(df['theta'])*np.cos(df['phi'])
    y = np.sin(df['theta'])*np.sin(df['phi'])
    z = np.cos(df['theta'])
    
    fig = plt.figure() # 创建一个画布figure，然后在这个画布上加各种元素。
    ax = plt.axes(projection ="3d")
    ax.scatter3D(x, y, z,s = 0.5)
    plt.show()
