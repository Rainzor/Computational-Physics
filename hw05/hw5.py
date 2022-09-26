#Using UTF-8 Code Python
#Runze Wang
#PB20020480

from Schrage_16807 import*
import matplotlib.pyplot as plt
from matplotlib import cm


if __name__=='__main__':
    s = seed_time()
    r = Schrage16807(seed=s)
    N = 3000

    #二维均匀随机数
    a = r.rand(N,2)
    a = a*2-1

    #选择圆内的点
    choose_in_circle = lambda x:((x[0]**2+x[1]**2<1))
    result = np.array([x for x in a if choose_in_circle(x)])
    u = result[:,0]
    v = result[:,1]
    r2 = u**2+v**2
    x = 2*u*np.sqrt(1-r2)
    y = 2*v*np.sqrt(1-r2)
    z = 1-2*r2

    #数据集
    # df = pd.DataFrame({"x":x,"y":y,"z":z})
    # df.to_csv("data_05.csv")

    #绘图
    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.scatter(x,y,z,s=0.5)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("3D spherical surface")

    ax = fig.add_subplot(1, 2, 2)
    ax.scatter(x,y,marker='x',s=30)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("2D projection on X-Y plane")

    plt.axis('equal')
    plt.show()
    

