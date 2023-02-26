主程序在“DLA.py”，“DBM.py”，“DBM_FAST.py”中
DLA.py: 含有 Growth_Model抽象类，DLA类继承了Growth_Model类
DBM.py: 调用了DLA.py，DBM类继承了Growth_Model类
DBM_FAST.py: 调用了DBM.py，DBM_FAST继承了DBM类

主要用到numpy、time、pandas、matplotlib包

DLA.gif是DLA方法生成的动画
DBM_FAST.gif是快速DBM法生成的动画

压缩包中提供了测试集数据“DLA.csv”,“DBM.csv”,“DBM_FAST.csv”
其中含有生长点的坐标

代码中没有直接调用数据，而是采取“时间-种子”方式，现场形成随机数，
可以直接调用程序即可画图

本实验对16807Schrage产生器进行单独封装，其代码在“Schrage_16807.py”中。
内部包含“seed_time”随机数产生种子函数，“Schage16807”类

