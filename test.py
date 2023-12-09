from vorticity import *

s1 = Surf('KB', fac=0.1, vi=3) # 'T' or 'KB' or 'PP'
ts = t_plot = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60] # 作图时间点
dt = 1e-2 # 时间步长
nu = 1e-3 # 耗散系数
s1.sim(dt, nu, ts) # 生成对应参数的文件夹，并在给定时间片作图