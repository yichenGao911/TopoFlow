from vorticity import *

s1 = Surf('Torus', fac=0.2, vi=1, pert=0.01, num=128) # 'Torus' or 'KB' or 'PP'
# ts = [0, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60] # 作图时间点
# ts += [65, 70, 75, 80, 90, 100]
# ts += [120, 150, 200]
ts = np.arange(0, 300, 5)
dt = 1e-2 # 时间步长
nu = 2e-4 # 耗散系数
s1.sim(dt, nu, ts) # 生成对应参数的文件夹，并在给定时间片作图