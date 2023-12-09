from vorticity import *

s1 = Surf('KB', fac=0.4)
ts = t_plot = [0, 10, 20, 40, 45, 50, 55, 60, 70, 80] # 作图时间点
dt = 1e-2
nu = 1e-3
s1.sim(dt, nu, ts)