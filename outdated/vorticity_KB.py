import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2
import numpy as np
import os

num = 64
Kmax = num//2
Jmax = num
jj = (0+1j)
ks = np.append(np.linspace(0, Kmax, Kmax+1), -np.linspace(Kmax-1, 1, Kmax-1))
js = np.append(np.linspace(0, Jmax, Jmax+1), -np.linspace(Jmax-1, 1, Jmax-1))
kjs = np.array(np.meshgrid(js, ks))
laplace = -(kjs[1]**2+kjs[0]**2)
laplace_0 = laplace.copy()
laplace_0[np.where(laplace_0==0)] = 1e9 # 消除0值
laplace_rev = 1/laplace_0

def rk4(x, dt, tend, arg):
    q1 = dt*tend(x, arg)
    q2 = dt*tend(x+q1/2, arg)
    q3 = dt*tend(x+q2/2, arg)
    q4 = dt*tend(x+q3, arg)
    return x+(q1+2*q2+2*q3+q4)/6

def ic(num, fac, intensity, pert):
    n_half = num//2
    n1 = int(n_half*(1-fac))
    n2 = n_half - n1

    vor_lat = np.concatenate([np.zeros(n1), -np.ones(n2), np.ones(n2), np.zeros(n1)])
    vor_lat *= intensity
    vor = np.expand_dims(vor_lat, 1).repeat(num, axis=1)

    vor[n_half][n_half] += pert

    return vor

def plot_it(t, lon, lat, vor):
    fname = "vorticity_"+str(int(t))+"_seconds_later.png"
    fig, ax = plt.subplots(figsize=(6,4))

    crange = np.linspace(-1.1, 1.1, 12)
    ax.set_aspect('equal', adjustable='box')
    ax.set_title(fname)
    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    cf = ax.contourf(lon, lat, vor, crange, cmap='bwr')
    cbar = fig.colorbar(cf)
    cbar.ax.set_ylabel("vorticity")

    plt.tight_layout()
    plt.savefig(fname, format='png')


def tend(vor, nu=0.1):
    vor_inv = -vor[::-1, :]
    vor_dbx = np.concatenate([vor_inv, vor], axis=1)
    f = fft2(vor_dbx)
    psi = f*laplace_rev
    fu = -psi*jj*kjs[1]
    fv = psi*jj*kjs[0]
    u = ifft2(fu).real
    # u -= u[0][0]*np.ones_like(u) # 增加一个匀速流场消除边界流
    v = ifft2(fv).real
    pvpx = ifft2(f*jj*kjs[0]).real
    pvpy = ifft2(f*jj*kjs[1]).real
    Lv = ifft2(f*laplace).real

    # print(np.max(pvpx), np.max(pvpy), np.max(Lv))
    tendency = -u*pvpx - v*pvpy + nu*Lv
    tendency = tendency[:, num:]

    return tendency

fac_ini = .2 # 非零初始涡度空间占比
vor_ini = 1 # 初始涡度
pert_ini = .01 # 初始扰动
t_plot = [0, 10, 20, 40, 45, 50, 55, 60, 70, 80] # 设置画图的时间点
dt = 1e-2 # 时间步长
nu = 1e-3
if ('KB' not in os.listdir()):
    os.mkdir('KB')
os.chdir('KB')
dirname = "fac="+str(fac_ini)+"_dt="+str(dt)+"_nu="+str(nu)+"_pert="+str(pert_ini)
if (dirname not in os.listdir()):
    os.mkdir(dirname)
os.chdir(dirname)
# dirname = ''
N = int(max(t_plot)/dt)

lon = np.delete(np.linspace(0, 2*np.pi, num+1), -1)
lat = np.delete(np.linspace(0, 2*np.pi, num+1), -1)

vor = ic(num, fac_ini, vor_ini, pert_ini)
# u, test=tend(vor)
# plt.contourf(lon, lat, test, cmap="bwr")
# plt.colorbar()

if 0 in t_plot:
    plot_it(0, lon, lat, vor)

for i in range(N):
    vor = rk4(vor, dt, tend, nu)
    t = dt*(i+1)
    # print(t, t_plot)
    if t in t_plot:
        plot_it(t, lon, lat, vor)

os.chdir('..')
os.chdir('..')
plt.show()


