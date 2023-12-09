import matplotlib.pyplot as plt
from numpy.fft import fft2, ifft2
import numpy as np
import os

class Surf:
    def __init__(self, kind, fac = 0.5, vi = 1, pert = 0.01, num = 128):
        self.kind = kind # 种类
        self.fac = fac # 初始非零涡度区占比
        self.vi = vi # 初始场强度
        self.pert = pert # 扰动大小
        self.num = num # 分辨率（格点数）
        self.vor = self.__get_initial_vor()
        self.laplace, self.laplace_rev = self.__get_spectral_util()
        self.lon = np.delete(np.linspace(0, 2*np.pi, num+1), -1)
        self.lat = np.delete(np.linspace(0, 2*np.pi, num+1), -1)

    def __get_spectral_util(self):
        Kmax = self.num//2
        Jmax = self.num//2
        if (self.kind == 'KB'):
            Jmax *= 2
        if (self.kind == 'PP'):
            Jmax *= 2
            Kmax *= 2
        self.jj = (0+1j)
        ks = np.append(np.linspace(0, Kmax, Kmax+1), -np.linspace(Kmax-1, 1, Kmax-1))
        js = np.append(np.linspace(0, Jmax, Jmax+1), -np.linspace(Jmax-1, 1, Jmax-1))
        self.kjs = np.array(np.meshgrid(js, ks))
        laplace = -(self.kjs[1]**2+self.kjs[0]**2)
        laplace_0 = laplace.copy()
        laplace_0[np.where(laplace_0==0)] = 1e9 # 消除0值
        laplace_rev = 1/laplace_0
        return laplace, laplace_rev

    def __get_initial_vor(self):
        n_half = self.num//2
        n1 = int(n_half*(1-self.fac))
        n2 = n_half - n1

        vor_lat = np.concatenate([np.zeros(n1), -np.ones(n2), np.ones(n2), np.zeros(n1)])
        vor_lat *= self.vi
        vor = np.expand_dims(vor_lat, 1).repeat(self.num, axis=1)

        vor[n_half][n_half] += self.pert
        return vor
    
    def rk4(self, dt, nu):
        q1 = dt*self.tend(self.vor, nu)
        q2 = dt*self.tend(self.vor+q1/2, nu)
        q3 = dt*self.tend(self.vor+q2/2, nu)
        q4 = dt*self.tend(self.vor+q3, nu)
        self.vor += (q1+2*q2+2*q3+q4)/6

    def tend(self, vor, nu):
        if (self.kind == 'KB'):
            vor_temp = np.concatenate([-vor[::-1, :], vor], axis=1)
        elif (self.kind == 'PP'):
            vor_invy = -vor[::-1, :]
            vor_invx = -vor[:, ::-1]
            vor_invxy = vor[::-1, ::-1]
            vor_dbx = np.concatenate([vor_invy, vor], axis=1)
            vor_dbx_rev = np.concatenate([vor_invxy, vor_invx], axis=1)
            vor_temp = np.concatenate([vor_dbx_rev, vor_dbx], axis=0) 
        elif (self.kind == 'T'):
            vor_temp = vor
        
        f = fft2(vor_temp)
        psi = f*self.laplace_rev
        fu = -psi*self.jj*self.kjs[1]
        fv = psi*self.jj*self.kjs[0]
        u = ifft2(fu).real
        v = ifft2(fv).real
        pvpx = ifft2(f*self.jj*self.kjs[0]).real
        pvpy = ifft2(f*self.jj*self.kjs[1]).real
        Lv = ifft2(f*self.laplace).real

        # print(np.max(pvpx), np.max(pvpy), np.max(Lv))
        tendency = -u*pvpx - v*pvpy + nu*Lv
        if (self.kind == 'KB'):
            tendency = tendency[:, self.num:]
        elif (self.kind == 'PP'):
            tendency = tendency[self.num:, self.num:]

        return tendency
    
    def goto_folder(self, dt, nu):
        dirname = "fac="+str(self.fac)+"_dt="+str(dt)+"_nu="+str(nu)+"_pert="+str(self.pert)+"_vi="+str(self.vi)
        if (dirname not in os.listdir()):
            os.mkdir(dirname)
        os.chdir(dirname)

    def plot_it(self, t):
        fname = "vorticity_"+str(int(t))+"_seconds_later.png"
        fig, ax = plt.subplots(figsize=(6,4))

        crange = np.linspace(-self.vi, self.vi, 10)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(fname)
        ax.set_xlabel('longitude')
        ax.set_ylabel('latitude')
        cf = ax.contourf(self.lon, self.lat, self.vor, crange, cmap='bwr')
        cbar = fig.colorbar(cf)
        cbar.ax.set_ylabel("vorticity")

        plt.tight_layout()
        plt.savefig(fname, format='png')

    def sim(self, dt, nu, ts, show=False): 
        if (self.kind not in os.listdir()):
            os.mkdir(self.kind)
        os.chdir(self.kind)

        self.goto_folder(dt, nu)
        N = int(max(ts)/dt)
        
        if 0 in ts:
            self.plot_it(0)

        for i in range(N):
            self.rk4(dt, nu)
            t = dt*(i+1)
            # print(t, t_plot)
            if t in ts:
                self.plot_it(t)

        os.chdir('..')
        os.chdir('..')
        if (show):
            plt.show()