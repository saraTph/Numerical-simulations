# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 12:02:25 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

n_values = [4.6e9, 4.554e9]
len_sim = 54

en = np.ones((len(n_values),len_sim))
size = np.ones((len(n_values),len_sim))
pop = np.ones((len(n_values),len_sim))
delta_values = np.ones((len(n_values),len_sim))

location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error"

file = os.path.join(location, 'sim_4_4.6n0.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_4_4.554n0.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


#%%

wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

L = np.sqrt(hbar/(m*wr))
E = hbar*wr
sigma_z = 14e-6

#%% size after tof

colors = plt.get_cmap('Set2').colors
lw = 1.7

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(8,5))

shift = [0.12, 0.122, 0.158, 0.257, 0.257, 0.376]
#for i in range(np.size(Omega_values)):

Den = (en[0,:]-en[1,:])  #J
DenDn = Den/(n_values[0]-n_values[1])
n_avg = (n_values[0]+n_values[1])/2
DnDs = -(n_avg/sigma_z)

N = np.sqrt(2*np.pi) * sigma_z * n_avg
F = - m*wz**2*sigma_z - 1/(N) * DenDn*DnDs
F_t = F[::-1]
d = delta_values[0,:]
d_t = d[::-1]

t = np.linspace(0,9e-3,54)
dt = np.diff(t)  # Time step differences
dt = np.append(dt, dt[-1])

delta_v = np.cumsum(F_t / m * dt)
Ekin = 0.5*m*delta_v**2

Hz_per_um_conversion = 1 / (hbar * 2 * np.pi) * 1e-6  # Convert J/m to Hz/um

#axS.plot(t *1e3, Ekin/(hbar*2*np.pi), lw=lw, color=colors[0], zorder =1)
axS.plot(d_t, Ekin/(hbar*2*np.pi), lw=lw, color=colors[0], zorder =1)

 
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
#axS.set_xlabel(r'$t (ms)$', fontsize=14)
axS.set_ylabel(r'$E_{kin} (Hz)$',fontsize = 14)  

axS.grid()

#%%
#figS.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\FORT2 off\Ekin_t.png', dpi = 300)
