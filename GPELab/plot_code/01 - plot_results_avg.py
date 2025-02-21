# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 17:46:41 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

location = r"C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\Export"
file = os.path.join(location, '4w_data.txt')
en_4w, size_4w, pop_4w, delta_values = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

location = r"C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\Export"
file = os.path.join(location, '8w_data.txt')
en_8w, size_8w, pop_8w, delta_values = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

pop = np.vstack((pop_4w, pop_8w))
en = np.vstack((en_4w, en_8w))
size = np.vstack((size_4w, size_8w))

data_exp = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\data_MF.mat') 
size_exp = data_exp.get('valY').squeeze()
scan_exp = data_exp.get('valX').squeeze()
#%%

wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz
Omega_values = [4*wr, 8*wr]

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

a_bohr = 0.52917721e-10
a11 = (85*a_bohr)
a22 = (33.4*a_bohr)
a12 = (-53.0*a_bohr)
g11 = (4*np.pi*hbar**2*a11*a_bohr)/m
g12 = (4*np.pi*hbar**2*a12*a_bohr)/m 
g22 = (4*np.pi*hbar**2*a22*a_bohr)/m
  
V = 1/(0.05e-6**3)

#%%
colors = plt.get_cmap('Set2').colors
lw = 2

#labels = [r'$\Omega = 4\omega_r$', r'$\Omega = 8\omega_r$', r'$\Omega \gg \omega_r$']
labels = [r'$\Omega$ = 676 Hz', r'$\Omega$ = 1.35 kHz', r'$\Omega$ = 3.38 kHz']

# Plotting the expressions
fig, ax = plt.subplots(1,2,constrained_layout=True, figsize=(9,4))


for i in range(np.size(Omega_values)):
    Om = Omega_values[i]
    
    ax[0].plot(delta_values.T, en[i,:], label=labels[i], lw=lw, color=colors[i])
    ax[1].plot(delta_values.T,pop[i,:], lw=lw, color=colors[i])

ax[0].set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax[0].set_ylabel(r'$E_{released}/\omega_r$',fontsize = 14)
ax[0].legend()

yticks = np.linspace(0,1,5)
ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
ax[1].set_yticks(yticks)
ax[1].set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax[1].set_ylabel(r'$P_{\uparrow\uparrow}$',fontsize = 14)
ax[1].grid()


#%% size after tof
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(5,4))

for i in range(np.size(Omega_values)):
    Om = Omega_values[i]
    
    axS.plot(delta_values.T, size[i,:]*10**3, label=labels[i], lw=lw, color=colors[i],zorder =1)

axS.plot(scan_exp+0.25, size_exp, label='exp', lw=lw, color='black' ,zorder =1)

xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)    
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)  
axS.set_ylim(0,0.26)  
axS.yaxis.set_major_locator(MultipleLocator(0.05))
axS.yaxis.set_major_formatter('{x:.2f}')
axS.yaxis.set_minor_locator(MultipleLocator(0.01))
#axS.spines['left'].set_position(('data', 0))
#ax.spines['left'].set_zorder(10)
axS.legend()
axS.grid()

#%%
# fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\Energy&Spin.png', dpi = 300)
# figS.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\Size_tof.png', dpi = 300)
