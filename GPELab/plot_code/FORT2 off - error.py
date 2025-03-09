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

n_values = [5, 4.6]
len_sim = 52

en = np.ones((len(n_values),len_sim))
size = np.ones((len(n_values),len_sim))
pop = np.ones((len(n_values),len_sim))
delta_values = np.ones((len(n_values),len_sim))

location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error"

file = os.path.join(location, 'sim_5_5n0.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_5_4.6n0.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


#%%

wr = 169*2*np.pi #Hz
# wz = 26*2*np.pi  #Hz

# u = 1.66053906660e-27
# m = 39*u
hbar = 1.054571818e-34 #J s

#%% size after tof

colors = plt.get_cmap('Set2').colors
lw = 1.7
labels = [r'$n_0$ = 5e9', r'$n_0$ = 4.6e9', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

from matplotlib.ticker import AutoMinorLocator, MultipleLocator

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(8,5))

shift = [0.12, 0.122, 0.158, 0.257, 0.257, 0.376]
#for i in range(np.size(Omega_values)):

DenDn = (en[0,:]-en[1,:])*hbar*wr
DnDs = -((n_values[0]-n_values[1])/2)/((size[0]-size[1])/2)
# devo calcolare la teglia del BEC for the two densities and for a 0 tof
axS.plot(delta_values[0,:], DenDn*DnDs, lw=lw, color=colors[0], zorder =1)


xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)    
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\Delta E/\Delta \sigma$',fontsize = 14)  
#axS.set_ylim(0.02,0.25)  
#axS.legend()
axS.grid()

#%%
#figS.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\MF_trueOmega.png', dpi = 300)
