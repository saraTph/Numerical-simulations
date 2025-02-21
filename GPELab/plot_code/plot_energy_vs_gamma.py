# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 21:25:51 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

#%%

mat = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\output_data.mat')


delta_values = mat.get('delta_values').squeeze()
Omega_values = mat.get('Omega_values').squeeze()
IE = mat.get('IE')
mu = mat.get('Chemical_potential')

wr = 169 #Hz
wz = 26  #Hz



a_bohr = 0.52917721e-10
a11 = (85*a_bohr)
a22 = (33.4*a_bohr)
a12 = (-53.0*a_bohr)

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

sigma_z = np.sqrt(hbar/(m*wr))
#n_values = np.linspace(1e8,1e10, 20);


gamma = np.zeros_like(Omega_values)


for i in range(len(Omega_values)):
    Om = Omega_values[i]
    n1D = 2.8e9
    T = 1/Om             # characteristic time
    L = 1/np.sqrt(wr)    # characteristic legth
    
    
    # define INTERACTION Adimensional parameters g22, g22, g12
    g11 = n1D * (4*np.pi* a11) * (T/L**2) 
    g22 = n1D * (4*np.pi* a22) * (T/L**2) 
    g12 = n1D * (4*np.pi* a12) * (T/L**2) 
    
    gamma[i] = (g11 + g22 - 0.5*g12)/4
    
gamma = gamma.squeeze()

#%%
colors = plt.get_cmap('Set3_r').colors
lw = 2
# Plotting the expressions
fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))

for i in range(len(delta_values)):
    d = delta_values[i]
    ax.plot(gamma, ((IE[:,i])*Om/wr)/(n1D*np.sqrt(2*np.pi)*sigma_z), label=fr'$\delta$= {d:.2g}', lw=lw, color=colors[i])

    

ax.set_xlabel(r'$\gamma$', fontsize=14)
ax.set_ylabel(r'$E_{int}/(N\hbar\omega_r)$',fontsize = 14)
ax.legend()

# yticks = np.linspace(0,1,5)
# ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
# ax[1].set_yticks(yticks)
# ax[1].set_xlabel(r'$\delta/\Omega$', fontsize=14)
# ax[1].set_ylabel(r'$P_{\uparrow\uparrow}$',fontsize = 14)
# ax[1].grid()
