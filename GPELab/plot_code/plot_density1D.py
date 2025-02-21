# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 10:35:40 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

#%%

mat = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\output_data.mat')


delta_values = mat.get('delta_values').squeeze()
Omega_values = mat.get('Omega_values').squeeze()
# P_down  = mat.get('P_down').squeeze()
# IE = mat.get('IE')
density1 = mat.get('densityProfile1D_comp1')
density2 = mat.get('densityProfile1D_comp2')

wr = 200 #Hz
wz = 50  #Hz
N = 5e6

a_bohr = 52.917721e-11
a11 = (85*a_bohr)
a22 = (33.4*a_bohr)
a12 = (-53.0*a_bohr)

u = 1.66053906660e-27
m = 39*u


#%%
colors = plt.get_cmap('Set3_r').colors
lw = 2

#labels = ['$\Omega = 2\omega_r$', '$\Omega = 4\omega_r$', '$\Omega \gg \omega_r$']
# Plotting the expressions
fig, ax = plt.subplots(1,2,constrained_layout=True, figsize=(9,4))

#for i in range(len(Omega_values)):
n_values = [2e7, 2e8, 2e9, 2e10]
for i in range(4):
    #Om = Omega_values[i]
    Om = n_values[i]
    #ax[0].plot(density1[i,:], lw=1,label=fr'$\Omega$= {Omega_values[i]:.3g} Hz')
    ax[0].plot(density1[i,:], lw=1,label=fr'n= {n_values[i]:.3g}')
    ax[1].plot(density2[i,:], lw=1)

ax[0].set_xlabel(r'x', fontsize=14)
ax[0].set_ylabel(r'$n_{1D}$ comp 1',fontsize = 14)
ax[0].legend()

#yticks = np.linspace(0,1,5)
#ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
#ax[1].set_yticks(yticks)
ax[1].set_xlabel(r'x', fontsize=14)
ax[1].set_ylabel(r'$n_{1D}$ comp 2',fontsize = 14)
ax[1].grid()

#%%
#fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\scan_gmax12\Energy&Spin.png', dpi = 300)

