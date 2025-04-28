# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 14:25:01 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io

#%%

mat = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\output_data.mat')


delta_values = mat.get('delta_values').squeeze()
Omega_values = mat.get('Omega_values').squeeze()
P_down  = mat.get('P_down').squeeze()
IE = mat.get('IE').squeeze()
KE = mat.get('KE').squeeze()
RE = mat.get('RE').squeeze()
PE = mat.get('PE').squeeze()
e_tot = mat.get('energy_tot').squeeze()

i = 0
P_down = P_down[i,:]
IE = IE[i,:]
KE = KE[i,:]
RE = RE[i,:]
PE = PE[i,:]
e_tot = e_tot[i,:]
Om = Omega_values[i]





wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz
#n1D = 2.3e9

# a_bohr = 0.52917721e-10
# a11 = (85*a_bohr)
# a22 = (33.4*a_bohr)
# a12 = (-53.0*a_bohr)

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s
  


#%%
colors = plt.get_cmap('Set3').colors
colors2 = plt.get_cmap('Dark2').colors
lw = 2

labels = [r'$\Omega = 4\omega_r$', r'$\Omega = 8\omega_r$', r'$\Omega \gg \omega_r$']
# Plotting the expressions
fig, ax = plt.subplots(1,2,constrained_layout=True, figsize=(9,4))


#ax[0].plot(delta_values.T, (IE[i,:]+PE[i,:]+KE[i,:])*Om/wr,label=labels[i], lw=lw, color=colors[i])
#e_rel = (e_tot[i,:]-RE[i,:]-hbar*wr)*Om/wr
#ax[0].plot(delta_values.T, e_rel, label=labels[i], lw=lw, color=colors[i])
ax[0].plot(delta_values.T, IE*Om/wr, label=r"$E_{int}$", lw=lw, color=colors2[1])
ax[0].plot(delta_values.T, KE*Om/wr, label=r'$E_{kin}$', lw=lw, color=colors2[2])
ax[0].plot(delta_values.T, PE*Om/wr, label=r'$E_{trap}$', lw=lw, color=colors2[3])
#ax[0].plot(delta_values.T, np.ones(len(IE))*(hbar*wr)/wr, label=r'$E_{w}$', lw=lw, color=colors2[4])
ax[0].plot(delta_values.T, (IE+KE+PE)*Om/wr - 1, label=r'$E_{tot}-E_{rabi}-\hbar\omega_r$', lw=lw, color=colors2[5])
ax[1].plot(delta_values.T,P_down, lw=lw, color='black')

ax[0].set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax[0].set_ylabel(r'$E/\omega_r$',fontsize = 14)
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



t_tof = 62.6e-3  #s
   
e_rel = (IE+KE+PE)*(hbar*Om) - hbar*wr  # energy tot - Rabi energy - hbar * wr
#e_rel = (KE[i,:]+IE[i,:]+PE[i,:])*(hbar*Om)
v = np.sqrt(2*e_rel/m)   
sizeBEC = v * t_tof
    
axS.plot(delta_values.T, sizeBEC*10**3, lw=lw, color=colors[6])
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)
axS.yaxis.set_major_locator(MultipleLocator(0.05))
axS.yaxis.set_major_formatter('{x:.2f}')
axS.yaxis.set_minor_locator(MultipleLocator(0.01))
#%%
# fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\EnergyDecomposition.png', dpi = 300)
