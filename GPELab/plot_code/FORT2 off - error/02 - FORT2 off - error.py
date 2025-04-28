
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 12:02:25 2025

@author: sarat

This program uploads the results from "01 - FORT 2 - error" where the kinetic energy accumulated during the sweep because of interaction has been calculated

Here we calculate at which time the FORT 2 trap should be turned off to do not excite the cloud at delta/Omega = 0.25 (min of interaction)

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

Omega_values = np.array([30400, 15200, 7600, 3800, 1900, 950])
len_sim = 64

delta_sweep = np.ones((len(Omega_values),len_sim))
F_sweep = np.ones((len(Omega_values),len_sim))
kinEn_sweep = np.ones((len(Omega_values),len_sim))

#location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"

file = os.path.join(location, 'kinEn_sweep_0.txt')
delta_sweep[0,:], F_sweep[0,:], kinEn_sweep[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'kinEn_sweep_1.txt')
delta_sweep[1,:], F_sweep[1,:], kinEn_sweep[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'kinEn_sweep_2.txt')
delta_sweep[2,:], F_sweep[2,:], kinEn_sweep[2,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'kinEn_sweep_3.txt')
delta_sweep[3,:], F_sweep[3,:], kinEn_sweep[3,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'kinEn_sweep_4.txt')
delta_sweep[4,:], F_sweep[4,:], kinEn_sweep[4,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'kinEn_sweep_5.txt')
delta_sweep[5,:], F_sweep[5,:], kinEn_sweep[5,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


#%%

wr = 169*2*np.pi # Hz
wz = 25.5*2*np.pi  # Hz
sigma_z = 10e-6

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s


#%% 

time_FORT2off = np.ones(((len(Omega_values))))
idx = int(np.where(delta_sweep[0,:] == 0.25)[0][0])

for i, Om in enumerate(Omega_values):
    
    Om = Om *2*np.pi
    time_FORT2off[i] = np.sqrt( kinEn_sweep[i,idx] * (2/m) * (wz**2*sigma_z)**(-2))
    #print(kinEn_sweep[i,idx]/(hbar*2*np.pi))
    print(time_FORT2off[i])

#%%
colors = plt.get_cmap('Set2').colors
lw = 1.7

fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 3.5))

ax.plot(Omega_values*1e-3, time_FORT2off*1e3, lw=lw, color=colors[0], zorder=1)
ax.set_xlabel(r'$\Omega$ (kHz)', fontsize=14)  # Time on bottom
ax.set_ylabel(r'$t$ (ms)', fontsize=14)


#%%

labels = [r'$\Omega$ = 30.4 kHz', r'$\Omega$ = 15.2 kHz', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

fig1, ax1 = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))

for i, Om in enumerate(Omega_values):
    ax1.scatter(delta_sweep[0,:], kinEn_sweep[i,:] /(hbar*2*np.pi), lw=lw, color=colors[i], zorder=1, label=labels[i])
    
ax1.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax1.set_ylabel(r'$E_{kin} (Hz)$', fontsize=14)
ax1.legend()
ax1.set_title(r'energy given by $\frac{\partial E}{\partial \sigma_z}$ in the adiabatic sweep and for each $\delta/\Omega$')

         

#%%
#figS.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\FORT2 off\Ekin_t.png', dpi = 300)

#location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\FORT2 off - error"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\FORT2 off - error"

# outarray = np.vstack((Omega_values, time_FORT2off)).T
# header = 'Rabi freq \t optimal time FORT2 off'
# np.savetxt(os.path.join(location, 'optimal_Time-FORT2off.txt'), outarray, header=header, delimiter='\t')