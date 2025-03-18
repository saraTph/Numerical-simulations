# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 18:00:10 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%%
location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\FORT2 off - error"

file = os.path.join(location, 'optimal_Time-FORT2off.txt')
Omega_values, time_FORT2off = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"
#location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"

Omega_values = np.array([30400, 15200, 7600, 3800, 1900, 950])
len_sim = 64

delta_sweep = np.ones((len(Omega_values),len_sim))
F_sweep = np.ones((len(Omega_values),len_sim))
kinEn_sweep = np.ones((len(Omega_values),len_sim))

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

#%% extimate real energy given to the system in the adiabatic sweep considering the switch off of FORT2

d = delta_sweep[0,:]

T = 9e-3
time_sweep = np.linspace(0,T,64)


F = np.ones((len(Omega_values),len_sim))

for j, Om in enumerate(Omega_values):
    for idx, t in enumerate(time_sweep):
        
        F[j,idx] = np.where(t < time_FORT2off[j], - m * wz**2 * sigma_z + 0*F_sweep[j, idx], - m * wz**2 * sigma_z +0*F_sweep[j, idx])
        
sweep = list(range(2, 66, 1))
energy_sweep = np.ones((len(Omega_values),len(sweep)))

for j, Om in enumerate(Omega_values):
    
    for idx, i in enumerate(sweep):
        
        d_scan = d[0:i]  # cut the scan of delta
        F_scan = F[j,0:i]  # cut the scan of the force
        
        #difine dt to perforn integation 
        delta_diff = np.abs(np.diff(d_scan))
        dt = (delta_diff / np.sum(delta_diff)) * T
            
        # integrate and evaluate kin energy
        delta_v = np.cumsum(F_scan[:-1] / m * dt)
        Ekin = 0.5*m*delta_v**2 
        #delta_sweep[j,idx] = d_scan[-1]
        energy_sweep[j,idx] = Ekin[-1]
        
#%%
lw=1.7
colors = plt.get_cmap('Set2').colors
labels = [r'$\Omega$ = 30.4 kHz', r'$\Omega$ = 15.2 kHz', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

fig1, ax1 = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))

for i, Om in enumerate(Omega_values):
    ax1.scatter(delta_sweep[0,:], energy_sweep[0,:] /(hbar*2*np.pi), lw=lw, color=colors[i], zorder=1, label=labels[i])
    #ax1.scatter(delta_sweep[0,:], F[0,:] /(hbar*2*np.pi), lw=lw, color=colors[i], zorder=1, label=labels[i])
    
ax1.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax1.set_ylabel(r'$E_{kin} (Hz)$', fontsize=14)
ax1.legend()
ax1.set_title(r'energy given by $\frac{\partial E}{\partial \sigma_z}$ in the adiabatic sweep and for each $\delta/\Omega$')

#%%
# fig1, ax1 = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))
# ax1.plot(F[0,:])