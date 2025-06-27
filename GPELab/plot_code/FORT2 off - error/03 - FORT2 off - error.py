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
#location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\FORT2 off - error"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\FORT2 off - error\Export\FORT2 off - error"
file = os.path.join(location, 'optimal_Time-FORT2off.txt')
Omega_values, time_FORT2off = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


#location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\FORT2 off - error\Export\kinEnergy_sweep"

Omega_values = np.array([30400, 15200, 7600, 3800, 1900, 950])
len_sim = 64

delta_sweep = np.ones((len(Omega_values),len_sim))
F_sweep = np.ones((len(Omega_values),len_sim))
kinEn_sweep = np.ones((len(Omega_values),len_sim))

for i in range(len(Omega_values)):
    file_path = os.path.join(location, f'kinEn_sweep_{i}.txt')
    
    delta_sweep[i, :], F_sweep[i, :], kinEn_sweep[i, :] = np.genfromtxt(
        file_path, delimiter='\t', skip_header=1, comments='#', unpack=True
    )


# store energy for each Omega
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\FORT2 off - error\Export\DeltaE_sim"
num_files = 6
num_columns = 65 
en = np.zeros((num_files, num_columns))
size = np.zeros((num_files, num_columns))
pop = np.zeros((num_files, num_columns))
delta_values = np.zeros((num_files, num_columns))

for N in range(num_files):
    file_path = os.path.join(location, f'sim_{N}_4.6n0.txt')

    en[N, :], size[N, :], pop[N, :], delta_values[N, :] = np.genfromtxt(
        file_path, delimiter='\t', skip_header=1, comments='#', unpack=True
    )
    
del size; del pop; del delta_values
 
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

        
sweep = list(range(2, 66, 1))
energy_sweep = np.ones((len(Omega_values),len(sweep)))

for j, Om in enumerate(Omega_values):
    
    for idx, i in enumerate(sweep):
        
        d_scan = d[0:i]                # cut the scan of delta
        F_sweep_cut = F_sweep[j,0:i]   # cut the force given by interactions 
        F_scan = np.ones((len(d_scan)))
        
        
        #define dt to perforn integation 
        delta_diff = np.abs(np.diff(d_scan))
        dt = (delta_diff / np.sum(delta_diff)) * T
        time_integration = np.concatenate([[0], np.cumsum(dt)])
        
        
        F_scan = np.where(time_integration < time_FORT2off[j], - m * wz**2 * sigma_z + F_sweep_cut, - 0*m * wz**2 * sigma_z + F_sweep_cut)
        
        #F_scan = F[j,0:i]  # cut the scan of the force
            
        # integrate and evaluate kin energy
        delta_v = np.cumsum(F_scan[:-1] / m * dt)
        Ekin = 0.5*m*delta_v**2 
        #delta_sweep[j,idx] = d_scan[-1]
        energy_sweep[j,idx] = Ekin[-1]
        

        
#%%

def moving_average(arr, window_size):
    return np.convolve(arr, np.ones(window_size) / window_size, mode='valid')

lw=1.7
colors = plt.get_cmap('Set2').colors
labels = [r'$\Omega$ = 30.4 kHz', r'$\Omega$ = 15.2 kHz', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

fig1, ax1 = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))

for i, Om in enumerate(Omega_values):
    ax1.plot(delta_sweep[0,:], (energy_sweep[i,:]/en[i,0:-1])*100, lw=lw, color=colors[i], zorder=1, label=labels[i])

    
ax1.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax1.set_ylabel(r'Energy relative error $(\%) $', fontsize=14)
ax1.legend()


#%%
fig1.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\FORT2 off - error\Figures', dpi = 300)