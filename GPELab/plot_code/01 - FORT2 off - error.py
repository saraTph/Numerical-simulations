# -*- coding: utf-8 -*-
"""
Created on Sun Mar  9 12:02:25 2025

@author: sarat

this program uploqds the results from the simulations GPELab runned for two different walues of densities
and calculates the Force generated during the adiabatic sweep fron delta/Omega = 8 to -8 

This force is calculated for each detuning value and considering a sweep duration of 9 ms (equal for each point)

From the force is calculated the Kinetic energy eccumulated during the sweep due to ONLY interaction term

Results are saved in txt files
 
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

n_values = [4.6e9, 4.554e9]
len_sim = 65

en = np.ones((len(n_values),len_sim))
size = np.ones((len(n_values),len_sim))
pop = np.ones((len(n_values),len_sim))
delta_values = np.ones((len(n_values),len_sim))

#location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error"
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\DeltaE_sim"

file = os.path.join(location, 'sim_0_4.6n0.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_0_4.554n0.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


#%%

wr = 169*2*np.pi # Hz
wz = 25.5*2*np.pi  # Hz

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

L = np.sqrt(hbar/(m*wr))
E = hbar*wr
sigma_z = 10e-6

#%%


Den = (en[0,:]-en[1,:])  # J
DenDn = Den/(n_values[0]-n_values[1])
n_avg = (n_values[0]+n_values[1])/2
DnDs = -(n_avg/sigma_z)

N = np.sqrt(2*np.pi) * sigma_z * n_avg
F = - 0*m*wz**2*sigma_z - DenDn*DnDs
# invert scans (from delta = 8 to delta = -8)
d = delta_values[0,:]
d = d[::-1]  # invert delta vector [+8 to -8]
F = F[::-1]  # invert Force vector [+8 to -8]

sweep = list(range(2, 66, 1))
delta_sweep = np.ones((len(sweep)))
energy_sweep = np.ones((len(sweep)))

for idx, i in enumerate(sweep):
    
    d_scan = d[0:i]  # cut the scan of delta 
    F_scan = F[0:i]  # cut the scan of the force

    #difine dt to perforn integation 
    delta_diff = np.abs(np.diff(d_scan))
    T = 9e-3 
    dt = (delta_diff / np.sum(delta_diff)) * T
    time_vector = np.concatenate([[0], np.cumsum(dt)])

    # integrate and evaluate kin energy
    delta_v = np.cumsum(F_scan[:-1] / m * dt)
    Ekin = 0.5*m*delta_v**2
    
    delta_sweep[idx] = d_scan[-1]
    energy_sweep[idx] = Ekin[-1]

#%%
# colors = plt.get_cmap('Set2').colors
# lw = 1.7

# figS, axT = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))

# axT.plot(time_vector[:-1] * 1e3, Ekin / (hbar * 2 * np.pi), lw=lw, color=colors[0], zorder=1)
# axT.set_xlabel(r'$t$ (ms)', fontsize=14)  # Time on bottom
# axT.set_ylabel(r'$E_{kin} (Hz)$', fontsize=14)


# delta_ticks = np.linspace(8, -8, num=17)  # Adjust `num` for more/less ticks
# time_ticks = np.interp(delta_ticks, d_scan[::-1], time_vector[::-1]) * 1e3  # Convert to ms
# axS = axT.secondary_xaxis('top')
# axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
# axS.set_xticks(time_ticks)  # Use interpolated time positions
# axS.set_xticklabels([f"{tick:.1f}" for tick in delta_ticks])  # Show delta values

# axT.grid(True)

#%%
fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(8, 5))
ax.scatter(delta_sweep,energy_sweep/(hbar * 2 * np.pi))
ax.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax.set_ylabel(r'$E_{kin} (Hz)$', fontsize=14)
ax.grid(True)

#%%

# location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\kinEnergy_sweep"
# outarray = np.vstack((delta_sweep, energy_sweep)).T
# header = 'delta_sweep \t kin_en_sweep'
# np.savetxt(os.path.join(location, 'kinEn_sweep_0.txt'), outarray, header=header, delimiter='\t')
