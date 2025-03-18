# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 10:20:13 2025

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os


#%%
#mat = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\output_data.mat')
mat = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\outputs\output_data.mat')

#import stimulaiton data
delta_values = mat.get('delta_values').squeeze()
Omega_values = mat.get('Omega_values').squeeze()
P_down  = mat.get('P_down')
IE = mat.get('IE')
KE = mat.get('KE')
RE = mat.get('RE')
PE = mat.get('PE')
e_tot = mat.get('energy_tot')
size = mat.get('sigma')


#%%
Om = Omega_values

wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

L = np.sqrt(hbar/(m*wr))
E = hbar*wr

#%% average over the TF density profile

#weights = np.array([0.1863, 0.1789, 0.1645, 0.1442, 0.1191, 0.0912, 0.0626, 0.0361, 0.0148, 0.0023]) 
weights = np.array([0.3651, 0.3087, 0.2103, 0.0987, 0.0171])
pop_avg = np.dot(weights, P_down)
e_tot_avg = np.dot(weights, e_tot)
RE_avg = np.dot(weights, RE)
IE_avg = np.dot(weights, IE)
PE_avg = np.dot(weights, PE)
KE_avg = np.dot(weights, KE)


e_rel = (e_tot_avg - RE_avg)*Om/wr - 1

#%% Extract energy and size 
colors = plt.get_cmap('Set3_r').colors
lw = 2
i = 2

t_tof = 62.6e-3  #s
e_rel_avg = (e_tot_avg-RE_avg)*(hbar*Om) - hbar*wr # energy tot - Rabi energy - hbar * wr

v = np.sqrt(2*e_rel_avg/m)   
sizeBEC = v * t_tof

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(5,4))
axS.plot(delta_values.T, sizeBEC*10**3, label='avg', lw=lw, color=colors[i],zorder =1)
 
xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)  
axS.legend()
axS.grid()


#%% Export
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\DeltaE_sim"
location = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\Export\FORT2 off-error\DeltaE_sim"

# outarray = np.vstack((e_rel_avg, sizeBEC, pop_avg, delta_values)).T
# header = 'energy \t size (m) \t pop \t delta_scan'
# np.savetxt(os.path.join(location, 'sim_4_4.6n0.txt'), outarray, header=header, delimiter='\t')