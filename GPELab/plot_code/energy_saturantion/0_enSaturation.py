# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 15:51:49 2025

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%%
wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz

# u = 1.66053906660e-27
# m = 39*u
# hbar = 1.054571818e-34 #J s

a_bohr = 0.52917721e-10
a11 = (86.4014*a_bohr) 
a22 = (33.2755*a_bohr)
a12 = (-53.1022*a_bohr)

# define INTERACTION Adimensional parameters g22, g22, g12
n1D = 4.3e9
g11 = (4*np.pi* a11)  
g22 = (4*np.pi* a22)
g12 = (4*np.pi* a12) 





#%%

N = 41
# folder_path = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\outputs\saturation_energy"
folder_path = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\saturation_energy"
files_name = [f"output_data_{n}.mat" for n in range(1, N)]

# import delta_values and Omega_values
path = os.path.join(folder_path, files_name[0])
mat = scipy.io.loadmat(path)
Omega_values = mat.get('Omega_values').squeeze()  #Rabi frequency
delta_values = mat.get('delta_values').squeeze()  #detuning

# define dimensions of data
n_Om = len(Omega_values)
n_delta = len(delta_values)

energy = np.ones((n_Om,n_delta))
for n, file_name in enumerate(files_name):
    
    path = os.path.join(folder_path, file_name)
    mat = scipy.io.loadmat(path)
    RE = mat.get('RE')
    e_tot = mat.get('energy_tot')

    # perform average over the Thomas-Fermi density profile
    weights = np.array([0.3651, 0.3087, 0.2103, 0.0987, 0.0171])
    RE_avg = np.dot(weights,RE)
    e_tot_avg = np.dot(weights,e_tot) 
    
    # evaluate energy
    energy[:,n] = (e_tot_avg - RE_avg)

#%%
#calcuate gamma vector and energy 2-body and 3-body


# calculate INTERACTION adimensional 2-body parameters
g_bar = (g11+g22-2*g12)/4 
g_inf = (g11*g22-g12**2)/(4*g_bar) 
k = (g11-g22)/(4*g_bar)
phi = np.linspace(0,np.pi,40)
g2 = g_inf+g_bar*(np.cos(phi)-k)**2


g3 = np.ones((n_Om,n_delta))
gamma_2b = np.ones((n_Om,n_delta))
gamma_3b = np.ones((n_Om,n_delta))

for idx, Om in enumerate(Omega_values):
    T = 1/Om             # characteristic time
    L = np.sqrt(1/wr)    # characteristic legth 
      
    gamma_2b[idx,:] = n1D * g2 *T/L**2
    
    g3[idx,:] = -3*g_bar**2/Om * (np.sin(phi))**3*(np.cos(phi)-k)**2
    gamma_3b[idx,:] =  n1D * g3[idx,:] *T/L**2
    
energy_2b = gamma_2b/2
energy_3b = gamma_3b/3*n1D


#%%

lw=2
ls = 14
# Plotting the expressions
fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))
delta = 30
ax.plot(gamma_2b[:,delta], energy[:,delta], '-', label=fr'energy; $\delta/\Omega$= {delta_values[0]:.2g}', lw=lw, color='r')
ax.plot(gamma_2b[:,delta], energy_2b[:,delta] + energy_3b[:,delta], label=fr'energy 2-body + 3-body; $\delta/\Omega$= {delta_values[0]:.2g}', lw=lw, color='k')
ax.set_xlabel(r'$\gamma$', fontsize=14)
ax.set_ylabel(r'$E_{int}/(\hbar\Omega)$', fontsize=14)
ax.legend()




tolerance = 0.1  # 10%
mask = np.abs(energy - (energy_2b+energy_3b)) < tolerance * energy

fig2D, ax2D = plt.subplots(1, 1, constrained_layout=True, figsize=(5, 4))
img = ax2D.imshow(energy, aspect='auto', origin='lower', cmap='viridis',
                  extent=[delta_values[0], delta_values[-1], Omega_values[0]/(2*np.pi)/1000, Omega_values[-1]/(2*np.pi)/1000])

# X, Y = np.meshgrid(delta_values, Omega_values)
# contour = ax2D.contour(X, Y, mask, levels=[0.5], colors='red', linewidths=2)

fig2D.colorbar(img, ax=ax2D, label=r'$E/N (\hbar\Omega)$')
ax2D.set_xlabel(r'$\delta$', size = ls)
ax2D.set_ylabel(r'$\Omega$ (kHz)', size = ls)
ax2D.set_title('Energy Map', size = ls)


#%%
# fig2D.savefig(r'C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\plot_code\energy_saturantion\Figures\energyMap', dpi = 300)