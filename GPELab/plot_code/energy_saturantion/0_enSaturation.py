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

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

a_bohr = 0.52917721e-10
a11 = (86.4014*a_bohr) 
a22 = (33.2755*a_bohr)
a12 = (-53.1022*a_bohr)

#%%

N = 41
#folder_path = r"C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\outputs\saturation_energy"
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

#calcuate gamma vector
gamma = np.ones(n_Om)
for idx, Om in enumerate(Omega_values):
    T = 1/Om             # characteristic time
    L = 1/np.sqrt(wr)    # characteristic legth
    L_z = 1/np.sqrt(wz) 
    
    # define INTERACTION Adimensional parameters g22, g22, g12
    n1D = 4.3e9
    g11 = (4*np.pi* a11)  
    g22 = (4*np.pi* a22)
    g12 = (4*np.pi* a12) 
    g_bar = (g11 + g22 - 2*g12)/4
    
    gamma[idx] =  g_bar*T/L**2*n1D*np.sqrt(1/wz)
    
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

lw = 2
# Plotting the expressions
# fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))
# ax.plot(gamma, energy[:,0], '-', label=fr'$\delta$= {delta_values[0]:.2g}', lw=lw, color='k')
# ax.plot(gamma, gamma/2, label=fr'$\delta$= {delta_values[0]:.2g}', lw=lw, color='r')
# ax.set_xlabel(r'$\gamma$', fontsize=14)
# ax.set_ylabel(r'$E_{int}/(\hbar\Omega)$', fontsize=14)

plt.imshow(energy, aspect='auto', origin='lower', cmap='Grays',
           extent=[delta_values[0], delta_values[-1], gamma[0], gamma[-1]])
plt.colorbar(label=r'$E/(\hbar\Omega)$')  # Adds color scale
plt.xlabel(r'$\delta$')
plt.ylabel(r'$\gamma$')
plt.title('Energy Map')
plt.show()

