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

#scan = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\exp_data\scan')
#data = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\exp_data\MFdata')
#mat = scipy.io.loadmat(r'C:\Users\sarat\OneDrive\Documenti\GitHub\Numerical-simulations\GPELab\outputs\output_data.mat')

#n = [2e9, 4e9,4.3e9,4.4e9,4.6e9,5e9]
#files = ["output_data_2.mat","output_data_4.mat", "output_data_4_3.mat", "output_data_4_4.mat",  "output_data_4_6.mat", "output_data_5.mat"]
files = "output_data.mat"
folder_path = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\saturation_energy"
path = os.path.join(folder_path, files)

mat = scipy.io.loadmat(path)
IE = mat.get('IE')
IE = IE[0:5,:]
Omega_values = mat.get('Omega_values').squeeze()  #Omega frequency
delta_values = mat.get('delta_values').squeeze()

# IE_avg = np.ones(len(Omega_values))
en_density = np.ones(len(Omega_values))
gamma = np.ones(len(Omega_values))

weights = np.array([0.3651, 0.3087, 0.2103, 0.0987, 0.0171])
en = np.dot(weights, IE)





for i in range(len(Omega_values)):
    
    
    
    Om = Omega_values[i]
    
    T = 1/Om             # characteristic time
    L = 1/np.sqrt(wr)    # characteristic legth
    L_z = 1/np.sqrt(wz) 
    # define INTERACTION Adimensional parameters g22, g22, g12
    n1D = 4.3e9
    g11 = n1D * (4*np.pi* a11) * (T/L**2) 
    g22 = n1D * (4*np.pi* a22) * (T/L**2) 
    g12 = n1D * (4*np.pi* a12) * (T/L**2) 
    
    gamma[i] = (g11 + g22 - 2*g12)/4* np.sqrt(2*np.pi)*L_z
    
    en_density[i] = en[i]



#%%

lw = 2
# Plotting the expressions
fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))
ax.plot(gamma, en_density, label=fr'$\delta$= {delta_values:.2g}', lw=lw, color='k')
ax.plot(gamma, gamma/2, label=fr'$\delta$= {delta_values:.2g}', lw=lw, color='r')
ax.set_xlabel('$\gamma$', fontsize=14)
ax.set_ylabel('$E_{int}/(\hbar\Omega)$', fontsize=14)
ax.set_xlim(0,1.2)
ax.set_ylim(0,0.8)



    
#plt.plot(gamma,en)
