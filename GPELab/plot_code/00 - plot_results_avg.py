# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 15:11:15 2025

@author: sarat
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os


#%%
#mat = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\output_data_4.mat')
mat = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\output_data.mat')

delta_values = mat.get('delta_values').squeeze()
Omega_values = mat.get('Omega_values').squeeze()
P_down  = mat.get('P_down')
IE = mat.get('IE')
KE = mat.get('KE')
RE = mat.get('RE')
PE = mat.get('PE')
e_tot = mat.get('energy_tot')

# data_exp5 = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\data_5.mat') 
# size_exp5 = data_exp5.get('valY').squeeze()
# scan_exp5 = data_exp5.get('valX').squeeze()



scan = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\scan') 
scan_exp0 = scan.get('valX').squeeze()
data = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\MFdata') 
size_exp0 = data.get('data5').squeeze()

#%%
Om = Omega_values

wr = 169*2*np.pi #Hz
wz = 26*2*np.pi  #Hz

u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

#%% average over the TF density profile

#weights = np.array([0.1863, 0.1789, 0.1645, 0.1442, 0.1191, 0.0912, 0.0626, 0.0361, 0.0148, 0.0023]) 
weights = np.array([0.3651, 0.3087, 0.2103, 0.0987, 0.0171])
pop_avg = np.dot(weights, P_down)
e_tot_avg = np.dot(weights, e_tot)
RE_avg = np.dot(weights, RE)
IE_avg = np.dot(weights, IE)
PE_avg = np.dot(weights, PE)
KE_avg = np.dot(weights, KE)

#%% Plots energy and pop

colors = plt.get_cmap('Set3_r').colors
lw = 2
i = 2

fig, ax = plt.subplots(1,2,constrained_layout=True, figsize=(9,4))

e_rel = (e_tot_avg - RE_avg)*Om/wr - 1
ax[0].plot(delta_values.T, e_rel, label='avg', lw=lw, color=colors[i])
ax[1].plot(delta_values.T, pop_avg, lw=lw, color=colors[i])

ax[0].set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax[0].set_ylabel(r'$E_{released}/\hbar\omega_r$',fontsize = 14)
ax[0].legend()
yticks = np.linspace(0,1,5)
ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
ax[1].set_yticks(yticks)
ax[1].set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax[1].set_ylabel(r'$P_{\uparrow\uparrow}$',fontsize = 14)
ax[1].grid()


#%% Plot size after tof
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
lw =2

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(5,4))

t_tof = 62.6e-3  #s

e_rel_avg = (e_tot_avg-RE_avg)*(hbar*Om) - hbar*wr # energy tot - Rabi energy - hbar * wr

v = np.sqrt(2*e_rel_avg/m)   
sizeBEC = v * t_tof


axS.plot(delta_values.T, sizeBEC*10**3, label='avg', lw=lw, color=colors[i],zorder =1)

#axS.scatter(scan_exp5+0.25, size_exp5, marker = '.' ,label='exp', lw=lw, color='black' ,zorder =1)
axS.scatter(scan_exp0 +  0.09, size_exp0, marker = '.' ,label='exp', lw=lw, color='black' ,zorder =1)
    
xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)  
axS.set_ylim(0,0.26)  
#axS.yaxis.set_major_locator(MultipleLocator(0.05))
#axS.yaxis.set_major_formatter('{x:.2f}')
#axS.yaxis.set_minor_locator(MultipleLocator(0.01))
#axS.spines['left'].set_position(('data', 0))
#axS.spines['left'].set_zorder(10)
axS.legend()
axS.grid()

#fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\energy_spin.png', dpi = 300)
#figS.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\size.png', dpi = 300)

#%% Export
location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export"

# outarray = np.vstack((e_rel, sizeBEC, pop_avg, delta_values)).T
# header = 'energy \t size (m) \t pop \t delta_scan'
# np.savetxt(os.path.join(location, 'sim_3.txt'), outarray, header=header, delimiter='\t')