# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 17:46:41 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

#load experimntal data and simulation results
data_energy = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\exp_data_energy') 
data_spin = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\exp_data_spin') 
loc_sim = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\GPE-output-analysis\Export\MF_tiemann_4.6n0_TrueOmega"

Omega_values = [30400, 15200, 7600, 3800, 1900, 950]

#initialize simulation results
len_sim = 52
len_energyData = 39
len_spinData = 25

en = np.ones((len(Omega_values),len_sim))
size = np.ones((len(Omega_values),len_sim))
pop = np.ones((len(Omega_values),len_sim))
delta_values = np.ones((len(Omega_values),len_sim))

#initialize data exp
size_exp = np.ones((len(Omega_values),len_energyData))
scan_energy = np.ones((len_energyData))
spin_exp = np.ones((len(Omega_values),len_spinData))
scan_spin = np.ones((len_spinData))


file = os.path.join(loc_sim, 'sim_0.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_1.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_2.txt')
en[2,:], size[2,:], pop[2,:], delta_values[2,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_3.txt')
en[3,:], size[3,:], pop[3,:], delta_values[3,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_4.txt')
en[4,:], size[4,:], pop[4,:], delta_values[4,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_5.txt')
en[5,:], size[5,:], pop[5,:], delta_values[5,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

# import experimental data
scan_energy = data_energy.get('valX').squeeze()
size_exp[5,:] = data_energy.get('data0') 
size_exp[4,:] = data_energy.get('data1')
size_exp[3,:] = data_energy.get('data2')
size_exp[2,:] = data_energy.get('data3') 
size_exp[1,:] = data_energy.get('data4') 
size_exp[0,:] = data_energy.get('data5') 

scan_spin = data_spin.get('scan')
spin_exp[5,:] = data_spin.get('pop5').squeeze()
spin_exp[4,:] = data_spin.get('pop4').squeeze()
spin_exp[3,:] = data_spin.get('pop3').squeeze()
spin_exp[2,:] = data_spin.get('pop2').squeeze()
spin_exp[1,:] = data_spin.get('pop1').squeeze()
spin_exp[0,:] = data_spin.get('pop0').squeeze()


#%%
colors = plt.get_cmap('Set2').colors
lw = 1.7
labels = [r'$\Omega$ = 30.4 kHz', r'$\Omega$ = 15.2 kHz', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

# Plotting the expressions
fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))


#for i in range(np.size(Omega_values)):
for i in range(6):
    Om = Omega_values[i] *2*np.pi
    
    ax.plot(delta_values[i,:], pop[i,:], lw=lw, color=colors[i], label=labels[i])
    ax.scatter(scan_spin, spin_exp[i,:], lw=2, marker = '2', color=colors[i])


yticks = np.linspace(0,1,5)
ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
ax.set_yticks(yticks)
ax.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax.set_ylabel(r'$P_{\uparrow\uparrow}$',fontsize = 14)
ax.legend()
ax.grid()


#%% size after tof

colors = plt.get_cmap('Set2').colors
lw = 1.7
labels = [r'$\Omega$ = 30.4 kHz', r'$\Omega$ = 15.2 kHz', r'$\Omega$ = 7.6 kHz', r'$\Omega$ = 3.8 kHz', r'$\Omega$ = 1.9 kHz', r'$\Omega$ = 950 Hz']

from matplotlib.ticker import AutoMinorLocator, MultipleLocator

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(10,7))

shift = [0.12, 0.122, 0.158, 0.257, 0.257, 0.376]
#for i in range(np.size(Omega_values)):
for i in range(0,6):
    Om = Omega_values[i]*2*np.pi
    
    axS.plot(delta_values[i,:], size[i,:]*10**3, label=labels[i], lw=lw, color=colors[i],zorder =1)
    axS.scatter(scan_energy + shift[i], size_exp[i,:], lw=2, marker = '2', color=colors[i])

xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)    
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)  
axS.set_ylim(0.02,0.25)  
axS.yaxis.set_major_locator(MultipleLocator(0.05))
axS.yaxis.set_major_formatter('{x:.2f}')
axS.yaxis.set_minor_locator(MultipleLocator(0.01))
#axS.spines['left'].set_position(('data', 0))
#ax.spines['left'].set_zorder(10)
axS.legend()
axS.grid()

#%%

#fig.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\Spin_DV.png', dpi = 300)
#figS.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\Energy_DV.png', dpi = 300)