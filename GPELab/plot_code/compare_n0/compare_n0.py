# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 15:41:29 2025

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

#load experimntal data and simulation results
data_energy = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\exp_data_energy') 
data_spin = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\exp_data_spin') 
loc_sim = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\GPE-output-analysis\Export\Omega30kHz_n0scan"

Om = 30400;
n0_values = [3.6, 4.6, 4.61, 4]

#initialize simulation results
len_sim = 52
len_energyData = 39
len_spinData = 25

en = np.ones((len(n0_values),len_sim))
size = np.ones((len(n0_values),len_sim))
pop = np.ones((len(n0_values),len_sim))
delta_values = np.ones((len(n0_values),len_sim))

#initialize data exp
size_exp = np.ones((len(n0_values),len_energyData))
scan_energy = np.ones((len_energyData))
spin_exp = np.ones((len(n0_values),len_spinData))
scan_spin = np.ones((len_spinData))


file = os.path.join(loc_sim, 'sim_3_6.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_4_6.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(loc_sim, 'sim_4_6_1.txt')
en[2,:], size[2,:], pop[2,:], delta_values[2,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)


file = os.path.join(loc_sim, 'sim_4.txt')
en[3,:], size[3,:], pop[3,:], delta_values[3,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

# import experimental data
scan_energy = data_energy.get('valX').squeeze()
# size_exp[5,:] = data_energy.get('data0') 
# size_exp[4,:] = data_energy.get('data1')
# size_exp[3,:] = data_energy.get('data2')
# size_exp[2,:] = data_energy.get('data3') 
# size_exp[1,:] = data_energy.get('data4') 
# size_exp[0,:] = data_energy.get('data5') 

#%% plots parameters

lw = 1.5
fs = 14
ls = 10

size_x = 5
size_y = 4

from matplotlib import cm

# Choose a sequential colormap (e.g., Greens)
colors = plt.get_cmap('Set1').colors
# colors = [cmap(x) for x in np.linspace(0.3, 0.9, 6)]  # Avoid very light and very dark ends
# colors = [cmap(x) for x in np.linspace(1, 0.4, 6)]

labels = [r'$n_0$ = 3.6e9', r'$n_0$ = 4.6e9', r'$n_0$ = 4.6_1e9', r'$n_0$ = 4']

#%% size after tof

from matplotlib.ticker import AutoMinorLocator, MultipleLocator

fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(size_x, size_y))

shift = [0.12, 0.122, 0.158, 0.257, 0.257, 0.376]
#for i in range(np.size(Omega_values)):
    
for i in range(0,4):
    
    ax.plot(delta_values[i,:], size[i,:]*10**3, lw=lw, color=colors[i],zorder =1,label=labels[i])
    # ax.scatter(scan_energy + shift[i], size_exp[i,:], lw=2, marker = '2', color=colors[i], label=labels[i],)

ax.set_xlim(-6.1,6)  
xticks = np.linspace(-6,6,7)
ax.set_xticks(xticks)    
ax.tick_params(labelsize=ls)
ax.set_xlabel(r'$\delta/\Omega$', fontsize=fs)
ax.set_ylabel(r'$\sigma(t_{TOF}= 62.6\ ms)$ (mm)',fontsize = fs)  
ax.set_ylim(-0.01,0.25)  
ax.yaxis.set_major_locator(MultipleLocator(0.05))
ax.yaxis.set_major_formatter('{x:.2f}')
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
#axS.spines['left'].set_position(('data', 0))
#ax.spines['left'].set_zorder(10)
legend = ax.legend( 
    loc='center left',
    fontsize=ls,              # Legend font size
    frameon=True,            # Show legend box
    fancybox=False,           # Rounded box corners
    framealpha=0.9,          # Slight transparency
    edgecolor='gray',        # Border color
    facecolor='white'        # Background color
)
ax.axhline(y=0, color='k', linestyle='--', linewidth=0.8)
ax.axvline(x=0, color='k', linestyle='--', linewidth=0.8)
