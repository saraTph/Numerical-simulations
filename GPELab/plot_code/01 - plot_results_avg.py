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

Omega_values = [22300, 11200, 5600, 2810, 1410, 706]

#initialize simulation results
en = np.ones((len(Omega_values),52))
size = np.ones((len(Omega_values),52))
pop = np.ones((len(Omega_values),52))
delta_values = np.ones((len(Omega_values),52))

#initialize data exp
size_exp = np.ones((len(Omega_values),39))
scan_exp = np.ones((39))

location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export"
#location = r"C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Export\MF_a11+2"
file = os.path.join(location, 'sim_0.txt')
en[0,:], size[0,:], pop[0,:], delta_values[0,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_1.txt')
en[1,:], size[1,:], pop[1,:], delta_values[1,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_2.txt')
en[2,:], size[2,:], pop[2,:], delta_values[2,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_3.txt')
en[3,:], size[3,:], pop[3,:], delta_values[3,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

file = os.path.join(location, 'sim_4.txt')
en[4,:], size[4,:], pop[4,:], delta_values[4,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

# file = os.path.join(location, 'sim_5.txt')
# en[5,:], size[5,:], pop[5,:], delta_values[5,:] = np.genfromtxt(file, delimiter='\t', skip_header=1, comments='#', unpack=True)

# import experimental data
scan = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\scan') 
scan_exp[:] = scan.get('valX').squeeze()

data = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\exp_data\MFdata') 
size_exp[5,:] = data.get('data0').squeeze()
size_exp[4,:] = data.get('data1').squeeze() 
size_exp[3,:] = data.get('data2').squeeze()
size_exp[2,:] = data.get('data3').squeeze()
size_exp[1,:] = data.get('data4').squeeze()
size_exp[0,:] = data.get('data5').squeeze()
#%%

# wr = 169*2*np.pi #Hz
# wz = 26*2*np.pi  #Hz

# u = 1.66053906660e-27
# m = 39*u
# hbar = 1.054571818e-34 #J s

# a_bohr = 0.52917721e-10
# a11 = (85*a_bohr)
# a22 = (33.4*a_bohr)
# a12 = (-53.0*a_bohr)
# g11 = (4*np.pi*hbar**2*a11*a_bohr)/m
# g12 = (4*np.pi*hbar**2*a12*a_bohr)/m 
# g22 = (4*np.pi*hbar**2*a22*a_bohr)/m
  
# V = 1/(0.05e-6**3)

#%%
colors = plt.get_cmap('Set2').colors
lw = 1

labels = [r'$\Omega$ = 22.3 kHz', r'$\Omega$ = 11.2 kHz', r'$\Omega$ = 5.6 kHz', r'$\Omega$ = 2.81 kHz', r'$\Omega$ = 1.41 kHz', r'$\Omega$ = 706 Hz']

# # Plotting the expressions
# fig, ax = plt.subplots(1,2,constrained_layout=True, figsize=(9,4))


# #for i in range(np.size(Omega_values)):
# for i in range(4):
#     Om = Omega_values[i] *2*np.pi
    
#     ax[0].plot(delta_values[i,:], en[i,:], label=labels[i], lw=lw, color=colors[i])
#     ax[1].plot(delta_values[i,:], pop[i,:], lw=lw, color=colors[i])

# ax[0].set_xlabel(r'$\delta/\Omega$', fontsize=14)
# ax[0].set_ylabel(r'$E_{released}/\omega_r$',fontsize = 14)
# ax[0].legend()

# yticks = np.linspace(0,1,5)
# ytick_labels = [r'$0$', r'$0.25$', r'$0.5$', r'$0.75$' , r'$1$']
# ax[1].set_yticks(yticks)
# ax[1].set_xlabel(r'$\delta/\Omega$', fontsize=14)
# ax[1].set_ylabel(r'$P_{\uparrow\uparrow}$',fontsize = 14)
# ax[1].grid()


#%% size after tof
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

figS, axS = plt.subplots(1,1,constrained_layout=True, figsize=(10,7))

shift = [0.12, 0.122, 0.158, 0.257, 0.257, 0.376]
#for i in range(np.size(Omega_values)):
for i in range(5):
    Om = Omega_values[i]*2*np.pi
    
    axS.plot(delta_values[i,:], size[i,:]*10**3, label=labels[i], lw=lw, color=colors[i],zorder =1)
    axS.scatter(scan_exp + shift[i], size_exp[i,:], lw=1.2, marker = '2', color=colors[i])
    #print(np.min(size_exp[i,:]))
    print(np.min( size[i,:]*10**3))

#axS.plot(scan_exp+0.25, size_exp, label='exp', lw=lw, color='black' ,zorder =1)

xticks = np.linspace(-8,8,17)
axS.set_xticks(xticks)    
axS.set_xlabel(r'$\delta/\Omega$', fontsize=14)
axS.set_ylabel(r'$\sigma$ (mm)',fontsize = 14)  
axS.set_ylim(0.025,0.25)  
axS.yaxis.set_major_locator(MultipleLocator(0.05))
axS.yaxis.set_major_formatter('{x:.2f}')
axS.yaxis.set_minor_locator(MultipleLocator(0.01))
#axS.spines['left'].set_position(('data', 0))
#ax.spines['left'].set_zorder(10)
axS.legend()
#axS.grid()

#%%
# fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\Energy&Spin.png', dpi = 300)
# figS.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\Size_tof.png', dpi = 300)
