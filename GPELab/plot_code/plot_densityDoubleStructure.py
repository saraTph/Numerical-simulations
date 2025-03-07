# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 18:26:24 2025

@author: Sarah
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import os

#%% Import

Omega_values = [22300, 11200, 5600, 2810, 1410, 706]

#initialize simulation results
len_sim = 129

matrix1 = np.ones((4, len_sim)) 
matrix2 = np.ones((4, len_sim)) 

densities1 = np.ones((6, len_sim)) #matrix for all the Rabi 
densities2 = np.ones((6, len_sim)) #matrix for all the Rabi

mat = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\doubleStructure\densities_5n0.mat')
#mat = scipy.io.loadmat(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs\doubleStructure\output_data_5n0.mat')
for i in range(6):
    # I do perform the mean over the densities 
    matrix1[0,:] = mat.get('d0_1')[i,:]  #density 0
    matrix1[1,:] = mat.get('d1_1')[i,:]  #density 1
    matrix1[2,:] = mat.get('d2_1')[i,:]  #density 2
    matrix1[3,:] = mat.get('d3_1')[i,:]  #density 3
    
    densities1[i,:] = np.mean(matrix1, axis=0)
    
for i in range(6):
    # I do perform the mean over the densities 
    matrix2[0,:] = mat.get('d0_2')[i,:]  #density 0
    matrix2[1,:] = mat.get('d1_2')[i,:]  #density 1
    matrix2[2,:] = mat.get('d2_2')[i,:]  #density 2
    matrix2[3,:] = mat.get('d3_2')[i,:]  #density 3
    
    densities2[i,:] = np.mean(matrix2, axis=0)
    
#     profile1 = np.mean()
# profile2 = mat.get('Density1_1D_2')


#%% size after tof
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

colors = plt.get_cmap('Set2').colors
lw = 1.7
titles = [r'$\Omega$ = 22.3 kHz', r'$\Omega$ = 11.2 kHz', r'$\Omega$ = 5.6 kHz', r'$\Omega$ = 2.81 kHz', r'$\Omega$ = 1.41 kHz', r'$\Omega$ = 706 Hz']
labels = [r'$\Psi_1$', r'$\Psi_2$']

figS, axS = plt.subplots(1,6,constrained_layout=True, figsize=(17,4))

#for i in range(np.size(Omega_values)):
for i in range(6):
    axS[i].plot(densities1[i,:], label=labels[0], lw=1.7, color=colors[0],zorder =1)
    axS[i].plot(densities2[i,:], label=labels[1], lw=1.7, color=colors[1],zorder =1)
    axS[i].set_title(titles[i])
    axS[i].legend()
   
   

#%%
# fig.savefig(r'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\25-02-17\Energy&Spin.png', dpi = 300)
# figS.savefig(r'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\plot_code\Figures\Double structure
#              .png', dpi = 300)
