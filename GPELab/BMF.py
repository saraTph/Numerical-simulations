# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 18:39:32 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt



u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

a_bohr = 0.52917721e-10
a11 = (86.4014*a_bohr); 
a22 = (33.2021*a_bohr);
a12 = ((-53.1907)*a_bohr);

g11 = (4*np.pi*hbar**2*a11)/m 
g22 = (4*np.pi*hbar**2*a22)/m
g12 = (4*np.pi*hbar**2*a12)/m  

theta = np.linspace(0,np.pi/2,50)
Omega_values = [706, 1410, 2810, 5600, 11200, 22300]
a_bmf = np.ones(len(Omega_values))


fig, ax = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))
labels = [r'$\Omega$ = 706 Hz', r'$\Omega$ = 1.41 kHz', r'$\Omega$ = 2.81 kHz', r'$\Omega$ = 5.60 kHz', r'$\Omega$ = 11.2 kHz',r'$\Omega$ = 22.23 kHz']
colors = plt.get_cmap('Set2').colors
lw =2


for i, Om in enumerate(Omega_values):
    Om = Om * 2*np.pi
    lOm = np.sqrt(hbar/m/Om)
    a_mm = 2**0.5*(np.sin(theta)*np.cos(theta))**4/lOm * (a11+a22-2*a12)**2
    ax.plot(theta,a_mm/a_bohr, label=labels[i], lw=lw, color=colors[i])
    a_bmf[i] = np.max(a_mm/a_bohr)
    print(np.max(a_mm/a_bohr))
    
ax.legend()
ax.set_xlabel(r'$\delta/\Omega$', fontsize=14)
ax.set_ylabel(r'$\tilde{a}/a_0$',fontsize = 14)
    

#%% MF 

fig2, ax2 = plt.subplots(1,1,constrained_layout=True, figsize=(9,4))
for i, Om in enumerate(Omega_values):
    Om = Om * 2*np.pi
    a_minus = a11*(np.sin(theta))**4 + a22*(np.cos(theta))**4 + 0.5*a12*(np.sin(2*theta))**2;
    a12 = a12 + 2*a_bmf[i]*a_bohr
    a_minus_mod = a11*(np.sin(theta))**4 + a22*(np.cos(theta))**4 + 0.5*a12*(np.sin(2*theta))**2;
    a = a_minus-a_minus_mod
    ax2.plot(theta,a/a_bohr, label=labels[i], lw=lw, color=colors[i])