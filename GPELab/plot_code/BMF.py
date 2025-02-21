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
a12 = (-53.1907*a_bohr);

g11 = (4*np.pi*hbar**2*a11)/m 
g22 = (4*np.pi*hbar**2*a22)/m
g12 = (4*np.pi*hbar**2*a12)/m  

#theta = np.linspace(0,np.pi/2,50)
Omega_values = [706, 1410, 2810, 5600, 11200, 22300]
a_bmf = np.ones(len(Omega_values))

d = np.linspace(-8,8,200)


fig, ax = plt.subplots(2,1,constrained_layout=True, figsize=(8,5))
labels = [r'$\Omega$ = 706 Hz', r'$\Omega$ = 1.41 kHz', r'$\Omega$ = 2.81 kHz', r'$\Omega$ = 5.60 kHz', r'$\Omega$ = 11.2 kHz',r'$\Omega$ = 22.23 kHz']
colors = plt.get_cmap('Set2').colors
lw =2


for i, Om in enumerate(Omega_values):
    Om = Om * 2*np.pi
    delta = d*Om
    Om_tilde = np.sqrt(delta**2+Om**2)
    lOm = np.sqrt(hbar/(m*Om_tilde))
    P = delta/Om_tilde
    theta = np.arccos(np.sqrt((1+P)/2))
    a_mm = 2**0.5*(np.sin(theta)*np.cos(theta))**4/lOm * (a11+a22-2*a12)**2
    a_bmf[i] = np.max(a_mm/a_bohr)
    
    ax[0].plot(theta,a_mm/a_bohr, label=labels[i], lw=lw, color=colors[i])
    print(a_bmf[i])
    
ax[0].legend()
ax[0].set_ylabel(r'$\tilde{a}/a_0$',fontsize = 14)
    

#%% MF 
for i, Om in enumerate(Omega_values):
    Om = Om * 2*np.pi
    delta = d*Om
    Om_tilde = np.sqrt(delta**2+Om**2)
    P = delta/Om_tilde
    theta = np.arccos(np.sqrt((1+P)/2))

    a_minus = a11*(np.sin(theta))**4 + a22*(np.cos(theta))**4 + 0.5*a12*(np.sin(2*theta))**2;
    a12_mod = a12 + 2*a_bmf[i]*a_bohr
    a_minus_mod = a11*(np.sin(theta))**4 + a22*(np.cos(theta))**4 + 0.5*a12_mod*(np.sin(2*theta))**2;
    
    a = a_minus_mod -a_minus
    ax[1].plot(theta,a/a_bohr, label=labels[i], lw=lw, color=colors[i])
    print(np.max(a/a_bohr))
ax[1].set_ylabel('$a_{new} - a_{MF}$', fontsize=14)
ax[1].set_xlabel(r'$\theta$', fontsize=14)