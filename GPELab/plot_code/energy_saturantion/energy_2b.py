# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 17:48:59 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s

a_bohr = 0.52917721e-10
a22 = (86.4014*a_bohr) 
a11 = (33.2755*a_bohr)
a12 = (-53.1022*a_bohr)

g11 = 4*hbar**2*np.pi*a11/4 
g22 = 4*hbar**2*np.pi*a22/4 
g12 = 4*hbar**2*np.pi*a12/4 
g_0 =  4*hbar**2*np.pi*a_bohr/4 

#%%

#1st notation
g_bar = (g11+g22-2*g12)/4
g_inf = (g11*g22-g12**2)/(4*g_bar)
k = (g11-g22)/(4*g_bar)

d = np.linspace(-4,4,100)
#cotg = 1/tg
#cotg(d) = theta
phi = np.linspace(0,np.pi,100)
g2 = g_inf+g_bar*(np.cos(phi)-k)**2

# 2nd notation
g = (g11+g22)
g_star = (g11-g12)/2
dg = g11-g22

#g2_new = g-g_star/(1+d**2)-dg/2/(1+d**2)**(1/2)

fig, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(4, 4))
ax.plot(d,g2/g_0)
#ax.plot(d,g2_new/g_0,'r')

