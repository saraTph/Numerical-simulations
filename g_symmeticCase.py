# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 12:17:10 2024

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt

# Define theta values
theta = np.linspace(0, np.pi / 2, 500)

# Define g values
u = 1.66053906660e-27
m = 39*u
hbar = 1.054571818e-34 #J s
a0 = 0.52917721e-10

g_up_up = 4*np.pi*hbar**2*(51.4809*a0)/m
g_down_down = 4*np.pi*hbar**2*(35.6377*a0)/m
g_up_down = 4*np.pi*hbar**2*(-53.7389*a0)/m
g0 = 4*np.pi*hbar**2*(a0)/m

# Define the expressions
g1 = g_down_down * np.cos(theta)**4 + g_up_up * np.sin(theta)**4 + 0.5 * g_up_down * np.sin(2*theta)**2
g2 = g_down_down * np.sin(theta)**4 + g_up_up * np.cos(theta)**4 + 0.5 * g_up_down * np.sin(2*theta)**2
g3 = (g_down_down + g_up_up) * np.sin(2*theta)**2 + 2 * g_up_down * np.cos(2*theta)**2
g4 = 1/4 * (g_down_down + g_up_up - 2 * g_up_down) * np.sin(2*theta)**2
g5 = np.sin(2*theta) * (g_up_up * np.sin(theta)**2 - g_down_down * np.cos(theta)**2 + g_up_down * np.cos(2*theta))
g6 = np.sin(2*theta) * (g_up_up * np.cos(theta)**2 - g_down_down * np.sin(theta)**2 - g_up_down * np.cos(2*theta))

Om=1
V=1
g2b = (2*g2/V + g4**2*np.sqrt(hbar*Om)/(V*np.sqrt(2)*np.pi))*1/4
#%%
lw = 1
fs = 14
colors = plt.get_cmap('Paired').colors

# Plotting the expressions
fig, ax = plt.subplots(1,1,constrained_layout=True)

plt.axhline(y=0, color='lightgrey', linestyle='-', linewidth=1, zorder=0)
ax.plot(theta, g2/g0, label=r'$g_2$', lw=lw, color=colors[5], linestyle='solid', zorder=0)
ax.plot(theta, g3/g0, label=r'$g_3$', lw=lw, color =colors[1], linestyle='dotted', zorder=0)
ax.plot(theta, g4/g0, label=r'$g_4$', lw=lw, color = colors[3], linestyle='-.', zorder=0)
ax.plot(theta, g6/g0, label=r'$g_6$', lw=lw, color =colors[7], linestyle='dashed', zorder=0)


theta_ticks = [0, np.pi / 4, np.pi / 2]  # Corresponding theta positions
xticks = [0, np.pi / 4, np.pi / 2]
xtick_labels = [r'$0$', r'$\pi/4$', r'$\pi/2$']
ax.set_xticks(xticks)
ax.set_xticklabels(xtick_labels, fontsize=fs)

ax.set_xlabel(r'$\theta$', fontsize=fs)
ax.set_ylabel(r'$g_i/g_0$', fontsize=fs)

#%% Save figure
#fig.savefig('g_i_color_dot.pdf', dpi=600)