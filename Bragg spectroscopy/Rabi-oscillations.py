# -*- coding: utf-8 -*-
"""
Created on Tue Apr 15 17:27:50 2025

@author: sarat
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters
Omega = 2 * np.pi * 500   # Coupling strength (Hz), related to lattice depth
Delta = 2 * np.pi * 0     # Energy detuning (Hz), 0 at BZ edge
t_max = 5e-3              # Total simulation time (s)
n_steps = 1000

# Time array
t = np.linspace(0, t_max, n_steps)

# Two-level system ODE
def band_dynamics(t, y):
    c0 = y[0] + 1j * y[1]
    c1 = y[2] + 1j * y[3]
    dc0dt = -1j * (Omega / 2) * c1
    dc1dt = -1j * (Omega / 2) * c0 - 1j * Delta * c1
    return [dc0dt.real, dc0dt.imag, dc1dt.real, dc1dt.imag]

# Initial state: all population in ground state
y0 = [1.0, 0.0, 0.0, 0.0]

# Solve the ODE
sol = solve_ivp(band_dynamics, [0, t_max], y0, t_eval=t)

# Extract populations
c0 = sol.y[0] + 1j * sol.y[1]
c1 = sol.y[2] + 1j * sol.y[3]
P0 = np.abs(c0)**2
P1 = np.abs(c1)**2

# Plot
plt.figure(figsize=(8, 5))
plt.plot(t * 1e3, P0, label='Ground band')
plt.plot(t * 1e3, P1, label='Excited band')
plt.xlabel("Time (ms)")
plt.ylabel("Population")
plt.title("Coherent Oscillations Between Bands at BZ Edge")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()
