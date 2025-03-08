from math import pi
import numpy as np
import matplotlib.pyplot as plt

# s'' = (1-2*n1D)*s^(-3)-s


#%% define functions

def scatteringLengthDressed(N,delta,a11,a22,a12):
    for i in range(N):
        theta = pi/2-np.arctan(delta)
        return a11*(np.sin(theta/2))**4 + a22*(np.cos(theta/2))**4 + 0.5*a12*(np.sin(2*theta/2))**2

def sweepLin(tpoints, N, di, df): #return the value of delta/Omega at time time
    # output: delta scan
    delta = np.zeros(N)
    Om = 1
    m=(df-di)/tpoints[-1]
    for j in range(N):
        t = tpoints[j]
        delta[j] = (di + m*t)/Om
    return delta


def FLin(u, v, t, n1D, a_s):
    if (t==0):
        print('Initial as =', a_s/a_bohr)
    return 1 * (1+2*n1D*a_s)/u**3 - u

#%% define variables 

hbar=6.62607015*10**(-34)/(2*np.pi)
u = 1.661*10**(-27)                      # [Kg] atomic mass units
m = 38.96370668*u                        # K39 mass
a_bohr = 0.52917721e-10;                 # [m] atomic Bohr radius
a11 = (86.4014*a_bohr); 
a22 = (33.2755*a_bohr);
a12 = (-53.1022*a_bohr);

#Om = 1                                  # omega Rabi (is set to 1 because do not matter (we care only about delta/Omega)
time_i = 0           
time_f = 9e-3                            # [s]
n1D = 1e9                                # [1/m] 1D density of BEC
wr = 169*(2*pi)                          # [Hz] radial frequency harmonic trap 
size  = np.sqrt(hbar/(m*wr))             # typical length harmonic trap


# adimensionalization 
T = 1/wr
L = np.sqrt(hbar/(m*wr))

# define sweep and dressed scattering lenght

N = 100                                  # number of steps
tf = time_f/T
ti = time_i/T
h = (tf-ti)/N
tpoints = np.arange(ti,tf,h)             # adimensional time vector

di = 8                                   # initial delta/Omega
df = -8                                  # final delta/Omega
delta = sweepLin(tpoints, N, di, df)
a_dressed = scatteringLengthDressed(N,delta,a11,a22,a12)
    

fig, ax = plt.subplots(2,1,constrained_layout=True, figsize=(10,8))
ax[0].plot(tpoints*T, delta, color = 'lightcoral',linewidth=3)
ax[0].set_xlabel('$t (s)$',size=20); 
ax[0].set_ylabel('$\delta/\Omega$',size=20)


ax[1].plot(tpoints*T , a_dressed/a_bohr ,color = 'lightcoral',linewidth=3)
ax[1].set_xlabel('$t (s)$',size=20); 
ax[1].set_ylabel('$a_{--}/a_0$',size=20)
fig.subplots_adjust(hspace=0.3)


#%% Size dynalics during the RD sweep (RK4 method)

u_eq11 = (hbar/(m*wr))**0.5*(1+2*n1D*a11)**(1/4)  #equilibrium solution
print('Equilibrium size a11= ',u_eq11)

u_eq22 = (hbar/(m*wr))**0.5*(1+2*n1D*a22)**(1/4)  #equilibrium solution
print('Equilibrium size a22= ',u_eq22)

print('')
print('Simulation trap frequency = 200 Hz begins...')
print('a_s changes during the linear sweep and equilibrium CI')
print('Trap frequency = ', wr, '(Hz)')
print('Trap length = ', size, '(m)')


u_sol = np.zeros(len(tpoints))
v_sol = np.zeros(len(tpoints))

#theta_i = pi/2-np.arctan(di)
u_eq = size*(1+2*n1D*a_dressed[0])**(1/4)      #equilibrium solution
u = (u_eq)/L
v = 0
print('initial conditions:')
print('u = ',u,' --> ',u*L,'[m]')
print('v = ',v,' --> ',v*(T/L),'[m/s]')
for i, t in enumerate(tpoints):
    u_sol[i] = u
    v_sol[i] = v
    a_s = a_dressed[i]

    m1 = h*v
    k1 = h*FLin(u, v, t, n1D, a_s)  #(s, s', t)

    m2 = h*(v + 0.5*k1)
    k2 = h*FLin(u+0.5*m1, v+0.5*k1, t+0.5*h, n1D, a_s)

    m3 = h*(v + 0.5*k2)
    k3 = h*FLin(u+0.5*m2, v+0.5*k2, t+0.5*h,n1D, a_s)

    m4 = h*(v + k3)
    k4 = h*FLin(u+m3, v+k3, t+h, n1D, a_s)

    u += (m1 + 2*m2 + 2*m3 + m4)/6
    v += (k1 + 2*k2 + 2*k3 + k4)/6

print('Size variation during the sweep= ', ((u-u_eq/L)*100)/(u_eq/L),'%')


en = hbar**2/(2*m* (u_sol*L)**2) *(1+2*n1D*a_dressed)+0.5*m*wr**2*(u_sol*L)**2

#%% plots 

figDyn, axDyn = plt.subplots(3,1,constrained_layout=True, figsize=(7,8))

axDyn[0].plot(delta, u_sol *L/size, label='$\sigma_{in}=\sigma_{eq}$: $f_{HO} = 200 Hz$, $a_s = a_{11}$')
#axDyn[0].set_xlabel('$t\omega_T$',size=15)
axDyn[0].set_xlabel('$\delta/\Omega$',size=15)
axDyn[0].set_ylabel('$\sigma/\sigma_{HO}$',size=15)

axDyn[1].plot(delta, u_sol*L * 1e6, label='$\sigma_{in}=\sigma_{eq}$: $f_{HO} = 200 Hz$, $a_s = a_{11}$')
#axDyn[1].set_xlabel('$t (ms)$',size=15)
axDyn[1].set_xlabel('$\delta/\Omega$',size=15)
axDyn[1].set_ylabel('$\sigma (\mu m)$',size=15)

axDyn[2].plot(delta, en/(hbar*wr), label='$\sigma_{in}=\sigma_{eq}$: $f_{HO} = 200 Hz$, Linear sweep')
#axDyn[2].set_xlabel('$t (ms)$',size=15)
axDyn[2].set_xlabel('$\delta/\Omega$',size=15)
axDyn[2].set_ylabel('$E/(\hbar\Omega_{r})$',size=15)