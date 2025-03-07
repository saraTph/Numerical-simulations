import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#gaussian function to fit the TF column density profile
# we consider only the x direction (y = 0)
def Gauss(x,n0,s_x,s_z):  #param = n0,s_z,s_z   
    return n0*np.exp(-0.5*(x/s_x)**2)*np.sqrt(2*np.pi*s_z**2) 

def TF(x,n0,R0):
    return n0/R0**2*4/3*(R0**2-(x**2))**(3/2)

def truncate(n, decimals=2):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

TF_colDen = np.zeros(1000)
x_axis = np.linspace(-3,3,1000)
# difine my TF column density vector
R0 = 1
n0 = 1

idx = 0
for i in x_axis:
    if (i>=-R0 and i<=R0):
        TF_colDen[idx] = TF(i,n0,R0)
    else:
        TF_colDen[idx] = 0
    idx = idx + 1

fig, ax = plt.subplots()
param0 = [1,1,1]
PARAM = np.zeros(3)
plt.plot(x_axis,TF_colDen,'-',color = "blue",label='TF column density')
param, cov = curve_fit(Gauss, x_axis, TF_colDen, p0=param0)
PARAM[0], PARAM[1], PARAM[2] = param
plt.plot(x_axis, Gauss(x_axis,PARAM[0], PARAM[1], PARAM[2]) , '-', label='Gauss fit',color = "red") 
PARAM[1] = truncate(PARAM[1],3)
ratio = PARAM[1]*100/R0
ratio =truncate(ratio,3)
textstr = '\n'.join((
    r'TF radius= ' +str(R0, ),
    r'gauss s_x= ' +str(PARAM[1], ),
    r'size ratio= ' + str(ratio,)+'%' ))


#textstr = '\n'.join((r'Gauss size vs TF radius=' +str(ratio)))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.95, 0.95, textstr,fontsize=12, bbox=props)
plt.legend()
plt.savefig('TFvsGauss.png', dpi=300)
plt.show()

print('Gauss size: Sigma=',PARAM[1])
print('TF Radius: R0=', R0)
print('size difference', ratio,'%')