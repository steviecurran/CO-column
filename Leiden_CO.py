#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import os 
import sys
from shutil import get_terminal_size
pd.set_option('display.width', get_terminal_size()[0]) 
pd.set_option('display.max_columns', None)

c = 299792458        # m/s
h = 6.62606957e-34      # J s
k = 1.3806488e-23       # J/K

df = pd.read_csv("Leiden_CO.dat", delim_whitespace=True); 

T = 568
C = 35.28

plt.rcParams.update({'font.size': 16})
ax = plt.gca();
plt.setp(ax.spines.values(), linewidth=2)
ax.tick_params(direction='in', length=6, width=1.5, which='major')
ax.tick_params(direction='in', length=3, width=1.5, which='minor')
ax.tick_params(axis='both', which='major', pad=5)

# IF USING DATA FROM https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat

df['N'] = (2*df['UP']+1)*np.exp(-(df['E_u(K)']/T) +C); #print(df)
df['N_total'] = round(df.N.cumsum(),3)

# OR COULD CALCULATE HERE 
stop = 41 # AS FAR AS LEIDEN GOES
freq = 115.2712018e9 # CO 1 -> 0
E_tot = 0
N_tot =0
data = []

for J in range(1,stop+1):
    E_J = J*h*freq/k
    E_tot = E_tot+E_J
    g = 2*J + 1
    N_J = g*np.exp(-E_tot/T + C)
    N_tot = N_tot + N_J
    #print(J, E_tot, N_J, N_tot)
    data.append(J)
    data.append(E_tot)
    data.append(N_tot)

data = np.reshape(data,(-1,3)); #print(data)

#############################################################
plt.rcParams.update({'font.size': 16})
ax = plt.gca();
plt.setp(ax.spines.values(), linewidth=2)
ax.tick_params(direction='in', length=6, width=1.5, which='major')
ax.tick_params(direction='in', length=3, width=1.5, which='minor')
ax.tick_params(axis='both', which='major', pad=5)

y = df['N_total']; x = df['UP']
print("Using Leiden up to J = 41,  N_total = %1.3e cm-2" %(df['N_total'].iloc[-1])) 

#x = data[:,1]; y = data[:,2]
print("Using my own summing up to J = %d, N_total = %1.3e cm-2" %(stop+1,data[-1,2])) 

ax.plot(x,y,'o', markersize=5, c = 'k')
plt.xlabel(r"Upper $J$ number of transition"); plt.ylabel("$\sum N_J$ [$10^{17}$ cm$^{-2}$]")

### SECONDARY AXES #####
from matplotlib.ticker import AutoMinorLocator
# TRICKY, BUT Leiden_CO-poly_fit.py  GIVES A DECENT 2ND ORDER POLY, SO TRY
A = 2.73975855; B = -2.13791121; CEE =  -2.93187467
J_test = 30; E_test = A*J_test**2 + B*J_test + CEE
print("J_test = %1.0f => E_total = %1.0f K" %(J_test, E_test)) # CHECK - SPOT ON

# def forward(x):
#     return A*(x**2) # A*x**2 + B*x + CEE  
# def inverse(x):
#     return A*(x**2)# A*x**2 + B*x + CEE  # NOT WORKING - BECAUSE NOT 1-TO-1 ?

x1,x2 = ax.get_xlim(); y1,y2 = ax.get_ylim()
ax.yaxis.get_offset_text().set_visible(False) # HIDING THE 1 E 17

xold = []
xnew = []
inc = 10

for j in range(int(x1*inc),int(x2*inc)): 
    J = float(j)#/inc
    E = A*J**2 + B*J + CEE 
    xold.append(j)
    xnew.append(E)

def forward(x):
    return np.interp(x, xold, xnew)

def inverse(x):
    return np.interp(x, xnew, xold)

secax = ax.secondary_xaxis('top', functions=(forward, inverse))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.tick_params(direction='in', length=6, width=1.5, which='major')
secax.tick_params(direction='in', length=3, width=1.5, which='minor')
secax.set_xlabel(r"$E_u$ [K]",labelpad=5)

plt.tight_layout()
out = 'Leiden_CO_N-T=%1.0f_C=%1.2f'  %(T,C)
plt.savefig("%s.eps" %(out), format='eps'); print('Written to %s' %(out))
plt.savefig("%s.png" %(out), format='png'); print('Written to %s' %(out))
plt.show()
plt.close()


