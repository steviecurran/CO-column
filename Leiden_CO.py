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
c = 35.28

plt.rcParams.update({'font.size': 16})
ax = plt.gca();
plt.setp(ax.spines.values(), linewidth=2)
ax.tick_params(direction='in', length=6, width=1.5, which='major')
ax.tick_params(direction='in', length=3, width=1.5, which='minor')
ax.tick_params(axis='both', which='major', pad=5)

x = df['UP']; y = df['E_u(K)']
ax.plot(x,y,'o', markersize=5, c = 'k')
fit = np.polyfit(x,y,2)
print(fit)
x1,x2 = ax.get_xlim(); y1,y2 = ax.get_ylim()
xplot = np.linspace(x1,x2,100)
yplot = fit[2]+fit[1]*xplot + fit[0]*xplot**2
plt.plot(xplot,yplot, c='red', linestyle ="dashed", lw=2.0)

plt.xlabel(r"Upper $J$ number of transition"); plt.ylabel("$E_u$ [K]")
plt.tight_layout()
out = 'Leiden_CO_E-T=%1.0f_c=%1.2f'  %(T,c)
plt.savefig("%s.eps" %(out), format='eps'); print('Written to %s' %(out))
plt.savefig("%s.png" %(out), format='png'); print('Written to %s' %(out))
plt.show()
plt.close()

############ SUMMING ##############
df['N'] = (2*df['UP']+1)*np.exp(-(df['E_u(K)']/T) +c);
df['N_total'] = round(df.N.cumsum(),3)

df['temp'] = df['E_u(K)'].shift(1) # PREVIOUS VALUE
df['temp2'] = h*df['FREQ(GHz)']*1e9 + df['temp']*k
df['E_calc'] = round(df['temp2']/k,3) # THAT'S IT!

del df['temp']; del df['temp2']

print(df)
plt.rcParams.update({'font.size': 16})
ax = plt.gca();
plt.setp(ax.spines.values(), linewidth=2)
ax.tick_params(direction='in', length=6, width=1.5, which='major')
ax.tick_params(direction='in', length=3, width=1.5, which='minor')
ax.tick_params(axis='both', which='major', pad=5)
x = df['UP']; y = df['N_total']
ax.plot(x,y,'o', markersize=5, c = 'k')

plt.xlabel(r"Upper $J$ number of transition"); plt.ylabel("$\sum N_J$ [cm$^{-2}$]")
plt.tight_layout()
out = 'Leiden_CO_N-T=%1.0f_c=%1.2f'  %(T,c)
plt.savefig("%s.eps" %(out), format='eps'); print('Written to %s' %(out))
plt.savefig("%s.png" %(out), format='png'); print('Written to %s' %(out))
plt.show()
plt.close()

#########  GET TOTAL N FROM FIT ############
df['E'] =  fit[0]*df['UP']**2 + fit[1]*df['UP'] + fit[2];
#print(df) # VALUES CLOSE, BUT BETTER WITH AN ARRAY
arr = []
stop = 10000
sum = 0
for J in range(1,stop+1):
    E = fit[0]*J**2 + fit[1]*J + fit[2];
    N = (2*J+1)*np.exp(-1*E/T +c)
    #print(N)
    sum = sum + N

print("Over first %d levels N = %1.3e" %(J,sum))
