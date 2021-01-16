# -*- coding: utf-8 -*-
"""
Created on Mon Aug 03 15:14:17 2015

@author: Dustin
"""

from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

v1 = 0.7
K1 = 1
n = 4
v2 = 0.35
K2 = 1
k3 = 0.7
v4 = 0.35
K4 = 1
k5 = 0.7
v6 = 0.35
K6 = 1
k7 = 0.35
v8 = 1
K8 = 1
vc = 0.4
Kc = 1
K = 0.5
L = 0

# solve the system dy/dt = f(y, t)
def f(y,t):
        Xi = y[0]
        Yi = y[1]
        Zi = y[2]
        Vi = y[3]
        F = y[4]
        # the model equations (see Munz et al. 2009)
        f0 = v1 * (K1 / (K1 + Zi)) - v2*(Xi/(K2+Xi)) + vc* ((K*F)/Kc + (K*F)) + L
        f1 = k3*Xi - v4*((Yi)/(K4+Yi))
        f2 = k5*Yi - v6*((Zi)/(K6+Zi))
        f3 = k7*Xi - v8*((Vi)/(K8+Vi))
        f4 = Vi 
        return [f0, f1, f2, f3, f4]
        
# initial conditions
Xi = 0       
Yi = 0        
Zi = 0
Vi = 0
F = 0   
y0 = [Xi, Yi, Zi, Vi, F]  
t  = np.linspace(0, 3000.)
        
# solve the DEs
soln = integrate.odeint(f, y0, t)
Xi = soln[:, 0]
Yi = soln[:, 1]
Zi = soln[:, 2]
Vi = soln[:, 3]
F = soln[:, 4]

plt.figure()
plt.plot(t, Xi, label='mRNA')
plt.plot(t, Yi, label='Protein')
plt.plot(t, Zi, label='Inhibitor')
plt.plot(t, Vi, label='Neuropeptide')
plt.plot(t, F, label='Mean Field')
plt.legend(loc=0)