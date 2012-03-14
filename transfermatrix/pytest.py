from pylab import *  
from numpy import * 
import tmatrix
import math

Lx = 4
Ly = 7

def G (E, Phi, Wrap):
    tl, c = tmatrix.gettrans (Lx, Ly, E, Phi, Wrap)
    tvals = array(tl)
    #print tvals
    g = tmatrix.conductance(tl)
    print E, Phi, g, c
    return g

figure()
Evals = arange (-4.0, 4.01, 0.01)
Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/400); 

X, Y = meshgrid (Phivals, Evals)
Z = zeros (shape(X))
for i in range (len (Evals)):
    for j in range(len (Phivals)):
        Z[i, j] = G(Y[i, j], X[i, j], False)
        
pcolor (X/2.0/math.pi, Y, Z)
xlim(0.0, 1.0)
ylim(-4.0, 4.0)
title ('Open boundary')
xlabel (r'Magnetic flux $\Phi / \Phi_0$')
ylabel (r'Energy $\epsilon/t$')
colorbar()

figure()
Evals = arange (-4.0, 4.01, 0.01)
Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/400); 

X, Y = meshgrid (Phivals, Evals)
Z = zeros (shape(X))
for i in range (len (Evals)):
    for j in range(len (Phivals)):
        Z[i, j] = G(Y[i, j], X[i, j], True)
        
pcolor (X/2.0/math.pi, Y, Z)
xlim(0.0, 1.0)
ylim(-4.0, 4.0)
title ('Wrapped boundary')
xlabel (r'Magnetic flux $\Phi / \Phi_0$')
ylabel (r'Energy $\epsilon/t$')
colorbar()

show ()
        
    
