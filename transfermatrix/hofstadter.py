from pylab import *  
from numpy import * 
import tmatrix
import math

Lx = 4
Ly = 8

def SaveData (data):
    import shelve
    f = shelve.open('hofstadter-%dx%d.dat' % (Lx, Ly))
    for k, v in data.items():
        f[k] = v
    f.close()
    

def CompareOldData (dict):
    s = 0.0
    try: 
        import shelve
        f = shelve.open('hofstadter-%dx%d.dat' % (Lx, Ly))
        for k, v in data.items(): 
            fv = f[k]
            d  = norm (v - fv)
            print "diff: ", k, d
            s += d
        f.close()
        return s
    except:
        pass
    print "Cannot read saved data to do comparison"
    return None
    
    
def G (E, Phi, current, gauge, wrap):
    tl, c = tmatrix.gettrans (current, gauge, Lx, Ly, E, Phi, wrap)
    tvals = array(tl)
    g = tmatrix.conductance(tl)
    if c > 1e-8: 
        import sys
        print E, Phi, g, c, wrap
        print "unitarity does not hold, exiting"
        sys.exit(-1)
    return g


def doScan(Evals, Phivals, wrap, dir, gauge): 
    X, Y = meshgrid (Phivals, Evals)
    Z = zeros (shape(X))
    for i in range (len (Evals)):
        print 'Wrap:', wrap, 'Dir:', dir, 'Gauge:', gauge, "E = ", Y[i, 0]
        for j in range(len (Phivals)):
            Z[i, j] = G(Y[i, j], X[i, j], dir, gauge, wrap)
    figure()     
    pcolor (X/2.0/math.pi, Y, Z)
    xlim(0.0, 1.0)
    ylim(-4.0, 4.0)
    tit = ''
    if wrap: 
       tit += 'Wrapped boundary, '
    else:
       tit += 'Open boundary, ' 
    tit += "current J || %s, " % dir
    tit += 'potential A || %s ' % gauge
    title (tit)
    xlabel (r'Magnetic flux $\Phi / \Phi_0$')
    ylabel (r'Energy $\epsilon/t$')
    colorbar()
    return X, Y, Z

def scanE (Evals):
    X = Evals
    Yx = zeros (shape(X))
    Yy = zeros (shape(X))
    for i, E in enumerate (X):
        Yx[i] = G (X[i], 0.0, 'X', 'X', False)
        Yy[i] = G (X[i], 0.0, 'Y', 'X', False)
    figure()
    plot (X, Yx, label='J || x')
    plot (X, Yy, label='J || y')
    legend()
    title ('Open boundary, zero flux')
    xlabel (r'Energy $\epsilon/t$')
    ylabel (r'Conductance $G(E)$')  
    return X, Yx

Evals = arange (-4.0, 4.01, 0.01)
Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/100); 

Xx, Yx        = scanE(Evals)
Xuy, Yuy, Zuy = doScan(Evals, Phivals, False, 'Y', 'Y')
Xux, Yux, Zux = doScan(Evals, Phivals, False, 'Y', 'X')
Xwy, Ywy, Zwy = doScan(Evals, Phivals, True,  'Y', 'Y')
Xwx, Ywx, Zwx = doScan(Evals, Phivals, True,  'Y', 'X')

#if CompareOldData (Xu, Yu, Zu, Xw, Yw, Zw) == None:
data = dict (Xux=Xux, Yux=Yux, Zux=Zux,
             Xuy=Xuy, Yuy=Yuy, Zuy=Zuy, 
             Xwx=Xwx, Ywx=Ywx, Zwx=Zwx, 
             Xwy=Xwy, Ywy=Ywy, Zwy=Zwx, 
             Xx=Xx, Yx=Yx)
             
if CompareOldData (data) == None: 
    print "Saving the data"
    SaveData (data)
#    SaveData (Xu, Yu, Zu, Xw, Yw, Zw)
print "Check gauge invariance: unwrapped, J||Y: ", norm(Zuy - Zux)             
    
show ()
        
    
