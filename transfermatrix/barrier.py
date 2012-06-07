from pylab import *  
from numpy import * 
import tmatrix
import math

Lx = 30
Ly = 60
W = 10

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
            if not f.has_key(k): 
               print "Saved ", k, "is not available, skip"
               continue
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
    
    
def T (E, U0, Delta):
    print E, U0, Delta
    V = zeros((Lx, Ly))
    Phi = 0.0
    wrap = True
    gauge = 'X'
    current = 'Y'
    V[0::2, 0::2] = -Delta
    V[1::2, 1::2] = -Delta
    V[0::2, 1::2] = Delta
    V[1::2, 0::2] = Delta
    y1 = Ly/2 - W/2
    y2 = Ly/2 + W/2
    V[:, y1:y2] += U0
  
    tl, c = tmatrix.gettrans (current, gauge, Lx, Ly, E + 1e-5, Phi, wrap, V)
    tvals = array(tl)
    g = tmatrix.conductance(tl)
    if c > 1e-8: 
        import sys
        print E, Phi, g, c, wrap
        print "unitarity does not hold, exiting"
        sys.exit(-1)
    return tl


def doScan(Evals, Uvals, Delta): 
    X, Y = meshgrid (Uvals, Evals)
    Z = zeros (shape(X))
    for i in range (len (Evals)):
        #print 'Wrap:', wrap, 'Dir:', dir, 'Gauge:', gauge, "E = ", Y[i, 0]
        print "E = ", Evals[i], "Delta = ", Delta
        for j in range(len (Uvals)):
            print Y[i, j], X[i, j], Delta
            Z[i, j] = T(Y[i, j], X[i, j], Delta)[0]
    figure()     
    pcolor (X, Y, Z)
    tit = 'Delta = %g W = %g' % (Delta, W)
    #xlim(0.0, 1.0)
    #ylim(-4.0, 4.0)
    #tit = ''
    #if wrap: 
    #   tit += 'Wrapped boundary, '
    #else:
    #   tit += 'Open boundary, ' 
    #tit += "current J || %s, " % dir
    #tit += 'potential A || %s ' % gauge
    title (tit)
    xlabel (r'Barrier height $U_0$')
    ylabel (r'Energy $\epsilon/t$')
    colorbar()
    #show ()
    return X, Y, Z

def doPlot (E, Delta, Uvals):
    figure()
    Tvals = [] 
    for U in Uvals: 
        Tvals.append (T(E, U, Delta))
    for i in [0, 1, 2, 3, 4, 5]:
       plot (Uvals, [t[i] for t in Tvals])
    title ('E = %g Delta = %g W = %g' % (E, Delta, W))


Evals = arange (-0.3, 0.3, 0.002)
Uvals = arange (-0.5, 0.5, 0.005)
#Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/100); 
#doScan(Evals, Uvals, 0.0)
#doScan(Evals, Uvals, 0.05)
#doScan(Evals, Uvals, 0.1)
#doScan(Evals, Uvals, 0.2)
#doScan(Evals, Uvals, 0.5)
#W = 10
#doPlot (-0.25, 0.0, Uvals)
#W = 2
#doPlot (-0.25, 0.0, Uvals)
#W = 1
#doPlot (-0.25, 0.0, Uvals)
#W = 10
#doPlot (0.11, 0.1, Uvals)
#W = 10
#doPlot (0.2, 0.1, Uvals)
W = 20
doPlot (0.033, 0.03, Uvals)
show ()

        
    
