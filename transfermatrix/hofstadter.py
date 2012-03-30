from pylab import *  
from numpy import * 
import tmatrix
import math

Lx = 4
Ly = 7

def SaveData (Xy, Yu, Zu, Xw, Yw, Zw):
    import shelve
    f = shelve.open('hofstadter-%dx%d.dat' % (Lx, Ly))
    f['Xw'] = Xw
    f['Yw'] = Yw
    f['Zw'] = Zw
    f['Xu'] = Xu
    f['Yu'] = Yu
    f['Zu'] = Zu
    f.close()
    

def CompareOldData (Xu, Yu, Zu, Xw, Yw, Zw):
    try: 
        import shelve
        f = shelve.open('hofstadter-%dx%d.dat' % (Lx, Ly))
        Xwold = f['Xw']
        Ywold = f['Yw']
        Zwold = f['Zw']
        Xuold = f['Xu']
        Yuold = f['Yu']
        Zuold = f['Zu']
        f.close()
        dxu = norm (Xu - Xuold)
        dyu = norm (Yu - Yuold)
        dzu = norm (Zu - Zuold)
        print "dxu: ", dxu, "dyu:", dyu, "dzw", dzu
        if (dxu >  1e-6) or (dyu > 1e-6):
            print "Grid has changed, comparison is meaningless"
            return None
    
        figure() 
        title ('Diff for unwrapped case')
        pcolor (Xu, Yu, Zu - Zuold)
        colorbar()
    
        dxw = norm (Xw - Xwold)
        dyw = norm (Yw - Ywold)
        dzw = norm (Zw - Zwold)
        print "dxw: ", dxw, "dyw:", dyw, "dzw:", dzw
        if (dxw >  1e-6) or (dyw > 1e-6):
            print "Grid has changed, comparison is meaningless"
            return None
        figure() 
        title ('Diff for wrapped case')
        pcolor (Xw, Yw, Zw - Zwold)
        colorbar()
        return dzu + dzw
    except:
        pass
    print "Cannot read saved data to do comparison"
    return None
    
    
def G (E, Phi, Wrap):
    tl, c = tmatrix.gettrans ('Y', 'Y', Lx, Ly, E, Phi, Wrap)
    tvals = array(tl)
    #print tvals
    g = tmatrix.conductance(tl)
    if c > 1e-8: 
        import sys
        print E, Phi, g, c, wrap
        print "unitarity does not hold, exiting"
        sys.exit(-1)
    return g

figure()
Evals = arange (-4.0, 4.01, 0.01)
Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/100); 

Xu, Yu = meshgrid (Phivals, Evals)
Zu = zeros (shape(Xu))
for i in range (len (Evals)):
    print "NOWRAP: E = ", Yu[i, 0]
    for j in range(len (Phivals)):
        Zu[i, j] = G(Yu[i, j], Xu[i, j], False)
        
pcolor (Xu/2.0/math.pi, Yu, Zu)
xlim(0.0, 1.0)
ylim(-4.0, 4.0)
title ('Open boundary')
xlabel (r'Magnetic flux $\Phi / \Phi_0$')
ylabel (r'Energy $\epsilon/t$')
colorbar()

figure()
Evals = arange (-4.0, 4.01, 0.01)
Phivals = arange (0.0, 2.0*math.pi + 0.0001, 2.0*math.pi/100); 

Xw, Yw = meshgrid (Phivals, Evals)
Zw = zeros (shape(Xw))
for i in range (len (Evals)):
    print "WRAP: E = ", Yw[i, 0]
    for j in range(len (Phivals)):
        Zw[i, j] = G(Yw[i, j], Xw[i, j], True)
        
pcolor (Xw/2.0/math.pi, Yw, Zw)
xlim(0.0, 1.0)
ylim(-4.0, 4.0)
title ('Wrapped boundary')
xlabel (r'Magnetic flux $\Phi / \Phi_0$')
ylabel (r'Energy $\epsilon/t$')
colorbar()

if CompareOldData (Xu, Yu, Zu, Xw, Yw, Zw) == None: 
    print "Saving the data"
    SaveData (Xu, Yu, Zu, Xw, Yw, Zw)
    
show ()
        
    
