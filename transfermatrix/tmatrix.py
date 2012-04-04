import _tmatrix
import numpy

def gettrans(current, gauge, Lx, Ly, E, phi, wrap):
    """
       This function calculates transmission
       eigenvalues
       
       Arguments:
           
           current -- current direction
           
           gauge   -- vector potential direction
           
           Lx   --- the size of the sample in x direction
                    (across the current flow)
           
           Ly   --- the size of the sample in y direction
                    (along the current)
           
           E    --- Energy of the particle
           
           Phi  --- Magnetic flux via one square
           
           wrap --- periodic b.c. along x, if true
           
      Return value: (tvals, check), where
      
           tvals --- the set of transmission coefficients

           check --- unitarity check
           
    """
    tvals, check = _tmatrix.gettrans(current.upper(), gauge.upper(), 
                                    Lx, Ly, E, phi, wrap)
    return tvals, check

def conductance (tvals):
    return sum (numpy.array(tvals)**2)

