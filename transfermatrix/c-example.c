#include <stdio.h>
#include <stdlib.h>

double gettrans_ (double *tvals, 
		  int    *limx,     // Fortran: pass by ref
		  int    *limy,     // ditto 
		  double *E,        
		  double *Phi, 
		  int    *wrap
		);

double G (double *tvals, int Lx) {
    double sum = 0.0; 
    int i = 0; 
   
    for (i = 0; i < Lx; i ++)
      sum += tvals[i] * tvals[i]; 
}
void main () {
     double E = -1.0; 
     double Phi = 0.0;
   
     int Lx = 2; 
     int Ly = 4;
     int wrap = 1; 
   
     double *tvals = NULL; 
     double check = 0.0;   
   
     tvals = (double *)malloc( Lx * sizeof(*tvals) );
     check = gettrans_(tvals, &Lx, &Ly, &E, &Phi, &wrap); 
     printf ("E = %g check = %g G = %g\n", E, check, G(tvals, Lx)); 
     free(tvals); 
}
