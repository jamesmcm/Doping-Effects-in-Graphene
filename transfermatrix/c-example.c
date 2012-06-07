#include <stdio.h>
#include <stdlib.h>

double gettrans_ (char   *current, 
		  char   *gauge, 
		  double *tvals, 
		  int    *ntvals, 
		  int    *limx,     // Fortran: pass by ref
		  int    *limy,     // ditto 
		  double *V, 
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
     int maxl = 10000, ntvals = 0; 
     int iv = 0, jv = 0; 
   
     double *tvals = NULL; 
     double check = 0.0;   
     double *V = NULL; 
     char current = 'Y'; 
     char gauge = 'Y'; 
   
     tvals = (double *)malloc( maxl * sizeof(*tvals) );
     V     = (double *)malloc( Lx * Ly * sizeof (*V)); 
     for (iv = 0; iv < Lx; iv ++) {
         for (jv = 0; jv < Ly; jv++) {
	     V[iv * Lx + jv] = 0.0; 
	 }
     }
     check = gettrans_(&current, &gauge, 
		       tvals, &ntvals, &Lx, &Ly, V, &E, &Phi, &wrap); 
     printf ("E = %g check = %g G = %g\n", E, check, G(tvals, ntvals)); 
     free(V); 
     free(tvals); 
     
}

