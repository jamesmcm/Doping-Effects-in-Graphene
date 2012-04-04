#include <stdio.h>
#include <stdlib.h>

double gettrans_ (char   *current, 
		  char   *gauge, 
		  double *tvals, 
		  int    *nt, 
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
     char cur[] = "Y"; 
     char gauge[] = "X";
     double tvals[10000]; 
     double check; 
   
     int Lx = 2; 
     int Ly = 4;
     int wrap = 1; 
     int Nt = 0; 
   
     check = gettrans_(cur, gauge, tvals, &Nt, &Lx, &Ly, &E, &Phi, &wrap); 
     printf ("E = %g check = %g G = %g\n", E, check, G(tvals, Nt)); 
}
