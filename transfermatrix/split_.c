#include <stdlib.h>
#include <string.h>

void split_ (const int *plimx, const double *mult, 
	     double *a, double *b, double *c, double *d) {
      // HACK! The arrays are in fact complex-valued, 
      // but here I interpret them as double arrays, 
      // by doubling their length
      int i; 
  
      int limx = *plimx;  
      int da   = 2 * limx; // 2 because of sizeof(complex) 
      int dm   = 2 * da;   // two halves of a col belong to
                           // two different matrices 
      int szcp = da * sizeof (double); // Chunk to copy
   
      double *pmult = mult;  // Set initial pointers
      double *pa = a; 
      double *pb = b;
      double *pc = c; 
      double *pd = d;
     
      // First limx cols of MULT are to be split into 
      // cols of A and B
      // First half goes into A, second half into B, 
      // and we increment the pointers accordingly.
      // We rely on Fortran array ordering convention:
      // A(1, 1), A(2, 1), A(3, 1), ... A(N, 1), A(1, 2), A(2, 2), etc
      for (i = 0; i < limx; i ++, pa += da, pb += da, pmult += dm) {
           memcpy (pa, pmult,      szcp); 
           memcpy (pb, pmult + da, szcp); 
      }
   
      // The remaining cols of MULT are to be split into 
      // cols of C and D
      pmult = mult + 2 * limx * da; 
      for (i = 0; i < limx; i ++, pc += da, pd += da, pmult += dm) {
           memcpy (pc, pmult,      szcp); 
           memcpy (pd, pmult + da, szcp); 
      }
}
