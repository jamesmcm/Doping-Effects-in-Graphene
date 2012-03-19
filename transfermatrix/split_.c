#include <stdlib.h>
#include <string.h>

void split_ (const int *plimx, const double *mult, 
	     double *a, double *b, double *c, double *d) {
      int i; 
   
      int limx = *plimx; 
      int da   = 2 * limx; 
      int dm   = 2 * da; 
      int szcp = da * sizeof (double);
   
      double *pmult = mult; 
      double *pa = a; 
      double *pb = b;
      double *pc = c; 
      double *pd = d;
     
      for (i = 0; i < limx; i ++, pa += da, pb += da, pmult += dm) {
	   //double *pmult = mult + 4 * limx * i; 
	   //double *pa = a + 2 * limx * i; 
	   //double *pb = b + 2 * limx * i; 
           memcpy (pa, pmult,      szcp); 
           memcpy (pb, pmult + da, szcp); 
      }
   
      pmult = mult + 2 * limx * da; 
      for (i = 0; i < limx; i ++, pc += da, pd += da, pmult += dm) {
//	   double *pmult = mult + 4 * limx * i + 4 * limx * limx; 
//	   double *pc = c + 2 * limx * i; 
//	   double *pd = d + 2 * limx * i; 
           memcpy (pc, pmult,      szcp); 
           memcpy (pd, pmult + da, szcp); 
      }
}
