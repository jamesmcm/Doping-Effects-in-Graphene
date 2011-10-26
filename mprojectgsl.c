#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
//#include <Accelerate/Accelerate.h> //Works on OS X ONLY

//Complex double: zheev.f Computes all eigenvalues and, optionally, eigenvectors of a complex, Hermitian matrix.
//

#define L 50 //colums
#define M 50 //rows
#define WRAPX 1
#define WRAPY 0

int main(){

  // create Hamiltonian
  static double H[L*M][L*M]={0}; //might want to use static here to ensure 0 setting in memory - this segfaults for large numbers of L,M
  //memset(H, 0, sizeof(double));
  int x=0;
  int y=0;
  int p=0;
  
  /* while(p<9){ */
  /*   printf("%i:%i\n",p, p+(((p%2)-!(p%2))*((p/L)%2))+(((-1*(p%2))+!(p%2))*!((p/L)%2)) ); */
  /*   p++; */
  /* } */

  while (y<M){
    while(x<L){
      p=(x+(y*L)); //maybe use macro for this
      if ((p+L)<(M*L)){
  	H[p][(p+L)]=1;
      }
      else {
  	H[p][(p+L)-(M*L)]=WRAPY;
      }
      if ((p-L)>=0){
  	H[p][(p-L)]=1;
      }
      else {
  	H[p][(p-L)+(M*L)]=WRAPY;
      }
      //printf("p:%i -> %i\n", p, (p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2))));
      if (((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))/L) != (p/L)){
  	H[p][(((((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))/L)-(p/L))<0)*((((p/L)+1)*L)-1))+(((((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))/L)-(p/L))>0)*((((p/L))*L)))]=WRAPX;
      }
      else {
  	H[p][(p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))]=1;
      }
      x++;
    }
    y++;
    x=0;
  }

  /* //Print Hamiltonian */
  /* int i=0; */
  /* int j=0; */
  /* while (j<(L*M)){ */
  /*   while (i<(L*M)){ */
  /*     printf("%.0f ", H[j][i]); */
  /*     i++; */
  /*   } */
  /*   i=0; */
  /*   j++; */
  /*   printf("\n"); */
  /* } */
  //double data[L*M*L*M];
  gsl_matrix_view m = gsl_matrix_view_array (*H, L*M, L*M);
     
  gsl_vector *eval = gsl_vector_alloc (L*M);
  //gsl_matrix *evec = gsl_matrix_alloc (L*M, L*M);
  gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (L*M); //Just worry about eigenvalues for now
  //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (L*M);
  gsl_eigen_symm (&m.matrix, eval, w);       
  //gsl_eigen_symmv (&m.matrix, eval, evec, w);
     
  gsl_eigen_symm_free (w);
     
  //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  int i=0;
  for (i = 0; i < L*M; i++)
    {
      double eval_i = gsl_vector_get (eval, i);
      //gsl_vector_view evec_i = gsl_matrix_column (evec, i);
      printf ("%g\n", eval_i); //eigenvalue
      //printf ("eigenvector = \n");
      //gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
           }
  gsl_vector_free (eval);
  //gsl_matrix_free (evec);
  return 0;
}
