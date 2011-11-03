#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h> //Initalise rand()

//Note ints limit 2147483647

#define L 60 //colums
#define M 60 //rows
#define WRAPX 1
#define WRAPY 0

int main(){

 int num=0;	
 int k=0; //Random function int
srand(time(NULL));

  // create Hamiltonian
  static double H[L*M][L*M]; //might want to use static here to ensure 0 setting in memory - this segfaults for large numbers of L,M
  //memset(H, 0, sizeof(double));
  long int x=0;
  long int y=0;
  long int p=0;
  
  /* while(p<9){ */
  /*   printf("%i:%i\n",p, p+(((p%2)-!(p%2))*((p/L)%2))+(((-1*(p%2))+!(p%2))*!((p/L)%2)) ); */
  /*   p++; */
  /* } */

  //array is always Hermitian or symmetric so should only need to run on half the values - fix this!
  while (y<M){
    while(x<L){
      p=(x+(y*L)); //maybe use macro for this
      int col=p%L;
      int row=p/L;
      int cond=(p+((((col)%2)-!((col)%2))*((row)%2))+(((-1*((col)%2))+!((col)%2))*!((row)%2))); //gives position of X link


      if ((p+L)<(M*L)){//link above
  	H[p][(p+L)]=1;
      }
      else {
  	H[p][(p+L)-(M*L)]=WRAPY;
      }

      if ((p-L)>=0){//link below
  	H[p][(p-L)]=1;
      }
      else {
  	H[p][(p-L)+(M*L)]=WRAPY;
      }

      if ((cond/L) != (row)){//xlink
  	H[p][((((cond/L)-(row))<0)*((((row)+1)*L)-1))+((((cond/L)-(row))>0)*((((row))*L)))]=WRAPX;
      }
      else {
  	H[p][cond]=1;
      }
      if(x==y){
		while(k<L*M){
  		   num = rand()%100;
    			if(num>=80){
        			H[k][k]=70;
    			}
     			++k;
 		}
	 }
      x++;
    }
    y++;
    x=0;
  }

  /* //Print Hamiltonian */
  /*int t=0; 
   int r=0; 
   while (r<(L*M)){ 
     while (t<(L*M)){ 
       printf("%.0f ", H[r][t]); 
       t++; 
     } 
     t=0; 
     r++; 
     printf("\n"); 
   } */

  //double data[L*M*L*M];
  gsl_matrix_view m = gsl_matrix_view_array (*H, L*M, L*M);
     
  gsl_vector *eval = gsl_vector_alloc (L*M);
  //gsl_matrix *evec = gsl_matrix_alloc (L*M, L*M);
  gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (L*M); //Just worry about eigenvalues for now
  //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (L*M);
  gsl_eigen_symm (&m.matrix, eval, w);
  //gsl_eigen_symmv (&m.matrix, eval, evec, w);
  //eigenvalues are symmetric so really only need half, don't think this can be fixed though
     
  gsl_eigen_symm_free (w);
     
  //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);
  long int i=0;
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
