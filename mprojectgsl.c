<<<<<<< HEAD
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

//Note ints limit 2147483647

#define L 3 //colums
#define M 4 //rows
#define WRAPX 1
#define WRAPY 0
#define RESTRICT(x, L) (((x) >= 0) ? ((x) % (L)) : ((x) + (L)))  // Restrict x to the interval [0:L-1]
                                                                 // This is needed for wrapping
                                                                 // Note: does not work when x < -L
#define SITE(x, y, L, M) (RESTRICT((x), (L)) + RESTRICT ((y), (M)) * (L) )  // Site index 

int main(){

  // create Hamiltonian
  static double H[L*M][L*M]={0}; //might want to use static here to ensure 0 setting in memory 
                                 // - this segfaults for large numbers of L,M
  //memset(H, 0, sizeof(double));
  int x=0;
  int y=0;


  //array is always Hermitian or symmetric so should only need 
  //to run on half the values - fix this!
  for (y = 0; y < M; y ++) {
    for (x = 0; x < L; x++) {
      int current     = SITE (x, y, L, M); 
      int col_stagger = (x % 2) ? 1 : -1; 
      int row_stagger = (y % 2) ? 1 : -1;
      int ab_stagger  = col_stagger * row_stagger; // +1 on A, -1 on B sublattice 
      int h_neigh     = SITE (x + ab_stagger, y, L, M); // horizontal neighbour,
                                                        // left or right, depending on AB
      int up_neigh    = SITE (x, y + 1, L, M);          // A neighbour upstairs
      int down_neigh  = SITE (x, y - 1, L, M);          // A neighbour downstairs
      
      if ((y < M - 1) || WRAPY ) { //link above
          H[current][up_neigh] = 1;
      }

      if ((y > 0) || WRAPY) {      //link below
          H[current][down_neigh] = 1;
      }

      if (((x > 0) && (x <  L - 1)) || WRAPX) 
           H[current][h_neigh] = 1; 
      }
     
  }
  int k;
  for(k=0; k<L*M; k++){
  int num = rand()%100;
  if(num>=80){
    H[k][k]=99;
  }
}

  //Print Hamiltonian
  int i=0;
  int j=0;
  while (j<(L*M)){
    while (i<(L*M)){
      printf("%.0f ", H[j][i]);
      i++;
    }
    i=0;
    j++;
    printf("\n");
  }
  //double data[L*M*L*M];
  /* gsl_matrix_view m = gsl_matrix_view_array (*H, L*M, L*M); */
     
  /* gsl_vector *eval = gsl_vector_alloc (L*M); */
  /* //gsl_matrix *evec = gsl_matrix_alloc (L*M, L*M); */
  /* gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (L*M); //Just worry about eigenvalues for now */
  /* //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (L*M); */
  /* gsl_eigen_symm (&m.matrix, eval, w); */
  /* //gsl_eigen_symmv (&m.matrix, eval, evec, w); */
  /* //eigenvalues are symmetric so really only need half, don't think this can be fixed though */
     
  /* gsl_eigen_symm_free (w); */
     
  /* //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC); */
  /* int i=0; */
  /* for (i = 0; i < L*M; i++) */
  /* { */
  /* double eval_i = gsl_vector_get (eval, i); */
  /* //gsl_vector_view evec_i = gsl_matrix_column (evec, i); */
  /* printf ("%g\n", eval_i); //eigenvalue */
  /* //printf ("eigenvector = \n"); */
  /* //gsl_vector_fprintf (stdout, &evec_i.vector, "%g"); */
  /* } */
  /* gsl_vector_free (eval); */
  /* //gsl_matrix_free (evec); */
  return 0;
}
=======
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

//Note ints limit 2147483647

#define L 3 //colums
#define M 4 //rows
#define WRAPX 1
#define WRAPY 0

int main(){

  // create Hamiltonian
  static double H[L*M][L*M]={0}; //might want to use static here to ensure 0 setting in memory - this segfaults for large numbers of L,M
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
      x++;
    }
    y++;
    x=0;
  }
  int k;
  for(k=0; k<L*M; k++){
  int num = rand()%100;
  if(num>=80){
    H[k][k]=99;
  }
}

  //Print Hamiltonian
  int i=0;
  int j=0;
  while (j<(L*M)){
    while (i<(L*M)){
      printf("%.0f ", H[j][i]);
      i++;
    }
    i=0;
    j++;
    printf("\n");
  }
  //double data[L*M*L*M];
  /* gsl_matrix_view m = gsl_matrix_view_array (*H, L*M, L*M); */
     
  /* gsl_vector *eval = gsl_vector_alloc (L*M); */
  /* //gsl_matrix *evec = gsl_matrix_alloc (L*M, L*M); */
  /* gsl_eigen_symm_workspace * w = gsl_eigen_symm_alloc (L*M); //Just worry about eigenvalues for now */
  /* //gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (L*M); */
  /* gsl_eigen_symm (&m.matrix, eval, w); */
  /* //gsl_eigen_symmv (&m.matrix, eval, evec, w); */
  /* //eigenvalues are symmetric so really only need half, don't think this can be fixed though */
     
  /* gsl_eigen_symm_free (w); */
     
  /* //gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC); */
  /* long int i=0; */
  /* for (i = 0; i < L*M; i++) */
  /*   { */
  /*     double eval_i = gsl_vector_get (eval, i); */
  /*     //gsl_vector_view evec_i = gsl_matrix_column (evec, i); */
  /*     printf ("%g\n", eval_i); //eigenvalue */
  /*     //printf ("eigenvector = \n"); */
  /*     //gsl_vector_fprintf (stdout, &evec_i.vector, "%g"); */
  /*          } */
  /* gsl_vector_free (eval); */
  /* //gsl_matrix_free (evec); */
  return 0;
}
>>>>>>> af45f6c95a2406730d51c5b5b2a7126e6e639083
