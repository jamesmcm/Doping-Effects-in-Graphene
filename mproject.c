#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h> //Works on OS X ONLY
#include <complex.h>

//Complex double: zheev.f Computes all eigenvalues and, optionally, eigenvectors of a complex, Hermitian matrix.
//

#define L 3 //colums
#define M 3 //rows
#define T 1
#define WRAPX 0
#define WRAPY 0


int main(){
  int i;
  // create Hamiltonian
  double complex *Hcomplex = malloc(L*M*L*M*sizeof(double complex));
  //Initialise array to zeroes
  for(i=0;i<(L*M*L*M);i++){
    *(Hcomplex+i)=0;
  }
  //*Hcomplex = 3.0 + (4.0*I);
  //*(Hcomplex+1)= 2.0+(1.0*I);

  //printf("%f + i%f\n", creal(*Hcomplex), cimag(*Hcomplex));
  //printf("%f + i%f\n", creal(*(Hcomplex+1)), cimag(*(Hcomplex+1)));

  //remember memory must be stored columns first, matrix must be hermitian
  int x=0;
  int y=0;
  int p=0;
  
  /* while(p<9){ */
  /*   printf("%i:%i\n",p, p+(((p%2)-!(p%2))*((p/L)%2))+(((-1*(p%2))+!(p%2))*!((p/L)%2)) ); */
  /*   p++; */
  /* } */

  //Hamiltonian must be Hermitian, so just generate up to diagnol and then flip
  while (y<M){
    while(x<L){
      p=(x+(y*L)); //maybe use macro for this - p is in non-transposed space
      ptrans=(y+(x*M));
      if (((ptrans+L)/M) != (ptrans/M)){ //Node above in real space brick model is always connected, when transposed the Y conditions become the complicated ones!
	//*(Hcomplex+(L*(p+L))+p)=T;
	*(Hcomplex+ptrans+L; // doesn't work
	//can use old conditions and swap for ptrans and M?q
      }
      else {
	//*(Hcomplex+(L*((p+L)-(M*L)))+p)=WRAPY*T;
      }
      if ((p-L)>=0){
	//*(Hcomplex+(L*(p-L))+p)=T;
      }
      else {
	//*(Hcomplex+(L*((p-L)+(M*L))+p))=WRAPY*T;
	
      }
      if (((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2)) +  (((-1*((p%L)%2))  +!((p%L)%2))*!((p/L)%2)))/L) != (p/L)){//connect xnode but this changes depending on position 
	//*(Hcomplex+p+(L*((((((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))/L)-(p/L))<0)*((((p/L)+1)*L)-1))+(((((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2)))/L)-(p/L))>0)*((((p/L))*L))))))=WRAPX*T;
      }
      else {
	//*(Hcomplex+p+(L*((p+((((p%L)%2)-!((p%L)%2))*((p/L)%2))+(((-1*((p%L)%2))+!((p%L)%2))*!((p/L)%2))))))=T;
      }
      x++;
    }
    y++;
    x=y;
  }

  //FLIP MATRIX

  //Print Hamiltonian
  i=0;
  int j=0;
  while (j<(L*M)){
    while (i<(L*M)){
      printf("%.0f ", creal(*(Hcomplex+j+(L*i))));
      i++;
    }
    i=0;
    j++;
    printf("\n");
  }

  return 0;
}
