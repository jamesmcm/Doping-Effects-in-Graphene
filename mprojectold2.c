#include <stdio.h>
#include <Accelerate/Accelerate.h> //Works on OS X ONLY

//Complex double: zheev.f Computes all eigenvalues and, optionally, eigenvectors of a complex, Hermitian matrix.
//

#define L 4 //colums
#define M 3 //rows
#define WRAPX 1
#define WRAPY 1

int main(){

  // create Hamiltonian
  int H[L*M][L*M]={0}; //might want to use static here to ensure 0 setting in memory

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

  //Print Hamiltonian
  int i=0;
  int j=0;
  while (j<(L*M)){
    while (i<(L*M)){
      printf("%i ", H[j][i]);
      i++;
    }
    i=0;
    j++;
    printf("\n");
  }

  return 0;
}
