#include <stdio.h>

#define L 10 //colums
#define M 10 //rows

int main(){

  int array[M][L];

  int curL=0;
  int curM=0;
  //fill array with Hexagon
  while (curM<M){
    while (curL<L){
      array[curM][curL]=(((curL%4)==1 || (curL%4)==2)*!(curM%2)) +(!((curL%4)==1 || (curL%4)==2)*(curM%2)) ;
      curL++;
    }
    curL=0;
    curM++;
  }


  //Debug print array
  curL=0;
  curM=0;
  while (curM<M){
    while (curL<L){
      printf("%i", array[curM][curL]);
      curL++;
    }
    printf("\n");
    curL=0;
    curM++;
  }


  //Generate Hamiltonian (when it is regular can use rules, with doping want to loop through array)
  //number sites from bottom left rightwards, then upwards - how do we only count valid sites?
  // since even with doping it is mostly regular, doing a full loop through is stupid

  return 0;
}
