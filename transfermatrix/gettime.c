#include <time.h>
#include <stdio.h>

/* int main(){ */

/*   time_t time1 = ((long)time(NULL)); */

/*   //printf("Current time is %li", time1); */

/*   return 0; */
/* } */

int gettime1_(){
  time_t time1 = ((int)time(NULL));

  return time1;
}
