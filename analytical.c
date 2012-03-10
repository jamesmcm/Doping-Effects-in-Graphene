#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define SIZE 2   //Range over which epsilon varies

//#define N 7    //Program is specific to a given N value

double comp(double lam,int n);
double real(double lam,int n);
double f_top(double t,double t_neg);
double f_down(double t,double t_neg);

int main (){

int k=0;
double i=0;
long double eps,lam,lam_neg;
long double a_n,b_n,c_n,d_n;
long double a_neg,b_neg,c_neg,d_neg;
long double t,t_neg;
long double top,down,t_d;

int n=0;
eps=0;
lam=0;
lam_neg=0;
top=0;
down=0;
      scanf("%d",&n);
  for(i=0;i<100000;i++){
	eps=((i/10000)-5.00);
	lam = (eps+1);
	lam_neg= (eps-1);	



		if((fabs(lam)<2)&&(fabs(lam_neg)<2)){  
    		t= comp(lam,n);
		t_neg=comp(lam_neg,n);

		//printf("%e %e %e\n",eps,t,t_neg);

		}
		if((fabs(lam)>2)&&(fabs(lam_neg)<2)){

		 t=real(lam,n);
		 t_neg=comp(lam_neg,n);
                
		//printf("%e %e %e\n",eps,t,t_neg);
		}

		if((fabs(lam)<2)&&(fabs(lam_neg)>2)){
		 t=comp(lam,n);  
    		 t_neg= real(lam_neg,n);
           
		 //printf("%e %e %e\n",eps,t,t_neg);

		}
		if((fabs(lam_neg)>2)&&(fabs(lam)>2)){
		 t=real(lam,n);
		 t_neg=real(lam_neg,n);
            
		//printf("%e %e %e\n",eps,t,t_neg);
		}
      
/*     
      down  = f_down(t,t_neg);
      top     = f_top(t,t_neg);
*/
      t_d = t+t_neg;
      
     

      printf("%.6e %.20e %.20e\n",eps,t_d);
   }
		


return 0;

}
//Function for complex lambda
double comp(double lam,int n){

double a_n=0,b_n=0,c_n=0,d_n=0,t=0,t_neg=0,arg1,arg;
double test=0;

arg = (atan2((sqrt(1-pow(lam,2)/4)),(lam/2)));

			a_n = (( -1*sin((n-1)*arg)) / (sin(arg)));
	
			b_n =  ((sin(n*(arg))) / (sin(arg)));

			c_n = ((-1*sin(n*arg)) / (sin(arg)));

			d_n = (sin((n+1)*arg) / (sin(arg)));   
	
			t = ((4)/(pow((a_n+d_n),2)+pow((c_n-b_n),2)));

return (t);

}
//Function for real lambda
double real(double lam,int n){
double a_n=0,b_n=0,c_n=0,d_n=0,t=0,xi=0,cons;

xi=((lam/2)+pow(( ((pow(lam,2))/(4))-1),0.5));
cons= (1/(sqrt( pow(lam,2)-4)));

			a_n= (-1*cons*((pow(xi,(n-1)))-(pow(xi,(-1*(n-1))))));

			b_n= (cons*((pow(xi,n))-(pow(xi,(-1*n)))));

			c_n= (-1*cons*((pow(xi,n))-(pow(xi,(-1*n)))));

			d_n= (cons*((pow(xi,(n+1)))-(pow(xi,(-1*(n+1))))));

			t = ((4)/(pow((a_n+d_n),2)+pow((c_n-b_n),2)));

return (t);
}

double f_top(double t,double t_neg){
 t   = (t>t_neg)? t:t_neg;


return (t);

}

double f_down(double t,double t_neg){

 t_neg =(t_neg<t)? t_neg:t;


return(t_neg);
}     





