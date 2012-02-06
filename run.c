#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define SIZE 2   //Range over which epsilon varies

#define N 7    //Program is specific to a given N value

double comp(double lam);
double real(double lam);

int main (){

int k=0;
double i=0;
double eps,lam,lam_neg;
double a_n,b_n,c_n,d_n;
double a_neg,b_neg,c_neg,d_neg;
double t,t_neg;

eps=0;
lam=0;
lam_neg=0;

  for(i=0;i<1000;i++){
	eps=((i/100)-5.001);
	lam = (eps+1);
	lam_neg= (eps-1);	



		if((fabs(lam)<2)&&(fabs(lam_neg)<2)){  
    		t= comp(lam);
		t_neg=comp(lam_neg);
		printf("%f %f %f\n",eps,t,t_neg);

		}
		if((fabs(lam)>2)&&(fabs(lam_neg)<2)){

		 t=real(lam);
		 t_neg=comp(lam_neg);
		printf("%f %f %f\n",eps,t,t_neg);
		}

		if((fabs(lam)<2)&&(fabs(lam_neg)>2)){
		 t=comp(lam);  
    		 t_neg= real(lam_neg);
		 printf("%f %f %f\n",eps,t,t_neg);

		}
		if((fabs(lam_neg)>2)&&(fabs(lam)>2)){
		 t=real(lam);
		 t_neg=real(lam_neg);
		printf("%f %f %f\n",eps,t,t_neg);
		}

   }
		


return 0;

}

double comp(double lam){

double a_n=0,b_n=0,c_n=0,d_n=0,t=0,t_neg=0,arg1,arg;
double test=0;

arg = (atan2((sqrt(1-pow(lam,2)/4)),(lam/2)));

			a_n = (( -1*sin((N-1)*arg)) / (sin(arg)));
	
			b_n =  ((sin(N*(arg))) / (sin(arg)));

			c_n = ((-1*sin(N*arg)) / (sin(arg)));

			d_n = (sin((N+1)*arg) / (sin(arg)));   
	
			t = ((4)/(pow((a_n+d_n),2)+pow((c_n-b_n),2)));

return (t);

}

double real(double lam){
double a_n=0,b_n=0,c_n=0,d_n=0,t=0,xi=0,cons;

xi=((lam/2)+pow(( ((pow(lam,2))/(4))-1),0.5));
cons= (1/(sqrt( pow(lam,2)-4)));

			a_n= (-1*cons*((pow(xi,(N-1)))-(pow(xi,(-1*(N-1))))));

			b_n= (cons*((pow(xi,N))-(pow(xi,(-1*N)))));

			c_n= (-1*cons*((pow(xi,N))-(pow(xi,(-1*N)))));

			d_n= (cons*((pow(xi,(N+1)))-(pow(xi,(-1*(N+1))))));

			t = ((4)/(pow((a_n+d_n),2)+pow((c_n-b_n),2)));

return (t);
}











