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
double transmission (double lam, int n); 

int main (){

  int k = 0;
  long int i = 0;  // double i is weird!
  long double eps; // lam not needed, really
                   // Ditto a, b, c, d
  long double t,t_neg;
  long double top,down,t_d;

  int n   = 0;
  eps     = 0;
  top     = 0;
  down    = 0;
   
  scanf("%d",&n);
  
  for(i = 0 ; i <= 100000L;  i++ ){
	eps = (i/10000.0)-5.00; // Never, ever divide integer by integer,
                                // if the result is supposed to be an integer
				// C would do the division
        t     = transmission (eps + 1.0, n); 
        t_neg = transmission (eps - 1.0, n); 

        // The logic with all four branches can be simplified
        t_d = t+t_neg;
      
        printf("%.20Le %.20Le\t%.20Le\t%.20Le\n",eps,t_d, t, t_neg);
  }
		
  return 0;

}

double transmission (double lam, int n) {
       if (fabs (lam) > 2.0) 
           return real(lam, n); 
       return comp (lam, n); 
  
}

//Function for complex lambda
double comp(double lam,int n){
   double a_n = 0,b_n   = 0,c_n = 0,d_n = 0,
          t   = 0,t_neg = 0, arg1, arg;
   double test = 0;
   
   arg = atan2( sqrt(1.0 - pow(lam, 2) / 4.0) , lam/2.0 );
   
   // a, b, c, d below are formally singular at arg == 0, 
   // however, they have a well-defined limit. 
   // To prevent floating point exceptions, we introduce 
   // regularisation: 
   // 
   if (fabs (arg) < 1e-8)
       arg = 1e-8; 
   
   if (fabs (arg - M_PI) < 1e-8)
       arg = M_PI - 1e-8; 
   
   a_n =  - sin( (n-1) * arg )  /  sin(arg) ;
	
   b_n =    sin( n * arg )  /  sin(arg) ;

   c_n =  - sin(n*arg)      /  sin(arg)  ;

   d_n =    sin( (n+1)*arg ) /  sin(arg)  ;   
	
   t =  4.0 / ( pow( a_n + d_n, 2) + pow( c_n - b_n, 2) );

   return t;

}
//Function for real lambda
double real(double lam,int n){
   double a_n = 0, b_n = 0, c_n = 0, d_n = 0, 
          t = 0, xi = 0, cons;
    
   if (fabs (lam - 2) < 1e-10)
       lam = 2 + 1e-10;
   
   if (fabs (lam + 2) < 1e-10)
       lam = -2 - 1e-10; 
   
   xi = lam / 2.0 + pow( pow(lam, 2) / 4.0 - 1.0, 0.5);
   
   cons = 1.0 / sqrt( pow(lam, 2) - 4.0);

   a_n =   -cons * ( pow( xi, n-1 ) - pow( xi, -(n-1) ) );

   b_n =    cons * ( pow( xi, n   ) - pow( xi, -n     ) );

   c_n=  -  cons * ( pow( xi, n   ) - pow( xi, -n     ) );

   d_n=     cons * ( pow( xi, n+1 ) - pow( xi, -(n+1)  ) );

   t = 4.0 / ( pow( a_n + d_n, 2) + pow( c_n - b_n, 2) );
   
   return t;
}

double f_top(double t,double t_neg){
       return (t > t_neg) ? t : t_neg; // You can return the value immediately
}

double f_down(double t,double t_neg){
       return (t_neg < t) ? t_neg : t;
}     





