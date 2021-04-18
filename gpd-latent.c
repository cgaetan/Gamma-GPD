#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include "memory.h"
#include <math.h> 
#include <R.h>
#include <R_ext/Arith.h>	/* NA handling */
#include <Rmath.h>	
#include <R_ext/Random.h>	/* ..RNGstate */
#include <R_ext/Applic.h>	/* NA handling */
#define SMALL 1e-320
#define MINF -1.0e15
#define max(x, y)  ( (x < y) ? y : x)


double dist_1(double a1, double b1, double a2, double b2 )
/*
evaluate the euclidean distance
*/
{

	double val, d1 , d2;
	
	d1 = a1 - a2;
	d2 = b1 - b2;
	val = fabs(d1)+fabs(d2);


	return(val);
}
double dist_max(double a1, double b1, double a2, double b2 )
/*
evaluate the euclidean distance
*/
{

	double val, d1 , d2;
	
	d1 = a1 - a2;
	d2 = b1 - b2;
	val = max(d1,d2);


	return(val);
}


double dist_2(double a1, double b1, double a2, double b2 )
/*
evaluate the euclidean distance
*/
{

	double val, d1 , d2;
	
	d1 = a1 - a2;
	d2 = b1 - b2;
	val = R_pow_di(d1, 2) + R_pow_di(d2, 2);
	val = sqrt(val);   

	return(val);
}




double ellipse_check_old(double x, double y, double x1, double y1 , double radiusx, double radiusy)
{
/*
check if 

$\frac{(x-x_1)^2}{r_x^2} + \frac{(y-y_1)^2}{r_y^2} \leq 1$
*/	
	int answer; 
	double val, dx , dy;
	

	dx = (x - x1)/radiusx;
	dy = (y - y1)/radiusy;
	val = R_pow_di(dx, 2) + R_pow_di(dy, 2);
	answer=0;
	if ( val  <= 1.0) 
		answer = 1;


	return(answer);
}


double ellipse_check(double x0, double y0, double x, double y , double a, double  b, double phi)
{
/*	
	 x0,y0 ellipse center coordinates
    phi  counterclockwise angle the major axis makes with respect to the  x-axis. 
    a,b be the semi-major and semi-minor axes, 
 	x,y  arbitrary point coordinates

*/	
	 double answer=0.0, X, Y, cosphi, sinphi, val;
	 
/* Translate and rotate coords.
 * X = (x-x0)*cos(t)+(y-y0)*sin(t); 
   Y = -(x-x0)*sin(t)+(y-y0)*cos(t);  to align with ellipse	 
*/
     cosphi = cos(phi);
     sinphi = sin(phi);
     X = (x-x0)*cosphi+(y-y0)*sinphi;
     Y = -(x-x0)*sinphi+(y-y0)*cosphi;
/*If X^2/a^2+Y^2/b^2
is less than 1, the point  lies inside the ellipse. If it equals 1, it is right on
 the ellipse. If it is greater than 1, the point is outside.
*/     
 
     X = X/a;
     Y = Y/b;
     val = R_pow_di(X, 2) + R_pow_di(Y, 2);
	 if ( val  <= 1.0) 
		answer = 1.0;
     return(answer);		

}
void integration2d(double *x, double *y, double *z,  int *n, double *a, double *b, double *phi, double *X,  double *Y, double *U, int *m)
{
/*
 This function calculates $\Gamma(I_A)=\sum_{j=1}^m u_jI_{\tau_j\in A}(\tau_j)$ in 2D
 */	
	double  val;
    int i, j;

	for (i = 0; i < *n; i++) {		
		z[i]=0.0;
		for (j = 0; j < *m; j++) { 
			
			val=ellipse_check(x[i],y[i], X[j], Y[j], *a, *b, *phi);
			z[i]+=val*U[j];
			
			}			
		}	
}
 

void integration3d(double *x, double *y, double *z, double *f, int *n, double *a, double *b, double *phi, double *duration,
                       double *X,  double *Y, double *Z,  double *U, double *V, int *m)
{
		double dz, xtilde, ytilde, val;
		int i, j;
/*
 This function calculates $\Gamma(I_A)=\sum_{j=1}^m u_jI_{\tau_j\in A}(\tau_j)$ in 3D
 */	
	
	for (i = 0; i < *n; i++) {		
		for (j = 0; j < *m; j++) { 
			if (Z[j] <= z[i]) {
				dz = z[i]-Z[j];
				if ( dz < *duration) {
					xtilde = x[i]-dz*V[0];
					ytilde = y[i]-dz*V[1];
					val = ellipse_check(xtilde,ytilde, X[j], Y[j], *a, *b, *phi);
					f[i] += val*U[j];			
				}	
			}							
		}			
	}	
}

