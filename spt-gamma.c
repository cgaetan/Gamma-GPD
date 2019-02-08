#include <stdio.h>
#include <stdlib.h> /* for malloc & free */
#include "memory.h"
#include <math.h> 
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <R_ext/Arith.h>	/* NA handling */
#include <R_ext/Random.h>	/* ..RNGstate */
#include <R_ext/Applic.h>	/* NA handling */

#define INF 1.0e32
#define MINF -1.0e32
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


double dist(double x, double y, double x1, double y1)
{
	double val;
	val = (x-x1)*(x-x1)+(y-y1)*(y-y1);
	val = sqrt(val);
	return(val);
} 


double ellipse_ellipse_overlap_gsl (double PHI_1, double A1, double B1, 
                                double H1, double K1, double PHI_2, 
                                double A2, double B2, double H2, double K2,
                                double X[4], double Y[4], int * NROOTS,
                                int *rtnCode, int choice);

void ellipse_ellipse_overlap_gsl_wrap (double *PHI, double *A, double *B, 
                                double *H1, double *K1, double *H2, double *K2, double *val) {
					double  x[4], y[4];									
					int nroots;
					int  rtnCode = 0; 
					int  choice = 1;
							
					*val =ellipse_ellipse_overlap_gsl(*PHI, *A, *B, *H1, *K1 , *PHI, *A, *B, *H2 , *K2, x, y,  &nroots,	
											&rtnCode, choice);									
}
                                
                                

double slated_elliptic_cylinders_intersection (double V[2], double D, double PHI, double A, double B, 
                                double H1, double K1, double T1, double H2, double K2, double T2)
                                {
									double dt, x[4], y[4], val, HTILDE, KTILDE, tmp;
									int nroots;
									int  rtnCode = 0; 
									int  choice = 1;
									dt = fabs(T1-T2);
									val = (D- dt); 
									if (val <= 0.0) 
									{
										val = 0.0;
										return(val);
									}
									if ( T2 > T1) {
										HTILDE = H2-dt*V[0];
										KTILDE = K2-dt*V[1];
										tmp = ellipse_ellipse_overlap_gsl(PHI, A, B, H1, K1 , PHI, A, B, HTILDE , KTILDE, x, y,  &nroots,	
											&rtnCode, choice);
										if (tmp > (M_PI*A*B)) {
											tmp = M_PI*A*B;
										}	
										
										val = val*tmp;
										
									} 
									else 
									{
										HTILDE = H1-dt*V[0];
										KTILDE = K1-dt*V[1];
										tmp = ellipse_ellipse_overlap_gsl(PHI, A, B, HTILDE, KTILDE , PHI, A, B, H2 , K2, x, y,  &nroots,	
											&rtnCode, choice);
										if (tmp > (M_PI*A*B)) {
											tmp = M_PI*A*B;
										}
										val = val*tmp;
									}
									
									return(val); 
                                    
								}
                                

double tr(double y, double xi, double sigma, double kappa)
{
	double val;
	val = (1.0+xi*y/sigma);
	val = R_pow(val,1.0/xi)-1.0;
	val = (kappa+1.0)*val;
						
	return(val);
} 

double dtr(double y, double xi, double sigma, double kappa)
{
	double val;
	val = (1.0+xi*y/sigma);
	val = R_pow(val,1.0/xi-1.0);
	val = (kappa+1.0)*val/sigma;
	return(val);
} 

double  test(double y, double xi, double sigma)
{
	double val;
	val = (1.0+xi*y/sigma);
						
	return(val);
}

double LT(double *a, double b, double x)
{

	double val;
	 val = R_pow(b,a[0])*R_pow((x+b),-a[0]);
	
	return(val);
}


double dLT(double *a, double b, double x)
{

	double val;
	 val = -a[0]* R_pow(b,a[0])*R_pow((x+b),-(a[0]+1));
	return(val);
}

double LT2(double *a, double b, double x, double y)
{
	double val;
       val = R_pow(b,(a[1]+a[2]+a[3]))*R_pow((x+b),-a[1])*R_pow((x+y+b),-a[2])*R_pow((y+b),-a[3]);
    return(val);
}

double DxLT2(double *a, double b, double x, double y){
	double val;
	 
	val = - R_pow(b,(a[1]+a[2]+a[3]))* (a[1]*R_pow(x+b,-(a[1]+1.0))*R_pow(x+y+b,-a[2])*R_pow(y+b,-a[3]) 
	            +a[2]*R_pow(x+b,-a[1])*R_pow(x+y+b,-(a[2]+1.0) )*R_pow(y+b,-a[3]));
	return(val);
}

double DyLT2(double *a, double b, double x, double y){
	double val;
	
	val = - R_pow(b,(a[1]+a[2]+a[3]))* (a[3]*R_pow(x+b,-a[1])*R_pow(x+y+b,-a[2])*R_pow(y+b,-(a[3]+1.0) ) 
	            +a[2]*R_pow(x+b,-a[1])*R_pow(x+y+b,-(a[2]+1.0))*R_pow(y+b,-a[3]));
	return(val);
}



double DxyLT2(double *a, double b, double x, double y) {
	
	double  val;
	
	val = R_pow(b,(a[1]+a[2]+a[3]))* (a[1]*a[2]*R_pow(x+b,-(a[1]+1.0))*R_pow(x+y+b,-(a[2]+1.0))*R_pow(y+b,-a[3]) 
	            + a[1]*a[3]*R_pow(x+b,-(a[1]+1.0))*R_pow(x+y+b,-a[2])*R_pow(y+b,-(a[3]+1.0))
	            + a[2]*(a[2]+1.0)*R_pow(x+b,-a[1])*R_pow(x+y+b,-(a[2]+2.0))*R_pow(y+b,-a[3])
	             +a[2]*a[3]*R_pow(x+b,-a[1])*R_pow(x+y+b,-(a[2]+1.0))*R_pow(y+b,-(a[3]+1.0)));
	return(val);
}


void chi(double *v, double *xcoords, double *ycoords, double *tcoords, double *param, double *val)
{
	double  alpha, beta, a[4],  duration,  kappa, phi, semiax1, semiax2, velocity[2], volume,  X[4], Y[4];
    
    semiax1 = param[2];
	semiax2 = param[3];
	phi = param[4];
	
	duration = param[5];
	velocity[0] = param[6]; // velocity
	velocity[1] = param[7];
	
	alpha = param[8];
	beta = param[9];
	kappa = param[10];
	volume = M_PI*semiax1*semiax2*duration; // volume of slated cylinder
	a[0] = alpha ;
	a[2] = alpha*slated_elliptic_cylinders_intersection(velocity, duration, phi, semiax1, semiax2, 
                                xcoords[0],ycoords[0],tcoords[0], xcoords[1],ycoords[1],tcoords[1])/volume;
     *val= R_pow(1+2*(*v+kappa)/beta,-a[2])*R_pow(1+(*v+kappa)/beta,2*a[2]-a[0]);
}
	


void  pwl_spatiotemporal(double *y, int *z,  int *n, double *theta, double *loglik, double *u, double *xcoords, double *ycoords, double *tcoords, double *delta)
{
	double  alpha, beta, a[4], distspatial, disttemporal , duration, indloglik, kappa, jacx, jacy, phi, semiax1, semiax2, sigma, tmp,  velocity[2], volume,xi, xx, yy, X[4], Y[4];
	int count = 0, i, j, ord, nroots, rtnCode, zcase;
/* this function calculates the pairwise likelihood for the exceedandes */	
	xi = theta[0];
    sigma = theta[1];

	semiax1 = theta[2];
	semiax2 = theta[3];
	
	phi = theta[4];
	duration = theta[5];
	velocity[0] = theta[6]; // velocity
	velocity[1] = theta[7];
	
	alpha = theta[8];
	beta = theta[9];
	kappa = theta[10];
	volume = M_PI*semiax1*semiax2*duration; // volume of slated cylinder
    
	
	if ( sigma <= 0) {
		*loglik = MINF*(1-sigma)*(1-sigma);		
		return;
	}
	if ( kappa <= 0) {
		*loglik = MINF*(1-kappa)*(1-kappa);		
		return;
	}
     
    if ( semiax1 <= 0) {
		*loglik = MINF*(1-semiax1)*(1-semiax1);		
		return;
	}
	 
   if ( semiax2 <= 0) {
		*loglik = MINF*(1-semiax2)*(1-semiax2);		
		return;
	} 
  
   if ( duration <= 0) {
		*loglik = MINF*(1-duration)*(1-duration);		
		return;
	}

	if ( (phi < 0) || phi >= M_PI_2)  {
		*loglik = MINF*(1-theta[4])*(1-theta[4]);		
		return;
	}
   	
	
	
    *loglik = 0;
			
	
    a[0] = alpha;
	
	for (i = 0; i < (*n-1); i++) {
		if (!ISNA(y[i])) {
				
		    tmp = test(y[i],xi,sigma); /* test if the argument is positive */
						if ( tmp <= 0) {
							*loglik = MINF*(1-tmp)*(1-tmp);		
						return;
						}						
			for (j = (i+1); j < (*n); j++){   
				if (!ISNA(y[j])) {
					distspatial = dist(xcoords[i],ycoords[i],xcoords[j],ycoords[j]);
					disttemporal = fabs(tcoords[i]-tcoords[j]);
					if ( (distspatial <= delta[0]) && (disttemporal <= delta[1])  ) {
						count+=1;
						a[2] = alpha*slated_elliptic_cylinders_intersection(velocity, duration, phi, semiax1, semiax2, 
                                xcoords[i],ycoords[i],tcoords[i], xcoords[j],ycoords[j],tcoords[j])/volume;
						a[1] = a[3] = (a[0]-a[2]); 
						zcase = z[i]  + 2 * z[j]; 
						switch(zcase) {
							case 0: indloglik = 1-2*LT(a,beta,kappa) + LT2(a,beta,kappa,kappa); 
									*loglik += log(indloglik);
									break;
							case 1: tmp = test(y[j],xi,sigma); /* # above -> below */
									if ( tmp <= 0) {
										*loglik = MINF*(1-tmp)*(1-tmp);		
									return;
									}					
									xx = tr(y[i],xi,sigma,kappa)+kappa;
									yy = kappa;
									indloglik = - dLT(a,beta,xx) + DxLT2(a,beta,xx,yy);
									jacx= dtr(y[i],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacx);
									break;
							case 2: tmp = test(y[j],xi,sigma); /* # below -> above  */								
									if ( tmp <= 0) {
										*loglik = MINF*(1-tmp)*(1-tmp);		
									return;
									}					
									xx = kappa;
									yy = tr(y[j],xi,sigma,kappa)+kappa;	
									indloglik = - dLT(a,beta,yy) + DyLT2(a,beta,xx,yy);
									jacy= dtr(y[j],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacy);
									break;
							case 3: tmp = test(y[j],xi,sigma);      /*# above -> above*/
									if ( tmp <= 0) {
										*loglik = MINF*(1<-tmp)*(1-tmp);		
									return;
									}					
									xx	= tr(y[i],xi,sigma,kappa)+kappa; /* tr the data */
									yy = tr(y[j],xi,sigma,kappa)+kappa;
									indloglik = DxyLT2(a,beta,xx,yy);	
									jacx = dtr(y[i],  xi,  sigma, kappa);
									jacy = dtr(y[j],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacx)+log(jacy);
									break;
							}
					}
				}		
			}
		}	
	}
	}


/* version for parallel computing */
void  pwl_spatiotemporal_parallel(double *y, int *z, int *nobs, int *nsites, int *nblocks, int *lag, double *theta, double *loglik, double *u, double *xcoords, 
                            double *ycoords, double *tcoords, double *delta, double *count)
{
	double  alpha, beta, a[4], distspatial,  duration, indloglik, kappa, jacx, jacy, phi, semiax1, semiax2, sigma, tmp,  velocity[2], volume,xi, xx, yy, X[4], Y[4];
	int end, i, j, k,  nroots, ntot, rtnCode, start, zcase;
/* this function calculates the pairwise likelihood for the exceedandes */	
	xi = theta[0];
    sigma = theta[1];

	semiax1 = theta[2];
	semiax2 = theta[3];
	phi = theta[4];
	
	duration = theta[5];
	velocity[0] = theta[6]; // velocity
	velocity[1] = theta[7];
	
	alpha = theta[8];
	beta = theta[9];
	kappa = theta[10];
	volume = M_PI*semiax1*semiax2*duration; // volume of slated cylinder

    *count=0.0;
	if ( sigma <= 0) {
		*loglik = MINF*(1-sigma)*(1-sigma);		
			
		return;
	}
	
	if ( kappa <= 0) {
		*loglik = MINF*(1-kappa)*(1-kappa);		
		
		return;
	}
    
    if ( semiax1 <= 0) {
		*loglik = MINF*(1-semiax1)*(1-semiax1);		
		
		return;
	}
	
   if ( semiax2 <= 0) {
		*loglik = MINF*(1-semiax2)*(1-semiax2);		
		
		return;
	} 
  
   if ( duration <= 0) {
		*loglik = MINF*(1-duration)*(1-duration);		
		
		return;
	}
	
	if ( (phi < 0) || (phi >= M_PI_2) ) {
		*loglik = MINF*(1-phi)*(1-phi);		
		
		return;
	} 
	
	
     
	*loglik = 0;
	

    a[0] = alpha;
	ntot = (*nsites)*(*lag);
	
	for (k = 0; k < (*nblocks); k ++) {
		start = k*(*nsites);
		end = (k+1)*(*nsites);
		/*if (end > *nobs )*/
		for (i = start; i < end; i++) {
			if (!ISNA(y[i])) {
		    tmp = test(y[i],xi,sigma); /* test if the argument is positive */
						if ( tmp <= 0) {
							*loglik = MINF*(1-tmp)*(1-tmp);		
						return;
						}					
			
			for (j = (i+1); j < MIN((end+ntot),*nobs); j++){ 
				if (!ISNA(y[j])) {
				distspatial = dist(xcoords[i],ycoords[i],xcoords[j],ycoords[j]);
				if ( (distspatial <= *delta)) {
	/*	double slated_elliptic_cylinders_intersection (double V[2], double D, double PHI, double A, double B, 
                                double H1, double K1, double T1, double H2, double K2, double T2)		
                                */						
					tmp = slated_elliptic_cylinders_intersection(velocity, duration, phi, semiax1, semiax2, 
                                xcoords[i],ycoords[i],tcoords[i], xcoords[j],ycoords[j],tcoords[j]);                     	
					a[2] = alpha*tmp/volume; 
                    a[1] = a[3] = (a[0]-a[2]); 
                    zcase = z[i]  + 2 * z[j];
					*count+=1.0;						
	                switch(zcase) {
							case 0: /* # below -> below  */
									indloglik = 1-2*LT(a,beta,kappa) + LT2(a,beta,kappa,kappa); 
									*loglik += log(indloglik);
									break;		
							case 1: /* # above -> below */
									tmp = test(y[j],xi,sigma);
									if ( tmp <= 0) {
										*loglik = MINF*(1-tmp)*(1-tmp);		
										return;
									}							
									xx = tr(y[i],xi,sigma,kappa)+kappa;
									yy = kappa;
									indloglik = - dLT(a,beta,xx) + DxLT2(a,beta,xx,yy);
									jacx = dtr(y[i],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacx);
									break;
							case 2: /* # below -> above  */
									tmp = test(y[j],xi,sigma);
									if ( tmp <= 0) {
										*loglik = MINF*(1-tmp)*(1-tmp);		
										return;
									}					
									xx = kappa;
									yy = tr(y[j],xi,sigma,kappa)+kappa;	
									indloglik = - dLT(a,beta,yy) + DyLT2(a,beta,xx,yy);
									jacy = dtr(y[j],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacy);
									break;
							case 3:        /*# above -> above*/
									tmp = test(y[j],xi,sigma);
									if ( tmp <= 0) {
										*loglik = MINF*(1-tmp)*(1-tmp);		
										return;
									}					
									xx = tr(y[i],xi,sigma,kappa)+kappa; /* transform the data */
									yy = tr(y[j],xi,sigma,kappa)+kappa;
									
									indloglik = DxyLT2(a,beta,xx,yy);	
									
									jacx = dtr(y[i],  xi,  sigma, kappa);
									jacy = dtr(y[j],  xi,  sigma, kappa);
									*loglik += log(indloglik)+log(jacx)+log(jacy);
									break;
								}
							}
					}
				}
			}
		}		
	}	
}

/* version for parallel computing which returns a vector of individual pairwise lik*/
void  pwl_spatiotemporal_parallel_vec(double *y, int *z, int *nobs, int *nsites, int *nblocks, int *lag, double *theta, double *loglik, double *u, double *xcoords, 
                            double *ycoords, double *tcoords, double *delta, int *count)
{
	double  alpha, beta, a[4], distspatial,  duration, indloglik, kappa, jacx, jacy, phi, semiax1, semiax2, sigma, tmp,  velocity[2], volume,xi, xx, yy, X[4], Y[4];
	int end, i, j, k,  nroots, ntot, rtnCode, start, zcase;
    
/* this function calculates the pairwise likelihood for the exceedandes */	
	xi = theta[0];
    sigma = theta[1];

	semiax1 = theta[2];
	semiax2 = theta[3];
	phi = theta[4];
	
	duration = theta[5];
	velocity[0] = theta[6]; // velocity
	velocity[1] = theta[7];
	
	alpha = theta[8];
	beta = theta[9];
	kappa = theta[10];
	volume = M_PI*semiax1*semiax2*duration; // volume of slated cylinder
    
   
    
    *count =0;
	if ( sigma <= 0) {
		loglik[*count] = MINF*(1-sigma)*(1-sigma);	
		*count +=1;	
		return;
	}
	if ( kappa <= 0) {
		loglik[*count] = MINF*(1-kappa)*(1-kappa);		
		*count +=1;
		return;
	}
     
    if ( semiax1 <= 0) {
		loglik[*count] = MINF*(1-semiax1)*(1-semiax1);		
		*count +=1;
		return;
	}
	 
   if ( semiax2 <= 0) {
		loglik[*count] = MINF*(1-semiax2)*(1-semiax2);		
		*count +=1;
		return;
	} 
  
   if ( duration <= 0) {
		loglik[*count] = MINF*(1-duration)*(1-duration);		
		*count +=1;
		return;
	} 
	
	if ( (theta[4] < 0) || (theta[4] >= M_PI_2) ) {
		loglik[*count] = MINF*(1-theta[4])*(1-theta[4]);		
		*count +=1;
		return;
	}
	
   	 
			

    a[0] = alpha;
	ntot = (*nsites)*(*lag);
	for (k = 0; k < (*nblocks); k ++) {
		start = k*(*nsites);
		end = (k+1)*(*nsites);
		/*if (end > *nobs )*/
		for (i = start; i < end; i++) {
			if (!ISNA(y[i])) {
		    tmp = test(y[i],xi,sigma); /* test if the argument is positive */
						if ( tmp <= 0) {
							loglik[*count] = MINF*(1-tmp)*(1-tmp);		
							*count +=1;
						return;
						}					
			 
			for (j = (i+1); j < MIN((end+ntot),*nobs); j++){
				if (!ISNA(y[j])) {   
				distspatial = dist(xcoords[i],ycoords[i],xcoords[j],ycoords[j]);
				if ( distspatial <= *delta ) {
					a[2] = alpha*slated_elliptic_cylinders_intersection(velocity, duration, phi, semiax1, semiax2, 
                                xcoords[i],ycoords[i],tcoords[i], xcoords[j],ycoords[j],tcoords[j])/volume;
                    a[1] = a[3] = (a[0]-a[2]); 
                    zcase = z[i]  + 2 * z[j];
                    
                    switch(zcase) {
							case 3: /*# above -> above*/
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
								*count +=1;		
								return;
								}					
								xx = tr(y[i],xi,sigma,kappa)+kappa; /* transform the data */
								yy = tr(y[j],xi,sigma,kappa)+kappa;
								indloglik = DxyLT2(a,beta,xx,yy);	
								jacx = dtr(y[i],  xi,  sigma, kappa);
								jacy = dtr(y[j],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacx)+log(jacy);
								*count +=1;
								break;
							case 1: /* # above -> below */
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
									*count +=1;		
								return;
								}							
								xx = tr(y[i],xi,sigma,kappa)+kappa;
								yy = kappa;
								indloglik = - dLT(a,beta,xx) + DxLT2(a,beta,xx,yy);
								jacx = dtr(y[i],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacx);
								*count +=1;
								break;
							case 2:  /* # below -> above  */
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
									*count +=1;		
								return;
								}					
								xx = kappa;
								yy = tr(y[j],xi,sigma,kappa)+kappa;	
								indloglik = - dLT(a,beta,yy) + DyLT2(a,beta,xx,yy);
								jacy= dtr(y[j],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacy);
								*count +=1;
								break;
							case 0:  /* # below -> below  */
								indloglik = 1-2*LT(a,beta,kappa) + LT2(a,beta,kappa,kappa); 
								loglik[*count] = log(indloglik);
								*count +=1;
					}
				}
			}	
			}	
		}
	}
	}
}	

/*  returns a vector of individual pairwise lik*/
void  pwl_spatiotemporal_vec(double *y, int *z, int *nobs, int *nsites, int *lag, double *theta, double *loglik, double *u, double *xcoords, 
                            double *ycoords, double *tcoords, double *delta, int *count)
{
	double  alpha, beta, a[4], distspatial,  duration, indloglik, kappa, jacx, jacy, phi, semiax1, semiax2, sigma, tmp,  velocity[2], volume,xi, xx, yy, X[4], Y[4];
	int end, i, j, k,  nroots, ntot, rtnCode, start, zcase;
    
/* this function calculates the pairwise likelihood for the exceedandes */	
	xi = theta[0];
    sigma = theta[1];

	semiax1 = theta[2];
	semiax2 = theta[3];
	phi = theta[4];
	
	duration = theta[5];
	velocity[0] = theta[6]; // velocity
	velocity[1] = theta[7];
	
	alpha = theta[8];
	beta = theta[9];
	kappa = theta[10];
	volume = M_PI*semiax1*semiax2*duration; // volume of slated cylinder
    
   
    
    *count =0;
	if ( sigma <= 0) {
		loglik[*count] = MINF*(1-sigma)*(1-sigma);	
		*count +=1;	
		return;
	}
	if ( kappa <= 0) {
		loglik[*count] = MINF*(1-kappa)*(1-kappa);		
		*count +=1;
		return;
	}
     
    if ( semiax1 <= 0) {
		loglik[*count] = MINF*(1-semiax1)*(1-semiax1);		
		*count +=1;
		return;
	}
	 
   if ( semiax2 <= 0) {
		loglik[*count] = MINF*(1-semiax2)*(1-semiax2);		
		*count +=1;
		return;
	} 
  
   if ( duration <= 0) {
		loglik[*count] = MINF*(1-duration)*(1-duration);		
		*count +=1;
		return;
	} 
	
	if ( (theta[4] < 0) || (theta[4] >= M_PI_2) ) {
		loglik[*count] = MINF*(1-theta[4])*(1-theta[4]);		
		*count +=1;
		return;
	}
	
   	 
			

    a[0] = alpha;
	ntot = (*nsites)*(*lag);
	
	for (i = 0; i < *nobs; i++) {
	    tmp = test(y[i],xi,sigma); /* test if the argument is positive */
		if ( tmp <= 0) {
			loglik[*count] = MINF*(1-tmp)*(1-tmp);		
			*count +=1;
			return;
			}					
			 
		for (j = (i+1); j < *nobs; j++){   
			distspatial = dist(xcoords[i],ycoords[i],xcoords[j],ycoords[j]);
			if ( distspatial <= *delta ) {
				a[2] = alpha*slated_elliptic_cylinders_intersection(velocity, duration, phi, semiax1, semiax2, 
                                xcoords[i],ycoords[i],tcoords[i], xcoords[j],ycoords[j],tcoords[j])/volume;
                    a[1] = a[3] = (a[0]-a[2]); 
                    zcase = z[i]  + 2 * z[j];
                    
                    switch(zcase) {
							case 3: /*# above -> above*/
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
								*count +=1;		
								return;
								}					
								xx = tr(y[i],xi,sigma,kappa)+kappa; /* transform the data */
								yy = tr(y[j],xi,sigma,kappa)+kappa;
								indloglik = DxyLT2(a,beta,xx,yy);	
								jacx = dtr(y[i],  xi,  sigma, kappa);
								jacy = dtr(y[j],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacx)+log(jacy);
								*count +=1;
								break;
							case 1: /* # above -> below */
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
									*count +=1;		
								return;
								}							
								xx = tr(y[i],xi,sigma,kappa)+kappa;
								yy = kappa;
								indloglik = - dLT(a,beta,xx) + DxLT2(a,beta,xx,yy);
								jacx = dtr(y[i],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacx);
								*count +=1;
								break;
							case 2:  /* # below -> above  */
								tmp = test(y[j],xi,sigma);
								if ( tmp <= 0) {
									loglik[*count] = MINF*(1-tmp)*(1-tmp);
									*count +=1;		
								return;
								}					
								xx = kappa;
								yy = tr(y[j],xi,sigma,kappa)+kappa;	
								indloglik = - dLT(a,beta,yy) + DyLT2(a,beta,xx,yy);
								jacy= dtr(y[j],  xi,  sigma, kappa);
								loglik[*count] = log(indloglik)+log(jacy);
								*count +=1;
								break;
							case 0:  /* # below -> below  */
								indloglik = 1-2*LT(a,beta,kappa) + LT2(a,beta,kappa,kappa); 
								loglik[*count] = log(indloglik);
								*count +=1;
					}
				}
			}	
		}
	
}	




