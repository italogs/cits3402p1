/* ilmsfxd.c: Quartic potential, Fixed boundary conditions */
/* This program can be continued from where it left off */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


#define chainlngth 2000
#define nsprngs (chainlngth+1)
#define halfchain 50 /* Ignore; not used here */
#define nmode 100  /* Ignore; not used here */
#define nprntstps 100 /* Number of output lines
							/* Also, t_end here because output is at t = 1, 2 etc. */ 
#define dt  0.001 /* Time step #CHANGED */
#define pi 3.14159 /* Ignore; not used here */
#define beta 0.7 /* Beware! beta is the nonlinear coefficient! */
		/* Usually alpha and beta appaears interchanged in literature */
#define alpha .16 /* alpha is the coefficient of the linear term! */
#define spcng 1.0 /* Just to make it explicit; not used in the program */

// An extra parameter was defined to define how many threads I'm using on this program 
#define nThreads 4

double accel(double *, double *);
double initialize(double *); /* Ignore; not used here */
double normalmode(int ,double *, double *, double *);/* Ignore; not used here */


int main() {
    struct timeval start, end;
    double delta;
	double a = alpha;
	double b = beta;
	printf("Alpha is:  %lf \n",  a);
	printf("Beta is :  %lf \n", b);


	FILE *fp, *fp1, *fp2, *fp3, *fp5, *fp6, *fp7, *fp8;
	int i, j, k, n, n1, prncnt, prntstps;
	double v[chainlngth], x[chainlngth], tke, tpe, te, ke1, omegak[nmode];
	double acc[chainlngth], ke[chainlngth], pe[nsprngs], fac, y[chainlngth];
	double t, t1, dx, hdt, hdt2,  twopi, twopisqr, alphaby4, cmass, cmom;

	prntstps = (int) (1.0 / dt);

	// The programmer might want to know the main parameters
	printf("INPUTS:\nchainlngth:%d\ndt:%.8lf\nnprntstps:%d\nnThreads:%d\nprntstps:%d\n",chainlngth,dt,nprntstps,nThreads,prntstps);

    gettimeofday(&start, NULL);

	alphaby4 = beta / 4.0;
	fp = fopen("updated/toten.dat","w");
	fp1 = fopen("updated/strsh.dat","w");
	fp3 = fopen("updated/velsh.dat","w");
	fp5 = fopen("updated/cmass.dat","w");
	fp2 = fopen("updated/restart.dat","w");
	fp6 = fopen("updated/ke.dat","w");
	fp7 = fopen("updated/acce.dat","w");
	fp8 = fopen("updated/pe.dat","w");



	/* Initialize the position, velocity, acceleration arrays */
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (i=0; i < chainlngth; i++) { 
			x[i] = acc[i] = v[i] = 0.0;
		}
	}

	/*
	Initial perturbation at the center of the chain and it is a double particle perturbation (EVEN Parity)
	*/
	x[50] = -0.98;
	x[51] = +0.98;	// Even Parity

	dx = tke = tpe = te = 0.0;
	for (i=0; i < chainlngth; i++) { 
		ke[i] = 0.5 * v[i] * v[i];
		fprintf(fp6,"%.10f\t", ke[i]);
		tke += ke[i];
		j = i-1;
		if (j == -1) { 
			dx = x[i];
		} else {
			dx = x[i] - x[j]; 
		}
		fac = dx * dx;
		pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;
		fprintf(fp8,"%.10f\t", pe[i]);
		tpe += pe[i];
	}
	dx = -x[chainlngth - 1];
	fac = dx * dx;
	pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;
	fprintf(fp8,"%.10f\n", pe[i]);
	tpe += pe[i];
	fprintf(fp6,"\n");
	te = tpe + tke; i = 0;
	fprintf(fp,"%d\t%.10f\n", i, te);


	#pragma omp parallel private(k)
	{
		/*
			As we want to keep the written data in the same order as the original output(procedural code output), 
			each thread should deal with only one file to avoid misplacement of data.
		*/
		#pragma omp sections 
		{
			#pragma omp section
			{
				for (k=0; k < chainlngth; k++) 
					fprintf(fp1,"%.10f\t", x[k]); 
				fprintf(fp1,"\n"); 
			}
			#pragma omp section
			{
				for (k=0; k < chainlngth; k++) 
					fprintf(fp3,"%.10f\t", v[k]); 
				fprintf(fp3,"\n"); 
			}
			#pragma omp section
			{
				for (k=0; k < chainlngth; k++) 
					fprintf(fp7,"%.10f\t", acc[k]); 
				fprintf(fp7,"\n"); 
			}
		}	
	}

	hdt = 0.5 * dt;
	hdt2 = dt * hdt; 
	n = 1; n1 = 1;

	while (n < nprntstps) { 

		/* new positions and mid-velocities; velocity-Verlet algorithm  */
		#pragma omp parallel
		{
			#pragma omp for schedule(static)
			for (j = 0; j < chainlngth; j++) {
				x[j] += dt * v[j] + hdt2 * acc[j];
				v[j] += hdt * acc[j];
			}
		}

		/* new accelerations */     
		accel(x, acc);


		/* new final velocities; Ignore the variables cmom and cmass */	    
		cmom = 0.0;
		#pragma omp parallel
		{
			#pragma omp for schedule(static) reduction(+:cmom)
			for (j = 0; j < chainlngth; j++) { 	
				v[j] += hdt * acc[j]; 
				cmom += v[j];
			}
			
		}
		cmom /= chainlngth;
		

		/* Kinetic energies */
		prncnt = n1 / prntstps;  //percent completion
		
		if (prncnt == 1) {  //if 100% done, print all
			tke = tpe = te = dx = 0.0;  //reset all variables
			cmass = 0.0;
			for (j = 0; j < chainlngth; j++) {
				ke[j] = 0.5 * v[j] * v[j]; 
				fprintf(fp6,"%.10f\t",ke[j]);
				tke += ke[j];
				k = j-1;
				if (k == -1) { 
					dx = x[j];
				} else {
					dx = x[j] - x[k]; 
				}
				fac = dx * dx;
				pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;
				fprintf(fp8,"%.10f\t", pe[i]);
				tpe += pe[i];
				cmass += x[j];
			}
			dx = -x[chainlngth - 1];
			fac = dx * dx;
			pe[i] = alpha * 0.5 * fac + alphaby4 * fac * fac;

			fprintf(fp8,"%.10f\n", pe[i]);
			tpe += pe[i];
			fprintf(fp5, "%d\t%.10f\n", i, cmass);
			cmass /= chainlngth;
			fprintf(fp6,"\n");
			te = tpe + tke;
			fprintf(fp,"%d\t%.10f\n", n, te);


			 #pragma omp parallel private(k) 
			{
				#pragma omp sections 
				{
					#pragma omp section
					{
						for (k=0; k < chainlngth; k++) {
							y[k] = x[k] - cmass;
							fprintf(fp1,"%.10f\t", y[k]); 
						}
						fprintf(fp1,"\n"); 
					}
					#pragma omp section 
					{
						for (k=0; k < chainlngth; k++) {
							fprintf(fp3,"%.10f\t", v[k]); 
						}
						fprintf(fp3,"\n"); 
					}
					#pragma omp section 
					{
						for (k=0; k < chainlngth; k++) {
							fprintf(fp7,"%.10f\t", acc[k]); 
						}
						fprintf(fp7,"\n");
					}
				}
			} 
			n++; n1 = 1;
		}
		n1++;
	}
	
	fclose(fp); fclose(fp1); fclose(fp3); fclose(fp5); 
	fclose(fp6); fclose(fp7); fclose(fp8);
  
	fprintf(fp2, "%d\t%d\n", i-1, n-1);
	for (j = 0; j < chainlngth; j++) 
		fprintf(fp2, "%.15f\t%.15f\t%.15f\n", x[j], v[j], acc[j]);
		fclose(fp2); 
        gettimeofday(&end, NULL);
        delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         end.tv_usec - start.tv_usec) / 1.e6;
        printf("Total time=%f seconds\n\n", delta);

        return 0;
}

double accel(double *x, double *acc) {

	double dxi, twodxi, dxipl1, dximn1, fac, fac1, fac2, fac13, fac23;
	int i, j, k;
	
	#pragma omp parallel private(k,j,dxi, twodxi, dxipl1, dximn1, fac, fac1, fac2, fac13, fac23)  num_threads(nThreads)
	{
		#pragma omp for schedule(auto)
		for (i = 0; i < chainlngth; i++) {
			j = i - 1;
			k = i + 1;
			if (j == -1) {
				dximn1 = 0.0;
			} else {
				dximn1 = x[j];
			}
			if (k == chainlngth) {
				dxipl1 = 0.0;
			} else {
				dxipl1 = x[k];
			}
			dxi = x[i];
			twodxi = 2.0 * dxi;
			fac = dxipl1 + dximn1 - twodxi;
			fac1 = dxipl1 - dxi;
			fac2 = dxi - dximn1;
			fac13 = fac1 * fac1 * fac1;
			fac23 = fac2 * fac2 * fac2;
			acc[i] = alpha * fac + beta * (fac13 - fac23);
		}
	}
}

