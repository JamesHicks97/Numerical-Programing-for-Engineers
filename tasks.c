/***************************************************************************
 *
 *   File        : tasks.c
 *   Student Id  : 833966
 *   Name        : James Hicks
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

#define MAXITERS 25 // Used in newton-raphson function. Usually raches a 
					// solution with 5 iterations so 25 is a generous number.
#define EPSILON 1e-6 // Accuracy factor to the 6th decimal place.
#define SQUARE 2 // Not gonna lie was just scared of magic numbers.
#define CHARBUFFER 30 // Allocates memory for headers.
#define DATABUFFER 6 // Allocates memory for data
#define PI 3.14159265359
#define XMAX 1.0
#define TIMEMAX 0.2

// Calculates value for f: Task 1
double f(double m, double theta, double gamma, double b);

// Calculates value for fprime: Task 1
double fprime(double m, double gamma, double b);

// Performs Newton-Raphson root finding method: Task 2
double  newtonraph(double m, double theta, double gamma, double b);

// converts from degrees to radians: Task 2
double degreetorad(double degree);

// converts for radians to degrees:  Task 2
double  radtodegree(double rad);

// Obtains intial guess for Newton-Raphson for lower root for a particular m:  Task 2
double lowguess(double m);

// Uses Thomas Algorithm for solving tridiagonal matrix:  Task 4 and 5
void solvetriag(double* A, double* B, double* C, double* Q, double* X, int n);

// Fills tridiagonal matrix needed to find C values in spline calculations.
void filltriag(double* TRIDIAG1, double* TRIDIAG2, double* TRIDIAG3, double* TRIDIAG4, double* H, double* X, double* A, int n);

// 2nd order Runge-kutta function
void rk2_C(double *F, double *F_INTER, double *F_NEXT, double delta_x, double delta_t, double c, int nx);

void rk2_U(double *F, double *F_INTER, double *F_NEXT, double delta_x, double delta_t, double c, int nx);

double central(int i, double delta_x, double *F, double c, int nx);

double upwind(int i, double delta_x, double *F, double c);

void initialf(double nx, double *F, double delta_x);
//
void shockwave(const char* q2_file)
{

	char header[CHARBUFFER];
	int thetaint;
	double m, theta, blguess, bl, buguess, bu, gamma;
	FILE *q2_read; 
	FILE *q2_write;
	q2_read = fopen(q2_file, "r");
	q2_write = fopen("out_shock.csv", "w");
	
	// Part A
	//Read first two lines used for part A.
	fscanf(q2_read, "%s", header);
	fscanf(q2_read, "%lf,%lf,%lf,%lf,%lf\n", 
		&m ,&theta ,&blguess ,&buguess ,&gamma);
	
	/*
	bl = blguess;
	bu = buguess; 
	
	bl = degreetorad(bl);
	bu = degreetorad(bu);
	theta = degreetorad(theta);
	
	bl = newtonraph(m, theta, gamma, bl);
    bu = newtonraph(m, theta, gamma, bu);
    
    bl = radtodegree(bl);
    bu = radtodegree(bu);
    
    printf("bl = %0.6f, bu = %0.6f \n", bl, bu);
    */
   
    //Part B.
    //Read two lines used for part B
    fscanf(q2_read, "%s", header);
	fscanf(q2_read, "%lf\n", &m);
	
	/*fprintf(q2write, "m, theta, bl, bu\n");
	
	blguess=radtodegree(lowguess(m));
	
	for (theta = 0; theta <= 90 ; theta++) {
		
		
		bl = blguess;
		bu = buguess; 
		
		bl = degreetorad(bl);
		bu = degreetorad(bu);
		theta = degreetorad(theta);
		
		bl = newtonraph(m, theta, gamma, bl);
		if (bl == EXIT_FAILURE) {
				break;
			}
		bu = newtonraph(m, theta, gamma, bu);
		
		bl = radtodegree(bl);
		bu = radtodegree(bu);
		theta = radtodegree(theta);
		
		fprintf(q2_write, " %lf, %lf, %lf\n",  theta, bl, bu);
		
	}
	*/
	
	// Part C
	// Read head for part c and print head for output file.
	fscanf(q2_read, "%c", header);
	fprintf(q2_write, "M, theta, beta_lower, beta_upper\n");
	
	// Take one M at a time and perform the necessary operations.
	while (fscanf(q2_read, "%lf", &m) != EOF) {
		
		// Use the lowest possible value of bl (when theta = 0) and set this
		// as the first guess for bl. Bu will always be 90 degrees.
		blguess = radtodegree(lowguess(m));
		
		//Cycle through the theta values we want to evaluate for.
		for (theta = 0 ; theta <=90 ; theta++) {
			
			// Convert values from degrees to radians to work in newton raphson
			// function.
			bl = degreetorad(blguess);
			bu = degreetorad(buguess);
			theta = degreetorad(theta);
			
			// Find a solution to bl.
			bl = newtonraph(m, theta, gamma, bl);
			// If no solution is found within the maximum iterations assume 
			// there is no solution and theta max had been reached. Break the 
			// theta loop and scan the next m value.
			if (bl == EXIT_FAILURE) {
				break;
			}
			// If a solution was found for bl find solution for bu.
			bu = newtonraph(m, theta, gamma, bu);
			
			// Convert bl, bu and theta back to degrees.
			bl = radtodegree(bl);
			bu = radtodegree(bu);
			theta = radtodegree(theta);
			thetaint = (int)theta;
			// Print to out file.
			fprintf(q2_write, "%0.6lf, %d, %0.6lf, %0.6lf\n", m, thetaint, bl, bu);
			
			}
			
		}
		
	
	// Close read and write files.
    fclose(q2_read);
    fclose(q2_write);
    
}



void linalgbsys(const char* q4_file)
{
   FILE *q4_read; 
   FILE *q4_write;
   q4_read = fopen(q4_file, "r");
   q4_write = fopen("out_linalsys.csv", "w");
   int i = 0,databuffer = DATABUFFER, n;
   double *A, *B, *C, *Q, *X;
   
   
   char header[CHARBUFFER];
   fscanf(q4_read, "%s", header);
   
   // Inital allocation of memory for each array.
   A = (double*)malloc(DATABUFFER*sizeof(double));
   B = (double*)malloc(DATABUFFER*sizeof(double));
   C = (double*)malloc(DATABUFFER*sizeof(double));
   Q = (double*)malloc(DATABUFFER*sizeof(double));
   X = (double*)malloc(DATABUFFER*sizeof(double));
   
   while (fscanf(q4_read, "%lf,%lf,%lf,%lf\n", &A[i], &B[i], &C[i], &Q[i]) != EOF) {
   	   
   	   i++;
   	   // Check if the allocated memory if full and if so reaccocates another databuffer.
   	   if (i == databuffer) {
   	   	   databuffer += DATABUFFER;	
   	   	   A = (double*)realloc(A, databuffer * sizeof(double));
   	   	   B = (double*)realloc(B, databuffer * sizeof(double));
   	   	   C = (double*)realloc(C, databuffer * sizeof(double));
   	   	   Q = (double*)realloc(Q, databuffer * sizeof(double));
   	   	   X = (double*)realloc(X, databuffer * sizeof(double));
   	   }
   	   
   }
   
   // Sets n to the number of data entries.
   n = i; 
   
   
   solvetriag(A, B, C, Q, X, n);
   
   //Cycle through the solution array and print solutions to the out file.
   for(i=0; i<n; i++) {
   	   fprintf(q4_write,"%.6lf\n", X[i]);
   	   
   }
   
   
   fclose(q4_read);
   fclose(q4_write);
   
   // Free the memory allocated to the arrays.
   free(A);
   free(B);
   free(C);
   free(Q);
   free(X);
}

void interp(const char* q5_file, const double xo)
{
    FILE *q5_read; 
	FILE *q5_write;
	q5_read = fopen(q5_file, "r");
	q5_write = fopen("out_interp.csv", "w");
	
	char header[CHARBUFFER];
	int i=0, databuffer=DATABUFFER, n;
	double *X, *A, *B, *C, *D, *H;
	double *TRIDIAG1, *TRIDIAG2, *TRIDIAG3, *TRIDIAG4;
	double *S;
	
	
	
	
	
	X = (double*)malloc(DATABUFFER*sizeof(double));
	A = (double*)malloc(DATABUFFER*sizeof(double));
	B = (double*)malloc(DATABUFFER*sizeof(double));
	C = (double*)malloc(DATABUFFER*sizeof(double));
	D = (double*)malloc(DATABUFFER*sizeof(double));
	H = (double*)malloc(DATABUFFER*sizeof(double));
	TRIDIAG1 = (double*)malloc(DATABUFFER*sizeof(double));
	TRIDIAG2 = (double*)malloc(DATABUFFER*sizeof(double));
	TRIDIAG3 = (double*)malloc(DATABUFFER*sizeof(double));
	TRIDIAG4 = (double*)malloc(DATABUFFER*sizeof(double));
	S = (double*)malloc(DATABUFFER*sizeof(double));
	
	
	fscanf(q5_read, "%s", header);
	while(fscanf(q5_read, "%lf,%lf", &X[i], &A[i]) != EOF) {
		i++;
   	   if (i == databuffer) {
   	   	   databuffer += DATABUFFER;
   	   	   X = (double*)realloc(X, databuffer * sizeof(double));
   	   	   A = (double*)realloc(A, databuffer * sizeof(double));
   	   	   
   	   	   
   	   }
   	}
   n = i;
   	
   	B = (double*)realloc(B, n * sizeof(double));
   	C = (double*)realloc(C, n * sizeof(double));
   	D = (double*)realloc(D, n * sizeof(double));
   	H = (double*)realloc(H, n * sizeof(double));
   	TRIDIAG1 = (double*)realloc(TRIDIAG1, n * sizeof(double));
   	TRIDIAG2 = (double*)realloc(TRIDIAG2, n * sizeof(double));
   	TRIDIAG3 = (double*)realloc(TRIDIAG3, n * sizeof(double));
   	TRIDIAG4 = (double*)realloc(TRIDIAG4, n * sizeof(double));
   	S = (double*)realloc(S, n * sizeof(double));
   	
   	
   	filltriag(TRIDIAG1, TRIDIAG2, TRIDIAG3, TRIDIAG4, H, X, A, n);
   	solvetriag(TRIDIAG2, TRIDIAG3, TRIDIAG1, TRIDIAG4, C, n);
   	
   	// Solve B and D arrays.
   	for (i=0; i<n-1 ; i++) {
   		B[i]=(1.0/H[i])*(A[i+1]-A[i])- (H[i]/3.0)*(2.0*C[i]+C[i+1]);
   		D[i]=(C[i+1]-C[i])/(3.0 * H[i]);
   		
   	}
   
   	// Code used to find spline graph
   	/*
   	for (i=0; i<n-1; i++) {
   		if (X[i]>X[i+1]){
   			for (x=X[i]; x > X[i+1]; x -= 0.00001){
   			
			S[i] = A[i] + B[i]*(x-X[i])+C[i]*pow((x-X[i]),2)+ D[i]*pow((x-X[i]),3);
			fprintf(q5_write, "%lf, %lf\n", x, S[i]);
			}
		}else if (X[i]<X[i+1]) {
			for (x=X[i]; x < X[i+1]; x += 0.00001){
   			
			S[i] = A[i] + B[i]*(x-X[i])+C[i]*pow((x-X[i]),2)+ D[i]*pow((x-X[i]),3);
			fprintf(q5_write, "%lf, %lf\n", x, S[i]);
			}
		}
	}
	*/
	fprintf(q5_write, "xo,f(xo)\n");
	// Cycle through the data points and finds what spline xo lies in
	for (i=0; i<n-1; i++) {
		//If xo is inbetween the point being assed and the next point, use
		// that points spline to find S[i]. 
		if ((xo<X[i] && xo> X[i+1]) || (xo>X[i] && xo< X[i+1])){
			fprintf(q5_write, "%lf,%lf\n", xo, A[i] + B[i]*(xo-X[i])+C[i]*pow((xo-X[i]),2)+ D[i]*pow((xo-X[i]),3));
		}
	}
	
	free(X); 
	free(A); 
	free(B); 
	free(C); 
	free(D); 
	free(H); 
	free(TRIDIAG1); 
	free(TRIDIAG2);
	free(TRIDIAG3);
	free(TRIDIAG4);
	free(S);
	
	fclose(q5_read);
	fclose(q5_write);	
	
}


void waveeqn(const char* q6_file)
{
    FILE *q6_read; 
	FILE *q6_write1;
	FILE *q6_write2;
	q6_read = fopen(q6_file, "r");
	q6_write1 = fopen("out_waveeqn_1U.csv", "w");
	q6_write2 = fopen("out_waveeqn_2C.csv", "w");
	
	char header[CHARBUFFER];
	int nx, out_iter, timesteps, i, k;
	double c, cfl, delta_x, delta_t, x /*,t*/;
	double *F, *F_NEXT, *F_INTER;

	
	fscanf(q6_read, "%s", header);
	fscanf(q6_read, "%lf,%d,%lf,%d", &c, &nx, &cfl, &out_iter);
	
	delta_x = XMAX/nx;
	delta_t = cfl*delta_x/c;
	timesteps = TIMEMAX/delta_t+1;
	
	F = (double*)malloc((nx+1)*sizeof(double));
	F_INTER = (double*)malloc((nx+1)*sizeof(double));
	F_NEXT = (double*)malloc((nx+1)*sizeof(double));
	
	//6.1
	//initial f values
	/*for (i=0; i<=nx+1; i++){
			x = (double)i*delta_x;
			F[i] = initialf(x);
			fprintf(q6_write, "%lf,%lf\n",i, delta_x, x, F[i]);
	}*/
	
	//Using upwind scheme
	initialf(nx, F, delta_x);
	
	for (k=0; k< out_iter; k++) {
		
	rk2_U(F, F_INTER, F_NEXT, delta_x, delta_t, c, nx);
			
	}
	fprintf(q6_write1, "x,f(x)\n");
	for (i=0; i<=nx; i++){
		x = (double)i*delta_x;
		fprintf(q6_write1, "%lf,%lf\n", x, F[i]);
	}
	
	//Using central scheme
	initialf(nx, F, delta_x);
	
	for (k=0; k< out_iter; k++) {
		
		rk2_C(F, F_INTER, F_NEXT, delta_x, delta_t, c, nx);		
	}
	
	fprintf(q6_write2, "x,f(x)\n");
	for (i=0; i<=nx; i++){
		x = (double)i*delta_x;
		fprintf(q6_write2, "%lf,%lf\n", x, F[i]);
	}
	/*
	// timeloop
	initialf(nx, F, delta_x);
	for (int n = 0; n <= timesteps; n ++) {     
         t = n*delta_t;
         if (t==0.05 || t==0.1 || t==0.15 || t==0.2){
         	 fprintf(q6_write2, "x,f(x), t=%lf, nx= %d, \n", t, nx);
         	 for (i=0; i<=nx; i++){
         	 	 x = (double)i*delta_x;
         	 	 fprintf(q6_write2, "%lf,%lf\n", x+(c*t), F[i]);
         	 }
         }
    }
	for (int n = 0; n <= timesteps; n ++) {     
         t = n*delta_t;
         rk2_U(F, F_INTER, F_NEXT, delta_x, delta_t, c, nx);
         if (t==0.05 || t==0.1 || t==0.15 || t==0.2){
         	 fprintf(q6_write2, "x,f(x), t=%lf, nx= %d, \n", t, nx);
         	 for (i=0; i<=nx; i++){
         	 	 x = (double)i*delta_x;
         	 	 fprintf(q6_write2, "%lf,%lf\n", x, F[i]);
         	 }
         }
    }
    */
	
	free(F);
	free(F_INTER);
	free(F_NEXT);
	fclose(q6_read);
	fclose(q6_write1);
	fclose(q6_write2);
}


double f(double m, double theta, double gamma, double b) {
	
	double fx;
	// Finds value of f(x) using equation (2) from assignment.
	fx = 2*(cos(b)/sin(b))*(pow(m,SQUARE)*pow(sin(b),SQUARE)-1)/
	(pow(m,SQUARE)*(gamma+cos(2*b))+2) - tan(theta);
	
	return fx;
}

double fprime(double m, double gamma, double b) {
	double fprimex, part1, part2, part3;
	
	// Finds value of fprime(x) using equation derived equation (2) from assignment.
	// Equation split into parts for ease of use.
	part1 = 4 * pow(m,SQUARE) * pow(cos(b),SQUARE) / 
			(pow(m,SQUARE) * (cos(2*b) + gamma) + 2);
			
	part2 = 4 * pow(m,SQUARE) * sin(2*b) * (cos(b)/sin(b)) * 
			(pow(m,SQUARE) * pow(sin(b),SQUARE) - 1 ) /
			pow(pow(m,SQUARE) * (gamma+cos(2*b)) + 2 ,SQUARE);
			
	part3 = 2 * pow(1/sin(b),SQUARE) * 
			(pow(m,SQUARE) * pow(sin(b),SQUARE)-1) /
			(pow(m,SQUARE) * (cos(2*b) + gamma) + 2);
			
	fprimex = part1 + part2 - part3; 
	return fprimex;
}

double newtonraph(double m, double theta, double gamma, double b) {
	int iters=0;
	
	// Checks if f(b) is within the error margin of zero
	while (fabs(f(m, theta, gamma, b)) > EPSILON ) {
			// If it is not it finds a new guess for b. increase the number of
			// iterations by one
			b = b - (f(m, theta, gamma, b)/fprime(m, gamma, b));
			iters++;
			// If the max iterations is reached there is most like no solution 
			// So newton-raphson fails and returns failure status.
			if (iters == MAXITERS) {
				return EXIT_FAILURE;
			}
	}
	// if f(b) is within the error margin of zero it returns b
	return b;
}


double  degreetorad(double degree) {
	
	double rad;
	
	rad = degree * PI / 180.0 ;
	
	return rad;
	
}

double radtodegree(double rad) {
	
	double degree;
	
	degree = rad * 180.0 / PI ; 
	
	return degree;
	
}

double lowguess( double m) {
	
	double bl;
	// The minimum bl occurs when theta = 0. Finds an intial low guess for an
	// m using eq (3) from assignment. bu guess will always be the same.
	bl= asin(1/m);
	return bl;
}


void solvetriag(double* A, double* B, double* C, double* Q, double* X, int n) {
	
	int i;
	
	//Cycles through the lines of the "matrix" to find new values using the 
	// THomas algorithm.
	for(i=0; i<n; i++) {
		//First line of matrix does not need to be changed.
   	   if (i == 0){
   	   	   continue;
   	   }
   	   else {
   	   	   A[i] = A[i] - B[i-1]*C[i]/A[i-1];
   	   	   Q[i] = Q[i] - Q[i-1]*C[i]/A[i-1];
   	   } 
   }
   // Backward solves to find the solution array.
   for(i=n-1; i>=0; i--) {
   	   if (i == n){
   	   	   X[i]= Q[i]/A[i];
   	   }
   	   else {
   	   	   X[i]= (Q[i]-B[i]*X[i+1])/A[i];
   	   }
   }
}

void filltriag(double* TRIDIAG1, double* TRIDIAG2, double* TRIDIAG3, double* TRIDIAG4, double* H, double* X, double* A, int n) {
	
	int i;
	// Fills the arrays needed to find the solutions to the C array.
	//TRIDIAG1 are the c values in solving a tridiagonal matrix, TRIDIAG2 are 
	//the a's, TRIDIAG3 are the b's and TRIDIAG4 are the q's
	for (i=0; i<n ; i++) {
   		if (i == 0) {
   			H[i]=X[i+1]-X[i];
   			TRIDIAG1[i] = 0;
   			TRIDIAG2[i] = 1;
   			TRIDIAG3[i] = 0;
   			TRIDIAG4[i] = 0;
   		} else if (i == n-1) {
   			H[i]=X[i+1]-X[i];
   			TRIDIAG1[i] = 0;
   			TRIDIAG2[i] = 1;
   			TRIDIAG3[i] = 0;
   			TRIDIAG4[i] = 0;
   		} else { 
   			H[i]=X[i+1]-X[i];
   			TRIDIAG1[i] = H[i-1];
   			TRIDIAG2[i] = 2.0*(H[i-1]+H[i]);
   			TRIDIAG3[i] = H[i];
   			TRIDIAG4[i]= (3.0/H[i])*(A[i+1]-A[i]) + (3.0/H[i-1])*(A[i-1]-A[i]);
   		}
   	}
 }
 

void rk2_U(double *F, double *F_INTER, double *F_NEXT, double delta_x, double delta_t, double c, int nx ) {
 	int i;
 	
	for (i=0; i<=nx; i++){
	F_INTER[i] = F[i]+(delta_t*upwind(i, delta_x, F, c));

	}
	for (i=0; i<=nx; i++){
	F_NEXT[i] = F[i] + ((1.0/2.0)*delta_t *
	(upwind(i, delta_x, F, c)+ upwind(i, delta_x, F_INTER, c)));
	}
	for (i=0; i<=nx; i++){
	F[i] = F_NEXT[i];
	}
}
 
void rk2_C(double *F, double *F_INTER, double *F_NEXT, double delta_x, double delta_t, double c, int nx) {
 	
	int i;
	
	for (i=0; i<=nx; i++){
	F_INTER[i] = F[i]+(delta_t*central(i, delta_x, F, c, nx));

	}
	for (i=0; i<=nx; i++){
	F_NEXT[i] = F[i] + ((1.0/2.0)*delta_t *
	(upwind(i, delta_x, F, c)+ central(i, delta_x, F_INTER, c, nx)));
	}
	for (i=0; i<=nx; i++){
	F[i] = F_NEXT[i];
	}
}
 	 
 
double upwind(int i, double delta_x, double *F, double c) {
	double rhs;
 	if (i == 0) {
 		rhs= -c*(F[1]-F[0])/delta_x;
 	} else {
 		rhs= -c*(F[i]-F[i-1])/delta_x;
 	}
 	
 	return rhs;
}
 
double central(int i, double delta_x, double *F, double c, int nx) {
 	double rhs;
 	 
 	if (i == 0) {
 	 	rhs= -c*(F[1]-F[0])/delta_x;
 	} else if (i == nx) {
 		rhs= -c*(F[nx]-F[nx-1])/delta_x;
	} else {
 	 	rhs= -c*(F[i+1]-F[i-1])/(2*delta_x);
 	}
 	return rhs;
}	 
 	 
 	 
void initialf(double nx, double *F, double delta_x) {
	int i;
	double x;
 	for (i=0; i<=nx; i++){
		x = (double)i*delta_x;
		if (x>=0.125 && x<=0.375){
			F[i]= 0.5*(1-cos(8*PI*(x-0.125)));				
		} else { 
			F[i]=0;
		}
 	}
}

 	 	 
 	 	 