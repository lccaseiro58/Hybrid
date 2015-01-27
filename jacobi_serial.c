/************************************************************
 * Program to solve a finite difference
 * discretization of the screened Poisson equation:
 * (d2/dx2)u + (d2/dy2)u - alpha u = f
 * with zero Dirichlet boundary condition using the iterative
 * Jacobi method with overrelaxation.
 *
 * RHS (source) function
 *   f(x,y) = -alpha*(1-x^2)(1-y^2)-2*[(1-x^2)+(1-y^2)]
 *
 * Analytical solution to the PDE
 *   u(x,y) = (1-x^2)(1-y^2)
 *
 * Current Version: Christian Iwainsky, RWTH Aachen University
 * MPI C Version: Christian Terboven, RWTH Aachen University, 2006
 * MPI Fortran Version: Dieter an Mey, RWTH Aachen University, 1999 - 2005
 * Modified: Sanjiv Shah,        Kuck and Associates, Inc. (KAI), 1998
 * Author:   Joseph Robicheaux,  Kuck and Associates, Inc. (KAI), 1998
 *
 * Unless READ_INPUT is defined, a meaningful input dataset is used (CT).
 *
 * Input : n     - grid dimension in x direction
 *         m     - grid dimension in y direction
 *         alpha - constant (always greater than 0.0)
 *         tol   - error tolerance for the iterative solver
 *         relax - Successice Overrelaxation parameter
 *         mits  - maximum iterations for the iterative solver
 *
 * On output
 *       : u(n,m)       - Dependent variable (solution)
 *       : f(n,m,alpha) - Right hand side function
 *
 *************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/*************************************************************
 * Performs one iteration of the Jacobi method and computes
 * the residual value.
 *
 * NOTE: u(0,*), u(maxXCount-1,*), u(*,0) and u(*,maxYCount-1)
 * are BOUNDARIES and therefore not part of the solution.
 *************************************************************/
double one_jacobi_iteration(double xStart, double yStart,
             int maxXCount, int maxYCount,
             double *src, double *dst,
             double deltaX, double deltaY,
             double alpha, double omega)
{
#define SRC(XX,YY) src[(YY)*maxXCount+(XX)]
#define DST(XX,YY) dst[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    // Coefficients
    double cx = 1.0/(deltaX*deltaX);
    double cy = 1.0/(deltaY*deltaY);
    double cc = -2.0*cx-2.0*cy-alpha;

    for (y = 1; y < (maxYCount-1); y++)
    {
   fY = yStart + (y-1)*deltaY;
   for (x = 1; x < (maxXCount-1); x++)
   {
       fX = xStart + (x-1)*deltaX;
       f = -alpha*(1.0-fX*fX)*(1.0-fY*fY)-2.0*(1.0-fX*fX)-2.0*(1.0-fY*fY);
       updateVal = ((SRC(x-1,y)+SRC(x+1,y))*cx +
          (SRC(x,y-1)+SRC(x,y+1))*cy +
          SRC(x,y)*cc - f)/cc;
       DST(x,y) = SRC(x,y) - omega*updateVal;
       error += updateVal*updateVal;
   }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}


/**********************************************************
 * Checks the error between numerical and exact solutions
 **********************************************************/
double checkSolution(double xStart, double yStart,
           int maxXCount, int maxYCount,
           double *u,
           double deltaX, double deltaY,
           double alpha)
{
#define U(XX,YY) u[(YY)*maxXCount+(XX)]
    int x, y;
    double fX, fY;
    double localError, error = 0.0;

    for (y = 1; y < (maxYCount-1); y++)
    {
   fY = yStart + (y-1)*deltaY;
   for (x = 1; x < (maxXCount-1); x++)
   {
       fX = xStart + (x-1)*deltaX;
       localError = U(x,y) - (1.0-fX*fX)*(1.0-fY*fY);
       error += localError*localError;
   }
    }
    return sqrt(error)/((maxXCount-2)*(maxYCount-2));
}


int main(int argc, char **argv)
{
    int n, m, mits;
    double alpha, tol, relax;
    double maxAcceptableError;
    double error;
    double *u, *u_old, *tmp;
    int allocCount;
    int iterationCount, maxIterationCount;

#ifdef READ_INPUT
    printf("Input n,m - grid dimension in x,y direction:\n");
    scanf("%d,%d", &n, &m);
    printf("Input alpha - Helmholts constant:\n");
    scanf("%lf", &alpha);
    printf("Input relax - Successive over-relaxation parameter:\n");
    scanf("%lf", &relax);
    printf("Input tol - error tolerance for iterrative solver:\n");
    scanf("%lf", &tol);
    printf("Input mits - Maximum iterations for solver:\n");
    scanf("%d", &mits);
#else
    n = 200;
    m = 200;
    alpha = 0.8;
    relax = 1.0;
    tol = 1e-10;
    mits = 10000;
#endif
    printf("-> %d, %d, %g, %g, %g, %d\n", n, m, alpha, relax, tol, mits);

    allocCount = (n+2)*(m+2);
    // Those two calls also zero the boundary elements
    u = (double*)calloc(sizeof(double), allocCount);
    u_old = (double*)calloc(sizeof(double), allocCount);
    if (u == NULL || u_old == NULL)
    {
   fprintf(stderr, "Not enough memory for two %ix%i matrices\n", n+2, m+2);
   exit(1);
    }
    maxIterationCount = mits;
    maxAcceptableError = tol;

    // Solve in [-1, 1] x [-1, 1]
    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    double deltaX = (xRight-xLeft)/(n-1);
    double deltaY = (yUp-yBottom)/(m-1);

    iterationCount = 0;
    error = HUGE_VAL;

    struct timeval start, end;
    gettimeofday(&start, NULL);
    
    /* Iterate as long as it takes to meet the convergence criterion */
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
   //printf("Iteration %i\n", iterationCount);
            error = one_jacobi_iteration(xLeft, yUp,
                 n+2, m+2,
                 u_old, u,
                 deltaX, deltaY,
                 alpha, relax);
   //printf("\tError %g\n", error);
   iterationCount++;
   // Swap the buffers
   tmp = u_old;
   u_old = u;
   u = tmp;
    }
    printf("Residual %g\n",error);
    double absoluteError;

    // u_old holds the solution after the most recent buffers swap
    absoluteError = checkSolution(xLeft, yUp, n+2, m+2,u_old,deltaX, deltaY,alpha);
    gettimeofday(&end, NULL);

    printf("The error of the iterative solution is %g\n", absoluteError);
    float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    printf("seconds: %f\n", delta);
    
    return 0;
}