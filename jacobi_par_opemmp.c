#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

typedef struct in{
        double alpha;
        double relax;
        double tol;
        int mits;
        int lines;
        int cols;
        double deltaX;
        double deltaY;
        int thread_count;
}input_data;

#define SRC(YY,XX) src[(YY)*inp.cols+(XX)]
#define DST(YY,XX) dst[(YY)*inp.cols+(XX)]
#define U(YY,XX) u[(YY)*inp.cols+(XX)]


double one_jacobi_iteration(double xStart, double yStart, input_data inp, double *src, double *dst)
{
    int i, j;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;

    // Coefficients
    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;


    # pragma omp parallel num_threads(inp.thread_count) default(none) \
                    shared(src,dst,yStart,xStart,inp,cx,cy,cc) private(i,j,fX,fY,f,updateVal) \
                    reduction(+:error)
    # pragma omp for collapse(2)
    for (i = 1; i < (inp.lines-1); i++)
    {
        for (j = 1; j < (inp.cols-1); j++)
        {
            fY = yStart + (i-1)*inp.deltaY;
            fX = xStart + (j-1)*inp.deltaX;
            f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) -2.0 * (1.0-fX*fX) -2.0 * (1.0-fY*fY);
            updateVal = ((SRC(i-1,j) + SRC(i+1,j)) *cx + (SRC(i,j-1) + SRC(i,j+1)) * cy + SRC(i,j)*cc - f) /cc;
            DST(i,j) = SRC(i,j) - inp.relax*updateVal;
            error += updateVal*updateVal;
        }
    }

    return sqrt(error)/((inp.cols-2)*(inp.lines-2));
}

double checkSolution(double xStart, double yStart, input_data inp, double *u)
{
    int j, i;
    double fX, fY;
    double localError, error = 0.0;

    # pragma omp parallel num_threads(inp.thread_count) default(none) \
                shared(u,inp,yStart,xStart) private(i,j,fX,fY,localError) \
                reduction(+:error)
    # pragma omp for collapse(2)
    for (i = 1; i < (inp.cols-1); ++i)
    {   
        for (j = 1; j < (inp.lines-1); ++j)
        {
            fY = yStart + (i-1)*inp.deltaY;
            fX = xStart + (j-1)*inp.deltaX;
            localError = U(i,j) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((inp.cols-2)*(inp.lines-2));
}

void update(double *dst, double *src, input_data inp)
{
    int i,j;
    # pragma omp parallel for num_threads(inp.thread_count) shared(src,dst) private(i,j) collapse(2)
    for(i=0;i<inp.lines;i++)
    {
        for(j=0;j<inp.cols;j++)
        {
            DST(i,j) = SRC(i,j);    
        }
    }
}


int main(int argc, char *argv[])
{
    input_data inp;

    FILE* f = fopen(argv[1],"r");
    if ( f == 0 )
    {
        printf( "Could not open file\n" );
        return 0;
    }
    else
    {
        fscanf(f, "%d %d %lf %lf %lf %d %d", &(inp.lines), &(inp.cols), &(inp.alpha),&(inp.relax),&(inp.tol),&(inp.mits),&(inp.thread_count));
        fclose(f);
        printf("-> %d, %d, %lf, %lf, %.12lf, %d, %d\n", inp.lines, inp.cols, inp.alpha, inp.relax, inp.tol, inp.mits, inp.thread_count);            
    }
    inp.lines+=2;
    inp.cols+=2;

    double* src = (double*)calloc(sizeof(double), inp.lines*inp.cols);
    double* dst = (double*)calloc(sizeof(double), inp.lines*inp.cols);

    double xLeft = -1.0, xRight = 1.0;
    double yBottom = -1.0, yUp = 1.0;

    inp.deltaX = (xRight-xLeft)/(inp.cols);
    inp.deltaY = (yUp-yBottom)/(inp.lines);

    int iterationCount, maxIterationCount;

    maxIterationCount = inp.mits;

    double maxAcceptableError = inp.tol;

    iterationCount = 0;

    double error = 1;


    double start = omp_get_wtime();
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i\n", iterationCount);
        error = one_jacobi_iteration(-1,1,inp,src,dst);
        iterationCount++;
        update(src,dst,inp);
    }
    printf("Residual %g\n",error);

    double absoluteError = checkSolution(-1,1,inp,src);
    double end = omp_get_wtime();
    printf("Iteration %i\n", iterationCount);
    printf("The error of the iterative solution is %g\n", absoluteError);
    printf("Time: \t %f \n", end-start);

    free(src);
    free(dst); 

  return 0;
}