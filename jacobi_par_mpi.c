#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>


#define TRUE 1
#define FALSE 0 
#define SRC(YY,XX) src[(YY)*inp.col_offset+(XX)]
#define DST(YY,XX) dst[(YY)*inp.col_offset+(XX)]
#define U(YY,XX) u[(YY)*inp.col_offset+(XX)]

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

typedef struct in{
        double alpha;
        double relax;
        double tol;
        int mits;
        int line_offset;
        int col_offset;
        double deltaX;
        double deltaY;
}input_data;

double NonCrit(double xStart, double yStart, double *src, double *dst, input_data inp)
{

    int i, j;
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;

    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;

    for (i = 2; i < inp.line_offset - 2 ; ++i)
    {
        fY = yStart - (i-1)*inp.deltaY;
        for (j = 2; j < inp.col_offset -2; ++j)
        {
            fX = xStart + (j-1)*inp.deltaX;
            f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);
            updateVal = ((SRC(i-1,j) + SRC(i+1,j))*cx + (SRC(i,j-1)+SRC(i,j+1))*cy + SRC(i,j)*cc - f) /cc;
            DST(i,j) = SRC(i,j) - inp.relax*updateVal;
            error += updateVal*updateVal;
        }
    }
    return sqrt(error)/((inp.line_offset-4)*(inp.col_offset-4));
}

double CritLeft(double xStart, double yStart,input_data inp,double *src, double *dst)
{

    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;

    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;

    int i;
    fX = xStart;
    for(i = 1; i < inp.line_offset - 1 ;++i)
    {
          fY = yStart - (i - 1)*inp.deltaY;
          f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);

          updateVal = ((SRC(1,i-1) + SRC(1,i+1))*cx + (SRC(0,i)+SRC(2,i))*cy + SRC(1,i)*cc - f) /cc;

          DST(i,0) = SRC(0,1) - inp.relax*updateVal;
          error += updateVal*updateVal;
    }
    return sqrt(error)/((inp.line_offset-2));
}

double CritRight(double xStart, double yStart,input_data inp,double *src, double *dst)
{
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;

    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;

    int i;
    fX = xStart + (inp.col_offset-2)*inp.deltaX;

    for(i = 1; i < inp.line_offset -1 ;++i)
    {
          fY = yStart - (i-1)*inp.deltaY;
          f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);

          updateVal = ((SRC(inp.col_offset - 2,i-1) + SRC(inp.col_offset - 2,i+1))*cx 
                        + (SRC(inp.col_offset-3,i)+SRC(inp.col_offset-1,i))*cy + SRC(inp.col_offset - 2,i)*cc - f) /cc;

          DST(inp.col_offset - 2,i) = SRC(inp.col_offset - 2,i) - inp.relax*updateVal;
          error += updateVal*updateVal;
    }
    return sqrt(error)/((inp.line_offset-2));
}

double CritUp(double xStart, double yStart,input_data inp,double *src, double *dst)
{
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;

    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;

    int j;
    
    fY = yStart;
    for(j = 2; j < inp.col_offset -2 ;++j)
    {
          fX = xStart + (j-1)*inp.deltaX;
          f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);

          updateVal = ((SRC(j,0) + SRC(j,2))*cx + (SRC(j-1,1)+SRC(j+1,1))*cy + SRC(j,1)*cc - f) /cc;

          DST(j,1) = SRC(j,1) - inp.relax*updateVal;
          error += updateVal*updateVal;
    }
    return sqrt(error)/((inp.col_offset-4));
}

double CritDown(double xStart, double yStart,input_data inp,double *src, double *dst)
{
    double fX, fY;
    double error = 0.0;
    double updateVal;
    double f;
    double cx = 1.0/(inp.deltaX*inp.deltaX);
    double cy = 1.0/(inp.deltaY*inp.deltaY);
    double cc = -2.0*cx-2.0*cy-inp.alpha;

    int j;
    
    fY = yStart - (inp.line_offset-2)*inp.deltaY;

    for(j = 2; j < inp.col_offset - 2 ;++j)
    {
          fX = xStart + (j-1)*inp.deltaX;
          f = -inp.alpha * (1.0-fX*fX) * (1.0-fY*fY) - 2.0*(1.0-fX*fX) - 2.0*(1.0-fY*fY);

          updateVal = ((SRC(j,inp.line_offset-3) + SRC(j,inp.line_offset-1))*cx + (SRC(j-1,inp.line_offset-2)
                        + SRC(j+1,inp.line_offset-2))*cy + SRC(j,inp.line_offset-2)*cc - f) /cc;

          DST(j,inp.line_offset-2) = SRC(j,inp.line_offset-2) - inp.relax*updateVal;
          error += updateVal*updateVal;
    }
    return sqrt(error)/((inp.col_offset-4));
}

 
double checkSolution(double xStart, double yStart,double *u, input_data inp)
{
    
    int i, j;
    double fX, fY;
    double localError, error = 0.0;

    for (i = 1; i < (inp.line_offset-1); ++i)
    {
        fY = yStart + (i-1)*inp.deltaY;
        for (j = 1; j < (inp.col_offset-1); ++j)
        {
            fX = xStart + (j-1)*inp.deltaX;
            localError = U(i,j) - (1.0-fX*fX)*(1.0-fY*fY);
            error += localError*localError;
        }
    }
    return sqrt(error)/((inp.line_offset-2)*(inp.col_offset-2));
}

void send_data(int recv, int direction,double* matrix, input_data inp, MPI_Datatype type, MPI_Request req[8], MPI_Comm my_grid)
{   

    if(recv == MPI_PROC_NULL || recv<0)
        return;

    int i;
    double* send = NULL;
    if(direction == UP)
    {
        send = (double*)calloc(sizeof(double), inp.col_offset);   
        for (i = 0; i < inp.col_offset; ++i)
        {
            send[i] = matrix[i];
        }
    }

    if(direction == DOWN)
    {
        send = (double*)calloc(sizeof(double), inp.col_offset);
        for (i = 0; i < inp.col_offset; ++i)
        {
            send[i] = matrix[i+(inp.col_offset * (inp.line_offset-1))];
        }
    }

    if(direction == LEFT)
    {
        send = (double*)calloc(sizeof(double), inp.line_offset);
        for (i = 0; i < inp.line_offset; ++i)
        {
            send[i] = matrix[i*inp.col_offset];
        }
    }

    if(direction == RIGHT)
    {
        send = (double*)calloc(sizeof(double), inp.line_offset);
        for (i = 0; i < inp.line_offset; ++i)
        {
            send[i] = matrix[(i*inp.col_offset) + (inp.line_offset-1)];
        }
    }        

    MPI_Isend(send, 1, type, recv, 0, my_grid,req);
}

void receive(int sender,int direction,double *matrix, input_data inp, MPI_Datatype type, MPI_Request req[8],MPI_Status stat[8],MPI_Comm my_grid)
{
    double* recv = NULL;
    int i;
    if(direction == UP)
    {
        recv = (double*)calloc(sizeof(double), inp.col_offset);
        MPI_Irecv(recv,1, type, sender, 0, my_grid, req);
        MPI_Wait(req,stat);
        for (i = 1; i < inp.col_offset - 1; ++i)
        {
            matrix[i] = recv[i];
        }
    }
    if(direction == DOWN)
    {
        recv = (double*)calloc(sizeof(double), inp.col_offset);      
        MPI_Irecv(recv,1, type, sender, 0, my_grid, req);
        MPI_Wait(req,stat);
        for (i = 1; i < inp.col_offset - 1; ++i)
        {
            matrix[(inp.col_offset*(inp.line_offset - 1)) + i] = recv[i];
        }        
    }
    if(direction == RIGHT)
    {
        recv = (double*)calloc(sizeof(double), inp.line_offset);
        MPI_Irecv(recv,1, type, sender, 0, my_grid, req);
        MPI_Wait(req,stat);
        for (i = 0; i < inp.line_offset; ++i)
        {
            matrix[(inp.line_offset - 1) + inp.col_offset*i] = recv[i];
        }
    }
    
    if(direction == LEFT)
    {
        recv = (double*)calloc(sizeof(double), inp.line_offset);
        MPI_Irecv(recv,1, type, sender, 0, my_grid, req);
        MPI_Wait(req,stat);
        for (i = 0; i < inp.line_offset; ++i)
        {
            matrix[inp.col_offset*i] = recv[i];
        }        
    }       
}

void update(double *dst, double *src, input_data inp)
{
    int i,j;
    for(i=0;i<inp.line_offset;i++)
    {
        for(j=0;j<inp.col_offset;j++)
        {
            DST(i,j) = SRC(i,j);    
        }
    }
}

int main(int argc, char **argv)
{
    int lines, cols;
    int rank,comm_size;
    MPI_Comm my_grid;
    int dim[2],period[2],reorder;
    int up,down,right,left;
    //int err=0;

    MPI_Init(&argc, &argv);
    double MPI_Wtime(void);
    double start, finish;
    start=MPI_Wtime();
    int tam = 0;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);

    tam = sqrt(comm_size);
    dim[0] = tam; 
    dim[1] = tam;
    dim[0] = 1*(comm_size-(dim[0]-dim[1]))/dim[1];

    period[0]=FALSE; 
    period[1]=FALSE;
    reorder=FALSE;

    MPI_Cart_create(MPI_COMM_WORLD,2,dim,period,reorder,&my_grid);


    MPI_Cart_shift(my_grid,0,1,&left,&right);
    MPI_Cart_shift(my_grid,1,1,&up,&down);

    printf("P:%d My neighbors are r: %d d:%d 1:%d u:%d\n",rank,right,down,left,up);

    input_data inp;

    int count = 8;
    int array_of_blocklengths[8] = {1, 1, 1, 1, 1, 1, 1, 1};

    MPI_Aint i_addr, alpha_addr, relax_addr, tol_addr, mits_addr, lines_addr, cols_addr, deltaX_addr, deltaY_addr;
    MPI_Aint array_of_displacements[8];

    MPI_Get_address(&inp,&i_addr);
    MPI_Get_address(&(inp.alpha),&alpha_addr);
    array_of_displacements[0] = alpha_addr - i_addr;
    MPI_Get_address(&(inp.relax),&relax_addr);
    array_of_displacements[1] = relax_addr - i_addr;
    MPI_Get_address(&(inp.tol),&tol_addr);
    array_of_displacements[2] = tol_addr - i_addr;
    MPI_Get_address(&(inp.mits),&mits_addr);
    array_of_displacements[3] = mits_addr - i_addr;
    MPI_Get_address(&(inp.line_offset),&lines_addr);
    array_of_displacements[4] = lines_addr - i_addr;
    MPI_Get_address(&(inp.col_offset),&cols_addr);
    array_of_displacements[5] = cols_addr - i_addr;
    MPI_Get_address(&(inp.deltaX),&deltaX_addr);
    array_of_displacements[6] = deltaX_addr - i_addr;
    MPI_Get_address(&(inp.deltaY),&deltaY_addr);
    array_of_displacements[7] = deltaY_addr - i_addr;

    MPI_Datatype array_of_types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Datatype input_type;
    MPI_Type_create_struct(count,array_of_blocklengths,array_of_displacements,array_of_types,&input_type);
    MPI_Type_commit(&input_type);

    if(rank ==0)
    {
        int flag = 1;
        if(flag ==1)
        {

    	   FILE* f = fopen(argv[1],"r");
            if ( f == 0 )
            {
                printf( "Could not open file\n" );
            }
            else
            {
                fscanf(f, "%d %d %lf %lf %lf %d", &lines, &cols, &(inp.alpha),&(inp.relax),&(inp.tol),&(inp.mits));
                fclose(f);
                printf("-> %d, %d, %lf, %lf, %.12lf, %d\n", lines, cols, inp.alpha, inp.relax, inp.tol, inp.mits);            
            }
        }   
        else
        {
            lines = 64;
            cols= 64;
            inp.alpha = 0.8;
            inp.relax = 1.0;
            inp.tol = 1e-7;
            inp.mits = 100;
        }
	    inp.line_offset = (lines / dim[0]) + 2;
        inp.col_offset = (cols / dim[1]) + 2 ;

        double xLeft = -1.0, xRight = 1.0;
        double yBottom = -1.0, yUp = 1.0;

        inp.deltaX = (xRight-xLeft)/(cols);
        inp.deltaY = (yUp-yBottom)/(lines);
        
    }
 

    MPI_Bcast(&(inp),1,input_type,0,my_grid);

    //if(rank == 0)
        //printf("-> %d, %d, %g, %g, %g, %d, %g, %g\n", inp.line_offset, inp.col_offset, inp.alpha, inp.relax, inp.tol, inp.mits, inp.deltaX, inp.deltaY);

    MPI_Datatype row;
    MPI_Type_contiguous(inp.col_offset, MPI_DOUBLE, &row);
    MPI_Type_commit(&row);

    MPI_Datatype col;
    MPI_Type_contiguous(inp.line_offset, MPI_DOUBLE, &col);
    MPI_Type_commit(&col);

    double* mlocal_old = (double*)calloc(sizeof(double), inp.line_offset*inp.col_offset);
    double* mlocal_new = (double*)calloc(sizeof(double), inp.line_offset*inp.col_offset);

    if (mlocal_old == NULL || mlocal_new == NULL)
    {
          fprintf(stderr, "Not enough memory for two %ix%i matrices\n", inp.line_offset, inp.col_offset);
          exit(1);
    }

    double errorNC=0,errorCD=0,errorCU=0,errorCL=0,errorCR=0;
    int iterationCount, maxIterationCount;

    maxIterationCount = inp.mits;

    double maxAcceptableError = inp.tol;

    iterationCount = 0;

    double error = 1;

    MPI_Request req[8];
    MPI_Status stat[8];

    int coords[2];
    MPI_Cart_coords(my_grid,rank,2,coords);
    //printf("My rand is %d and my coords are x = %d and y = %d\n",rank, coords[0],coords[1]);


    // Iterate as long as it takes to meet the convergence criterion
    while (iterationCount < maxIterationCount && error > maxAcceptableError)
    {
        //printf("Iteration %i\n", iterationCount);

        if(up >= 0 && up != MPI_PROC_NULL)
        {
            send_data(up, UP, mlocal_old, inp, row, req,my_grid);
            //printf("Data sent from %d to %d\n", rank,up);
        }
        if (down >= 0 && down != MPI_PROC_NULL)
        {
            send_data(down, DOWN, mlocal_old, inp, row, req,my_grid);
            //printf("Data sent from %d to %d\n", rank,down);

        }
        if (left >= 0 && left != MPI_PROC_NULL)
        {
            send_data(left, LEFT, mlocal_old, inp, col, req,my_grid);
            //printf("Data sent from %d to %d\n", rank,left);

        }

        if (right >= 0 && right != MPI_PROC_NULL)
        {
            send_data(right, RIGHT, mlocal_old, inp, col, req,my_grid);
            //printf("Data sent from %d to %d\n", rank,right);

        }

        errorNC = NonCrit(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY), mlocal_old , mlocal_new,inp);
        //printf("%g\n", errorNC);
        if(down >= 0 && down != MPI_PROC_NULL)
        {

            receive(down,DOWN,mlocal_old,inp,row,req,stat,my_grid);
            //printf("Passou rank: %d\n", rank);
            errorCD=CritDown(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY),inp, mlocal_old , mlocal_new);
            //printf("%g\n", errorCD);
        }

        if(right >= 0 && right != MPI_PROC_NULL)
        {
            receive(right,RIGHT,mlocal_old,inp,col,req,stat,my_grid);
            errorCR=CritRight(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY),inp,mlocal_old, mlocal_new);
            //printf("%g\n", errorCR);
        }

        if(up >= 0 && up != MPI_PROC_NULL)
        {
            receive(up,UP,mlocal_old,inp,row,req,stat,my_grid);
            errorCU=CritUp(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY),inp,mlocal_old,mlocal_new);
            //printf("%g\n", errorCU);
        }

        if(left >= 0 && left != MPI_PROC_NULL)
        {
            receive(left,LEFT,mlocal_old,inp,col,req,stat,my_grid);
            errorCL=CritLeft(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY),inp,mlocal_old,mlocal_new);
            //printf("%g\n", errorCL);
        }

        errorNC += errorCU + errorCR + errorCL + errorCD;
        iterationCount++;
        //printf("Passou rank: %d\n", rank);
        update(mlocal_old,mlocal_new,inp);

        MPI_Allreduce(&errorNC,&error,1,MPI_DOUBLE,MPI_SUM,my_grid);
    }

    double check_error = 0;
    check_error = checkSolution(-1 + (coords[0]*inp.col_offset*inp.deltaX), 1 - (coords[1]*inp.line_offset*inp.deltaY),mlocal_new,inp);
    double final_error = 0;
    MPI_Reduce(&check_error,&final_error,1,MPI_DOUBLE,MPI_SUM,0,my_grid);
    finish=MPI_Wtime();

    if(rank == 0)
    {
        printf("tempo decorrido: %f\n", finish - start);
        printf("Iteration %i\n", iterationCount);
        printf("Residual %g, rank: %d\n",error,rank);
        printf("The error of the iterative solution is %g, rank %d\n", final_error,rank);
    }
    /*
    free(mlocal_old);
    free(mlocal_new);
    */
    MPI_Finalize();

    return 0;
}