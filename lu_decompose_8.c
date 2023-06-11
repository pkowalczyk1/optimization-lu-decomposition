#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define IDX(i, j, n) (((j)+ (i)*(n)))
static double gtod_ref_time_sec = 0.0;
static char * result_file_name = "output_O3.m";



int LUPDecompose(double *A, int N) {

    int i, j, k; 

    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            A[IDX(j,i,N)] /= A[IDX(i,i,N)];

            for (k = i + 1; k < N; k++)
                A[IDX(j,k,N)] -= A[IDX(j,i,N)] * A[IDX(i,k,N)];
        }
    }

    return 0;  //decomposition done 
}

double dclock()
{
    double the_time, norm_sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    if (gtod_ref_time_sec == 0.0)
        gtod_ref_time_sec = (double)tv.tv_sec;
    norm_sec = (double)tv.tv_sec - gtod_ref_time_sec;
    the_time = norm_sec + tv.tv_usec * 1.0e-6;
    return the_time;
}

double calculate_gflops(int size) {
    double flops = 0;
    for (int i=0; i<size; i++) {
        flops += (size - i - 1) + 2 * (size - i - 1) * (size - i - 1);
    }
    return flops * 1.0e-09;
}

int main(){
    int i,j, iret, reps = 5;

    double dtime, dtime_best;

    double* matrix;

    FILE * output = fopen(result_file_name, "w");
    fprintf(output, "version = 'LU_Decompose_O3';");
    fprintf(output, "LU_Decompose = [\n");

    for (int size = 100; size <= 1000; size += 40)
    {
        matrix = malloc(size*size*sizeof(double));

        srand(1);

        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                matrix[IDX(i, j, size)] = rand();
            }
        }

        for (int rep=0; rep<reps; rep++) {
            dtime = dclock();
            iret = LUPDecompose(matrix, size);
            dtime = dclock() - dtime;
            if (rep == 0) {
                dtime_best = dtime;
            } else {
                dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
            }
        }
        printf("Time: %le \n", dtime_best);
        fprintf(output, "%d %le 0.1 \n", size, calculate_gflops(size) / dtime_best);

        double check = 0.0;
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                check = check + matrix[IDX(i,j,size)];
            }
        }
        printf("Check:%le \n", check);
        fflush(stdout);
        free(matrix);
    }
    fprintf(output, "];\n");
    fclose(output);
}