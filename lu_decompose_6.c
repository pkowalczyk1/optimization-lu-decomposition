#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define IDX(i, j, n) (((j)+ (i)*(n)))
static double gtod_ref_time_sec = 0.0;
static char * result_file_name = "output_6.m";


int LUPDecompose(double *A, int N) {

    register int i, j, k, size;
    register double division_var, temp;
    size = N;

    for (i = 0; i < size; i++) {
        division_var = A[IDX(i, i, size)];
        for (j = i + 1; j < size; j++) {
            A[IDX(j,i,size)] /= division_var;
            temp = A[IDX(j,i,size)];

            for (k = i + 1; k < size;) {
                if (k + 15 < size) {
                    temp = A[IDX(j, i, size)];
                    A[IDX(j, k, size)] -= temp * A[IDX(i, k, size)];
                    A[IDX(j, k+1, size)] -= temp * A[IDX(i, k+1, size)];
                    A[IDX(j, k+2, size)] -= temp * A[IDX(i, k+2, size)];
                    A[IDX(j, k+3, size)] -= temp * A[IDX(i, k+3, size)];
                    A[IDX(j, k+4, size)] -= temp * A[IDX(i, k+4, size)];
                    A[IDX(j, k+5, size)] -= temp * A[IDX(i, k+5, size)];
                    A[IDX(j, k+6, size)] -= temp * A[IDX(i, k+6, size)];
                    A[IDX(j, k+7, size)] -= temp * A[IDX(i, k+7, size)];
                    A[IDX(j, k+8, size)] -= temp * A[IDX(i, k+8, size)];
                    A[IDX(j, k+9, size)] -= temp * A[IDX(i, k+9, size)];
                    A[IDX(j, k+10, size)] -= temp * A[IDX(i, k+10, size)];
                    A[IDX(j, k+11, size)] -= temp * A[IDX(i, k+11, size)];
                    A[IDX(j, k+12, size)] -= temp * A[IDX(i, k+12, size)];
                    A[IDX(j, k+13, size)] -= temp * A[IDX(i, k+13, size)];
                    A[IDX(j, k+14, size)] -= temp * A[IDX(i, k+14, size)];
                    A[IDX(j, k+15, size)] -= temp * A[IDX(i, k+15, size)];
                    k += 16;
                } else {
                    A[IDX(j, k, size)] -= A[IDX(j, i, size)] * A[IDX(i, k, size)];
                    k++;
                }
            }
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
    fprintf(output, "version = 'LU_Decompose_6';");
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