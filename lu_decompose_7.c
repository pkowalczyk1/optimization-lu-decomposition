#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <immintrin.h>

#define IDX(i, j, n) (((j)+ (i)*(n)))
static double gtod_ref_time_sec = 0.0;
static char * result_file_name = "output_7.m";


int LUPDecompose(double *A, int N) {

    register int i, j, k, size;
    register double division_var;
    register __m256d tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, mult;
    size = N;

    for (i = 0; i < size; i++) {
        division_var = A[IDX(i, i, size)];
        for (j = i + 1; j < size; j++) {
            A[IDX(j,i,size)] /= division_var;
            mult = _mm256_set1_pd(A[IDX(j,i,size)]);

            for (k = i + 1; k < size;) {
                if (k + 15 < size) {
                    tmp0 = _mm256_loadu_pd(A+IDX(i, k, size));
                    tmp1 = _mm256_loadu_pd(A+IDX(i, k+4, size));
                    tmp2 = _mm256_loadu_pd(A+IDX(i, k+8, size));
                    tmp3 = _mm256_loadu_pd(A+IDX(i, k+12, size));

                    tmp4 = _mm256_loadu_pd(A+IDX(j, k, size));
                    tmp5 = _mm256_loadu_pd(A+IDX(j, k+4, size));
                    tmp6 = _mm256_loadu_pd(A+IDX(j, k+8, size));
                    tmp7 = _mm256_loadu_pd(A+IDX(j, k+12, size));

                    tmp0 = _mm256_mul_pd(mult, tmp0);
                    tmp1 = _mm256_mul_pd(mult, tmp1);
                    tmp2 = _mm256_mul_pd(mult, tmp2);
                    tmp3 = _mm256_mul_pd(mult, tmp3);

                    tmp0 = _mm256_sub_pd(tmp4, tmp0);
                    tmp1 = _mm256_sub_pd(tmp5, tmp1);
                    tmp2 = _mm256_sub_pd(tmp6, tmp2);
                    tmp3 = _mm256_sub_pd(tmp7, tmp3);

                    A[IDX(j, k, size)] = tmp0[0];
                    A[IDX(j, k+1, size)] = tmp0[1];
                    A[IDX(j, k+2, size)] = tmp0[2];
                    A[IDX(j, k+3, size)] = tmp0[3];
                    A[IDX(j, k+4, size)] = tmp1[0];
                    A[IDX(j, k+5, size)] = tmp1[1];
                    A[IDX(j, k+6, size)] = tmp1[2];
                    A[IDX(j, k+7, size)] = tmp1[3];
                    A[IDX(j, k+8, size)] = tmp2[0];
                    A[IDX(j, k+9, size)] = tmp2[1];
                    A[IDX(j, k+10, size)] = tmp2[2];
                    A[IDX(j, k+11, size)] = tmp2[3];
                    A[IDX(j, k+12, size)] = tmp3[0];
                    A[IDX(j, k+13, size)] = tmp3[1];
                    A[IDX(j, k+14, size)] = tmp3[2];
                    A[IDX(j, k+15, size)] = tmp3[3];

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
    fprintf(output, "version = 'LU_Decompose_7';");
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