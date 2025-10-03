#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Gaussian.h"


/* (B) 拷贝式：不修改原 A，把 2D 拷贝到连续缓冲区，再求解 */
GAUSSIAN_Err gauss_pp_from_2d_copy(size_t n, const double A[][n+1], double x[]) {
    size_t lda = n + 1;
    double *buf = (double*)malloc(n * lda * sizeof(double));
    if (!buf) return -1;

    /* 二维→扁平化：逐行 memcpy（紧致：lda = n+1） */
    for (size_t i = 0; i < n; ++i)
        memcpy(&buf[AIDX(i,0,lda)], &A[i][0], (n+1) * sizeof(double));

    GAUSSIAN_Err ret = gauss_pp_core(n, buf, lda, x);
    free(buf);
    return ret;
}

GAUSSIAN_Err gauss_jordan_from_2d_copy(size_t n, const double A[][n+1], double x[]) {
    size_t lda = n + 1;
    double *buf = (double*)malloc(n * lda * sizeof(double));
    if (!buf) return -1;

    /* 二维→扁平化：逐行 memcpy（紧致：lda = n+1） */
    for (size_t i = 0; i < n; ++i)
        memcpy(&buf[AIDX(i,0,lda)], &A[i][0], (n+1) * sizeof(double));

    GAUSSIAN_Err ret = gauss_jordan_solve(n, buf, lda, x);
    free(buf);
    return ret;
}

/* ------------------ 演示：用 5x6 增广矩阵测试 ------------------ */
static void print_vec(const char *name, const double *x, size_t n) {
    printf("%s = [", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%s%.10g", (i? ", ":""), x[i]);
    }
    puts("]");
}

/* 简单打印工具 */
void print_mat(const char *name, const double *A, size_t n, size_t lda) {
    printf("%s =\n", name);
    for (size_t i = 0; i < n; ++i) {
        printf("  ");
        for (size_t j = 0; j < n; ++j) {
            printf("%10.6f%s", A[IDX(i,j,lda)], (j+1==n? "":" "));
        }
        puts("");
    }
}

int main(void) {
    /* 与我们之前的 5x6 案例一致（5 个未知数 + 常数列） */
    double A[5][6] = {
        { 2,  3, -1,  1,  2,  4},
        { 1, -1,  2, -2,  1, -1},
        { 3,  2,  3,  1,  4, 10},
        { 2,  1,  1,  1, -1,  5},
        { 1,  4, -2,  2,  3,  7}
    };
    double x[5];
    double x1[5];

    double B[3][4] = {
        {2, 1, -1, 8},
        {-3, -1, 2, -11},
        {-2, 1, 2, -3}
    };
    double y[3];
    double y1[3];

    /* 方式二：拷贝式（不修改原 A） */
    {
        GAUSSIAN_Err ret = gauss_pp_from_2d_copy(5, A, x);
        GAUSSIAN_Err ret1 = gauss_jordan_from_2d_copy(5,A,x1);
        if (ret != GAUSSIAN_SUCCESS || ret1 != GAUSSIAN_SUCCESS) {
            printf("copy solver failed: %d %d\n", ret, ret1);
            return 1;
        }
        print_vec("x (copy)   ", x, 5);
        print_vec("x1 (copy)   ", x1, 5);
    }

    /* 期望解（有理数）：[-46/15, 86/15, 13/3, -1/5, -19/15]
       数值约   ：[-3.066666..., 5.733333..., 4.333333..., -0.2, -1.266666...] */

    {
        GAUSSIAN_Err ret = gauss_pp_from_2d_copy(3, B, y);
        GAUSSIAN_Err ret1 = gauss_jordan_from_2d_copy(3,B,y1);
        if (ret != GAUSSIAN_SUCCESS || ret1 != GAUSSIAN_SUCCESS) {
            printf("copy solver failed: %d %d\n", ret, ret1);
            return 1;
        }
        print_vec("y (copy)   ", y, 3);
        print_vec("y1 (copy)   ", y1, 3);
    }
    /* 期望解（有理数）：[2, 3, -1]
       数值约   ：[2, 3, -1] */

    {
        const size_t n = 3; const size_t lda = n;
        double C[3][3] = {
            {2, 1, -1},
            {4, 5, -5},
            {-2, -5, 7}
        };

        size_t piv[3];
        GAUSSIAN_Err ret2 = lu_decompose_pp(n, &C[0][0], lda, piv);
        if (ret2 != GAUSSIAN_SUCCESS) {
            fprintf(stderr, "LU failed: %d\n", ret2);
            return 1;
        }

        double L[3][3], U[3][3];
        lu_extract(n, &C[0][0], lda, &L[0][0], n, &U[0][0], n);

        print_mat("L", &L[0][0], n, n);
        print_mat("U", &U[0][0], n, n);

        /* 如需验证：构造置换矩阵 P（由 piv 给出），检查 P*A0 ≈ L*U */
        /* 此处略去验证代码，专注于分解与显示 */
    }

    return 0;
}