#ifndef NUMERICAL_ANALYSIS_GAUSSIAN_H
#define NUMERICAL_ANALYSIS_GAUSSIAN_H
#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>
#define AIDX(i,j,lda) ((i)*(lda) + (j))
#define IDX(i,j,lda) ((i)*(lda) + (j))

typedef enum Gaussian{
    GAUSSIAN_SUCCESS = 0,
    GAUSSIAN_BAD_MATRIX = 1,
    GAUSSIAN_INVALID_INPUT = 2
}GAUSSIAN_Err;


// 非API设计
GAUSSIAN_Err gauss_pp_core(size_t n, double *A, size_t lda, double *x);
GAUSSIAN_Err gauss_jordan_solve(size_t n, double *A, size_t lda, double *x);
GAUSSIAN_Err lu_decompose_pp(size_t n, double *A, size_t lda, size_t *piv);
void lu_extract(size_t n, const double *LU, size_t lda,
                double *L, size_t ldl,
                double *U, size_t ldu);



#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_GAUSSIAN_H