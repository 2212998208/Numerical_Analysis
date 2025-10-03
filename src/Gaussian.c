#include "Gaussian.h"
#include <stdlib.h>
#include <math.h>



/* ------------------ 核心：扁平化 + 行跨度 lda ------------------ */
/* 高斯消元（部分选主元），A 为 n x (n+1) 的增广矩阵（就地修改）
 * 返回 0 成功；1 奇异/病态；2 参数非法 */
GAUSSIAN_Err gauss_pp_core(size_t n, double *A, size_t lda, double *x) {
    if (!A || !x || lda < n+1) return GAUSSIAN_INVALID_INPUT;
    const double EPS = 1e-12;  /* 根据数据量级可调 */

    for (size_t k = 0; k < n; ++k) {
        /* 1) 选主元（列 k 的 k..n-1 中 |a_ik| 最大） */
        size_t piv = k;
        double maxv = fabs(A[AIDX(k, k, lda)]);
        for (size_t i = k + 1; i < n; ++i) {
            double v = fabs(A[AIDX(i, k, lda)]);
            if (v > maxv) { maxv = v; piv = i; }
        }
        if (maxv < EPS) return GAUSSIAN_BAD_MATRIX; /* 奇异或严重病态 */

        /* 2) 行交换（含右端列 j=n） */
        if (piv != k) {
            for (size_t j = k; j <= n; ++j) {
                double tmp = A[AIDX(k, j, lda)];
                A[AIDX(k, j, lda)]   = A[AIDX(piv, j, lda)];
                A[AIDX(piv, j, lda)] = tmp;
            }
        }

        /* 3) 归一化主元行（让 a_kk = 1，便于回代） */
        double akk = A[AIDX(k, k, lda)];
        for (size_t j = k; j <= n; ++j)
            A[AIDX(k, j, lda)] /= akk;

        /* 4) 用第 k 行消掉其下方元素 */
        for (size_t i = k + 1; i < n; ++i) {
            double lik = A[AIDX(i, k, lda)];
            if (fabs(lik) < EPS) continue;
            for (size_t j = k; j <= n; ++j) {
                A[AIDX(i, j, lda)] -= lik * A[AIDX(k, j, lda)];
            }
        }
    }

    /* 5) 回代（对角已归一） */
    for (ptrdiff_t i = (ptrdiff_t)n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (size_t j = (size_t)i + 1; j < n; ++j)
            sum += A[AIDX((size_t)i, j, lda)] * x[j];
        x[i] = A[AIDX((size_t)i, n, lda)] - sum;
    }
    return GAUSSIAN_SUCCESS;
}


// 高斯-乔丹法
/* Gauss–Jordan elimination with partial pivoting.
 * Solves A_left * x = b by operating in-place on the augmented matrix A (n x (n+1)):
 * - Input : A   -> pointer to the first element of augmented matrix [A_left | b]
 *           n   -> number of unknowns (rows = n, cols = n+1)
 *           lda -> row stride (must be >= n+1)
 * - Output: x   -> solution vector of length n
 * Returns: 0 on success; 1 if (near-)singular; -1 on invalid args.
 *
 * This routine transforms A to [I | x] by eliminating both above and below each pivot.
 * It modifies A in-place.
 */
GAUSSIAN_Err gauss_jordan_solve(size_t n, double *A, size_t lda, double *x)
{
    if (!A || !x || lda < n + 1) return GAUSSIAN_INVALID_INPUT;
    const double EPS = 1e-12;

    for (size_t k = 0; k < n; ++k) {
        /* 1) Partial pivoting: find pivot row with max |A[i,k]|, i>=k */
        size_t piv = k;
        double maxv = A[k*lda + k]; if (maxv < 0) maxv = -maxv;
        for (size_t i = k + 1; i < n; ++i) {
            double v = A[i*lda + k]; if (v < 0) v = -v;
            if (v > maxv) { maxv = v; piv = i; }
        }
        if (maxv < EPS) return GAUSSIAN_BAD_MATRIX; /* singular or ill-conditioned */

        /* 2) Swap current row k with pivot row piv (all columns 0..n) */
        if (piv != k) {
            for (size_t j = 0; j <= n; ++j) {
                double tmp = A[k*lda + j];
                A[k*lda + j]   = A[piv*lda + j];
                A[piv*lda + j] = tmp;
            }
        }

        /* 3) Normalize pivot row so that A[k,k] = 1 */
        double akk = A[k*lda + k];
        for (size_t j = 0; j <= n; ++j) A[k*lda + j] /= akk;

        /* 4) Eliminate column k in all other rows (above and below) */
        for (size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            double factor = A[i*lda + k];
            if (factor == 0.0) continue;
            for (size_t j = 0; j <= n; ++j) {
                A[i*lda + j] -= factor * A[k*lda + j];
            }
        }
    }

    /* 5) Read solution: left block is I, rightmost column is x */
    for (size_t i = 0; i < n; ++i) x[i] = A[i*lda + n];
    return GAUSSIAN_SUCCESS;
}

// LU分解法
/* ============================================================
 * LU decomposition with partial pivoting (Doolittle form)
 *   PA = LU
 * - A: in-place stores L (strictly lower) and U (upper incl. diag)
 * - L has unit diagonal (not stored in A; extract when needed)
 * - piv: permutation vector, size n; row i swapped with piv[i] in step i
 * - Returns: 0 ok; 1 near-singular; -1 invalid args
 * ============================================================ */
GAUSSIAN_Err lu_decompose_pp(size_t n, double *A, size_t lda, size_t *piv) {
    if (!A || !piv || lda < n) return GAUSSIAN_INVALID_INPUT;
    const double EPS = 1e-12;

    for (size_t i = 0; i < n; ++i) piv[i] = i;

    for (size_t k = 0; k < n; ++k) {
        /* --- choose pivot row r with max |A[r,k]|, r >= k --- */
        size_t r = k;
        double maxv = fabs(A[IDX(k,k,lda)]);
        for (size_t i = k + 1; i < n; ++i) {
            double v = fabs(A[IDX(i,k,lda)]);
            if (v > maxv) { maxv = v; r = i; }
        }
        if (maxv < EPS) return GAUSSIAN_BAD_MATRIX; /* singular or ill-conditioned */

        /* --- swap rows k <-> r (all columns 0..n-1) --- */
        if (r != k) {
            for (size_t j = 0; j < n; ++j) {
                double tmp = A[IDX(k,j,lda)];
                A[IDX(k,j,lda)] = A[IDX(r,j,lda)];
                A[IDX(r,j,lda)] = tmp;
            }
            size_t tp = piv[k]; piv[k] = piv[r]; piv[r] = tp;
        }

        /* --- factorization step: eliminate below pivot --- */
        double akk = A[IDX(k,k,lda)];
        for (size_t i = k + 1; i < n; ++i) {
            A[IDX(i,k,lda)] /= akk;                /* L(i,k) */
            double lik = A[IDX(i,k,lda)];
            for (size_t j = k + 1; j < n; ++j) {
                A[IDX(i,j,lda)] -= lik * A[IDX(k,j,lda)];
            }
        }
    }
    return GAUSSIAN_SUCCESS;
}

/* 显式提取 L 与 U（便于显示/验证）：L 为单位下三角，U 为上三角 */
void lu_extract(size_t n, const double *LU, size_t lda,
                double *L, size_t ldl,
                double *U, size_t ldu) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i > j) {      /* 严格下三角 -> L */
                L[IDX(i,j,ldl)] = LU[IDX(i,j,lda)];
                U[IDX(i,j,ldu)] = 0.0;
            } else if (i == j) { /* 对角：L=1, U=diag */
                L[IDX(i,j,ldl)] = 1.0;
                U[IDX(i,j,ldu)] = LU[IDX(i,j,lda)];
            } else {          /* 上三角 -> U */
                L[IDX(i,j,ldl)] = 0.0;
                U[IDX(i,j,ldu)] = LU[IDX(i,j,lda)];
            }
        }
    }
}
