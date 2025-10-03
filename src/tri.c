#include "tri.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef TRI_DEFAULT_TOL
#define TRI_DEFAULT_TOL 1e-12
#endif

struct tri_solver {
    size_t n;
    double *a, *b, *c, *d; /* 内部拷贝的系数与右端 */
    int has_coeffs;
    int has_rhs;
    double tol;
};

/* 将浮点判断为“近零” */
static int tri_is_near_zero(double v, double tol) {
    /* 使用相对-绝对混合准则：|v| <= tol * (1 + |v|) 近似等价于 |v| <= tol */
    return fabs(v) <= tol;
}

const char* tri_strerror(int err) {
    switch (err) {
        case TRI_SUCCESS:        return "success";
        case TRI_ERR_NULL_PTR:   return "null pointer";
        case TRI_ERR_ALLOC:      return "memory allocation failure";
        case TRI_ERR_DIM:        return "invalid dimension";
        case TRI_ERR_BAD_STATE:  return "bad state";
        case TRI_ERR_ZERO_PIVOT: return "zero (or near-zero) pivot encountered";
        case TRI_ERR_SINGULAR:   return "matrix singular or ill-conditioned";
        case TRI_ERR_INTERNAL:   return "internal error";
        default:                 return "unknown error";
    }
}

/* Thomas 前向消元 + 回代（不修改输入数组），需要临时工作区 */
static int tri_thomas_core(size_t n,
                           const double* a,
                           const double* b,
                           const double* c,
                           const double* d,
                           double* x_out,
                           double tol)
{
    if (n == 0) return TRI_ERR_DIM;
    if (!b || !d || !x_out) return TRI_ERR_NULL_PTR;
    if ((n > 1) && (!a || !c)) return TRI_ERR_NULL_PTR;

    /* n==1 特判 */
    if (n == 1) {
        if (tri_is_near_zero(b[0], tol)) return TRI_ERR_ZERO_PIVOT;
        x_out[0] = d[0] / b[0];
        return TRI_SUCCESS;
    }

    /* 工作数组：c'(0..n-2)，d'(0..n-1)
     * 说明：传统公式中 c'[n-1] 不用；这里单独开辟 cprime(n-1)，dprime(n)
     */
    double* cprime = (double*)malloc((n - 1) * sizeof(double));
    double* dprime = (double*)malloc(n * sizeof(double));
    if (!cprime || !dprime) {
        free(cprime); free(dprime);
        return TRI_ERR_ALLOC;
    }

    /* i=0: 归一化第一行 */
    if (tri_is_near_zero(b[0], tol)) {
        free(cprime); free(dprime);
        return TRI_ERR_ZERO_PIVOT;
    }
    cprime[0] = c[0] / b[0];
    dprime[0] = d[0] / b[0];

    /* i=1..n-2: 递推 c', d'；i=n-1 时只递推 d' */
    for (size_t i = 1; i < n; ++i) {
        double denom = b[i] - a[i - 1] * cprime[i - 1];
        if (tri_is_near_zero(denom, tol)) {
            free(cprime); free(dprime);
            return TRI_ERR_ZERO_PIVOT; /* 或判为病态/奇异 */
        }
        if (i < n - 1)
            cprime[i] = c[i] / denom;
        dprime[i] = (d[i] - a[i - 1] * dprime[i - 1]) / denom;
    }

    /* 回代 */
    x_out[n - 1] = dprime[n - 1];
    for (size_t k = n - 1; k-- > 0; ) {
        x_out[k] = dprime[k] - cprime[k] * x_out[k + 1];
    }

    free(cprime);
    free(dprime);
    return TRI_SUCCESS;
}

/*===========================
 * 面向对象 API 实现
 *===========================*/
int tri_create(size_t n, tri_solver_t** out) {
    if (!out) return TRI_ERR_NULL_PTR;
    *out = NULL;
    if (n == 0) return TRI_ERR_DIM;

    tri_solver_t* s = (tri_solver_t*)calloc(1, sizeof(*s));
    if (!s) return TRI_ERR_ALLOC;
    s->n = n;
    s->tol = TRI_DEFAULT_TOL;
    s->has_coeffs = 0;
    s->has_rhs = 0;
    *out = s;
    return TRI_SUCCESS;
}

int tri_set_tolerance(tri_solver_t* solver, double tol) {
    if (!solver) return TRI_ERR_NULL_PTR;
    if (tol <= 0) tol = TRI_DEFAULT_TOL;
    solver->tol = tol;
    return TRI_SUCCESS;
}

int tri_get_tolerance(const tri_solver_t* solver, double* tol_out) {
    if (!solver || !tol_out) return TRI_ERR_NULL_PTR;
    *tol_out = solver->tol;
    return TRI_SUCCESS;
}

int tri_set_coeffs(tri_solver_t* solver,
                   const double* a, const double* b, const double* c) {
    if (!solver) return TRI_ERR_NULL_PTR;
    size_t n = solver->n;
    if (!b || (n > 1 && (!a || !c))) return TRI_ERR_NULL_PTR;

    /* 释放旧内存，重新拷贝 */
    free(solver->a); free(solver->b); free(solver->c);
    solver->a = solver->b = solver->c = NULL;

    if (n > 1) {
        solver->a = (double*)malloc((n - 1) * sizeof(double));
        solver->c = (double*)malloc((n - 1) * sizeof(double));
        if (!solver->a || !solver->c) {
            free(solver->a); free(solver->c);
            solver->a = solver->c = NULL;
            return TRI_ERR_ALLOC;
        }
        memcpy(solver->a, a, (n - 1) * sizeof(double));
        memcpy(solver->c, c, (n - 1) * sizeof(double));
    }
    solver->b = (double*)malloc(n * sizeof(double));
    if (!solver->b) {
        free(solver->a); free(solver->c);
        solver->a = solver->c = NULL;
        return TRI_ERR_ALLOC;
    }
    memcpy(solver->b, b, n * sizeof(double));

    solver->has_coeffs = 1;
    return TRI_SUCCESS;
}

int tri_set_rhs(tri_solver_t* solver, const double* d) {
    if (!solver || !d) return TRI_ERR_NULL_PTR;
    size_t n = solver->n;

    free(solver->d);
    solver->d = (double*)malloc(n * sizeof(double));
    if (!solver->d) {
        solver->has_rhs = 0;
        return TRI_ERR_ALLOC;
    }
    memcpy(solver->d, d, n * sizeof(double));
    solver->has_rhs = 1;
    return TRI_SUCCESS;
}

int tri_solve(const tri_solver_t* solver, double* x_out) {
    if (!solver || !x_out) return TRI_ERR_NULL_PTR;
    if (!solver->has_coeffs || !solver->has_rhs) return TRI_ERR_BAD_STATE;

    return tri_thomas_core(solver->n,
                           solver->a, solver->b, solver->c, solver->d,
                           x_out, solver->tol);
}

int tri_get_n(const tri_solver_t* solver, size_t* n_out) {
    if (!solver || !n_out) return TRI_ERR_NULL_PTR;
    *n_out = solver->n;
    return TRI_SUCCESS;
}

int tri_destroy(tri_solver_t** psolver) {
    if (!psolver) return TRI_ERR_NULL_PTR;
    tri_solver_t* s = *psolver;
    if (!s) { *psolver = NULL; return TRI_SUCCESS; }
    free(s->a); free(s->b); free(s->c); free(s->d);
    s->a = s->b = s->c = s->d = NULL;
    free(s);
    *psolver = NULL;
    return TRI_SUCCESS;
}

/*===========================
 * 直接数组 API
 *===========================*/
int tri_solve_raw(size_t n,
                  const double* a,
                  const double* b,
                  const double* c,
                  const double* d,
                  double* x_out,
                  double tol)
{
    if (tol <= 0) tol = TRI_DEFAULT_TOL;
    return tri_thomas_core(n, a, b, c, d, x_out, tol);
}
