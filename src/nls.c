//
// Created by ninico on 2025/10/18.
//

#include "nls.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef NLS_DEFAULT_MAXIT
#define NLS_DEFAULT_MAXIT 50
#endif
#ifndef NLS_DEFAULT_TOL_RES
#define NLS_DEFAULT_TOL_RES 1e-10
#endif
#ifndef NLS_DEFAULT_TOL_STEP
#define NLS_DEFAULT_TOL_STEP 1e-12
#endif
#ifndef NLS_DEFAULT_LS_C
#define NLS_DEFAULT_LS_C 1e-4
#endif
#ifndef NLS_DEFAULT_LS_BETA
#define NLS_DEFAULT_LS_BETA 0.5
#endif
#ifndef NLS_DEFAULT_LS_MAXBT
#define NLS_DEFAULT_LS_MAXBT 20
#endif

struct nls_solver {
    size_t n;
    nls_fun_F F;
    nls_fun_J J;
    void* user;

    /* 选项 */
    int    max_iter;
    double tol_res, tol_step, fd_eps, ls_c, ls_beta;
    int    ls_max_backtrack;

    /* 工作区与统计 */
    double *f, *Jm, *dx, *x_try, *f_try;
    int iters;
    double last_res_inf, last_step_inf;
};

/*===========================
 * 工具函数
 *===========================*/
const char* nls_strerror(int err) {
    switch (err) {
        case NLS_SUCCESS:       return "success";
        case NLS_ERR_NULL_PTR:  return "null pointer";
        case NLS_ERR_ALLOC:     return "memory allocation failure";
        case NLS_ERR_DIM:       return "invalid dimension";
        case NLS_ERR_BAD_STATE: return "bad state";
        case NLS_ERR_SINGULAR:  return "singular or ill-conditioned linear system";
        case NLS_ERR_MAXIT:     return "max iterations reached (no convergence)";
        case NLS_ERR_JACOBIAN:  return "jacobian computation failed";
        case NLS_ERR_INTERNAL:  return "internal error";
        default:                return "unknown error";
    }
}

static double vnorm_inf(const double* x, size_t n) {
    double v = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double a = fabs(x[i]); if (a > v) v = a;
    }
    return v;
}

/* LU with partial pivoting, in-place on A (row-major n×n) */
static int lu_decomp(double* A, int n, int* piv) {
    for (int i = 0; i < n; ++i) piv[i] = i;
    for (int k = 0; k < n; ++k) {
        int p = k;
        double maxv = fabs(A[k*n + k]);
        for (int i = k+1; i < n; ++i) {
            double v = fabs(A[i*n + k]);
            if (v > maxv) { maxv = v; p = i; }
        }
        if (maxv <= DBL_EPSILON) return 0;
        if (p != k) {
            for (int j = 0; j < n; ++j) {
                double tmp = A[k*n + j]; A[k*n + j] = A[p*n + j]; A[p*n + j] = tmp;
            }
            int tp = piv[k]; piv[k] = piv[p]; piv[p] = tp;
        }
        for (int i = k+1; i < n; ++i) {
            A[i*n + k] /= A[k*n + k];
            double lik = A[i*n + k];
            for (int j = k+1; j < n; ++j) A[i*n + j] -= lik * A[k*n + j];
        }
    }
    return 1;
}
static void lu_solve(const double* LU, int n, const int* piv, double* b) {
    /* PB */
    double* y = (double*)malloc(n*sizeof(double));
    for (int i = 0; i < n; ++i) y[i] = b[piv[i]];
    /* Ly=y */
    for (int i = 1; i < n; ++i) {
        double s = y[i];
        for (int j = 0; j < i; ++j) s -= LU[i*n + j]*y[j];
        y[i] = s;
    }
    /* Ux=y */
    for (int i = n-1; i >= 0; --i) {
        double s = y[i];
        for (int j = i+1; j < n; ++j) s -= LU[i*n + j]*b[j];
        b[i] = s / LU[i*n + i];
    }
    free(y);
}

/* 数值差分雅可比（前向差分，列更新） */
static int fd_jacobian(nls_solver_t* s, const double* x, const double* f0, double* J) {
    size_t n = s->n;
    double* xt = s->x_try;
    double* ft = s->f_try;
    memcpy(xt, x, n*sizeof(double));
    for (size_t j = 0; j < n; ++j) {
        double h = s->fd_eps * (1.0 + fabs(x[j]));
        xt[j] += h;
        if (s->F(xt, ft, s->user) != 0) return NLS_ERR_JACOBIAN;
        for (size_t i = 0; i < n; ++i) {
            J[i*n + j] = (ft[i] - f0[i]) / h;
        }
        xt[j] -= h;
    }
    return NLS_SUCCESS;
}

/* 线搜索目标采用 φ(x)=0.5*||F(x)||^2 的单调下降准则
   使用 Armijo-like 条件：||F(x+αΔx)|| ≤ (1 - c α) ||F(x)|| */
static int backtracking(nls_solver_t* s, const double* x, const double* f, double fnorm,
                        const double* dx, double* x_out, double* f_out, double* fnorm_out,
                        double* alpha_used)
{
    size_t n = s->n;
    double alpha = 1.0;
    double best_norm = INFINITY;
    double* xt = s->x_try;
    double* ft = s->f_try;

    for (int k = 0; k < s->ls_max_backtrack; ++k) {
        for (size_t i = 0; i < n; ++i) xt[i] = x[i] + alpha * dx[i];
        if (s->F(xt, ft, s->user) != 0) return NLS_ERR_INTERNAL;
        double fn = vnorm_inf(ft, n);
        if (fn < best_norm) { best_norm = fn; memcpy(x_out, xt, n*sizeof(double)); memcpy(f_out, ft, n*sizeof(double)); }
        if (fn <= (1.0 - s->ls_c * alpha) * fnorm) {
            *fnorm_out = fn; *alpha_used = alpha; return NLS_SUCCESS;
        }
        alpha *= s->ls_beta;
    }
    /* 未满足强条件：取最佳者（保证单调下降） */
    *fnorm_out = best_norm;
    memcpy(x_out, xt, n*sizeof(double)); /* xt 当前是最后一次计算值，但 best_norm 可能不是最后一次 */
    memcpy(f_out, s->f_try, n*sizeof(double)); /* 上面已复制 best 时刻，故这两行仅作为兜底 */
    *alpha_used = alpha;
    return NLS_SUCCESS;
}

/*===========================
 * 对外 API
 *===========================*/
int nls_create(size_t n, nls_solver_t** out) {
    if (!out) return NLS_ERR_NULL_PTR;
    *out = NULL;
    if (n == 0) return NLS_ERR_DIM;
    nls_solver_t* s = (nls_solver_t*)calloc(1, sizeof(*s));
    if (!s) return NLS_ERR_ALLOC;
    s->n = n;
    s->max_iter = NLS_DEFAULT_MAXIT;
    s->tol_res  = NLS_DEFAULT_TOL_RES;
    s->tol_step = NLS_DEFAULT_TOL_STEP;
    s->fd_eps   = sqrt(DBL_EPSILON);
    s->ls_c     = NLS_DEFAULT_LS_C;
    s->ls_beta  = NLS_DEFAULT_LS_BETA;
    s->ls_max_backtrack = NLS_DEFAULT_LS_MAXBT;

    s->f     = (double*)malloc(n*sizeof(double));
    s->Jm    = (double*)malloc(n*n*sizeof(double));
    s->dx    = (double*)malloc(n*sizeof(double));
    s->x_try = (double*)malloc(n*sizeof(double));
    s->f_try = (double*)malloc(n*sizeof(double));
    if (!s->f || !s->Jm || !s->dx || !s->x_try || !s->f_try) {
        free(s->f); free(s->Jm); free(s->dx); free(s->x_try); free(s->f_try); free(s);
        return NLS_ERR_ALLOC;
    }
    *out = s;
    return NLS_SUCCESS;
}

int nls_set_function(nls_solver_t* s, nls_fun_F F, void* user) {
    if (!s) return NLS_ERR_NULL_PTR;
    s->F = F; s->user = user;
    return NLS_SUCCESS;
}
int nls_set_jacobian(nls_solver_t* s, nls_fun_J J) {
    if (!s) return NLS_ERR_NULL_PTR;
    s->J = J;
    return NLS_SUCCESS;
}
int nls_set_options(nls_solver_t* s, const nls_options_t* opt) {
    if (!s || !opt) return NLS_ERR_NULL_PTR;
    if (opt->max_iter > 0) s->max_iter = opt->max_iter;
    if (opt->tol_res  > 0) s->tol_res  = opt->tol_res;
    if (opt->tol_step > 0) s->tol_step = opt->tol_step;
    if (opt->fd_eps   > 0) s->fd_eps   = opt->fd_eps;
    if (opt->ls_c     > 0) s->ls_c     = opt->ls_c;
    if (opt->ls_beta  > 0 && opt->ls_beta < 1) s->ls_beta = opt->ls_beta;
    if (opt->ls_max_backtrack > 0) s->ls_max_backtrack = opt->ls_max_backtrack;
    return NLS_SUCCESS;
}

int nls_get_stats(const nls_solver_t* s, int* iters, double* res_inf, double* step_inf) {
    if (!s) return NLS_ERR_NULL_PTR;
    if (iters) *iters = s->iters;
    if (res_inf) *res_inf = s->last_res_inf;
    if (step_inf) *step_inf = s->last_step_inf;
    return NLS_SUCCESS;
}

int nls_destroy(nls_solver_t** ps) {
    if (!ps) return NLS_ERR_NULL_PTR;
    nls_solver_t* s = *ps;
    if (!s) { *ps = NULL; return NLS_SUCCESS; }
    free(s->f); free(s->Jm); free(s->dx); free(s->x_try); free(s->f_try);
    free(s);
    *ps = NULL;
    return NLS_SUCCESS;
}

int nls_solve(nls_solver_t* s, double* x) {
    if (!s || !x) return NLS_ERR_NULL_PTR;
    if (!s->F) return NLS_ERR_BAD_STATE;

    const size_t n = s->n;
    int* piv = (int*)malloc(n*sizeof(int));
    if (!piv) return NLS_ERR_ALLOC;

    /* 初始 F */
    if (s->F(x, s->f, s->user) != 0) { free(piv); return NLS_ERR_INTERNAL; }
    double fnorm = vnorm_inf(s->f, n);
    s->iters = 0;

    for (int it = 0; it < s->max_iter; ++it) {
        /* 收敛检查 */
        if (fnorm < s->tol_res) { s->iters = it; s->last_res_inf = fnorm; s->last_step_inf = 0.0; free(piv); return NLS_SUCCESS; }

        /* 雅可比 */
        int ret;
        if (s->J) {
            ret = s->J(x, s->Jm, s->user);
            if (ret) { free(piv); return NLS_ERR_JACOBIAN; }
        } else {
            ret = fd_jacobian(s, x, s->f, s->Jm);
            if (ret) { free(piv); return ret; }
        }

        /* 解 J dx = -F */
        for (size_t i = 0; i < n; ++i) s->dx[i] = -s->f[i];
        double* LU = s->Jm; /* 复用 J 存 LU */
        double* A  = (double*)malloc(n*n*sizeof(double));
        if (!A) { free(piv); return NLS_ERR_ALLOC; }
        memcpy(A, LU, n*n*sizeof(double));
        if (!lu_decomp(A, (int)n, piv)) { free(A); free(piv); return NLS_ERR_SINGULAR; }
        lu_solve(A, (int)n, piv, s->dx);
        free(A);

        double step_inf = vnorm_inf(s->dx, n);
        if (step_inf < s->tol_step) { s->iters = it; s->last_res_inf = fnorm; s->last_step_inf = step_inf; free(piv); return NLS_SUCCESS; }

        /* 回溯线搜索 */
        double alpha_used=0.0, fnorm_new=0.0;
        double* xnew = s->x_try;
        double* fnew = s->f_try;
        ret = backtracking(s, x, s->f, fnorm, s->dx, xnew, fnew, &fnorm_new, &alpha_used);
        if (ret) { free(piv); return ret; }

        /* 接受并迭代 */
        memcpy(x, xnew, n*sizeof(double));
        memcpy(s->f, fnew, n*sizeof(double));
        fnorm = fnorm_new;

        s->last_res_inf = fnorm;
        s->last_step_inf = step_inf;
        s->iters = it + 1;
    }
    free(piv);
    return NLS_ERR_MAXIT;
}

int nls_solve_raw(size_t n,
                  nls_fun_F F, nls_fun_J J, void* user,
                  const nls_options_t* opt,
                  double* x_inout)
{
    if (!F || !x_inout) return NLS_ERR_NULL_PTR;
    nls_solver_t* s = NULL;
    int err = nls_create(n, &s);
    if (err) return err;
    nls_set_function(s, F, user);
    nls_set_jacobian(s, J);
    if (opt) nls_set_options(s, opt);
    err = nls_solve(s, x_inout);
    nls_destroy(&s);
    return err;
}
