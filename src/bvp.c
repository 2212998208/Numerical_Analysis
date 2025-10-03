//
// Created by ninico on 2025/10/18.
//
#include "bvp.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#ifndef BVP_DEFAULT_MAXIT
#define BVP_DEFAULT_MAXIT 50
#endif
#ifndef BVP_DEFAULT_TOL_RES
#define BVP_DEFAULT_TOL_RES 1e-8
#endif
#ifndef BVP_DEFAULT_TOL_STEP
#define BVP_DEFAULT_TOL_STEP 1e-10
#endif

struct bvp_solver {
    size_t m, N;     /* 维度与内点数 */
    double a, b, h;  /* 区间与步长 */
    /* 回调 */
    bvp_fun_G G;
    bvp_fun_J J_y;
    bvp_fun_J J_yp;
    void* user;

    /* 数据：边界与内点当前解 */
    double* ya;      /* 长度 m */
    double* yb;      /* 长度 m */
    double* y;       /* 内点解，长度 N*m，按节点块存储 */

    /* 选项 */
    int max_iter;
    double tol_res, tol_step;

    int has_boundary;
    int has_init;
};

/*===========================
 * 文本错误信息
 *===========================*/
const char* bvp_strerror(int err) {
    switch (err) {
        case BVP_SUCCESS:       return "success";
        case BVP_ERR_NULL_PTR:  return "null pointer";
        case BVP_ERR_ALLOC:     return "memory allocation failure";
        case BVP_ERR_DIM:       return "invalid dimension";
        case BVP_ERR_BAD_STATE: return "bad state";
        case BVP_ERR_SINGULAR:  return "singular or ill-conditioned linear system";
        case BVP_ERR_MAXIT:     return "max iterations reached (no convergence)";
        case BVP_ERR_JACOBIAN:  return "jacobian computation failed";
        case BVP_ERR_INTERNAL:  return "internal error";
        default:                return "unknown error";
    }
}

/*===========================
 * 小型矩阵运算（m×m），带 LU 主元
 *===========================*/
static int lu_decomp(double* A, int n, int* piv) { /* A in-place -> LU */
    for (int i = 0; i < n; ++i) piv[i] = i;
    for (int k = 0; k < n; ++k) {
        /* 选主元 */
        int p = k;
        double maxv = fabs(A[k*n + k]);
        for (int i = k+1; i < n; ++i) {
            double v = fabs(A[i*n + k]);
            if (v > maxv) { maxv = v; p = i; }
        }
        if (maxv <= DBL_EPSILON) return 0; /* singular */
        if (p != k) {
            for (int j = 0; j < n; ++j) {
                double tmp = A[k*n + j];
                A[k*n + j] = A[p*n + j];
                A[p*n + j] = tmp;
            }
            int tp = piv[k]; piv[k] = piv[p]; piv[p] = tp;
        }
        /* 消元 */
        for (int i = k+1; i < n; ++i) {
            A[i*n + k] /= A[k*n + k];
            double lik = A[i*n + k];
            for (int j = k+1; j < n; ++j) {
                A[i*n + j] -= lik * A[k*n + j];
            }
        }
    }
    return 1;
}

static void lu_solve(const double* LU, int n, const int* piv, double* b, int nrhs) {
    /* 行置换 */
    for (int k = 0; k < n; ++k) {
        int pk = piv[k];
        if (pk != k) {
            for (int j = 0; j < nrhs; ++j) {
                double tmp = b[k*nrhs + j];
                b[k*nrhs + j] = b[pk*nrhs + j];
                b[pk*nrhs + j] = tmp;
            }
        }
    }
    /* Ly = Pb */
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double lij = LU[i*n + j];
            for (int k = 0; k < nrhs; ++k) {
                b[i*nrhs + k] -= lij * b[j*nrhs + k];
            }
        }
    }
    /* Ux = y */
    for (int i = n-1; i >= 0; --i) {
        double uii = LU[i*n + i];
        for (int k = 0; k < nrhs; ++k) {
            b[i*nrhs + k] /= uii;
        }
        for (int j = 0; j < i; ++j) {
            double uji = LU[j*n + i];
            for (int k = 0; k < nrhs; ++k) {
                b[j*nrhs + k] -= uji * b[i*nrhs + k];
            }
        }
    }
}

static void mat_eye(double* M, size_t m) {
    memset(M, 0, m*m*sizeof(double));
    for (size_t i = 0; i < m; ++i) M[i*m + i] = 1.0;
}

static void mat_copy(double* dst, const double* src, size_t m) {
    memcpy(dst, src, m*m*sizeof(double));
}
static void mat_add_inplace(double* A, const double* B, size_t m, double alpha) {
    size_t L = m*m;
    for (size_t i = 0; i < L; ++i) A[i] += alpha * B[i];
}
static void mat_scale(double* A, size_t m, double s) {
    size_t L = m*m;
    for (size_t i = 0; i < L; ++i) A[i] *= s;
}
static void mat_mul(const double* A, const double* B, double* C, size_t m) {
    /* C = A B */
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < m; ++k) sum += A[i*m+k]*B[k*m+j];
            C[i*m+j] = sum;
        }
    }
}
static void mat_vec(const double* A, const double* x, double* y, size_t m) {
    for (size_t i = 0; i < m; ++i) {
        double s = 0.0;
        for (size_t j = 0; j < m; ++j) s += A[i*m+j]*x[j];
        y[i] = s;
    }
}
static void vec_axpy(double* y, const double* x, size_t m, double alpha) {
    for (size_t i = 0; i < m; ++i) y[i] += alpha * x[i];
}
static double vec_norm_inf(const double* x, size_t L) {
    double v = 0.0;
    for (size_t i = 0; i < L; ++i) {
        double a = fabs(x[i]);
        if (a > v) v = a;
    }
    return v;
}

/*===========================
 * 数值差分雅可比（前向差分）
 *===========================*/
static int fd_J(size_t m, double eps_base,
                bvp_fun_G G, void* user,
                double x, const double* y, const double* yp,
                int need_Jy, double* Jy,
                int need_Jyp, double* Jyp)
{
    int ret;
    double* g0 = (double*)malloc(m*sizeof(double));
    double* g1 = (double*)malloc(m*sizeof(double));
    double* ytmp = (double*)malloc(m*sizeof(double));
    double* yptmp= (double*)malloc(m*sizeof(double));
    if (!g0 || !g1 || !ytmp || !yptmp) {
        free(g0); free(g1); free(ytmp); free(yptmp);
        return BVP_ERR_ALLOC;
    }
    memcpy(ytmp, y, m*sizeof(double));
    memcpy(yptmp, yp, m*sizeof(double));
    ret = G(x, ytmp, yptmp, g0, user);
    if (ret != 0) { free(g0); free(g1); free(ytmp); free(yptmp); return BVP_ERR_JACOBIAN; }

    if (need_Jy) {
        for (size_t k = 0; k < m; ++k) {
            double sk = eps_base * (1.0 + fabs(y[k]));
            ytmp[k] += sk;
            ret = G(x, ytmp, yptmp, g1, user);
            if (ret != 0) { free(g0); free(g1); free(ytmp); free(yptmp); return BVP_ERR_JACOBIAN; }
            for (size_t i = 0; i < m; ++i) {
                Jy[i*m + k] = (g1[i] - g0[i]) / sk; /* 第 k 列 */
            }
            ytmp[k] -= sk;
        }
    }
    if (need_Jyp) {
        for (size_t k = 0; k < m; ++k) {
            double sk = eps_base * (1.0 + fabs(yp[k]));
            yptmp[k] += sk;
            ret = G(x, ytmp, yptmp, g1, user);
            if (ret != 0) { free(g0); free(g1); free(ytmp); free(yptmp); return BVP_ERR_JACOBIAN; }
            for (size_t i = 0; i < m; ++i) {
                Jyp[i*m + k] = (g1[i] - g0[i]) / sk;
            }
            yptmp[k] -= sk;
        }
    }
    free(g0); free(g1); free(ytmp); free(yptmp);
    return BVP_SUCCESS;
}

/*===========================
 * 构建残差与块三对角雅可比
 * r_i = y_{i-1} - 2 y_i + y_{i+1} - h^2 G(x_i, y_i, y'_i)
 * A_i = I + (h/2) J_yp,  B_i = -2I - h^2 J_y,  C_i = I - (h/2) J_yp
 * 注意：未知仅为内点 y_1..y_N，故 i=0 无 A_0，i=N-1 无 C_{N-1}
 *===========================*/
static int build_residual_and_blocks(bvp_solver_t* s,
                                     const double* y_in,
                                     double* r,       /* 长度 N*m */
                                     double* Ablk,    /* (N-1)×(m×m) */
                                     double* Bblk,    /*  N   ×(m×m) */
                                     double* Cblk)    /* (N-1)×(m×m) */
{
    size_t m = s->m, N = s->N;
    double h = s->h;
    double *yi   = (double*)malloc(m*sizeof(double));
    double *yim1 = (double*)malloc(m*sizeof(double));
    double *yip1 = (double*)malloc(m*sizeof(double));
    double *ypi  = (double*)malloc(m*sizeof(double));
    double *Gi   = (double*)malloc(m*sizeof(double));
    double *Jy   = (double*)malloc(m*m*sizeof(double));
    double *Jyp  = (double*)malloc(m*m*sizeof(double));
    if (!yi || !yim1 || !yip1 || !ypi || !Gi || !Jy || !Jyp) {
        free(yi); free(yim1); free(yip1); free(ypi); free(Gi); free(Jy); free(Jyp);
        return BVP_ERR_ALLOC;
    }
    for (size_t i = 0; i < N; ++i) {
        /* 取三点与端点 */
        const double xi = s->a + (i+1)*h;

        const double* y_i = &y_in[i*m];
        memcpy(yi, y_i, m*sizeof(double));

        if (i == 0) memcpy(yim1, s->ya, m*sizeof(double));
        else        memcpy(yim1, &y_in[(i-1)*m], m*sizeof(double));
        if (i+1 < N) memcpy(yip1, &y_in[(i+1)*m], m*sizeof(double));
        else         memcpy(yip1, s->yb, m*sizeof(double));

        for (size_t k = 0; k < m; ++k) {
            ypi[k] = (yip1[k] - yim1[k]) / (2.0*h);
        }

        /* G 与雅可比 */
        int ret = s->G(xi, yi, ypi, Gi, s->user);
        if (ret) { free(yi); free(yim1); free(yip1); free(ypi); free(Gi); free(Jy); free(Jyp); return BVP_ERR_INTERNAL; }

        int have_Jy  = (s->J_y  != NULL);
        int have_Jyp = (s->J_yp != NULL);

        if (have_Jy)  { ret = s->J_y (xi, yi, ypi, Jy,  s->user); if (ret) return BVP_ERR_JACOBIAN; }
        if (have_Jyp) { ret = s->J_yp(xi, yi, ypi, Jyp, s->user); if (ret) return BVP_ERR_JACOBIAN; }

        if (!have_Jy || !have_Jyp) {
            double eps_base = sqrt(DBL_EPSILON);
            ret = fd_J(m, eps_base, s->G, s->user, xi, yi, ypi,
                       !have_Jy,  Jy,
                       !have_Jyp, Jyp);
            if (ret) { free(yi); free(yim1); free(yip1); free(ypi); free(Gi); free(Jy); free(Jyp); return ret; }
        }

        /* 残差：r_i */
        for (size_t k = 0; k < m; ++k) {
            r[i*m + k] = (yim1[k] - 2.0*yi[k] + yip1[k]) - h*h * Gi[k];
        }

        /* 组块：A_i, B_i, C_i */
        double *Bi = &Bblk[i*m*m];
        for (size_t t = 0; t < m*m; ++t) Bi[t] = - h*h * Jy[t];
        for (size_t k = 0; k < m; ++k) Bi[k*m + k] += -2.0;

        if (i > 0) {
            double *Ai = &Ablk[(i-1)*m*m];
            memcpy(Ai, Jyp, m*m*sizeof(double));
            mat_scale(Ai, m, (h/2.0));
            for (size_t k = 0; k < m; ++k) Ai[k*m + k] += 1.0;
        }
        if (i+1 < N) {
            double *Ci = &Cblk[i*m*m];
            memcpy(Ci, Jyp, m*m*sizeof(double));
            mat_scale(Ci, m, -(h/2.0));
            for (size_t k = 0; k < m; ++k) Ci[k*m + k] += 1.0;
        }
    }
    free(yi); free(yim1); free(yip1); free(ypi); free(Gi); free(Jy); free(Jyp);
    return BVP_SUCCESS;
}

/*===========================
 * 块三对角 Thomas（矩阵版）
 * 解： (A_i) x_{i-1} + (B_i) x_i + (C_i) x_{i+1} = d_i
 * i=0..N-1，其中 A_0, C_{N-1} 不存在
 *===========================*/
static int block_tridiag_solve(size_t N, size_t m,
                               const double* Ablk, const double* Bblk, const double* Cblk,
                               double* d, /* inout: RHS -> 解（每块向量长 m） */
                               double* workLU, int* workPiv, double* workMat)
{
    /* 预分配：
       workLU: N * (m×m)
       workPiv: N * m
       workMat: max(2, m) ×(m×m)，临时 */
    /* i=0 */
    memcpy(&workLU[0], &Bblk[0], m*m*sizeof(double));
    if (!lu_decomp(&workLU[0], (int)m, &workPiv[0])) return BVP_ERR_SINGULAR;

    /* C'_0 = solve(B0, C0), d'_0 = solve(B0, d0) */
    if (N > 1) {
        memcpy(workMat, &Cblk[0], m*m*sizeof(double));
        lu_solve(&workLU[0], (int)m, &workPiv[0], workMat, (int)m); /* m 右端列 */
        memcpy((double*)&Cblk[0], workMat, m*m*sizeof(double)); /* 把 C' 写回 Cblk (就地复用) */
    }
    lu_solve(&workLU[0], (int)m, &workPiv[0], &d[0], 1);

    for (size_t i = 1; i < N; ++i) {
        /* S_i = B_i - A_i * C'_{i-1} */
        const double* Ai = &Ablk[(i-1)*m*m];
        const double* Bi = &Bblk[i*m*m];
        const double* Cim1p = &Cblk[(i-1)*m*m];
        double* Si = &workLU[i*m*m];

        /* Si = Bi - Ai*C'_{i-1} */
        mat_mul(Ai, Cim1p, workMat, m);
        memcpy(Si, Bi, m*m*sizeof(double));
        mat_add_inplace(Si, workMat, m, -1.0);

        if (!lu_decomp(Si, (int)m, &workPiv[i*m])) return BVP_ERR_SINGULAR;

        /* d_i := d_i - A_i*d'_{i-1} */
        double* di = &d[i*m];
        double* dim1p = &d[(i-1)*m];
        mat_vec(Ai, dim1p, workMat, m);
        for (size_t k = 0; k < m; ++k) di[k] -= workMat[k];

        /* C'_i = solve(S_i, C_i) */
        if (i+1 < N) {
            memcpy(workMat, &Cblk[i*m*m], m*m*sizeof(double));
            lu_solve(Si, (int)m, &workPiv[i*m], workMat, (int)m);
            memcpy((double*)&Cblk[i*m*m], workMat, m*m*sizeof(double));
        }
        /* d'_i = solve(S_i, d_i) */
        lu_solve(Si, (int)m, &workPiv[i*m], di, 1);
    }

    /* 回代：x_{N-1} 已在 d 中（即 d'_{N-1}） */
    for (long i = (long)N - 2; i >= 0; --i) {
        /* d_i := d'_i - C'_i * d'_{i+1} */
        double* di = &d[i*m];
        double* dip1 = &d[(i+1)*m];
        mat_vec(&Cblk[i*m*m], dip1, workMat, m);
        for (size_t k = 0; k < m; ++k) di[k] -= workMat[k];
    }
    return BVP_SUCCESS;
}

/*===========================
 * 组装并执行一次牛顿步（阻尼线搜索）
 *===========================*/
static int newton_solve(bvp_solver_t* s, double* y_interior) {
    size_t m = s->m, N = s->N;
    double h = s->h;

    /* 工作区 */
    double* r    = (double*)malloc(N*m*sizeof(double));
    double* Ablk = (double*)malloc((N>1 ? (N-1) : 1)*m*m*sizeof(double));
    double* Bblk = (double*)malloc( N   *m*m*sizeof(double));
    double* Cblk = (double*)malloc((N>1 ? (N-1) : 1)*m*m*sizeof(double));
    double* d    = (double*)malloc(N*m*sizeof(double));   /* -r */
    double* ynew = (double*)malloc(N*m*sizeof(double));
    double* LU   = (double*)malloc(N*m*m*sizeof(double));
    int*    piv  = (int*)   malloc(N*m*sizeof(int));
    double* MAT  = (double*)malloc(m*m*sizeof(double));
    if (!r||!Ablk||!Bblk||!Cblk||!d||!ynew||!LU||!piv||!MAT) {
        free(r); free(Ablk); free(Bblk); free(Cblk); free(d);
        free(ynew); free(LU); free(piv); free(MAT);
        return BVP_ERR_ALLOC;
    }

    int it = 0;
    double last_res = INFINITY;
    while (1) {
        int ret = build_residual_and_blocks(s, y_interior, r, Ablk, Bblk, Cblk);
        if (ret) { free(r); free(Ablk); free(Bblk); free(Cblk); free(d); free(ynew); free(LU); free(piv); free(MAT); return ret; }

        for (size_t i = 0; i < N*m; ++i) d[i] = -r[i];

        ret = block_tridiag_solve(N, m, Ablk, Bblk, Cblk, d, LU, piv, MAT);
        if (ret) { free(r); free(Ablk); free(Bblk); free(Cblk); free(d); free(ynew); free(LU); free(piv); free(MAT); return ret; }

        double res_norm = vec_norm_inf(r, N*m);
        double step_norm = vec_norm_inf(d, N*m);

        if (res_norm < s->tol_res || step_norm < s->tol_step) {
            break; /* 收敛 */
        }
        if (it >= s->max_iter) {
            free(r); free(Ablk); free(Bblk); free(Cblk); free(d); free(ynew); free(LU); free(piv); free(MAT);
            return BVP_ERR_MAXIT;
        }

        /* 阻尼线搜索：尝试 alpha = 1, 1/2, 1/4, ... */
        double alpha = 1.0;
        double best_res = INFINITY;
        double best_alpha = 0.0;
        int improved = 0;
        for (int ls = 0; ls < 8; ++ls) {
            for (size_t i = 0; i < N*m; ++i) ynew[i] = y_interior[i] + alpha * d[i];

            /* 计算新残差范数 */
            int ret2 = build_residual_and_blocks(s, ynew, r, Ablk, Bblk, Cblk);
            if (ret2) { free(r); free(Ablk); free(Bblk); free(Cblk); free(d); free(ynew); free(LU); free(piv); free(MAT); return ret2; }
            double rn = vec_norm_inf(r, N*m);
            if (rn < best_res) { best_res = rn; best_alpha = alpha; }
            if (rn < res_norm * 0.8) { improved = 1; break; }
            alpha *= 0.5;
        }
        if (!improved && best_alpha == 0.0) best_alpha = 0.5; /* 兜底 */

        for (size_t i = 0; i < N*m; ++i) y_interior[i] += best_alpha * d[i];

        last_res = res_norm;
        ++it;
        (void)last_res; (void)h;
    }

    free(r); free(Ablk); free(Bblk); free(Cblk); free(d); free(ynew); free(LU); free(piv); free(MAT);
    return BVP_SUCCESS;
}

/*===========================
 * 对外 API
 *===========================*/
int bvp_create(size_t m, size_t N, double a, double b, bvp_solver_t** out) {
    if (!out) return BVP_ERR_NULL_PTR;
    *out = NULL;
    if (m == 0 || N == 0) return BVP_ERR_DIM;
    if (!(b > a)) return BVP_ERR_DIM;

    bvp_solver_t* s = (bvp_solver_t*)calloc(1, sizeof(*s));
    if (!s) return BVP_ERR_ALLOC;
    s->m = m; s->N = N; s->a = a; s->b = b; s->h = (b-a)/(N+1);
    s->max_iter = BVP_DEFAULT_MAXIT;
    s->tol_res  = BVP_DEFAULT_TOL_RES;
    s->tol_step = BVP_DEFAULT_TOL_STEP;
    s->ya = (double*)malloc(m*sizeof(double));
    s->yb = (double*)malloc(m*sizeof(double));
    s->y  = (double*)malloc(N*m*sizeof(double));
    if (!s->ya || !s->yb || !s->y) {
        free(s->ya); free(s->yb); free(s->y); free(s);
        return BVP_ERR_ALLOC;
    }
    s->has_boundary = 0;
    s->has_init = 0;
    *out = s;
    return BVP_SUCCESS;
}

int bvp_set_function(bvp_solver_t* s, bvp_fun_G G, void* user) {
    if (!s) return BVP_ERR_NULL_PTR;
    s->G = G; s->user = user;
    return BVP_SUCCESS;
}
int bvp_set_jacobians(bvp_solver_t* s, bvp_fun_J J_y, bvp_fun_J J_yp) {
    if (!s) return BVP_ERR_NULL_PTR;
    s->J_y = J_y; s->J_yp = J_yp;
    return BVP_SUCCESS;
}
int bvp_set_boundary(bvp_solver_t* s, const double* ya, const double* yb) {
    if (!s || !ya || !yb) return BVP_ERR_NULL_PTR;
    memcpy(s->ya, ya, s->m*sizeof(double));
    memcpy(s->yb, yb, s->m*sizeof(double));
    s->has_boundary = 1;
    return BVP_SUCCESS;
}
int bvp_set_initial_guess(bvp_solver_t* s, const double* y_interior) {
    if (!s) return BVP_ERR_NULL_PTR;
    if (!y_interior) return BVP_ERR_NULL_PTR;
    memcpy(s->y, y_interior, s->N*s->m*sizeof(double));
    s->has_init = 1;
    return BVP_SUCCESS;
}
int bvp_set_options(bvp_solver_t* s, int max_iter, double tol_res, double tol_step) {
    if (!s) return BVP_ERR_NULL_PTR;
    if (max_iter > 0) s->max_iter = max_iter;
    if (tol_res  > 0) s->tol_res  = tol_res;
    if (tol_step > 0) s->tol_step = tol_step;
    return BVP_SUCCESS;
}
int bvp_get_grid(const bvp_solver_t* s, double* x_all_out) {
    if (!s || !x_all_out) return BVP_ERR_NULL_PTR;
    for (size_t i = 0; i < s->N+2; ++i) x_all_out[i] = s->a + i*s->h;
    return BVP_SUCCESS;
}
int bvp_get_meta(const bvp_solver_t* s, size_t* m, size_t* N, double* a, double* b) {
    if (!s) return BVP_ERR_NULL_PTR;
    if (m) *m = s->m; if (N) *N = s->N; if (a) *a = s->a; if (b) *b = s->b;
    return BVP_SUCCESS;
}
int bvp_destroy(bvp_solver_t** ps) {
    if (!ps) return BVP_ERR_NULL_PTR;
    bvp_solver_t* s = *ps;
    if (!s) { *ps = NULL; return BVP_SUCCESS; }
    free(s->ya); free(s->yb); free(s->y);
    free(s);
    *ps = NULL;
    return BVP_SUCCESS;
}

int bvp_solve(bvp_solver_t* s, double* y_all_out) {
    if (!s || !y_all_out) return BVP_ERR_NULL_PTR;
    if (!s->G) return BVP_ERR_BAD_STATE;
    if (!s->has_boundary) return BVP_ERR_BAD_STATE;

    /* 初值：若用户未提供，则线性插值 */
    if (!s->has_init) {
        for (size_t i = 0; i < s->N; ++i) {
            double t = (double)(i+1)/(s->N+1);
            for (size_t k = 0; k < s->m; ++k) {
                s->y[i*s->m + k] = (1.0 - t)*s->ya[k] + t*s->yb[k];
            }
        }
    }
    int ret = newton_solve(s, s->y);
    if (ret) return ret;

    /* 输出全节点：端点 + 内点 */
    memcpy(&y_all_out[0], s->ya, s->m*sizeof(double));
    for (size_t i = 0; i < s->N; ++i)
        memcpy(&y_all_out[(i+1)*s->m], &s->y[i*s->m], s->m*sizeof(double));
    memcpy(&y_all_out[(s->N+1)*s->m], s->yb, s->m*sizeof(double));
    return BVP_SUCCESS;
}

int bvp_solve_raw(size_t m, size_t N, double a, double b,
                  const double* ya, const double* yb,
                  bvp_fun_G G, bvp_fun_J J_y, bvp_fun_J J_yp, void* user,
                  const bvp_options_t* opt,
                  double* y_all_out)
{
    if (!ya || !yb || !G || !y_all_out) return BVP_ERR_NULL_PTR;
    bvp_solver_t* s = NULL;
    int ret = bvp_create(m, N, a, b, &s);
    if (ret) return ret;
    bvp_set_function(s, G, user);
    bvp_set_jacobians(s, J_y, J_yp);
    bvp_set_boundary(s, ya, yb);
    if (opt) bvp_set_options(s, opt->max_iter, opt->tol_res, opt->tol_step);
    ret = bvp_solve(s, y_all_out);
    bvp_destroy(&s);
    return ret;
}
