#include "successive_approximation.h"

#include <math.h>
#include <stdlib.h>

struct FixedPointFunction {
    double (*g)(double x);
};


struct NonlinearSA {
    struct FixedPointFunction gfunc;
    double x0;
    double tol;
    size_t max_iter;
    char *name;
};


// 收敛性定理（仅验证选择的初始点,设置Δx为1e-3）
static Succesive_Err contraction_mapping(double (*g)(double x), const double x0) {
    const double x1 = x0 + 1e-3;
    const double y0 = g(x0);
    const double y1 = g(x1);
    const long double dy2dx_neg = (y1 - y0) / (x1 - x0);
    const long double dy2dy_pos = (y0 - y1) / (x0 - x1);

    const long double Dy2dx = (dy2dy_pos + dy2dx_neg) / 2.0L;
    if (fabsl(Dy2dx) >= 0 && fabsl(Dy2dx) < 1) {
        return SA_OK;
    }
    return SA_ERR_NOAPPROXIMATION;
}


// 逐次逼近法（不动点迭代法）
Succesive_Err nonlinear_sa_solve(const NonlinearSA *outSA, double *outRoot) {
    if (!*outSA || (*outSA)->tol == 0 || (*outSA)->max_iter == 0) {
        return SA_ERR_INVAL;
    }
    double (*g)(double x) = (*outSA)->gfunc.g;
    double *x0 = &(*outSA)->x0;
    double *x1 = (double *)malloc(sizeof(double));

    for (size_t i = 0; i < (*outSA)->max_iter; i++) {
        *x1 = g(*x0);
        Succesive_Err err = contraction_mapping(g, *x1);
        if (!*x1) {
            free(x1);
            return SA_ERR_INVAL;
        }
        if (err == SA_ERR_NOAPPROXIMATION) {
            free(x1);
            return err;
        }
        if (fabsl(*x1 - *x0) < (*outSA)->tol) {
            if (fabs(*x1) < 1e-7) { // 根接近于0
                *outRoot = 0;
                free(x1);
                return SA_OK;
            }
            *outRoot = *x1;
            free(x1);
            return SA_OK;
        }
        *x0 = *x1;
    }
    free(x1);
    return SA_ERR_MAXITER;
}






// 构造函数
Succesive_Err NonlinearSA_create(double (*g)(double x), const double x0,
    const double tol, const size_t max_iter,
    NonlinearSA *outSA,
    const char *name) {
    if (!outSA || !g || tol <= 1e-65 || max_iter == 0 || !name) {
        return SA_ERR_INVAL;
    }
    const Succesive_Err err = contraction_mapping(g, x0);
    if (err != SA_OK) {
        return err;
    }
    const NonlinearSA sa = (NonlinearSA)malloc(sizeof(struct NonlinearSA));
    if (!sa) {
        return SA_ERR_NOMEM;
    }
    sa->gfunc.g = g;
    sa->x0 = x0;
    sa->tol = tol;
    sa->max_iter = max_iter;
    sa->name = (char *)name;
    *outSA = sa;
    return SA_OK;
}







// 析构函数
Succesive_Err NonlinearSA_destroy(const NonlinearSA *inSA) {
    if (!*inSA) {
        return SA_ERR_INVAL;
    }
    free(*inSA);
    inSA = NULL;
    return SA_OK;
}






// 公共API
const SAAPI SA = {
    .nonlinear_sa_solve = nonlinear_sa_solve,
    .NonlinearSA_create = NonlinearSA_create,
    .NonlinearSA_destroy = NonlinearSA_destroy
};



