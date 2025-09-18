#include "newton_raphson.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>


struct NonLinearFunction {
    double (*f)(double x);
};

struct NonLinearRange {
    struct NonLinearFunction func;
    double x0;           // 初始猜测值
    double tol;         // 容差
    size_t max_iter;    // 最大迭代次数
    char *name;         // 方程名称
};



// 计算点导数（设置Δx为1e-3）
static NewtonRaphson_Err derivative(double (*f)(double x), const double x0, long double *outDy2dx) {
    const double x1 = x0 + 1e-3;
    const double y1 = f(x1);
    const double y0 = f(x0);
    const long double dy2dx_neg = (y1 - y0) / (x1 - x0);
    const long double dy2dx_pos = (y0 - y1) / (x0 - x1);
    if (fabsl(dy2dx_neg) < 1e-9 || fabsl(dy2dx_pos) < 1e-9) {
        return NEWTON_RAPHSON_ERR_DERIVATIVE_ZERO;
    }

    if (fabsl(dy2dx_neg - dy2dx_pos) > 1e-9) {
        return NEWTON_RAPHSON_ERR_DERIVATIVE_UNSTABLE;
    }
    *outDy2dx = (dy2dx_neg + dy2dx_pos) / 2.0L;
    return NEWTON_RAPHSON_OK;
}


// 牛顿-拉夫森法求解非线性方程
NewtonRaphson_Err newton_raphson_solve(const NonLinearRange *outRange, double *outRoot) {
    if (!*outRange || (*outRange)->tol == 0 || (*outRange)->max_iter == 0) {
        return NEWTON_RAPHSON_ERR_INVALID;
    }

    double (*f)(double x) = (*outRange)->func.f;
    double *x0 = &(*outRange)->x0;
    double *x1 = (double *)malloc(sizeof(double));
    long double dy2dx = 0.0;

    for (size_t i = 0; i < (*outRange)->max_iter; i++) {

        NewtonRaphson_Err err = derivative(f, *x0, &dy2dx);
        if (err == NEWTON_RAPHSON_ERR_DERIVATIVE_ZERO || err == NEWTON_RAPHSON_ERR_DERIVATIVE_UNSTABLE) {
            free(x1);
            return err;
        }
        *x1 = *x0 - f(*x0) / dy2dx;
        if (!*x1) {
            free(x1);
            return NEWTON_RAPHSON_ERR_INVALID;
        }
        if (fabsl(*x1 - *x0) < (*outRange)->tol || fabs(f(*x1)) < (*outRange)->tol) {
            *outRoot = *x1;
            free(x1);
            return NEWTON_RAPHSON_OK;
        }
        *x0 = *x1;
    }
    free(x1);
    return NEWTON_RAPHSON_ERR_MAXITER;

}


// 构造函数
NewtonRaphson_Err NonLinearRange_create(double (*f)(double x), const double x0,
    const double tol, const size_t max_iter,
    NonLinearRange *outRange,
    const char *name) {
    if (!outRange || !f || tol <= 1e-15 || max_iter == 0 || !name) {
        return NEWTON_RAPHSON_ERR_INVALID;
    }
    NonLinearRange range = (NonLinearRange)malloc(sizeof(struct NonLinearRange));
    if (!range) {
        return NEWTON_RAPHSON_ERR_NOMEM;
    }
    range->func.f = f;
    range->x0 = x0;
    range->tol = tol;
    range->max_iter = max_iter;
    range->name = (char *)name;
    *outRange = range;
    return NEWTON_RAPHSON_OK;

}


// 析构函数
NewtonRaphson_Err NonLinearRange_destroy(const NonLinearRange *inRange) {
    if (!inRange || !*inRange) {
        return NEWTON_RAPHSON_ERR_INVALID;
    }
    free(*inRange);
    return NEWTON_RAPHSON_OK;
}


// 公共API
const NewtonRaphsonAPI NewtonRaphson = {
    .NonLinearRange_create = NonLinearRange_create,
    .NonLinearRange_destroy = NonLinearRange_destroy,
    .newton_raphson_solve = newton_raphson_solve
};