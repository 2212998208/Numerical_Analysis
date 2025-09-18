#include "secant.h"

#include <math.h>
#include <stdlib.h>


struct NonLinearFunction {
    double (*f)(double x);
};

struct NonLinearScant {
    struct NonLinearFunction func;
    double x0;           // 初始猜测值
    double x1;           // 第二个猜测值
    double tol;         // 容差
    size_t max_iter;    // 最大迭代次数
    char *name;         // 方程名称
};


// 计算割线斜率
static Secant_Err secant_slope(double (*f)(double x), const double x0, const double x1, long double *outDy2dx) {
    if (fabs(x1 - x0) < 1e-9) {
        return SECANT_ERR_DIVIDE_BY_ZERO;
    }

    const double y1 = f(x1);
    const double y0 = f(x0);

    const double dy2dx = (y1 - y0) / (x1 - x0);
    if (fabs(dy2dx) < 1e-9) {
        return SECANT_ERR_DIVIDE_BY_ZERO;
    }

    *outDy2dx = dy2dx;
    return SECANT_OK;
}

// 割线法求解非线性方程
Secant_Err secant_solve(const NonLinearScant *outScant, double *outRoot) {
    if (!*outScant || (*outScant)->tol == 0 || (*outScant)->max_iter == 0) {
        return SECANT_ERR_INVALID;
    }
    double (*f)(double x) = (*outScant)->func.f;
    double *x0 = &(*outScant)->x0;
    double *x1 = &(*outScant)->x1;
    double *x2 = (double *)malloc(sizeof(double));
    long double dy2dx = 0.0;

    for (size_t i = 0; i < (*outScant)->max_iter; i++) {

        Secant_Err err = secant_slope(f, *x0, *x1, &dy2dx);
        if (err == SECANT_ERR_DIVIDE_BY_ZERO) {
            free(x2);
            return err;
        }
        *x2 = *x1 - f(*x1) / dy2dx;
        if (!*x2) {
            free(x2);
            return SECANT_ERR_INVALID;
        }
        if (fabsl(*x2 - *x1) < (*outScant)->tol || fabs(f(*x2)) < (*outScant)->tol) {
            if (fabs(*x2) < 2e-8) { // 根接近于0
                *outRoot = 0;
                free(x2);
                return SECANT_OK;
            }
            *outRoot = *x2;
            free(x2);
            return SECANT_OK;
        }
        *x0 = *x1;
        *x1 = *x2;
    }
    free(x2);
    return SECANT_ERR_MAXITER;
}

// 构造函数
Secant_Err NonLinearScant_create(double (*f)(double x), const double x0, const double x1,
    const double tol, const size_t max_iter,
    NonLinearScant *outScant,
    const char *name) {
    if (!f || tol <= 1e-65 || max_iter == 0 || !outScant || !name) {
        return SECANT_ERR_INVALID;
    }
    NonLinearScant scant = (NonLinearScant)malloc(sizeof(struct NonLinearScant));
    if (!scant) {
        return SECANT_ERR_NOMEM;
    }
    scant->func.f = f;
    scant->x0 = x0;
    scant->x1 = x1;
    scant->tol = tol;
    scant->max_iter = max_iter;
    scant->name = (char *)name;
    *outScant = scant;
    return SECANT_OK;
}


// 析构函数
Secant_Err NonLinearScant_destroy(const NonLinearScant *inScant) {
    if (!inScant || !*inScant) {
        return SECANT_ERR_INVALID;
    }
    free(*inScant);
    return SECANT_OK;
}


// 公共API
const SecantAPI Secant = {
    .secant_solve = secant_solve,
    .NonLinearScant_create = NonLinearScant_create,
    .NonLinearScant_destroy = NonLinearScant_destroy
};