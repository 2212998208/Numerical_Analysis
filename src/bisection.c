#include "bisection.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

struct BisectionFunction {
    double (*f)(double x);
};

struct BisectionRange {
    struct BisectionFunction func;
    double a;
    double b;
    double tol;
    size_t max_iter;
    char *name;
};


// 区间迭代次数计算
static Bisection_Err bisection_compute_iterations(const BisectionRange *outRange) {
    if (!*outRange || (*outRange)->tol <= 1e-15) {
        return BISECTION_ERR_INVALID;
    }

    const double a = (*outRange)->a;
    const double b = (*outRange)->b;
    const double tol = (*outRange)->tol;
    size_t max_iter = 0;

    for (double power = 1;fabs(b-a) / power > tol; max_iter++) {
        power *= 2;
    }

    (*outRange)->max_iter = max_iter;

    return BISECTION_OK;
}


// 介值定理
static Bisection_Err bisection_compute_roots(const BisectionRange *outRange) {
    if (!*outRange || (*outRange)->tol <= 1e-15 || (*outRange)->max_iter == 0) {
        return BISECTION_ERR_INVALID;
    }

    double *a = &((*outRange)->a);
    double *b = &((*outRange)->b);

    double (*f)(double x) = (*outRange)->func.f;

    if (f(*a) * f(*b) > 0) {
        return BISECTION_ERR_INVALID;
    }
    if (f(*a) * f(*b) < 0) {
    for (size_t i = 0; i < (*outRange)->max_iter; i++) {
        double M = (*a + (*b)) / 2;
        if (f(*a) * f(M) < 0) {
            *b = M;
        } else if (f(*b) * f(M) < 0) {
            *a = M;
        }
    }
    }
    return BISECTION_OK;
}


// 非线性方程数值解
Bisection_Err bisection_solve(BisectionRange *outRange) {
    if (!*outRange || (*outRange)->tol <= 1e-15) {
        return BISECTION_ERR_INVALID;
    }

    Bisection_Err err;

    err = bisection_compute_iterations(outRange);
    if (err != BISECTION_OK) {
        return err;
    }

    err = bisection_compute_roots(outRange);
    if (err != BISECTION_OK) {
        return err;
    }

    // 数值解
    printf("[左区间]a=%.15f [右区间]b=%.15f [迭代次数]max_iter=%zu [精度]epsilon=%e [非线性方程]f(x)=%s \n", (*outRange)->a,
        (*outRange)->b,
        (*outRange)->max_iter,
        (*outRange)->tol,
        (*outRange)->name);

    return BISECTION_OK;
}


// 构造函数
Bisection_Err BisectionRange_create(double (*f)(double x), double a,
    double b, double tol, BisectionRange *outRange, const char *name) {
    if (!outRange || !f || tol <= 1e-15 || a >= b) {
        return BISECTION_ERR_INVALID;
    }
    BisectionRange range = (BisectionRange)malloc(sizeof(struct BisectionRange));
    if (!range) {
        return BISECTION_ERR_NOMEM;
    }
    range->func.f = f;
    range->a = a;
    range->b = b;
    range->tol = tol;
    range->max_iter = 0;
    range->name = (char *)name;
    *outRange = range;
    return BISECTION_OK;

}

// 析构函数
Bisection_Err BisectionRange_destroy(BisectionRange *outRange) {
    if (!outRange || !*outRange) {
        return BISECTION_ERR_INVALID;
    }
    free(*outRange);
    *outRange = NULL;
    return BISECTION_OK;
}

// 获取区间中点
double BisectionRange_get_midpoint(const BisectionRange *outRange) {
    if (!*outRange) {
        return NAN;
    }
    return ((*outRange)->a + (*outRange)->b) / 2;
}


// 公共API
const BisectionAPI Bisection = {
    .bisection_solve = bisection_solve,
    .bisection_create = BisectionRange_create,
    .bisection_destroy = BisectionRange_destroy,
    .bisection_get_midpoint = BisectionRange_get_midpoint
};