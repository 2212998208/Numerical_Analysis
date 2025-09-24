#include "RK4.h"

#include <stdlib.h>


struct RK4Data {
    double (*dx2dt) (double x, double t);
    double t0;
    double x0;
    double h;
    size_t max_iter;
};


// 构造函数
RK4_Err rk4_create(rk *outRK, double (*f)(double x, double t),
    const double t0, const double x0, const double h, size_t max_iter) {
    if (!outRK || !f || h <= 0|| max_iter == 0) {
        return RK4_INVALID;
    }

    const rk newRK = malloc(sizeof(struct RK4Data));
    if (!newRK) {
        return RK4_POINTER_ERROR;
    }

    newRK->dx2dt = f;
    newRK->t0 = t0;
    newRK->x0 = x0;
    newRK->h = h;
    newRK->max_iter = max_iter;

    *outRK = newRK;

    return RK4_OK;

}


// 析构函数
RK4_Err rk4_destroy(rk *outRK) {
    if (!outRK || !*outRK) {
        return RK4_POINTER_ERROR;
    }

    free(*outRK);
    *outRK = NULL;

    return RK4_OK;
}


// 四阶龙格-库塔法求解常微分方程初值问题
RK4_Err rk4_solve(const rk *inRk, double *outY) {
    if (!*inRk || !outY) {
        return RK4_POINTER_ERROR;
    }


    double t0 = (*inRk)->t0;
    double x0 = (*inRk)->x0;
    double h = (*inRk)->h;
    size_t max_iter = (*inRk)->max_iter;
    double (*f)(double x, double t) = (*inRk)->dx2dt;

    // 计算核心
    for (size_t i = 0; i < max_iter; ++i) {
        double k1 = h * f(x0, t0);
        double k2 = h * f(x0 + 0.5 * k1, t0 + 0.5 * h);
        double k3 = h * f(x0 + 0.5 * k2, t0 + 0.5 * h);
        double k4 = h * f(x0 + k3, t0 + h);

        x0 += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
        t0 += h;
    }

    *outY = x0;

    return RK4_OK;
}


// API设计
const RK4_API RK4 = {
    .create = rk4_create,
    .destroy = rk4_destroy,
    .solve = rk4_solve
};