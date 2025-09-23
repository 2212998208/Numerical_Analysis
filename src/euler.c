#include "euler.h"

#include <stdlib.h>


struct EulerDifferential {
    double (*dx2dt)(double x, double t);
    double x0;
    double t0;
    double Δt;
};

// 构造函数
Euler_Err euler_create(Euler *outEuler, double (*dx2dt)(double x, double t), double x0, double t0, double Δt) {
    if (outEuler == NULL || dx2dt == NULL || Δt <= t0) {
        return EULER_INVALID;
    }
    Euler euler = (Euler)malloc(sizeof(struct EulerDifferential));
    if (euler == NULL) {
        return EULER_INVALID;
    }
    euler->dx2dt = dx2dt;
    euler->x0 = x0;
    euler->t0 = t0;
    euler->Δt = Δt;
    *outEuler = euler;
    return EULER_OK;
}


// 析构函数
Euler_Err euler_destroy(Euler *euler) {
    if (euler == NULL || *euler == NULL) {
        return EULER_INVALID;
    }
    free(*euler);
    *euler = NULL;
    return EULER_OK;
}

// 前向欧拉法解决OEDS的初值问题
Euler_Err euler_solve(const Euler *inEuler, size_t max_iter, double *xn) {
    if (inEuler == NULL || *inEuler == NULL || xn == NULL) {
        return EULER_INVALID;
    }

    double x0 = (*inEuler)->x0;
    double t0 = (*inEuler)->t0;

    double delta_t = (*inEuler)->Δt;
    double (*dx2dt)(double x, double t) = (*inEuler)->dx2dt;

    for (size_t i = 0; i < max_iter; i++) {
        *xn = x0 + delta_t * dx2dt(x0, t0);
        x0 = *xn;
        t0 += delta_t;
    }

    return EULER_OK;

}


// 公共API
const EulerAPI EULER  = {
    .euler_create = euler_create,
    .euler_destroy = euler_destroy,
    .euler_solve = euler_solve
};