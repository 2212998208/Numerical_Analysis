#include "Trapezoidal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>


struct IntegrationApproximation {
    double (*f)(double x);
    double a;
    double b;
    size_t max_iter;
    char *name;
};



// 梯形法近似积分
Trapezoidal_Err trapezoidal_integration(const Trapezoidal *inTrap, double *outApproxIntegral) {
    if (!inTrap || !*inTrap || !outApproxIntegral) {
        return TRAP_ERR_INVAL;
    }

    double (*f)(double x) = (*inTrap)->f;
    double a = (*inTrap)->a;
    double b = (*inTrap)->b;
    size_t max_iter = (*inTrap)->max_iter;

    if (a >= b || max_iter == 0) {
        return TRAP_ERR_INVAL;
    }

    // 初始步长
    double delta = b - a;
    double h = delta / (int)max_iter;
    double integral = 0.0;

    // 迭代计算积分值
    for (size_t i = 0; i < max_iter; i++) {

        // 计算当前步长下的积分值
        integral += 0.5 * (f(a + ((int)i) * h) + f(a + ((int)(i + 1)) * h)) * h;
    }

    *outApproxIntegral = integral;

    return TRAP_OK;
}


// 构造函数
Trapezoidal_Err trapezoidal_create(double (*f)(double x), double a, double b, size_t max_iter, Trapezoidal *outTrap, const char *name) {
    if (!f || !outTrap) {
        return TRAP_ERR_INVAL;
    }
    if (a >= b || max_iter == 0) {
        return TRAP_ERR_INVAL;
    }

    const Trapezoidal trap = (Trapezoidal)malloc(sizeof(struct IntegrationApproximation));
    if (!trap) {
        return TRAP_ERR_NOMEM;
    }

    trap->f = f;
    trap->a = a;
    trap->b = b;
    trap->max_iter = max_iter;

    if (name) {
        size_t name_len = strlen(name);
        trap->name = (char *)malloc(name_len + 1);
        if (!trap->name) {
            free(trap);
            return TRAP_ERR_NOMEM;
        }
        strcpy(trap->name, name);
    } else {
        trap->name = NULL;
    }

    *outTrap = trap;
    return TRAP_OK;
}


// 析构函数
Trapezoidal_Err trapezoidal_destroy(Trapezoidal *inTrap) {
    if (inTrap && *inTrap) {
        if ((*inTrap)->name) {
            free((*inTrap)->name);
        }
        free(*inTrap);
        *inTrap = NULL;
        return TRAP_OK;
    }

    return TRAP_ERR_INVAL;
}

const TrapezoidalIntegrationAPI TI = {
    .trapezoidal_integration = trapezoidal_integration,
    .trapezoidal_create = trapezoidal_create,
    .trapezoidal_destroy = trapezoidal_destroy
};