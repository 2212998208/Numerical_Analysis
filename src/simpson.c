#include "simpson.h"

#include <stdlib.h>
#include <string.h>


struct IntegrationApproximation {
    double (*f)(double x);
    double a;
    double b;
    size_t max_iter;
    char *name;
};


// 辛普森法近似积分
Simpson_Err simpson_integration(const Simpson *inSimpson, double *outApproxIntegral) {
    if (!inSimpson || !*inSimpson || !outApproxIntegral) {
        return SIMPSON_ERR_INVAL;
    }

    double (*f)(double x) = (*inSimpson)->f;
    const double a = (*inSimpson)->a;
    const double b = (*inSimpson)->b;
    const size_t max_iter = (*inSimpson)->max_iter;

    if (a >= b || max_iter == 0 || max_iter % 2 != 0) {
        return SIMPSON_ERR_INVAL;
    }

    // 初始步长
    const double delta = b - a;
    const double h = delta / (int)max_iter;
    double integral = (f(a) + f(b)) * h / 3.0;

    // 迭代计算积分值
    for (int i = 1; i < max_iter; i++) {
        double x = a + i * h;
        // 计算当前步长下的积分值
        integral += (i % 2 == 0) ? 2.0 * f(x) * h / 3.0 : 4.0 * f(x) * h / 3.0;
    }

    *outApproxIntegral = integral;

    return SIMPSON_OK;
}




// 构造函数
Simpson_Err simpson_create(double (*f)(double x), double a, double b, size_t max_iter, Simpson *outSimpson, const char *name) {
    if (!f || !outSimpson) {
        return SIMPSON_ERR_INVAL;
    }
    if (a >= b || max_iter == 0 || max_iter % 2 != 0) {
        return SIMPSON_ERR_INVAL;
    }

    Simpson simpson = (Simpson)malloc(sizeof(struct IntegrationApproximation));
    if (!simpson) {
        return SIMPSON_ERR_NOMEM;
    }

    simpson->f = f;
    simpson->a = a;
    simpson->b = b;
    simpson->max_iter = max_iter;

    if (name) {
        const size_t name_len = strlen(name);
        simpson->name = (char *)malloc(name_len + 1);
        if (!simpson->name) {
            free(simpson);
            return SIMPSON_ERR_NOMEM;
        }
        strcpy(simpson->name, name);
    } else {
        simpson->name = NULL;
    }

    *outSimpson = simpson;
    return SIMPSON_OK;
}




// 析构函数
Simpson_Err simpson_destroy(Simpson *inSimpson) {
    if (!inSimpson || !*inSimpson) {
        return SIMPSON_ERR_INVAL;
    }

    if ((*inSimpson)->name) {
        free((*inSimpson)->name);
    }
    free(*inSimpson);
    *inSimpson = NULL;

    return SIMPSON_OK;
}

const SimpsonAPI SI = {
    .simpson_create = simpson_create,
    .simpson_integration = simpson_integration,
    .simpson_destroy = simpson_destroy,
};