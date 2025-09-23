#include "double_simpson.h"

#include <stdlib.h>
#include <string.h>


struct DoubleIntegrationApproximation {
    double (*f)(double x, double y);
    double x_a;
    double x_b;
    double y_c;
    double y_d;
    size_t n;
    size_t m;
    char *name;
};


// 权重系数
static int get_simpson_weight_1d(int index, int max_index) {
    if (index == 0 || index == max_index) {
        // 端点权重为 1
        return 1;
    }
    if (index % 2 == 1) {
        // 奇数点权重为 4
        return 4;
    }
    // 偶数内点权重为 2
    return 2;
}


// 二重辛普森法近似积分
Simpson_Err double_simpson_integration(const DoubleSimpson *inDoubleSimpson, double *outApproxIntegral) {
    if (!inDoubleSimpson || !*inDoubleSimpson || !outApproxIntegral) {
        return SIMPSON_ERR_INVAL;
    }

    double (*f)(double x, double y) = (*inDoubleSimpson)->f;
    const double x_a = (*inDoubleSimpson)->x_a;
    const double x_b = (*inDoubleSimpson)->x_b;
    const double y_c = (*inDoubleSimpson)->y_c;
    const double y_d = (*inDoubleSimpson)->y_d;
    const size_t n = (*inDoubleSimpson)->n;
    const size_t m = (*inDoubleSimpson)->m;

    if (x_a >= x_b || y_c >= y_d || n == 0 || m == 0 || n % 2 != 0 || m % 2 != 0) {
        return SIMPSON_ERR_INVAL;
    }

    // --- 步骤 2: 计算步长 ---
    double h = (x_b - x_a) / n;
    double k = (y_d - y_c) / m;
    double total_sum = 0.0;

    // --- 步骤 3: 遍历所有网格点，计算加权和 ---
    // 外层循环遍历 y 方向 (从 j=0 到 m)
    for (int j = 0; j <= m; ++j) {
        // 内层循环遍历 x 方向 (从 i=0 到 n)
        for (int i = 0; i <= n; ++i) {
            // 计算当前点的坐标
            double x_i = x_a + i * h;
            double y_j = y_c + j * k;

            // 获取 x 和 y 方向各自的权重
            int weight_x = get_simpson_weight_1d(i, n);
            int weight_y = get_simpson_weight_1d(j, m);

            // 二维权重是两个一维权重的乘积
            int weight_ij = weight_x * weight_y;

            // 将加权后的函数值累加到总和中
            total_sum += weight_ij * f(x_i, y_j);
        }
    }

    // --- 步骤 4: 应用最终公式 ---
    // 最终结果 = (h*k/9) * 加权总和
     *outApproxIntegral = (h * k / 9.0) * total_sum;

    return SIMPSON_OK;
}


// 构造函数
Simpson_Err double_simpson_create(double (*f)(double x, double y),
                                 const double x_a, const double x_b,
                                 const double y_c, const double y_d,
                                 const size_t n, const size_t m,
                                 DoubleSimpson *outDoubleSimpson,
                                 const char *name) {
    if (!f || !outDoubleSimpson) {
        return SIMPSON_ERR_INVAL;
    }
    if (x_a >= x_b || y_c >= y_d || n == 0 || m == 0 || n % 2 != 0 || m % 2 != 0) {
        return SIMPSON_ERR_INVAL;
    }

    DoubleSimpson doubleSimpson = (DoubleSimpson)malloc(sizeof(struct DoubleIntegrationApproximation));
    if (!doubleSimpson) {
        return SIMPSON_ERR_NOMEM;
    }

    doubleSimpson->f = f;
    doubleSimpson->x_a = x_a;
    doubleSimpson->x_b = x_b;
    doubleSimpson->y_c = y_c;
    doubleSimpson->y_d = y_d;
    doubleSimpson->n = n;
    doubleSimpson->m = m;

    if (name) {
        doubleSimpson->name = (char *)malloc(strlen(name) + 1);
        if (!doubleSimpson->name) {
            free(doubleSimpson);
            return SIMPSON_ERR_NOMEM;
        }
        strcpy(doubleSimpson->name, name);
    } else {
        doubleSimpson->name = NULL;
    }

    *outDoubleSimpson = doubleSimpson;
    return SIMPSON_OK;
}


// 析构函数
Simpson_Err double_simpson_destroy(DoubleSimpson *inDoubleSimpson) {
    if (!inDoubleSimpson || !*inDoubleSimpson) {
        return SIMPSON_ERR_INVAL;
    }

    if ((*inDoubleSimpson)->name) {
        free((*inDoubleSimpson)->name);
    }
    free(*inDoubleSimpson);
    *inDoubleSimpson = NULL;

    return SIMPSON_OK;
}

const DoubleSimpsonAPI DS = {
    .double_simpson_create = double_simpson_create,
    .double_simpson_destroy = double_simpson_destroy,
    .double_simpson_integrate = double_simpson_integration
};

