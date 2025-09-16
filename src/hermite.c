

#include "hermite.h"

#include <math.h>
#include <stdlib.h>

// 函数前向声明
static Hermite_Err hermite_create_interpolator_impl(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator);
static Hermite_Err hermite_destroy_interpolator_impl(HermiteInterpolator *inInterpolator);
static HermiteDataset create_hermite_dataset_impl(size_t size, const double *x, const double *y, const double *dy);
static Hermite_Err destroy_hermite_dataset_impl(HermiteDataset *inDataset);
static Hermite_Err hermite_evaluate_impl(const HermiteInterpolator *inInterpolator, double x, double *outY, double *outDy);
static Hermite_Err hermite_compute_coefficients(HermiteInterpolator *interpolator, const HermiteDataset inDataset);
static struct HermitePoint make_hermite_point(double x, double y, double dy);

/** @brief Hermite 插值数据点，包含函数值和一阶导数值。*/
struct HermitePoint {
    double x;  /**< x 坐标 */
    double y;  /**< y 坐标 (函数值 f(x)) */
    double dy; /**< 一阶导数值 f'(x) */
};

/** @brief Hermite 数据集，包含一组 HermitePoint。*/
struct HermiteDataset {
    struct HermitePoint *points; /**< 指向数据点数组的指针 */
    size_t size;                 /**< 数据点的数量 (n) */
};

/** @brief Hermite 插值器，存储预计算的系数。*/
struct HermiteInterpolator {
    double *coefficients; /**< 差商表（或多项式系数）*/
    double *z_nodes;      /**< 重复节点序列 z_0, z_1, ..., z_{2n-1} */
    size_t size;          /**< 插值器的总阶数 (2n) */
};

/**
 * @brief 创建 Hermite 插值器的实现。
 *
 * 这是核心的准备步骤，它构建重复节点序列，并调用 `hermite_compute_coefficients`
 * 来计算差商表，最后将所有内容封装到 HermiteInterpolator 结构中。
 */
Hermite_Err hermite_create_interpolator_impl(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator) {
    if (!(*inDataset) || (*inDataset)->size == 0 || !outInterpolator || *outInterpolator != NULL) {
        return HERMITE_ERR_INVALID;
    }

    const size_t n = (*inDataset)->size;
    const size_t N = 2 * n; // 总节点数

    // 分配插值器结构
    HermiteInterpolator interpolator = (HermiteInterpolator)calloc(1, sizeof(struct HermiteInterpolator));
    if (!interpolator) return HERMITE_ERR_NOMEM;

    interpolator->z_nodes = (double *)malloc(N * sizeof(double));
    interpolator->coefficients = (double *)malloc(N * sizeof(double));
    if (!interpolator->z_nodes || !interpolator->coefficients) {
        free(interpolator->z_nodes);
        free(interpolator->coefficients);
        free(interpolator);
        return HERMITE_ERR_NOMEM;
    }
    interpolator->size = N;

    // 构造重复节点序列 z 和零阶差商
    for (size_t i = 0; i < n; i++) {
        interpolator->z_nodes[2 * i] = (*inDataset)->points[i].x;
        interpolator->z_nodes[2 * i + 1] = (*inDataset)->points[i].x;
        interpolator->coefficients[2 * i] = (*inDataset)->points[i].y;
        interpolator->coefficients[2 * i + 1] = (*inDataset)->points[i].y;
    }

    // 计算完整的差商表
    Hermite_Err err = hermite_compute_coefficients(&interpolator, *inDataset);
    if (err != HERMITE_OK) {
        hermite_destroy_interpolator_impl(&interpolator); // 清理
        return err;
    }

    *outInterpolator = interpolator;
    return HERMITE_OK;
}

/**
 * @brief 销毁 Hermite 插值器并释放内存的实现。
 */
Hermite_Err hermite_destroy_interpolator_impl(HermiteInterpolator *inInterpolator) {
    if (inInterpolator && *inInterpolator) {
        free((*inInterpolator)->coefficients);
        free((*inInterpolator)->z_nodes);
        free(*inInterpolator);
        *inInterpolator = NULL;
    }
    return HERMITE_OK;
}

/**
 * @brief 计算 Hermite 插值的差商表。
 *
 * 此函数是 Hermite 插值算法的核心。它建立在一个扩展的节点序列上，
 * 并对差商的定义进行了扩展：
 * - f[z_i] = f(z_i)
 * - f[z_i, z_{i+1}] = f'(z_i)  如果 z_i = z_{i+1}
 * - f[z_i, ..., z_j] = (f[z_{i+1}, ..., z_j] - f[z_i, ..., z_{j-1}]) / (z_j - z_i)
 */
static Hermite_Err hermite_compute_coefficients(HermiteInterpolator *interpolator, const HermiteDataset *inDataset) {
    const size_t N = (*interpolator)->size;
    double *c = (*interpolator)->coefficients; // 使用短别名 c
    const double *z = (*interpolator)->z_nodes; // 使用短别名 z

    // 恢复初始值
    for (size_t i = 0; i < inDataset->size; i++) {
        c[2 * i] = inDataset->points[i].y;
        c[2 * i + 1] = inDataset->points[i].y;
    }

    // 计算一阶差商
    for (size_t i = N - 1; i > 0; i--) {
        if (i % 2 != 0) { // i 是奇数, e.g. i=1,3,5...
            c[i] = inDataset->points[i/2].dy;
        } else { // i 是偶数
            c[i] = (c[i] - c[i-1]) / (z[i] - z[i-1]);
        }
    }

    // 计算更高阶差商
    for (size_t j = 2; j < N; j++) {
        for (size_t i = N - 1; i >= j; i--) {
            double denominator = z[i] - z[i-j];
             if (fabs(denominator) < 1e-9) {
                return HERMITE_ERR_DIVBYZERO;
            }
            c[i] = (c[i] - c[i-1]) / denominator;
        }
    }

    return HERMITE_OK;
}

/**
 * @brief 创建 Hermite 数据集的实现。
 */
HermiteDataset create_hermite_dataset_impl(size_t size, const double *x, const double *y, const double *dy) {
    if (size == 0 || !x || !y || !dy) {
        return NULL;
    }
    HermiteDataset dataset = (HermiteDataset)malloc(sizeof(struct HermiteDataset));
    if (!dataset) return NULL;

    dataset->points = (struct HermitePoint *)malloc(size * sizeof(struct HermitePoint));
    if (!dataset->points) {
        free(dataset);
        return NULL;
    }
    for (size_t i = 0; i < size; i++) {
        dataset->points[i] = make_hermite_point(x[i], y[i], dy[i]);
    }
    dataset->size = size;
    return dataset;
}

/**
 * @brief 销毁 Hermite 数据集的实现。
 */
Hermite_Err destroy_hermite_dataset_impl(HermiteDataset *inDataset) {
    if (inDataset && *inDataset) {
        free((*inDataset)->points);
        free(*inDataset);
        *inDataset = NULL;
    }
    return HERMITE_OK;
}

/**
 * @brief 使用 Hermite 插值器在指定点求值（包括函数值和导数值）的实现。
 *
 * 采用扩展的秦九韶算法（Horner's method）同时计算 P(x) 和 P'(x)。
 */
Hermite_Err hermite_evaluate_impl(const HermiteInterpolator *inInterpolator, const double x, double *outY, double *outDy) {
    if (!(*inInterpolator) || (*inInterpolator)->size == 0) {
        return HERMITE_ERR_INVALID;
    }

    const size_t N = (*inInterpolator)->size;
    const double *c = (*inInterpolator)->coefficients;
    const double *z = (*inInterpolator)->z_nodes;

    // 从最高阶系数开始
    double result = c[N - 1];
    double derivative = 0.0; // 导数初值

    // 嵌套计算
    for (int i = (int)N - 2; i >= 0; i--) {
        derivative = derivative * (x - z[i]) + result;
        result = result * (x - z[i]) + c[i];
    }

    if (outY) *outY = result;
    if (outDy) *outDy = derivative;

    return HERMITE_OK;
}

/**
 * @brief 创建一个 HermitePoint 对象的便捷函数。
 */
static struct HermitePoint make_hermite_point(const double x, const double y, const double dy) {
    struct HermitePoint p = { x, y, dy };
    return p;
}


// 全局API实例
const HermiteAPI Hermite = {
    .hermite_create_interpolator = hermite_create_interpolator_impl,
    .hermite_destroy_interpolator = hermite_destroy_interpolator_impl,
    .create_hermite_dataset = create_hermite_dataset_impl,
    .destroy_hermite_dataset = destroy_hermite_dataset_impl,
    .hermite_evaluate = hermite_evaluate_impl,
};