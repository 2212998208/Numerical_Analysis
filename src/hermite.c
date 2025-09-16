

#include "hermite.h"

#include <math.h>
#include <stdlib.h>

static Hermite_Err hermite_compute_coefficients(HermiteInterpolator *interpolator, HermiteDataset *inDataset);
static struct HermitePoint hermite_point(double x, double y, double dy);

struct HermitePoint {
    double x;
    double y;
    double dy;
};

struct HermiteDataset {
    struct HermitePoint *points;
    size_t size;
};

struct HermiteInterpolator {
    double *coefficients;
    double *z_node;
    size_t size;
};


// 创建插值器
Hermite_Err hermite_create_interpolator(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator) {
    if (!*inDataset || (*inDataset)->size == 0) {
        return HERMITE_ERR_INVALID;
    }

    const size_t n = (*inDataset)->size;
    const size_t N = 2 * n;

    double *z_nodes = malloc(N * sizeof(double));
    if (!z_nodes) {
        return HERMITE_ERR_NOMEM;
    }

    double *coefficients = malloc(N * sizeof(double));
    if (!coefficients) {
        free(z_nodes);
        return HERMITE_ERR_NOMEM;
    }

    for (size_t i = 0; i < n; i++) {
        z_nodes[2 * i] = (*inDataset)->points[i].x;
        z_nodes[2 * i + 1] = (*inDataset)->points[i].x;
        coefficients[2 * i] = (*inDataset)->points[i].y;
        coefficients[2 * i + 1] = (*inDataset)->points[i].y;
    }

    HermiteInterpolator interpolator = malloc(sizeof(struct HermiteInterpolator));
    if (!interpolator) {
        free(z_nodes);
        free(coefficients);
        return HERMITE_ERR_NOMEM;
    }

    interpolator->coefficients = coefficients;
    interpolator->z_node = z_nodes;
    interpolator->size = N;

    Hermite_Err err = hermite_compute_coefficients(&interpolator, inDataset);
    if (err != HERMITE_OK) {
        free(z_nodes);
        free(coefficients);
        free(interpolator);
        return err;
    }

    *outInterpolator = interpolator;
    return HERMITE_OK;
}


// 销毁插值器
Hermite_Err hermite_destroy_interpolator(HermiteInterpolator *InInterpolator) {
    if (!*InInterpolator) {
        return HERMITE_ERR_INVALID;
    }
    free((*InInterpolator)->coefficients);
    free((*InInterpolator)->z_node);
    free(*InInterpolator);
    *InInterpolator = NULL;
    return HERMITE_OK;
}






// 均差系数计算
static Hermite_Err hermite_compute_coefficients(HermiteInterpolator *interpolator, HermiteDataset *inDataset) {
    if (!*interpolator || !*inDataset || (*inDataset)->size == 0 ||
        (*interpolator)->size != 2 * (*inDataset)->size ||
        (*interpolator)->size == 0) {
        return HERMITE_ERR_INVALID;
    }

    const size_t loop_size = (*interpolator)->size;
    double *coefficients = (*interpolator)->coefficients;
    const double *z_nodes = (*interpolator)->z_node;

    for (size_t i = 1; i < loop_size; i++) {
        for (size_t j = loop_size - 1; j >= i; j--) {
            double numerator = coefficients[j] - coefficients[j - 1];
            double denominator = z_nodes[j] - z_nodes[j - i];
            if (fabs(denominator) < 1e-9 || fabs(denominator) == 0) {
                // 处理重根的情况
                size_t original_index = j / 2;
                coefficients[j] = (*inDataset)->points[original_index].dy;
            } else {
                coefficients[j] = numerator / denominator;
            }

        }
    }

    return HERMITE_OK;
}






// 数据集的创建
HermiteDataset create_hermite_dataset(size_t size, const double *x, const double *y, const double *dy) {
    if (size == 0) {
        return NULL;
    }
    HermiteDataset dataset = malloc(sizeof(struct HermiteDataset));
    if (!dataset) {
        return NULL;
    }
    dataset->points = (struct HermitePoint *)malloc(size * sizeof(struct HermitePoint));
    if (!dataset->points) {
        free(dataset);
        return NULL;
    }
    for (size_t i = 0; i < size; i++) {
        dataset->points[i] = hermite_point(x[i], y[i], dy[i]);
    }
    dataset->size = size;

    return dataset;
}






// 数据集的销毁
Hermite_Err destroy_hermite_dataset(HermiteDataset *inDataset) {
    if (!*inDataset) {
        return HERMITE_ERR_INVALID;
    }
    free((*inDataset)->points);
    free(*inDataset);
    *inDataset = NULL;
    return HERMITE_OK;
}


// 求值函数
Hermite_Err hermite_evaluate(const HermiteInterpolator *InInterpolator, const double x, double *outY, double *outDy) {
    if (!*InInterpolator || (*InInterpolator)->size == 0) {
        return HERMITE_ERR_INVALID;
    }

    const size_t N = (*InInterpolator)->size;
    const double *coefficients = (*InInterpolator)->coefficients;
    const double *z_nodes = (*InInterpolator)->z_node;

    double result = coefficients[N - 1];
    double derivative = 0.0;

    for (size_t i = N - 1; i > 0; i--) {
        derivative = derivative * (x - z_nodes[i - 1]) + result;
        result = result * (x - z_nodes[i - 1]) + coefficients[i - 1];
    }

    if (outY) {
        *outY = result;
    }
    if (outDy) {
        *outDy = derivative;
    }

    return HERMITE_OK;
}

static struct HermitePoint hermite_point(const double x, const double y, const double dy) {
    struct HermitePoint p = { x, y, dy };
    return p;
}


// 全局API实例
const HermiteAPI Hermite = {
    hermite_create_interpolator,
    hermite_destroy_interpolator,
    create_hermite_dataset,
    destroy_hermite_dataset,
    hermite_evaluate,
};