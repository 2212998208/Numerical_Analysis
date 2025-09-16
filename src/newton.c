#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "newton.h"
#include "lagrange.h"


// 函数前向声明
static Newton_Err newton_interpolate_impl(NewtonDataSet **inNewtonDataSet, DataSet **inDataset, double x, double *outY);
static Newton_Err create_newton_dataset_impl(NewtonDataSet **outDataset, size_t size);
static Newton_Err destroy_dataset_impl(NewtonDataSet **outDataset);
static Newton_Err print_newton_dataset_impl(const NewtonDataSet **outDataset);
static Newton_Err compute_divided_differences(NewtonDataSet *dataset, const Point *points, size_t size);
static Newton_Err newton_evaluate(const NewtonDataSet *dataset, double x, double *outY);

/**
 * @brief 牛顿插值所需的数据集结构。
 *
 * 存储差商表和对应的 x 节点。
 */
struct NewtonDataSet {
    double *divided_differences; /**< 差商表，存储 f[x_0], f[x_0,x_1], ... */
    double *x_nodes;             /**< 插值节点的 x 坐标 */
    size_t x_nodes_size;         /**< 节点的数量 */
};

/**
 * @brief 牛顿插值法多项式计算的实现。
 *
 * 此函数是一个高级封装，它首先根据输入的数据点计算差商表（如果尚未计算），
 * 然后使用该表来评估在点 x 处的插值。
 */
Newton_Err newton_interpolate_impl(NewtonDataSet **inNewtonDataset, DataSet **inDataset, double x, double *outY) {
    if (!inNewtonDataset || !*inNewtonDataset || !inDataset || !*inDataset || !outY) {
        return NEWTON_ERR_INVALID;
    }

    size_t size = (*inNewtonDataset)->x_nodes_size;
    Point *points = Lagrange.get_points(inDataset);

    if (size == 0 || !points) {
        return NEWTON_ERR_INVALID;
    }

    if (size == 1) {
        *outY = points[0].y; // 只有一个点，直接返回y值
        return NEWTON_OK;
    }

    // 如果差商表尚未计算，则进行计算
    if ((*inNewtonDataset)->divided_differences == NULL) {
        Newton_Err err = compute_divided_differences(*inNewtonDataset, points, size);
        if (err != NEWTON_OK) {
            return err;
        }
    }

    // 使用计算好的差商表进行求值
    return newton_evaluate(*inNewtonDataset, x, outY);
}

/**
 * @brief 创建一个空的牛顿数据集的实现。
 *
 * 仅分配结构本身和设置大小，差商表和 x 节点将在 `compute_divided_differences` 中分配。
 */
Newton_Err create_newton_dataset_impl(NewtonDataSet **outDataset, const size_t size) {
    if (!outDataset || *outDataset != NULL) { // 不应覆盖已存在的指针
        return NEWTON_ERR_INVALID;
    }
    NewtonDataSet *dataset = (NewtonDataSet *)calloc(1, sizeof(NewtonDataSet));
    if (!dataset) {
        return NEWTON_ERR_NOMEM;
    }
    dataset->x_nodes_size = size;
    *outDataset = dataset;
    return NEWTON_OK;
}

/**
 * @brief 销毁牛顿数据集并释放所有相关内存的实现。
 */
Newton_Err destroy_dataset_impl(NewtonDataSet **outDataset) {
    if (outDataset && *outDataset) {
        free((*outDataset)->divided_differences);
        free((*outDataset)->x_nodes);
        free(*outDataset);
        *outDataset = NULL;
    }
    return NEWTON_OK;
}

/**
 * @brief 使用预先计算的差商表来评估牛顿插值多项式。
 *
 * 采用秦九韶算法（Horner's method）进行高效计算。
 * P(x) = f[x_0] + f[x_0,x_1](x-x_0) + ...
 *      = f[x_0] + (x-x_0)(f[x_0,x_1] + (x-x_1)(...))
 *
 * @param dataset 包含差商表和 x 节点的数据集。
 * @param x       要求值的点。
 * @param outY    用于存储结果的指针。
 * @return        操作状态。
 */
static Newton_Err newton_evaluate(const NewtonDataSet *dataset, double x, double *outY) {
    if (!dataset || !dataset->divided_differences || !dataset->x_nodes) {
        return NEWTON_ERR_INVALID;
    }

    const size_t size = dataset->x_nodes_size;
    // 从最高阶的差商开始
    double result = dataset->divided_differences[size - 1];

    // 嵌套乘法
    for (int i = (int)size - 2; i >= 0; i--) {
        result = result * (x - dataset->x_nodes[i]) + dataset->divided_differences[i];
    }

    *outY = result;
    return NEWTON_OK;
}

/**
 * @brief 打印牛顿数据集（差商表和节点）以供调试的实现。
 */
Newton_Err print_newton_dataset_impl(const NewtonDataSet **outDataset) {
    if (!outDataset || !*outDataset) {
        printf("NewtonDataSet pointer is NULL.\n");
        return NEWTON_ERR_INVALID;
    }
    const NewtonDataSet* dataset = *outDataset;
    if (!dataset->divided_differences || !dataset->x_nodes) {
        printf("NewtonDataSet is not fully initialized (missing data arrays).\n");
        return NEWTON_ERR_INVALID;
    }

    printf("Newton Dataset (size: %zu):\n", dataset->x_nodes_size);
    printf("  X Nodes:           ");
    for (size_t i = 0; i < dataset->x_nodes_size; i++) {
        printf("%-8.3f ", dataset->x_nodes[i]);
    }
    printf("\n");
    printf("  Divided Differences: ");
    for (size_t i = 0; i < dataset->x_nodes_size; i++) {
        printf("%-8.3f ", dataset->divided_differences[i]);
    }
    printf("\n");

    return NEWTON_OK;
}

/**
 * @brief 计算差商表。
 *
 * 这是牛顿插值的核心计算步骤。它填充 NewtonDataSet 中的
 * `divided_differences` 和 `x_nodes` 数组。
 *
 * @param dataset 指向要填充的牛顿数据集。
 * @param points  原始数据点。
 * @param size    数据点的数量。
 * @return        操作状态。
 */
static Newton_Err compute_divided_differences(NewtonDataSet *dataset, const Point *points, const size_t size) {
    if (!dataset || !points || size == 0) {
        return NEWTON_ERR_INVALID;
    }

    // 分配内存
    dataset->divided_differences = (double *)malloc(size * sizeof(double));
    dataset->x_nodes = (double *)malloc(size * sizeof(double));
    if (!dataset->divided_differences || !dataset->x_nodes) {
        free(dataset->divided_differences);
        free(dataset->x_nodes);
        dataset->divided_differences = NULL;
        dataset->x_nodes = NULL;
        return NEWTON_ERR_NOMEM;
    }
    dataset->x_nodes_size = size;

    // 初始化零阶差商（即 y 值）和 x 节点
    for (size_t i = 0; i < size; i++) {
        dataset->divided_differences[i] = points[i].y;
        dataset->x_nodes[i] = points[i].x;
    }

    // 逐阶计算差商
    // j 是阶数
    for (size_t j = 1; j < size; j++) {
        // i 是当前计算的差商在数组中的索引
        for (size_t i = size - 1; i >= j; i--) {
            double numerator = dataset->divided_differences[i] - dataset->divided_differences[i - 1];
            double denominator = dataset->x_nodes[i] - dataset->x_nodes[i - j];
            if (fabs(denominator) < 1e-9) {
                // 遇到重复的x节点，无法计算差商
                return NEWTON_ERR_DIVBYZERO;
            }
            dataset->divided_differences[i] = numerator / denominator;
        }
    }

    return NEWTON_OK;
}


// 公共API
const NewtonAPI Newton = {
    .newton_interpolate = newton_interpolate_impl,
    .create_newton_dataset = create_newton_dataset_impl,
    .destroy_dataset = destroy_dataset_impl,
    .print_newton_dataset = print_newton_dataset_impl
};