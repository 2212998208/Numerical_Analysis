#include "lagrange.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


// 函数前向声明
static double calculate_basis_polynomial(const DataSet *dataset, int k, double x);
static Lagrange_Err create_dataset_impl(DataSet **outDataset, Point *points, size_t size);
static void destroy_dataset_impl(DataSet *dataset);
static double lagrange_interpolate_impl(DataSet *dataset, double x);
static Lagrange_Err create_dataset_from_points_impl(DataSet **outDataset, size_t size);
static DataSet* empty_dataset_impl(void);
static Point* get_points_impl(DataSet **outDataset);


/**
 * @brief 拉格朗日插值所需的数据集结构。
 *
 * 内部定义，对用户透明。
 */
struct DataSet {
    Point *points; /**< 指向数据点数组的指针 */
    size_t size;   /**< 数据点的数量 */
};

/**
 * @brief 拉格朗日插值法多项式计算的实现。
 *
 * P(x) = sum_{i=0 to n-1} (y_i * L_i(x))
 * 其中 L_i(x) 是拉格朗日基函数。
 */
double lagrange_interpolate_impl(DataSet *dataset, const double x) {
    if (!dataset || dataset->size == 0 || !dataset->points) {
        // 返回一个特定的错误值，或者可以设置一个错误状态
        return NAN; // 使用 NAN 表示无效操作
    }

    if (dataset->size == 1) {
        return (dataset->points)[0].y; // 只有一个点，插值结果就是该点的y值
    }

    double total_sum = 0.0;
    for (size_t i = 0; i < dataset->size; ++i) {
        double y_i = (dataset->points)[i].y;

        // 计算第 i 个基函数 L_i(x)
        double l_i = calculate_basis_polynomial(dataset, (int)i, x);
        if (isnan(l_i)) {
            // 如果基函数计算中出现除零（即 x 与某个数据点 x_j 重合），
            // 并且 x_j 不是当前基函数对应的 x_i，这是一个错误。
            // 但如果 x == x_i，则插值结果应为 y_i。
            // 这里为了简化，我们检查 x 是否与任何一个数据点 x_j 接近。
            for (size_t j = 0; j < dataset->size; ++j) {
                if (fabs(x - dataset->points[j].x) < 1e-9) {
                    return dataset->points[j].y;
                }
            }
            return NAN; // 如果不是精确匹配，则返回错误
        }
        total_sum += y_i * l_i;
    }

    return total_sum;
}

/**
 * @brief 计算拉格朗日基函数 L_k(x) 的值。
 *
 * L_k(x) = product_{j=0 to n-1, j!=k} ( (x - x_j) / (x_k - x_j) )
 */
static double calculate_basis_polynomial(const DataSet *dataset, const int k, const double x) {
    double product = 1.0;
    double x_k = (dataset->points)[k].x;

    for (size_t i = 0; i < dataset->size; ++i) {
        if ((int)i == k) continue;

        double x_i = (dataset->points)[i].x;
        double denominator = x_k - x_i;

        // 检查分母是否为零，这表示数据点中有重复的 x 值，这是无效的。
        if (fabs(denominator) < 1e-9) {
            return NAN; // 返回 NAN 表示计算失败
        }
        double numerator = x - x_i;
        product *= (numerator / denominator);
    }

    return product;
}

/**
 * @brief 从一个 Point 数组创建并初始化数据集的实现。
 *
 * 函数会为数据集和其内部的 points 数组分配新的内存，并复制数据。
 */
Lagrange_Err create_dataset_impl(DataSet **outDataset, Point *points, size_t size) {
    if (!outDataset || !points || size == 0) {
        return LAGRANGE_ERR_INVALID;
    }
    // 分配 DataSet 结构本身的内存
    DataSet *dataset = (DataSet *)malloc(sizeof(DataSet));
    if (!dataset) {
        return LAGRANGE_ERR_NOMEM;
    }
    // 分配存储 Point 的数组内存
    dataset->points = (Point *)malloc(size * sizeof(Point));
    if (!dataset->points) {
        free(dataset); // 分配失败，清理已分配的内存
        return LAGRANGE_ERR_NOMEM;
    }
    // 复制数据
    for (size_t i = 0; i < size; ++i) {
        (dataset->points)[i] = points[i];
    }
    dataset->size = size;
    *outDataset = dataset;
    return LAGRANGE_OK;
}

/**
 * @brief 销毁一个数据集并释放其占用的所有内存的实现。
 */
void destroy_dataset_impl(DataSet *dataset) {
    if (dataset) {
        free(dataset->points); // 先释放内部数组
        free(dataset);         // 再释放结构体本身
    }
}

/**
 * @brief (不推荐) 从标准输入创建数据集的实现。
 * @note  此函数具有副作用（打印到 stdout，从 stdin 读取），
 *        通常不应在库函数中出现。仅用于测试或特定应用。
 */
Lagrange_Err create_dataset_from_points_impl(DataSet **outDataset, size_t size) {
    if (outDataset == NULL || *outDataset != NULL || size == 0) {
        return LAGRANGE_ERR_INVALID;
    }

    Point* points = (Point*)malloc(size * sizeof(Point));
    if (!points) {
        return LAGRANGE_ERR_NOMEM;
    }

    for (size_t i = 0; i < size; ++i) {
        printf("请输入第 %zu 个点 (格式为 x,y): ", i + 1);
        fflush(stdout);

        if (scanf("%lf,%lf", &points[i].x, &points[i].y) != 2) {
            // 清理输入缓冲区，以防用户输入错误
            int ch;
            while ((ch = getchar()) != '\n' && ch != EOF);
            free(points);
            return LAGRANGE_ERR_INVALID; // 输入格式错误
        }
    }

    Lagrange_Err err = create_dataset_impl(outDataset, points, size);
    free(points); // create_dataset_impl 已复制数据，所以可以释放临时数组
    return err;
}


/** @deprecated  此函数设计存在缺陷，建议废弃。*/
DataSet *empty_dataset_impl(void) {
    return (DataSet *)calloc(1, sizeof(DataSet));
}


/** @deprecated  此函数破坏封装，建议废弃。*/
Point *get_points_impl(DataSet **outDataset) {
    if (!outDataset || !*outDataset) {
        return NULL;
    }
    return (*outDataset)->points;
}
// 公共API
const LagrangeAPI Lagrange = {
    .lagrange_interpolate = lagrange_interpolate_impl,
    .empty_dataset = empty_dataset_impl,
    .create_dataset = create_dataset_impl,
    .create_dataset_from_points = create_dataset_from_points_impl,
    .destroy_dataset = destroy_dataset_impl,
    .get_points = get_points_impl
};
