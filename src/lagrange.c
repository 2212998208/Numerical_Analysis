#include "lagrange.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


static double calculate_basis_polynomial(const DataSet *dataset, int k, double x);

// 数值分析Lagrange的必要的数据结构
struct DataSet {
    Point *points;
    size_t size;
};


// Lagrange插值法多项式计算API
double lagrange_interpolate(DataSet *dataset, const double x) {
    if (!dataset || dataset->size == 0 || !dataset->points) {
        return LAGRANGE_ERR_INVALID;
    }

    if (dataset->size == 1) {
        return (dataset->points)[0].y; // 只有一个点，直接返回y值
    }

    // Lagrange插值多项式计算核心算法
    double total_sum =0.0;
    for (size_t i = 0; i < dataset->size; ++i) {
        double y_i = (dataset->points)[i].y;

        double L_k = calculate_basis_polynomial(dataset, (int)i, x);
        if ((Lagrange_Err)L_k == LAGRANGE_ERR_DIVBYZERO) {
            return LAGRANGE_ERR_DIVBYZERO; // 遇到除零错误，返回错误码
        }

        total_sum += y_i * L_k;
    }

    return total_sum;
}

// 计算Lagrange基函数
static double calculate_basis_polynomial(const DataSet *dataset, const int k, const double x) {
    double product = 1.0;
    double x_k = (dataset->points)[k].x;

    for (size_t i = 0; i < dataset->size; ++i) {
        if ((int)i == k) continue;

        double x_i = (dataset->points)[i].x;

        double denominator = x_k - x_i;
        if (fabs(denominator) < 1e-9) {
            return LAGRANGE_ERR_DIVBYZERO;
        }
        double numerator = x - x_i;
        product = product * (numerator / denominator);
    }

    return product;
}

// 数据集构建函数
Lagrange_Err create_dataset(DataSet **outDataset, Point *points, size_t size) {
    if (!outDataset || !points || size == 0) {
        return LAGRANGE_ERR_INVALID;
    }
    DataSet *dataset = (DataSet *)malloc(sizeof(DataSet));
    if (!dataset) {
        return LAGRANGE_ERR_NOMEM;
    }
    dataset->points = (struct Point *)malloc(size * sizeof(struct Point));
    if (!dataset->points) {
        free(dataset);
        return LAGRANGE_ERR_NOMEM;
    }
    for (size_t i = 0; i < size; ++i) {
        (dataset->points)[i] = points[i];
    }
    dataset->size = size;
    *outDataset = dataset;
    return LAGRANGE_OK;
}

// 数据集销毁函数
void destroy_dataset(DataSet *dataset) {
    if (dataset) {
        free(dataset->points);
        free(dataset);
    }
}

// 从点集的创建数据集
Lagrange_Err create_dataset_from_points(DataSet **outDataset, size_t size) {
    double x,y;


    if (!*outDataset || size == 0) {
        return LAGRANGE_ERR_INVALID;
    }
    DataSet *dataset = (DataSet *)malloc(sizeof(DataSet));
    if (!dataset) {
        return LAGRANGE_ERR_NOMEM;
    }
    dataset->points = (struct Point *)malloc(size * sizeof(struct Point));
    if (!dataset->points) {
        free(dataset);
        return LAGRANGE_ERR_NOMEM;
    }
    for (size_t i = 0; i < size; ++i) {
        // 获取你的点集数据
        printf("请输入第%zu个点格式为 x,y:", i + 1);
        fflush(stdout); // 确保提示信息立即输出


        if (scanf("%lf,%lf", &x, &y) != 2) {
            int ch;
            while ((ch = getchar()) != '\n' && ch != EOF);
            destroy_dataset(dataset);
            return LAGRANGE_ERR_INVALID;
        }
        (dataset->points)[i] = point_make(x,y); // 初始化y为0.0
    }
    dataset->size = size;
    *outDataset = dataset;
    return LAGRANGE_OK;
}

// 构建空数据集函数
DataSet *empty_dataset(void) {
    return (DataSet *)malloc(sizeof(DataSet));
}
// 公共API
const LagrangeAPI Lagrange = {
    lagrange_interpolate,
    empty_dataset,
    create_dataset,
    create_dataset_from_points,
    destroy_dataset
};
