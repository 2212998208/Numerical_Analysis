#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "newton.h"
#include "lagrange.h"


static Newton_Err newton_create_interpolator(NewtonDataSet **outDataset, DataSet **inDataset, size_t size);
static Newton_Err destroy_dataset(NewtonDataSet **outDataset);
static Newton_Err compute_divided_differences(NewtonDataSet *dataset, Point *points, size_t size);
static Newton_Err newton_evaluate(NewtonDataSet *dataset, double x, double *outY);

// 数值分析Newton的必要的数据结构
struct NewtonDataSet {
    double *divided_differences;
    double *x_nodes;
    size_t x_nodes_size;
};


// Newton插值法多项式计算API
Newton_Err newton_interpolate(NewtonDataSet **inNewtonDataset, DataSet **inDataset,double x, double *outY) {
    if (!*inNewtonDataset || !outY ||
        (*inNewtonDataset)->x_nodes_size == 0 ||
        !*inDataset || !(Lagrange.get_points(inDataset))) {
        return NEWTON_ERR_INVALID;
    }

    if (((*inNewtonDataset)->x_nodes_size) == 1) {
        *outY = Lagrange.get_points(inDataset)[0].y; // 只有一个点，直接返回y值
        return NEWTON_OK;
    }

    Newton_Err err1 = newton_create_interpolator(inNewtonDataset, inDataset, (*inNewtonDataset)->x_nodes_size);
    if (err1 != NEWTON_OK) {
        return err1;
    }

    Newton_Err err2 = newton_evaluate(*inNewtonDataset, x, outY);
    if (err2 != NEWTON_OK) {
        return err2;
    }
    return NEWTON_OK;
}



// 创建插值器
static Newton_Err newton_create_interpolator(NewtonDataSet **outDataset, DataSet **inDataset, const size_t size) {
    if (!*outDataset || !*inDataset || size == 0) {
        return NEWTON_ERR_INVALID;
    }


    return compute_divided_differences(*outDataset, Lagrange.get_points(inDataset), size);

}


// 销毁插值器
Newton_Err destroy_dataset(NewtonDataSet **outDataset) {
    if (!*outDataset) {
        return NEWTON_ERR_INVALID;
    }
    if (!(*outDataset)->divided_differences || !(*outDataset)->x_nodes) {
        free(*outDataset);
        *outDataset = NULL;
        return NEWTON_OK;
    }

    free((*outDataset)->divided_differences);
    free((*outDataset)->x_nodes);
    free(*outDataset);
    *outDataset = NULL;
    return NEWTON_OK;

}




// 计算插值结果
static Newton_Err newton_evaluate(NewtonDataSet *dataset, double x, double *outY) {
    if (!dataset) {
        return NEWTON_ERR_INVALID;
    }
    const int size = (int)dataset->x_nodes_size;
    double result = dataset->divided_differences[size - 1];

    for (int i = size - 1; i > 0; i--) {
        result = result * (x - dataset->x_nodes[i-1]) + dataset->divided_differences[i-1];
    }

    *outY = result;
    return NEWTON_OK;
}



// 打印插值器信息
Newton_Err print_newton_dataset(const NewtonDataSet **outDataset) {
    if (!*outDataset || !(*outDataset)->divided_differences || !(*outDataset)->x_nodes ||
        (*outDataset)->x_nodes_size == 0 || (*outDataset)->x_nodes_size == 1) {
        printf("Dataset is NULL\n");
        return NEWTON_ERR_INVALID;
    }
    printf("Newton Dataset:\n");
    printf("Size: %zu\n", (*outDataset)->x_nodes_size);
    printf("X Nodes: ");
    for (size_t i = 0; i < (*outDataset)->x_nodes_size; i++) {
        printf("%f ", (*outDataset)->x_nodes[i]);
    }
    printf("\nDivided Differences: ");
    for (size_t i = 0; i < (*outDataset)->x_nodes_size; i++) {
        printf("%f ", (*outDataset)->divided_differences[i]);
    }
    printf("\n");

    return NEWTON_OK;
}



// 数据集的构建函数
Newton_Err create_newton_dataset(NewtonDataSet **outDataset, const size_t size) {
    if (!outDataset) {
        return NEWTON_ERR_INVALID;
    }
    NewtonDataSet *dataset = (NewtonDataSet *)malloc(sizeof(NewtonDataSet));
    if (!dataset) {
        return NEWTON_ERR_NOMEM;
    }
    dataset->divided_differences = NULL;
    dataset->x_nodes = NULL;
    dataset->x_nodes_size = size;

    *outDataset = dataset;
    return NEWTON_OK;
}

// 计算均差系数
static Newton_Err compute_divided_differences(NewtonDataSet *dataset, Point *points, const size_t size) {
    if (!dataset || !points || size == 0) {
        return NEWTON_ERR_INVALID;
    }

    // 分配均差系数数组
    dataset->divided_differences = (double *)malloc(size * sizeof(double));
    if (!dataset->divided_differences) {
        return NEWTON_ERR_NOMEM;
    }
    // 分配x节点数组
    dataset->x_nodes = (double *)malloc(size * sizeof(double));
    if (!dataset->x_nodes) {
        free(dataset->divided_differences);
        return NEWTON_ERR_NOMEM;
    }
    dataset->x_nodes_size = size;

    // 初始化均差系数和x节点
    for (size_t i = 0; i < size; i++) {
        dataset->divided_differences[i] = points[i].y; // 零阶均差系数即为y值
        dataset->x_nodes[i] = points[i].x;
    }

    // 计算均差系数
    for (size_t j = 1; j < size; j++) { // j表示当前计算的均差阶数
        for (size_t i = size - 1; i >= j; i--) { // 从后向前计算，避免覆盖前面的值，因为dd[i]依赖于dd[i-1]。
            double numerator = dataset->divided_differences[i] - dataset->divided_differences[i - 1];
            double denominator = dataset->x_nodes[i] - dataset->x_nodes[i - j];
            if (fabs(denominator) < 1e-9 || fabs(denominator) == 0) {
                free(dataset->divided_differences);
                free(dataset->x_nodes);
                return NEWTON_ERR_DIVBYZERO;
            }
            dataset->divided_differences[i] = numerator / denominator;
        }
    }

    return NEWTON_OK;
}


// 公共API
const NewtonAPI Newton = {
    newton_interpolate,
    create_newton_dataset,
    destroy_dataset,
    print_newton_dataset
};