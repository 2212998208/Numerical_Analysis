#ifndef NUMERICAL_ANALYSIS_LAGRANGE_H
#define NUMERICAL_ANALYSIS_LAGRANGE_H
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 数据集结构的前向声明。
 *
 * 具体的定义在 lagrange.c 中，此处仅为声明，以实现信息隐藏。
 * 用户通过指针与此结构交互。
 */
typedef struct DataSet DataSet;

/**
 * @brief 表示二维平面上的一个点 (x, y)。
 *
 * 这是拉格朗日插值所需的基本数据单位。
 */
typedef struct Point {
    double x; /**< 点的 x 坐标 */
    double y; /**< 点的 y 坐标 */
} Point;

/**
 * @brief 拉格朗日插值模块的错误码。
 *
 * 用于标识函数执行过程中可能发生的各种错误。
 */
typedef enum Lagrange_Err {
    LAGRANGE_OK = 0,          /**< 操作成功 */
    LAGRANGE_ERR_NOMEM = 1,     /**< 内存分配失败 */
    LAGRANGE_ERR_INVALID = 2, /**< 无效的输入参数（例如，空指针或大小为零） */
    LAGRANGE_ERR_DIVBYZERO = 3  /**< 发生除零错误（例如，插值点与某个数据点相同） */
} Lagrange_Err;

/**
 * @brief 创建一个 Point 对象的便捷函数。
 *
 * @param x 点的 x 坐标。
 * @param y 点的 y 坐标。
 * @return  一个初始化后的 Point 对象。
 */
static inline Point point_make(const double x, const double y) {
    const Point p = {x, y};
    return p;
}

/**
 * @brief 拉格朗日插值算法的 API 集合。
 *
 * 封装了与拉格朗日插值相关的所有操作，包括数据集的创建、销毁和插值计算。
 */
typedef struct {
    /**
     * @brief 在给定的数据集上执行拉格朗日插值。
     *
     * @param dataset 指向一个已初始化的数据集的指针。
     * @param x       需要计算插值结果的点的 x 坐标。
     * @return        在点 x 处的插值结果 y。如果数据点中存在 x 坐标相同的点，行为未定义。
     */
    double (*lagrange_interpolate)(DataSet *dataset, double x);

    /**
     * @brief 创建一个空的数据集。
     * @deprecated 此函数可能已废弃，建议使用 create_dataset。
     * @return 指向新创建的空数据集的指针，如果失败则返回 NULL。
     */
    DataSet *(*empty_dataset)(void);

    /**
     * @brief 从一个 Point 数组创建并初始化数据集。
     *
     * @param[out] dataset 指向一个指针，该指针将用于存储新创建的数据集的地址。
     * @param points      一个包含插值数据点的数组。
     * @param size        数组中的数据点数量。
     * @return            如果成功，返回 LAGRANGE_OK；否则返回相应的错误码。
     */
    Lagrange_Err (*create_dataset)(DataSet **dataset, Point *points, size_t size);

    /**
     * @brief 根据指定大小创建一个数据集，但不填充具体数据点。
     * @deprecated 此函数可能已废弃或用途特殊，建议使用 create_dataset。
     * @param[out] dataset 指向一个指针，用于存储新数据集的地址。
     * @param size        要为数据集预分配空间的数据点数量。
     * @return            如果成功，返回 LAGRANGE_OK；否则返回错误码。
     */
    Lagrange_Err (*create_dataset_from_points)(DataSet **dataset, size_t size);

    /**
     * @brief 销毁一个数据集并释放其占用的所有内存。
     *
     * @param dataset 指向要销毁的数据集的指针。如果为 NULL，则不执行任何操作。
     */
    void (*destroy_dataset)(DataSet *dataset);

    /**
     * @brief 获取数据集内部的原始数据点数组。
     * @deprecated 此函数可能已废弃或设计不佳，因为它破坏了封装性。
     * @param outDataset 指向数据集的指针。
     * @return 指向内部 Point 数组的指针。
     */
    Point *(*get_points)(DataSet **outDataset);
} LagrangeAPI;


extern const LagrangeAPI Lagrange;

#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_LAGRANGE_H