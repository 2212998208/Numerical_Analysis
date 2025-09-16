#ifndef NUMERICAL_ANALYSIS_NEWTON_H
#define NUMERICAL_ANALYSIS_NEWTON_H
#include "lagrange.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 牛顿插值数据集结构的前向声明。
 *
 * 具体的定义在 newton.c 中，对用户隐藏实现细节。
 */
typedef struct NewtonDataSet NewtonDataSet;

/**
 * @brief 牛顿插值模块的错误码。
 */
typedef enum Newton_Err {
    NEWTON_OK = 0,          /**< 操作成功 */
    NEWTON_ERR_NOMEM = 1,     /**< 内存分配失败 */
    NEWTON_ERR_INVALID = 2, /**< 无效的输入参数 */
    NEWTON_ERR_DIVBYZERO = 3  /**< 发生除零错误 */
} Newton_Err;

/**
 * @brief 牛顿插值算法的 API 集合。
 *
 * 封装了与牛顿插值相关的所有操作。
 */
typedef struct {
    /**
     * @brief 执行牛顿插值计算。
     *
     * 该函数首先根据输入的数据点（intData）构建差商表（存储在 inNewtonDataSet 中），
     * 然后用该表来计算在点 x 处的插值 y 值。
     *
     * @param[in, out] inNewtonDataSet 指向牛顿数据集的指针。它将被用来存储和重用差商表。
     * @param[in]      intData         指向包含原始数据点 (x, y) 的拉格朗日数据集的指针。
     * @param[in]      x               需要进行插值的点的 x 坐标。
     * @param[out]     outY            指向一个 double 变量的指针，用于存储计算出的插值结果 y。
     * @return                         如果成功，返回 NEWTON_OK；否则返回错误码。
     */
    Newton_Err (*newton_interpolate)(NewtonDataSet **inNewtonDataSet, DataSet **intData, double x, double *outY);

    /**
     * @brief 创建一个用于牛顿插值的数据集。
     *
     * @param[out] outDataset 指向一个指针，该指针将存储新创建的牛顿数据集的地址。
     * @param[in]  size       数据点的数量，用于为差商表分配内存。
     * @return                如果成功，返回 NEWTON_OK；否则返回错误码。
     */
    Newton_Err (*create_newton_dataset)(NewtonDataSet **outDataset, size_t size);

    /**
     * @brief 销毁一个牛顿数据集并释放相关内存。
     *
     * @param[in, out] outDataset 指向要销毁的数据集指针的指针。成功后，指针将被设为 NULL。
     * @return                    如果成功，返回 NEWTON_OK。
     */
    Newton_Err (*destroy_dataset)(NewtonDataSet **outDataset);

    /**
     * @brief 打印（或记录）牛顿数据集的内容，通常是差商表。
     * @note 此函数主要用于调试。
     *
     * @param[in] outDataset 指向要打印的数据集的指针。
     * @return               如果成功，返回 NEWTON_OK。
     */
    Newton_Err (*print_newton_dataset)(const NewtonDataSet **outDataset);
} NewtonAPI;


extern const NewtonAPI Newton;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_NEWTON_H