
#ifndef NUMERICAL_ANALYSIS_HERMITE_H
#define NUMERICAL_ANALYSIS_HERMITE_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Hermite 插值数据点结构的前向声明。
 *
 * 包含一个点的 x 坐标、y 坐标和一阶导数值。
 */
typedef struct HermitePoint* HermitePoint;

/**
 * @brief Hermite 插值数据集结构的前向声明。
 *
 * 包含一组 HermitePoint。
 */
typedef struct HermiteDataset* HermiteDataset;

/**
 * @brief Hermite 插值器结构的前向声明。
 *
 * 存储插值多项式的系数或其他预计算结果，用于高效求值。
 */
typedef struct HermiteInterpolator* HermiteInterpolator;

/**
 * @brief Hermite 插值模块的错误码。
 */
typedef enum Hermite_Err {
    HERMITE_OK = 0,          /**< 操作成功 */
    HERMITE_ERR_NOMEM = 1,     /**< 内存分配失败 */
    HERMITE_ERR_INVALID = 2, /**< 无效的输入参数 */
    HERMITE_ERR_DIVBYZERO = 3  /**< 发生除零错误 */
} Hermite_Err;

/**
 * @brief Hermite 插值算法的 API 集合。
 */
typedef struct {
    /**
     * @brief 从数据集创建 Hermite 插值器。
     *
     * 此函数预计算插值多项式所需的系数（如差商），并存储在插值器中。
     *
     * @param[in]  inDataset       指向包含 Hermite 数据点的数据集的指针。
     * @param[out] outInterpolator 指向一个指针，用于存储新创建的插值器。
     * @return                     如果成功，返回 HERMITE_OK；否则返回错误码。
     */
    Hermite_Err (*hermite_create_interpolator)(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator);

    /**
     * @brief 销毁 Hermite 插值器并释放内存。
     *
     * @param[in, out] inInterpolator 指向要销毁的插值器指针的指针。成功后，指针将被设为 NULL。
     * @return                        如果成功，返回 HERMITE_OK。
     */
    Hermite_Err (*hermite_destroy_interpolator)(HermiteInterpolator *inInterpolator);

    /**
     * @brief 从原始数据数组创建 Hermite 数据集。
     *
     * @param[in] size 数据点的数量。
     * @param[in] x    指向包含 x 坐标的数组的指针。
     * @param[in] y    指向包含 y 坐标（函数值）的数组的指针。
     * @param[in] dy   指向包含一阶导数值的数组的指针。
     * @return         如果成功，返回指向新创建的数据集的指针；否则返回 NULL。
     */
    HermiteDataset (*create_hermite_dataset)(size_t size, const double *x, const double *y, const double *dy);

    /**
     * @brief 销毁 Hermite 数据集并释放内存。
     *
     * @param[in, out] inDataset 指向要销毁的数据集指针的指针。成功后，指针将被设为 NULL。
     * @return                   如果成功，返回 HERMITE_OK。
     */
    Hermite_Err (*destroy_hermite_dataset)(HermiteDataset *inDataset);

    /**
     * @brief 使用插值器在指定点进行求值。
     *
     * @param[in]  inInterpolator 指向已创建的 Hermite 插值器的指针。
     * @param[in]  x              要求值的点的 x 坐标。
     * @param[out] outY           指向一个 double 变量的指针，用于存储插值函数值。
     * @param[out] outDy          指向一个 double 变量的指针，用于存储插值导数值。
     * @return                    如果成功，返回 HERMITE_OK；否则返回错误码。
     */
    Hermite_Err (*hermite_evaluate)(const HermiteInterpolator *inInterpolator, double x, double *outY, double *outDy);
} HermiteAPI;







extern const HermiteAPI Hermite;
#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_HERMITE_H