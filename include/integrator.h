#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 通用被积函数指针类型。
 *
 * 定义了一个函数原型，用于表示被积函数 f(x)。
 * 这种设计允许用户传入任何匹配此签名的函数进行积分。
 *
 * @param x           自变量，积分点。
 * @param user_data   一个可选的、指向用户自定义数据的指针。
 *                    当被积函数需要额外参数时，可通过此指针传递，避免使用全局变量。
 *                    如果不需要，可以为 NULL。
 * @return            函数在点 x 处的值 f(x)。
 */
typedef double (*IntegrandFn)(double x, void *user_data);

/**
 * @brief 自适应积分过程的状态码。
 *
 * 用于标识自适应积分算法的执行结果。
 */
typedef enum {
    /** 积分成功完成。 */
    INTEGRATOR_OK = 0,
    /** 达到最大迭代次数，积分可能未收敛到期望精度。 */
    INTEGRATOR_MAX_STEPS_REACHED = 1
} IntegratorStatus;

/**
 * @brief 自适应积分算法的配置参数。
 *
 * 存储控制自适应积分过程的精度和迭代限制。
 */
typedef struct {
    /** 绝对误差容限 (Absolute Tolerance)。积分的绝对误差应小于此值。 */
    double abs_tol;
    /** 相对误差容限 (Relative Tolerance)。积分的相对误差应小于此值。 */
    double rel_tol;
    /**
     * @brief 最大迭代（细分）次数。
     *
     * 为防止无限递归或计算时间过长，限制了区间细分的最大深度。
     * 例如，如果初始区间数为 N，最大迭代次数为 M，则最多会计算 N * 2^M 个子区间。
     */
    int max_iterations;
} AdaptiveConfig;

/**
 * @brief 数值积分算法的 API 集合。
 *
 * 使用结构体封装一组函数指针，提供一个统一的、可扩展的积分器接口。
 * 所有与积分相关的操作都通过这个结构的实例（如全局的 Integrator 对象）进行调用。
 */
typedef struct {
    /**
     * @brief 使用固定步长的四阶龙格-库塔 (RK4) 方法计算定积分。
     *
     * @param f       被积函数，符合 IntegrandFn 类型。
     * @param user    传递给被积函数的自定义数据指针，可为 NULL。
     * @param a       积分下限。
     * @param b       积分上限。
     * @param steps   积分区间的总步数，必须为正数。步长 h = (b - a) / steps。
     * @return        定积分的近似值。
     */
    double (*rk4_fixed)(IntegrandFn f, void *user, double a, double b, int steps);

    /**
     * @brief 使用自适应步长的四阶龙ge-库塔 (RK4) 方法计算定积分。
     *
     * 该方法会根据误差估计自动调整步长，以在满足用户指定的精度容限下完成积分。
     * 它采用双网格（步长 h 和 h/2）的 Richardson 误差估计来控制精度。
     *
     * @param f       被积函数，符合 IntegrandFn 类型。
     * @param user    传递给被积函数的自定义数据指针，可为 NULL。
     * @param a       积分下限。
     * @param b       积分上限。
     * @param cfg     自适应算法的配置，包括误差容限和最大迭代次数。
     * @param status  一个指向 IntegratorStatus 枚举的指针，用于返回积分过程的状态。
     *                如果积分成功，其值将被设为 INTEGRATOR_OK；
     *                如果达到最大迭代次数仍未收敛，则为 INTEGRATOR_MAX_STEPS_REACHED。
     *                此参数不能为 NULL。
     * @return        定积分的近似值。如果无法完成，可能返回 NAN 或部分计算结果。
     */
    double (*rk4_adaptive)(IntegrandFn f, void *user, double a, double b,
                           AdaptiveConfig cfg, IntegratorStatus *status);
} IntegratorAPI;

// 全局只读实例
extern const IntegratorAPI Integrator;

#ifdef __cplusplus
}
#endif

#endif

