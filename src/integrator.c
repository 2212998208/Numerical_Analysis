#include "integrator.h"
#include <math.h>
// 为避免引入 <math.h> 带来的依赖，定义本地的绝对值和最大值函数。
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }

// 定义跨平台的强制内联宏，以优化性能。
#ifndef INTEGRATOR_INLINE
#if defined(_MSC_VER)
#define INTEGRATOR_INLINE __forceinline
#else
#define INTEGRATOR_INLINE inline __attribute__((always_inline))
#endif
#endif

/**
 * @brief 执行单步四阶龙格-库塔 (RK4) 计算。
 *
 * 这是求解常微分方程 y' = f(x, y) 的核心步骤。对于定积分 y' = f(x)，
 * 此函数计算从 x 到 x+h 的增量。
 * @note 这里的 k2 和 k3 相同，因为 f(x) 不依赖于 y。
 *       如果求解 f(x,y)，则 k3 = f(x + 0.5 * h, y + 0.5 * h * k2)。
 */
static INTEGRATOR_INLINE double rk4_step(IntegrandFn f, void *user,
                                         double x, double y, double h) {
    // 计算四个斜率估计值
    double k1 = f(x, user);
    double k2 = f(x + 0.5 * h, user);
    double k3 = f(x + 0.5 * h, user); // 对于 y'=f(x)，k3 与 k2 相同
    double k4 = f(x + h, user);
    // 根据加权平均更新 y 值
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

/**
 * @brief 固定步长 RK4 积分的内部实现。
 *
 * 将积分区间 [a, b] 分为 `steps` 个相等的部分，并累加每一步的 RK4 计算结果。
 * 积分问题被看作是求解初始值问题：y' = f(x)，y(a) = 0。
 */
static double rk4_fixed_impl(IntegrandFn f, void *user, double a, double b, int steps) {
    if (steps <= 0) return 0.0; // 防御性编程：步数为非正数时返回0
    double h = (b - a) / (double)steps;
    double x = a;
    double y = 0.0; // 积分初始值为0
    for (int i = 0; i < steps; ++i) {
        y = rk4_step(f, user, x, y, h);
        x += h;
    }
    return y;
}

/**
 * @brief 自适应步长 RK4 积分的内部实现。
 *
 * 采用 Richardson 外推法进行误差估计。通过比较使用步长 h 和 h/2（即双倍步数）
 * 计算出的两个积分值，来估计截断误差。如果误差小于容限，则接受结果；
 * 否则，加倍步数，继续迭代，直到满足精度或达到最大迭代次数。
 *
 * 误差公式 E ≈ |I_refined - I_prev| / (2^p - 1)，其中 p=4 是 RK4 的阶数。
 * 因此分母为 15。
 */
static double rk4_adaptive_impl(IntegrandFn f, void *user, double a, double b,
                                AdaptiveConfig cfg, IntegratorStatus *status) {
    // 初始化状态和参数
    if (status) *status = INTEGRATOR_OK;
    if (a == b) return 0.0;
    if (cfg.max_iterations <= 0) cfg.max_iterations = 20; // 默认最大迭代次数
    if (cfg.abs_tol <= 0.0) cfg.abs_tol = 1e-9;           // 默认绝对容限
    if (cfg.rel_tol <= 0.0) cfg.rel_tol = 1e-9;           // 默认相对容限

    int steps = 8; // 初始步数，可以根据问题调整
    double integral_prev = rk4_fixed_impl(f, user, a, b, steps);

    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        steps *= 2; // 步数加倍，步长减半
        double integral_refined = rk4_fixed_impl(f, user, a, b, steps);

        // Richardson 误差估计
        double error_est = NA_ABS(integral_refined - integral_prev) / 15.0;
        // 确定容限标度
        double scale = NA_MAX(cfg.abs_tol, NA_ABS(integral_refined) * cfg.rel_tol);

        // 检查是否满足精度要求
        if (error_est <= scale) {
            return integral_refined; // 结果收敛，返回更精确的值
        }
        integral_prev = integral_refined; // 未收敛，准备下一次迭代
    }

    // 达到最大迭代次数，标记状态并返回当前最佳结果
    if (status) *status = INTEGRATOR_MAX_STEPS_REACHED;
    return integral_prev;
}

// 公共 API
const IntegratorAPI Integrator = {
    rk4_fixed_impl,
    rk4_adaptive_impl
};
