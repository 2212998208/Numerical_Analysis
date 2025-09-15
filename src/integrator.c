#include "integrator.h"
#include <math.h>
// 本地替代 (不依赖 fabs / fmax)
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }

#ifndef INTEGRATOR_INLINE
#if defined(_MSC_VER)
#define INTEGRATOR_INLINE __forceinline
#else
#define INTEGRATOR_INLINE inline __attribute__((always_inline))
#endif
#endif

// 单步 RK4
static INTEGRATOR_INLINE double rk4_step(IntegrandFn f, void *user,
                                         double x, double y, double h) {
    double k1 = f(x, user);
    double k2 = f(x + 0.5 * h, user);
    double k3 = k2; // 可复用：若函数评估昂贵可再算，这里为示例保持对称性简单
    double k4 = f(x + h, user);
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// 固定步长 RK4 积分 (把积分视作 y' = f(x), y(a)=0)
static double rk4_fixed_impl(IntegrandFn f, void *user, double a, double b, int steps) {
    if (steps <= 0) return 0.0;
    double h = (b - a) / (double)steps;
    double x = a;
    double y = 0.0;
    for (int i = 0; i < steps; ++i) {
        y = rk4_step(f, user, x, y, h);
        x += h;
    }
    return y;
}

// 自适应 RK4 (双步比较): 误差估计 E ≈ (I_{h/2} - I_h) / 15 (阶数4的 Richardson)
static double rk4_adaptive_impl(IntegrandFn f, void *user, double a, double b,
                                AdaptiveConfig cfg, IntegratorStatus *status) {
    if (status) *status = INTEGRATOR_OK;
    if (a == b) return 0.0;
    if (cfg.max_iterations <= 0) cfg.max_iterations = 20;
    if (cfg.abs_tol <= 0.0) cfg.abs_tol = 1e-9;
    if (cfg.rel_tol <= 0.0) cfg.rel_tol = 1e-9;

    int steps = 8; // 初始分段
    double integral_prev = rk4_fixed_impl(f, user, a, b, steps);
    for (int iter = 0; iter < cfg.max_iterations; ++iter) {
        steps *= 2;
        double integral_refined = rk4_fixed_impl(f, user, a, b, steps);
        double error_est = NA_ABS(integral_refined - integral_prev) / 15.0;
        double scale = NA_MAX(cfg.abs_tol, NA_ABS(integral_refined) * cfg.rel_tol);
        if (error_est <= scale) {
            return integral_refined;
        }
        integral_prev = integral_refined;
    }
    if (status) *status = INTEGRATOR_MAX_STEPS_REACHED;
    return integral_prev;
}

// 公共 API
const IntegratorAPI Integrator = {
    rk4_fixed_impl,
    rk4_adaptive_impl
};
