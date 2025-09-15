#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// 通用被积函数: x -> f(x); user_data 便于传递参数 (可为 NULL)
typedef double (*IntegrandFn)(double x, void *user_data);

// 自适应积分状态
typedef enum {
    INTEGRATOR_OK = 0,
    INTEGRATOR_MAX_STEPS_REACHED = 1
} IntegratorStatus;

// 自适应配置
typedef struct {
    double abs_tol;        // 绝对误差容限
    double rel_tol;        // 相对误差容限
    int    max_iterations; // 细分最大轮次 (指数加深)
} AdaptiveConfig;

// API 结构 (可扩展更多算法)
typedef struct {
    // 固定步长 RK4: steps 为区间总步数
    double (*rk4_fixed)(IntegrandFn f, void *user, double a, double b, int steps);
    // 自适应 RK4: 返回积分值; status 输出状态
    double (*rk4_adaptive)(IntegrandFn f, void *user, double a, double b,
                           AdaptiveConfig cfg, IntegratorStatus *status);
} IntegratorAPI;

// 全局只读实例
extern const IntegratorAPI Integrator;

#ifdef __cplusplus
}
#endif

#endif

