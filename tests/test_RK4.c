#include "RK4.h"
#include <math.h>
#include <stdio.h>
#include <windows.h>


#define EXP_1 2.71828182845904523536
#define EXP_N1 0.36787944117144232159


// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))


// ODES定义
double ode1(double x, double t) {
    return x;  // dx/dt = x, 解为 x(t) = x0 * exp(t - t0)
}

double ode2(double x, double t) {
  return -2 * t * x;  // dx/dt = -2tx, 解为 x(t) = x0 * exp(-t^2 + t0^2)
}

double ode3(double x, double t) {
    return exp(-x) - sin(x) + sqrt(x); // dx/dt = exp(-x) - sin(x) + sqrt(-x), 仅在x>=0时有意义
}
// 用例描述
typedef struct {
    const char *name;
    double (*dx2dt)(double x, double t);
    double x0;
    double t0;
    double h;
    size_t max_iter;
    double expected_xn; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    RK4_Err expected_err; // 期望的错误码
}RK4TestCase;

static RK4TestCase cases [] = {
    {"ode1: dx/dt=x", ode1, 1.0, 0.0, 0.01, 100, EXP_1, 0, RK4_OK}, // 期望结果 e^1
    {"ode2: dx/dt=-2tx", ode2, 1.0, 0.0, 0.01, 100, EXP_N1, 0, RK4_OK}, // 期望结果 e^-1
    {"ode3: dx/dt=exp(-x)-sin(x)+sqrt(x)", ode3, 1.0, 0.0, 0.01, 100, 1.4665701099, 0, RK4_OK}, // 期望结果 -1
    {"Invalid max_iter", ode1, 1.0, 0.0, 0.1, 0, 0.0, 1, RK4_INVALID}, // 错误的最大迭代次数
    {"Invalid h", ode1, 1.0, 0.0, -0.1, 10, 0.0, 1, RK4_INVALID}, // 错误的Δt
    {"Null function", NULL, 1.0, 0.0, 0.1, 10, 0.0, 1, RK4_INVALID}, // 错误的函数指针
};


int test_rk4() {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const RK4TestCase *tc = &cases[i];


        // 构造函数
        rk rk4 = NULL;
        RK4_Err err = RK4.create(&rk4, tc->dx2dt, tc->t0, tc->x0, tc->h, tc->max_iter);
        if (err != RK4_OK) {
            if (tc->expect_error && err == tc->expected_err && err == RK4_INVALID) {
                printf("[TEST] %s: 构造函数参数错误: err=%d  PASSED \n", tc->name, err);
                passed++;
                continue;
            }
            printf("[TEST] %s: 构造函数错误: err=%d, expected_err=%d  FAILED\n", tc->name, err, tc->expected_err);
            failed++;
            continue;
        }

        // 求解
        double xn = 0.0;
        err = RK4.solve(&rk4, &xn);
        if (err != RK4_OK) {
            if (tc->expect_error && err == tc->expected_err && err == RK4_INVALID) {
                printf("[TEST] %s: 求解函数参数错误: err=%d  PASSED\n", tc->name, err);
                passed++;
                RK4.destroy(&rk4);
                continue;
            }
            printf("[TEST] %s: 求解函数错误: err=%d, expected_err=%d  FAILED\n", tc->name, err, tc->expected_err);
            failed++;
            RK4.destroy(&rk4);
            continue;
        }

        // 验证结果
        if (TEST_ABS_REL_CLOSE(xn, tc->expected_xn, 1e-8, 1e-8)) {
            printf("[TEST] %-25s 计算结果: xn=%.10f, 期望结果: %.10f  CALCULATED\n", tc->name, xn, tc->expected_xn);
            passed++;
        }
        else {
            printf("[TEST] %-25s 计算结果错误: xn=%.10f, 期望结果: %.10f  FAILED\n", tc->name, xn, tc->expected_xn);
            failed++;
        }

        // 析构函数
        RK4.destroy(&rk4);
    }

    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return failed == 0 ? 0 : 1;
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_rk4();
}