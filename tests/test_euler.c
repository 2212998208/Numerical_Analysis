#include <math.h>
#include <stdio.h>
#include <windows.h>

#include "euler.h"
#define EXP_1 2.71828182845904523536
#define EXP_N1 0.36787944117144232159

// ODES定义
double ode1(double x, double t) {
    return x;  // dx/dt = x, 解为 x(t) = x0 * exp(t - t0)
}

double ode2(double x, double t) {
  return -2 * t * x;  // dx/dt = -2tx, 解为 x(t) = x0 * exp(-t^2 + t0^2)
}

// 用例描述
typedef struct {
    const char *name;
    double (*dx2dt)(double x, double t);
    double x0;
    double t0;
    double Δt;
    size_t max_iter;
    double expected_xn; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    Euler_Err expected_err; // 期望的错误码
}EulerTestCase;

static EulerTestCase cases [] = {
    {"ode1: dx/dt=x", ode1, 1.0, 0.0, 0.1, 10, EXP_1, 0, EULER_OK}, // 期望结果 e^1
    {"ode2: dx/dt=-2tx", ode2, 1.0, 0.0, 0.1, 10, EXP_N1, 0, EULER_OK}, // 期望结果 e^-1
    {"Invalid Δt", ode1, 1.0, 0.0, -0.1, 10, 0.0, 1, EULER_INVALID}, // 错误的Δt
    {"Null function", NULL, 1.0, 0.0, 0.1, 10, 0.0, 1, EULER_INVALID}, // 错误的函数指针
};


int test_euler() {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const EulerTestCase *tc = &cases[i];


        // 构造函数
        Euler euler = NULL;
        Euler modified_euler = NULL;
        Euler_Err err = EULER.euler_create(&euler, tc->dx2dt, tc->x0, tc->t0, tc->Δt);
        Euler_Err modified_err = EULER.euler_create(&modified_euler, tc->dx2dt, tc->x0, tc->t0, tc->Δt);
        if (err != EULER_OK || modified_err != EULER_OK) {
            if (tc->expect_error && err == tc->expected_err && err == EULER_INVALID || modified_err == EULER_INVALID
                && modified_err == tc->expected_err && tc->expect_error) {
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
        double modified_xn = 0.0;
        err = EULER.euler_solve(&euler, tc->max_iter, &xn);
        modified_err = EULER.modified_euler_solve(&modified_euler, tc->max_iter, &modified_xn);
        if (err != EULER_OK || modified_err != EULER_OK) {
            if (tc->expect_error && err == tc->expected_err && err == EULER_INVALID
                || tc->expect_error && modified_err == EULER_INVALID && modified_err == tc->expected_err) {
                printf("[TEST] %s: 求解函数参数错误: err=%d  PASSED\n", tc->name, err);
                passed++;
                EULER.euler_destroy(&euler);
                EULER.euler_destroy(&modified_euler);
                continue;
            }
            printf("[TEST] %s: 求解函数错误: err=%d, expected_err=%d  FAILED\n", tc->name, err, tc->expected_err);
            failed++;
            EULER.euler_destroy(&euler);
            EULER.euler_destroy(&modified_euler);
            continue;
        }

        // 验证结果
        printf("[TEST] %-25s 计算结果: xn=%.10f, 期望结果: %.10f  CALCULATED\n", tc->name, xn, tc->expected_xn);
        printf("[TEST] %-25s 修正结果: modified_xn=%.10f, 期望结果: %.10f  CALCULATED\n", tc->name, modified_xn, tc->expected_xn);
        passed++;

        // 析构函数
        EULER.euler_destroy(&euler);
        EULER.euler_destroy(&modified_euler);
    }

    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return failed == 0 ? 0 : 1;
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_euler();
}