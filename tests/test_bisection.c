#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <windows.h>

#include "bisection.h"

// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))


// 函数定义
double func1(double x) {return x*x - 4;}  // x=+-2
double func2(double x) {return x*x*x - x - 1;} // x=1.324717957244746
double func3(double x) {return cos(x) - x;} // x=0.7390851332151607

// 用例描述
typedef struct {
    const char *name;
    double (*func)(double x);
    double a;
    double b;
    double tol;
    double expected_root; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    Bisection_Err expected_err; // 期望的错误码
}BisectionTestCase;


static BisectionTestCase cases [] = {
    {"x^2 - 4", func1, 0, 3, 1e-9, 2.0, 0, BISECTION_OK},
    {"x^2 - 4", func1, -3, 0, 1e-9, -2.0, 0, BISECTION_OK},
    {"x^3 - x - 1", func2, 1, 2, 1e-9, 1.324717957244746, 0, BISECTION_OK},
    {"cos(x) - x", func3, 0, 1, 1e-9, 0.7390851332151607, 0, BISECTION_OK},
    //{"无根区间", func1, 3, 4, 1e-9, 0.0, 1, BISECTION_ERR_INVALID},
    //{"无效区间", func1, 4, 3, 1e-9, 0.0, 1, BISECTION_ERR_INVALID},
    //{"过小精度", func1, 0, 3, 1e-20, 0.0, 1, BISECTION_ERR_INVALID}
};

int test_bisection(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const BisectionTestCase *tc = &cases[i];

        // 构造函数
        BisectionRange range = NULL;
        Bisection_Err err = Bisection.bisection_create(tc->func, tc->a, tc->b, tc->tol, &range, tc->name);
        if (err != BISECTION_OK) {
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-8s 构建区间期望失败: err=%d PASS\n", tc->name, err);
                passed++;
            } else {
                printf("[TEST] %-8s 构建区间失败: err=%d FAIL\n", tc->name, err);
                failed++;
            }
        }

        // 数值解
        err = Bisection.bisection_solve(&range);
        if (err != BISECTION_OK) {
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-8s 求解期望失败: err=%d PASS\n", tc->name, err);
                passed++;
            } else {
                printf("[TEST] %-8s 求解失败: err=%d FAIL\n", tc->name, err);
                failed++;
            }
        } else {
            if (TEST_ABS_REL_CLOSE(Bisection.bisection_get_midpoint(&range), tc->expected_root, tc->tol, tc->tol)) {
                printf("[TEST] %-8s 求解成功: root=%.15f PASS\n", tc->name, Bisection.bisection_get_midpoint(&range));
                passed++;
            } else {
                printf("[TEST] %-8s 求解结果错误: root=%.15f expected=%.15f FAIL\n", tc->name,
                       Bisection.bisection_get_midpoint(&range), tc->expected_root);
                failed++;
            }
        }

        // 析构函数
        err = Bisection.bisection_destroy(&range);
        if (err != BISECTION_OK) {
            printf("[TEST] %-8s 析构失败: err=%d FAIL\n", tc->name, err);
            failed++;
        }
    }

    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return (passed == (passed + failed)) ? 0 : 1;
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_bisection();
}
