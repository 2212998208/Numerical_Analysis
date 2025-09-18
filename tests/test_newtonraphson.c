#include <math.h>
#include <stdio.h>
#include <windows.h>

#include "newton_raphson.h"



// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))

// 函数定义
double func1(double x) {return x*x - 4;}  // x=+-2
double func2(double x) {return x*x*x - x - 1;} // x=1.324717957244746
double func3(double x) {return cos(x) - x;} // x=0.7390851332151607
double func4(double x) {return exp(x) - x - 1;} // x=0
double func5(double x) {return exp(2*x/M_PI) - exp(sin(x));} // x=0
// 用例描述
typedef struct {
    const char *name;
    double (*func)(double x);
    double x0;
    double tol;
    size_t max_iter;
    double expected_root; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    NewtonRaphson_Err expected_err; // 期望的错误码
}BisectionTestCase;


static BisectionTestCase cases [] = {
    {"x^2 - 4", func1, 0.1, 1e-2, 8, 2.0, 0, NEWTON_RAPHSON_OK},
    {"x^2 - 4", func1, -3, 1e-2, 8, -2.0, 0, NEWTON_RAPHSON_OK},
    {"x^3 - x - 1", func2, 1, 1e-2, 8, 1.324717957244746, 0, NEWTON_RAPHSON_OK},
    {"cos(x) - x", func3, 0.1, 1e-2, 8, 0.7390851332151607, 0, NEWTON_RAPHSON_OK},
    {"e^x - x - 1", func4, 1, 1e-9, 1024, 0.0, 0, NEWTON_RAPHSON_OK},  // 敏感边界
    {"e^(2x/pi) - esin(x)",func5, 1.2, 1e-5, 1024, M_PI/2, 0, NEWTON_RAPHSON_OK}
};


int test_newtonraphson(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const BisectionTestCase *tc = &cases[i];

        // 构造函数
        NonLinearRange range = NULL;
        NewtonRaphson_Err err = NewtonRaphson.NonLinearRange_create(tc->func,tc->x0,tc->tol,tc->max_iter,&range,tc->name);
        if (err != NEWTON_RAPHSON_OK) {
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-8s 构建区间期望失败: err=%d PASS\n", tc->name, err);
                passed++;
            } else {
                printf("[TEST] %-8s 构建区间失败: err=%d FAIL\n", tc->name, err);
                failed++;
            }
        }

        // 数值解
        double root = 0.0;
        err = NewtonRaphson.newton_raphson_solve(&range, &root);
        if (err != NEWTON_RAPHSON_OK) {
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-8s 求解期望失败: err=%d PASS\n", tc->name, err);
                passed++;
            } else {
                printf("[TEST] %-8s 求解失败: err=%d FAIL\n", tc->name, err);
                failed++;
            }
        } else {
            if (tc->expect_error) {
                printf("[TEST] %-8s 求解期望失败但成功了: root=%g FAIL\n", tc->name, root);
                failed++;
            } else if (TEST_ABS_REL_CLOSE(root, tc->expected_root, tc->tol, tc->tol)) {
                printf("[TEST] %-8s 求解成功: root=%g PASS\n", tc->name, root);
                passed++;
            } else {
                printf("[TEST] %-8s 求解结果错误: root=%g expected=%g FAIL\n", tc->name, root, tc->expected_root);
                failed++;
            }
        }

        // 析构函数
        err = NewtonRaphson.NonLinearRange_destroy(&range);
        if (err != NEWTON_RAPHSON_OK) {
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
    return test_newtonraphson();
}