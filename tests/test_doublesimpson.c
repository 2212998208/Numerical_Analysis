#include <math.h>

#include "double_simpson.h"
#include <stdio.h>
#include <windows.h>

#define EXP_1 1.7182818284590452353602874713527
// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))


// 被积函数定义
double f1(double x, double y) {
    return pow(x, 2) + pow(y, 3);
}

double f2(double x, double y) {
    return sin(x) + cos(y);
}

double f3(double x, double y) {
    return sqrt(x) + sqrt(y);
}

double f4(double x, double y) {
    return exp(x) * exp(y);
}




// 用例描述
typedef struct {
    const char *name;
    double (*f)(double x, double y);
    double x_a;
    double x_b;
    double y_c;
    double y_d;
    size_t n;
    size_t m;
    double expected_integral; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    Simpson_Err expected_err; // 期望的错误码
}DoubleSimpsonTestCase;

static DoubleSimpsonTestCase cases [] = {
    {"f1(x,y)=x^2+y^3", f1, 0.0, 1.0, 1.0, 2.0, 2, 2, 49.0/12.0, 0, SIMPSON_OK},
    {"f2(x,y)=sin(x)+cos(y)", f2, 0.0, M_PI, 0.0, M_PI, 200, 200, 2.0 * M_PI, 0, SIMPSON_OK},
    {"f3(x,y)=sqrt(x)+sqrt(y)", f3, 0.0, 1.0, 0.0, 1.0, 10000, 10000, 4.0/3.0, 0, SIMPSON_OK},
    {"f4(x,y)=exp(x)*exp(y)", f4, 0.0, 1.0, 0.0, 1.0, 100, 100, EXP_1*EXP_1, 0, SIMPSON_OK},
    {"f1(x,y)=x^2+y^3", f1, 1.0, 0.0, 1.0, 2.0, 2, 2, 49.0/12.0, 1, SIMPSON_ERR_INVAL}, // 错误的积分区间
    {"f1(x,y)=x^2+y^3", f1, 0.0, 1.0, 1.0, 2.0, 3, 2, 49.0/12.0, 1, SIMPSON_ERR_INVAL}, // 错误的n
    {"f1(x,y)=x^2+y^3", f1, 0.0, 1.0, 1.0, 2.0, 2, 3, 49.0/12.0, 1, SIMPSON_ERR_INVAL}, // 错误的m
};


int test_doublesimpson(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const DoubleSimpsonTestCase *tc = &cases[i];

        // 构造函数
        DoubleSimpson DoubleSimpson = NULL;
        Simpson_Err err = DS.double_simpson_create(tc->f, tc->x_a, tc->x_b, tc->y_c, tc->y_d, tc->n, tc->m, &DoubleSimpson, tc->name);
        if (err != SIMPSON_OK) {
            if (tc->expect_error && err == tc->expected_err && err == SIMPSON_ERR_INVAL) {
                printf("[TEST] %-25s 构造函数参数错误: err=%d PASS\n", tc->name, err);
                passed++;
                continue;
            }
            printf("[TEST] %-25s 构造函数失败: err=%d FAIL\n", tc->name, err);
            failed++;
            continue;
        }

        // 执行积分
        double result = 0.0;
        err = DS.double_simpson_integrate(&DoubleSimpson, &result);
        if (err != SIMPSON_OK) {
            if (tc->expect_error && err == tc->expected_err && err == SIMPSON_ERR_INVAL) {
                printf("[TEST] %-25s 积分函数参数错误: err=%d PASS\n", tc->name, err);
                DS.double_simpson_destroy(&DoubleSimpson);
                passed++;
                continue;
            }
            printf("[TEST] %-25s 积分函数失败: err=%d FAIL\n", tc->name, err);
            DS.double_simpson_destroy(&DoubleSimpson);
            failed++;
            continue;
        }

        // 验证结果
        if (tc->expect_error) {
            printf("[TEST] %-25s 期望错误但未发生 FAIL\n", tc->name);
            DS.double_simpson_destroy(&DoubleSimpson);
            failed++;
            continue;
        }

        if (TEST_ABS_REL_CLOSE(result, tc->expected_integral, 1e-6, 1e-6)) {
            printf("[TEST] %-25s 积分结果正确: result=%.10f PASS\n", tc->name, result);
            passed++;
        } else {
            printf("[TEST] %-25s 积分结果错误: result=%.10f expected=%.10f FAIL\n", tc->name, result, tc->expected_integral);
            failed++;
        }
        // 析构函数
        DS.double_simpson_destroy(&DoubleSimpson);
    }
    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return failed == 0 ? 0 : 1;
}

int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_doublesimpson();
}