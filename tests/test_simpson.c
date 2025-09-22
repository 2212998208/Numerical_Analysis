#include <math.h>
#include <stdio.h>
#include <windows.h>

#include "simpson.h"
#define Log_2 0.693147180559945323846
// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))

// 被积函数定义
double f1(double x) {return x * x;} // 积分区间[0,1], 积分值=1/3
double f2(double x) {return sin(x);} // 积分区间[0,π], 积分值=2
double f3(double x) {return exp(x);} // 积分区间[0,1], 积分值=e-1
double f4(double x) {return 1.0 / (1.0 + x * x);} // 积分区间[0,1], 积分值=π/4
double f5(double x) {return log(x + 1);} // 积分区间[0,1], 积分值=2ln(2)-1
double f6(double x) {return sqrt(x);} // 积分区间[0,1], 积分值=2/3

// 用例描述
typedef struct {
    const char *name;
    double (*f)(double x);
    double a;
    double b;
    size_t max_iter;
    double expected_integral; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    Simpson_Err expected_err; // 期望的错误码
}TrapTestCase;


static TrapTestCase cases [] = {
    {"f1(x)=x^2", f1, 0.0, 1.0, 4, 1.0/3.0, 0, SIMPSON_OK},
    {"f2(x)=sin(x)", f2, 0.0, M_PI, 100, 2.0, 0, SIMPSON_OK},
    {"f3(x)=exp(x)", f3, 0.0, 1.0, 100, M_E - 1.0, 0, SIMPSON_OK},
    {"f4(x)=1/(1+x^2)", f4, 0.0, 1.0, 10, M_PI / 4.0, 0, SIMPSON_OK},
    {"f5(x)=log(x+1)", f5, 0.0, 1.0, 100, 2 * Log_2 - 1.0, 0, SIMPSON_OK},
    {"f6(x)=sqrt(x)", f6, 0.0, 1.0, 200000, 2.0 / 3.0, 0, SIMPSON_OK},
    {"f1(x)=x^2", f1, 1.0, 0.0, 10, 1.0/3.0, 1, SIMPSON_ERR_INVAL}, // 错误的积分区间
    {"f1(x)=x^2", f1, 0.0, 1.0, 0, 1.0/3.0, 1, SIMPSON_ERR_INVAL}, // 错误的最大迭代次数
};

int test_simpson(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const TrapTestCase *tc = &cases[i];

        // 构造函数
        Simpson Simpson = NULL;
        Simpson_Err err = SI.simpson_create(tc->f, tc->a, tc->b, tc->max_iter, &Simpson, tc->name);
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

        // 计算积分
        double integral = 0.0;
        err = SI.simpson_integration(&Simpson, &integral);
        if (err != SIMPSON_OK) {
            if (tc->expect_error && err == tc->expected_err && err == SIMPSON_ERR_INVAL) {
                printf("[TEST] %-25s 积分计算参数错误: err=%d PASS\n", tc->name, err);
                passed++;
                SI.simpson_destroy(&Simpson);
                continue;
            } else {
                printf("[TEST] %-25s 积分计算失败: err=%d FAIL\n", tc->name, err);
                failed++;
                SI.simpson_destroy(&Simpson);
                continue;
            }
        }

        // 检查结果
        if (tc->expect_error) {
            printf("[TEST] %-25s 期望错误但计算成功: integral=%.10f FAIL\n", tc->name, integral);
            failed++;
        } else if (TEST_ABS_REL_CLOSE(integral, tc->expected_integral, 1e-8, 1e-8)) {
            printf("[TEST] %-25s 积分结果正确: integral=%.10f PASS\n", tc->name, integral);
            passed++;
        } else {
            printf("[TEST] %-25s 积分结果错误: integral=%.10f expected=%.10f FAIL\n", tc->name, integral, tc->expected_integral);
            failed++;
        }
        // 析构函数
        SI.simpson_destroy(&Simpson);
    }
    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return failed == 0 ? 0 : 1;
}


int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_simpson();
}