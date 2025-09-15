#include <stdio.h>
#include <math.h>
#include <windows.h>

// 本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }

#include "integrator.h"
#include "test_integrator.h"

// 替换原 TEST_ABS_REL_CLOSE
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
    (NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))

typedef struct {
    const char *name;
    double (*func)(double, void *);
    double a, b;
    double expected;
} TestCase;

static double f_x2(double x, void *u)  { (void)u; return x * x; }
static double f_sin(double x, void *u) { (void)u; return sin(x); }
static double f_exp(double x, void *u) { (void)u; return exp(x); }

int run_all_tests(void) {
    TestCase cases[] = {
        {"x^2",  f_x2,  0.0, 2.0, (8.0/3.0)},
        {"sin",  f_sin, 0.0, M_PI, 2.0},
        {"exp",  f_exp, 0.0, 1.0, (M_E - 1.0)}
    };
    const int N = (int)(sizeof(cases)/sizeof(cases[0]));
    int passed = 0;

    for (int i = 0; i < N; ++i) {
        const TestCase *tc = &cases[i];
        double fixed_val = Integrator.rk4_fixed(tc->func, NULL, tc->a, tc->b, 2000);
        IntegratorStatus st;
        double adapt_val = Integrator.rk4_adaptive(tc->func, NULL, tc->a, tc->b,
            (AdaptiveConfig){1e-10, 1e-10, 24}, &st);

        int ok_fixed = TEST_ABS_REL_CLOSE(fixed_val, tc->expected, 1e-9, 1e-9);
        int ok_adapt = TEST_ABS_REL_CLOSE(adapt_val, tc->expected, 5e-10, 5e-10) && st == INTEGRATOR_OK;

        printf("[TEST] %-6s fixed=%.12f adaptive=%.12f ref=%.12f  %s/%s\n",
               tc->name, fixed_val, adapt_val, tc->expected,
               ok_fixed ? "OK" : "FAIL",
               ok_adapt ? "OK" : "FAIL");

        if (ok_fixed && ok_adapt) ++passed;
    }
    printf("测试通过: %d / %d\n", passed, N);
    return (passed == N) ? 0 : 1;
}

int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    // 运行测试
    int test_ret = run_all_tests();
    if (test_ret != 0) {
        puts("部分测试失败，仍继续示例...");
    }
}