#include <stdio.h>
#include <locale.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "integrator.h"
#include "test_integrator.h"

// 控制台 UTF-8 设置
static void set_console_utf8(void) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");
#else
    setlocale(LC_ALL, "C.UTF-8");
#endif
}

// 示例函数
static double fx2(double x, void *u) { (void)u; return x * x; }

// 本地绝对值 (避免使用 fabs)
static inline double na_abs(double v) { return v < 0 ? -v : v; }

int main(int argc, char *argv[]) {
    (void)argc; (void)argv;
    set_console_utf8();

    puts("=== 自适应 RK4 数值积分示例 ===");

    double a = 0.0, b = 2.0;
    double fixed_result = Integrator.rk4_fixed(fx2, NULL, a, b, 4000);
    IntegratorStatus st;
    double adaptive_result = Integrator.rk4_adaptive(fx2, NULL, a, b,
        (AdaptiveConfig){1e-12, 1e-12, 28}, &st);
    double exact = 8.0 / 3.0;

    printf("固定步长结果: %.15f\n", fixed_result);
    printf("自适应结果  : %.15f (status=%d)\n", adaptive_result, st);
    printf("解析解      : %.15f\n", exact);
    printf("相对误差(自适应): %.3e\n",
           na_abs(adaptive_result - exact) / exact);

    return 0;
}