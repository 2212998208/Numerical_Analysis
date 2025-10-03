//
// Created by ninico on 2025/10/18.
//
#include "nls.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <windows.h>

/*===========================
 * 测试用例描述结构
 *===========================*/
typedef int (*test_run_fn)(void);
typedef struct {
    const char* name;
    const char* desc;
    int expect_success;
    test_run_fn run;
} nls_test_case;

/* 工具 */
static double vnorm_inf(const double* x, size_t n) {
    double v=0.0; for (size_t i=0;i<n;++i){ double a=fabs(x[i]); if(a>v)v=a;} return v;
}

/*==========================================================
 * 用例1：Rosenbrock 两元系统（解析雅可比）
 * F1 = 10(x2 - x1^2), F2 = 1 - x1，根 (1,1)
 *==========================================================*/
static int F_rosen(const double* x, double* f, void* user) {
    (void)user;
    f[0] = 10.0*(x[1] - x[0]*x[0]);
    f[1] = 1.0 - x[0];
    return 0;
}
static int J_rosen(const double* x, double* J, void* user) {
    (void)user;
    /* 行主序 2×2 */
    J[0*2+0] = -20.0*x[0];  J[0*2+1] = 10.0;
    J[1*2+0] = -1.0;        J[1*2+1] = 0.0;
    return 0;
}
static int run_case1(void) {
    const size_t n=2;
    double x[2] = {-1.2, 1.0}; /* 经典初值 */
    nls_options_t opt = { .max_iter=50, .tol_res=1e-12, .tol_step=1e-14,
                          .fd_eps=0, .ls_c=1e-4, .ls_beta=0.5, .ls_max_backtrack=20 };
    int err = nls_solve_raw(n, F_rosen, J_rosen, NULL, &opt, x);
    if (err) { printf("  err=%s\n", nls_strerror(err)); return 0; }
    printf("  sol=(%.10f, %.10f)\n", x[0], x[1]);
    return (fabs(x[0]-1.0) < 1e-10 && fabs(x[1]-1.0) < 1e-10);
}

/*==========================================================
 * 用例2：圆+直线（不提供雅可比，测试数值差分）
 * F1 = x + y - 1
 * F2 = x^2 + y^2 - 1
 * 期望根：(1,0) 或 (0,1)
 *==========================================================*/
typedef struct { int dummy; } ctx2_t;
static int F_circline(const double* x, double* f, void* user) {
    (void)user;
    f[0] = x[0] + x[1] - 1.0;
    f[1] = x[0]*x[0] + x[1]*x[1] - 1.0;
    return 0;
}
static int run_case2(void) {
    const size_t n=2;
    double x[2] = {0.5, 0.5};
    nls_options_t opt = { .max_iter=5000, .tol_res=1e-12, .tol_step=1e-12,
                          .fd_eps=0, .ls_c=1e-4, .ls_beta=0.5, .ls_max_backtrack=25 };
    int err = nls_solve_raw(n, F_circline, NULL, NULL, &opt, x);
    if (err) { printf("  err=%s\n", nls_strerror(err)); return 0; }
    printf("  sol≈(%.10f, %.10f)\n", x[0], x[1]);
    int ok1 = (fabs(x[0]-1.0) < 1e-9 && fabs(x[1]-0.0) < 1e-9);
    int ok2 = (fabs(x[0]-0.0) < 1e-9 && fabs(x[1]-1.0) < 1e-9);
    return ok1 || ok2;
}

/*==========================================================
 * 用例3：奇异雅可比触发错误
 * F = [x1^2, x1*x2]，在 (0,0) 处 J=0，期望返回 NLS_ERR_SINGULAR
 *==========================================================*/
static int F_sing(const double* x, double* f, void* user) {
    (void)user;
    f[0] = x[0]*x[0];
    f[1] = x[0]*x[1];
    return 0;
}
static int run_case3(void) {
    const size_t n=2;
    double x[2] = {0.0, 0.0};
    nls_options_t opt = { .max_iter=30, .tol_res=1e-14, .tol_step=1e-14,
                          .fd_eps=0, .ls_c=1e-4, .ls_beta=0.5, .ls_max_backtrack=10 };
    int err = nls_solve_raw(n, F_sing, NULL, NULL, &opt, x);
    printf("  expect SINGULAR, got err=%s\n", nls_strerror(err));
    return (err == NLS_ERR_SINGULAR);
}

/*==========================================================
 * 用例4：制造解（中等规模、轻耦合；不提供雅可比）
 * 设真解 x* 随机固定，令
 *   F_i(x) = sin(x_i) + γ * Σ_j x_j^3 - b_i
 *   其中 b_i = sin(x*_i) + γ * Σ_j x*_j^3
 * γ=0.05 以保证良性耦合与收敛
 *==========================================================*/
typedef struct {
    size_t n;
    double gamma;
    double* b;
} ctx4_t;

static int F_manu(const double* x, double* f, void* user) {
    ctx4_t* C = (ctx4_t*)user;
    size_t n = C->n;
    double S=0.0; for (size_t j=0;j<n;++j) S += x[j]*x[j]*x[j];
    for (size_t i=0;i<n;++i) f[i] = sin(x[i]) + C->gamma*S - C->b[i];
    return 0;
}

static int run_case4(void) {
    const size_t n=6;
    const double gamma = 0.05;
    double x_true[6] = { -0.9, 0.3, 1.1, -0.7, 0.2, 0.8 };
    double S=0.0; for(size_t j=0;j<n;++j) S += x_true[j]*x_true[j]*x_true[j];
    double* b = (double*)malloc(n*sizeof(double));
    for (size_t i=0;i<n;++i) b[i] = sin(x_true[i]) + gamma*S;
    ctx4_t C = { .n=n, .gamma=gamma, .b=b };

    double x[6];
    for (size_t i=0;i<n;++i) x[i] = x_true[i] + ((i%2)? 0.05 : -0.05); /* 轻微扰动作为初值 */

    nls_options_t opt = { .max_iter=80, .tol_res=1e-12, .tol_step=1e-12,
                          .fd_eps=0, .ls_c=1e-4, .ls_beta=0.5, .ls_max_backtrack=25 };
    int err = nls_solve_raw(n, F_manu, NULL, &C, &opt, x);
    if (err) { printf("  err=%s\n", nls_strerror(err)); free(b); return 0; }

    double err_inf=0.0; for(size_t i=0;i<n;++i){ double d=fabs(x[i]-x_true[i]); if(d>err_inf)err_inf=d; }
    printf("  max|x - x*|_inf = %.3e\n", err_inf);
    free(b);
    return (err_inf < 1e-9);
}

/*===========================
 * 主函数：运行全部用例
 *===========================*/
int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    const nls_test_case cases[] = {
        { "Rosenbrock-2D", "解析雅可比，根在 (1,1)", 1, run_case1 },
        { "Circle+Line",   "不提供雅可比，数值差分验证", 1, run_case2 },
        { "Singular-J",    "奇异雅可比触发错误", 1, run_case3 },
        { "Manufactured-6D", "制造解，轻耦合中等规模", 1, run_case4 },
    };
    const size_t ncases = sizeof(cases)/sizeof(cases[0]);

    int allpass = 1;
    for (size_t i = 0; i < ncases; ++i) {
        printf("=== [%s] %s ===\n", cases[i].name, cases[i].desc);
        int ok = cases[i].run();
        if (cases[i].expect_success && !ok) {
            printf("  -> FAILED\n");
        } else if (!cases[i].expect_success && ok) {
            printf("  -> FAILED (unexpected pass)\n");
        } else {
            printf("  -> PASSED\n");
        }
        allpass = allpass && ok;
    }
    printf("\nOVERALL: %s\n", allpass ? "PASS" : "FAIL");
    return allpass ? 0 : 1;
}
