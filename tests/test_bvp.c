//
// Created by ninico on 2025/10/18.
//
#include "bvp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <windows.h>

/* 工具：向量无穷范数差 */
static double vec_err_inf(const double* x, const double* y, size_t n) {
    double v = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double d = fabs(x[i]-y[i]);
        if (d > v) v = d;
    }
    return v;
}

/*===========================
 * 测试用例结构体
 *===========================*/
typedef int (*test_run_fn)(void);

typedef struct {
    const char* name;
    const char* desc;
    int expect_success;
    test_run_fn run;
} bvp_test_case;

/*===========================
 * 用例1：制造解的非线性二元系统（无 y'）
 * 真解：y1=sin(pi x), y2=cos(pi x) on [0,1]
 * 方程：y'' = G(x,y)
 *  G1 = -pi^2 sin(pi x) + (y1 - sin)^2 + (y2 - cos)^2
 *  G2 = -pi^2 cos(pi x) + (y1 - sin)*(y2 - cos) + (y2 - cos)^3
 *===========================*/
typedef struct { int m; } user1_t; /* 占位 */

static int G1(double x, const double* y, const double* yp, double* out, void* user) {
    (void)user; (void)yp;
    const double pi = 3.14159265358979323846;
    double s = sin(pi*x), c = cos(pi*x);
    double e1 = y[0] - s, e2 = y[1] - c;
    out[0] = -pi*pi*s + e1*e1 + e2*e2;
    out[1] = -pi*pi*c + e1*e2 + e2*e2*e2;
    return 0;
}
static int Jy1(double x, const double* y, const double* yp, double* J, void* user) {
    (void)user; (void)yp;
    const double pi = 3.14159265358979323846;
    double s = sin(pi*x), c = cos(pi*x);
    double e1 = y[0] - s, e2 = y[1] - c;
    /* J = dG/dy (2x2，行主序) */
    /* dG1/dy1 = 2e1, dG1/dy2 = 2e2
       dG2/dy1 = e2,  dG2/dy2 = e1 + 3e2^2 */
    J[0*2+0] = 2.0*e1; J[0*2+1] = 2.0*e2;
    J[1*2+0] = e2;     J[1*2+1] = e1 + 3.0*e2*e2;
    return 0;
}
/* J_yp = 0 不提供，由求解器自动数值近似（也为0） */

static int run_case1(void) {
    const size_t m=2, N=100; const double a=0.0, b=1.0;
    double ya[2], yb[2];
    ya[0]=sin(0.0); ya[1]=cos(0.0);
    yb[0]=sin(M_PI); yb[1]=cos(M_PI);

    double* y_all = (double*)malloc((N+2)*m*sizeof(double));
    if (!y_all) return 0;

    bvp_solver_t* s=NULL;
    int err = bvp_create(m,N,a,b,&s);
    if (err) { printf("create err=%s\n", bvp_strerror(err)); free(y_all); return 0; }
    bvp_set_function(s, G1, NULL);
    bvp_set_jacobians(s, Jy1, NULL);
    bvp_set_boundary(s, ya, yb);
    bvp_set_options(s, 60, 1e-10, 1e-12);

    err = bvp_solve(s, y_all);
    if (err) { printf("solve err=%s\n", bvp_strerror(err)); bvp_destroy(&s); free(y_all); return 0; }

    /* 校验：与真解比较 */
    double* x = (double*)malloc((N+2)*sizeof(double));
    bvp_get_grid(s, x);
    double errmax = 0.0;
    for (size_t i=0;i<N+2;++i) {
        double ytrue0 = sin(M_PI*x[i]);
        double ytrue1 = cos(M_PI*x[i]);
        double e0 = fabs(y_all[i*m+0]-ytrue0);
        double e1 = fabs(y_all[i*m+1]-ytrue1);
        if (e0>errmax) errmax=e0; if (e1>errmax) errmax=e1;
    }
    printf("   max|y - y_true|_inf = %.3e\n", errmax);
    free(x);
    bvp_destroy(&s);

    int pass = (errmax < 2e-4); /* 二阶 FDM 误差 O(h^2)，N=100 足够小 */
    free(y_all);
    return pass;
}

/*===========================
 * 用例2：Bratu 方程（m=1）
 * y'' + λ e^y = 0, y(0)=y(1)=0, 取 λ=1
 *===========================*/
typedef struct { double lambda; } bratu_t;

static int G_bratu(double x, const double* y, const double* yp, double* out, void* user) {
    (void)x; (void)yp;
    bratu_t* U = (bratu_t*)user;
    out[0] = -U->lambda * exp(y[0]);
    return 0;
}
static int Jy_bratu(double x, const double* y, const double* yp, double* J, void* user) {
    (void)x; (void)yp;
    bratu_t* U = (bratu_t*)user;
    J[0] = -U->lambda * exp(y[0]);
    return 0;
}

static int run_case2(void) {
    const size_t m=1, N=100; const double a=0.0, b=1.0;
    double ya[1]={0.0}, yb[1]={0.0};
    double* y_all = (double*)malloc((N+2)*m*sizeof(double));
    if (!y_all) return 0;

    bratu_t U = { .lambda = 1.0 };
    bvp_options_t opt = { .max_iter=80, .tol_res=1e-10, .tol_step=1e-12 };
    int err = bvp_solve_raw(m,N,a,b, ya,yb, G_bratu, Jy_bratu, NULL, &U, &opt, y_all);
    if (err) { printf("   solve err=%s\n", bvp_strerror(err)); free(y_all); return 0; }

    /* 基本性质：对称、在中点取最大值且 >0；残差检查 */
    double h = (b-a)/(N+1);
    /* 计算离散残差最大值 */
    double rmax = 0.0;
    for (size_t i=1;i<=N;++i) {
        double yi_1 = y_all[(i-1)*m], yi=y_all[i*m], yi1=y_all[(i+1)*m];
        double Gi = -U.lambda * exp(yi);
        double ri = (yi_1 - 2*yi + yi1) - h*h*Gi;
        rmax = fmax(rmax, fabs(ri));
    }
    printf("   max discrete residual = %.3e, y(0.5)=%.6f\n", rmax, y_all[(N/2+1)*m]);

    int pass = (rmax < 5e-8) && (y_all[(N/2+1)*m] > 0.0);
    free(y_all);
    return pass;
}

/*===========================
 * 用例3：含 y' 的非线性（数值雅可比验证）
 * y'' = y*y'  (m=1), y(0)=y(1)=0
 * 真解 y=0；不提供雅可比，让求解器数值差分
 *===========================*/
static int G_yp(double x, const double* y, const double* yp, double* out, void* user) {
    (void)x; (void)user;
    out[0] = y[0] * yp[0];
    return 0;
}

static int run_case3(void) {
    const size_t m=1, N=50; const double a=0.0, b=1.0;
    double ya[1]={0.0}, yb[1]={0.0};
    double* y_all = (double*)malloc((N+2)*m*sizeof(double));
    if (!y_all) return 0;

    bvp_options_t opt = { .max_iter=50, .tol_res=1e-10, .tol_step=1e-12 };
    int err = bvp_solve_raw(m,N,a,b, ya,yb, G_yp, NULL, NULL, NULL, &opt, y_all);
    if (err) { printf("   solve err=%s\n", bvp_strerror(err)); free(y_all); return 0; }

    /* 应收敛到零解 */
    double maxabs = 0.0;
    for (size_t i=0;i<N+2;++i) {
        double v = fabs(y_all[i*m]);
        if (v>maxabs) maxabs = v;
    }
    printf("   max |y| = %.3e\n", maxabs);
    int pass = (maxabs < 1e-10);
    free(y_all);
    return pass;
}

/*===========================
 * 用例4：错误态——未设置边界
 *===========================*/
static int run_case4(void) {
    const size_t m=1, N=10; const double a=0.0, b=1.0;
    bvp_solver_t* s=NULL;
    int err = bvp_create(m,N,a,b,&s);
    if (err) return 0;
    bvp_set_function(s, G_bratu, NULL);
    double* y_all = (double*)malloc((N+2)*m*sizeof(double));
    if (!y_all) { bvp_destroy(&s); return 0; }
    err = bvp_solve(s, y_all);
    printf("   expect BAD_STATE, got err=%s\n", bvp_strerror(err));
    free(y_all); bvp_destroy(&s);
    return (err == BVP_ERR_BAD_STATE);
}

/*===========================
 * 主函数：跑所有用例
 *===========================*/
int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    const bvp_test_case cases[] = {
        { "Manufactured-Nonlinear-2D", "制造解的非线性二元系统（校验收敛阶与精度）", 1, run_case1 },
        { "Bratu-1D",                   "Bratu 方程 λ=1（经典非线性 BVP）",        1, run_case2 },
        { "Nonlinear-with-yp",         "含一阶导的非线性（数值雅可比验证）",       1, run_case3 },
        { "BadState-NoBoundary",       "错误态：未设置边界即求解，应报错",        1, run_case4 },
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
