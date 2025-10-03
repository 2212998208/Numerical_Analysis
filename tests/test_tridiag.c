#include "tri.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <windows.h>

/* 判等：||x - y||_inf <= eps */
static int approx_equal_vec(const double* x, const double* y, size_t n, double eps) {
    double maxdiff = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = fabs(x[i] - y[i]);
        if (diff > maxdiff) maxdiff = diff;
        if (diff > eps) return 0;
    }
    (void)maxdiff; /* 可用于调试 */
    return 1;
}

/* 生成三对角矩阵-向量乘积 y = A x
   a[0..n-2], b[0..n-1], c[0..n-2] */
static void tri_matvec(size_t n,
                       const double* a,
                       const double* b,
                       const double* c,
                       const double* x,
                       double* y)
{
    if (n == 0) return;
    if (n == 1) { y[0] = b[0] * x[0]; return; }

    y[0] = b[0] * x[0] + c[0] * x[1];
    for (size_t i = 1; i + 1 < n; ++i) {
        y[i] = a[i - 1] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1];
    }
    y[n - 1] = a[n - 2] * x[n - 2] + b[n - 1] * x[n - 1];
}

/*===========================
 * 测试用例结构体
 *===========================*/
typedef struct {
    const char* name;        /* 用例名称 */
    const char* desc;        /* 用例描述 */
    int expect_success;      /* 是否期望成功 */
    size_t n;                /* 规模 */
    /* 生成 a,b,c,d 和期望解 x_true 的函数指针 */
    void (*build)(size_t n,
                  double* a, double* b, double* c,
                  double* d, double* x_true);
} tri_test_case;

/* 用例1：n=3，解析可解（Poisson 1D 的经典三对角：2,-1,-1） */
static void build_case_poisson_3(size_t n,
                                 double* a, double* b, double* c,
                                 double* d, double* x_true)
{
    (void)n;
    /* A = [[ 2 -1  0]
            [-1  2 -1]
            [ 0 -1  2]]，令真解 x_true = [1,2,3]^T */
    x_true[0] = 1.0; x_true[1] = 2.0; x_true[2] = 3.0;
    b[0] = 2.0; b[1] = 2.0; b[2] = 2.0;
    a[0] = -1.0; a[1] = -1.0;
    c[0] = -1.0; c[1] = -1.0;
    tri_matvec(3, a, b, c, x_true, d);
}

/* 用例2：n=1 的边界情况 */
static void build_case_n1(size_t n,
                          double* a, double* b, double* c,
                          double* d, double* x_true)
{
    (void)a; (void)c; (void)n;
    b[0] = 5.0;
    x_true[0] = 2.0;
    d[0] = b[0] * x_true[0];
}

/* 用例3：构造零主元（b0=0）导致失败 */
static void build_case_zero_pivot(size_t n,
                                  double* a, double* b, double* c,
                                  double* d, double* x_true)
{
    (void)x_true;
    /* n=3，b0=0, c0=1, 其余随意，d 任意 */
    b[0] = 0.0; b[1] = 2.0; b[2] = 2.0;
    a[0] = -1.0; a[1] = -1.0;
    c[0] = 1.0;  c[1] = -1.0;
    d[0] = 1.0;  d[1] = 0.0; d[2] = 0.0;
}

/* 用例4：大规模随机严格对角占优（保证稳健性），检验精度 */
static void build_case_random_diag_dom(size_t n,
                                       double* a, double* b, double* c,
                                       double* d, double* x_true)
{
    /* 固定随机种子，保证可复现 */
    srand(12345u);
    for (size_t i = 0; i < n; ++i) {
        x_true[i] = (double)rand() / RAND_MAX * 2.0 - 1.0; /* 真解 in [-1,1] */
    }
    if (n == 1) {
        a[0] = 0.0; c[0] = 0.0; /* 忽略 */
        b[0] = 3.0;
        tri_matvec(1, a, b, c, x_true, d);
        return;
    }
    for (size_t i = 0; i < n - 1; ++i) {
        a[i] = (double)rand() / RAND_MAX - 0.5;
        c[i] = (double)rand() / RAND_MAX - 0.5;
    }
    for (size_t i = 0; i < n; ++i) {
        /* 严格对角占优：b_i = |a_{i-1}| + |c_i| + 1 */
        double left = (i > 0)     ? fabs(a[i - 1]) : 0.0;
        double right= (i < n - 1) ? fabs(c[i])     : 0.0;
        b[i] = left + right + 1.0;
    }
    tri_matvec(n, a, b, c, x_true, d);
}

static void run_one_case(const tri_test_case* tc, double atol) {
    printf("=== [%s] %s ===\n", tc->name, tc->desc);

    size_t n = tc->n;
    double *a=NULL, *b=NULL, *c=NULL, *d=NULL, *x_true=NULL, *x_raw=NULL, *x_obj=NULL;
    if (n > 1) {
        a = (double*)malloc((n - 1) * sizeof(double));
        c = (double*)malloc((n - 1) * sizeof(double));
    } else {
        a = (double*)malloc(sizeof(double));
        c = (double*)malloc(sizeof(double));
    }
    b = (double*)malloc(n * sizeof(double));
    d = (double*)malloc(n * sizeof(double));
    x_true = (double*)malloc(n * sizeof(double));
    x_raw  = (double*)malloc(n * sizeof(double));
    x_obj  = (double*)malloc(n * sizeof(double));

    if (!a || !b || !c || !d || !x_true || !x_raw || !x_obj) {
        fprintf(stderr, "malloc failed\n");
        exit(1);
    }

    tc->build(n, a, b, c, d, x_true);

    /* 1) 直接数组 API */
    int err_raw = tri_solve_raw(n, a, b, c, d, x_raw, 1e-12);

    /* 2) 对象式 API */
    tri_solver_t* solver = NULL;
    int err_obj = tri_create(n, &solver);
    if (err_obj == TRI_SUCCESS) err_obj = tri_set_coeffs(solver, a, b, c);
    if (err_obj == TRI_SUCCESS) err_obj = tri_set_rhs(solver, d);
    if (err_obj == TRI_SUCCESS) err_obj = tri_solve(solver, x_obj);

    /* 期望与断言 */
    if (tc->expect_success) {
        if (err_raw != TRI_SUCCESS || err_obj != TRI_SUCCESS) {
            printf("  -> FAILED (unexpected error). raw=%s, obj=%s\n",
                    tri_strerror(err_raw), tri_strerror(err_obj));
        } else {
            int ok1 = approx_equal_vec(x_raw, x_true, n, atol);
            int ok2 = approx_equal_vec(x_obj, x_true, n, atol);
            if (ok1 && ok2) {
                printf("  -> PASSED (|x-x_true|_inf <= %.1e)\n", atol);
            } else {
                printf("  -> FAILED (solution mismatch). max |x-x_true|_inf > %.1e\n", atol);
            }
        }
    } else { /* 预期失败 */
        if (err_raw == TRI_SUCCESS || err_obj == TRI_SUCCESS) {
            printf("  -> FAILED (expected failure but solved)\n");
        } else {
            printf("  -> PASSED (caught error). raw=%s, obj=%s\n",
                    tri_strerror(err_raw), tri_strerror(err_obj));
        }
    }

    tri_destroy(&solver);
    free(a); free(b); free(c); free(d); free(x_true); free(x_raw); free(x_obj);
}

int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    const tri_test_case cases[] = {
        {
            .name = "Poisson-3",
            .desc = "经典离散泊松三对角（n=3），验证解析解",
            .expect_success = 1,
            .n = 3,
            .build = build_case_poisson_3
        },
        {
            .name = "N=1",
            .desc = "边界情形 n=1",
            .expect_success = 1,
            .n = 1,
            .build = build_case_n1
        },
        {
            .name = "ZeroPivot",
            .desc = "构造 b0=0 的零主元例子，应报错",
            .expect_success = 0,
            .n = 3,
            .build = build_case_zero_pivot
        },
        {
            .name = "RandDiagDom-1000",
            .desc = "n=1000 的严格对角占优随机例子，检验稳健性与精度",
            .expect_success = 1,
            .n = 1000,
            .build = build_case_random_diag_dom
        },
    };

    const size_t ncases = sizeof(cases) / sizeof(cases[0]);
    const double atol = 1e-9;

    for (size_t i = 0; i < ncases; ++i) {
        run_one_case(&cases[i], atol);
    }
    return 0;
}
