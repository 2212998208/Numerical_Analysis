#include <math.h>
#include <stdio.h>
#include <windows.h>

#include "successive_approximation.h"

// 误差计算本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))

// 函数定义
double gfunc1(double x) {return sqrt(2 * x + 3);} // x=3
double gfunc2(double x) {return -sqrt(2 * x + 3);} // x=-1
double gfunc3(double x) {return cos(x);} // x=0.7390851332151607
double gfunc4(double x) {return exp(-x);} // x=0.5671432904097838
double gfun5(double x) {return 0.5 * x * x - 0.5 * x + 1;} // x=1.0和x=2.0
double gfun6(double x) {return M_PI + 1.05 * sin(x);} // x= M_PI

// 用例描述
typedef struct {
    const char *name;
    double (*gfunc)(double x);
    double x0;
    double tol;
    size_t max_iter;
    double expected_root; // 若 expect_error=1 则忽略
    int expect_error;       // 是否期望插值时出现错误
    Succesive_Err expected_err; // 期望的错误码
}SATestCase;

static SATestCase cases [] = {
    {"g1(x)=sqrt(2x+3)", gfunc1, 0.0, 1e-8, 256, 3.0, 0, SA_OK},
    {"g2(x)=-sqrt(2x+3)", gfunc2, -2, 1e-8, 256, -1.0, 1, SA_ERR_NOAPPROXIMATION}, // 不动点迭代法不收敛
    {"g3(x)=cos(x)", gfunc3, 0.5, 1e-8, 256, 0.7390851332151607, 0, SA_OK},
    {"g4(x)=exp(-x)", gfunc4, 1.0, 1e-8, 512, 0.5671432904097838, 0, SA_OK},
    {"g4(x)=exp(-x)", gfunc4, -1.0, 1e-8, 1024, 0.5671432904097838, 1, SA_ERR_NOAPPROXIMATION}, // 不动点迭代法不收敛
    {"g4(x)=exp(-x)", gfunc4, 1.0, 1e-64, 2048, 0.5671432904097838, 0, SA_OK},
    {"g5(x)=0.5x^2-0.5x+1", gfun5, 0.0, 1e-8, 256, 1.0, 0, SA_OK}, // 收敛到x=1
    {"g5(x)=0.5x^2-0.5x+1", gfun5, 1.5, 1e-8, 256, 2.0, 1, SA_ERR_NOAPPROXIMATION}, // 不动点迭代法不收敛
    {"g5(x)=0.5x^2-0.5x+1", gfun5, -0.5, 1e-8, 256, 1.0, 0, SA_OK}, // 负区间收敛到x=1
    {"g5(x)=0.5x^2-0.5x+1", gfun5, 2.1, 1e-8, 256, 2.0, 1, SA_ERR_NOAPPROXIMATION}, // 不动点迭代法不收敛
    // {"g6(x)=π+sin(x)", gfun6, M_PI/1.5, 1e-8, 10248, M_PI, 1, SA_ERR_NOAPPROXIMATION}, // 过程不收敛
};


int test_SA(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;
    int failed = 0;
    for (int i = 0; i < N; ++i) {
        const SATestCase *tc = &cases[i];

        // 构造函数
        NonlinearSA sa = NULL;
        Succesive_Err err = SA.NonlinearSA_create(tc->gfunc,tc->x0,tc->tol,tc->max_iter,&sa,tc->name);
        if (err != SA_OK) {
            if (tc->expect_error && err == tc->expected_err && err == SA_ERR_NOAPPROXIMATION) {
                printf("[TEST] %-25s 初始点不满足收敛定理: err=%d PASS\n", tc->name, err);
                passed++;
                continue;
            }
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-25s 构建区间期望失败: err=%d PASS\n", tc->name, err);
                passed++;
                continue;
            }
            printf("[TEST] %-25s 构建区间失败: err=%d FAIL\n", tc->name, err);
            failed++;
            continue;
        }

        // 不动点迭代
        double root = 0.0;
        err = SA.nonlinear_sa_solve(&sa, &root);
        if (err != SA_OK) {
            if (tc->expect_error && err == tc->expected_err && err == SA_ERR_NOAPPROXIMATION) {
                printf("[TEST] %-25s 存在点不满足收敛定理: err=%d PASS\n", tc->name, err);
                passed++;
            }
            if (tc->expect_error && err == tc->expected_err) {
                printf("[TEST] %-25s 不动点迭代期望失败: err=%d PASS\n", tc->name, err);
                passed++;
            } else {
                printf("[TEST] %-25s 不动点迭代失败: err=%d FAIL\n", tc->name, err);
                failed++;
            }
        } else {
            if (tc->expect_error) {
                printf("[TEST] %-25s 不动点迭代期望失败但实际成功: root=%.15f FAIL\n", tc->name, root);
                failed++;
            } else if (TEST_ABS_REL_CLOSE(root, tc->expected_root, tc->tol, tc->tol)) {
                printf("[TEST] %-25s 不动点迭代成功: root=%.15f PASS\n", tc->name, root);
                passed++;
            } else {
                printf("[TEST] %-25s 不动点迭代结果错误: root=%.15f expected=%.15f FAIL\n", tc->name, root, tc->expected_root);
                failed++;

            }
        }

        // 析构函数
        err = SA.NonlinearSA_destroy(&sa);
        if (err != SA_OK) {
            printf("[TEST] %-25s 析构失败: err=%d FAIL\n", tc->name, err);
            failed++;
        }
    }
    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return (passed == (passed + failed)) ? 0 : 1;
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_SA();
}