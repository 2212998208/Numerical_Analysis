#include <stdio.h>
#include <windows.h>

#include <newton.h>
#include <lagrange.h>
#include <test_newton.h>

// 本地工具
static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
    (NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))


// 查询点结果期望
typedef struct {
    double x;
    double expected; // 若 expect_error=1 则忽略
} Query;

// 用例描述
typedef struct {
    const char *name;
    Point *points;          // 数据点集合
    size_t point_count;     // 点数
    Query *queries;         // 查询集合
    size_t query_count;     // 查询数量
    double abs_tol;         // 绝对误差容限
    double rel_tol;         // 相对误差容限
    int expect_error;       // 是否期望插值时出现错误
    Newton_Err expected_err; // 期望的错误码
} NewtonTestCase;

// ------------------------- 测试数据定义 -------------------------
// 线性: y = 2x + 1 (点: (1,3),(3,7))
static Point linear_points[] = { {1,3}, {3,7} };
static Query linear_queries[] = { {2.0,5.0}, {1.0,3.0} };

// 二次: y = x^2 (点: (1,1),(2,4),(3,9))
static Point quadratic_points[] = { {1,1},{2,4},{3,9} };
static Query quadratic_queries[] = { {2.5,6.25}, {3.0,9.0}, {7.0, 49.0}};

// 常数: y = 10 (点: (5,10)) 任意x -> 10 (算法当前实现返回单点y)
static Point constant_points[] = { {5,10} };
static Query constant_queries[] = { {0.0,10.0}, {5.0,10.0} };

// 重复X: 触发除零错误 (点: (1,2),(2,5),(1,8))
static Point duplicate_points[] = { {1,2}, {2,5}, {1,8} };
static Query duplicate_queries[] = { {1.5, 0.0} }; // expected 忽略

// 乱序节点 (与 quadratic 相同集合不同顺序)
static Point unordered_points[] = { {3,9}, {1,1}, {2,4} };
static Query unordered_queries[] = { {2.5,6.25}, {2.0,4.0}, {5.0, 25.0}, {7.0, 49.0}};

static NewtonTestCase cases[] = {
    {"线性插值",      linear_points,    2, linear_queries,    2, 1e-9, 1e-9, 0, NEWTON_OK},
    {"二次插值",      quadratic_points, 3, quadratic_queries, 3, 1e-9, 1e-9, 0, NEWTON_OK},
    {"常数插值",      constant_points,  1, constant_queries,  2, 1e-9, 1e-9, 0, NEWTON_OK},
    {"重复X值",       duplicate_points, 3, duplicate_queries, 1, 1e-9, 1e-9, 1, NEWTON_ERR_DIVBYZERO},
    {"乱序节点",      unordered_points, 3, unordered_queries, 4, 1e-9, 1e-9, 0, NEWTON_OK},
};

// ------------------------- 测试执行 -------------------------
int test_newton(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int passed = 0;

    for (int i = 0; i < N; i++) {
        NewtonTestCase *tc = &cases[i];

        // 创建数据集
        NewtonDataSet *dataset = NULL;
        DataSet *inDataset = NULL;
        Newton_Err newton_err = Newton.create_newton_dataset(&dataset, tc->point_count);
        Lagrange_Err lagrange_err = Lagrange.create_dataset(&inDataset, tc->points, tc->point_count);
        if (newton_err != NEWTON_OK || lagrange_err != LAGRANGE_OK) {
            printf("[TEST] %-8s 构建数据集失败: err=%d FAIL\n", tc->name, newton_err);
            continue;
        }


        // 计算Newton插值
        int case_ok = 1;
        if (tc->expect_error) {
            // 触发一次插值以验证错误码
            double v = 0.0;
            Newton_Err err= Newton.newton_interpolate(&dataset,&inDataset,tc->queries[0].x, &v);
            if (err != tc->expected_err) {
                printf("[TEST] %-8s 期望错误码 %d 实际 %d FAIL\n", tc->name, tc->expected_err, err);
                case_ok = 0;
                continue;
            }
            printf("[TEST] %-8s 期望错误码 %d 实际 %d PASS\n", tc->name, tc->expected_err, err);
            continue;
        }
        for (size_t q = 0; q < tc->query_count; ++q) {
            double outY = 0.0;
            Newton_Err err = Newton.newton_interpolate(&dataset, &inDataset, tc->queries[q].x, &outY);
            if (err != NEWTON_OK) {
                printf("[TEST] %-8s 插值失败 err=%d FAIL\n", tc->name, err);
                case_ok = 0;
                break;
            }
            if (!TEST_ABS_REL_CLOSE(outY, tc->queries[q].expected, tc->abs_tol, tc->rel_tol)) {
                printf("[TEST] %-8s x=%.6f 期望 %.6f 实际 %.6f FAIL\n",
                       tc->name, tc->queries[q].x, tc->queries[q].expected, outY);
                case_ok = 0;
            } else {
                printf("[TEST] %-8s x=%.6f 期望 %.6f 实际 %.6f PASS\n",
                       tc->name, tc->queries[q].x, tc->queries[q].expected, outY);
            }
        }
        if (case_ok) {
            ++passed;
        }
        Newton.print_newton_dataset(&dataset);
        Newton.destroy_dataset(&dataset);
        if (!inDataset) {
        }else {
            Lagrange.destroy_dataset(inDataset);
        }
    }
    ++passed;

    printf("测试通过用例: %d / %d\n", passed, N);
    return (passed == N) ? 0 : 1;
}

int main(void) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_newton();
}