
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

#include "hermite.h"


static inline double NA_ABS(double v) { return v < 0 ? -v : v; }
static inline double NA_MAX(double a, double b) { return (a > b) ? a : b; }
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
(NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))

// 查询点结果期望
typedef struct {
    double x;
    double y_expected; // 若 expect_error=1 则忽略
    double dy_expected; // 若 expect_error=1 则忽略
} Query;

// 数据点
typedef struct {
    double x;
    double y;
    double dy; // 导数
} Point;

// 用例描述
typedef struct {
    const char *name;
    Point *points;   // 数据点集合
    size_t point_count;     // 点数
    Query *queries;         // 查询集合
    size_t query_count;     // 查询数量
    double abs_tol;         // 绝对误差容限
    double rel_tol;         // 相对误差容限
    int expect_error;       // 是否期望插值时出现错误
    Hermite_Err expected_err; // 期望的错误码
} HermiteTestCase;

// ------------------------- 测试数据定义 -------------------------
// 线性: y = 2x + 1 (点: (1,3),(3,7), dy=2)
static Point linear_points[] = { {1,3,2}, {3,7,2} };
static Query linear_queries[] = { {2.0,5.0, 2.0}, {1.0,3.0, 2.0} };

// 二次：y=x^2 (点: (1,1),(2,4),(3,9), dy=2x)
static Point quadratic_points[] = { {1,1,2}, {2,4,4}, {3,9,6} };
static Query quadratic_queries[] = { {2.5,6.25,5.0}, {3.0,9.0,6.0}, {7.0, 49.0,14.0}};

static HermiteTestCase cases[] = {
    {"线性插值", linear_points, 2, linear_queries, 2, 1e-9, 1e-9, 0, HERMITE_OK},
    {"二次插值", quadratic_points, 3, quadratic_queries, 3, 1e-9, 1e-9, 0, HERMITE_OK}
};


// ------------------------- 测试执行 -------------------------
int test_hermite(void) {
    const int N = (int)(sizeof(cases) / sizeof(cases[0]));
    int failed = 0;
    int passed = 0;

    for (int i = 0; i < N; ++i) {
        HermiteTestCase *tc = &cases[i];

        // 创建数据集
        HermiteDataset dataset = NULL;
        dataset = Hermite.create_hermite_dataset(tc->point_count,
                                                        (const double[]){tc->points[0].x, tc->points[1].x},
                                                        (const double[]){tc->points[0].y, tc->points[1].y},
                                                        (const double[]){tc->points[0].dy, tc->points[1].dy});
        if (!dataset) {
            printf("[TEST] %-8s 构建数据集失败: err=%d FAIL\n", tc->name, HERMITE_ERR_INVALID);
            failed++;
        }

        // 创建插值器
        HermiteInterpolator interpolator = NULL;
        Hermite_Err hermite_err = Hermite.hermite_create_interpolator(&dataset, &interpolator);
        if (hermite_err != HERMITE_OK) {
            if (tc->expect_error && hermite_err == tc->expected_err) {
                printf("[TEST] %-8s 创建插值器期望失败: err=%d PASS\n", tc->name, hermite_err);
                passed++;
            } else {
                printf("[TEST] %-8s 创建插值器失败: err=%d FAIL\n", tc->name, hermite_err);
                failed++;
            }
            Hermite.destroy_hermite_dataset(&dataset);
        }

        // 计算插值
        for (size_t q = 0; q < tc->query_count; ++q) {
            const Query *query = &(tc->queries[q]);
            double y_out = 0.0;
            double dy_out = 0.0;
            hermite_err = Hermite.hermite_evaluate(&interpolator, query->x, &y_out, &dy_out);
            if (hermite_err != HERMITE_OK) {
                if (tc->expect_error && hermite_err == tc->expected_err) {
                    printf("[TEST] %-8s 查询(%.3f)期望失败: err=%d PASS\n", tc->name, query->x, hermite_err);
                    passed++;
                } else {
                    printf("[TEST] %-8s 查询(%.3f)失败: err=%d FAIL\n", tc->name, query->x, hermite_err);
                    failed++;
                }
                Hermite.hermite_destroy_interpolator(&interpolator);
                Hermite.destroy_hermite_dataset(&dataset);
                continue;
            }

            // 验证结果
            if (TEST_ABS_REL_CLOSE(y_out, query->y_expected, tc->abs_tol, tc->rel_tol) &&
                TEST_ABS_REL_CLOSE(dy_out, query->dy_expected, tc->abs_tol, tc->rel_tol)) {
                printf("[TEST] %-8s x=%.3f 结果正确: y=%.6f dy=%.6f PASS\n", tc->name, query->x, y_out, dy_out);
                passed++;
            } else {
                printf("[TEST] %-8s x=%.3f 结果错误: got y=%.6f dy=%.6f, expected y=%.6f dy=%.6f FAIL\n",
                       tc->name, query->x, y_out, dy_out, query->y_expected, query->dy_expected);
                failed++;
            }
        }

        Hermite.hermite_destroy_interpolator(&interpolator);
        Hermite.destroy_hermite_dataset(&dataset);
    }
    printf("[TEST] 通过 %d / %d 个用例\n", passed, passed + failed);
    return (passed == (passed+failed)) ? 0 : 1;
}

int main(int argc, char *argv[]) {
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    return test_hermite();
}

