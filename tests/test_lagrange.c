#include "lagrange.h"
#include "test_lagrange.h"
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

// 容差计算工具
static inline double NA_ABS(const double v) { return v < 0 ? -v : v; };
static inline double NA_MAX(const double a, const double b) { return (a > b) ? a : b; };


// 宏函数定义误差判断
#define TEST_ABS_REL_CLOSE(val, ref, abs_tol, rel_tol) \
    (NA_ABS((val) - (ref)) <= NA_MAX((abs_tol), NA_ABS(ref) * (rel_tol)))


// 测试结构体
typedef struct {
    DataSet **datasets;
    double *x;
    double *expected;
    double abs_tol;  // 绝对误差容限
    double rel_tol;  // 相对误差容限
    const int num_points;
    const char *name;
}TestCase;

// 创建一个线性插值检测（2个点的线性插值）{y=2x+1}
TestCase linear_case = {
    NULL,
    (double []) {2.0, 1.0},
     (double []) {5.0, 3.0},
    1e-9,
    1e-9,
    2,
    "线性插值"
};

// 创建一个二次插值检测（3个点的二次插值）{y=x^2}
TestCase quadratic_case = {
    NULL,
     (double []) {2.5, 3.0},
     (double []) {6.25, 9.0},
    1e-9,
    1e-9,
    2,
    "二次插值"
};

// 创建一个边界情况检测（1个点的常数插值）{y=10}
TestCase constant_case = {
    NULL,
     (double []) {0.0, 10.0},
     (double []) {10, 10},
    1e-9,
    1e-9,
    2,
    "常数插值"
};

// 创建错误处理-重复的X值检测
TestCase cubic_case = {
    NULL,
    (double []) {1.5},
    (double []) {LAGRANGE_ERR_DIVBYZERO},
    1e-9,
    1e-9,
    1,
    "重复X值"
};

// 创建乱序节点，Lagrange插值不受顺序影响
TestCase unordered_case = {
    NULL,
    (double []) {2.5, 3.0},
    (double []) {6.25, 9.0},
    1e-9,
    1e-9,
    2,
    "乱序节点"
};

int run_all_tests(void) {
    const int N1 = linear_case.num_points;
    const int N2 = quadratic_case.num_points;
    const int N3 = constant_case.num_points;
    const int N4 = cubic_case.num_points;
    const int N5 = unordered_case.num_points;

    Lagrange_Err calculated_err = LAGRANGE_OK;

    int passed = 0;

    Point linear_points[] = {
        point_make(1,3),
        point_make(3,7)
    };

    Point quadratic_points[] = {
        point_make(1,1),
        point_make(2,4),
        point_make(3,9),
    };

    Point constant_points[] = {
        point_make(5,10)
    };

    Point cubic_points[] = {
        point_make(1,2),
        point_make(2,5),
        point_make(1,8)
    };

    Point unordered_points[] = {
        point_make(3,9),
        point_make(1,1),
        point_make(2,4)
    };

    linear_case.datasets = Lagrange.empty_dataset();
    quadratic_case.datasets = Lagrange.empty_dataset();
    constant_case.datasets = Lagrange.empty_dataset();
    cubic_case.datasets = Lagrange.empty_dataset();
    unordered_case.datasets = Lagrange.empty_dataset();

    const Lagrange_Err err1 = Lagrange.create_dataset(linear_case.datasets, linear_points, 2);
    const Lagrange_Err err2 = Lagrange.create_dataset(quadratic_case.datasets, quadratic_points, 3);
    const Lagrange_Err err3 = Lagrange.create_dataset(constant_case.datasets, constant_points, 1);
    const Lagrange_Err err4 = Lagrange.create_dataset(cubic_case.datasets, cubic_points, 3);
    const Lagrange_Err err5 = Lagrange.create_dataset(unordered_case.datasets, unordered_points, 3);


    if (err1 != LAGRANGE_OK) {
        return -1;
    }
    if (err2 != LAGRANGE_OK) {
        return -1;
    }
    if (err3 != LAGRANGE_OK) {
        return -1;
    }
    if (err4 != LAGRANGE_OK) {
        return -1;
    }
    if (err5 != LAGRANGE_OK) {
        return -1;
    }

    printf("CREATE LAGRANGE TESTS OK\n");

    for (size_t i = 0; i < N1; ++i) {
        double Lagrange_calculated = Lagrange.lagrange_interpolate(*linear_case.datasets, linear_case.x[i]);
        int ok = TEST_ABS_REL_CLOSE(Lagrange_calculated, linear_case.expected[i], linear_case.abs_tol, linear_case.rel_tol);
        printf("[TEST] %-6s x=%.2f calculated=%.12f expected=%.12f  %s\n",
               linear_case.name, linear_case.x[i], Lagrange_calculated, linear_case.expected[i],
               ok ? "OK" : "FAIL");
    }
    ++passed;

    for (size_t i = 0; i < N2; ++i) {
        double Lagrange_calculated = Lagrange.lagrange_interpolate(*quadratic_case.datasets, quadratic_case.x[i]);
        int ok = TEST_ABS_REL_CLOSE(Lagrange_calculated, quadratic_case.expected[i], quadratic_case.abs_tol, quadratic_case.rel_tol);
        printf("[TEST] %-6s x=%.2f calculated=%.12f expected=%.12f  %s\n",
               quadratic_case.name, quadratic_case.x[i], Lagrange_calculated, quadratic_case.expected[i],
               ok ? "OK" : "FAIL");
    }
    ++passed;

    for (size_t i = 0; i < N3; ++i) {
        double Lagrange_calculated = Lagrange.lagrange_interpolate(*constant_case.datasets, constant_case.x[i]);
        int ok = TEST_ABS_REL_CLOSE(Lagrange_calculated, constant_case.expected[i], constant_case.abs_tol, constant_case.rel_tol);
        printf("[TEST] %-6s x=%.2f calculated=%.12f expected=%.12f  %s\n",
               constant_case.name, constant_case.x[i], Lagrange_calculated, constant_case.expected[i],
               ok ? "OK" : "FAIL");
    }
    ++passed;

    for (size_t i = 0; i < N4; ++i) {
        calculated_err = Lagrange.lagrange_interpolate(*cubic_case.datasets, cubic_case.x[i]);
    }

    if (calculated_err == LAGRANGE_ERR_DIVBYZERO) {
        printf("[TEST] %-6s  重复X值正确返回错误码 %d  OK\n", cubic_case.name, err4);
        ++passed;
    } else {
        printf("[TEST] %-6s  重复X值未正确返回错误码 %d  FAIL\n", cubic_case.name, err4);
    }

    for (size_t i = 0; i < N5; ++i) {
        double Lagrange_calculated = Lagrange.lagrange_interpolate(*unordered_case.datasets, unordered_case.x[i]);
        int ok = TEST_ABS_REL_CLOSE(Lagrange_calculated, unordered_case.expected[i], unordered_case.abs_tol, unordered_case.rel_tol);
        printf("[TEST] %-6s x=%.2f calculated=%.12f expected=%.12f  %s\n",
               unordered_case.name, unordered_case.x[i], Lagrange_calculated, unordered_case.expected[i],
               ok ? "OK" : "FAIL");
    }
    ++passed;


    // 清理
    Lagrange.destroy_dataset(*linear_case.datasets);
    Lagrange.destroy_dataset(*quadratic_case.datasets);
    Lagrange.destroy_dataset(*constant_case.datasets);
    Lagrange.destroy_dataset(*cubic_case.datasets);
    Lagrange.destroy_dataset(*unordered_case.datasets);
    free(linear_case.datasets);
    free(quadratic_case.datasets);
    free(constant_case.datasets);
    free(cubic_case.datasets);
    free(unordered_case.datasets);


    printf("PASSED\n");
    return (passed == 5) ? 0 : 1;

}

int main(const int argc, char **argv) {
    (void)argc; (void)argv;
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);

    /*
    (void)argc; (void)argv;


    DataSet **dataset = Lagrange.empty_dataset();
    DataSet **dataset2 = Lagrange.empty_dataset();
    Point points[] = {
        point_make(2,4),
        point_make(3,9),
        point_make(5,25)
    };

    const Lagrange_Err err = Lagrange.create_dataset_from_points(dataset, 3);
    const Lagrange_Err err2 = Lagrange.create_dataset(dataset2, points, 3);
    if (err != LAGRANGE_OK) {
        return -1;
    }
    if (err2 != LAGRANGE_OK) {
        return -1;
    }
    printf("LAGRANGE OK\n");


    // 测试插值
    const double x = 4.0;
    const double y = Lagrange.lagrange_interpolate(*dataset, x);
    const double y2 = Lagrange.lagrange_interpolate(*dataset2, x);
    printf("插值结果: f(%.2f) = %.2f\n", x, y);
    printf("插值结果: f(%.2f) = %.2f\n", x, y2);


    // 清理
    Lagrange.destroy_dataset(*dataset);
    Lagrange.destroy_dataset(*dataset2);
    free(dataset);
    free(dataset2);
    */

    run_all_tests();
    return 0;
}
