#include "lagrange.h"
#include "test_lagrange.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv) {
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

    return 0;
}
