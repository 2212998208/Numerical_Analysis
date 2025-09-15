#ifndef NUMERICAL_ANALYSIS_LAGRANGE_H
#define NUMERICAL_ANALYSIS_LAGRANGE_H
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

typedef struct DataSet DataSet;
typedef struct Point {
        double x;
        double y;
} Point;

typedef enum Lagrange_Err {
    LAGRANGE_OK = 0,
    LAGRANGE_ERR_NOMEM = 1,
    LAGRANGE_ERR_INVALID = 2,
    LAGRANGE_ERR_DIVBYZERO = 3
} Lagrange_Err;


// 便捷Point构造器
static inline Point point_make(const double x, const double y) {
    const Point p = {x, y};
    return p;
}


// API 结构(可扩展更多算法)
typedef struct {
    double (*lagrange_interpolate)(DataSet *dataset, double x);
    DataSet *(*empty_dataset)(void);
    Lagrange_Err (*create_dataset)(DataSet **dataset, Point *points, size_t size);
    Lagrange_Err (*create_dataset_from_points)(DataSet **dataset, size_t size);
    void (*destroy_dataset)(DataSet *dataset);
    Point *(*get_points)(DataSet **outDataset);
} LagrangeAPI;


extern const LagrangeAPI Lagrange;

#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_LAGRANGE_H