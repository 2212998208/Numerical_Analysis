#ifndef NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H
#define NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef enum Simpson_Err {
    SIMPSON_OK = 0,
    SIMPSON_ERR_NOMEM = 1,
    SIMPSON_ERR_INVAL = 2,
    SIMPSON_ERR_MAXITER = 3,
    SIMPSON_ERR_DIVIDE_BY_ZERO = 4,
}Simpson_Err;

typedef struct DoubleIntegrationApproximation *DoubleSimpson;


typedef struct DoubleSimpson {
    Simpson_Err (*double_simpson_create)(double (*f)(double x, double y),
                                  const double x_a, const double x_b,
                                  const double y_c, const double y_d,
                                  const size_t n, const size_t m,
                                  DoubleSimpson *outDoubleSimpson,
                                  const char *name);
    Simpson_Err (*double_simpson_destroy)(DoubleSimpson *inDoubleSimpson);
    Simpson_Err (*double_simpson_integrate)(const DoubleSimpson *inDoubleSimpson, double *outApproxIntegral);
} DoubleSimpsonAPI;


extern const DoubleSimpsonAPI DS;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H