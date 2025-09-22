

#ifndef NUMERICAL_ANALYSIS_SIMPSON_H
#define NUMERICAL_ANALYSIS_SIMPSON_H
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


typedef struct IntegrationApproximation *Simpson;





typedef struct Simpson {
    Simpson_Err (*simpson_create)(double (*f)(double x), double a, double b, size_t max_iter, Simpson *outSimpson, const char *name);
    Simpson_Err (*simpson_integration)(const Simpson *inSimpson, double *outApproxIntegral);
    Simpson_Err (*simpson_destroy)(Simpson *inSimpson);
}SimpsonAPI;

extern const SimpsonAPI SI;


#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_SIMPSON_H