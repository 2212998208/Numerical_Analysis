#ifndef NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H
#define NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif






typedef enum NewtonRaphson_Err {
    NEWTON_RAPHSON_OK = 0,
    NEWTON_RAPHSON_ERR_NOMEM = 1,
    NEWTON_RAPHSON_ERR_INVALID = 2,
    NEWTON_RAPHSON_ERR_MAXITER = 3,
    NEWTON_RAPHSON_ERR_DERIVATIVE_ZERO = 4,
    NEWTON_RAPHSON_ERR_DERIVATIVE_UNSTABLE = 5
} NewtonRaphson_Err;




typedef struct NonLinearRange *NonLinearRange;


typedef struct NewtonRaphson {
    NewtonRaphson_Err (*NonLinearRange_create)(double (*f)(double x), const double x0,
        const double tol, const size_t max_iter,
        NonLinearRange *outRange,
        const char *name);
    NewtonRaphson_Err (*NonLinearRange_destroy)(const NonLinearRange *inRange);
    NewtonRaphson_Err (*newton_raphson_solve)(const NonLinearRange *outRange, double *outRoot);
}NewtonRaphsonAPI;


extern const NewtonRaphsonAPI NewtonRaphson;
#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H