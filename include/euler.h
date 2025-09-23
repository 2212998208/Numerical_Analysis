#ifndef NUMERICAL_ANALYSIS_EULER_H
#define NUMERICAL_ANALYSIS_EULER_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct EulerDifferential *Euler;


typedef enum Euler_Err {
    EULER_OK = 0,
    EULER_INVALID = 1,
    EULER_DIVISION_BY_ZERO = 2,
}Euler_Err;

typedef struct Euler {
    Euler_Err (*euler_create) (Euler *outEuler, double (*dx2dt)(double x, double t), double x0, double t0, double Δt);
    Euler_Err (*euler_destroy) (Euler *euler);
    Euler_Err (*euler_solve) (const Euler *inEuler, size_t max_iter, double *xn);
}EulerAPI;











extern const EulerAPI EULER;



#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_EULER_H