#ifndef NUMERICAL_ANALYSIS_RK4_H
#define NUMERICAL_ANALYSIS_RK4_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stddef.h>


typedef struct RK4Data *rk;

typedef enum RK4 {
    RK4_OK = 0,
    RK4_INVALID = 1,
    RK4_POINTER_ERROR = 2,
}RK4_Err;


// 不采用API设计
RK4_Err rk4_create(rk *outRK, double (*f)(double x, double t),
    const double t0, const double x0, const double h, size_t max_iter);

RK4_Err rk4_destroy(rk *outRK);

RK4_Err rk4_solve(const rk *inRk, double *outY);


// 采用API设计
typedef struct RK4_API {
    RK4_Err (*create)(rk *outRK, double (*f)(double x, double t),
        const double t0, const double x0, const double h, size_t max_iter);
    RK4_Err (*destroy)(rk *outRK);
    RK4_Err (*solve)(const rk *inRk, double *outY);
}RK4_API;


extern const RK4_API RK4;
#ifdef __cplusplus
}
#endif


#endif //NUMERICAL_ANALYSIS_RK4_H