#ifndef NUMERICAL_ANALYSIS_BISECTION_H
#define NUMERICAL_ANALYSIS_BISECTION_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BisectionRange *BisectionRange;



typedef enum Bisection_Err {
    BISECTION_OK = 0,
    BISECTION_ERR_NOMEM = 1,
    BISECTION_ERR_INVALID = 2,
    BISECTION_ERR_MAXITER = 3
} Bisection_Err;

typedef struct Bisection {
    Bisection_Err (*bisection_solve)(BisectionRange *outRange);
    Bisection_Err (*bisection_create)(double (*f)(double x), double a,
        double b, double tol,
        BisectionRange *outRange,
        const char *name);
    Bisection_Err (*bisection_destroy)(BisectionRange *outRange);
    double (*bisection_get_midpoint)(const BisectionRange *outRange);
}BisectionAPI;


extern const BisectionAPI Bisection;

#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_BISECTION_H