#ifndef NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H
#define NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct NonlinearSA *NonlinearSA;

typedef enum Succesive_Err {
    SA_OK = 0,
    SA_ERR_NOMEM = 1,
    SA_ERR_INVAL = 2,
    SA_ERR_MAXITER = 3,
    SA_ERR_DIVIDE_BY_ZERO = 4,
    SA_ERR_NOAPPROXIMATION = 5,
}Succesive_Err;


typedef struct Succesive_Approximation {
    Succesive_Err (*nonlinear_sa_solve)(const NonlinearSA *outSA, double *outRoot);
    Succesive_Err (*NonlinearSA_create)(double (*g)(double x), const double x0,
        const double tol, const size_t max_iter,
        NonlinearSA *outSA,
        const char *name);
    Succesive_Err (*NonlinearSA_destroy)(const NonlinearSA *inSA);
}SAAPI;

extern const SAAPI SA;
#endif //NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H