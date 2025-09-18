#ifndef NUMERICAL_ANALYSIS_SECANT_H
#define NUMERICAL_ANALYSIS_SECANT_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct NonLinearScant *NonLinearScant;


typedef enum Secant_Err {
    SECANT_OK = 0,
    SECANT_ERR_NOMEM = 1,
    SECANT_ERR_INVALID = 2,
    SECANT_ERR_MAXITER = 3,
    SECANT_ERR_DIVIDE_BY_ZERO = 4
} Secant_Err;

typedef struct Secant {
    Secant_Err (*secant_solve)(const NonLinearScant *outScant, double *outRoot);
    Secant_Err (*NonLinearScant_create)(double (*f)(double x), const double x0, const double x1,
        const double tol, const size_t max_iter,
        NonLinearScant *outScant,
        const char *name);
    Secant_Err (*NonLinearScant_destroy)(const NonLinearScant *inScant);
}SecantAPI;

extern const SecantAPI Secant;


#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_SECANT_H