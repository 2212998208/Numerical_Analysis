#ifndef NUMERICAL_ANALYSIS_TRAPEZOIDAL_H
#define NUMERICAL_ANALYSIS_TRAPEZOIDAL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif



typedef enum Trapezoidal_Err {
    TRAP_OK = 0,
    TRAP_ERR_NOMEM = 1,
    TRAP_ERR_INVAL = 2,
    TRAP_ERR_MAXITER = 3,
    TRAP_ERR_DIVIDE_BY_ZERO = 4,
}Trapezoidal_Err;


typedef struct IntegrationApproximation *Trapezoidal;









typedef struct TrapezoidalIntegration {
    Trapezoidal_Err (*trapezoidal_integration)(const Trapezoidal *inTrap, double *outApproxIntegral);
    Trapezoidal_Err (*trapezoidal_create)(double (*f)(double x), double a, double b, size_t max_iter, Trapezoidal *outTrap, const char *name);
    Trapezoidal_Err (*trapezoidal_destroy)(Trapezoidal *inTrap);
}TrapezoidalIntegrationAPI;



extern const TrapezoidalIntegrationAPI TI;























#ifdef __cplusplus
}
#endif







#endif //NUMERICAL_ANALYSIS_TRAPEZOIDAL_H