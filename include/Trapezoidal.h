#ifndef NUMERICAL_ANALYSIS_TRAPEZOIDAL_H
#define NUMERICAL_ANALYSIS_TRAPEZOIDAL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif



/**
 * @brief Represents the error codes for the trapezoidal integration module.
 */
typedef enum Trapezoidal_Err {
    TRAP_OK = 0, /**< No error. */
    TRAP_ERR_NOMEM = 1, /**< Memory allocation failed. */
    TRAP_ERR_INVAL = 2, /**< Invalid input. */
    TRAP_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    TRAP_ERR_DIVIDE_BY_ZERO = 4, /**< Division by zero. */
}Trapezoidal_Err;


/**
 * @brief An opaque handle to an integration approximation.
 */
typedef struct IntegrationApproximation *Trapezoidal;

/**
 * @brief Defines the API for the trapezoidal integration module.
 */
typedef struct TrapezoidalIntegration {
    /**
     * @brief Performs trapezoidal integration.
     * @param inTrap A pointer to the integration approximation handle.
     * @param outApproxIntegral A pointer to a variable to store the approximate integral.
     * @return A Trapezoidal_Err error code.
     */
    Trapezoidal_Err (*trapezoidal_integration)(const Trapezoidal *inTrap, double *outApproxIntegral);

    /**
     * @brief Creates a trapezoidal integration handle.
     * @param f The function to integrate.
     * @param a The lower limit of integration.
     * @param b The upper limit of integration.
     * @param max_iter The maximum number of iterations.
     * @param outTrap A pointer to a pointer to the integration approximation handle to create.
     * @param name A name for the integration.
     * @return A Trapezoidal_Err error code.
     */
    Trapezoidal_Err (*trapezoidal_create)(double (*f)(double x), double a, double b, size_t max_iter, Trapezoidal *outTrap, const char *name);

    /**
     * @brief Destroys a trapezoidal integration handle.
     * @param inTrap A pointer to the integration approximation handle to destroy.
     * @return A Trapezoidal_Err error code.
     */
    Trapezoidal_Err (*trapezoidal_destroy)(Trapezoidal *inTrap);
}TrapezoidalIntegrationAPI;



/**
 * @brief The global instance of the trapezoidal integration API.
 */
extern const TrapezoidalIntegrationAPI TI;

#ifdef __cplusplus
}
#endif







#endif //NUMERICAL_ANALYSIS_TRAPEZOIDAL_H
