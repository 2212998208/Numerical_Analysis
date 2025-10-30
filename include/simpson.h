
#ifndef NUMERICAL_ANALYSIS_SIMPSON_H
#define NUMERICAL_ANALYSIS_SIMPSON_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents the error codes for the Simpson's rule integration module.
 */
typedef enum Simpson_Err {
    SIMPSON_OK = 0, /**< No error. */
    SIMPSON_ERR_NOMEM = 1, /**< Memory allocation failed. */
    SIMPSON_ERR_INVAL = 2, /**< Invalid input. */
    SIMPSON_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    SIMPSON_ERR_DIVIDE_BY_ZERO = 4, /**< Division by zero. */
}Simpson_Err;


/**
 * @brief An opaque handle to an integration approximation.
 */
typedef struct IntegrationApproximation *Simpson;

/**
 * @brief Defines the API for the Simpson's rule integration module.
 */
typedef struct Simpson {
    /**
     * @brief Creates a Simpson's rule integration handle.
     * @param f The function to integrate.
     * @param a The lower limit of integration.
     * @param b The upper limit of integration.
     * @param max_iter The maximum number of iterations.
     * @param outSimpson A pointer to a pointer to the integration approximation handle to create.
     * @param name A name for the integration.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*simpson_create)(double (*f)(double x), double a, double b, size_t max_iter, Simpson *outSimpson, const char *name);

    /**
     * @brief Performs Simpson's rule integration.
     * @param inSimpson A pointer to the integration approximation handle.
     * @param outApproxIntegral A pointer to a variable to store the approximate integral.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*simpson_integration)(const Simpson *inSimpson, double *outApproxIntegral);

    /**
     * @brief Destroys a Simpson's rule integration handle.
     * @param inSimpson A pointer to the integration approximation handle to destroy.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*simpson_destroy)(Simpson *inSimpson);
}SimpsonAPI;

/**
 * @brief The global instance of the Simpson's rule integration API.
 */
extern const SimpsonAPI SI;


#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_SIMPSON_H
