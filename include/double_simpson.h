#ifndef NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H
#define NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents the error codes for the double Simpson's rule integration module.
 * @note This enum is identical to the one in simpson.h.
 */
typedef enum Simpson_Err {
    DOUBLE_SIMPSON_OK = 0, /**< No error. */
    DOUBLE_SIMPSON_ERR_NOMEM = 1, /**< Memory allocation failed. */
    DOUBLE_SIMPSON_ERR_INVAL = 2, /**< Invalid input. */
    DOUBLE_SIMPSON_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    DOUBLE_SIMPSON_ERR_DIVIDE_BY_ZERO = 4, /**< Division by zero. */
}Double_Simpson_Err;

/**
 * @brief An opaque handle to a double integration approximation.
 */
typedef struct DoubleIntegrationApproximation *DoubleSimpson;


/**
 * @brief Defines the API for the double Simpson's rule integration module.
 */
typedef struct DoubleSimpson {
    /**
     * @brief Creates a double Simpson's rule integration handle.
     * @param f The function to integrate.
     * @param x_a The lower limit of integration for x.
     * @param x_b The upper limit of integration for x.
     * @param y_c The lower limit of integration for y.
     * @param y_d The upper limit of integration for y.
     * @param n The number of intervals for x.
     * @param m The number of intervals for y.
     * @param outDoubleSimpson A pointer to a pointer to the double integration approximation handle to create.
     * @param name A name for the integration.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*double_simpson_create)(double (*f)(double x, double y),
                                  const double x_a, const double x_b,
                                  const double y_c, const double y_d,
                                  const size_t n, const size_t m,
                                  DoubleSimpson *outDoubleSimpson,
                                  const char *name);
    /**
     * @brief Destroys a double Simpson's rule integration handle.
     * @param inDoubleSimpson A pointer to the double integration approximation handle to destroy.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*double_simpson_destroy)(DoubleSimpson *inDoubleSimpson);

    /**
     * @brief Performs double Simpson's rule integration.
     * @param inDoubleSimpson A pointer to the double integration approximation handle.
     * @param outApproxIntegral A pointer to a variable to store the approximate integral.
     * @return A Simpson_Err error code.
     */
    Simpson_Err (*double_simpson_integrate)(const DoubleSimpson *inDoubleSimpson, double *outApproxIntegral);
} DoubleSimpsonAPI;


/**
 * @brief The global instance of the double Simpson's rule integration API.
 */
extern const DoubleSimpsonAPI DS;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_DOUBLE_SIMPSON_H
