#ifndef NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H
#define NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents the error codes for the Newton-Raphson module.
 */
typedef enum NewtonRaphson_Err {
    NEWTON_RAPHSON_OK = 0, /**< No error. */
    NEWTON_RAPHSON_ERR_NOMEM = 1, /**< Memory allocation failed. */
    NEWTON_RAPHSON_ERR_INVALID = 2, /**< Invalid input. */
    NEWTON_RAPHSON_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    NEWTON_RAPHSON_ERR_DERIVATIVE_ZERO = 4, /**< The derivative is zero. */
    NEWTON_RAPHSON_ERR_DERIVATIVE_UNSTABLE = 5 /**< The derivative is unstable. */
} NewtonRaphson_Err;

/**
 * @brief An opaque handle to a non-linear range.
 */
typedef struct NonLinearRange *NonLinearRange;

/**
 * @brief Defines the API for the Newton-Raphson module.
 */
typedef struct NewtonRaphson {
    /**
     * @brief Creates a non-linear range for the Newton-Raphson solver.
     * @param f The function for which to find the root.
     * @param x0 The initial guess.
     * @param tol The tolerance.
     * @param max_iter The maximum number of iterations.
     * @param outRange A pointer to a pointer to the non-linear range to create.
     * @param name A name for the non-linear range.
     * @return A NewtonRaphson_Err error code.
     */
    NewtonRaphson_Err (*NonLinearRange_create)(double (*f)(double x), const double x0,
        const double tol, const size_t max_iter,
        NonLinearRange *outRange,
        const char *name);

    /**
     * @brief Destroys a non-linear range.
     * @param inRange A pointer to the non-linear range to destroy.
     * @return A NewtonRaphson_Err error code.
     */
    NewtonRaphson_Err (*NonLinearRange_destroy)(const NonLinearRange *inRange);

    /**
     * @brief Solves for the root of a function using the Newton-Raphson method.
     * @param outRange A pointer to the non-linear range.
     * @param outRoot A pointer to a variable to store the root.
     * @return A NewtonRaphson_Err error code.
     */
    NewtonRaphson_Err (*newton_raphson_solve)(const NonLinearRange *outRange, double *outRoot);
}NewtonRaphsonAPI;


/**
 * @brief The global instance of the Newton-Raphson API.
 */
extern const NewtonRaphsonAPI NewtonRaphson;
#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_NEWTON_RAPHSON_H
