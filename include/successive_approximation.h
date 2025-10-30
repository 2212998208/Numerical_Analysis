#ifndef NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H
#define NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief An opaque handle to a non-linear successive approximation solver.
 */
typedef struct NonlinearSA *NonlinearSA;

/**
 * @brief Represents the error codes for the successive approximation module.
 */
typedef enum Succesive_Err {
    SA_OK = 0, /**< No error. */
    SA_ERR_NOMEM = 1, /**< Memory allocation failed. */
    SA_ERR_INVAL = 2, /**< Invalid input. */
    SA_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    SA_ERR_DIVIDE_BY_ZERO = 4, /**< Division by zero. */
    SA_ERR_NOAPPROXIMATION = 5, /**< No approximation could be found. */
}Succesive_Err;


/**
 * @brief Defines the API for the successive approximation module.
 */
typedef struct Succesive_Approximation {
    /**
     * @brief Solves for the root of a function using the successive approximation method.
     * @param outSA A pointer to the non-linear successive approximation solver handle.
     * @param outRoot A pointer to a variable to store the root.
     * @return A Succesive_Err error code.
     */
    Succesive_Err (*nonlinear_sa_solve)(const NonlinearSA *outSA, double *outRoot);

    /**
     * @brief Creates a non-linear successive approximation solver handle.
     * @param g The function g(x) for the fixed-point iteration x = g(x).
     * @param x0 The initial guess.
     * @param tol The tolerance.
     * @param max_iter The maximum number of iterations.
     * @param outSA A pointer to a pointer to the solver handle to create.
     * @param name A name for the solver.
     * @return A Succesive_Err error code.
     */
    Succesive_Err (*NonlinearSA_create)(double (*g)(double x), const double x0,
        const double tol, const size_t max_iter,
        NonlinearSA *outSA,
        const char *name);

    /**
     * @brief Destroys a non-linear successive approximation solver handle.
     * @param inSA A pointer to the solver handle to destroy.
     * @return A Succesive_Err error code.
     */
    Succesive_Err (*NonlinearSA_destroy)(const NonlinearSA *inSA);
}SAAPI;

/**
 * @brief The global instance of the successive approximation API.
 */
extern const SAAPI SA;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_SUCCESSIVE_APPROXIMATION_H
