#ifndef NUMERICAL_ANALYSIS_SECANT_H
#define NUMERICAL_ANALYSIS_SECANT_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief An opaque handle to a non-linear scant.
 */
typedef struct NonLinearScant *NonLinearScant;


/**
 * @brief Represents the error codes for the secant module.
 */
typedef enum Secant_Err {
    SECANT_OK = 0, /**< No error. */
    SECANT_ERR_NOMEM = 1, /**< Memory allocation failed. */
    SECANT_ERR_INVALID = 2, /**< Invalid input. */
    SECANT_ERR_MAXITER = 3, /**< Maximum iterations reached. */
    SECANT_ERR_DIVIDE_BY_ZERO = 4 /**< Division by zero. */
} Secant_Err;

/**
 * @brief Defines the API for the secant module.
 */
typedef struct Secant {
    /**
     * @brief Solves for the root of a function using the secant method.
     * @param outScant A pointer to the non-linear scant.
     * @param outRoot A pointer to a variable to store the root.
     * @return A Secant_Err error code.
     */
    Secant_Err (*secant_solve)(const NonLinearScant *outScant, double *outRoot);

    /**
     * @brief Creates a non-linear scant for the secant solver.
     * @param f The function for which to find the root.
     * @param x0 The first initial guess.
     * @param x1 The second initial guess.
     * @param tol The tolerance.
     * @param max_iter The maximum number of iterations.
     * @param outScant A pointer to a pointer to the non-linear scant to create.
     * @param name A name for the non-linear scant.
     * @return A Secant_Err error code.
     */
    Secant_Err (*NonLinearScant_create)(double (*f)(double x), const double x0, const double x1,
        const double tol, const size_t max_iter,
        NonLinearScant *outScant,
        const char *name);

    /**
     * @brief Destroys a non-linear scant.
     * @param inScant A pointer to the non-linear scant to destroy.
     * @return A Secant_Err error code.
     */
    Secant_Err (*NonLinearScant_destroy)(const NonLinearScant *inScant);
}SecantAPI;

/**
 * @brief The global instance of the secant API.
 */
extern const SecantAPI Secant;


#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_SECANT_H
