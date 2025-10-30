#ifndef NUMERICAL_ANALYSIS_BISECTION_H
#define NUMERICAL_ANALYSIS_BISECTION_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief An opaque handle to a bisection range.
 */
typedef struct BisectionRange *BisectionRange;



/**
 * @brief Represents the error codes for the bisection module.
 */
typedef enum Bisection_Err {
    BISECTION_OK = 0, /**< No error. */
    BISECTION_ERR_NOMEM = 1, /**< Memory allocation failed. */
    BISECTION_ERR_INVALID = 2, /**< Invalid input. */
    BISECTION_ERR_MAXITER = 3 /**< Maximum iterations reached. */
} Bisection_Err;

/**
 * @brief Defines the API for the bisection module.
 */
typedef struct Bisection {
    /**
     * @brief Solves for the root of a function within a given range using the bisection method.
     * @param outRange A pointer to the bisection range handle.
     * @return A Bisection_Err error code.
     */
    Bisection_Err (*bisection_solve)(BisectionRange *outRange);

    /**
     * @brief Creates a bisection range.
     * @param f The function for which to find the root.
     * @param a The lower bound of the range.
     * @param b The upper bound of the range.
     * @param tol The tolerance for the root.
     * @param outRange A pointer to a pointer to the bisection range to create.
     * @param name A name for the bisection range.
     * @return A Bisection_Err error code.
     */
    Bisection_Err (*bisection_create)(double (*f)(double x), double a,
        double b, double tol,
        BisectionRange *outRange,
        const char *name);

    /**
     * @brief Destroys a bisection range.
     * @param outRange A pointer to the bisection range to destroy.
     * @return A Bisection_Err error code.
     */
    Bisection_Err (*bisection_destroy)(BisectionRange *outRange);

    /**
     * @brief Gets the midpoint of the current bisection range.
     * @param outRange A pointer to the bisection range.
     * @return The midpoint of the range.
     */
    double (*bisection_get_midpoint)(const BisectionRange *outRange);
}BisectionAPI;


/**
 * @brief The global instance of the bisection API.
 */
extern const BisectionAPI Bisection;

#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_BISECTION_H
