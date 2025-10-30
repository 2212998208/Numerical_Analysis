#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#pragma once
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Defines a function pointer for the integrand function.
 * @param x The independent variable.
 * @param user_data A pointer to user-defined data that can be passed to the function.
 * @return The value of the function at x.
 */
typedef double (*IntegrandFn)(double x, void *user_data);

/**
 * @brief Represents the status of the adaptive integration process.
 */
typedef enum {
    INTEGRATOR_OK = 0, /**< Integration completed successfully. */
    INTEGRATOR_MAX_STEPS_REACHED = 1 /**< The maximum number of iterations was reached. */
} IntegratorStatus;

/**
 * @brief Configuration for the adaptive integrator.
 */
typedef struct {
    double abs_tol;        /**< The absolute error tolerance. */
    double rel_tol;        /**< The relative error tolerance. */
    int    max_iterations; /**< The maximum number of iterations. */
} AdaptiveConfig;

/**
 * @brief Defines the API for the integrator module.
 */
typedef struct {
    /**
     * @brief Performs fixed-step RK4 integration.
     * @param f The integrand function.
     * @param user The user data to pass to the integrand function.
     * @param a The lower limit of integration.
     * @param b The upper limit of integration.
     * @param steps The number of steps to use.
     * @return The value of the integral.
     */
    double (*rk4_fixed)(IntegrandFn f, void *user, double a, double b, int steps);

    /**
     * @brief Performs adaptive RK4 integration.
     * @param f The integrand function.
     * @param user The user data to pass to the integrand function.
     * @param a The lower limit of integration.
     * @param b The upper limit of integration.
     * @param cfg The adaptive configuration.
     * @param status A pointer to a variable to store the status of the integration.
     * @return The value of the integral.
     */
    double (*rk4_adaptive)(IntegrandFn f, void *user, double a, double b,
                           AdaptiveConfig cfg, IntegratorStatus *status);
} IntegratorAPI;

/**
 * @brief The global instance of the integrator API.
 */
extern const IntegratorAPI Integrator;

#ifdef __cplusplus
}
#endif

#endif
