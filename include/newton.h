#ifndef NUMERICAL_ANALYSIS_NEWTON_H
#define NUMERICAL_ANALYSIS_NEWTON_H
#include "lagrange.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents a dataset for Newton interpolation.
 */
typedef struct NewtonDataSet NewtonDataSet;

/**
 * @brief Represents the error codes for the Newton interpolation module.
 */
typedef enum Newton_Err {
    NEWTON_OK = 0, /**< No error. */
    NEWTON_ERR_NOMEM = 1, /**< Memory allocation failed. */
    NEWTON_ERR_INVALID = 2, /**< Invalid input. */
    NEWTON_ERR_DIVBYZERO = 3 /**< Division by zero. */
} Newton_Err;

/**
 * @brief Defines the API for the Newton interpolation module.
 */
typedef struct {
    /**
     * @brief Performs Newton interpolation.
     * @param inNewtonDataSet A pointer to a pointer to the Newton dataset.
     * @param intData A pointer to a pointer to the Lagrange dataset.
     * @param x The x-coordinate at which to interpolate.
     * @param outY A pointer to a variable to store the interpolated y-coordinate.
     * @return A Newton_Err error code.
     */
    Newton_Err (*newton_interpolate)(NewtonDataSet **inNewtonDataSet, DataSet **intData, double x, double *outY);

    /**
     * @brief Creates a Newton dataset.
     * @param outDataset A pointer to a pointer to the dataset to create.
     * @param size The number of points in the dataset.
     * @return A Newton_Err error code.
     */
    Newton_Err (*create_newton_dataset)(NewtonDataSet **outDataset, size_t size);

    /**
     * @brief Destroys a Newton dataset.
     * @param outDataset A pointer to a pointer to the dataset to destroy.
     * @return A Newton_Err error code.
     */
    Newton_Err (*destroy_dataset)(NewtonDataSet **outDataset);

    /**
     * @brief Prints a Newton dataset.
     * @param outDataset A pointer to a pointer to the dataset to print.
     * @return A Newton_Err error code.
     */
    Newton_Err (*print_newton_dataset)(const NewtonDataSet **outDataset);
}NewtonAPI;


/**
 * @brief The global instance of the Newton interpolation API.
 */
extern const NewtonAPI Newton;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_NEWTON_H
