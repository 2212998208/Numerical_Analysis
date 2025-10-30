
#ifndef NUMERICAL_ANALYSIS_HERMITE_H
#define NUMERICAL_ANALYSIS_HERMITE_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief An opaque handle to a Hermite point.
 */
typedef struct HermitePoint* HermitePoint;

/**
 * @brief An opaque handle to a Hermite dataset.
 */
typedef struct HermiteDataset* HermiteDataset;

/**
 * @brief An opaque handle to a Hermite interpolator.
 */
typedef struct HermiteInterpolator* HermiteInterpolator;


/**
 * @brief Represents the error codes for the Hermite interpolation module.
 */
typedef enum Hermite_Err {
    HERMITE_OK = 0, /**< No error. */
    HERMITE_ERR_NOMEM = 1, /**< Memory allocation failed. */
    HERMITE_ERR_INVALID = 2, /**< Invalid input. */
    HERMITE_ERR_DIVBYZERO = 3 /**< Division by zero. */
} Hermite_Err;


/**
 * @brief Defines the API for the Hermite interpolation module.
 */
typedef struct {
    /**
     * @brief Creates a Hermite interpolator from a dataset.
     * @param inDataset A pointer to the Hermite dataset.
     * @param outInterpolator A pointer to a pointer to the interpolator to create.
     * @return A Hermite_Err error code.
     */
    Hermite_Err (*hermite_create_interpolator)(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator);

    /**
     * @brief Destroys a Hermite interpolator.
     * @param inInterpolator A pointer to the interpolator to destroy.
     * @return A Hermite_Err error code.
     */
    Hermite_Err (*hermite_destroy_interpolator)(HermiteInterpolator *inInterpolator);

    /**
     * @brief Creates a Hermite dataset.
     * @param size The number of points in the dataset.
     * @param x An array of x-coordinates.
     * @param y An array of y-coordinates.
     * @param dy An array of derivatives.
     * @return A pointer to the created dataset.
     */
    HermiteDataset (*create_hermite_dataset)(size_t size, const double *x, const double *y, const double *dy);

    /**
     * @brief Destroys a Hermite dataset.
     * @param inDataset A pointer to the dataset to destroy.
     * @return A Hermite_Err error code.
     */
    Hermite_Err (*destroy_hermite_dataset)(HermiteDataset *inDataset);

    /**
     * @brief Evaluates the Hermite interpolator at a given point.
     * @param inInterpolator A pointer to the interpolator.
     * @param x The x-coordinate at which to evaluate.
     * @param outY A pointer to a variable to store the interpolated y-coordinate.
     * @param outDy A pointer to a variable to store the interpolated derivative.
     * @return A Hermite_Err error code.
     */
    Hermite_Err (*hermite_evaluate)(const HermiteInterpolator *inInterpolator, double x, double *outY, double *outDy);
}HermiteAPI;







/**
 * @brief The global instance of the Hermite interpolation API.
 */
extern const HermiteAPI Hermite;
#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_HERMITE_H
