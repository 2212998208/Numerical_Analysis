#ifndef NUMERICAL_ANALYSIS_LAGRANGE_H
#define NUMERICAL_ANALYSIS_LAGRANGE_H
#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Represents a dataset for Lagrange interpolation.
 */
typedef struct DataSet DataSet;

/**
 * @brief Represents a point with x and y coordinates.
 */
typedef struct Point {
        double x; /**< The x-coordinate of the point. */
        double y; /**< The y-coordinate of the point. */
} Point;

/**
 * @brief Represents the error codes for the Lagrange interpolation module.
 */
typedef enum Lagrange_Err {
    LAGRANGE_OK = 0, /**< No error. */
    LAGRANGE_ERR_NOMEM = 1, /**< Memory allocation failed. */
    LAGRANGE_ERR_INVALID = 2, /**< Invalid input. */
    LAGRANGE_ERR_DIVBYZERO = 3 /**< Division by zero. */
} Lagrange_Err;


/**
 * @brief Creates a Point object.
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @return The created Point object.
 */
static inline Point point_make(const double x, const double y) {
    const Point p = {x, y};
    return p;
}


/**
 * @brief Defines the API for the Lagrange interpolation module.
 */
typedef struct {
    /**
     * @brief Performs Lagrange interpolation.
     * @param dataset The dataset to use for interpolation.
     * @param x The x-coordinate at which to interpolate.
     * @return The interpolated y-coordinate.
     */
    double (*lagrange_interpolate)(DataSet *dataset, double x);

    /**
     * @brief Creates an empty dataset.
     * @return A pointer to the created dataset.
     */
    DataSet *(*empty_dataset)(void);

    /**
     * @brief Creates a dataset from an array of points.
     * @param dataset A pointer to a pointer to the dataset to create.
     * @param points An array of points.
     * @param size The number of points in the array.
     * @return A Lagrange_Err error code.
     */
    Lagrange_Err (*create_dataset)(DataSet **dataset, Point *points, size_t size);

    /**
     * @brief Creates a dataset from a given size.
     * @param dataset A pointer to a pointer to the dataset to create.
     * @param size The number of points in the dataset.
     * @return A Lagrange_Err error code.
     */
    Lagrange_Err (*create_dataset_from_points)(DataSet **dataset, size_t size);

    /**
     * @brief Destroys a dataset.
     * @param dataset A pointer to the dataset to destroy.
     */
    void (*destroy_dataset)(DataSet *dataset);

    /**
     * @brief Gets the points from a dataset.
     * @param outDataset A pointer to a pointer to the dataset.
     * @return A pointer to the array of points.
     */
    Point *(*get_points)(DataSet **outDataset);
} LagrangeAPI;


/**
 * @brief The global instance of the Lagrange interpolation API.
 */
extern const LagrangeAPI Lagrange;

#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_LAGRANGE_H
