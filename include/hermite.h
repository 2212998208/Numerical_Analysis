
#ifndef NUMERICAL_ANALYSIS_HERMITE_H
#define NUMERICAL_ANALYSIS_HERMITE_H
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef struct HermitePoint* HermitePoint;
typedef struct HermiteDataset* HermiteDataset;
typedef struct HermiteInterpolator* HermiteInterpolator;


typedef enum Hermite_Err {
    HERMITE_OK = 0,
    HERMITE_ERR_NOMEM = 1,
    HERMITE_ERR_INVALID = 2,
    HERMITE_ERR_DIVBYZERO = 3
} Hermite_Err;


typedef struct {
    Hermite_Err (*hermite_create_interpolator)(HermiteDataset *inDataset, HermiteInterpolator *outInterpolator);
    Hermite_Err (*hermite_destroy_interpolator)(HermiteInterpolator *inInterpolator);
    HermiteDataset (*create_hermite_dataset)(size_t size, const double *x, const double *y, const double *dy);
    Hermite_Err (*destroy_hermite_dataset)(HermiteDataset *inDataset);
    Hermite_Err (*hermite_evaluate)(const HermiteInterpolator *inInterpolator, double x, double *outY, double *outDy);
}HermiteAPI;







extern const HermiteAPI Hermite;
#ifdef __cplusplus
}
#endif

#endif //NUMERICAL_ANALYSIS_HERMITE_H