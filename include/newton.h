#ifndef NUMERICAL_ANALYSIS_NEWTON_H
#define NUMERICAL_ANALYSIS_NEWTON_H
#include "lagrange.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct NewtonDataSet NewtonDataSet;

typedef enum Newton_Err {
    NEWTON_OK = 0,
    NEWTON_ERR_NOMEM = 1,
    NEWTON_ERR_INVALID = 2,
    NEWTON_ERR_DIVBYZERO = 3
} Newton_Err;

typedef struct {
    Newton_Err (*newton_interpolate)(NewtonDataSet **inNewtonDataSet, DataSet **intData, double x, double *outY);
    Newton_Err (*create_newton_dataset)(NewtonDataSet **outDataset, size_t size);
    Newton_Err (*destroy_dataset)(NewtonDataSet **outDataset);
    Newton_Err (*print_newton_dataset)(const NewtonDataSet **outDataset);
}NewtonAPI;


extern const NewtonAPI Newton;
#ifdef __cplusplus
}
#endif
#endif //NUMERICAL_ANALYSIS_NEWTON_H