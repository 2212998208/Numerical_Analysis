//
// Created by ninico on 2025/10/18.
//

#ifndef NUMERICAL_ANALYSIS_NLX_H
#define NUMERICAL_ANALYSIS_NLX_H
#ifndef NLS_H
#define NLS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/*===========================
 * 错误码（0 为成功）
 *===========================*/
typedef enum {
    NLS_SUCCESS        = 0,  /* 成功 */
    NLS_ERR_NULL_PTR   = 1,  /* 空指针 */
    NLS_ERR_ALLOC      = 2,  /* 内存分配失败 */
    NLS_ERR_DIM        = 3,  /* 维度非法 */
    NLS_ERR_BAD_STATE  = 4,  /* 状态错误（未设置函数等） */
    NLS_ERR_SINGULAR   = 5,  /* 线性系统奇异/病态 */
    NLS_ERR_MAXIT      = 6,  /* 迭代未收敛 */
    NLS_ERR_JACOBIAN   = 7,  /* 雅可比计算失败 */
    NLS_ERR_INTERNAL   = 8   /* 内部错误 */
} nls_err_t;

/* 文本化错误 */
const char* nls_strerror(int err);

/*===========================
 * 不透明句柄（隐式指针）
 *===========================*/
typedef struct nls_solver nls_solver_t;

/*===========================
 * 回调类型
 *===========================*/
/* F: ℝ^n → ℝ^n */
typedef int (*nls_fun_F)(const double* x, double* f, void* user);
/* J: 雅可比，行主序 n×n；可传 NULL 则用数值差分 */
typedef int (*nls_fun_J)(const double* x, double* J, void* user);

/*===========================
 * 选项
 *===========================*/
typedef struct {
    int    max_iter;           /* 默认 50 */
    double tol_res;            /* 残差阈值，默认 1e-10 (||F||_inf) */
    double tol_step;           /* 步长阈值，默认 1e-12 (||Δx||_inf) */
    double fd_eps;             /* 数值差分基准步长，默认 sqrt(DBL_EPSILON) */
    double ls_c;               /* 线搜索 Armijo 系数，默认 1e-4 */
    double ls_beta;            /* 回溯因子 ∈(0,1)，默认 0.5 */
    int    ls_max_backtrack;   /* 最大回溯次数，默认 20 */
} nls_options_t;

/*===========================
 * 面向对象风格 API
 *===========================*/
int nls_create(size_t n, nls_solver_t** out);
int nls_set_function(nls_solver_t* s, nls_fun_F F, void* user);
int nls_set_jacobian(nls_solver_t* s, nls_fun_J J);
int nls_set_options(nls_solver_t* s, const nls_options_t* opt);

/* 解方程 F(x)=0
 *  输入：x_inout 作为初值
 *  输出：x_inout 被原地更新为解
 */
int nls_solve(nls_solver_t* s, double* x_inout);

/* 统计量：迭代步数、最后一次 ||F||_inf 与 ||Δx||_inf */
int nls_get_stats(const nls_solver_t* s, int* iters, double* res_inf, double* step_inf);

int nls_destroy(nls_solver_t** ps);

/*===========================
 * 便捷一次性 API
 *===========================*/
int nls_solve_raw(size_t n,
                  nls_fun_F F, nls_fun_J J, void* user,
                  const nls_options_t* opt,
                  double* x_inout);

#ifdef __cplusplus
}
#endif
#endif /* NLS_H */

#endif //NUMERICAL_ANALYSIS_NLX_H