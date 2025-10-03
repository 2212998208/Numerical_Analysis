#ifndef NUMERICAL_ANALYSIS_TRI_H
#define NUMERICAL_ANALYSIS_TRI_H
#ifndef TRI_H
#define TRI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h> /* size_t */

/*===========================
 * 错误码（0 为成功）
 *===========================*/
typedef enum {
    TRI_SUCCESS          = 0,  /* 成功 */
    TRI_ERR_NULL_PTR     = 1,  /* 空指针 */
    TRI_ERR_ALLOC        = 2,  /* 内存分配失败 */
    TRI_ERR_DIM          = 3,  /* 维度非法（n==0等） */
    TRI_ERR_BAD_STATE    = 4,  /* 状态错误（未设置系数/右端等） */
    TRI_ERR_ZERO_PIVOT   = 5,  /* 消元遇到零主元/近零主元 */
    TRI_ERR_SINGULAR     = 6,  /* 矩阵奇异或严重病态 */
    TRI_ERR_INTERNAL     = 7   /* 其他内部错误 */
} tri_err_t;

/* 将错误码转为可读字符串，便于日志或测试输出 */
const char* tri_strerror(int err);

/*===========================
 * 不透明句柄（隐式指针封装）
 *===========================*/
typedef struct tri_solver tri_solver_t;

/*===========================
 * 面向对象风格 API
 *===========================*/

/* 创建求解器：仅记录规模 n，不分配系数数组（待 set 接口拷贝） */
int tri_create(size_t n, tri_solver_t** out);

/* 设置三对角系数（会**拷贝**到内部缓冲区）
   a: 下对角 a[0..n-2]
   b: 主对角 b[0..n-1]
   c: 上对角 c[0..n-2]
*/
int tri_set_coeffs(tri_solver_t* solver,
                   const double* a,
                   const double* b,
                   const double* c);

/* 设置右端项 d（会**拷贝**到内部缓冲区）d[0..n-1] */
int tri_set_rhs(tri_solver_t* solver, const double* d);

/* 允许设置（或查询）数值容差（用于“近零主元”判断），默认 1e-12 */
int tri_set_tolerance(tri_solver_t* solver, double tol);
int tri_get_tolerance(const tri_solver_t* solver, double* tol_out);

/* 求解：输出 x[0..n-1]，调用者提供外部缓冲区 */
int tri_solve(const tri_solver_t* solver, double* x_out);

/* 获取规模 n */
int tri_get_n(const tri_solver_t* solver, size_t* n_out);

/* 释放资源（置 *solver=NULL） */
int tri_destroy(tri_solver_t** solver);

/*===========================
 * 直接数组风格的便捷 API
 *===========================*/
/* 基于原始数组的一次性求解（不会保存状态）
   a[0..n-2], b[0..n-1], c[0..n-2], d[0..n-1] -> x[0..n-1]
   tol 为主元容差（传 <=0 则使用默认 1e-12）
*/
int tri_solve_raw(size_t n,
                  const double* a,
                  const double* b,
                  const double* c,
                  const double* d,
                  double* x_out,
                  double tol);

#ifdef __cplusplus
}
#endif

#endif /* TRI_H */

#endif //NUMERICAL_ANALYSIS_TRI_H