//
// Created by ninico on 2025/10/18.
//

#ifndef NUMERICAL_ANALYSIS_BVP_H
#define NUMERICAL_ANALYSIS_BVP_H
#ifndef BVP_H
#define BVP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h> /* size_t */

/*===========================
 * 错误码（0 为成功）
 *===========================*/
typedef enum {
    BVP_SUCCESS        = 0,  /* 成功 */
    BVP_ERR_NULL_PTR   = 1,  /* 空指针 */
    BVP_ERR_ALLOC      = 2,  /* 内存分配失败 */
    BVP_ERR_DIM        = 3,  /* 维度非法（m,N 等） */
    BVP_ERR_BAD_STATE  = 4,  /* 状态错误（未设置函数/边界等） */
    BVP_ERR_SINGULAR   = 5,  /* 线性系统奇异/严重病态 */
    BVP_ERR_MAXIT      = 6,  /* 迭代未收敛（达到最大步数） */
    BVP_ERR_JACOBIAN   = 7,  /* 雅可比计算失败 */
    BVP_ERR_INTERNAL   = 8   /* 内部错误 */
} bvp_err_t;

const char* bvp_strerror(int err);

/*===========================
 * 不透明句柄（隐式指针）
 *===========================*/
typedef struct bvp_solver bvp_solver_t;

/*===========================
 * 回调类型
 * G: y'' = G(x, y, y')
 * J_y / J_yp：可选解析雅可比（若未提供，将用数值差分近似）
 *===========================*/
typedef int (*bvp_fun_G)(
    double x, const double* y, const double* yp,
    double* out,                /* out: 长度 m */
    void* user
);

typedef int (*bvp_fun_J)(
    double x, const double* y, const double* yp,
    double* J /* out: 行主序 m×m */,
    void* user
);

/*===========================
 * 创建 / 销毁 与设置
 *===========================*/
int bvp_create(size_t m, size_t N, double a, double b, bvp_solver_t** out);
/* 设置强制项函数与可选雅可比 */
int bvp_set_function(bvp_solver_t* s, bvp_fun_G G, void* user);
int bvp_set_jacobians(bvp_solver_t* s, bvp_fun_J J_y, bvp_fun_J J_yp);

/* Dirichlet 边界：y(a)=ya, y(b)=yb （向量，长度 m） */
int bvp_set_boundary(bvp_solver_t* s, const double* ya, const double* yb);

/* 初值：内点共 N 个节点（不含端点），长度 N*m；若不设置，采用端点线性插值作为初值 */
int bvp_set_initial_guess(bvp_solver_t* s, const double* y_interior);

/* 选项：最大迭代、残差阈值、步长阈值（<=0 则使用默认） */
int bvp_set_options(bvp_solver_t* s, int max_iter, double tol_res, double tol_step);

/* 求解：输出所有 N+2 个节点（含端点）上的解，长度 (N+2)*m，按节点快排布：i=0..N+1，每节点 m 个分量 */
int bvp_solve(bvp_solver_t* s, double* y_all_out);

/* 获取网格坐标（长度 N+2），等距 [a,b]；调用方提供缓冲区 */
int bvp_get_grid(const bvp_solver_t* s, double* x_all_out);

/* 获取 m、N、区间 [a,b] */
int bvp_get_meta(const bvp_solver_t* s, size_t* m, size_t* N, double* a, double* b);

/* 销毁 */
int bvp_destroy(bvp_solver_t** ps);

/*===========================
 * 便捷一次性 API
 *===========================*/
typedef struct {
    int max_iter;     /* 默认 50 */
    double tol_res;   /* 默认 1e-8 */
    double tol_step;  /* 默认 1e-10 */
} bvp_options_t;

/* 若 J_y/J_yp 传 NULL，将自动用数值差分近似 */
int bvp_solve_raw(size_t m, size_t N, double a, double b,
                  const double* ya, const double* yb,
                  bvp_fun_G G, bvp_fun_J J_y, bvp_fun_J J_yp, void* user,
                  const bvp_options_t* opt,
                  double* y_all_out);

#ifdef __cplusplus
}
#endif
#endif /* BVP_H */

#endif //NUMERICAL_ANALYSIS_BVP_H