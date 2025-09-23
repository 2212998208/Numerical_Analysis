# Numerical_Analysis 项目说明（完整版）

本仓库实现常见数值分析算法：插值（Lagrange、Newton、Hermite）、方程求根（Bisection、Newton–Raphson、Secant、逐次逼近）、数值积分（梯形、Simpson、双重Simpson、基于 RK4 的积分工具）。项目采用 C 语言 + CMake 组织，接口统一以“API 结构体 + 函数指针”的方式暴露。

---
## 目录
- 功能概览
- 代码组织结构
- CMake 构建与目标组织
- API 设计风格（统一说明）
- 各模块 API 与数据类型
  - Integrator（RK4 积分工具）
  - Lagrange 插值
  - Newton 插值
  - Hermite 插值
  - Bisection 二分法
  - Newton–Raphson 牛顿法
  - Secant 割线法
  - Trapezoidal 梯形求积
  - Simpson 单变量 Simpson
  - Double Simpson 二重积分 Simpson
  - Successive Approximation 逐次逼近
- 使用示例（选摘）
- 测试说明
- 开发建议与代码注释规范

---
## 功能概览
- 插值算法：
  - Lagrange 多项式插值
  - Newton 多项式插值
  - Hermite 插值（带导数信息）
- 非线性方程求根：
  - Bisection（区间二分）
  - Newton–Raphson（切线法）
  - Secant（割线法）
  - Successive Approximation（逐次逼近/不动点迭代）
- 数值积分：
  - Trapezoidal（梯形法）
  - Simpson（单变量）
  - Double Simpson（二重积分）
  - Integrator（独立的 RK4 工具，含固定步长与自适应配置）

---
## 代码组织结构
```
Numerical_Analysis/
├─ CMakeLists.txt              # CMake 构建配置
├─ main.c                      # 示例/入口（与各 src 组合构成主程序）
├─ include/                    # 头文件（公共 API）
│  ├─ integrator.h             # RK4 积分工具 API
│  ├─ lagrange.h               # Lagrange 插值 API
│  ├─ newton.h                 # Newton 插值 API
│  ├─ hermite.h                # Hermite 插值 API
│  ├─ bisection.h              # 二分法 API
│  ├─ newton_raphson.h         # 牛顿法 API
│  ├─ secant.h                 # 割线法 API
│  ├─ Trapezoidal.h            # 梯形积分 API
│  ├─ simpson.h                # Simpson 积分 API
│  ├─ double_simpson.h         # 双重 Simpson 积分 API
│  ├─ successive_approximation.h # 逐次逼近 API
│  ├─ test_integrator.h        # 积分测试声明（run_all_tests）
│  ├─ test_lagrange.h          # Lagrange 测试声明（run_all_tests）
│  └─ test_newton.h            # Newton 测试声明（run_all_tests）
├─ src/                        # 源码实现
│  ├─ integrator.c, lagrange.c, newton.c, hermite.c
│  ├─ bisection.c, newton_raphson.c, secant.c
│  ├─ Trapezoidal.c, simpson.c, double_simpson.c, successive_approximation.c
├─ tests/                      # 各模块测试
│  ├─ test_integrator.c, test_lagrange.c, test_newton.c, test_hermite.c, ...
├─ docs/
│  ├─ Integrator_MindMap.svg
│  └─ lagrange.md              # Lagrange 模块改进建议（详见文档）
└─ cmake-build-*/              # 构建产物（由 IDE/CMake 生成）
```

---
## CMake 构建与目标组织
CMake 以“模块源 + 测试源”组合生成多个可执行文件：
- 主程序：`Numerical_Analysis`，由 `main.c` 与各模块 `src/*.c` 组合。
- 测试程序：为每个模块生成独立的测试可执行文件，例如：
  - `Numerical_Analysis_tests_lagrange` 由 `src/lagrange.c + tests/test_lagrange.c` 组成
  - `Numerical_Analysis_tests_integrator` 由 `src/integrator.c + tests/test_integrator.c` 组成
  - 其余模块同理（bisection/newton/…），并通过 `enable_testing()` + `add_test()` 注册到 CTest

构建与运行（Windows，命令提示符）：
```
mkdir build
cd build
cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release

# 运行所有测试（可选）
ctest -C Release -V

# 或直接运行某个测试可执行文件
./Numerical_Analysis_tests_lagrange.exe
```
提示：若使用 MSVC/CLion，IDE 将自动生成并管理上述目标与测试。

---
## API 设计风格（统一说明）
- 每个功能模块以一个“API 结构体”的只读全局实例暴露，例如：
  - `extern const LagrangeAPI Lagrange;`
  - `extern const IntegratorAPI Integrator;`
- 结构体字段为函数指针，指向对应实现；源文件中通过静态内部函数完成初始化：
  - 这样可隐藏实现细节（内部函数 `static` 仅在本翻译单元可见），对外仅暴露单一入口对象。
  - 初始化常见写法：
    - 位置初始化：`const HermiteAPI Hermite = { hermite_create_interpolator, ... };`
    - 指定字段初始化（推荐可读）：`const HermiteAPI Hermite = { .hermite_evaluate = hermite_evaluate_impl, ... };`
- Opaque Handle（不透明句柄）
  - 多模块使用 `typedef struct Foo* Foo;` 作为句柄类型，配套 `create/destroy` 负责资源生命周期。
- 错误码
  - 各模块定义独立枚举，例如 `Lagrange_Err/Simpson_Err/...`，统一以 0 为 OK，其余表示异常。

---
## 各模块 API 与数据类型
以下内容摘自头文件，便于总览（详细参数含义请参考具体头文件注释与测试用例）。

### Integrator（include/integrator.h）
- 数据类型
  - `typedef double (*IntegrandFn)(double x, void *user_data);`
  - `typedef enum IntegratorStatus { INTEGRATOR_OK, INTEGRATOR_MAX_STEPS_REACHED }`;
  - `typedef struct AdaptiveConfig { double abs_tol; double rel_tol; int max_iterations; }`;
- API
  - `extern const IntegratorAPI Integrator;`
  - 成员：
    - `double (*rk4_fixed)(IntegrandFn f, void *user, double a, double b, int steps);`
    - `double (*rk4_adaptive)(IntegrandFn f, void *user, double a, double b, AdaptiveConfig cfg, IntegratorStatus *status);`

### Lagrange（include/lagrange.h）
- 数据类型
  - `typedef struct DataSet DataSet;`（不透明数据集）
  - `typedef struct Point { double x; double y; } Point;`
  - `typedef enum Lagrange_Err { LAGRANGE_OK, LAGRANGE_ERR_NOMEM, LAGRANGE_ERR_INVALID, LAGRANGE_ERR_DIVBYZERO }`;
  - 内联辅助：`static inline Point point_make(double x, double y);`
- API
  - `extern const LagrangeAPI Lagrange;`
  - 成员：
    - `double (*lagrange_interpolate)(DataSet *dataset, double x);`
    - `DataSet *(*empty_dataset)(void);`
    - `Lagrange_Err (*create_dataset)(DataSet **dataset, Point *points, size_t size);`
    - `Lagrange_Err (*create_dataset_from_points)(DataSet **dataset, size_t size);`
    - `void (*destroy_dataset)(DataSet *dataset);`
    - `Point *(*get_points)(DataSet **outDataset);`

### Newton（include/newton.h）
- 数据类型
  - `typedef struct NewtonDataSet NewtonDataSet;`
  - `typedef enum Newton_Err { NEWTON_OK, NEWTON_ERR_NOMEM, NEWTON_ERR_INVALID, NEWTON_ERR_DIVBYZERO }`;
- API
  - `extern const NewtonAPI Newton;`
  - 成员：
    - `Newton_Err (*newton_interpolate)(NewtonDataSet **inNewtonDataSet, DataSet **intData, double x, double *outY);`
    - `Newton_Err (*create_newton_dataset)(NewtonDataSet **outDataset, size_t size);`
    - `Newton_Err (*destroy_dataset)(NewtonDataSet **outDataset);`
    - `Newton_Err (*print_newton_dataset)(const NewtonDataSet **outDataset);`

### Hermite（include/hermite.h）
- 数据类型（不透明句柄）
  - `typedef struct HermitePoint* HermitePoint;`
  - `typedef struct HermiteDataset* HermiteDataset;`
  - `typedef struct HermiteInterpolator* HermiteInterpolator;`
  - `typedef enum Hermite_Err { HERMITE_OK, HERMITE_ERR_NOMEM, HERMITE_ERR_INVALID, HERMITE_ERR_DIVBYZERO }`;
- API
  - `extern const HermiteAPI Hermite;`
  - 成员：
    - `Hermite_Err (*hermite_create_interpolator)(HermiteDataset*, HermiteInterpolator*);`
    - `Hermite_Err (*hermite_destroy_interpolator)(HermiteInterpolator*);`
    - `HermiteDataset (*create_hermite_dataset)(size_t size, const double *x, const double *y, const double *dy);`
    - `Hermite_Err (*destroy_hermite_dataset)(HermiteDataset*);`
    - `Hermite_Err (*hermite_evaluate)(const HermiteInterpolator*, double x, double *outY, double *outDy);`

### Bisection（include/bisection.h）
- 句柄与错误码
  - `typedef struct BisectionRange *BisectionRange;`
  - `typedef enum Bisection_Err { BISECTION_OK, BISECTION_ERR_NOMEM, BISECTION_ERR_INVALID, BISECTION_ERR_MAXITER }`;
- API：`extern const BisectionAPI Bisection;`
  - `Bisection_Err (*bisection_create)(double (*f)(double), double a, double b, double tol, BisectionRange *outRange, const char *name);`
  - `Bisection_Err (*bisection_solve)(BisectionRange *outRange);`
  - `Bisection_Err (*bisection_destroy)(BisectionRange *outRange);`
  - `double (*bisection_get_midpoint)(const BisectionRange *outRange);`

### Newton–Raphson（include/newton_raphson.h）
- 句柄与错误码
  - `typedef struct NonLinearRange *NonLinearRange;`
  - `typedef enum NewtonRaphson_Err { NEWTON_RAPHSON_OK, ... }`;
- API：`extern const NewtonRaphsonAPI NewtonRaphson;`
  - `NewtonRaphson_Err (*NonLinearRange_create)(double (*f)(double), double x0, double tol, size_t max_iter, NonLinearRange *outRange, const char *name);`
  - `NewtonRaphson_Err (*NonLinearRange_destroy)(const NonLinearRange *inRange);`
  - `NewtonRaphson_Err (*newton_raphson_solve)(const NonLinearRange *outRange, double *outRoot);`

### Secant（include/secant.h）
- 句柄与错误码
  - `typedef struct NonLinearScant *NonLinearScant;`
  - `typedef enum Secant_Err { SECANT_OK, ... }`;
- API：`extern const SecantAPI Secant;`
  - `Secant_Err (*NonLinearScant_create)(double (*f)(double), double x0, double x1, double tol, size_t max_iter, NonLinearScant *outScant, const char *name);`
  - `Secant_Err (*NonLinearScant_destroy)(const NonLinearScant *inScant);`
  - `Secant_Err (*secant_solve)(const NonLinearScant *outScant, double *outRoot);`

### Trapezoidal（include/Trapezoidal.h）
- 句柄与错误码
  - `typedef struct IntegrationApproximation *Trapezoidal;`
  - `typedef enum Trapezoidal_Err { TRAP_OK, ... }`;
- API：`extern const TrapezoidalIntegrationAPI TI;`
  - `Trapezoidal_Err (*trapezoidal_create)(double (*f)(double), double a, double b, size_t max_iter, Trapezoidal *outTrap, const char *name);`
  - `Trapezoidal_Err (*trapezoidal_integration)(const Trapezoidal *inTrap, double *outApproxIntegral);`
  - `Trapezoidal_Err (*trapezoidal_destroy)(Trapezoidal *inTrap);`

### Simpson（include/simpson.h）
- 句柄与错误码
  - `typedef struct IntegrationApproximation *Simpson;`
  - `typedef enum Simpson_Err { SIMPSON_OK, ... }`;
- API：`extern const SimpsonAPI SI;`
  - `Simpson_Err (*simpson_create)(double (*f)(double), double a, double b, size_t max_iter, Simpson *outSimpson, const char *name);`
  - `Simpson_Err (*simpson_integration)(const Simpson *inSimpson, double *outApproxIntegral);`
  - `Simpson_Err (*simpson_destroy)(Simpson *inSimpson);`

### Double Simpson（include/double_simpson.h）
- 句柄与错误码
  - `typedef struct DoubleIntegrationApproximation *DoubleSimpson;`
  - `typedef enum Simpson_Err { SIMPSON_OK, ... }`（与单变量定义名相同，注意命名冲突的潜在风险）
- API：`extern const DoubleSimpsonAPI DS;`
  - `Simpson_Err (*double_simpson_create)(double (*f)(double, double), double x_a, double x_b, double y_c, double y_d, size_t n, size_t m, DoubleSimpson *outDoubleSimpson, const char *name);`
  - `Simpson_Err (*double_simpson_integrate)(const DoubleSimpson *inDoubleSimpson, double *outApproxIntegral);`
  - `Simpson_Err (*double_simpson_destroy)(DoubleSimpson *inDoubleSimpson);`

### Successive Approximation（include/successive_approximation.h）
- 句柄与错误码
  - `typedef struct NonlinearSA *NonlinearSA;`
  - `typedef enum Succesive_Err { SA_OK, ... }`;
- API：`extern const SAAPI SA;`
  - `Succesive_Err (*NonlinearSA_create)(double (*g)(double), double x0, double tol, size_t max_iter, NonlinearSA *outSA, const char *name);`
  - `Succesive_Err (*nonlinear_sa_solve)(const NonlinearSA *outSA, double *outRoot);`
  - `Succesive_Err (*NonlinearSA_destroy)(const NonlinearSA *inSA);`

---
## 使用示例（选摘）
以下伪代码示意，具体可参见 tests/ 目录中的对应用例：
- Lagrange 插值
  1) 准备点集 `Point pts[] = { {x0,y0}, {x1,y1}, ... }`；
  2) `DataSet *ds = NULL; Lagrange.create_dataset(&ds, pts, n);`
  3) `double y = Lagrange.lagrange_interpolate(ds, xq);`
  4) `Lagrange.destroy_dataset(ds);`
- RK4 自适应积分
  1) 定义被积函数 `double f(double x, void *u)`；
  2) 配置 `AdaptiveConfig cfg = {1e-8, 1e-8, 20};`
  3) `IntegratorStatus st; double I = Integrator.rk4_adaptive(f, NULL, a, b, cfg, &st);`
- Newton–Raphson 求根
  1) `NonLinearRange range; NewtonRaphson.NonLinearRange_create(f, x0, tol, iters, &range, "case");`
  2) `double root; NewtonRaphson.newton_raphson_solve(range, &root);`
  3) `NewtonRaphson.NonLinearRange_destroy(range);`

---
## 测试说明
- 每个模块均有独立测试源 `tests/test_*.c`，与实现 `src/*.c` 一起编译为单独测试可执行文件，例如：
  - `Numerical_Analysis_tests_lagrange.exe`
  - `Numerical_Analysis_tests_integrator.exe`
  - 其它模块类似（bisection/newton/secant/simpson/...）。
- 运行方式：
  - 通过 CTest 统一运行：`ctest -C Debug -V`
  - 或直接双击/命令行运行对应的 `*.exe`
- 测试风格：
  - 每个 `tests/test_xxx.c` 通常包含 `run_all_tests()` 聚合全部测试用例；
  - 对常规、边界与异常场景进行覆盖（例如重复结点、无效区间、最大迭代等）。

---
## 开发建议与代码注释规范
- 头文件注释建议使用 Doxygen 风格，明确：
  - 功能描述、参数/返回值、错误码、所有权（资源创建/释放者）、前置条件/后置条件。
- 一致的 API 设计：
  - 统一采用 `extern const <Module>API <ModuleName>;` + 内部 `static` 实现进行初始化；
  - 推荐使用“指定字段初始化”提升可读性：`.field = impl_fn`；
  - 资源句柄类型统一为不透明指针，配套 `create/destroy`；
  - 错误码从 0 开始表示成功；
  - 注意跨模块的类型命名冲突（例如 Simpson_Err 在单/双 Simpson 中重复定义）。
- 构建与平台：
  - Windows 下建议开启 `/W4 /utf-8`（MSVC）或 `-Wall -Wextra`（GCC/Clang）并保持 `C11` 标准；
  - `-D_USE_MATH_DEFINES` 使 MSVC 提供 `M_PI/M_E`。

---
## 备注
- Lagrange 模块的具体改进建议参见 `docs/lagrange.md`（已在文档中列出不足与修改方向）。
- 本 README 内容依据仓库头文件与 CMake 配置整理，若后续接口有所调整，请同步更新本说明与相应测试。
