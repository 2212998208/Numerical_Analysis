# Numerical_Analysis 项目说明

> 本仓库：数值分析常用算法实现，包括 Runge–Kutta 四阶 (RK4) 数值积分、Lagrange/牛顿/埃尔米特插值等，支持模块化测试与扩展。

---
## 1. 功能简介
- **数值积分**：基于 RK4 方法，支持固定步长与自适应步长（双网格 Richardson 误差估计）。
- **插值法**：实现 Lagrange、Newton、Hermite 多项式插值，支持异常处理与多种边界情况。
- **测试框架**：每个模块均有独立测试入口，便于单元测试与回归。

---
## 2. 文件组织结构
```
Numerical_Analysis/
├── CMakeLists.txt           # CMake 构建配置
├── main.c                   # 程序入口
├── include/                 # 头文件目录
│   ├── integrator.h         # 积分 API
│   ├── lagrange.h           # Lagrange 插值 API
│   ├── newton.h             # Newton 插值 API
│   ├── hermite.h            # Hermite 插值 API
│   ├── test_integrator.h    # 积分测试声明
│   ├── test_lagrange.h      # Lagrange 测试声明
│   └── test_newton.h        # Newton 测试声明
├── src/                     # 源码实现
│   ├── integrator.c         # 积分实现
│   ├── lagrange.c           # Lagrange 实现
│   ├── newton.c             # Newton 实现
│   └── hermite.c            # Hermite 实现
├── tests/                   # 各模块测试实现
│   ├── test_integrator.c    # 积分测试
│   ├── test_lagrange.c      # Lagrange 测试
│   └── test_newton.c        # Newton 测试
├── docs/                    # 文档
│   ├── Integrator_MindMap.svg
│   └── lagrange.md          # Lagrange 模块改进建议
```

---
## 3. CMakeLists.txt 组织方式
- 采用 CMake 管理多模块编译，支持单独或全部测试可执行文件生成。
- 头文件统一放在 include/，源码在 src/，测试在 tests/。
- 每个测试文件（如 test_lagrange.c）可独立编译运行。
- 示例：
  - `add_executable(Numerical_Analysis_tests_lagrange tests/test_lagrange.c src/lagrange.c)`
  - `add_executable(Numerical_Analysis main.c src/integrator.c ...)`

---
## 4. 主要 .h/.c 文件说明与用法
### integrator.h / integrator.c
- 提供 RK4 积分 API：
  - `double Integrator.rk4_fixed(...)` 固定步长
  - `double Integrator.rk4_adaptive(...)` 自适应步长
- 结构体：`AdaptiveConfig`, `IntegratorStatus`
- 用法见 test_integrator.c

### lagrange.h / lagrange.c
- 提供 Lagrange 插值 API：
  - `Lagrange.create_dataset`, `Lagrange.lagrange_interpolate`, `Lagrange.destroy_dataset` 等
- 结构体：`Point`, `DataSet`, `Lagrange_Err`
- 用法见 test_lagrange.c
- 改进建议详见 docs/lagrange.md

### newton.h / newton.c
- 提供 Newton 插值法相关接口
- 用法见 test_newton.c

### hermite.h / hermite.c
- 提供 Hermite 插值法相关接口

### test_xxx.h / test_xxx.c
- 每个模块有独立测试头文件（声明 int run_all_tests(void)）和实现文件
- 测试实现文件包含 main()，可独立运行
- 测试用例覆盖常规、边界、异常情况

---
## 5. 测试说明
- **运行方式**：
  - 通过 CMake 生成的可执行文件（如 Numerical_Analysis_tests_lagrange.exe）直接运行
  - Windows 下建议命令行运行，自动设置 UTF-8 控制台编码
- **测试内容**：
  - 积分模块：验证 x^2、sin、exp 等函数的积分精度
  - Lagrange 插值：线性、二次、常数、重复节点、乱序节点等多种情况
  - 其他插值法：见对应 test_xxx.c
- **测试输出**：
  - 每个用例输出计算值、期望值、误差、通过/失败标记
  - 统计通过用例数

---
## 6. 注释与代码风格建议
- 头文件建议为每个结构体、函数添加 Doxygen 风格注释，说明参数、返回值、所有权、错误码等
- 结构体命名简洁明了，接口参数尽量 const、简洁
- 内存管理接口需健壮，避免未初始化和重复释放

---
## 7. 参考与扩展
- Lagrange 插值模块改进建议详见 docs/lagrange.md
- 可扩展更多数值分析算法与测试
