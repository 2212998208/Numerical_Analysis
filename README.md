# 数值分析算法库

> **项目简介**: 本项目是一个用 C 语言实现的数值分析算法库，专注于提供清晰、模块化且易于测试的数值计算工具。目前已实现的核心功能包括四阶龙格-库塔 (RK4) 数值积分、拉格朗日插值、牛顿插值和埃尔米特插值等。

---

## 1. 功能特性

- **数值积分**:
  - **固定步长 RK4**: 简单高效，适用于对计算速度要求高的场景。
  - **自适应步长 RK4**: 通过 Richardson 外推法进行误差估计，自动调整步长以达到预设精度，兼顾效率与准确性。

- **多项式插值**:
  - **拉格朗日插值**: 理论形式简单，易于理解。
  - **牛顿插值**: 具有承袭性，方便在增加新节点时重用之前的计算结果（已实现差商表）。
  - **埃尔米特插值**: 同时匹配函数值和导数值，提供更高阶的接触，插值结果更平滑。

- **模块化设计**:
  - 每个算法模块（积分、插值等）都通过统一的 `API` 结构体提供接口，代码解耦，易于扩展。
  - 采用信息隐藏，将数据结构 (`struct`) 的具体定义放在 `.c` 文件中，用户只能通过头文件中的指针和 API 函数进行操作。

- **独立的测试框架**:
  - 每个模块都有对应的测试文件 (`tests/test_*.c`)，可独立编译和运行，方便进行单元测试和回归测试。

---

## 2. 文件组织结构

```
Numerical_Analysis/
├── CMakeLists.txt           # CMake 构建配置
├── README.md                # 本文档
├── include/                 # 头文件目录 (公开API)
│   ├── integrator.h         # 数值积分 API
│   ├── lagrange.h           # 拉格朗日插值 API
│   ├── newton.h             # 牛顿插值 API
│   └── hermite.h            # 埃尔米特插值 API
├── src/                     # 源码实现
│   ├── integrator.c
│   ├── lagrange.c
│   ├── newton.c
│   └── hermite.c
├── tests/                   # 测试代码
│   ├── test_integrator.c
│   ├── test_lagrange.c
│   └── test_newton.c
│   └── test_hermite.c
└── docs/                    # 项目文档
    ├── Integrator_MindMap.svg
    └── lagrange.md
```

---

## 3. 如何编译与运行

本项目使用 `CMake` 进行构建管理。

### 编译步骤

1.  **确保已安装 CMake 和 C 编译器** (如 GCC, Clang, MSVC)。
2.  在项目根目录下创建一个构建目录：
    ```bash
    mkdir build
    cd build
    ```
3.  运行 CMake 生成构建系统（例如，为 Makefiles）：
    ```bash
    cmake ..
    ```
4.  编译项目：
    ```bash
    make
    ```
    或者在 Windows + Visual Studio 环境下，使用 `cmake --build .`。

### 运行测试

编译后，`build/bin` 目录下会生成各个模块的独立测试程序。例如，要运行拉格朗日插值的测试：

```bash
./bin/Numerical_Analysis_tests_lagrange
```

每个测试程序会输出详细的用例执行情况，包括计算值、期望值、误差和通过/失败状态。

---

## 4. 使用示例

所有公开接口都通过全局的 `const` API 结构体（如 `Integrator`, `Lagrange`）调用。

### 示例 1: 数值积分

计算函数 `f(x) = x^2` 在 `[0, 1]` 上的定积分。

```c
#include "integrator.h"
#include <stdio.h>

// 定义被积函数
double square_func(double x, void *user_data) {
    (void)user_data; // 本示例中未使用
    return x * x;
}

int main() {
    // 1. 使用固定步长积分
    double result_fixed = Integrator.rk4_fixed(square_func, NULL, 0.0, 1.0, 100);
    printf("固定步长积分结果: %f (期望值: ~0.333)\n", result_fixed);

    // 2. 使用自适应步长积分
    AdaptiveConfig config = { .abs_tol = 1e-7, .rel_tol = 1e-7, .max_iterations = 20 };
    IntegratorStatus status;
    double result_adaptive = Integrator.rk4_adaptive(square_func, NULL, 0.0, 1.0, config, &status);

    if (status == INTEGRATOR_OK) {
        printf("自适应积分结果: %f (期望值: ~0.333)\n", result_adaptive);
    } else {
        printf("自适应积分未达到指定精度。\n");
    }

    return 0;
}
```

### 示例 2: 拉格朗日插值

使用数据点 `(1,1), (2,4), (3,9)` 进行插值，计算 `x=2.5` 时的值。

```c
#include "lagrange.h"
#include <stdio.h>

int main() {
    Point points[] = { {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0} };
    size_t num_points = sizeof(points) / sizeof(points[0]);

    DataSet *dataset = NULL;
    Lagrange_Err err = Lagrange.create_dataset(&dataset, points, num_points);

    if (err == LAGRANGE_OK) {
        double x_interp = 2.5;
        double y_interp = Lagrange.lagrange_interpolate(dataset, x_interp);
        printf("在 x=%.2f 处的拉格朗日插值结果: %f (期望值: 6.25)\n", x_interp, y_interp);

        Lagrange.destroy_dataset(dataset);
    } else {
        fprintf(stderr, "创建数据集失败，错误码: %d\n", err);
    }

    return 0;
}
```

### 示例 3: 牛顿插值

同样使用 `(1,1), (2,4), (3,9)` 数据点，计算 `x=2.5` 时的值。

```c
#include "newton.h"
#include "lagrange.h" // 牛顿插值依赖于拉格朗日的数据结构
#include <stdio.h>

int main() {
    Point points[] = { {1.0, 1.0}, {2.0, 4.0}, {3.0, 9.0} };
    size_t num_points = sizeof(points) / sizeof(points[0]);

    // 1. 创建拉格朗日数据集（用于提供原始数据）
    DataSet *lagrange_data = NULL;
    Lagrange.create_dataset(&lagrange_data, points, num_points);

    // 2. 创建牛顿数据集（用于存储差商表）
    NewtonDataSet *newton_data = NULL;
    Newton.create_newton_dataset(&newton_data, num_points);

    if (lagrange_data && newton_data) {
        double x_interp = 2.5;
        double y_interp;

        // 3. 执行插值（此函数会内部计算差商表）
        Newton_Err err = Newton.newton_interpolate(&newton_data, &lagrange_data, x_interp, &y_interp);

        if (err == NEWTON_OK) {
            printf("在 x=%.2f 处的牛顿插值结果: %f (期望值: 6.25)\n", x_interp, y_interp);
        } else {
             fprintf(stderr, "牛顿插值失败，错误码: %d\n", err);
        }

        // 4. 清理
        Lagrange.destroy_dataset(lagrange_data);
        Newton.destroy_dataset(&newton_data);
    }

    return 0;
}
```

### 示例 4: 埃尔米特插值

使用点 `(1,1)` 且导数 `f'(1)=2`，和点 `(2,4)` 且导数 `f'(2)=4` 进行插值。

```c
#include "hermite.h"
#include <stdio.h>

int main() {
    // 数据点: x, y, y'
    const double x[] = {1.0, 2.0};
    const double y[] = {1.0, 4.0};
    const double dy[] = {2.0, 4.0}; // 对应 f(x)=x^2 的导数
    size_t num_points = 2;

    // 1. 创建数据集
    HermiteDataset dataset = Hermite.create_hermite_dataset(num_points, x, y, dy);

    // 2. 创建插值器（预计算系数）
    HermiteInterpolator interpolator = NULL;
    Hermite_Err err = Hermite.hermite_create_interpolator(&dataset, &interpolator);

    if (err == HERMITE_OK) {
        double x_eval = 1.5;
        double y_eval, dy_eval;

        // 3. 在 x=1.5 处求值
        Hermite.hermite_evaluate(&interpolator, x_eval, &y_eval, &dy_eval);

        printf("在 x=%.2f 处的 Hermite 插值:\n", x_eval);
        printf("  - 函数值: %f (期望值: 2.25)\n", y_eval);
        printf("  - 导数值: %f (期望值: 3.0)\n", dy_eval);
    } else {
        fprintf(stderr, "Hermite 插值器创建失败, 错误码: %d\n", err);
    }

    // 4. 清理
    Hermite.destroy_hermite_dataset(&dataset);
    Hermite.hermite_destroy_interpolator(&interpolator);

    return 0;
}
```

---

## 5. 代码风格与文档

- **注释**: 所有公开的头文件 (`.h`) 都遵循 **Doxygen** 风格进行注释，清晰地描述了每个函数、结构体、参数和返回值的用途。实现文件 (`.c`) 中也添加了必要的注释来解释关键算法和实现细节。
- **编码**: 所有代码文件和注释均采用 **UTF-8** 编码。
- **命名**: API 函数和结构体遵循模块化命名（如 `Lagrange.create_dataset`），实现函数使用 `_impl` 后缀以区分。

---

## 6. 未来扩展

- **增加更多算法**: 可以方便地扩展新的数值分析方法，如其他插值算法（样条插值）、线性代数求解器等。
- **完善错误处理**: 引入更详细的错误码和错误信息获取机制。
- **性能优化**: 对计算密集型部分进行性能分析和优化。
