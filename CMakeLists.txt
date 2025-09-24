cmake_minimum_required(VERSION 3.15)
project(Numerical_Analysis C)

set(CMAKE_C_STANDARD 11)

set(INTEGRATOR
    main.c
    src/integrator.c
)
set(LAGRANGE
    main.c
    src/lagrange.c
)
set(NEWTON
    main.c
    src/newton.c
)
set(HERMITE
    main.c
    src/hermite.c
)
set(BISECTION
    main.c
    src/bisection.c
)
set(NEWTONRAPHSON
    main.c
    src/newton_raphson.c
)
set(SECANT
    main.c
    src/secant.c
)
set(SA
    main.c
    src/successive_approximation.c
)
set(TRAPEZOIDAL
    src/Trapezoidal.c
    main.c
)
set(SI
    main.c
    src/simpson.c
)
set(DS
    main.c
    src/double_simpson.c
)
set(EULER
    main.c
    src/euler.c
)
set(RK4
    main.c
    src/RK4.c
)
set(TESTS_INTEGRATOR
    src/integrator.c
    tests/test_integrator.c
)
set(TESTS_LAGRANGE
    src/lagrange.c
    tests/test_lagrange.c
)
set(TESTS_NEWTON
    src/lagrange.c
    src/newton.c
    tests/test_newton.c
)
set(TESTS_HERMITE
    src/hermite.c
    tests/test_hermite.c
)
set(TESTS_BISECTION
    src/bisection.c
    tests/test_bisection.c
)
set(TESTS_NEWTONRAPHSON
    src/newton_raphson.c
    tests/test_newtonraphson.c
)
set(TESTS_SECANT
    src/secant.c
    tests/test_secant.c
)
set(TESTS_SA
    src/successive_approximation.c
    tests/test_successiveapproximation.c
)
set(TESTS_TI
    src/Trapezoidal.c
    tests/test_trapezoidal.c
)
set(TESTS_SI
    src/simpson.c
    tests/test_simpson.c
)
set(TESTS_DS
    src/double_simpson.c
    tests/test_doublesimpson.c
)
set(TESTS_EULER
    src/euler.c
    tests/test_euler.c
)
set(TESTS_RK4
    src/RK4.c
    tests/test_RK4.c
)
#===================================================================
add_executable(Numerical_Analysis
               ${INTEGRATOR}
               ${LAGRANGE}
               ${NEWTON}
               ${HERMITE}
               ${BISECTION}
               ${NEWTONRAPHSON}
               ${SECANT}
               ${SA}
               ${TRAPEZOIDAL}
               ${SI}
               ${DS}
               ${EULER}
               ${RK4}
)

target_include_directories(Numerical_Analysis PRIVATE include)

#===================================================================

# 使 MSVC 获得 M_PI / M_E
add_definitions(-D_USE_MATH_DEFINES)

if (MSVC)
    target_compile_options(Numerical_Analysis PRIVATE /W4 /utf-8 /O2)
else()
    target_compile_options(Numerical_Analysis PRIVATE -Wall -Wextra -O3)
endif()

# 生成目录
set_target_properties(Numerical_Analysis PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

enable_testing()
#===================================================================
# 测试integrator
add_executable(Numerical_Analysis_tests_integrator
               ${TESTS_INTEGRATOR})
target_include_directories(Numerical_Analysis_tests_integrator PRIVATE include)

add_test(NAME Numerical_Analysis_tests_integrator COMMAND Numerical_Analysis_tests_integrator)

#===================================================================

#===================================================================
# 测试lagrange
add_executable(Numerical_Analysis_tests_lagrange
        ${TESTS_LAGRANGE})
target_include_directories(Numerical_Analysis_tests_lagrange PRIVATE include)

add_test(NAME Numerical_Analysis_tests_lagrange COMMAND Numerical_Analysis_tests_lagrange)

#===================================================================

#===================================================================
# 测试newton
add_executable(Numerical_Analysis_tests_newton
        ${TESTS_NEWTON})
target_include_directories(Numerical_Analysis_tests_newton PRIVATE include)
add_test(NAME Numerical_Analysis_tests_newton COMMAND Numerical_Analysis_tests_newton)
#===================================================================

#===================================================================
# 测试hermite
add_executable(Numerical_Analysis_tests_hermite
            ${TESTS_HERMITE})
target_include_directories(Numerical_Analysis_tests_hermite PRIVATE include)
add_test(NAME Numerical_Analysis_tests_hermite COMMAND Numerical_Analysis_tests_hermite)
#===================================================================

#===================================================================
# 测试bisection
add_executable(Numerical_Analysis_tests_bisection
            ${TESTS_BISECTION})
target_include_directories(Numerical_Analysis_tests_bisection PRIVATE include)
add_test(NAME Numerical_Analysis_tests_bisection COMMAND Numerical_Analysis_tests_bisection)
#===================================================================

#===================================================================
# 测试newtonraphson
add_executable(Numerical_Analysis_tests_newtonraphson
            ${TESTS_NEWTONRAPHSON})
target_include_directories(Numerical_Analysis_tests_newtonraphson PRIVATE include)
add_test(NAME Numerical_Analysis_tests_newtonraphson COMMAND Numerical_Analysis_tests_newtonraphson)
#===================================================================

#===================================================================
# 测试secant
add_executable(Numerical_Analysis_tests_secant
            ${TESTS_SECANT})
target_include_directories(Numerical_Analysis_tests_secant PRIVATE include)
add_test(NAME Numerical_Analysis_tests_secant COMMAND Numerical_Analysis_tests_secant)
#===================================================================

#===================================================================
# 测试successive approximation
add_executable(Numerical_Analysis_tests_sa
            ${TESTS_SA})
target_include_directories(Numerical_Analysis_tests_sa PRIVATE include)
add_test(NAME Numerical_Analysis_tests_sa COMMAND Numerical_Analysis_tests_sa)
#===================================================================

#===================================================================
# 测试trapezoidal
add_executable(Numerical_Analysis_tests_ti
            ${TESTS_TI})
target_include_directories(Numerical_Analysis_tests_ti PRIVATE include)
add_test(NAME Numerical_Analysis_tests_ti COMMAND Numerical_Analysis_tests_ti)
#===================================================================

#===================================================================
# 测试simpson
add_executable(Numerical_Analysis_tests_si
            ${TESTS_SI})
target_include_directories(Numerical_Analysis_tests_si PRIVATE include)
add_test(NAME Numerical_Analysis_tests_si COMMAND Numerical_Analysis_tests_si)
#===================================================================

#===================================================================
# 测试double simpson
add_executable(Numerical_Analysis_tests_ds
            ${TESTS_DS})
target_include_directories(Numerical_Analysis_tests_ds PRIVATE include)
add_test(NAME Numerical_Analysis_tests_ds COMMAND Numerical_Analysis_tests_ds)
#===================================================================

#===================================================================
# 测试euler
add_executable(Numerical_Analysis_tests_euler
            ${TESTS_EULER})
target_include_directories(Numerical_Analysis_tests_euler PRIVATE include)
add_test(NAME Numerical_Analysis_tests_euler COMMAND Numerical_Analysis_tests_euler)
#===================================================================

#===================================================================
# 测试RK4
add_executable(Numerical_Analysis_tests_rk4
            ${TESTS_RK4})
target_include_directories(Numerical_Analysis_tests_rk4 PRIVATE include)
add_test(NAME Numerical_Analysis_tests_rk4 COMMAND Numerical_Analysis_tests_rk4)
#===================================================================

