#include <stdio.h>
#include <locale.h>
#ifdef _WIN32
#include <windows.h>
#endif

// 控制台 UTF-8 设置
static void set_console_utf8(void) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    setlocale(LC_ALL, ".UTF-8");
#else
    setlocale(LC_ALL, "C.UTF-8");
#endif
}


/* ======================= Demo main ============================= */
int main(void) {
    set_console_utf8();

    const double A[3][4] = {
        {2, 1, -1, 8},
        {-3, -1, 2, -11},
        {-2, 1, 2, -3}
    };
    const double *buf = &A[0][0];

    printf("A[1][1] = %d", (int)buf[5]);
};
