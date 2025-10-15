#include "Tools.h"
#include <vector>

void KEEPTIME(const char *str, clock_t &stamp){
    // eprintf("%s %lf\n", str, 1.0 * (clock() - stamp) / CLOCKS_PER_SEC);
    stamp = clock();
    return;
}