#include <chrono>
#include <ctime>
#include <functional>
#include "Rinternals.h"

using namespace std;
clock_t START_TIMER;

void tic()
{
    START_TIMER = clock();
}

void toc()
{
    auto end = std::clock();
    Rprintf("Elapsed time: %.6g ns\n", (end - START_TIMER) / (double)CLOCKS_PER_SEC * 1000);
}
