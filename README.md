
# highPrecisionSteck

The goal of highPrecisionSteck is to use the high precision library to compute the steck's determinant and probability.

## Installation instruction
### Dependencies
Before installing the package, you must have GNU Multiple Precision Arithmetic Library(gmp) and Floating-Point Reliable Library(mpfr) installed. The source code can be found at 
[gmp](https://gmplib.org/) and [mpfr](https://www.mpfr.org/). If you are using Windows, you can also install the precompiled version from [msys2](https://www.msys2.org/).

### Build the package
The package using the variable `MPFR_LIB` and `GMP_LIB` to find gmp and mpfr library respectively. For setting these variables, you can either manually edit them in the file `src/Makevars`, or you can add them to your global environment. If you are on Windows and using `msys2`, you can set these variables to "path-to-msys2/mingw32" or "path-to-msys2/mingw64" depending on your system.

After configuring these variables, the package can be installed as usual.


## Simulation
This is the code for running the simulation in the paper:
```r
library(highPrecisionSteck)
N_list <- c(50,100,200,500,1000,2000,5000,10000,15000)
sim_num <- 100

record = list()
for(i in seq_along(N_list)){
    n <-  N_list[i]
    message("Simulating steck for n = ",n)
    record[[i]]  <- run_simulation(n, sim_num, progress = TRUE)
}
```
For showing the result, the code below compute the average precision, time for each algorithm.
``` r
final <- lapply(record, function(x) colMeans(x))
df <- data.frame(matrix(unlist(final), nrow = length(final), byrow = T))
colnames(df) <- colnames(record[[1]])
df[c("n", "high_steck_prec", "high_row_prec", "high_steck_error", "high_row_error", 
     "time_rational_steck", "time_rational_row", "timing_high_steck", 
     "timing_high_row")]
```

