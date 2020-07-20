
# highPrecisionSteck

<!-- badges: start -->
<!-- badges: end -->

The goal of highPrecisionSteck is to ...

## Installation

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
For showing the result
``` r
final <- lapply(record, function(x) colMeans(x))
df <- data.frame(matrix(unlist(final), nrow = length(final), byrow = T))
colnames(df) <- colnames(record[[1]])
df
```

