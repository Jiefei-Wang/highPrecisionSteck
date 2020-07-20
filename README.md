# Installation instruction
## Dependencies
Before installing the package, you must have GNU Multiple Precision Arithmetic Library(gmp) and Floating-Point Reliable Library(mpfr) installed. The source code can be found at 
[gmp](https://gmplib.org/) and [mpfr](https://www.mpfr.org/). If you are using Windows, you can also install the precompiled version from [msys2](https://www.msys2.org/).

## Build the package
The package using the variable `MPFR_LIB` and `GMP_LIB` to find gmp and mpfr library respectively. For setting these variables, you can either manually edit them in the file `src/Makevars`, or you can add them to your global environment. If you are on Windows and using `msys2`, you can set these variables to "path-to-msys2/mingw32" or "path-to-msys2/mingw64" depending on your system.

After configuring these variables, the package can be installed as usual.


# Simulation
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

# Example
Suppose your working directory is at the repository root. Here is the code for generating the empirical CDF.
```
library(sROC)
mydata <- read.csv("data\\Yang_YoungOnsetHypertension_Illumina550_SBAS_rawpv.csv")
theGene <- "NAV2"
gene_p <- mydata[mydata$Symbol==theGene,"CLR_N_BMI_pv"]
x <- seq(0,1,by=0.01)
est_CDF <- kCDF(gene_p,xgrid = x)$Fhat

plot(x,est_CDF, xlab ="", ylab = "", type="l", 
     axes =FALSE,xlim = c(0,1),ylim=c(-0.05,1),col ="gray",lwd=2)
axis(1,at=c(0,0.5,1),labels=c("0", "0.5", "1"), cex.axis = 1)
axis(2,at=c(0,0.5,1),labels=c("0", "0.5", "1"), cex.axis = 1)
axis(1,at=gene_p,labels =FALSE,tck = 0.05)
lines(x, x, lty= "dashed",lwd=2)
```

For computing the Berk-Jone statistic and its p-value.
```
mydata <- read.csv("data\\Yang_YoungOnsetHypertension_Illumina550_SBAS_rawpv.csv")
theGene <- "NAV2"
gene_p <- mydata[mydata$Symbol==theGene, "CLR_N_BMI_pv"]
n <- length(gene_p)
stat <- BJStat(gene_p,indexU = NULL)

lc <- BJLocalCritical(stat,n)
l <- lc$l
h <- rep(1, n)
pvalue <- compute_high_steck(l,h)$pval
```
