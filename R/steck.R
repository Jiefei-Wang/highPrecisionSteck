compute_high_prec <- function(func,l,h,prec){
    n <- length(l)
    upperbound <- compute_upper_bound(l,h,50)
    new_precision = prec
    old_res <- 0
    new_res <- 1
    ## if an out-of-bound result has been found,
    ## new_res will be -1.
    while (new_res==-1 ||
           abs(new_res-old_res)>0.001) {
        prec <- new_precision
        old_res <- new_res
        pval <- func(l,h,prec,upperbound)
        new_res <- pval[n+1]
        new_precision <- new_precision * 2
        if(prec > 1024*8){
            stop()
        }
    }
    list(pval = 1 - new_res, prec=prec)
}

compute_high_steck <- function(l,h,prec = 128){
    compute_high_prec(high_steck,l,h,prec)
}
compute_high_row <- function(l,h,prec = 128){
    compute_high_prec(high_row,l,h,prec)
}
compute_rational_steck <- function(l,h){
    1 - rational_steck(l,h)
}
compute_rational_row <- function(l,h){
    1 - rational_row(l,h)
}





