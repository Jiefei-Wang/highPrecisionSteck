bound_func <- function(x, k){
    pbeta(x[seq_len(k)],1:k,k:1)
}

compute_upper_bound_k <- function(l,h,k){
    min(bound_func(h,k) - bound_func(l,k)) 
}

## compute the upper bound
compute_upper_bound <- function(l,h,check_number){
    n <- length(l)
    res <- rep(1,n)
    check_point <- round(seq(from = 1, to= n, length.out = check_number))
    for(k in check_point){
        ## We add a small value to handle the inaccuracy of double-precision number
        res[k] <- compute_upper_bound_k(l,h,k) + 10^-6
    }
    res
}