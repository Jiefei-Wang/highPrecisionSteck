#' @param x data
#' @param indexL lower side index
#' @param indexU upper side index
BJStat<-function(x,indexL=seq_along(x),indexU=seq_along(x)){
    n <- length(x)
    sx <- sort(x)
    sx[sx == 0] <- min(10 ^ -6, sx[sx != 0])
    sx[sx == 1] <- max(1 - 10 ^ -6, sx[sx != 1])
    BJLevel(
        n = n,
        x = x,
        sx = sx,
        indexU = indexU ,
        indexL = indexL
    )
}

BJPlusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    vapply(seq_along(x)[index], function(x)
        pbeta(sx[x], x, n - x + 1),numeric(1))
}
BJMinusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    1 - BJPlusLevel(n, x, sx, index)
}
BJLevel <- function(n, x, sx, indexL, indexU) {
    BJPlus <- min(BJPlusLevel(n, x, sx,indexL),1)
    BJMinus <- min(BJMinusLevel(n, x, sx,indexU),1)
    min(BJPlus, BJMinus)
}



## given the statistics, compute the local critical value
BJLocalCritical<-function(stat,n){
    l=vapply(seq_len(n),function(x)qbeta(stat,x,n-x+1),numeric(1))
    h=vapply(seq_len(n),function(x)qbeta(1 - stat,x,n-x+1),numeric(1))
    list(l =l,h= h)
}