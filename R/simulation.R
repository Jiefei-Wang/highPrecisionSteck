
run_simulation <- function(N, rep=100, precision = 128, progress = FALSE){
    record <- c()
    if(progress)
        pb = txtProgressBar(min = 0, max = rep, initial = 0, style = 3) 
    for(i in 1:rep){
        x <-runif(N)
        stat<- BJStat(x)
        lc <- BJLocalCritical(stat,N)
        l <- lc$l
        h <- lc$h
        
        timing_rational1 <- system.time({
            rational_steck <- compute_rational_steck(l,h)[N+1]
        })[3]
        timing_rational2 <- system.time({
            rational_row <- compute_rational_row(l,h)[N+1]
        })[3]
        
        if(rational_steck!=rational_row){
            stop("Rational results is not consistent")
        }
        
        timing_high1 <- system.time({
            high_steck <- compute_high_steck(l,h,precision)
        })[3]
        timing_high2 <- system.time({
            high_row <- compute_high_row(l,h,precision)
        })[3]
        
        diff_high <- abs(high_steck$pval - rational_steck)
        diff_row <- abs(high_row$pval - rational_steck)
        
        record<-rbind(record,
                      c(N,stat,
                        rational_steck,
                        high_steck$pval,high_row$pval,
                        high_steck$prec,high_row$prec,
                        diff_high,diff_row,
                        timing_rational1,timing_rational2,
                        timing_high1,timing_high2))
        if(progress)
            setTxtProgressBar(pb,i)
    }
    if(progress)
        close(pb)
    colnames(record) <- c("n","stat",
                          "true_pval",
                          "high_steck_pval","high_row_pval",
                          "high_steck_prec","high_row_prec",
                          "high_steck_error","high_row_error",
                          "time_rational_steck","time_rational_row",
                          "timing_high_steck","timing_high_row")
    record
}