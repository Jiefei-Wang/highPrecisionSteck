#include <vector>
#include "tools.h"
#include <algorithm>  
#include <cmath>
using std::vector;

vector<double> log_prod;
// [[Rcpp::export]]
void set_max_binomial(int n)
{
    if (log_prod.size() == 0)
    {
        log_prod.push_back(0);
    }
    if (log_prod.size() < (size_t)n + 1)
    {
        log_prod.reserve(n + 1);
        for (int i = log_prod.size(); i < n + 1; i++)
        {
            log_prod.push_back(log_prod[i - 1] + log2(i));
        }
    }
}
// [[Rcpp::export]]
std::vector<double> get_log_prod()
{
    return (log_prod);
}
// [[Rcpp::export]]
double C_log_prod(int n){
    return log_prod[n];
}
// [[Rcpp::export]]
double C_log_beta(int a,int b){
    return C_log_prod(a-1)+C_log_prod(b-1)-C_log_prod(a+b-1);
}


// [[Rcpp::export]]
double C_lchoose(int n, int k)
{
    return log_prod[n] - log_prod[k] - log_prod[n - k];
}
