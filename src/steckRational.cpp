#include <Rcpp.h>
#include "mpfr.h"
#include "tools.h"
#include <algorithm>  
#include <vector>
using namespace Rcpp;
using std::vector;

vector<mpq_ptr> tmp_vector;
void mpq_power(mpq_ptr res,mpq_ptr op, int i){
    if(i==0){
        mpq_set_d(res, 1);
    }else{
        if(i==1){
           mpq_set(res,op);
        }else{

        }
    }
    switch(i){
        case 0:
            mpq_set_d(res, 1);
            return;
        case 1:
            mpq_set(res,op);
            return;
        case 2:
            mpq_mul(res,op,op);
            return;
    }
    double n = floor(log2(i)) + 1;
    tmp_vector.reserve(n);
    for(int j = tmp_vector.size(); j<n;j++){
        mpq_ptr a = new __mpq_struct;
        mpq_init(a);
        tmp_vector.push_back(a);
    }
    mpq_set(tmp_vector[0],op);
    for(int j = 1; j<n;j++){
        mpq_mul(tmp_vector[j],tmp_vector[j-1],tmp_vector[j-1]);
    }
    mpq_set_d(res,1);
    while(i!=0){
        int k= floor(log2(i));
        mpq_mul(res,res,tmp_vector[k]);
        i = i - (int)pow(2,k);
    }
}
// [[Rcpp::export]]
NumericVector rational_steck(NumericVector l, NumericVector h)
{
    R_xlen_t n = LENGTH(l);
    set_max_binomial(n);
    NumericVector res(n + 1);
    res(0) = 1;
    mpq_ptr rl = new  __mpq_struct[n];
    mpq_ptr rh = new  __mpq_struct[n];
    mpq_ptr P = new __mpq_struct[n + 1];
    mpq_init(P);
    mpq_set_d(P, 1);
    for (int i = 1; i <= n; i++)
    {
        mpq_init(P + i);
        mpq_set_d(P + i, 0);
        int _i = i-1;
        mpq_init(rl + _i);
        mpq_init(rh + _i);
        mpq_set_d(rl + _i, l(_i));
        mpq_set_d(rh + _i, h(_i));
        //Rprintf("%f,%f\n",mpq_get_d(rl+_i),mpq_get_d(rh+_i));
    }


    mpq_t tmp;
    mpq_init(tmp);
    mpq_t binom;
    mpq_init(binom);
    for (unsigned int k = 1; k <= n; k++)
    {
        //Rprintf("%d\n",k);
        mpq_ptr cur_p = P + k;
        unsigned int _k = k - 1;
        // initial to k choose 1
        mpq_set_d(binom, k);

        bool pos = true;
        for (unsigned int i = 1; i <= k; i++)
        {
            unsigned int _i = i - 1;
            if (h(_k - _i) - l(_k) <= 0)
            {
                break;
            }
            // v-u
            mpq_sub(tmp, rh+_k - _i,rl+_k);
            // tmp = (v-u)^i
            mpq_power(tmp,tmp,i);
            // tmp = binom * tmp = (k choose i) * (v-u)^i
            mpq_mul(tmp, tmp, binom);
            // tmp = binom * tmp * P_{k-i} = (k choose i) * (v-u)^i * P_{k-i}
            mpq_mul(tmp, tmp, P + k - i);
            if (pos)
            {
                mpq_add(cur_p, cur_p, tmp);
            }
            else
            {
                mpq_sub(cur_p, cur_p, tmp);
            }
            if (i != k)
            {   
                mpq_set_d(tmp,k-i);
                mpq_mul(binom, binom, tmp);
                mpq_set_d(tmp,i + 1);
                mpq_div(binom, binom, tmp);
            }
            //DEBUG(Rprintf("res: %f\n", mpfr_get_d(cur_p, MPFR_RNDN)));
            pos = !pos;
        }
        res(k) = mpq_get_d(cur_p);
    }
    mpq_clear(binom);
    mpq_clear(tmp);
    for (int i = 0; i < n + 1; i++)
    {
        mpq_clear(P + i);
    }
    delete[] P;
    return res;
}



// [[Rcpp::export]]
NumericVector rational_row(NumericVector l, NumericVector h)
{
    R_xlen_t n = LENGTH(l);
    set_max_binomial(n);
    NumericVector res(n+1);
    mpq_ptr mpq_l = new  __mpq_struct[n];
    mpq_ptr mpq_h = new  __mpq_struct[n];
    for(int i=0;i<n;i++){
        mpq_init(mpq_l + i);
        mpq_init(mpq_h + i);
        mpq_set_d(mpq_l + i,l[i]);
        mpq_set_d(mpq_h + i,h[i]);
    }
    mpq_ptr Q = new __mpq_struct[n];
    for (int i = 0; i < n; i++)
    {
        mpq_ptr curQ = Q + i;
        mpq_init(curQ);
        if(h(0) - l(i)>0){
            //curQ = h(0) - l(i);
            mpq_sub(curQ,mpq_h,mpq_l+i);
            mpq_power(curQ, curQ, i+1);
        }
        else{
            mpq_set_d(curQ, 0);
        }
    }
    mpq_t tmp;
    mpq_init(tmp);

    mpq_t binom;
    mpq_init(binom);
    //int count=0;
    res(0) = 1;
    res(1) = mpq_get_d(Q);
    for (unsigned int i = 1; i <= n-1; i++)
    {
        mpq_set_d(binom, i + 1);
        //int count1=0;
        unsigned int _i = i - 1;
        mpq_ptr Q_i = Q + _i;
        // initial to k choose 1
        for (unsigned int k = i+1; k <= n; k++)
        {
            unsigned int _k = k-1;
            mpq_ptr Q_k = Q + _k;
            if (h(_i + 1) - l(_k) <= 0)
            {
                break;
            }
            // v-u
            //tmp = h(_i + 1) - l(_k)
            mpq_sub(tmp, mpq_h+_i + 1, mpq_l+_k);
            // tmp = (v-u)^(k-i)
            mpq_power(tmp, tmp, k-i);
            // tmp = tmp * Q_i = b * (v-u)^(k-i) * Q_i
            mpq_mul(tmp, tmp, Q_i);
            // tmp = binom * tmp = b * (v-u)^(k-i) * Q_i
            mpq_mul(tmp, tmp, binom);
            // Q_k = Q_k - t
            mpq_sub(Q_k, Q_k, tmp);

            
            mpq_set_d(tmp,k + 1);
            mpq_mul(binom, binom, tmp);
            mpq_set_d(tmp,k + 1 - i);
            mpq_div(binom, binom, tmp);
        }
        res(i+1) =abs(mpq_get_d(Q + i));
    }

    mpq_clear(binom);
    mpq_clear(tmp);
    for (int i = 0; i < n; i++)
    {
        mpq_clear(Q + i);
        mpq_clear(mpq_l + i);
        mpq_clear(mpq_h + i);
    }
    delete[] Q;
    delete[] mpq_l;
    delete[] mpq_h;
    return res;
}