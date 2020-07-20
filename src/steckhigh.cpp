#include <Rcpp.h>
#include "mpfr.h"
#include "tools.h"
#include <algorithm>  
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector high_steck(NumericVector l, NumericVector h, int prec,NumericVector upperBound)
{
    R_xlen_t n = LENGTH(l);
    set_max_binomial(n);
    //NumericVector upperBound = computeUpperBound2(l,h);
    NumericVector res(n + 1);
    IntegerVector P_precision(n);
    res(0) = 1;
    P_precision(0) = prec;
    double prec_diff = ((double)prec-53)/n;
    prec_diff=0;
    if(prec_diff<0)
        prec_diff=0;
    for (int i = 1; i < n; i++)
    {
        P_precision(i)=prec - (int)floor(prec_diff*i);
    }
    mpfr_ptr mpfr_l = new  __mpfr_struct[n];
    mpfr_ptr mpfr_h = new  __mpfr_struct[n];
    mpfr_ptr P = new __mpfr_struct[n + 1];
    mpfr_init2(P, 2);
    mpfr_set_d(P, 1, MPFR_RNDN);
    int cur_precision = prec;
    for (int i = 1; i <= n; i++)
    {
        mpfr_init2(P + i, cur_precision);
        mpfr_set_d(P + i, 0, MPFR_RNDN);
        cur_precision=P_precision(i-1);
        int _i = i-1;
        mpfr_init2(mpfr_l + _i, prec);
        mpfr_init2(mpfr_h + _i, prec);
        mpfr_set_d(mpfr_l + _i,l[_i],MPFR_RNDN);
        mpfr_set_d(mpfr_h + _i,h[_i],MPFR_RNDN);
    }
    mpfr_t tmp;
    mpfr_init2(tmp, prec+n);

    mpz_t binom;
    mpz_init(binom);
    cur_precision = prec;
    //int count=0;
    for (unsigned int k = 1; k <= n; k++)
    {
        //int count1=0;
        cur_precision = P_precision(k-1);
        //Rprintf("%d,%d\n",k,cur_precision);
        mpfr_ptr cur_p = P + k;
        unsigned int _k = k - 1;
        // initial to k choose 1
        mpz_set_ui(binom, k);

        bool pos = true;
        for (unsigned int i = 1; i <= k; i++)
        {
            unsigned int _i = i - 1;
            if (h(_k - _i) - l(_k) <= 0)
            {
                break;
            }
            //count1++;
           /* double cur_log_a = C_lchoose(k, i);
            double cur_log_b = i * log2(h(_k - _i) - l(_k));
            double cur_log_ab = cur_log_a + cur_log_b;
            int tmp_prec = std::max(cur_precision + (int)ceil(cur_log_ab),53);*/
           mpfr_set_prec(tmp, cur_precision);
            // v-u
            //mpfr_set_d(tmp_raw, h[_k - _i] - l[_k], MPFR_RNDN);
            mpfr_sub(tmp, mpfr_h+_k - _i, mpfr_l+_k, MPFR_RNDN);
            // tmp = (v-u)^i
            mpfr_pow_si(tmp, tmp, (long)i, MPFR_RNDN);
            // tmp = binom * tmp = (k choose i) * (v-u)^i
            mpfr_mul_z(tmp, tmp, binom, MPFR_RNDN);
            // tmp = binom * tmp * P_{k-i} = (k choose i) * (v-u)^i * P_{k-i}
            mpfr_mul(tmp, tmp, P + k - i, MPFR_RNDN);

            //Rprintf("tmp:%f\n",mpfr_get_d(tmp, MPFR_RNDN));
            //mpfr_frac(tmp, tmp, MPFR_RNDN);
            //Rprintf("tmp:%f\n",mpfr_get_d(tmp, MPFR_RNDN));
            // P_k = P_k + tmp
            if (pos)
            {
                mpfr_add(cur_p, cur_p, tmp, MPFR_RNDN);
            }
            else
            {
                mpfr_sub(cur_p, cur_p, tmp, MPFR_RNDN);
            }
            //mpfr_frac(cur_p, cur_p, MPFR_RNDN);
            if (i != k)
            {
                mpz_mul_ui(binom, binom, k - i);
                mpz_fdiv_q_ui(binom, binom, i + 1);
            }
            DEBUG(Rprintf("res: %f\n", mpfr_get_d(cur_p, MPFR_RNDN)));
            pos = !pos;
        }
        //count = count1>count?count1:count;
        //Rprintf("%d,%d\n",count1,count);
        //Rprintf("%f\n",mpfr_get_d(cur_p, MPFR_RNDN));
        //mpfr_frac(cur_p, cur_p, MPFR_RNDN);
        /*if (mpfr_sgn(cur_p) <= 0)
        {
            mpfr_add_ui(cur_p, cur_p, 1, MPFR_RNDN);
        }*/
        res(k) = mpfr_get_d(cur_p, MPFR_RNDN);
        if(res(k)>upperBound(k-1)){
            //Rf_warning("Upperbound, Stop at k=%d", k);
            res(n) = -1;
            res(0) = k;
            break;
        }
        if(res(k)<res(k-1)*R::fmax2(h(k-1)-l(k-1),0)){
            //Rf_warning("lowerbound, Stop at k=%d", k);
            res(n) = -1;
            res(0) = k;
            break;
        }
        //if (log2(res[k]) < -50)
        //    break;
        // || (res[k] - res[k-1])>0
        /*if (res[k] < 0 || res[k] > 1)
        {
           Rf_warning("Stop at k=%d", k);
           break;
        }*/
        
    }
    mpz_clear(binom);
    mpfr_clear(tmp);
    for (int i = 0; i < n + 1; i++)
    {
        mpfr_clear(P + i);
    }
    delete[] P;

    for (int i = 0; i < n; i++)
    {
        mpfr_clear(mpfr_l + i);
        mpfr_clear(mpfr_h + i);
    }
    delete[] mpfr_l;
    delete[] mpfr_h;
    //List L = List::create(Named("p") = res);
    return res;
}


// [[Rcpp::export]]
NumericVector high_row(NumericVector l, NumericVector h, int prec, NumericVector upperBound)
{
    R_xlen_t n = LENGTH(l);
    set_max_binomial(n);
    //NumericVector upperBound = computeUpperBound2(l,h);
    NumericVector res(n+1);
    IntegerVector Q_precision(n);
    mpfr_ptr mpfr_l = new  __mpfr_struct[n];
    mpfr_ptr mpfr_h = new  __mpfr_struct[n];
    for(int i=0;i<n;i++){
        mpfr_init2(mpfr_l + i, prec);
        mpfr_init2(mpfr_h + i, prec);
        mpfr_set_d(mpfr_l + i,l[i],MPFR_RNDN);
        mpfr_set_d(mpfr_h + i,h[i],MPFR_RNDN);
    }
    mpfr_ptr Q = new __mpfr_struct[n];
    double prec_diff = ((double)prec-53)/n;
    prec_diff=0;
    if(prec_diff<0)
        prec_diff=0;
    for (int i = 0; i < n; i++)
    {
        Q_precision(i)= prec - (int)floor(prec_diff*i);
        mpfr_ptr curQ = Q + i;
        mpfr_init2(curQ, Q_precision(i));
        if(h(0) - l(i)>0){
            //mpfr_set_d(curQ, h(0) - l(i), MPFR_RNDN);
            mpfr_sub(curQ,mpfr_h,mpfr_l+i, MPFR_RNDN);
            mpfr_pow_si(curQ, curQ, (long)(i+1), MPFR_RNDN);
        }
        else{
            mpfr_set_d(curQ, 0, MPFR_RNDN);
        }
    }
    mpfr_t tmp;
    mpfr_init2(tmp, prec);

    mpz_t binom;
    mpz_init(binom);
    //int count=0;
    res(0) = 1;
    res(1) = mpfr_get_d(Q, MPFR_RNDN);
    for (unsigned int i = 1; i <= n-1; i++)
    {
        mpz_set_ui(binom, i + 1);
        //int count1=0;
        unsigned int _i = i - 1;
        mpfr_ptr Q_i = Q + _i;
        // initial to k choose 1
        for (unsigned int k = i+1; k <= n; k++)
        {
            unsigned int _k = k-1;
            mpfr_ptr Q_k = Q + _k;
            int cur_precision = Q_precision(_k);
            mpfr_set_prec(tmp, cur_precision);
            if (h(_i + 1) - l(_k) <= 0)
            {
                break;
            }
            // v-u
            //mpfr_set_d(tmp_raw, h(_i + 1) - l(_k), MPFR_RNDN);
            mpfr_sub(tmp, mpfr_h+_i + 1, mpfr_l+_k, MPFR_RNDN);
            // tmp = (v-u)^(k-i)
            mpfr_pow_si(tmp, tmp, (long)(k-i), MPFR_RNDN);
            // tmp = tmp * Q_i = b * (v-u)^(k-i) * Q_i
            mpfr_mul(tmp, tmp, Q_i, MPFR_RNDN);
            // tmp = binom * tmp = b * (v-u)^(k-i) * Q_i
            mpfr_mul_z(tmp, tmp, binom, MPFR_RNDN);
            // Q_k = Q_k - t
            mpfr_sub(Q_k, Q_k, tmp,MPFR_RNDN);

            
            mpz_mul_ui(binom, binom, k + 1);
            mpz_fdiv_q_ui(binom, binom, k + 1 - i);
        }
        int index = _i+2;
        //Rprintf("%d\n",index);
        res(index) =abs(mpfr_get_d(Q + _i +1, MPFR_RNDN));
        if(res(index)>upperBound(index-1)){
           // Rf_warning("Upperbound, Stop at i=%d", i);
            res(n) = -1;
            res(0) = i;
            break;
        }
        if(res(index)<res(index-1)*R::fmax2(h(index-1)-l(index-1),0)){
            //Rf_warning("lowerbound, Stop at i=%d", index);
            res(n) = -1;
            res(0) = index;
            break;
        }
        //if (log2(res[k]) < -50)
        //    break;
        // || (res[k] - res[k-1])>0
        /*if (res[k] < 0 || res[k] > 1)
        {
           Rf_warning("Stop at k=%d", k);
           break;
        }*/
        
    }

    mpz_clear(binom);
    mpfr_clear(tmp);
    for (int i = 0; i < n; i++)
    {
        mpfr_clear(Q + i);
        mpfr_clear(mpfr_l + i);
        mpfr_clear(mpfr_h + i);
    }
    delete[] Q;
    delete[] mpfr_l;
    delete[] mpfr_h;
    //List L = List::create(Named("p") = res);
    return res;
}