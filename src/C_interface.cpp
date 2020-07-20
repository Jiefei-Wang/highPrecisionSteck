#include <Rcpp.h>
#include <stdio.h>
#include <stdarg.h>
#include "mpfr.h"
#include <vector>
#include "tools.h"
//#include <gmp.h>
using namespace Rcpp;
using std::vector;



// [[Rcpp::export]]
NumericVector steck(int prec, NumericVector h, NumericVector l) {
	int n = LENGTH(h);
	//vector<__mpfr_struct*> P(n + 1);
	__mpfr_struct* P = new __mpfr_struct[n + 1];


	NumericVector res(n+1);

	__mpfr_struct* p = P;
	mpfr_t tmp;
	mpfr_inits2(prec, p, tmp, (mpfr_ptr)0);
	mpfr_set_d(p, 1, MPFR_RNDN);
	mpz_t binom;
	mpz_init(binom);
	DEBUG(Rprintf("Check\n"));
	for (int k = 1; k <= n; k++) {
		//__mpfr_struct* cur_p = new __mpfr_struct[1];
		__mpfr_struct* cur_p = P + k;
		mpfr_init2(cur_p, prec);
		mpfr_set_d(cur_p, 0, MPFR_RNDN);
		bool pos = true;
		for (int i = 1; i <= k; i++) {
			DEBUG(Rprintf("iter:%d,%d\n", k,i));
			DEBUG(Rprintf("%d,%d,%f,%f\n", k - i,k-1,h[k - i], l[k-1]));
			if (h[k - i] <= l[k-1]) {
				break;
			}
			// v-u
			mpfr_set_d(tmp, h[k-i]-l[k-1], MPFR_RNDN);
			DEBUG(Rprintf("res: %f\n", mpfr_get_d(tmp, MPFR_RNDN)));
			// tmp = (v-u)^i
			mpfr_pow_si(tmp, tmp, (long)i, MPFR_RNDN);
			DEBUG(Rprintf("res: %f\n", mpfr_get_d(tmp, MPFR_RNDN)));
			// binom = k choose i
			mpz_bin_uiui(binom, k, i);
			// tmp = binom * tmp = (k choose i) * (v-u)^i
			mpfr_mul_z(tmp, tmp, binom, MPFR_RNDN);
			DEBUG(Rprintf("res: %f\n", mpfr_get_d(tmp, MPFR_RNDN)));
			// tmp = binom * tmp * P_{k-i} = (k choose i) * (v-u)^i * P_{k-i}
			mpfr_mul(tmp, tmp, P+k - i, MPFR_RNDN);
			DEBUG(Rprintf("res: %f\n", mpfr_get_d(tmp, MPFR_RNDN)));
			// P_k = P_k + tmp
			if (pos) {
				mpfr_add(cur_p, cur_p, tmp, MPFR_RNDN);
			}
			else {
				mpfr_sub(cur_p, cur_p, tmp, MPFR_RNDN);
			}
			DEBUG(Rprintf("res: %f\n", mpfr_get_d(cur_p, MPFR_RNDN)));
			pos = !pos;
		}
		res[k] = mpfr_get_d(cur_p, MPFR_RNDN);
	}
	mpz_clear(binom);
	mpfr_clear(tmp);
	for (int i = 0; i < n + 1; i++) {
		mpfr_clear(P+i);
	}
	delete[] P;
  return res;
}


#include "timing.h"
// [[Rcpp::export]]
double test_pow(NumericVector x, int power, int prec) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
	__mpfr_struct* P = new __mpfr_struct[n];
	mpfr_t tmp;
	mpfr_init2(tmp, prec);
	for (int i = 0; i < n; i++) {
		mpfr_init2(P + i, prec);
		mpfr_set_d(P + i, x[i], MPFR_RNDN);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpfr_pow_si(tmp, P + i, power, MPFR_RNDN);
		res += mpfr_get_d(tmp, MPFR_RNDN);
	}
	toc();

	mpfr_clear(tmp);
	for (int i = 0; i < n; i++) {
		mpfr_clear(P + i);
	}
	delete[] P;

	tic();
	for (int i = 0; i < n; i++) {
		res -= pow(x[i],power);
	}
	toc();

	return res;
}



// [[Rcpp::export]]
double test_add(NumericVector x, int prec) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
	__mpfr_struct* P = new __mpfr_struct[n];
	mpfr_t tmp;
	mpfr_init2(tmp, prec);
	for (int i = 0; i < n; i++) {
		mpfr_init2(P + i, prec);
		mpfr_set_d(P + i, x[i], MPFR_RNDN);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpfr_add(tmp, P + i, P + i, MPFR_RNDN);
		res += mpfr_get_d(tmp, MPFR_RNDN);
	}
	toc();
	tic();
	for (int i = 0; i < n; i++) {
		res -= x[i] + x[i];
	}
	toc();


	mpfr_clear(tmp);
	for (int i = 0; i < n; i++) {
		mpfr_clear(P + i);
	}
	delete[] P;


	return res;
}



// [[Rcpp::export]]
double test_sub(NumericVector x, int prec) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
	__mpfr_struct* P = new __mpfr_struct[n];
	mpfr_t tmp;
	mpfr_init2(tmp, prec);
	for (int i = 0; i < n; i++) {
		mpfr_init2(P + i, prec);
		mpfr_set_d(P + i, x[i], MPFR_RNDN);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpfr_sub(tmp, P + i, P + i, MPFR_RNDN);
		res += mpfr_get_d(tmp, MPFR_RNDN);
	}
	toc();

	tic();
	for (int i = 0; i < n; i++) {
		res -= x[i] - x[i];
	}
	toc();


	mpfr_clear(tmp);
	for (int i = 0; i < n; i++) {
		mpfr_clear(P + i);
	}
	delete[] P;


	return res;
}

// [[Rcpp::export]]
double test_mul(NumericVector x, int prec) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
	__mpfr_struct* P = new __mpfr_struct[n];
	mpfr_t tmp;
	mpfr_init2(tmp, prec);
	for (int i = 0; i < n; i++) {
		mpfr_init2(P + i, prec);
		mpfr_set_d(P + i, x[i], MPFR_RNDN);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpfr_mul(tmp, P + i, P + i, MPFR_RNDN);
		res += mpfr_get_d(tmp, MPFR_RNDN);
	}
	toc();

	tic();
	for (int i = 0; i < n; i++) {
		res -= x[i]*x[i];
	}
	toc();

	mpfr_clear(tmp);
	for (int i = 0; i < n; i++) {
		mpfr_clear(P + i);
	}
	delete[] P;


	return res;
}


// [[Rcpp::export]]
double test_bin(NumericVector x, NumericVector y) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
		__mpz_struct* P = new __mpz_struct[n];
	for (int i = 0; i < n; i++) {
		mpz_init(P + i);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpz_bin_uiui(P+i, x[i], y[i]);
		res += mpz_get_d(P + i);
	}
	toc();


	for (int i = 0; i < n; i++) {
		mpz_clear(P + i);
	}
	delete[] P;


	return res;
}



// [[Rcpp::export]]
double test_mpz_add(NumericVector x) {
	int n = LENGTH(x);
	double res = 0;
	//vector<__mpfr_struct*> P(n + 1);
	__mpz_struct* P = new __mpz_struct[n];
	for (int i = 0; i < n; i++) {
		mpz_init(P + i);
		mpz_set_si(P + i, 2);
		mpz_pow_ui(P + i, P + i, (int)x[i]);
	}
	tic();
	for (int i = 0; i < n; i++) {
		mpz_mul_si(P + i, P + i, (int)x[i]);
		res += mpz_get_d(P + i);
	}
	toc();

	for (int i = 0; i < n; i++) {
		mpz_clear(P + i);
	}
	delete[] P;


	return res;
}