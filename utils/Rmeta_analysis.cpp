#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix getMetaP(NumericMatrix x) {
    int nrow = x.nrow(); int ncol = x.ncol();
    int npop = ncol / 3; // BETA, SE, and P
    NumericMatrix meta(nrow, 3);

    for (int i = 0; i < nrow; i++) {
        double B = 0;
        double S2 = 0;

        for (int j = 0; j < npop; j++) {
            double pj = x(i, 3*j+2);
            if (R_IsNA(pj)) {
                if (npop == 2) {
                    S2 = 0;
                    break;
                }
                continue;
            }
            double bj = x(i, 3*j);
            double s2j = x(i, 3*j+1) * x(i, 3*j+1);
            B += bj / s2j;
            S2 += 1 / s2j;
        }
        if (S2 == 0) {
            for (int k = 0; k < 3; k++) {meta(i, k) = NA_REAL;}
            continue;
        }

        meta(i, 0) = B / S2;
        meta(i, 1) = sqrt(1 / S2);
        double Z = meta(i, 0) / meta(i, 1);
        meta(i, 2) = R::pchisq(Z*Z, 1, 0, 0);
    }

    return meta;
}

Function quantile("quantile");

// [[Rcpp::export]]
double getLambda(NumericMatrix x) {
    int nrow = x.nrow();
    NumericVector chi(nrow);
    double lambda;

    for (int i = 0; i < nrow; i++) {
        double bi = x(i, 0);
        double s2i = x(i, 1) * x(i, 1);
        chi(i) = (bi / s2i) / sqrt(1 / s2i);
    }
    lambda = as<double>(quantile(chi*chi, 0.5, 1, 1, 7)) / 0.4549364;
    return lambda;
}


// [[Rcpp::export]]
NumericVector getMinP_meta(NumericMatrix x) {
    int nrow = x.nrow(); int ncol = x.ncol();
    int npop = ncol / 3; // BETA, SE, and P
    NumericVector metaZ_GC(nrow);
    NumericVector minP((npop+1) * 3 + 1); // (pops) + meta + (pops_GC) + meta_GC + meta_GC2 + (pops_lambda) + meta_lambda
    NumericVector minZ(npop+1);
    NumericVector lambda(npop+1);

    for (int j = 0; j <= npop; j++) {
        minP(j) = R_PosInf;
        if (j < npop) {
            int jj = (npop+1)*2 + j + 1;
            minP(jj) = getLambda(x(_, Range(3*j, 3*j+2)));
            lambda(j) = minP(jj) > 1 ? minP(jj) : 1;
        }
    }
    minP(npop*2 + 1) = R_PosInf; // meta_GC

    for (int i = 0; i < nrow; i++) {
        double B = 0;
        double S2 = 0;
        double S2_GC = 0;
        // double t = 0, u = 0;
        for (int j = 0; j < npop; j++) {
            int jj = 3*j;
            double pj = x(i, jj+2);
            if (R_IsNA(pj)) {continue;}
            double bj = x(i, jj);
            double s2j = x(i, jj+1) * x(i, jj+1);
            double zj = (bj / s2j) / sqrt(1 / s2j);
            if (pj < minP(j) && 0 < pj) {
                minP(j) = pj;
                minZ(j) = zj;
            }
            B += bj / s2j;
            S2 += 1 / s2j;
            S2_GC += (1 / s2j) / lambda(j);
            // t += bj*(1/x(i, jj+1))*(1/x(i, jj+1));
            // u += (1/x(i, jj+1))*(1/x(i, jj+1));
        }
        if (S2 == 0) {continue;} // If all p-vals are NA, ignore them.

        double Z = (B / S2) / sqrt(1 / S2);
        double Z_GC = (B / S2) / sqrt(1 / S2_GC);
        // double P = 2 * R::pnorm(fabs(Z), 0.0, 1.0, 0, 0);
        // double P_GC = 2 * R::pnorm(fabs(Z_GC), 0.0, 1.0, 0, 0);
        double P = R::pchisq(Z*Z, 1, 0, 0);
        double P_GC = R::pchisq(Z_GC*Z_GC, 1, 0, 0);
        metaZ_GC(i) = Z_GC;

        if (P < minP(npop) && 0 < P) {
            minP(npop) = P;
        }
        if (P_GC < minP(npop*2 + 1) && 0 < P_GC) {
            minP(npop*2 + 1) = P_GC;
            minZ(npop) = Z_GC;
        }
    }

    for (int j = 0; j < npop; j++) {
        double chi = minZ(j);
        int jj = (npop+1) + j;
        minP(jj) = R::pchisq(chi*chi / lambda(j), 1, 0, 0);
    }

    // double GC correction
    minP((npop+1) * 3) = as<double>(quantile(metaZ_GC*metaZ_GC, 0.5, 1, 1, 7)) / 0.4549364;
    lambda(npop) = minP((npop+1) * 3) > 1 ? minP((npop+1) * 3) : 1;
    double chi_GC = minZ(npop);
    minP((npop+1) * 2) = R::pchisq(chi_GC*chi_GC / lambda(npop), 1, 0, 0);

    return minP;
}


// [[Rcpp::export]]
NumericVector getMinP_genome(NumericMatrix x, int start, int end) {
    int nrow = x.nrow(); int ncol = x.ncol();
    int N = end - start + 1;
    int nchr = nrow / N;

    NumericMatrix minP(N, ncol);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < ncol; j++) {
            minP(i, j) = R_PosInf;
            for (int c = 0; c < nchr; c++) {
               int ii = c * N + i;
               if (x(ii, j) < minP(i, j) && 0 < x(ii, j)) {
                   minP(i, j) = x(ii, j);
               }
            }
        }
    }

    return minP;
}

