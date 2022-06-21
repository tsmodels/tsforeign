#include "filters.h"
#include <Rcpp.h>
using namespace Rcpp;


inline arma::mat injectQ(arma::mat &Q, const arma::uvec indexQ, const arma::vec q) {
    if(indexQ.n_elem > 0) {
        Q.elem(indexQ) = q;
    }
    return(Q);
}

inline arma::mat injectG(arma::mat &G, const arma::uvec indexG, const arma::vec g) {
    if(indexG.n_elem > 0) {
        G.elem(indexG) = g;
    }
    return(G);
}

// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::export]]
Rcpp::List bsts_posterior_predict(Rcpp::NumericVector model, arma::rowvec Z, arma::mat G, arma::mat Q, arma::mat R, arma::mat init_state, 
                                  arma::mat X, arma::mat B, arma::mat gpriors, arma::mat qpriors, arma::vec vpriors, arma::mat U)
{
    // model: time[t] series[n]
    try {
        int t = static_cast<int>(model[0]);
        int draws = static_cast<int>(model[1]);
        int n_states = static_cast<int>(model[2]);
        int j,i;
        arma::uvec indexQ = arma::find_nonfinite(Q);
        arma::uvec indexG = arma::find_nonfinite(G);
        arma::cube pstates(n_states, t, draws);
        int nq = Q.n_cols;
        arma::mat y(draws, t + 1);
        arma::mat a(n_states, t + 1);
        for(j = 0; j<draws; j++){
            Q = injectQ(Q, indexQ, qpriors.col(j));
            G = injectG(G, indexG, gpriors.col(j));
            a.col(0) = init_state.col(j);
            arma::vec obs_noise = U.col(j) * vpriors(j);
            arma::mat state_noise = arma::trans(rmvnorm(t+1, arma::zeros<arma::vec>(nq), Q));
            for(i = 1; i<=t; i++){
                a.col(i) = (R * state_noise.col(i-1)) + G * a.col(i-1);
                y(j,i) = arma::as_scalar(Z * a.col(i) + X.row(i) * B.col(j) + obs_noise(i));
            }
            pstates.slice(j) = a.cols(1,t);
        }
        
        Rcpp::List output = Rcpp::List::create(Rcpp::Named("y") = y.cols(1, t),
                                               Rcpp::Named("states") = pstates,
                                               Rcpp::Named("Q") = Q, Rcpp::Named("G") = G);
        return(output);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "tsforeign--> bsts_posterior_predict exception (unknown reason)" );
    }
    return R_NilValue;
}

