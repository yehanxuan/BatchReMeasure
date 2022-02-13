#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::uvec my_setdiff2(arma::uvec& x, const arma::uvec& y){

  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}

// [[Rcpp::export]]
double Update_rho1_S2(arma::mat Zc1, arma::mat Zc3, arma::vec Yc1, arma::vec Yc3, double a3H, arma::vec betaH,
                      double sigma1H, double sigma3H, arma::uvec Index_C) {
  int nc3 = Zc3.n_rows;
  arma::vec mean1 = Zc1 * betaH;
  arma::vec mean3 = a3H + Zc3 * betaH;
  Index_C = Index_C - 1;

  arma::vec Ys1 = Yc1.elem(Index_C);
  double W1p = sum( (Ys1 - mean1.elem(Index_C))%(Ys1 - mean1.elem(Index_C)));
  double W3p = sum( (Yc3 - mean3) % (Yc3 - mean3 ));
  double W13p = sum( (Yc1.elem(Index_C) - mean1.elem(Index_C)) % (Yc3 - mean3) );
  double c3 = 1, c2 = -W13p/(nc3*sigma1H*sigma3H);
  double c1 = W1p/(nc3*sigma1H*sigma1H) + W3p/(nc3*sigma3H*sigma3H) - 1,
    c0 = -W13p/(nc3*sigma1H*sigma3H);
  arma::vec P(4);
  P(0) = c3, P(1) = c2, P(2) = c1, P(3) = c0;
  arma::cx_vec R = roots(P);
  arma::uword i = abs(imag(R)).index_min();// real number index
  double rho1H = real(R(i));
  if (rho1H > 0.99)
    rho1H = 0.95;
  else if (rho1H < -0.99 )
    rho1H = -0.95;
  return rho1H;
}

// [[Rcpp::export]]
double Update_rho2_S2(arma::mat Zt2, arma::mat Zt3, arma::vec Yt2, arma::vec Yt3,
                      double a0H, double a1H, double a3H,
                      arma::vec betaH, double sigma2H, double sigma3H, arma::uvec Index_T) {
  int nt3 = Zt3.n_rows;
  arma::vec mean2 = a0H + a1H + Zt2 * betaH;
  arma::vec mean4 = a0H + a3H + Zt3 * betaH;
  Index_T = Index_T - 1;
  double W2p = sum( (Yt2.elem(Index_T) - mean2.elem(Index_T)) %
                    (Yt2.elem(Index_T) - mean2.elem(Index_T)) );
  double W4p = sum( (Yt3 - mean4)% (Yt3 - mean4) );
  double W24p = sum( (Yt2.elem(Index_T) - mean2.elem(Index_T)) % (Yt3 - mean4) );
  double c3 = 1, c2 = -W24p/(nt3*sigma2H*sigma3H);
  double c1 = W2p/(nt3*sigma2H*sigma2H) + W4p/(nt3*sigma3H*sigma3H) - 1,
    c0 = -W24p/(nt3*sigma2H*sigma3H);

  arma::vec P(4);
  P(0) = c3, P(1) = c2, P(2) = c1, P(3) = c0;
  arma::cx_vec R = roots(P);
  arma::uword i = abs(imag(R)).index_min();// real number index
  double rho2H = real(R(i));
  if (rho2H > 0.99)
    rho2H = 0.95;
  else if (rho2H < -0.99 )
    rho2H = -0.95;
  return rho2H;
}

// [[Rcpp::export]]
double Update_sigma1_S2(arma::mat Zc1, arma::mat Zc3, arma::vec Yc1, arma::vec Yc3, double a3H, arma::vec betaH,
                        double rho1H, double sigma3H, arma::uvec Index_C) {
  int nc1 = Zc1.n_rows, nc3 = Zc3.n_rows;
  arma::vec mean1 = Zc1 * betaH;
  arma::vec mean3 = a3H + Zc3 * betaH;
  Index_C = Index_C - 1;
  arma::uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff2(allid, Index_C);

  double W1p = sum( (Yc1.elem(Index_C) - mean1.elem(Index_C) ) %
                    (Yc1.elem(Index_C) - mean1.elem(Index_C) ));
  double W1s;
  if (nc3 == nc1)
    W1s = 0;
  else
    W1s = sum( (Yc1.elem(idx_c) - mean1.elem(idx_c)) %
      (Yc1.elem(idx_c) - mean1.elem(idx_c)));
  double W13p = sum( (Yc1.elem(Index_C) - mean1.elem(Index_C)) % (Yc3 - mean3) );
  if (rho1H > 0.99)
    rho1H = 0.95;
  else if (rho1H < -0.99)
    rho1H = -0.95;
  double c0 = ( -W1p - (1 - rho1H*rho1H)*W1s )/nc1, c1 = rho1H*W13p/(nc1*sigma3H),
    c2 = (1 - rho1H*rho1H);
  arma::vec P = {c2, c1, c0};
  arma::cx_vec R = roots(P);
  arma::vec RT = real(R);
  double sigma1H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma1H;
}

// [[Rcpp::export]]
double Update_sigma2_S2(arma::mat Zt2, arma::mat Zt3,arma::vec Yt2, arma::vec Yt3,
                        double a0H, double a1H, double a3H,
                        arma::vec betaH, double rho2H, double sigma3H, arma::uvec Index_T) {
  int nt2 = Zt2.n_rows, nt3 = Zt3.n_rows;
  arma::vec mean2 = a0H + a1H + Zt2 * betaH;
  arma::vec mean4 = a0H + a3H + Zt3 * betaH;
  Index_T = Index_T - 1;
  double W2p = sum( (Yt2.elem(Index_T) - mean2.elem(Index_T)) %
                    (Yt2.elem(Index_T) - mean2.elem(Index_T)) );
  arma::uvec allid = regspace<uvec>(0, 1, nt2 - 1);
  arma::uvec idx_t = my_setdiff2(allid, Index_T);
  double W2s;
  if (nt3 == nt2)
    W2s = 0;
  else
    W2s = sum( (Yt2.elem(idx_t) - mean2.elem(idx_t)) %
      (Yt2.elem(idx_t) - mean2.elem(idx_t)) );
  double W24p = sum( (Yt2.elem(Index_T) - mean2.elem(Index_T)) % (Yt3 - mean4) );
  if (rho2H > 0.99)
    rho2H = 0.95;
  else if (rho2H < -0.99)
    rho2H = -0.95;
  double c0 = -W2p/nt2 - (1 - rho2H*rho2H)*W2s/nt2, c1 = rho2H*W24p/(nt2*sigma3H),
    c2 = (1 - rho2H*rho2H);
  arma::vec P = {c2, c1, c0};
  arma::cx_vec R=roots(P);
  arma::vec RT = real(R);
  double sigma2H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma2H;
}

// [[Rcpp::export]]
double Update_sigma3_S2(arma::mat Zc1, arma::mat Zt2, arma::mat Zc3, arma::mat Zt3, arma::vec Yc1, arma::vec Yt2, arma::vec Yc3, arma::vec Yt3,
                        double a0H, double a1H, double a3H, arma::vec betaH, double rho1H, double rho2H,
                        double sigma1H, double sigma2H, arma::uvec Index_C, arma::uvec Index_T) {
  int nc3 = Zc3.n_rows, nt3 = Zt3.n_rows;
  int N3 = nc3 + nt3;
  arma::vec mean1 = Zc1 * betaH;
  arma::vec mean2 = a0H + a1H + Zt2 * betaH;
  arma::vec mean3 = a3H + Zc3 * betaH;
  arma::vec mean4 = a0H + a3H + Zt3 * betaH;
  Index_C = Index_C - 1;
  Index_T = Index_T - 1;
  double W3p = sum( (Yc3 - mean3) % (Yc3 - mean3));
  double W4p = sum( (Yt3 - mean4) %  (Yt3 - mean4) );
  double W13p = sum( (Yc1.elem(Index_C) - mean1.elem(Index_C))%(Yc3 - mean3) );
  double W24p = sum( (Yt2.elem(Index_T) - mean2.elem(Index_T) )%(Yt3 - mean4) );
  if (rho1H > 0.99)
    rho1H = 0.95;
  else if (rho1H < -0.99)
    rho1H = -0.95;
  if (rho2H > 0.99)
    rho2H = 0.95;
  else if (rho2H < -0.99)
    rho2H = -0.95;
  double c0 = -(1 - rho2H*rho2H)*W3p/N3 - (1 - rho1H*rho1H)*W4p/N3;
  double c1 = (1 - rho2H*rho2H)*rho1H*W13p/sigma1H  + (1 - rho1H*rho1H)*rho2H*W24p/sigma2H;
  c1 = c1/N3;
  double c2 = (1 - rho1H*rho1H)*(1 - rho2H*rho2H);
  arma::vec P = {c2, c1, c0};
  arma::cx_vec R = roots(P);
  arma::vec RT = real(R);
  double sigma3H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma3H;
}

// [[Rcpp::export]]
Rcpp::List Update_a_S2(arma::mat Zc1, arma::mat Zt2, arma::mat Zc3, arma::mat Zt3,
                       arma::vec Yc1, arma::vec Yt2, arma::vec Yc3, arma::vec Yt3,
                       arma::vec betaH, double rho1H, double rho2H,
                      double sigma1H, double sigma2H, double sigma3H, arma::uvec Index_C, arma::uvec Index_T) {
  int nt2 = Zt2.n_rows, nt3 = Zt3.n_rows;
  Index_C = Index_C - 1;
  Index_T = Index_T - 1;
  double R1p = mean( Yc1.elem(Index_C) - Zc1.rows(Index_C) * betaH );
  double R3p = mean(Yc3 - Zc3 * betaH);
  double R2p = mean( Yt2.elem(Index_T) - Zt2.rows(Index_T) * betaH );
  double R4p = mean( Yt3 - Zt3 * betaH);
  double R2s;
  arma::uvec allid_t = regspace<uvec>(0, 1, nt2 - 1);
  arma::uvec idx_t = my_setdiff2(allid_t, Index_T);
  if (nt3 == nt2)
    R2s = 0;
  else
    R2s = mean( Yt2.elem(idx_t) - Zt2.rows(idx_t)*betaH );
  double ratioT = static_cast<double>(nt3)/static_cast<double>(nt2);

  double a3H = R3p - rho1H*sigma3H*R1p/sigma1H;
  double a0H = R4p - a3H - rho2H*sigma3H/sigma2H*
    (1 - ratioT)*(R2p-R2s);
  double a1H = a3H +
    (ratioT + (1-ratioT)*rho2H*sigma3H/sigma2H )*R2p - R4p +
    (1 - ratioT)*(1 - rho2H*sigma3H/sigma2H)*R2s;
  return Rcpp::List::create(Rcpp::Named("a0H") = a0H,
                            Rcpp::Named("a1H") = a1H,
                            Rcpp::Named("a3H") = a3H);
}

// [[Rcpp::export]]
arma::vec Update_beta_S2(arma::mat Zc1, arma::mat Zt2, arma::mat Zc3, arma::mat Zt3,
                         arma::vec Yc1, arma::vec Yt2, arma::vec Yc3, arma::vec Yt3, double rho1H, double rho2H,
                         double sigma1H, double sigma2H, double sigma3H, arma::uvec Index_C, arma::uvec Index_T) {
  if (sigma3H < 0.01)
    sigma3H = 0.1;
  int nc1 = Zc1.n_rows, nt2 = Zt2.n_rows, nc3 = Zc3.n_rows, nt3 = Zt3.n_rows;
  Index_C = Index_C - 1;
  Index_T = Index_T - 1;

  arma::uvec allid_c = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff2(allid_c, Index_C);
  arma::uvec allid_t = regspace<uvec>(0, 1, nt2 - 1);
  arma::uvec idx_t = my_setdiff2(allid_t, Index_T);
  double ratioT = static_cast<double>(nt3)/static_cast<double>(nt2);
  double R1 = rho1H*sigma3H/sigma1H, R2 = rho2H*sigma3H/sigma2H;
  arma::rowvec Zc3mean = mean(Zc3, 0);
  arma::mat Zc_scale = Zc3.each_row() - (1 - R1)*Zc3mean;
  arma::vec Yc_scale = Yc3 - mean(Yc3) + R1*mean(Yc1.elem(Index_C));
  arma::mat Cov1 = ( Zc3.t() * Zc3/(sigma1H*sigma1H) - 2*Zc3.t()*Zc_scale*rho1H/(sigma1H*sigma3H) +
    Zc_scale.t() * Zc_scale/(sigma3H*sigma3H) )/(1 - rho1H*rho1H);
  arma::vec Cor1 = ( Zc3.t() * Yc1.elem(Index_C)/(sigma1H*sigma1H) - Zc3.t()*Yc_scale*rho1H/(sigma1H*sigma3H) -
    Zc_scale.t() * Yc1.elem(Index_C)*rho1H/(sigma1H*sigma3H) + Zc_scale.t()*Yc_scale/(sigma3H*sigma3H) )/(1 - rho1H*rho1H);
  if (nc3 != nc1) {
    arma::mat Zc1cs = Zc1.rows(idx_c);
    Cov1 += Zc1cs.t() * Zc1cs/(sigma1H*sigma1H);
    Cor1 += Zc1cs.t() * Yc1.elem(idx_c)/(sigma1H*sigma1H);
  }

  arma::mat Zt3mean = mean(Zt3, 0);
  arma::mat Zt2_scale = Zt3.each_row() - ratioT*Zt3mean;
  arma::vec Yt2_scale = Yt2.elem(Index_T) - ratioT*mean( Yt2.elem(Index_T) );
  arma::mat Zt3_scale = Zt3.each_row()- Zt3mean;
  arma::vec Yt3_scale = Yt3 - mean(Yt3);
  //mat Zt2ts, tmp, Cov4;
  arma::mat Zt2ts, tmp;
  if (nt3 != nt2) {
    Zt2_scale = Zt2_scale.each_row() - (1-ratioT)*mean(Zt2.rows(idx_t), 0);
    Yt2_scale = Yt2_scale - (1 - ratioT)*mean(Yt2.elem(idx_t));
    Zt3_scale = Zt3_scale.each_row() +
      R2*(1 - ratioT)*( mean(Zt3, 0) - mean(Zt2.rows(idx_t), 0) );
    Yt3_scale = Yt3_scale + R2*(1-ratioT)*(mean(Yt2.elem(Index_T)) - mean(Yt2.elem(idx_t)) );
    Zt2ts = Zt2.rows(idx_t);
    tmp = Zt2ts.each_row() - (ratioT*mean(Zt3, 0) + (1-ratioT)*mean(Zt2ts, 0) );

    Cov1 += tmp.t()*tmp/(sigma2H*sigma2H);
    Cor1 += tmp.t() * (Yt2.elem(idx_t) - ratioT*mean(Yt2.elem(Index_T)) -
      (1 - ratioT)*mean(Yt2.elem(idx_t)) )/(sigma2H*sigma2H);
  }
  arma::mat Cov3 = ( Zt2_scale.t() * Zt2_scale/(sigma2H*sigma2H) - 2*rho2H*Zt2_scale.t()*Zt3_scale/(sigma2H*sigma3H) +
    Zt3_scale.t()*Zt3_scale/(sigma3H*sigma3H))/(1 - rho2H*rho2H);
  arma::mat Cov = Cov1 + Cov3;

  arma::vec Cor3 = ( Zt2_scale.t()*Yt2_scale/(sigma2H*sigma2H) + Zt3_scale.t()*Yt3_scale/(sigma3H*sigma3H) -
    Zt2_scale.t() * Yt3_scale*rho2H/(sigma2H*sigma3H) -
    Zt3_scale.t() * Yt2_scale*rho2H/(sigma2H*sigma3H) )/(1 - rho2H*rho2H);
  arma::vec Cor = Cor1  + Cor3;
  arma::vec betaH = solve(Cov, Cor);
  return betaH;
}

// [[Rcpp::export]]
double Rcpp_Objective_S2(arma::mat Zc1, arma::mat Zt2, arma::mat Zc3, arma::mat Zt3, arma::vec Yc1, arma::vec Yt2, arma::vec Yc3, arma::vec Yt3, arma::vec betaH,
                         double a0H, double a1H, double a3H, double rho1H, double rho2H,
                         double sigma1H, double sigma2H, double sigma3H, arma::uvec Index_C, arma::uvec Index_T) {
  int nc1 = Zc1.n_rows, nt2 = Zt2.n_rows, nc3 = Zc3.n_rows, nt3 = Zt3.n_rows;
  Index_C = Index_C - 1;
  Index_T = Index_T - 1;
  arma::uvec allid_c = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff2(allid_c, Index_C);
  arma::uvec allid_t = regspace<uvec>(0, 1, nt2 - 1);
  arma::uvec idx_t = my_setdiff2(allid_t, Index_T);
  arma::vec Yc1Scale = (Yc1 - Zc1 * betaH)/sigma1H;
  arma::vec Yt2Scale = (Yt2 - a0H - a1H - Zt2*betaH)/sigma2H;
  arma::vec Yc3Scale = (Yc3 - a3H - Zc3 * betaH)/sigma3H;
  arma::vec Yt3Scale = (Yt3 - a0H - a3H - Zt3*betaH)/sigma3H;

  double part1 = as_scalar( nc3*log(sigma1H*sigma3H) + nc3*log(1 - rho1H*rho1H)/2 +
    ( Yc1Scale.elem(Index_C).t()* Yc1Scale.elem(Index_C) - 2*rho1H*Yc1Scale.elem(Index_C).t()*Yc3Scale +
    Yc3Scale.t()*Yc3Scale )/(2*(1 - rho1H*rho1H)) );

  if (nc3 != nc1)
    part1 += as_scalar( (nc1 - nc3)*log(sigma1H) + Yc1Scale.elem(idx_c).t()*Yc1Scale.elem(idx_c)/2);

  double part3 = as_scalar( nt3*log(sigma2H*sigma3H) + nt3*log(1-rho2H*rho2H)/2 +
    ( Yt2Scale.elem(Index_T).t()*Yt2Scale.elem(Index_T) - 2*rho2H*Yt2Scale.elem(Index_T).t()*Yt3Scale +
    Yt3Scale.t()*Yt3Scale )/(2*(1 - rho2H*rho2H)) );

  if (nt3 != nt2)
    part3 += as_scalar( (nt2-nt3)*log(sigma2H) + Yt2Scale.elem(idx_t).t()*Yt2Scale.elem(idx_t)/2 );

  double obj = part1+part3;
  return obj;
}







