#include <RcppArmadillo.h>
//using namespace Rcpp;

using namespace arma;


// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
arma::uvec my_setdiff1(arma::uvec& x, const arma::uvec& y){

  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}

// [[Rcpp::export]]
double Update_rho_S1(arma::mat Zc1, arma::mat Zc2, arma::vec Yc1, arma::vec Yc2, double a0H, double a1H,
                     arma::vec betaH, double sigma1H, double sigma2H, arma::uvec Index) {
  int nc2 = Zc2.n_rows;
  arma::vec mean1 = Zc1 * betaH;
  arma::vec mean3 = a1H + Zc2 * betaH;
  Index = Index - 1;  // Index different from R

  arma::vec Ys1 = Yc1.elem(Index);
  double W1p = sum( (Ys1 - mean1.elem(Index))%(Ys1 - mean1.elem(Index)));
  double W3p = sum( (Yc2 - mean3)% (Yc2 - mean3));
  double W13p = sum( (Ys1 - mean1.elem(Index))%(Yc2 - mean3) );

  double c3 = 1, c2 = -W13p/(nc2*sigma1H*sigma2H);
  double c1 = (W1p/(nc2*sigma1H*sigma1H) +  W3p/(nc2*sigma2H*sigma2H) - 1),
    c0 = -W13p/(nc2*sigma1H*sigma2H);
  arma::vec P(4);
  P(0) = c3, P(1) = c2, P(2) = c1, P(3) = c0;
  arma::cx_vec R = roots(P);
  double rhoH = real(R[2]); // real number
  if (rhoH > 0.99)
    rhoH = 0.97;
  else if (rhoH <-0.99)
    rhoH = -0.97;
  return rhoH;
}


// [[Rcpp::export]]
double Update_sigma1_S1(arma::mat Zc1, arma::mat Zc2, arma::vec Yc1, arma::vec Yc2, double a0H, double a1H,
                        arma::vec betaH, double rhoH, double sigma2H, arma::uvec Index) {
  int nc1 = Yc1.size(), nc2 = Yc2.size();
  arma::vec mean1 = Zc1 * betaH;
  arma::vec mean3 = a1H + Zc2 * betaH;
  Index = Index - 1;
  arma::uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff1(allid, Index);

  arma::vec Ys1 = Yc1.elem(Index);
  double W1p = sum( (Ys1 - mean1.elem(Index))%(Ys1 - mean1.elem(Index)));
  double W1s;
  if (nc2 == nc1)
    W1s = 0;
  else
    W1s = sum( (Yc1.elem(idx_c) - mean1.elem(idx_c)) %
      (Yc1.elem(idx_c) - mean1.elem(idx_c)));
  double W13p = sum( ( Yc1.elem(Index) - mean1.elem(Index) )%(Yc2 - mean3));

  if (rhoH > 0.99)
    rhoH = 0.97;
  else if (rhoH < -0.99)
    rhoH = -0.97;
  double c0 = -W1p/nc1 - (1 - pow(rhoH,2))*W1s/nc1, c1 = rhoH*W13p/(nc1*sigma2H),
    c2 = (1 - pow(rhoH, 2));
  arma::vec P = {c2, c1, c0};
  arma::cx_vec R=roots(P);
  arma::vec RT = real(R);
  double sigma1H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma1H;
}

// [[Rcpp::export]]
double Update_sigma2_S1(arma::mat Zc1, arma::mat Zt2, arma::mat Zc2, arma::vec Yc1, arma::vec Yt2, arma::vec Yc2,
                        double a0H, double a1H,
                        arma::vec betaH, double rhoH, double sigma1H, arma::uvec Index) {
  int nc2 = Zc2.n_rows;
  int nt2 = Zt2.n_rows;
  int N2 = nc2 + nt2;
  arma::vec mean1 = Zc1*betaH;
  arma::vec mean2 = a0H + a1H + Zt2 * betaH;
  arma::vec mean3 = a1H + Zc2 * betaH;

  Index = Index - 1;
  double W2 = sum( (Yt2 - mean2)% (Yt2 - mean2));
  double W3p = sum((Yc2 - mean3)%(Yc2 - mean3) );
  double W13p = sum( (Yc1.elem(Index) - mean1.elem(Index)) % (Yc2 - mean3) );
  if (rhoH > 0.99)
    rhoH = 0.97;
  else if (rhoH < -0.99)
    rhoH = -0.97;

  double c2 = (1 - pow(rhoH, 2) ), c1 = rhoH*W13p/(N2*sigma1H),
    c0 = -W3p/N2 -(1 - pow(rhoH, 2) )*W2/N2;

  arma::vec P = {c2, c1, c0};
  arma::cx_vec R=roots(P);
  arma::vec RT = real(R);
  double sigma2H = as_scalar(  RT.elem(find(RT > 0 )) );
  return sigma2H;
}

// [[Rcpp::export]]
arma::vec Update_beta_S1(arma::mat Zc1, arma::mat Zt2, arma::mat Zc2, arma::vec Yc1, arma::vec Yt2, arma::vec Yc2,
                   double a0H, double a1H, double rhoH, double sigma1H, double sigma2H, arma::uvec Index) {
  int nc1 = Yc1.size();
  int nc2 = Yc2.size();
  double R = rhoH*sigma2H/sigma1H;
  Index = Index - 1;
  arma::uvec total_vec = linspace<uvec>(0, nc1-1, nc1);
  arma::mat Cov1 = Zc2.t() * Zc2/( pow(sigma1H,2)*(1- pow(rhoH, 2)));
  arma::mat Cor1 = Zc2.t() * Yc1.elem(Index)/( pow(sigma1H,2)*(1- pow(rhoH, 2)));

  arma::uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff1(allid, Index);

  if (nc2 != nc1) {
    arma::mat Zc1cs = Zc1.rows(idx_c);
    Cov1 = Cov1 + Zc1cs.t() * Zc1cs/(sigma1H*sigma1H);
    Cor1 = Cor1 + Zc1cs.t() * Yc1.elem(idx_c)/(sigma1H*sigma1H);
  }

  arma::rowvec Zc2mean = mean(Zc2, 0);
  arma::mat Ztilde = Zc2.each_row() - (1-R)*Zc2mean;
  arma::mat Zt2_ct = Zt2.each_row() - mean(Zt2, 0);

  arma::mat  Cov2 = ( -2*rhoH*Zc2.t()*Ztilde/(sigma1H*sigma2H) +
    Ztilde.t() * Ztilde/(sigma2H*sigma2H) )/(1 - rhoH*rhoH) + Zt2_ct.t() * Zt2_ct/(sigma2H*sigma2H);

  arma::vec Ytilde = Yc2 - mean(Yc2) + R*mean(Yc1.elem(Index));

  arma::vec Cor2 = ( -Zc2.t()*Ytilde - Ztilde.t()*Yc1.elem(Index) )*rhoH/(sigma1H*sigma2H*(1-rhoH*rhoH)) +
    Ztilde.t()*Ytilde/(sigma2H*sigma2H*(1-rhoH*rhoH)) + Zt2_ct.t() * (Yt2 - mean(Yt2) )/(sigma2H*sigma2H);

  arma::mat S = Cov1 + Cov2;
  arma::vec t = Cor1 + Cor2;
  arma::vec betaH = solve(S, t);
  return betaH;
}

// [[Rcpp::export]]
arma::vec Update_a_S1(arma::mat Zc1, arma::mat Zt2, arma::mat Zc2, arma::vec Yc1, arma::vec Yt2, arma::vec Yc2,
                double a0H, double a1H, arma::vec betaH, double rhoH, double sigma1H, double sigma2H, arma::uvec Index) {

  double R = rhoH*sigma2H/sigma1H;
  Index = Index - 1;
  a1H = mean( Yc2 - R*Yc1.elem(Index) - (Zc2 - R*Zc1.rows(Index))*betaH );
  a0H = mean(Yt2 - Zt2*betaH) - a1H;
  arma::vec a(2);
  a(0) = a1H; a(1) = a0H;
  return a;
}

// [[Rcpp::export]]
double Rcpp_Objective_S1(arma::mat Zc1, arma::mat Zt2, arma::mat Zc2, arma::vec Yc1, arma::vec Yt2, arma::vec Yc2,
                         double a0H, double a1H, arma::vec betaH, double rhoH, double sigma1H, double sigma2H, arma::uvec Index) {
  Index = Index - 1;
  int nc1 = Yc1.size();
  int nc2 = Yc2.size();
  int nt2 = Yt2.size();
  arma::uvec allid = regspace<uvec>(0, 1, nc1 - 1);
  arma::uvec idx_c = my_setdiff1(allid, Index);

  arma::vec Yc1Scale = (Yc1 - Zc1*betaH)/sigma1H;
  arma::vec Yt2Scale = (Yt2 - a0H - a1H - Zt2*betaH)/sigma2H;
  arma::vec Yc2Scale = (Yc2 - a1H - Zc2*betaH)/sigma2H;

  double part1 = as_scalar( nc2*log(sigma1H) + nc2*log(sigma2H) + nc2*log(1 - pow(rhoH, 2))/2 +
                            1/(2*(1-pow(rhoH, 2))) * ( Yc1Scale.elem(Index).t()*Yc1Scale.elem(Index) -
                            2*rhoH*Yc1Scale.elem(Index).t() * Yc2Scale +  Yc2Scale.t() * Yc2Scale
                            ) + nt2*log(sigma2H) + Yt2Scale.t()*Yt2Scale/2 );
  double part2;
  if (nc2 == nc1)
    part2 = 0;
  else
    part2 = as_scalar( (nc1 - nc2)*log(sigma1H) + Yc1Scale.elem(idx_c).t() * Yc1Scale.elem(idx_c)/2 );
  return part1 + part2;
}


