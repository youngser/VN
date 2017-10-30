// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Rcpp::as;

//' Multiply two matrices using RcppArmadillo
//'
//' @param A B two matrices
//' @export
// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
    arma::mat C = A * B;

    return Rcpp::wrap(C);
}

//' Transpose of a matrix using RcppArmadillo
//'
//' @param A matrix
//' @export
// [[Rcpp::export]]
SEXP armaMatT(arma::mat A){
    arma::mat C = A.t();

    return Rcpp::wrap(C);
}

//' Multiply two matrices using RcppEigen
//'
//' @param A B two matrices
//' @export
// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

//' Transpose of a matrix using RcppEigen
//'
//' @param A matrix
//' @export
// [[Rcpp::export]]
SEXP eigenMatT(Eigen::MatrixXd A){
    Eigen::MatrixXd C = A.transpose();

    return Rcpp::wrap(C);
}

//' Multiply two matrices using RcppEigen::Map
//'
//' @param A B two matrices
//' @export
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}

//' Transpose matrix using RcppEigen::Map
//'
//' @param A matrix
//' @export
// [[Rcpp::export]]
SEXP eigenMapMatT(const Eigen::Map<Eigen::MatrixXd> A){
    Eigen::MatrixXd At(A.transpose());

    return Rcpp::wrap(At);
}
