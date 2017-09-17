#' Generates a matrix \eqn{B} which has correlation \eqn{rho}
#'          with the input matrix
#'
#' @param A adjacency matrix of the original graph \eqn{G}
#' @param P matrix of probabilities, same size as A
#' @param corr correlation of edge in B v. A (should be positive)
#' @return a graph of the same size as \code{A} such that the correlation coefficient between the entries of the two adjacency matrices is \code{corr}.
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export

#adjcorrH <- function(A,P,corr,permutation){
adjcorrH <- function(A,P,corr){
    Q <- P+corr*(1-P)
    n <- nrow(A)
    B <- A
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            if(A[i,j]==1 && runif(1) > Q[i,j]){
                B[i,j]=0
                B[j,i]=0
            }
            else if(A[i,j]==0 && runif(1) < P[i,j]*(1-corr)){
                B[i,j]=1
                B[j,i]=1
            }
        }
    }
#    B<-permutation%*%B%*%t(permutation)
#    B<-as.matrix(B)
    return(B)
}


