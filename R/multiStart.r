#'
#' Runs SGM algorithm starting with rsp (rand start point)
#' \code{R} times and generates an average probability matrix over the \code{R}
#' restarts
#'
#' @param A adjacency matrix for \eqn{G_1}
#' @param B adjacency for correlated matrix \eqn{G_2}
#' @param R number of restarts
#' @param s number of seeds
#' @param g gamma (upper bound on alpha) should be in (0,1), max tolerance for alpha, how far away from the barycenter user is willing to go for
#' the initialization of \code{sgm} on any given iteration
#' @param pad a scalar value for padding for sgm
#' @param iter number of iterations allowed to be used in sgm step
#'
#' @return \eqn{n \times n} average probability matrix over all restarts.
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export
### S.D.G
#


#sgm3 <- function(A,B,m,start,iteration){

multiStart <- function(A,B,R,s,g,pad=0,iter=20){
#     source("~Hspace/sgmvn/ForTopK/Functions/sgm3.R")
    n <- max(nrow(A),nrow(B))
    P <- matrix(0,n,n)

    for(i in 1:R){
        M <- rsp(n-s,g)
#        if (n-s==1) {
#            Ps <- sgm.igraph(A,B,s,start=M,iteration=iter)
#        } else {
            Ps <- sgm.ordered(A,B,s,start=M,pad=pad,iteration=iter)
#        }
        P <- P + Ps$P
    } ## END FOR

    D <- (1/R)*P

    return(D)
} ## END FUNCTION



#   Time:
##  Working status:
### Comments:
####Soli Deo Gloria
