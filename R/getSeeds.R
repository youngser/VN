#'
#' Finds a set of seeds
#'
#' @param W vector of indices of the shared vertices of \eqn{G_1} \eqn{G_2}
#' @param x vector of indices for vertices of interest (voi) in \eqn{G_1}
#' @param s number of seeds
#'
#' @return a vector of \eqn{s} seeds for \eqn{sgm}
#'
#' @author Youngser Park <youngser@jhu.edu>
#' @export
#

getSeeds <- function(W,x,s){
#    W <- intersect(V(g1),V(g2))
    W <- setdiff(W,x) # exclude x from W
    (maxseed <- min(length(W),s))
    (S <- sort(sample(W,maxseed)))

    return(S)
} ## END FUNCTION

getSeeds2 <- function(g1,g2,x,s,h){
    W <- intersect(V(g1),V(g2))
    W <- setdiff(W,x) # exclude x from W
    N1 <- unlist(ego(g1,h,nodes=x,mindist=1)) # mindist=0: close, 1: open, don't have x
    s1 <- sample(intersect(W,N1),1)

    Ns1 <- unlist(ego(g1,h,nodes=s1,mindist=1)) # mindist=0: close, 1: open
    Ns2 <- unlist(ego(g2,h,nodes=s1,mindist=1)) # mindist=0: close, 1: open
    Ns <- setdiff(intersect(Ns1,Ns2),x) # don't have x
    maxseed <- min(length(Ns),s-1)
    S <- c(s1,sample(Ns,maxseed))

    return(S)
} ## END FUNCTION
