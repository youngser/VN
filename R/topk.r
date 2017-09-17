#' Creates top \code{k} nomination list of vertices in \eqn{G_2} for each vertex
#'      in \eqn{G_1}
#'
#' @param k number of vertices we will consider from \eqn{G_2}
#' @param M multistart matrix (doubly stochastic)
#'
#' @return  probability matrix based on likelihood each vertex in
#'                 \eqn{G2} is matched to corresponding vertex of
#'                 \eqn{G_1} i.e. the \eqn{P[i,:]} gives the
#'                 probabilities of the top \code{k}
#'                 most likely candidates from \eqn{G_2} for vertex
#'                 \eqn{i} in \eqn{G_1}
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export
### S.D.G
#


topk <- function(k,M){
#    M <- multiStart(n,g,R)
    noms = M
    n <- nrow(M)
    for(i in 1:n){
        q <- sort(M[i,],decreasing=TRUE,index.return=TRUE)
        q$x[(k+1):n] = 0;

        p = rep(0,n)
        for(i in 1:n){
            p[q$ix[i]] = q$x[i]
        } ## END FOR

        noms[i,]=p
    } ## END FOR

    return(nominations=noms)
} ## END FUNCTION




#   Time:
##  Working status:
### Comments:
####Soli Deo Gloria
