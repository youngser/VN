#' Generates a pair of graphs from a \eqn{rho}-correlated
#'              SBM
#'
#' @param  B block model matrix of probabilities
#' @param  n vector of sizes of blocks in main matrix A
#' @param  m vector of sizes of corresponding matrix created from A
#' @param  P matrix of probabilities size of sum(n)
#' @param  corr correlation of edge in B v. A (should be positive)
#' @param  setseed True, then set seed. False, then doesn't
#'
#' @return  A  the primary adjacency matrix A
#' @return  Bcor the corresponding matrix from A
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export

randmatcorr <- function(B,n,m,P,corr,setseed=FALSE){
    library(igraph)
    #source("~/Hspace/sgmvn/FromVince/adjcorr2.R")
    #source("~/Hspace/sgmvn/FromVince/adjcorr.R")
    #source("~/Hspace/sgmvn/ForTopK/Functions/randmat.r")
    #source("~/Hspace/sgmvn/ForTopK/Functions/adjcorrH.r")
    if(setseed){
        set.seed(setseed)
        #set.seed(15345)
    }
    A <- randmat(B,n,m)
    np<-c(0,cumsum(n))
    samps <- NULL
    for(i in 1:length(m)){ # get m[i] samples from each block i
        samps <- c(samps,sample((np[i]+1):np[i+1], m[i]))
    }
    samps <- sort(samps)
    A1 <- A[samps,samps]
    nsamps <- sort(setdiff(c(1:nrow(A)),samps))
    A <- A[c(samps,nsamps),c(samps,nsamps)] # put m "samps" first followed by the rest of n-m "nsamps"
    #A <- A[c(nsamps,samps),c(nsamps,samps)]
    P1 <- P[samps,samps]
#    Bcorr <- adjcorrH(A1,P=P1,corr,permutation)
    Bcorr <- adjcorrH(A1,P=P1,corr) # generate a correlated graph of the first m vertices of A
#    Q<-P+corr*(1-P)
#    Bcorr <- adjcorr2(A1,P=P1,Q,corr)
    L <- list(A=A,Bcorr=Bcorr)
    return(L)
} # END FUNCTION

randmat <- function(B,n,m){
    l <- nrow(B)
    #set.seed(4328437)
    #set.seed(123456)
    if (l == length(n) && length(m) == length(n)){
    G1<-sbm.game(sum(n),B,n,directed=FALSE,loops=FALSE)
    #G2<-sbm.game(sum(m),B,m,directed=FALSE,loops=FALSE)
    A<-as.matrix(get.adjacency(G1))
    return(A)
    #C<-as.matrix(get.adjacency(G2))
    #L <- list(A=A,Acor=C)
    #return(L)
    }
    else{
        print("ERROR: Dimensions must agree")
        L <- NULL
        return(L)
    }
} ## END FUNCTION


# if(FALSE){
# B <- matrix(c(.7,.3,.4,.3,.7,.3,.4,.3,.7),3,3)
# n <- c(100,100,100)
# m <- c(75,75,75)
# np <- c(0,n)
# mp <- c(0,m)
# P <- matrix(0,sum(n),sum(n))
# for(i in 1:length(n)){
#     for(j in 1:length(n)){
#         P[(sum(np[1:i])+1):sum(np[1:i+1]),
#           (sum(np[1:j])+1):sum(np[1:j+1])] <-
#             B[i,j]*matrix(1,np[i+1],np[j+1])
#     }
# }
# corr <- 0.8
# sm <- sample(sum(m))
# I <- diag(sum(m))
# permutation <- I[sm,]
# }

#   Time:
##  Working status:
### Comments:
####Soli Deo Gloria
