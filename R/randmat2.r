#' Generates a pair of graphs from a \eqn{rho}-correlated
#'              SBM. 
#' Assuming block 1 of the first graph is to have j vertices and block 1
#' of the second graph is to have l vertices, the upper lefthand portion
#' of the first block corresponding to vertices 1:min(j,l) in the two
#' graphs will have the specified correlation.
#'
#' @param  B block model matrix of probabilities
#' @param  n vector of sizes of blocks in main matrix A
#' @param  m vector of sizes of corresponding matrix created from A
#' @param  P matrix of probabilities size of sum(n)
#' @param  corr correlation of edge in B v. A 
#' @param  setseed True, then set seed. False, then doesn't
#'
#' @return  A  the primary adjacency matrix A
#' @return  Bcor the corresponding matrix from A
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export

randmatcorr2 <- function(B,n,m,corr,setseed=FALSE){
    library(igraph)
    if(setseed){
        set.seed(setseed)
    }

    l <- nrow(B)
    if(length(n)==l && length(n) == length(m)){
        maxvals <- rep(NA,nrow(B))
        for(j in 1:l){
            maxvals[j] <- max(n[j],m[j])
        }
        
        if(maxvals==n){
            A <- randmat(B,n)
            np<-c(0,cumsum(n))

            np1 <- c(0,n)
            P <- matrix(0,sum(n),sum(n))
            for(ip in 1:length(n)){
                for(j in 1:length(n)){
                    P[(sum(np1[1:ip])+1):sum(np1[1:ip+1]),
                      (sum(np1[1:j])+1):sum(np1[1:j+1])] <-
                        B[ip,j]*matrix(1,np1[ip+1],np1[j+1])
                }
            }

            samps <- NULL
            for(i in 1:length(m)){ # get m[i] samples from each block i
                samps <- c(samps,sample((np[i]+1):np[i+1], m[i]))
            }
            samps <- sort(samps)
            A1 <- A[samps,samps]
            nsamps <- sort(setdiff(c(1:nrow(A)),samps))
            A <- A[c(samps,nsamps),c(samps,nsamps)] # put m "samps" first followed by the rest of n-m "nsamps"
            P1 <- P[samps,samps]
            Bcorr <- adjcorrH(A1,P=P1,corr) # generate a correlated graph of the first m vertices of A
            L <- list(A=A,Bcorr=Bcorr)
            return(L)
        }else{
            Amax <- randmat(B,maxvals)
            mvp<-c(0,cumsum(maxvals))

            np1 <- c(0,maxvals)
            P <- matrix(0,sum(maxvals),sum(maxvals))
            for(ip in 1:l){
                for(j in 1:l){
                    P[(sum(np1[1:ip])+1):sum(np1[1:ip+1]),
                      (sum(np1[1:j])+1):sum(np1[1:j+1])] <-
                        B[ip,j]*matrix(1,np1[ip+1],np1[j+1])
                }
            }

            Bmax <- adjcorrH(Amax,P,corr)


            np <- c(0,cumsum(n))
            mp <- c(0,cumsum(m))
            Mp <- c(0,cumsum(maxvals))

            Aidx <- c() #rep(NA,sum(n))
            Bidx <- c() #rep(NA,sum(m))
            for(i in 1:l){
                lb <- Mp[i] + 1
                difn <- np[i+1] - np[i]
                ubn <- lb + difn - 1
                Aidx <- c(Aidx, lb:ubn)
                difm <- mp[i+1] - mp[i]
                ubm <- lb + difm - 1
                Bidx <- c(Bidx, lb:ubm)
            }
                
            A <- Amax[Aidx,Aidx]
            Bcorr <- Bmax[Bidx,Bidx]
            L <- list(A=A,Bcorr=Bcorr)
            return(L)

        }
    }else{
        print("ERROR: Dimensions must agree")
        L <- NULL
        return(L)
    }
} # END FUNCTION


#' Generates a random undirected, loop-free graph from an SBM
#'
#' @param  B block model matrix of probabilities
#' @param  n vector of sizes of each block
#'
#' @return  A the adjacency matrix for a graph realized from the
#' specified SBM
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export


randmat2 <- function(B,n){
    l <- nrow(B)
    if (l == length(n)){
        G1<-sbm.game(sum(n),B,n,directed=FALSE,loops=FALSE)
        A<-as.matrix(get.adjacency(G1))
        return(A)
    }
    else{
        print("ERROR: Dimensions must agree")
        L <- NULL
        return(L)
    }
} ## END FUNCTION



#   Time:
##  Working status:
### Comments:
####Soli Deo Gloria
