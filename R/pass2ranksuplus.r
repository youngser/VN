#############################################################################
#' Generates a matrix whose \eqn{i,j-th} entry represents the location of
#' element \eqn{i,j} when the entries in row \eqn{i} are ordered from
#'  greatest to least. ties broken by averaging.
#'
#' @param M A matrix that we wish to pass to ranks on.
#' @return ranked matrix (by row)
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export

### Used this one in the final versions of the code.
pass2ranksuplus <- function(M){
    rankM <- matrix(NA,nrow(M),ncol(M))
    for(ro in 1:nrow(M)){
        #rankM[ro,] <- rank(-M[ro,],ties.method="min")
        rankM[ro,] <- rank(-M[ro,],ties.method="average")
    }
    return(rankM)
}

pass2ranks <- function(Pst){
    p2r <- t(apply(-t(Pst),1,rank,ties.method='first'))
    return(p2r)
}

#############################################################################
#' Generates a matrix whose \eqn{i,j-th} entry represents the location of
#' element \eqn{i,j} when the entries in row \eqn{i} are ordered from
#'  greatest to least. ties broken by averaging.
#'
#' @param M A matrix that we wish to pass to ranks on.
#' @return ranked matrix (element-wise with rank 1 as smallest number in
#' matrix) and unique entries sorted
#'
#' @author Heather Gaddy Patsolic <hgaddy1@jhu.edu>
#' @export

pass2ranksu <- function(M){
    SortM <- sort(unique(as.vector(M)))
    Mprime <- matrix(NA,nrow(M),ncol(M))
    for(u in 1:length(SortM)){
        for(v in 1:nrow(M)){
            for(c in 1:ncol(M)){
                if(M[v,c] == SortM[u]){
                    Mprime[v,c] <- u
                }
            }
        }
    }
    rankList <- list(mat=Mprime,label=SortM) #returns the ranked matrix and the unique entries sorted
    return(rankList)
}


