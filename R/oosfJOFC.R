#
# OOS-embedding with weighted raw-stress 
#
# This code uses some of the same input/output data massaging structures as the Smacof R package
# De Leeuw, Jan, and Patrick Mair. "Multidimensional scaling using majorization: SMACOF in R." Department of Statistics, UCLA (2011).
# Smacof package available at https://CRAN.R-project.org/package=smacof

oosfJOFC <- function(D, X, w,
                  verbose  = TRUE,
                  itmax    = 1000,
                  eps      = 1e-6) {
  ## Input:
  ##       D   : list of m dissimilarities each of size n+1 (n in sample, 1 oos)
  ##       X   : list of m within-sample embeddings, each pf dimension n x d
  ##       w   : weighting in fJOFC for oos commensurability
  
  ## N= n*m where n is the number of (in-sample) objects
  n <- nrow(X[[1]])
  m <-length(X)
  N <-n*m
  d <- ncol(X[[1]])
  
  ## NN is the total number of objects * m
  NN<- N+m
  nn <- n+1
  
  ## Normalize D,
  # make the delta's into diss objects, same subroutine used in SmacofSym from the Smacof package
  
  diss<-D
  for(i in 1:m){
    if ((is.matrix(diss[[i]])) || (is.data.frame(diss[[i]]))) {
      dx<-as.vector(diss[[i]][lower.tri(diss[[i]])])
      diss[[i]] <- structure(dx, Size = n+1, call = quote(as.dist.default(m=b)), 
                             class = "dist", Diag = FALSE, Upper = FALSE)
      attr(diss[[i]], "Labels") <- c((nn*(i-1)+1):(nn*i))
    }
  }
  wgthsd<-as.dist(matrix(1,nn,nn)-diag(nn))
  dhat<-list()
  for(i in 1:m){
    dhat[[i]] <- nDiss(diss[[i]], wgthsd)*sqrt(choose(NN,2)/(m*choose(nn,2)))  #normalize each dissimilarity to nvert(nvert-1)/2
  }
  del<-list()
  for(i in 1:m){dd<-as.matrix(dhat[[i]])
  del[[i]]<-dd[1:n,n+1]}
  
  
  ## randomly initialize y the oos embedding
  y<-list()
  for(i in 1:m){
    y[[i]] <- t(matrix(mvrnorm(1, c(rep(0,d)), diag(apply(X[[i]],2,var),d,d))))
  } 
  
  ## compute baseline stress
  oosdist <- list()
  for(i in 1:m){
  oosdist[[i]] <- rect.dist.sqrt(X[[i]],y[[i]])
  }
  stress <-sum( Reduce("+",sapply( mapply("-",oosdist,del,SIMPLIFY = FALSE) , function(x) x^2,simplify=FALSE ) ) )
  for(i in 2:m){for(j in 1:(i-1)){ stress<-stress+w*vec_sq(y[[i]]-y[[j]])}}
  stress<-stress/(n*m+choose(m,2))
  
  itel<-0 #iteration number
  if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                   " Stress (normalized):", formatC(c(stress),digits=8,width=12,format="f"),"\n")
  itel<- 1 
  
  # iterative majorization
  repeat{
    sold<-stress
    xi<-list()
    for(i in 1:m){
      xi[[i]]<- t( -del[[i]]/oosdist[[i]]+c(rep(1,n)) )%*%X[[i]]
    }
    phi<-list()
    for(i in 1:m){
      phi[[i]]<- sum( del[[i]]/oosdist[[i]])
    }
    z<-list()
   for(i in 1:m){z[[i]]<-1/(n+m*w)*xi[[i]]+w/(n*(n+m*w))*Reduce("+",xi)
   z[[i]]<-z[[i]]+(phi[[i]])/(n+m*w)*y[[i]]+w/(n*(n+m*w))*Reduce("+",mapply("*",phi,y,SIMPLIFY = FALSE))
   }
    # recompute stress
    oosdist <- list()
    for(i in 1:m){
      oosdist[[i]] <- rect.dist.sqrt(X[[i]],z[[i]])
    }
    stress <-sum( Reduce("+",sapply( mapply("-",oosdist,del,SIMPLIFY = FALSE) , function(x) x^2,simplify=FALSE ) ) )
    for(i in 2:m){for(j in 1:(i-1)){ stress<-stress+w*vec_sq(z[[i]]-z[[j]])}}
    stress<-stress/(n*m+choose(m,2))
    
    #print out intermediate stress
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                     " Stress (normalized):", formatC(c(stress),digits=8,width=12,format="f"),
                     " Difference: ", formatC(sold-stress,digits=8,width=12,format="f"),"\n")
  
    
    if (((sold-stress)<eps) || (itel == itmax)) break()
    itel<-itel+1
    y<-z
  }
  result<-list(oosDist=oosdist,Y=y,stress=stress)
  result}


vec_sq <- function(x){sum(x^2)}


procrustes <- function(X,Y, type = "1"){
  if(type == "c"){
    Y <- Y*norm(X, type = "F")/norm(Y, type = "F")
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }
  
  tmp <- t(X)%*%Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  return(list(error = norm(X%*%W - Y, type = "F"), W = W))
}

nDiss <-function(d,w)
{
  l <- length(d)
  nn <- d/sqrt(sum(w*d^2))*sqrt(l)
  return(nn)
}
