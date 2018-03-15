# This code uses the same input/output data massaging structures as the Smacof R package
# De Leeuw, Jan, and Patrick Mair. "Multidimensional scaling using majorization: SMACOF in R." Department of Statistics, UCLA (2011).
# Smacof package available at https://CRAN.R-project.org/package=smacof

fJOFC <- function(w1,delta, ndim = 2,verbose = FALSE, init=NULL,
                        itmax = 1000, eps = 1e-6)  
  {
  require(shapes)
  require(smacof)
  # w1 is the value such that W is then of the form:
  #       ee^T-I on diagonal blocks
  #       w1 I off diagonal
  # delta ... list of dissmilarity matrices to be jointly embedded
  #           vertices are assumed matched
  m<-length(delta)
  nvert<-nrow(delta[[1]])
  # ndim ... number of dimensions
  # init ... matrix with starting values of dimension nvert*m \times ndim
  # itmax ... maximum number of iterations
  # eps ... change in loss function


  # make the delta's into diss objects, same subroutine used in SmacofSym from the Smacof package

  diss<-delta
  for(i in 1:m){
  if ((is.matrix(diss[[i]])) || (is.data.frame(diss[[i]]))) {
    dx<-as.vector(diss[[i]][lower.tri(diss[[i]])])
    diss[[i]] <- structure(dx, Size = nvert, call = quote(as.dist.default(m=b)), 
                           class = "dist", Diag = FALSE, Upper = FALSE) 
    attr(diss[[i]], "Labels") <- c((nvert*(i-1)+1):(nvert*i))
  }
  }
  
  p <- ndim                                     
  n <- nvert*m
  NN<-n*(n-1)/2
  if (p > (n - 1)) stop("Maximum number of dimensions is nvert*m-1!")
  
  nn <- nvert*(nvert-1)/2
  m <- length(diss)
  
  wgthsd<-as.dist(matrix(1,nvert,nvert)-diag(nvert))
  wgthso<-w1*diag(nvert)
  
  # can parallelize here
  dhat<-list()
  for(i in 1:m){
    qqq<-length(diss[[i]])
  dhat[[i]] <- diss[[i]]/sqrt(sum(wgthsd*diss[[i]]^2))*sqrt(qqq)*sqrt(NN/(m*nn))  #normalize each dissimilarity to nvert(nvert-1)/2
  }
  if (is.null(init)){   
  # to initialize, instead of embedding everything at once using torgerson
  # which would be inefficient speed-wise, let's embed each separately
  # and procrustes fit them to each other and call that the initialization
      x<-list()
      xmean <- torgerson(sqrt(Reduce('+',diss)), p=p)
      for(i in 1:m){
          x[[i]] <- torgerson(sqrt(diss[[i]]), p=p) # intialize each x
          xx     <- procOPA(xmean,x[[i]],scale=FALSE,reflect=TRUE)
          x[[i]] <- xx$Bhat
      } 
  } else{
      x <- list()
      for(i in 1:m){x[[i]]<-as.matrix(init)[((i-1)*nvert+1):(nvert*i),]}
  }
  
  itel <- 1  #iteration number
  # make the distance matrix list
  d    <- list()
  for (i in 1:m){
      d[[i]]<-list()
      for (j in i:m){
          d[[i]][[j]]<-rect.dist.sqrt(x[[i]],x[[j]])
          if (i==j){d[[i]][[j]]<-as.dist(d[[i]][[j]])}
      }
  }
  numer   <-0
  denom   <-0
  for(i in 1:m){numer <-numer+sum(d[[i]][[i]]*dhat[[i]])}
  for(i in 1:m){
    for(j in i:m){
      if(i==j){denom<-denom+sum(wgthsd*d[[i]][[i]]^2)}
      else{denom<-denom+sum(w1*diag(d[[i]][[j]]^2))}
     }
  }
  
  lb   <- numer/denom      #normalization tr(X'VX); 
  for(i in 1:m){
      x[[i]]    <- lb*x[[i]]      #modify x with lb-factor
      for(j in i:m){
          d[[i]][[j]]    <- lb*d[[i]][[j]]  #modify d with lb-factor
      }
  }
  c1   <- (nvert+w1)/(nvert*(nvert+m*w1))
  c2   <- w1/(nvert*(nvert+m*w1))
  
  # calculate stress (can be parallelized)
  sold<-0
  for(i in 1:m){
    for(j in i:m){
      if(i==j){sold<-sold+sum(wgthsd*(dhat[[i]]-d[[i]][[j]])^2)}
      else{sold<-sold+sum(w1*diag(d[[i]][[j]]^2))}
    }
  }
  sold <- sold/(NN)         #stress (to be minimized in repeat loop)
  
  repeat {    
    # make B*x's
    bx <- list()
    for(i in 1:m){
        qqq<-ifelse(d[[i]][[i]]<1e-12,1,0)
        bbb<-as.matrix( (wgthsd*dhat[[i]]*(1-qqq)) /(d[[i]][[i]]+qqq))
        rrr<-rowSums(bbb)
        b <- diag(rrr)-bbb
        bx[[i]] <- b%*%x[[i]]
    }
    # make updated embedding
    y <- list()
    ysum<-c2*Reduce("+",bx)
    for(i in 1:m){
      y[[i]]<-ysum+(c1-c2)*bx[[i]]
    }
    # make distance matrices for y
    e    <- list()
    for (i in 1:m){
      e[[i]]<-list()
      for (j in i:m){
        e[[i]][[j]]<-rect.dist.sqrt(y[[i]],y[[j]])
        if (i==j){e[[i]][[j]]<-as.dist(e[[i]][[j]])}
      }
    }

    
    # calculate new stress (can be parrallelized)
    snon<-0
    for(i in 1:m){
      for(j in i:m){
        if(i==j){snon<-snon+sum(wgthsd*(dhat[[i]]-e[[i]][[j]])^2)}
        else{snon<-snon+sum(w1*diag(e[[i]][[j]])^2)}
      }
    }
    snon <- snon/NN
    
    #print out intermediate stress, same output structure as in SmacofSym to allow for ease of comparison
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d"),
                     " Stress (normalized):", formatC(c(snon),digits=8,width=12,format="f"),
                     " Difference: ", formatC(sold-snon,digits=8,width=12,format="f"),"\n")
    
    if (((sold-snon)<eps) || (itel == itmax)) break()
    x <- y                           #update configurations
    d <- e                           #update configuration distances
    sold <- snon                     #update stress
    itel <- itel+1	                 #increase iterations
  }
  
  stress <- sqrt(snon)               #stress normalization
  
  confdiss<-e
  for(i in 1:m){
    for(j in i:m){
      if(i==j){
        qqq<-length(e[[i]][[j]])
        confdiss[[i]][[j]] <- e[[i]][[j]]/sqrt(sum(wgthsd*e[[i]][[j]]^2))*sqrt(qqq)}
      else{
        qqq<-length(e[[i]][[j]])
        confdiss[[i]][[j]] <- e[[i]][[j]]/sqrt(sum(wgthso*e[[i]][[j]]^2))*sqrt(qqq)}
    }
  }
  
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!") 
  y <- do.call(rbind,y)
  
  #return configurations, configuration distances, normalized observed distances 
  result <- list(delta = diss, obsdiss = dhat, confdiss = confdiss, conf=y,stress=stress,niter = itel)
  result 
}

rect.dist.sqrt <- function(X,Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  m <- nrow(Y)
  tmp1 <- X%*%t(Y)
  tmp2 <- outer(rep(1, n), rowSums(Y^2))
  tmp3 <- outer(rowSums(X^2), rep(1,m))
  
  D <- tmp2 - 2*tmp1 + tmp3
  D[D<0]<-0
  return(sqrt(D))
}

