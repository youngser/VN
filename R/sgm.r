#' @useDynLib VN
#' @importFrom Rcpp sourceCpp
NULL

#' Matches Graphs given a seeding of vertex correspondences
#'
#' Given two adjacency matrices \code{A} and \code{B} of the same size, match
#' the two graphs with the help of \code{m} seed vertex pairs which correspond
#' to \code{m} rows (and columns) of the adjacency matrices.
#'
#' The approximate graph matching problem is to find a bijection between the
#' vertices of two graphs , such that the number of edge disagreements between
#' the corresponding vertex pairs is minimized. For seeded graph matching, part
#' of the bijection that consist of known correspondences (the seeds) is known
#' and the problem task is to complete the bijection by estimating the
#' permutation matrix that permutes the rows and columns of the adjacency
#' matrix of the second graph.
#'
#' It is assumed that for the two supplied adjacency matrices \code{A} and
#' \code{B}, both of size \eqn{n\times n}{n*n}, the first \eqn{m} rows(and
#' columns) of \code{A} and \code{B} correspond to the same vertices in both
#' graphs. That is, the \eqn{n \times n}{n*n} permutation matrix that defines
#' the bijection is \eqn{I_{m} \bigoplus P} for a \eqn{(n-m)\times
#' (n-m)}{(n-m)*(n-m)} permutation matrix \eqn{P} and \eqn{m} times \eqn{m}
#' identity matrix \eqn{I_{m}}. The function \code{match_vertices} estimates
#' the permutation matrix \eqn{P} via an optimization algorithm based on the
#' Frank-Wolfe algorithm.
#'
#' See references for further details.
#'
# @aliases match_vertices seeded.graph.match
#' @param A a numeric matrix, the adjacency matrix of the first graph
#' @param B a numeric matrix, the adjacency matrix of the second graph
#' @param seeds a numeric matrix, the number of seeds x 2 matching vertex table.
#' If \code{S} is \code{NULL}, then it is using a \eqn{soft} seeding algorithm.
#' @param hard a bloolean, TRUE for hard seeding, FALSE for soft seeding.
#' @param pad a scalar value for padding
#' @param maxiter The number of maxiters for the Frank-Wolfe algorithm
#' @return A numeric matrix which is the permutation matrix that determines the
#' bijection between the graphs of \code{A} and \code{B}
#' @author Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}
#' @references Vogelstein, J. T., Conroy, J. M., Podrazik, L. J., Kratzer, S.
#' G., Harley, E. T., Fishkind, D. E.,Vogelstein, R. J., Priebe, C. E. (2011).
#' Fast Approximate Quadratic Programming for Large (Brain) Graph Matching.
#' Online: \url{http://arxiv.org/abs/1112.5507}
#'
#' Fishkind, D. E., Adali, S., Priebe, C. E. (2012). Seeded Graph Matching
#' Online: \url{http://arxiv.org/abs/1209.0367}
#'
#' @export

sgm <- function (A,B,seeds,hard=TRUE,pad=0,start="barycenter",maxiter=20){
    gamma <- 0.1
    nv1<-nrow(A)
    nv2<-nrow(B)
    nv<-max(nv1,nv2)

    if(is.null(seeds)) { # no seed!
        m=0
        if (start=="barycenter") {
            S<-matrix(1/nv,nv,nv)
        } else {
            S <- rsp(nv,gamma)
        }
        AA <- A
        BB <- B
    }else {
        A.ind <- c(seeds[,1], setdiff(1:nv1, seeds[,1]))
        B.ind <- c(seeds[,2], setdiff(1:nv2, seeds[,2]))
        AA <- A[A.ind, A.ind]
        BB <- B[B.ind, B.ind]

        m <- nrow(seeds)
        if (hard==TRUE) {
            n <- nv-m
            if (start=="barycenter") {
                S <- matrix(1/n,n,n)
            } else {
                S <- rsp(n,gamma)
            }
        } else {
            s <- m
            m <- 0
            if (start=="barycenter") {
                diag1 <- diag(s)
                diag2 <- matrix(1/(nv-s),nv-s,nv-s)
                offdiag <- matrix(0,s,nv-s)
                S <- rbind(cbind(diag1,offdiag), cbind(t(offdiag),diag2))
            } else {
                M <- rsp(nv-s,gamma)
                S <- diag(nv);
                S[(s+1):nv,(s+1):nv] <- M
            }
        }
    }

    P <- sgm.ordered(AA,BB,m,S,pad,maxiter)
    return(P)
}

sgm2 <- function (A,B,seeds,hard=TRUE,pad=0,start="barycenter",maxiter=20){
    gamma <- 0.1
    if(is.null(seeds)){ # no seed!
        m=0
        nv1<-nrow(A)
        nv2<-nrow(B)
        nv<-max(nv1,nv2)
        if (start=="barycenter") {
            S<-matrix(1/nv,nv,nv)
        } else {
            S <- rsp(nv,gamma)
        }
        AA <- A
        BB <- B
    }else{
        nv1<-nrow(A)
        nv2<-nrow(B)
        nv<-max(nv1,nv2)

        s1<-seeds[,1]
        temp<-s1
        s1<-rep(FALSE,nv1)
        s1[temp]<- TRUE
        ns1 <- !s1

        s2<-seeds[,2]
        temp<-s2
        s2<-rep(FALSE,nv2)
        s2[temp]<- TRUE
        ns2 <- !s2

        A11<-A[s1,s1]
        A12<-A[s1,ns1]
        A22<-A[ns1,ns1]
        A21<-A[ns1,s1]

        AA<-rbind(cbind(A11,A12),cbind(A21,A22))

        B11<-B[s2,s2]
        B12<-B[s2,ns2]
        B22<-B[ns2,ns2]
        B21<-B[ns2,s2]

        BB<-rbind(cbind(B11,B12),cbind(B21,B22))

        if(hard==TRUE){
            m <- length(temp)
            n <- nv-m
            if (start=="barycenter") {
                S <- matrix(1/n,n,n)
            } else {
                S <- rsp(n,gamma)
            }
        }else{
            m <- 0
            s <- length(temp)
            if (start=="barycenter") {
                S <- rbind(cbind(diag(s),matrix(0,s,nv-s)),cbind(matrix(0,nv-s,s),matrix(1/(nv-s),nv-s,nv-s)))
            } else {
                M <- rsp(nv-s,gamma)
                S <- diag(nv);
                S[(s+1):nv,(s+1):nv] <- M
            }
        }
    }

    P <- sgm.ordered(AA,BB,m,S,pad,maxiter)
    return(P)
}


#' Matches Graphs given a seeding of vertex correspondences
#'
#' Given two adjacency matrices \code{A} and \code{B} of the same size, match
#' the two graphs with the help of \code{m} seed vertex pairs which correspond
#' to the first \code{m} rows (and columns) of the adjacency matrices.
#'
#' The approximate graph matching problem is to find a bijection between the
#' vertices of two graphs , such that the number of edge disagreements between
#' the corresponding vertex pairs is minimized. For seeded graph matching, part
#' of the bijection that consist of known correspondences (the seeds) is known
#' and the problem task is to complete the bijection by estimating the
#' permutation matrix that permutes the rows and columns of the adjacency
#' matrix of the second graph.
#'
#' It is assumed that for the two supplied adjacency matrices \code{A} and
#' \code{B}, both of size \eqn{n\times n}{n*n}, the first \eqn{m} rows(and
#' columns) of \code{A} and \code{B} correspond to the same vertices in both
#' graphs. That is, the \eqn{n \times n}{n*n} permutation matrix that defines
#' the bijection is \eqn{I_{m} \bigoplus P} for a \eqn{(n-m)\times
#' (n-m)}{(n-m)*(n-m)} permutation matrix \eqn{P} and \eqn{m} times \eqn{m}
#' identity matrix \eqn{I_{m}}. The function \code{match_vertices} estimates
#' the permutation matrix \eqn{P} via an optimization algorithm based on the
#' Frank-Wolfe algorithm.
#'
#' See references for further details.
#'
# @aliases match_vertices seeded.graph.match
#' @param A a numeric matrix, the adjacency matrix of the first graph
#' @param B a numeric matrix, the adjacency matrix of the second graph
#' @param m The number of seeds. The first \code{m} vertices of both graphs are
#' matched.
#' @param start a numeric matrix, the permutation matrix estimate is
#' initialized with \code{start}
#' @param pad a scalar value for padding
#' @param maxiter The number of maxiters for the Frank-Wolfe algorithm
#' @param LAP a character either "exact" or "approx"
#' @return A numeric matrix which is the permutation matrix that determines the
#' bijection between the graphs of \code{A} and \code{B}
#' @author Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}
#' @references Vogelstein, J. T., Conroy, J. M., Podrazik, L. J., Kratzer, S.
#' G., Harley, E. T., Fishkind, D. E.,Vogelstein, R. J., Priebe, C. E. (2011).
#' Fast Approximate Quadratic Programming for Large (Brain) Graph Matching.
#' Online: \url{http://arxiv.org/abs/1112.5507}
#'
#' Fishkind, D. E., Adali, S., Priebe, C. E. (2012). Seeded Graph Matching
#' Online: \url{http://arxiv.org/abs/1209.0367}
#'
#' @export
sgm.ordered <- function(A,B,m,start,pad=0,maxiter=20,LAP="exact",verbose=FALSE){
    #seeds are assumed to be vertices 1:m in both graphs
#    suppressMessages(library(clue))
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
#        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
        B <- cbind(B, matrix(pad, nrow(B), diff))
        B <- rbind(B, matrix(pad, diff, ncol(B)))
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
#        for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
        A <- cbind(A, matrix(pad, nrow(A), diff))
        A <- rbind(A, matrix(pad, diff, ncol(A)))
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m==0){
        A12 <- A21 <- B12 <- B21 <- matrix(0,n,n)
    } else {
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (n==1) {
        A12 <- A21 <- B12 <- B21 <- t(A12)
    }

    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    tA22 <- t(A22)
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tB22 <- t(B22)
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    x<- A21 %*% t(B21)
    y<- t(A12) %*% B12
    xy <- x + y
    while (toggle==1 & iter<maxiter)
    {
        iter<-iter+1
        z <- A22 %*% P %*% tB22
        w <- tA22 %*% P %*% B22
        Grad <- xy+z+w;

        if (LAP=="exact") {
            mm <- max(abs(Grad))
            ind<-matrix(clue::solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        } else { # approx
            temp <- matrix(0, n, n)
            Grad1 <- rbind(cbind(temp, t(Grad)),cbind(Grad, temp))
            ind <- parallelMatch(Grad1)
        }

        Pdir <- diag(n)
        Pdir <- Pdir[ind,]
        tPdir <- t(Pdir)
        tP <- t(P)
        wt <- tA22 %*% Pdir %*% B22
        c <- sum(diag(w %*% tP))
        d <- sum(diag(wt %*% tP)) + sum(diag(w %*% tPdir))
        e <- sum(diag(wt %*% tPdir))
        u <- sum(diag(tP %*% x + tP %*% y))
        v <- sum(diag(tPdir %*% x + tPdir %*% y))
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha <- 0
        }else{
            alpha <- -(d-2*e+u-v)/(2*(c-d+e))}
        f0 <- 0
        f1 <- c-e+u-v
        falpha <- (c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P <- alpha*P+(1-alpha)*Pdir
        }else if(f0 > f1){
            P <- Pdir
        }else{
            toggle<-0}
        if (verbose) cat("iter = ", iter, "\n")
    }
    D<-P

    if (LAP=="exact") {
        corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
    } else {
        PP <- rbind(cbind(temp, t(P)),cbind(P, temp))
        corr<- t(parallelMatch(PP))
    }
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D, iter=iter))
}

#' @export
sgm.ordered.rk1<-function(A,B,m,start,pad=0,maxiter=20,symmetric=TRUE){
    #seeds are assumed to be vertices 1:m in both graphs
#    require('clue')
    totv1<-ncol(A)
    totv2<-ncol(B)
    X<-stfp(A,1)$X
    Y<-stfp(B,1)$X
    if(symmetric==FALSE){
        XX<-stfp(t(A),1)$X
        X<-rbind(X,XX)
        YY<-stfp(t(B),1)$X
        Y<-rbind(Y,YY)
        A=A-t(X)%*%X
        B=B-t(Y)%*%Y
    }else{
        A=A-X%*%t(X)
        B=B-Y%*%t(Y)}
    if(totv1>totv2){
        diff<-totv1-totv2
#        for (j in 1:diff){B<-cbind(rbind(B,0),0)}
        B <- cbind(B, matrix(pad, nrow(B), diff))
        B <- rbind(B, matrix(pad, diff, ncol(B)))
    }else if(totv1<totv2){
        diff<-totv2-totv1
#        for (j in 1:diff){A<-cbind(rbind(A,0),0)}
        A <- cbind(A, matrix(pad, nrow(A), diff))
        A <- rbind(A, matrix(pad, diff, ncol(A)))
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m==0){
        A12<-matrix(0,n,n)
        A21<-matrix(0,n,n)
        B12<-matrix(0,n,n)
        B21<-matrix(0,n,n)
    } else {
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (n==1){
        A12=t(A12)
        A21=t(A21)
        B12=t(B12)
        B21=t(B21)
    }
    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    tA22 <- t(A22)
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tB22 <- t(B22)
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    x<- A21%*%t(B21)
    y<- t(A12)%*%B12
    while (toggle==1 & iter<maxiter)
    {
        iter<-iter+1
        z<- A22 %*% P %*% tB22
        w<- tA22 %*% P %*% B22
        Grad<-x+y+z+w;
        mm=max(abs(Grad))
        ind<-matrix(clue::solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        T<-diag(n)
        T<-T[ind,]
        tT <- t(T)
        wt<-tA22 %*% T %*% B22
        tP <- t(P)
        c<-sum(diag(w%*%tP))
        d<-sum(diag(wt%*%tP)) + sum(diag(w%*%tT))
        e<-sum(diag(wt%*%tT))
        u<-sum(diag(tP%*%x + tP%*%y))
        v<-sum(diag(tT%*%x + tT%*%y))
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha<-0
        }else{
            alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
        f0<-0
        f1<- c-e+u-v
        falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P<- alpha*P+(1-alpha)*T
        }else if(f0 > f1){
            P<-T
        }else{
            toggle<-0}
    }
    D<-P
    corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D, iter=iter))
}


#' @export
sgm.ordered.cross <- function(A,B,m,start,pad=0,maxiter=20){
    #seeds are assumed to be vertices 1:m in both graphs
    #    suppressMessages(library(clue))
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
#        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
        B <- cbind(B, matrix(pad, nrow(B), diff))
        B <- rbind(B, matrix(pad, diff, ncol(B)))
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
#        for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
        A <- cbind(A, matrix(pad, nrow(A), diff))
        A <- rbind(A, matrix(pad, diff, ncol(A)))
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m==0){
        A12<-matrix(0,n,n)
        A21<-matrix(0,n,n)
        B12<-matrix(0,n,n)
        B21<-matrix(0,n,n)
    } else {
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (n==1) {
        A12=t(A12)
        A21=t(A21)
        B12=t(B12)
        B21=t(B21)
    }

    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    x<- tcrossprod(A21, B21) # A21 %*% t(B21)
    y<- crossprod(A12, B12) #t(A12) %*% B12
    while (toggle==1 & iter<maxiter)
    {
        iter <- iter+1
        z <- tcrossprod(crossprod(t(A22),P), B22) #     A22 %*% P %*% t(B22)
        w <- crossprod(t(crossprod(A22, P)), B22) #     t(A22) %*% P %*% B22
        Grad<-x+y+z+w;
        mm=max(abs(Grad))
        ind<-matrix(clue::solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        T<-diag(n)
        T<-T[ind,]
        wt <- crossprod(t(crossprod(A22, T)), B22) #     t(A22) %*% T %*% B22
        c <- sum(diag(tcrossprod(w, P)))  # w %*% t(P)
        d <- sum(diag(tcrossprod(wt, P))) + sum(diag(tcrossprod(w, T))) # wt %*% t(P), w %*% t(T)
        e<-sum(diag(tcrossprod(wt, T))) # wt %*% t(T)
        u<-sum(diag(crossprod(P, x) + crossprod(P, y))) # t(P) %*% x, t(P) %*% y
        v<-sum(diag(crossprod(T, x) + crossprod(T, y))) # t(T) %*% x, t(T) %*% y
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha<-0
        }else{
            alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
        f0<-0
        f1<- c-e+u-v
        falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P<- alpha*P+(1-alpha)*T
        }else if(f0 > f1){
            P<-T
        }else{
            toggle<-0}
    }
    D<-P
    corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D, iter=iter))
}

#' @export
sgm.ordered.rcpp <- function(A,B,m,start,pad=0,maxiter=20){
    #seeds are assumed to be vertices 1:m in both graphs
#    suppressMessages(library(clue))
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
        B <- cbind(B, matrix(pad, nrow(B), diff))
        B <- rbind(B, matrix(pad, diff, ncol(B)))
#        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
        A <- cbind(A, matrix(pad, nrow(A), diff))
        A <- rbind(A, matrix(pad, diff, ncol(A)))
#        for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m==0){
        A12<-matrix(0,n,n)
        A21<-matrix(0,n,n)
        B12<-matrix(0,n,n)
        B21<-matrix(0,n,n)
    } else {
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (n==1) {
        A12=t(A12)
        A21=t(A21)
        B12=t(B12)
        B21=t(B21)
    }

    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    x<- armaMatMult(A21, t(B21)) #A21 %*% t(B21)
    y<- armaMatMult(t(A12), B12) #t(A12) %*% B12
    while (toggle==1 & iter<maxiter)
    {
        iter<-iter+1
        z<- armaMatMult(armaMatMult(A22, P), t(B22))  #A22 %*% P %*% t(B22)
        w<- armaMatMult(armaMatMult(t(A22), P), B22)  # t(A22) %*% P %*% B22
        Grad<-x+y+z+w;
        mm=max(abs(Grad))
        ind<-matrix(clue::solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        T<-diag(n)
        T<-T[ind,]
        wt<- armaMatMult(armaMatMult(t(A22), T), B22)    #t(A22) %*% T %*% B22
        c<-sum(diag(armaMatMult(w, t(P))))
        d<-sum(diag(armaMatMult(wt, t(P)))) + sum(diag(armaMatMult(w, t(T))))
        e<-sum(diag(armaMatMult(wt, t(T))))
        u<-sum(diag(armaMatMult(t(P), x) + armaMatMult(t(P),y)))
        v<-sum(diag(armaMatMult(t(T), x) + armaMatMult(t(T),y)))
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha<-0
        }else{
            alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
        f0<-0
        f1<- c-e+u-v
        falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P<- alpha*P+(1-alpha)*T
        }else if(f0 > f1){
            P<-T
        }else{
            toggle<-0}
    }
    D<-P
    corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D, iter=iter))
}

#' @export
sgm.ordered.gpu <- function(A,B,m,start,pad=0,maxiter=20){
    #seeds are assumed to be vertices 1:m in both graphs
#    suppressMessages(library(clue))
    suppressMessages(library(gpuR))
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
        B <- cbind(B, matrix(pad, nrow(B), diff))
        B <- rbind(B, matrix(pad, diff, ncol(B)))
        #        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
        A <- cbind(A, matrix(pad, nrow(A), diff))
        A <- rbind(A, matrix(pad, diff, ncol(A)))
        #        for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m==0){
        A12<-matrix(0,n,n)
        A21<-matrix(0,n,n)
        B12<-matrix(0,n,n)
        B21<-matrix(0,n,n)
    } else {
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (n==1) {
        A12=t(A12)
        A21=t(A21)
        B12=t(B12)
        B21=t(B21)
    }

    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    A21.g <- gpuMatrix(A21, type="float")
    tB21.g <- gpuMatrix(t(B21), type="float")
    x.g <- A21.g %*% tB21.g; x <- as.matrix(x.g)
    tA12.g <- gpuMatrix(t(A12), type="float")
    B12.g <- gpuMatrix(B12, type="float")
    y.g <- tA12.g %*% B12.g; y <- as.matrix(y.g)

    A22.g <- gpuMatrix(A22, type="float")
    B22.g <- gpuMatrix(B22, type="float")
    tA22.g <- gpuMatrix(t(A22), type="float")
    tB22.g <- gpuMatrix(t(B22), type="float")

    while (toggle==1 & iter<maxiter)
    {
        iter<-iter+1
        P.g <- gpuMatrix(P, type="float")
        tP.g <- gpuMatrix(t(P), type="float")
        z<- A22.g %*% P.g %*% tB22.g; z <- as.matrix(z)
        w.g <- tA22.g %*% P.g %*% B22.g; w <- as.matrix(w.g)
        Grad<-x+y+z+w;
        mm=max(abs(Grad))
        ind<-matrix(clue::solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        T<-diag(n)
        T<-T[ind,]
        T.g <- gpuMatrix(T, type="float")
        tT.g <- gpuMatrix(t(T), type="float")
        wt.g <- tA22.g %*% T.g %*% B22.g; wt <- as.matrix(wt.g)
        c<-sum(diag(as.matrix(w.g %*% tP.g)))
        d<-sum(diag(as.matrix(wt.g %*% tP.g))) + sum(diag(as.matrix(w.g %*% tT.g)))
        e<-sum(diag(as.matrix(wt.g %*% tT.g)))
        u<-sum(diag(as.matrix(tP.g %*% x.g) + as.matrix(tP.g %*% y.g)))
        v<-sum(diag(as.matrix(tT.g %*% x.g) + as.matrix(tT.g %*% y.g)))
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha<-0
        }else{
            alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
        f0<-0
        f1<- c-e+u-v
        falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P<- alpha*P+(1-alpha)*T
        }else if(f0 > f1){
            P<-T
        }else{
            toggle<-0}
    }
    D<-P
    corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D, iter=iter))
}


sgm.ordered.sparse <- function(A,B,m,start,pad=0,maxiter=20){
  #seeds are assumed to be vertices 1:m in both graphs
  totv1<-ncol(A)
  totv2<-ncol(B)
  if(totv1>totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv1-totv2
#    for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    B <- cbind(B, matrix(pad, nrow(B), diff))
    B <- rbind(B, matrix(pad, diff, ncol(B)))
  }else if(totv1<totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv2-totv1
#    for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
    A <- cbind(A, matrix(pad, nrow(A), diff))
    A <- rbind(A, matrix(pad, diff, ncol(A)))
  }
  totv<-max(totv1,totv2)
  n<-totv-m
  if (m==0){
      A12<-matrix(0,n,n)
      A21<-matrix(0,n,n)
      B12<-matrix(0,n,n)
      B21<-matrix(0,n,n)
  } else {
      A12<-rbind(A[1:m,(m+1):(m+n)])
      A21<-cbind(A[(m+1):(m+n),1:m])
      B12<-rbind(B[1:m,(m+1):(m+n)])
      B21<-cbind(B[(m+1):(m+n),1:m])
  }
  if (n==1) {
      A12=Matrix::t(A12)
      A21=Matrix::t(A21)
      B12=Matrix::t(B12)
      B21=Matrix::t(B21)
  }

  A22<-A[(m+1):(m+n),(m+1):(m+n)]
  B22<-B[(m+1):(m+n),(m+1):(m+n)]
  tol<-1
  P <- start
  toggle<-1
  iter<-0
  x<- A21%*%Matrix::t(B21)
  y<- Matrix::t(A12)%*%B12
  while (toggle==1 & iter<maxiter)
  {
    iter<-iter+1
    z<- A22 %*% P %*% Matrix::t(B22)
    w<- Matrix::t(A22) %*% P %*% B22
    Grad<-x+y+z+w;
    mm=max(abs(Grad))
    ind<-matrix(clue::solve_LSAP(as.matrix(Grad)+matrix(mm,totv-m,totv-m), maximum =TRUE))
    T<-diag(n)
    T<-T[ind,]
    wt<-Matrix::t(A22) %*% T %*% B22
    c<-sum(diag(w%*%Matrix::t(P)))
    d<-sum(diag(wt%*%Matrix::t(P)))+sum(diag(w%*%Matrix::t(T)))
    e<-sum(diag(wt%*%Matrix::t(T)))
    u<-sum(diag(Matrix::t(P)%*%x+Matrix::t(P)%*%y))
    v<-sum(diag(Matrix::t(T)%*%x+Matrix::t(T)%*%y))
    if( c-d+e==0 && d-2*e+u-v==0){
      alpha<-0
    }else{
      alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
    f0<-0
    f1<- c-e+u-v
    falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
    if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
      P<- alpha*P+(1-alpha)*T
    }else if(f0 > f1){
      P<-T
    }else{
      toggle<-0}
  }
  D<-P
  corr<-matrix(clue::solve_LSAP(P, maximum = TRUE))
  P=diag(n)
  P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
  corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
  return(list(corr=corr[,2], P=P, D=D, iter=iter))
}


#' Matches Graphs given a seeding of vertex correspondences (igraph implmentation)
#'
#' Given two adjacency matrices \code{A} and \code{B} of the same size, match
#' the two graphs with the help of \code{m} seed vertex pairs which correspond
#' to the first \code{m} rows (and columns) of the adjacency matrices.
#'
#' The approximate graph matching problem is to find a bijection between the
#' vertices of two graphs , such that the number of edge disagreements between
#' the corresponding vertex pairs is minimized. For seeded graph matching, part
#' of the bijection that consist of known correspondences (the seeds) is known
#' and the problem task is to complete the bijection by estimating the
#' permutation matrix that permutes the rows and columns of the adjacency
#' matrix of the second graph.
#'
#' It is assumed that for the two supplied adjacency matrices \code{A} and
#' \code{B}, both of size \eqn{n\times n}{n*n}, the first \eqn{m} rows(and
#' columns) of \code{A} and \code{B} correspond to the same vertices in both
#' graphs. That is, the \eqn{n \times n}{n*n} permutation matrix that defines
#' the bijection is \eqn{I_{m} \bigoplus P} for a \eqn{(n-m)\times
#' (n-m)}{(n-m)*(n-m)} permutation matrix \eqn{P} and \eqn{m} times \eqn{m}
#' identity matrix \eqn{I_{m}}. The function \code{match_vertices} estimates
#' the permutation matrix \eqn{P} via an optimization algorithm based on the
#' Frank-Wolfe algorithm.
#'
#' See references for further details.
#'
# @aliases match_vertices seeded.graph.match
#' @param A a numeric matrix, the adjacency matrix of the first graph
#' @param B a numeric matrix, the adjacency matrix of the second graph
#' @param m The number of seeds. The first \code{m} vertices of both graphs are
#' matched.
#' @param start a numeric matrix, the permutation matrix estimate is
#' initialized with \code{start}
#' @param maxiter The number of maxiters for the Frank-Wolfe algorithm
#' @return A numeric matrix which is the permutation matrix that determines the
#' bijection between the graphs of \code{A} and \code{B}
#' @author Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}
#' @references Vogelstein, J. T., Conroy, J. M., Podrazik, L. J., Kratzer, S.
#' G., Harley, E. T., Fishkind, D. E.,Vogelstein, R. J., Priebe, C. E. (2011).
#' Fast Approximate Quadratic Programming for Large (Brain) Graph Matching.
#' Online: \url{http://arxiv.org/abs/1112.5507}
#'
#' Fishkind, D. E., Adali, S., Priebe, C. E. (2012). Seeded Graph Matching
#' Online: \url{http://arxiv.org/abs/1209.0367}
#'
#' @export

sgm.igraph<-function(A,B,m,start,iteration){
    require(igraph)
    return(match_vertices(A,B,m,start,iteration))
}



#' Matches Graphs given a seeding of vertex correspondences
#'
#' Given two adjacency matrices \code{A} and \code{B} of the same size, match
#' the two graphs with the help of \code{m} seed vertex pairs which correspond
#' to the first \code{m} rows (and columns) of the adjacency matrices.
#'
#' The approximate graph matching problem is to find a bijection between the
#' vertices of two graphs , such that the number of edge disagreements between
#' the corresponding vertex pairs is minimized. For seeded graph matching, part
#' of the bijection that consist of known correspondences (the seeds) is known
#' and the problem task is to complete the bijection by estimating the
#' permutation matrix that permutes the rows and columns of the adjacency
#' matrix of the second graph.
#'
#' It is assumed that for the two supplied adjacency matrices \code{A} and
#' \code{B}, both of size \eqn{n\times n}{n*n}, the first \eqn{m} rows(and
#' columns) of \code{A} and \code{B} correspond to the same vertices in both
#' graphs. That is, the \eqn{n \times n}{n*n} permutation matrix that defines
#' the bijection is \eqn{I_{m} \bigoplus P} for a \eqn{(n-m)\times
#' (n-m)}{(n-m)*(n-m)} permutation matrix \eqn{P} and \eqn{m} times \eqn{m}
#' identity matrix \eqn{I_{m}}. The function \code{match_vertices} estimates
#' the permutation matrix \eqn{P} via an optimization algorithm based on the
#' Frank-Wolfe algorithm.
#'
#' See references for further details.
#'
# @aliases match_vertices seeded.graph.match
#' @param A a numeric matrix, the adjacency matrix of the first graph
#' @param B a numeric matrix, the adjacency matrix of the second graph
#' @param S a numeric matrix, the number of seeds x 2 matching vertex table.
#' If \code{S} is \code{NULL}, then it is using a \eqn{soft} seeding algorithm,
#' otherwise a \eqn{hard} seeding is used.
#' @param start a numeric matrix, the permutation matrix estimate is
#' initialized with \code{start}
#' @param pad a scalar value for padding
#' @param iteration The number of iterations for the Frank-Wolfe algorithm
#' @return A numeric matrix which is the permutation matrix that determines the
#' bijection between the graphs of \code{A} and \code{B}
#' @author Vince Lyzinski \url{http://www.ams.jhu.edu/~lyzinski/}
#' @references Vogelstein, J. T., Conroy, J. M., Podrazik, L. J., Kratzer, S.
#' G., Harley, E. T., Fishkind, D. E.,Vogelstein, R. J., Priebe, C. E. (2011).
#' Fast Approximate Quadratic Programming for Large (Brain) Graph Matching.
#' Online: \url{http://arxiv.org/abs/1112.5507}
#'
#' Fishkind, D. E., Adali, S., Priebe, C. E. (2012). Seeded Graph Matching
#' Online: \url{http://arxiv.org/abs/1209.0367}
#'
#' @export

sgm.ori <- function(A,B,start,S=NULL,pad=0,iteration=20){
    # forcing seeds to be first m vertices in both graphs
    m <- nrow(S)
    if(is.null(m)){
        m <- 0
    }else{
        Aseeds <- S[,1]
        Anotseeds <- setdiff(1:nrow(A),Aseeds)
        Aorder <- c(Aseeds,Anotseeds)
        Anew <- A[Aorder,Aorder]
        A <- Anew

        Bseeds <- S[,2]
        Bnotseeds <- setdiff(1:nrow(B),Bseeds)
        Border <- c(Bseeds,Bnotseeds)
        Bnew <- B[Border,Border]
        B <- Bnew
    }

    #seeds are assumed to be vertices 1:m in both graphs
    require('clue')
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
        for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
    }
    totv<-max(totv1,totv2)
    n<-totv-m
    if (m != 0){
        A12<-rbind(A[1:m,(m+1):(m+n)])
        A21<-cbind(A[(m+1):(m+n),1:m])
        B12<-rbind(B[1:m,(m+1):(m+n)])
        B21<-cbind(B[(m+1):(m+n),1:m])
    }
    if (m==0){
        A12<-matrix(0,n,n)
        A21<-matrix(0,n,n)
        B12<-matrix(0,n,n)
        B21<-matrix(0,n,n)
    }
    if (n==1){ A12=t(A12)
    A21=t(A21)
    B12=t(B12)
    B21=t(B21)}

    A22<-A[(m+1):(m+n),(m+1):(m+n)]
    B22<-B[(m+1):(m+n),(m+1):(m+n)]
    tol<-1
    P<-start
    toggle<-1
    iter<-0
    x<- A21%*%t(B21)
    y<- t(A12)%*%B12
    while (toggle==1 & iter<maxiter)
    {
        iter<-iter+1
        z<- A22%*%P%*%t(B22)
        w<- t(A22)%*%P%*%B22
        Grad<-x+y+z+w;
        mm=max(abs(Grad))
        ind<-matrix(solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
        T<-diag(n)
        T<-T[ind,]
        wt<-t(A22)%*%T%*%B22
        c<-sum(diag(w%*%t(P)))
        d<-sum(diag(wt%*%t(P)))+sum(diag(w%*%t(T)))
        e<-sum(diag(wt%*%t(T)))
        u<-sum(diag(t(P)%*%x+t(P)%*%y))
        v<-sum(diag(t(T)%*%x+t(T)%*%y))
        if( c-d+e==0 && d-2*e+u-v==0){
            alpha<-0
        }else{
            alpha<- -(d-2*e+u-v)/(2*(c-d+e))}
        f0<-0
        f1<- c-e+u-v
        falpha<-(c-d+e)*alpha^2+(d-2*e+u-v)*alpha
        if(alpha < tol && alpha > 0 && falpha > f0 && falpha > f1){
            P<- alpha*P+(1-alpha)*T
        }else if(f0 > f1){
            P<-T
        }else{
            toggle<-0}
    }
    D<-P
    corr<-matrix(solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(corr=corr[,2], P=P, D=D))
}



stfp <- function(A, dim = NULL, scaling = TRUE)
{
    suppressMessages(require(irlba))

    if (dim(A)[1]==dim(A)[2]) {
        L <- nonpsd.laplacian(A) ## diagonal augmentation
    } else {
        L <- A
    }

    if(is.null(dim)) {
        L.svd <- irlba(L)
        #        L.svd <- svd(L)
        dim <- getElbows(L.svd$d,3,plot=FALSE)[2] # ignore this, use "getElbows" instead
    } else {
        L.svd <- irlba(L,dim,dim)
        #        L.svd <- svd(L)
    }

    L.svd.values <- L.svd$d[1:dim]
    L.svd.vectors <- L.svd$v[,1:dim]

    if(scaling == TRUE) { # projecting to sphere
        if(dim == 1)
            L.coords <- sqrt(L.svd.values) * L.svd.vectors
        else
            L.coords <- L.svd.vectors %*% diag(sqrt(L.svd.values))
    }
    else {
        L.coords <- L.svd.vectors
    }

    return(list(X=L.coords,D=L.svd.values))
}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = F) {
    if (is.matrix(dat))
        d <- sort(apply(dat,2,sd), decreasing=T)
    else
        d <- sort(dat,decreasing=TRUE)

    if (!is.logical(threshold))
        d <- d[d > threshold]

    p <- length(d)
    if (p == 0)
        stop(paste("d must have elements that are larger than the threshold ",
                   threshold), "!", sep="")

    lq <- rep(0.0, p)                     # log likelihood, function of q
    for (q in 1:p) {
        mu1 <- mean(d[1:q])
        mu2 <- mean(d[-(1:q)])              # = NaN when q = p
        sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) / (p - 1 - (q < p))
        lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
            sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
    }

    q <- which.max(lq)
    if (n > 1 && q < p) {
        q <- c(q, q + getElbows(d[(q+1):p], n=n-1, plot=FALSE))
    }

    if (plot==TRUE) {
        if (is.matrix(dat)) {
            sdv <- d # apply(dat,2,sd)
            plot(sdv,type="b",xlab="dim",ylab="stdev")
            points(q,sdv[q],col=2,pch=19)
        } else {
            plot(dat, type="b")
            points(q,dat[q],col=2,pch=19)
        }
    }

    return(q)
}

nonpsd.laplacian <- function(A)
{
    A <- as.matrix(A)
    n <- nrow(A)
    s <- rowSums(A)/(n-1)
    diag(A) <- diag(A)+s
    return(A)
}

## for aLAP
# Help method of parrallelMatch
matchVertex <- function(s, candidate, mate, Q) {
    if (candidate[candidate[s]] == s) {
        mate[s] <- candidate[s];
        mate[candidate[s]] <- s;
        Q <- c(Q, s, candidate[s]);
    }
    return(list(Q, mate));
}

findMate <- function(s, graph, mate) {
    # Initialization
    max_wt <- -Inf;
    max_wt_id <- NaN;

    # Find the locally dominant vertices for one single vertex
    if (s <= dim(graph)[2] / 2) {
        for (i in (dim(graph)[2] / 2 + 1) : dim(graph)[2]) {
            if (is.nan(mate[i]) && max_wt < graph[s, i]) {
                max_wt <- graph[s, i];
                max_wt_id <- i;
            }
        }
    } else {
        for (i in 1 : (dim(graph)[2] / 2)) {
            if (is.nan(mate[i]) && max_wt < graph[s, i]) {
                max_wt <- graph[s, i];
                max_wt_id <- i;
            }
        }
    }
    return(max_wt_id);
}


# Function to solve Linear Assignment Problem
# Using algorithm in A multithreaded algorithm for network alignment
# via approximate matching
# graph has the form: [0 C; C' 0] where C is the cost function
# Output mate is a permutation
# Written by Ao Sun and Lingyao Meng
parallelMatch <- function(graph) {
    ##require('foreach')
    ##source('findMate.R');
    ##source('matchVertex.R');
    ## adjacency matrix corresponding to graph G
    numVertices <- nrow(graph)

    ## Initialize the matching vector with NaN
    mate <- rep(NaN, numVertices)

    ## Initialize the other parameters
    candidate <- rep(0, numVertices);
    qC <- vector(); qN <- vector();

    #--------------------------------- Phase 1 ------------------------------------------
    candidate <- sapply(seq_len(numVertices), function(x) findMate(x,graph,mate))

    for (j in 1 : numVertices) {
        temp <- matchVertex(j, candidate, mate, qC);
        qC <- temp[[1]]; mate <- temp[[2]];
    }
    #--------------------------------- Phase 2 ------------------------------------------
    repeat {
        for (k in 1 : length(qC)) {
            ## Return to the index of adjacent Vertex
            if (qC[k] <= (numVertices / 2)) {

                for (h in (numVertices / 2 + 1) : numVertices) {
                    if ((candidate[h] == qC[k]) && (h != mate[qC[k]])) {
                        candidate[h] <- findMate(h, graph, mate);
                        temp2 <- matchVertex(h, candidate, mate, qN);
                        qN <- unlist(temp2[1]); mate <- unlist(temp2[2]);
                    }
                }
            } else {
                for (h in 1 : (numVertices / 2)) {
                    if ((candidate[h] == qC[k]) && h != mate[qC[k]]) {
                        candidate[h] <- findMate(h, graph, mate);
                        temp3 <- matchVertex(h, candidate, mate, qN);
                        qN <- unlist(temp3[1]); mate <- unlist(temp3[2]);
                    }
                }
            }
        }

        qC <- qN;
        qN <- vector();
        if (length(qC) == 0) {
            break;
        }
    }
    mate <- mate[(numVertices / 2 + 1) : length(mate)];
    return(mate)
}

