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

sgm <- function (A,B,seeds,hard=TRUE,start="barycenter",maxiter=20){
    gamma <- 0.1
    if(is.null(seeds)){
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

    P <- sgm.ordered(AA,BB,m,S,pad=0,maxiter)
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

sgm.ordered <- function(A,B,m,start,pad=0,maxiter=20){
    #seeds are assumed to be vertices 1:m in both graphs
    suppressMessages(library(clue))
    totv1<-ncol(A)
    totv2<-ncol(B)
    if(totv1>totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv1-totv2
        B.org <- B
        B <- A
        B[1:totv2,1:totv2] <- B.org
        B[-(1:totv2),-(1:totv2)] <- pad
#        for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    }else if(totv1<totv2){
        A[A==0]<- -1
        B[B==0]<- -1
        diff<-totv2-totv1
        A.org <- A
        A <- B
        A[1:totv1,1:totv1] <- A.org
        A[-(1:totv1),-(1:totv1)] <- pad
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
        ind<-matrix(solve_LSAP(Grad+matrix(mm,totv-m,totv-m), maximum =TRUE))
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
    corr<-matrix(solve_LSAP(P, maximum = TRUE))
    P=diag(n)
    P=rbind(cbind(diag(m),matrix(0,m,n)),cbind(matrix(0,n,m),P[corr,]))
    corr<-cbind(matrix((m+1):totv, n),matrix(m+corr,n))
    return(list(A=A22, B=B22, corr=corr[,2], P=P, D=D, iter=iter))
}

#' @export
sgm.ordered.ori <- function(A,B,m,start,pad=0,maxiter=20){
    #seeds are assumed to be vertices 1:m in both graphs
    suppressMessages(library(clue))
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
    return(list(A=A22, B=B22, corr=corr[,2], P=P, D=D, iter=iter))
}

#' @export
sgm.ordered.sparse <- function(A,B,m,start,pad=0,maxiter=20){
  #seeds are assumed to be vertices 1:m in both graphs
  totv1<-ncol(A)
  totv2<-ncol(B)
  if(totv1>totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv1-totv2
#    for (j in 1:diff){B<-cbind(rbind(B,pad),pad)}
    B.org <- B
    B <- A
    B[1:totv2,1:totv2] <- B.org
    B[-(1:totv2),-(1:totv2)] <- pad
  }else if(totv1<totv2){
    A[A==0]<- -1
    B[B==0]<- -1
    diff<-totv2-totv1
#    for (j in 1:diff){A<-cbind(rbind(A,pad),pad)}
    A.org <- A
    A <- B
    A[1:totv1,1:totv1] <- A.org
    A[-(1:totv1),-(1:totv1)] <- pad
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
  return(list(A=A22, B=B22, corr=corr[,2], P=P, D=D, iter=iter))
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

sgm2 <- function(A,B,start,S=NULL,pad=0,iteration=20){
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
