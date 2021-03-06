% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rdpg.R
\name{rdpg.sample.correlated}
\alias{rdpg.sample.correlated}
\title{Samples a pair of correlated \eqn{G(X)} random dot product graphs}
\usage{
rdpg.sample.correlated(X, rho)
}
\arguments{
\item{X}{\eqn{n \times d} random matrix}

\item{rho}{Numeric scalar in the unit interval, the target Pearson correlation between
the adjacency matrieces of the original and the generated graph.}
}
\value{
A list of two igraph objects, named \code{A} and \code{B}, which are
two graphs whose adjacency matrix entries are correlated with \code{rho}.
}
\description{
Samples a pair of correlated \eqn{G(X)} random dot product graphs
}
\author{
Minh Tang <scent.of.time@gmail.com>
}
