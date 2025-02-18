\name{T_dcor}
\alias{T_dcor}
\title{
Proximity adjusted independence test based on distance correlation coefficient.
}
\description{
Carries out a proximity adjusted independence test based on distance correlation coefficient.
}
\usage{
T_dcor(x, y, B = 200)
}
\arguments{
  \item{x}{
A numeric vector containing the values from the random variable X.
}
  \item{y}{
A numeric vector containing the values from the random variable Y.
}
  \item{B}{
Number of permutations.
}
}

\value{
A numeric value which is the p-value for this test.
}
\examples{
x = rnorm(10)
y = rnorm(10)
T_dcor(x, y, B = 200)
}
