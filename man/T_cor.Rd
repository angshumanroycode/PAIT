\name{T_cor}
\alias{T_cor}
\title{
Proximity adjusted independence test based on correlation coefficient.
}
\description{
Carries out a proximity adjusted independence test based on correlation coefficient.
}
\usage{
T_cor(x, y, B = 200)
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
T_cor(x, y, B = 200)
}
