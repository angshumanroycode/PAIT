\name{T_phi}
\alias{T_phi}
\title{
Proximity adjusted independence test based on phi coefficient.
}
\description{
Carries out a proximity adjusted independence test based on phi coefficient.
}
\usage{
T_phi(x, y, B = 200)
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
T_phi(x, y, B = 200)
}
