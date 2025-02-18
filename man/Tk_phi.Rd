\name{Tk_phi}
\alias{Tk_phi}
\title{
Vector of average phi coefficients for different sized neighborhoods.
}
\description{
Provides a vector of average phi coefficients for different sized neighborhoods.
}
\usage{
Tk_phi(x, y)
}
\arguments{
  \item{x}{
A numeric vector containing the values from the random variable X.
}
  \item{y}{
A numeric vector containing the values from the random variable Y.
}
}
\value{
A numeric vector containing average value of phi coefficients for different sized neighborhoods.
}
\examples{
x = rnorm(10)
y = rnorm(10)
Tk_phi(x, y)
}
