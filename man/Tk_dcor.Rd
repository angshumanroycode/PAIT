\name{Tk_dcor}
\alias{Tk_dcor}
\title{
Vector of average distance correlation coefficients for different sized neighborhoods.
}
\description{
Provides a vector of average distance correlation coefficients for different sized neighborhoods.
}
\usage{
Tk_dcor(x, y)
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
A numeric vector containing average value of distance correlation coefficients for different sized neighborhoods.
}
\examples{
x = rnorm(10)
y = rnorm(10)
Tk_dcor(x, y)
}
