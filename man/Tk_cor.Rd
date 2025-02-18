\name{Tk_cor}
\alias{Tk_cor}
\title{
Vector of average correlation coefficients for different sized neighborhoods.
}
\description{
Provides a vector of average correlation coefficients for different sized neighborhoods.
}
\usage{
Tk_cor(x, y)
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
A numeric vector containing average value of correlation coefficients for different sized neighborhoods.
}
\examples{
x = rnorm(10)
y = rnorm(10)
Tk_cor(x, y)
}
