\name{PTA-package}
\alias{PTA-package}
\encoding{UTF-8}
\docType{package}
\title{
Parsimonious temporal aggregation
}
\description{
Efficiently perform aggregation on ranged data based on various metrics.
}
\details{
\tabular{ll}{
Package: \tab PTA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-01-28\cr
License: \tab BSD\cr
}

The original method described in the paper is meant for simplifying temporal
ranged data, but this package supports other metrics for merging, like
correlation and absolute value.
}
\author{
Kiprianas Spiridonovas

Maintainer: Kiprianas Spiridonovas <k.spiridonovas@gmail.com>
}
\references{
Gordevičius, Juozas, Johann Gamper, and Michael Böhlen. "Parsimonious temporal aggregation." The VLDB Journal 21.3 (2012): 309-332.
}
\keyword{ package }
\seealso{
\link{IRanges}, \link{GenomicRanges}
}
\examples{
set.seed(123)
PTA(RangedData(ranges=IRanges(start=seq(0, 99, 10), width=11), score=rnorm(10)),
    error.bound=0.1)
}
