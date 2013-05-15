\name{adapt.ranged}
\alias{adapt.ranged}

\title{adapt.ranged}
\description{Adapt ranged data to a different set of ranges.}
\usage{
adapt.ranged(ranged.data, dest.ranges, na.rm=FALSE, add.error=FALSE)
}
\arguments{
    \item{ranged.data}{either a RangedData or GRanges instance with a \code{score} column.}
    \item{dest.ranges}{destination ranges.}
    \item{na.rm}{whether to remove ranges which don't have any relevant sources for score, or leave them as NAs.}
    \item{add.error}{add an \code{error} column which contains the sum squared error of each range's adaptation.}
}
\value{
Returns a RangedData instance with the adapted data.
}
\seealso{
\link[PTA]{PTA-package}
}
