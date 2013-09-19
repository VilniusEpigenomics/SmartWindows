\name{adapt.ranged}
\alias{adapt.ranged}
\alias{adapt.ranged.raw}

\title{adapt.ranged}
\description{Adapt ranged data to a different set of ranges.}
\usage{
adapt.ranged(ranged.data, dest.ranges, na.rm=FALSE, add.error=FALSE)
adapt.ranged.raw(src.start, src.end, src.score, dest.start, dest.end)
}
\arguments{
    \item{ranged.data}{either a RangedData or GRanges instance with a \code{score} column.}
    \item{src.start, src.end,  src.score}{source data.}
    \item{dest.ranges, dest.start, dest.end}{destination ranges.}
    \item{na.rm}{whether to remove ranges which don't have any relevant sources for score, or leave them as NAs.}
    \item{add.error}{add an \code{error} column which contains the sum squared error of each range's adaptation.}
}
\value{
\code{adapt.ranged} returns a RangedData instance with the adapted data.

\code{adapt.ranged.raw} returns a list with items \code{start}, \code{end} and \code{score}.
}
\seealso{
\link[PTA]{PTA-package}
}
