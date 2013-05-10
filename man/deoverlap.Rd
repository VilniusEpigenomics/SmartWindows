\name{deoverelap}
\alias{deoverlap}
\alias{deoverlap.raw}

\title{deoverlap}
\description{Remove simple overlaps from ranged data.}
\usage{
deoverlap(x)

deoverlap.raw(start, end)
}
\arguments{
    \item{x}{ranged data object.}
    \item{start, end}{vectors of beginnings and ends of ranges.}
}
\value{
\code{deoverlap} finds overlaps between pairs of adjacent ranges. It takes the
overlapping ends and sets them both to the mean of their position. This is
useful for passing mTAG-chip probes to PTA. It returns \code{x} with the
modified ranges

\code{deoverlap.raw} takes vectors describing the ranges and returns a list with
items \code{start} and \code{end} containing the result.
}
\seealso{
\link[PTA]{PTA-package}
}
