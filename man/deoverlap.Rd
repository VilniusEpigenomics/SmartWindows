\name{deoverelap}
\alias{deoverlap}

\title{deoverlap}
\description{Remove simple overlaps from ranged data.}
\usage{
deoverlapRanges(x)

deoverlap(start, end)
}
\arguments{
    \item{x}{a GRanges object.}
    \item{start, end}{vectors of beginnings and ends of ranges.}
}
\value{
\code{deoverlapRanges} finds overlaps between pairs of adjacent ranges. It takes the
overlapping ends and sets them both to the mean of their position. This is
useful for passing mTAG-chip probes to SmartWindows. It returns \code{x} with the
modified ranges

\code{deoverlap} takes vectors describing the ranges and returns a list with
items \code{start} and \code{end} containing the result.
}
\seealso{
\link[SmartWindows]{SmartWindows-package}
}
