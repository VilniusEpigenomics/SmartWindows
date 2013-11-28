\name{adaptRanged}
\alias{adaptRanged}

\title{adaptRanged}
\description{Adapt ranged data to a different set of ranges.}
\usage{
adaptRanged(srcStart, srcEnd, srcScore, destStart, destEnd)
}
\arguments{
    \item{srcStart, srcEnd,  srcScore}{source data. \code{srcScore} can be a matrix where rows correspond to ranges.}
    \item{destRanges, destStart, destEnd}{destination ranges.}
}
\value{
\code{adaptRanged} returns a list with items \code{start}, \code{end} and \code{score}.
}
\seealso{
\link[PTA]{PTA-package}
}
