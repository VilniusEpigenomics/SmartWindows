\name{PTA}
\alias{PTA}

\title{PTA}
\description{Aggregate ranged data in various ways}
\usage{
    PTA(granges, ...)

    PTA(data, ...)
    
    PTA(start, end, scores, chr=NULL,
        countBound=1, errorBound=Inf, cumulativeErrorBound=Inf,
        adjacencyThreshold=1, skip=0, mode=c("normal", "correlation"),
        correlationBound=-1, correlationSpearman=FALSE,
        correlationAbsolute=TRUE,
        individualParameter=numeric(), individualParameterWeight=0.5)
}
\arguments{
    \item{granges}{a GRanges object.}
    \item{start, end, chr}{vectors of beginnings, ends and possibly chromosome
        (or other kind of group) of ranges. The ranges should be sorted and
        non-overlapping (see \code{\link{deoverlap}}).}
    \item{scores}{matrix of values, the i-th row of which is a vector
        corresponding to the i-th range.}
    \item{data}{a list containing elements named start, end, scores and chr
        corresponding to the arguments.}
    \item{adjacencyThreshold}{ignore gaps of this size when considering whether
        ranges are adjacent.}
    \item{skip}{allow skipping over a given amount of ranges, ignoring the
        error they introduce. Useful if data has outliers that need to be ignored.
        Skipping still happens no further than adjacency threshold, which
        should be increased together with this parameter.}
    \item{countBound}{stop processing when this count of ranges is reached.}
    \item{mode}{either 'normal' (default) which uses sum-squared-error as
        dissimilarity metric, or 'correlation' which uses correlation as
        similarity}
    \item{errorBound}{stop processing when per-merge sum-squared-error reaches
        this value. Only used in normal mode.}
    \item{cumulativeErrorBound}{stop processing when the cumulative error
        reaches this value, expressed as a fraction of approximate maximum
        merge error for the given data. Only used in normal mode.}
    \item{correlationBound}{stop processing when the maximum correlation
        between nodes considered for merging is lower than this number.}
    \item{correlationSpearman}{use Spearman correlation.}
    \item{correlationAbsolute}{if this is true, values are considered similar
        whether they are directly or inversely correlated.}
    \item{individualParameter, individualParameterWeight}{if this is given,
        correlation is calculated not in the usual way, but using the formula 
        \deqn{\mathrm{cor}(A, B)^{1-w} (\mathrm{cor}(A, p)^2 + \mathrm{cor}(B, p)^2)^{w / 2}}{%
              cor(A, B)^(1-w) (cor(A, p)^2 + cor(B, p)^2)^(w / 2)}
        where p is the parameter vector, A and B are scores for the regions that are being compared and w is the weight.

        This takes into account the ranges' correlation with an external parameter, constant for each individual.}
}
\value{
    \code{PTA} returns a list. The items \code{start}, \code{end}, \code{chr}
    and \code{scores} which contain processed data.

    The item \code{groups} contains the index (1-based) of the merged range
    that each original range belongs to. An index of 0 means that the range was
    skipped.

    In normal mode the list also contains the cumulative error reached
    (excluding nodes skipped in merges when \code{skip} > 0) in item
    \code{cumulativeError}.

    In correlation mode, the list also contains \code{coefficients} and
    \code{intercepts} describing approximately the transforms from
    original-value-space to merged-value-space.
}
\seealso{
    \link[PTA]{PTA-package}
}
\examples{
library(GenomicRanges)
set.seed(123)

# generate data
x <- rnorm(20)
y <- x + 0.5 * rnorm(20)
z <- 0.3 * x - 0.3 * y + 0.4 * rnorm(20)
d <- GRanges(seqnames="chr1", ranges=IRanges(start=(1:20)*10, width=6), x=x, y=y, z=z)

# normal mode
PTA(d, adjacencyThreshold=10, cumulativeErrorBound=0.1)

# correlation mode
PTA(d, adjacencyThreshold=10, mode="correlation", correlationBound=0.8)
}
