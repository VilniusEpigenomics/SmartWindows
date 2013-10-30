\name{PTA}
\alias{PTA}
\alias{PTA.raw}

\title{PTA}
\description{Aggregate ranged data in various ways}
\usage{
    PTA(data, space=1, ...)
    PTA.raw(start, end, scores, count.bound=1,
            error.bound=Inf, cumulative.error.bound=Inf,
            adjacency.threshold=1, skip=0, mode=c("normal", "correlation"),
            correlation.bound=-1, correlation.spearman=FALSE,
            correlation.absolute=TRUE)
}
\arguments{
    \item{data}{either a RangedData or GRanges instance with the data to be
        processed.}
    \item{space}{the space or seqlevel (in case of GRanges) to process. The
        default is the first one.}
    \item{...}{these arguments are shared with PTA.raw.}
    \item{start, end}{vectors of beginnings and ends of ranges. The ranges
        should be sorted and non-overlapping (see \code{\link{deoverlap}}).}
    \item{scores}{matrix of values, the i-th row of which is a vector
        corresponding to the i-th range}
    \item{adjacency.threshold}{ignore gaps of this size when considering whether
        ranges are adjacent.}
    \item{skip}{allow skipping over a given amount of ranges, ignoring the
        error they introduce. Useful if data has outliers that need to be ignored.
        Skipping still happens no further than adjacency threshold, which
        should be increased together with this parameter.}
    \item{count.bound}{stop processing when this count of ranges is reached.}
    \item{mode}{either 'normal' (default) which use sum-squared-error as
        dissimilarity metric, or 'correlation' which uses correlation as
        similarity}
    \item{error.bound}{stop processing when per-merge sum-squared-error reaches
        this value. Only used in normal mode.}
    \item{cumulative.error.bound}{stop processing when the cumulative error
        reaches this value, expressed as a fraction of maximum merge error for the
        given data. Only used in normal mode.}
    \item{correlation.bound}{stop processing when the maximum correlation
        between nodes considered for merging is lower than this number.}
    \item{correlation.spearman}{use Spearman correlation.}
    \item{correlation.absolute}{if this is true, values are considered similar
        whether they are directly or inversely correlated.}
}
\value{
    \code{PTA} returns either RangedData or GRanges, depending on what was
    passed, containing processed data and only the space given with
    \code{space}.

    \code{PTA.raw} returns a list with items \code{start}, \code{end} and
    \code{scores} which contain processed data. The list also contains the
    cumulative error reached (excluding nodes skipped in merges when
    \code{skip} > 0) in item \code{cumulative.error}.
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
PTA(d, adjacency.threshold=10, cumulative.error.bound=0.1)

# correlation mode
PTA(d, adjacency.threshold=10, mode="correlation", correlation.bound=0.8)
}
