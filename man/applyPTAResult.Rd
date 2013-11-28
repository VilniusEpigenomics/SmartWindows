\name{applyPTAResult}
\alias{applyPTAResult}

\title{applyPTAResult}
\description{Apply ranges found with PTA to a different score matrix}
\usage{
    applyPTAResult(result, start, end, scores)
}
\arguments{
    \item{result}{result of previous PTA application.}
    \item{start, end}{vectors of original ranges (should be the same that were
        passed to PTA).}
    \item{scores}{matrix of scores with rows corresponding to ranges.}
}
\value{
    Returns a score matrix with the new ranges applied.
}
\seealso{
    \link[PTA]{PTA-package}
}
