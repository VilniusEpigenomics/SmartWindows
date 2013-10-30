\name{apply.PTA.result}
\alias{apply.PTA.result}

\title{apply.PTA.result}
\description{Apply ranges found with PTA to a different score matrix}
\usage{
    apply.PTA.result(result, start, end, scores)
}
\arguments{
    \item{result}{result of previous PTA.raw application.}
    \item{start, end}{vectors of original ranges (should be the same that were
        passed to PTA.raw).}
    \item{scores}{matrix of scores with rows corresponding to ranges.}
}
\value{
    Returns a score matrix with the new ranges applied.
}
\seealso{
    \link[PTA]{PTA-package}
}
