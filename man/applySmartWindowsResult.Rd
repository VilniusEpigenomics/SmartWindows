\name{applySmartWindowsResult}
\alias{applySmartWindowsResult}

\title{applySmartWindowsResult}
\description{Apply ranges found with SmartWindows to a different score matrix}
\usage{
    applySmartWindowsResult(result, start, end, scores, chr=NULL)
}
\arguments{
    \item{result}{result of previous SmartWindows application.}
    \item{start, end, chr}{vectors of original ranges (should be the same that were
        passed to SmartWindows).}
    \item{scores}{matrix of scores with rows corresponding to ranges.}
}
\value{
    Returns a score matrix with the new ranges applied.
}
\seealso{
    \link[SmartWindows]{SmartWindows-package}
}
