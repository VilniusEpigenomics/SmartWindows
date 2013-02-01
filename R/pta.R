library(IRanges)

PTA <- function(data, count=0, error=Inf){
    if ((count == 0) && (error == Inf)) {
        stop("Either count or error bound should be specified.")
    }
    result <- .Call("PTA", start(data), end(data), data[["score"]], count, error, PACKAGE = "pta")
    RangedData(IRanges(start=result$start, end=result$end), score=result$score)
}

gPTA <- function(data, count=0, error=Inf){
    if ((count == 0) && (error == Inf)) {
        stop("Either count or error bound should be specified.")
    }
    result <- .Call("gPTA", start(data), end(data), data[["score"]], count, error, PACKAGE = "pta")
    RangedData(IRanges(start=result$start, end=result$end), score=result$score)
}

