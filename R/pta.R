library(IRanges)

PTA <- function(data, count=1, error=Inf, adjacency.treshold=1){
    result <- .Call("PTA", start(data), end(data), score(data), count, error, adjacency.treshold, PACKAGE = "pta")
    RangedData(IRanges(start=result$start, end=result$end), score=result$score)
}

gPTA <- function(data, count=1, error=Inf, adjacency.treshold=1){
    result <- .Call("gPTA", start(data), end(data), score(data), count, error, adjacency.treshold, PACKAGE = "pta")
    RangedData(IRanges(start=result$start, end=result$end), score=result$score)
}

