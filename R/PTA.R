library(IRanges)

PTA <- function(data,
                count=1, error=Inf, adjacency.treshold=1, skip=0, space=1, mode=c("normal", "correlation"), correlation.bound=0) {
    d.start <- start(data)
    d.end <- end(data)

    if (class(data) == "GRanges") {
        if (is.numeric(space)) {
            space <- as.vector(seqlevels(data)[space])
        }
        df <- values(data[seqnames(data) == space])
    } else if (class(data) == "RangedData") {
        df <- values(data)[[space]]
    } else {
        error("Unsupported data type")
    }

    d.scores <- matrix(nrow=nrow(df), ncol=ncol(df))
    colnames(d.scores) <- colnames(df)
    for (i in 1:ncol(df)) {
        d.scores[, i] <- df[, i]
    }
    rm(df)

    result <- PTA.raw(d.start, d.end, d.scores, count, error, adjacency.treshold, skip, space, mode, correlation.bound)

    ranges <- IRanges(start=result$start, end=result$end)
    if (class(data) == "GRanges") {
        GRanges(ranges=ranges, seqinfo=seqinfo(data), seqnames=space, mcols=as.data.frame(result$scores))
    } else if (class(data) == "RangedData") {
        r <- RangedData(ranges=ranges, space=space)
        for (col in colnames(result$scores)) {
            values(r)[[1]][[col]] <- result$scores[, col]
        }
        r
    }
}

PTA.raw <- function(start, end, scores,
                    count=1, error=Inf, adjacency.treshold=1, skip=0, space=1, mode=c("normal", "correlation"), correlation.bound=0) {
    mode <- match.arg(mode)
    mode.int <- switch(mode, normal=0, correlation=1)

    result <- .Call("PTA",
                    start, end, scores,
                    count, error, adjacency.treshold, skip, mode.int, correlation.bound,
                    PACKAGE="PTA")

    colnames(result$scores) <- colnames(scores)

    result
}
