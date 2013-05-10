library(IRanges)

PTA <- function(data, space=1, ...) {
    d.start <- start(data)
    d.end <- end(data)

    if (class(data) == "GRanges") {
        if (is.numeric(space)) {
            space <- as.vector(seqlevels(data)[space])
        }
        df <- mcols(data[seqnames(data) == space])
    } else if (class(data) == "RangedData") {
        df <- values(data)[[space]]
    } else {
        stop("Unsupported data type")
    }

    d.scores <- matrix(nrow=nrow(df), ncol=ncol(df))
    colnames(d.scores) <- colnames(df)
    for (i in 1:ncol(df)) {
        d.scores[, i] <- df[, i]
    }
    rm(df)

    result <- PTA.raw(d.start, d.end, d.scores, ...)

    ranges <- IRanges(start=result$start, end=result$end)
    if (class(data) == "GRanges") {
        r <- GRanges(ranges=ranges, seqinfo=seqinfo(data), seqnames=space)
        for (col in colnames(result$scores)) {
            values(r)[[col]] <- result$scores[, col]
        }
        r
    } else if (class(data) == "RangedData") {
        r <- RangedData(ranges=ranges, space=space)
        for (col in colnames(result$scores)) {
            values(r)[[1]][[col]] <- result$scores[, col]
        }
        r
    }
}

PTA.raw <- function(start, end, scores,
                    count.bound=1, error.bound=Inf, adjacency.treshold=1, skip=0, mode=c("normal", "correlation"), correlation.bound=0) {
    mode <- match.arg(mode)
    mode.int <- switch(mode, normal=0, correlation=1)

    result <- .Call("PTA",
                    start, end, scores,
                    count.bound, error.bound, adjacency.treshold, skip, mode.int, correlation.bound,
                    PACKAGE="PTA")

    colnames(result$scores) <- colnames(scores)

    result
}

deoverlap <- function(x) {
    result <- deoverlap.raw(start(x), end(x))
    start(x) <- result$start
    end(x) <- result$end
    x
}

deoverlap.raw <- function(start, end) {
    .Call("deoverlap", start, end)
}
