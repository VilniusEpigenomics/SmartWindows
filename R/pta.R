library(IRanges)

PTA <- function(data, count=1, error=Inf, adjacency.treshold=1, skip=0, space=1, mode=c("normal", "correlation")) {
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
    for (i in 1:ncol(df)) {
        d.scores[, i] <- df[, i]
    }

    mode <- match.arg(mode)
    mode.int <- switch(mode, normal=0, correlation=1)

    result <- .Call("PTA",
                    d.start, d.end, d.scores,
                    count, error, adjacency.treshold, skip, mode.int,
                    PACKAGE="pta")
    result
}
