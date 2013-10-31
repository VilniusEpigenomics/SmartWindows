library(IRanges)

#' @useDynLib PTA

#' @export
PTA <- function(data, space=1, ...) {
    d.start <- start(data)
    d.end <- end(data)

    if (class(data) == "GRanges") {
        if (is.numeric(space)) {
            space <- as.vector(GenomicRanges::seqlevels(data)[space])
        }
        df <- GenomicRanges::mcols(data[seqnames(data) == space])
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
        r <- GRanges(ranges=ranges, seqinfo=GenomicRanges::seqinfo(data), seqnames=space)
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

#' @export
PTA.raw <- function(start, end, scores,
                    count.bound=1, error.bound=Inf, cumulative.error.bound=Inf,
                    adjacency.threshold=1, skip=0, mode=c("normal", "correlation"),
                    correlation.bound=-1, correlation.spearman=FALSE, correlation.absolute=TRUE,
                    sample.parameter=numeric(), sample.parameter.weight=0.5) {
    mode <- match.arg(mode)
    mode <- switch(mode, normal=0, correlation=1)

    if (is.vector(scores)) {
        scores <- matrix(scores)
    }

    if (length(start) != length(end)) {
        stop("Start and end counts differ.")
    }

    if (length(start) != nrow(scores)) {
        stop("Range and score counts differ")
    }

    sample.parameter.given <- length(sample.parameter) > 0
    if (sample.parameter.given && length(sample.parameter) != ncol(scores)) {
        stop("There should be as many sample.parameter values as there are samples.")
    }

    arguments <- as.list(environment())

    result <- .Call("PTA", arguments, PACKAGE="PTA")

    colnames(result$scores) <- colnames(scores)

    result
}

#' @export
apply.PTA.result <- function(result, start, end, scores) {
    stopifnot(length(result$groups) == nrow(scores))
    len <- end - start
    s <- len * if (!is.null(result$coefficients)) {
        result$coefficients * scores + result$intercept
    } else {
        scores
    }

    filter <- result$groups != -1
    s <- s[filter,]
    groups <- result$groups[filter]

    grouplen <- tapply(len, groups, sum)

    x <- apply(s, 2,
               function(col) {
                   tapply(col, groups, sum) / grouplen
                })
    x
}

#' @export
deoverlap <- function(x) {
    result <- deoverlap.raw(start(x), end(x))
    start(x) <- result$start
    end(x) <- result$end
    x
}

#' @export
deoverlap.raw <- function(start, end) {
    .Call("deoverlap", start, end)
}

#' @export
adapt.ranged <- function(ranged.data, dest.ranges, na.rm=FALSE, add.error=FALSE) {
    destcount <- nrow(dest.ranges)
    if (is.null(destcount)) destcount <- length(dest.ranges)
    srccount <- nrow(ranged.data)
    if (is.null(srccount)) srccount <- length(ranged.data)

    destscore <- if (is.null(score(dest.ranges))) numeric() else score(dest.ranges)
    sum <- rep(0, destcount)
    len <- rep(0, destcount)
    found <- rep(FALSE, destcount)
    error <- rep(0, 0)
    .Call("adapt_ranged",
          sum, len, found, error, FALSE,
          destcount, start(dest.ranges), end(dest.ranges),
          destscore,
          srccount, start(ranged.data), end(ranged.data),
          if (is.null(score(ranged.data))) numeric() else score(ranged.data))
    destscore <- sum / len

    if (add.error) {
        sum <- rep(0, destcount)
        len <- rep(0, destcount)
        found <- rep(FALSE, destcount)
        error <- rep(0, destcount)
        .Call("adapt_ranged",
              sum, len, found, error, TRUE,
              destcount, start(dest.ranges), end(dest.ranges),
              destscore,
              srccount, start(ranged.data), end(ranged.data),
              if (is.null(score(ranged.data))) numeric() else score(ranged.data))
    }

    ranges <- IRanges(start=start(dest.ranges), end=end(dest.ranges))
    newscore <- sum / len
    if (na.rm) {
        newscore <- newscore[found]
        error <- error[found]
        ranges <- ranges[found]
    } else {
        newscore[!found] <- NA
    }
    if (add.error) {
        RangedData(ranges=ranges, score=newscore, error=error)
    } else {
        RangedData(ranges=ranges, score=newscore)
    }
}

#' @export
adapt.ranged.raw <- function(src.start, src.end, src.score, dest.start, dest.end) {
    dest.count <- length(dest.start)
    src.count <- length(src.start)
    if (is.matrix(src.score)) {
        newscore <- apply(src.score, 2,
                          function(score) {
                              sum <- rep(0, dest.count)
                              len <- rep(0, dest.count)
                              found <- rep(FALSE, dest.count)
                              error <- rep(0, 0)
                              .Call("adapt_ranged", sum, len, found, error, FALSE,
                                    dest.count, dest.start, dest.end, rep(0, dest.count),
                                    src.count, src.start, src.end, score)
                              ifelse(found, sum / len, NA)
                          })
    } else {
        sum <- rep(0, dest.count)
        len <- rep(0, dest.count)
        found <- rep(FALSE, dest.count)
        error <- rep(0, 0)
        .Call("adapt_ranged", sum, len, found, error, FALSE,
              dest.count, dest.start, dest.end, rep(0, dest.count),
              src.count, src.start, src.end, src.score)
        newscore <- ifelse(found, sum / len, NA)
    }
    list(start=dest.start,
         end=dest.end,
         score=newscore)
}
