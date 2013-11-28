#' @useDynLib PTA

extractGRanges <- function(data, seqname=NULL) {
    if (is.null(seqname)) {
        seqname <- GenomicRanges::seqlevels(data)[[1]]
    }
    d <- list()
    d$start <- GenomicRanges::start(data)
    d$end <- GenomicRanges::end(data)
    s <- GenomicRanges::mcols(data[seqnames(data) == seqname])
    d$scores <- as.matrix(as.data.frame(s))
    d
}

#' @export
PTA <- function(start=NULL, end=NULL, scores=NULL, data=NULL, chr=NULL,
                countBound=1, errorBound=Inf, cumulativeErrorBound=Inf,
                adjacencyThreshold=1, skip=0, mode=c("normal", "correlation"),
                correlationBound=-1, correlationSpearman=FALSE, correlationAbsolute=TRUE,
                individualParameter=numeric(), individualParameterWeight=0.5) {
    if (class(start) == "GRanges") {
        data <- extractGRanges(start, chr)
    }

    if (!is.null(data)) {
        start <- data$start
        end <- data$end
        scores <- data$scores
    }

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

    individualParameterGiven <- length(individualParameter) > 0
    if (individualParameterGiven && length(individualParameter) != ncol(scores)) {
        stop("There should be as many individualParameter values as there are samples.")
    }

    arguments <- as.list(environment())

    result <- .Call("PTA", arguments, PACKAGE="PTA")

    colnames(result$scores) <- colnames(scores)

    result$groups <- result$groups + 1

    result
}

#' @export
applyPTAResult <- function(result, start, end, scores) {
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
    len <- len[filter]

    grouplen <- tapply(len, groups, sum)

    x <- apply(s, 2,
               function(col) {
                   tapply(col, groups, sum) / grouplen
                })
    x
}

#' @export
deoverlapRanges <- function(x) {
    result <- deoverlap(IRanges::start(x), IRanges::end(x))
    IRanges::start(x) <- result$start
    IRanges::end(x) <- result$end
    x
}

#' @export
deoverlap <- function(start, end) {
    .Call("deoverlap", start, end)
}

#' @export
adaptRanged <- function(srcStart, srcEnd, srcScore, destStart, destEnd) {
    destCount <- length(destStart)
    srcCount <- length(srcStart)
    if (is.matrix(srcScore)) {
        newscore <- apply(srcScore, 2,
                          function(score) {
                              sum <- rep(0, destCount)
                              len <- rep(0, destCount)
                              found <- rep(FALSE, destCount)
                              error <- rep(0, 0)
                              .Call("adapt_ranged", sum, len, found, error, FALSE,
                                    destCount, destStart, destEnd, rep(0, destCount),
                                    srcCount, srcStart, srcEnd, score)
                              ifelse(found, sum / len, NA)
                          })
    } else {
        sum <- rep(0, destCount)
        len <- rep(0, destCount)
        found <- rep(FALSE, destCount)
        error <- rep(0, 0)
        .Call("adapt_ranged", sum, len, found, error, FALSE,
              destCount, destStart, destEnd, rep(0, destCount),
              srcCount, srcStart, srcEnd, srcScore)
        newscore <- ifelse(found, sum / len, NA)
    }
    list(start=destStart,
         end=destEnd,
         score=newscore)
}
