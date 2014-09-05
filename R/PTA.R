#' @useDynLib PTA

extractGRanges <- function(data) {
    d <- list()
    d$start <- GenomicRanges::start(data)
    d$end <- GenomicRanges::end(data)
    d$chr <- BiocGenerics::as.vector(GenomicRanges::seqnames(data))
    d$scores <- IRanges::as.matrix(GenomicRanges::mcols(data))
    d
}

#' @export
PTA <- function(start=NULL, end=NULL, scores=NULL, chr=NULL, data=NULL, ...) {
    if (class(data) == "GRanges") {
        data <- extractGRanges(data)
    }

    if (class(start) == "GRanges") {
        data <- extractGRanges(start)
    }

    if (!is.null(data)) {
        start <- data$start
        end <- data$end
        scores <- data$scores
        chr <- data$chr
    }

    if (is.null(chr) || length(unique(chr)) == 1) {
        PTA1(start, end, scores, ...)
    } else {
        result <- list(chr=c(),
                       start=c(),
                       end=c(),
                       scores=matrix(nrow=0, ncol=ncol(scores)),
                       coefficients=c(),
                       intercepts=c())
        for (chr1 in sort(unique(chr))) {
            pta <- PTA1(start[chr == chr1],
                        end[chr == chr1],
                        scores[chr == chr1, ],
                        ...)
            result$chr <- c(result$chr, rep(chr1, length(pta$start)))
            result$start <- c(result$start, pta$start)
            result$end <- c(result$end, pta$end)
            result$scores <- rbind(result$scores, pta$scores)
            result$coefficients <- c(result$coefficients, pta$coefficients)
            result$intercepts <- c(result$intercepts, pta$intercepts)
        }
        result$chr <- factor(result$chr)
        result
    }
}

PTA1 <- function(start, end, scores,
                 countBound=1, cumulativeErrorBound=Inf,
                 adjacencyThreshold=1, skip=0,
                 mode=c("normal", "correlation", "correlationSpearman"),
                 correlationBound=-1)
{
    mode <- match.arg(mode)
    mode <- switch(mode, normal=0, correlation=1, correlationSpearman=2)

    if (is.vector(scores)) {
        scores <- matrix(scores)
    }

    if (length(start) != length(end)) {
        stop("Start and end counts differ.")
    }

    if (length(start) != nrow(scores)) {
        stop("Range and score counts differ")
    }

    arguments <- as.list(environment())

    result <- .Call("PTA", arguments, PACKAGE="PTA")

    colnames(result$scores) <- colnames(scores)

    result$groups <- result$groups + 1

    result
}

#' @export
intersectionAggregate <- function(group, start, end, scores)
{
    if (length(group) != length(start)) { stop("Group and start counts differ.") }
    if (length(start) != length(end)) { stop("Start and end counts differ.") }
    if (is.vector(scores)) { scores <- as.matrix(scores) }
    if (length(start) != nrow(scores)) { stop("Range and score counts differ") }
    arguments <- as.list(environment())
    result <- .Call("intersectionAggregate", arguments, PACKAGE="PTA")
    colnames(result$scores) <- colnames(scores)
    result
}

#' @export
spanAggregate <- function(start=NULL, end=NULL, scores=NULL, chr=NULL, data=NULL, span) {
    if (class(data) == "GRanges") {
        data <- extractGRanges(data)
    }

    if (class(start) == "GRanges") {
        data <- extractGRanges(start)
    }

    if (!is.null(data)) {
        start <- data$start
        end <- data$end
        scores <- data$scores
        chr <- data$chr
    }

    if (is.null(chr) || length(unique(chr)) == 1) {
        r <- .Call("spanAggregate",
                   list(start=start, end=end, scores=scores, span=span),
                   PACKAGE="PTA")
        colnames(r$scores) <- colnames(scores)
        r
    } else {
        result <- list(chr=c(),
                       start=c(),
                       end=c(),
                       scores=matrix(nrow=0, ncol=ncol(scores)))
        colnames(result$scores) <- colnames(scores)
        for (chr1 in sort(unique(chr))) {
            r <- .Call("spanAggregate",
                       list(start=start[chr == chr1],
                            end=end[chr == chr1],
                            scores=as.matrix(scores[chr == chr1, ]),
                            span=span),
                       PACKAGE="PTA")
            result$chr <- c(result$chr, rep(chr1, length(r$start)))
            result$start <- c(result$start, r$start)
            result$end <- c(result$end, r$end)
            result$scores <- rbind(result$scores, r$scores)
        }
        result$chr <- factor(result$chr)
        result
    }
}

#' @export
applyPTAResult <- function(result,
                           start=NULL, end=NULL, scores=NULL, chr=NULL, data=NULL)
{
    if (class(data) == "GRanges") {
        data <- extractGRanges(data)
    }

    if (class(start) == "GRanges") {
        data <- extractGRanges(start)
    }

    if (!is.null(data)) {
        start <- data$start
        end <- data$end
        scores <- data$scores
        chr <- data$chr
    }

    if (is.null(chr) || length(unique(chr)) == 1) {
        applyPTAResult1(result, start, end, scores)
    } else {
        outscores <- matrix(nrow=nrow(result$scores), ncol=ncol(scores))
        for (chr1 in sort(unique(chr))) {
            subpta <- list(
                groups = result$groups[result$chr == chr1],
                coefficients = result$coefficients[result$chr == chr1],
                intercept = result$intercept[result$chr == chr1]
            )
            newscores <- applyPTAResult1(subpta,
                                         start[chr == chr1],
                                         end[chr == chr1],
                                         scores[chr == chr1, ])
            outscores[result$chr == chr1,] <- newscores
        }
        outscores
    }
}

applyPTAResult1 <- function(result, start, end, scores) {
    stopifnot(length(result$groups) == nrow(scores))
    len <- end - start
    s <- len * if (!is.null(result$coefficients)) {
        result$coefficients * scores + result$intercept
    } else {
        scores
    }

    filter <- result$groups != 0
    s <- s[filter,]
    groups <- result$groups[filter]
    len <- len[filter]

    grouplen <- tapply(len, groups, sum)

    x <- apply(as.matrix(s), 2,
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
    newscore <- if (is.matrix(srcScore)) {
        apply(srcScore, 2,
              function(scoreCol) {
                  .Call("adapt_ranged",
                        as.numeric(srcStart), as.numeric(srcEnd), as.numeric(scoreCol),
                        as.numeric(destStart), as.numeric(destEnd),
                        FALSE)
              })
    } else {
        .Call("adapt_ranged",
              as.numeric(srcStart), as.numeric(srcEnd), as.numeric(srcScore),
              as.numeric(destStart), as.numeric(destEnd),
              FALSE)
    }
    list(start=destStart,
         end=destEnd,
         score=newscore)
}
