context("PTA")

library(GenomicRanges)

set.seed(123)
n <- 100
x <- rnorm(n)
y <- x + 0.5 * rnorm(n)
z <- 0.3 * x - 0.3 * y + 0.4 * rnorm(n)
dStart <- (1:n)*10
dEnd <- dStart + 6
scores <- cbind(x, y, z)
d <- GRanges(seqnames="chr0", ranges=IRanges(start=dStart, end=dEnd), x=x, y=y, z=z)

test_that("normal mode with errorBound works", {
    p <- PTA(d, adjacencyThreshold=10, cumulativeErrorBound=0.1)
    expect_true(nrow(p$scores) < length(d))

    p <- PTA(d, adjacencyThreshold=10, errorBound=100)
    expect_true(nrow(p$scores) < length(d))
})

test_that("skip works", {
    p <- PTA(d, skip=1, adjacencyThreshold=20, cumulativeErrorBound=0.1)
    expect_true(nrow(p$scores) < length(d))
})

test_that("normal mode with countBound works", {
    p <- PTA(d, adjacencyThreshold=10, countBound=10)
    expect_equal(nrow(p$scores), 10)
})

test_that("correlation mode works", {
    p <- PTA(d, adjacencyThreshold=10, mode="correlation", correlationBound=0.8)
    expect_true(nrow(p$scores) < length(d))
})

test_that("spearman correlation works", {
    p <- PTA(d, adjacencyThreshold=10, mode="correlation", correlationBound=0.8, correlationSpearman=TRUE)
    expect_true(nrow(p$scores) < length(d))
})


test_that("PTA works on raw data", {
    p <- PTA(dStart, dEnd, scores, adjacencyThreshold=10, countBound=10)
    expect_true(10 == nrow(p$scores))
    expect_true(10 == length(p$start))
    expect_true(10 == length(p$end))
    expect_true(ncol(p$scores) == 3)
    expect_true(p$cumulativeError > 0)
})

test_that("PTA doesn't change arguments", {
    start <- as.integer(1:10)
    end <- as.integer(start + 1)
    score <- as.double(1:10)
    p <- PTA(start, end, score, adjacencyThreshold=2)
    expect_equal(start, 1:10)
    expect_equal(end, 1:10 + 1)
    expect_equal(score, 1:10) 
})

test_that("skipped nodes have group 0", {
    p <- PTA(dStart, dEnd, scores, adjacencyThreshold=20, skip=1, mode="correlation", correlationBound=0.8)
    expect_equal(length(unique(p$groups)), 1 + max(p$groups))
    expect_true(0 %in% p$groups)

    diffSkipped <- diff(p$groups == 0)
    for (i in 1:length(diffSkipped - 1)) {
        if (diffSkipped[i] == 1) {
            expect_equal(diffSkipped[i + 1], -1)
        }
    }
})

test_that("NaN and infinite values return an error", {
    scoresInf <- scores
    scoresInf[1, 1] <- 1.0/0
    expect_that(PTA(dStart, dEnd, scoresInf), throws_error())
    scoresInf[1, 1] <- -1.0/0
    expect_that(PTA(dStart, dEnd, scoresInf), throws_error())
    scoresInf[1, 1] <- 0.0/0
    expect_that(PTA(dStart, dEnd, scoresInf), throws_error())
})


test_that("applyPTAResult works", {
    scores2 <- 2 * scores + rnorm(n)
    p <- PTA(dStart, dEnd, scores, adjacencyThreshold=10, mode="correlation", countBound=10)
    p2Scores <- applyPTAResult(p, dStart, dEnd, scores2)
    expect_true(all(dim(p$scores) == dim(p2Scores)))
})

test_that("applyPTAResult merges same data with small difference", {
    p <- PTA(dStart, dEnd, scores, mode="correlation", countBound=10)
    p2Scores <- applyPTAResult(p, dStart, dEnd, scores)
    err <- abs(p$scores - p2Scores)
    expect_true(sum(err) / sum(scores) < 1e-15)

    p <- PTA(dStart, dEnd, scores, mode="correlation", countBound=10, skip=1)
    p2Scores <- applyPTAResult(p, dStart, dEnd, scores)
    err <- abs(p$scores - p2Scores)
    expect_true(sum(err) / sum(scores) < 1e-15)
})
