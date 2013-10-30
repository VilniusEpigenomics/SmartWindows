context("PTA")

library(GenomicRanges)

set.seed(123)
x <- rnorm(20)
y <- x + 0.5 * rnorm(20)
z <- 0.3 * x - 0.3 * y + 0.4 * rnorm(20)
d.start <- (1:20)*10
d.end <- d.start + 6
scores <- cbind(x, y, z)
d <- RangedData(ranges=IRanges(start=d.start, end=d.end), x=x, y=y, z=z)
d.gr <- GRanges(seqnames="chr0", ranges=IRanges(start=d.start, end=d.end), x=x, y=y, z=z)

test_that("normal mode with error bounds works", {
    p <- PTA(d, adjacency.threshold=10, cumulative.error.bound=0.1)
    expect_true(nrow(p) < nrow(d))
    p <- PTA(d.gr, adjacency.threshold=10, cumulative.error.bound=0.1)
    expect_true(length(p) < length(d.gr))

    p <- PTA(d, adjacency.threshold=10, error.bound=100)
    expect_true(nrow(p) < nrow(d))
    p <- PTA(d.gr, adjacency.threshold=10, error.bound=100)
    expect_true(length(p) < length(d.gr))
})

test_that("skip works", {
    p <- PTA(d, skip=1, adjacency.threshold=20, cumulative.error.bound=0.1)
    expect_true(nrow(p) < nrow(d))
})

test_that("normal mode with count.bound works", {
    p <- PTA(d, adjacency.threshold=10, count.bound=10)
    expect_true(nrow(p) == 10)
    p <- PTA(d.gr, adjacency.threshold=10, count.bound=10)
    expect_true(length(p) == 10)
})

test_that("correlation mode works", {
    p <- PTA(d, adjacency.threshold=10, mode="correlation", correlation.bound=0.8)
    expect_true(nrow(p) < nrow(d))
})

test_that("spearman correlation works", {
    p <- PTA(d, adjacency.threshold=10, mode="correlation", correlation.bound=0.8, correlation.spearman=TRUE)
    expect_true(nrow(p) < nrow(d))
})

test_that("PTA.raw works", {
    p <- PTA.raw(d.start, d.end, scores, adjacency.threshold=10, count.bound=10)
    expect_true(10 == nrow(p$scores))
    expect_true(10 == length(p$start))
    expect_true(10 == length(p$end))
    expect_true(ncol(p$scores) == 3)
    expect_true(p$cumulative.error > 0)
})

test_that("PTA.raw doesn't change arguments", {
    start <- as.integer(1:10)
    end <- as.integer(start + 1)
    score <- as.double(1:10)
    p <- PTA.raw(start, end, score, adjacency.threshold=2)
    expect_equal(start, 1:10)
    expect_equal(end, 1:10 + 1)
    expect_equal(score, 1:10) 
})

test_that("apply.PTA.result works", {
    scores2 <- 2 * scores + rnorm(20)
    p <- PTA.raw(d.start, d.end, scores, adjacency.threshold=10, mode="correlation", count.bound=10)
    p2.scores <- apply.PTA.result(p, d.start, d.end, scores2)
    expect_true(all(dim(p$scores) == dim(p2.scores)))
})

test_that("apply.PTA.result merges same data with small difference", {
    p <- PTA.raw(d.start, d.end, scores, mode="correlation", count.bound=10)
    p2.scores <- apply.PTA.result(p, d.start, d.end, scores)
    err <- abs(p$scores - p2.scores)
    expect_true(sum(err) / sum(scores) < 1e-16)

    p <- PTA.raw(d.start, d.end, scores, mode="correlation", count.bound=10, skip=1)
    p2.scores <- apply.PTA.result(p, d.start, d.end, scores)
    err <- abs(p$scores - p2.scores)
    expect_true(sum(err) / sum(scores) < 1e-16)
})
