context("PTA")

set.seed(123)
x <- rnorm(20)
y <- x + 0.5 * rnorm(20)
z <- 0.3 * x - 0.3 * y + 0.4 * rnorm(20)
d <- RangedData(ranges=IRanges(start=(1:20)*10, width=6), x=x, y=y, z=z)

test_that("normal mode with error.bound works", {
    p <- PTA(d, adjacency.threshold=10, error.bound=0.1)
    expect_true(nrow(p) < nrow(d))
})

test_that("skip works", {
    p <- PTA(d, skip=1, adjacency.threshold=20, error.bound=0.1)
    expect_true(nrow(p) < nrow(d))
})

test_that("normal mode with count.bound works", {
    p <- PTA(d, adjacency.threshold=10, count.bound=10)
    expect_true(nrow(p) == 10)
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
    p <- PTA.raw(1:10, 1:10+1, 1:10)
    expect_true(nrow(p$scores) == length(p$start))
    expect_true(nrow(p$scores) == length(p$end))
    expect_true(ncol(p$scores) == 1)
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
