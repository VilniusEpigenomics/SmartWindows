context("spanAggregate")

library(GenomicRanges)

set.seed(123)
n <- 100
x <- rnorm(n)
y <- x + 0.5 * rnorm(n)
z <- 0.3 * x - 0.3 * y + 0.4 * rnorm(n)
dStart <- (1:n)*10
dEnd <- dStart + 6
scores <- cbind(x, y, z)
oneChr <- rep("chr1", 100)
twoChr <- rep(c("chr1", "chr2"), each=50)
gr <- GRanges(seqnames=oneChr, ranges=IRanges(start=dStart, end=dEnd), x=x, y=y, z=z)
gr2 <- GRanges(seqnames=twoChr, ranges=IRanges(start=dStart, end=dEnd), x=x, y=y, z=z)

colTApply <- function(mtx, index, fun) {
    apply(mtx, 2, function(col) tapply(col, index, fun))
}

test_that("correctly merges to one large span", {
    r <- spanAggregate(gr, span=1 + tail(dEnd, 1))
    expect_true(all(r$scores - colMeans(scores) < .Machine$double.eps))
    r <- spanAggregate(gr2, span=1 + tail(dEnd, 1))
    realScores <- colTApply(scores, rep(c(1, 2), each=50), mean)
    expect_true(all(r$scores - realScores < .Machine$double.eps))
})

test_that("correctly merges to ten spans", {
    r <- spanAggregate(gr, span=100)
    realScores <- apply(scores, 2, function(col) tapply(col, rep(1:10, each=10), mean))
    expect_true(all(r$scores - realScores < .Machine$double.eps))
    r <- spanAggregate(gr2, span=100)
    realScores <- apply(scores, 2, function(col) tapply(col, rep(1:10, each=10), mean))
    expect_true(all(r$scores - realScores < .Machine$double.eps))
})
