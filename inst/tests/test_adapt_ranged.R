context("adapt.ranged")

empty <- RangedData(ranges=IRanges(), score=numeric())
d1 <- RangedData(ranges=IRanges(start=c(1, 3, 6), end=c(3, 6, 9)), score=1:3)
d2 <- RangedData(ranges=IRanges(start=c(1, 2, 7, 8), end=c(2, 7, 8, 10)), score=1:4)
d3 <- RangedData(ranges=IRanges(start=c(1, 7, 8), end=c(7, 8, 10)), score=1:3)
d4 <- RangedData(ranges=IRanges(start=c(1, 7, 8), end=c(4, 8, 10)), score=1:3)

test_that("empty ranges work", {
    expect_equal(adapt.ranged(d2, empty), empty)
    expect_equal(adapt.ranged(empty, d2, na.rm=TRUE), empty)

    expect_equal(length(adapt.ranged.raw(start(d2), end(d2), score(d2), start(empty), end(empty))$score), 0)
    expect_equal(adapt.ranged.raw(integer(), integer(), integer(), start(d2), end(d2))$score, rep(NA, 4))
})

test_that("doesn't change identical intervals", {
    expect_equal(adapt.ranged(d1, d1), d1)
    expect_equal(adapt.ranged(d2, d2), d2)

    expect_equal(adapt.ranged.raw(start(d1), end(d1), score(d1), start(d1), end(d1)),
                 list(start=start(d1), end=end(d1), score=score(d1)))
    expect_equal(adapt.ranged.raw(start(d2), end(d2), score(d2), start(d2), end(d2)),
                 list(start=start(d2), end=end(d2), score=score(d2)))
})

test_that("multiple src intervals per dest work", {
    expect_equal(adapt.ranged(d2, d3),
                 RangedData(ranges=ranges(d3), score=c((1+2*5)/6, 3, 4)))

    expect_equal(adapt.ranged.raw(start(d2), end(d2), score(d2), start(d3), end(d3)),
                 list(start=start(d3), end=end(d3), score=c((1+2*5)/6, 3, 4)))
})

test_that("merges intersecting intervals", {
    expect_equal(adapt.ranged(d2, d1),
                 RangedData(ranges=ranges(d1), score=c(1.5, 2, 3)))

    expect_equal(adapt.ranged.raw(start(d2), end(d2), score(d2), start(d1), end(d1)),
                 list(start=start(d1), end=end(d1), score=c(1.5, 2, 3)))
})

test_that("works with spaces between intervals", {
    expect_equal(adapt.ranged(d2, d4),
                 RangedData(ranges=ranges(d4), score=c((1+2*2)/3, 3, 4)))

    expect_equal(adapt.ranged.raw(start(d2), end(d2), score(d2), start(d4), end(d4)),
                 list(start=start(d4), end=end(d4), score=c((1+2*2)/3, 3, 4)))
})

test_that("error is calculated correctly", {
    expect_equal(adapt.ranged(d1, d1, add.error=TRUE)$error, c(0, 0, 0))
    expect_equal(adapt.ranged(d3, d4, add.error=TRUE)$error, c(0, 0, 0))
    expect_true(sum(adapt.ranged(d1, d4, add.error=TRUE)$error) > 0)
})
