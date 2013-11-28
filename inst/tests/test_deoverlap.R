context("deoverlap")

test_that("deoverlaps simple data", {
    initial <- IRanges(start=c(1, 3, 10), end=c(5, 12, 15))
    expect_equal(deoverlapRanges(initial),
                 IRanges(start=c(1, 4, 11), end=c(4, 11, 15)))
})
