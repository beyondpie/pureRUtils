test_that("sampleBasedOnDepth2 works", {
  log10UMI <- rnorm(100, mean = 1000, sd = 10)
  expect_equal(50, length(sampleBasedOnDepth2(log10UMI, n = 50)))
})
