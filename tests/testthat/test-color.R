test_that("getMultipleColors works", {
  expect_vector(
    getMultipleColors(1:10), ptype = character(), size = 10
  )
})

test_that("getContinuousColors works", {
  expect_equal(
    nchar(getContinuousColors()(0.1)), 9
  )
})
