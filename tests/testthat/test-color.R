test_that("getMultipleColors works", {
  expect_vector(
    getMultipleColors(1:10), ptype = character(), size = 10
  )
})
