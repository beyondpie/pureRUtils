test_that("loadRData works", {
  expect_equal(loadRData(fileRData = system.file("test/a.RData", package = "pureRUtils")), c(1,2,3))
})
