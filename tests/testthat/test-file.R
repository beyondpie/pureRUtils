test_that("loadRData works", {
  expect_equal(loadRData(fileRData = test_path("testdata", "a.RData")), c(1,2,3))
})
