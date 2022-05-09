test_that("snapGmat2Seurat works", {
  s <- readRDS(test_path("testdata/", "testSnap.rds"))
  ss <- snapGmat2Seurat(snap = s, eigDims = 1:50,
                        assay = "gmat", pcaPrefix = "SnapATAC_",
                        useSnapATACEmbed = TRUE)
  expect_s4_class(ss, "Seurat")
})
