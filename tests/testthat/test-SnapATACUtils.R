test_that("snapGmat2Seurat works", {
  s <- readRDS(test_path("testdata/", "testSnap.rds"))
  ss <- snapGmat2Seurat(snap = s, eigDims = 1:50,
                        assay = "gmat", pcaPrefix = "SnapATAC_",
                        useSnapATACEmbed = TRUE)
  expect_s4_class(ss, "Seurat")
})

test_that("plotConsensus works", {
  withr::with_tempfile(
    new = "tf",
    fileext = "pdf",
    code = {
      plotConsensus(test_path("testdata/",
                              "all_consensus_stat.txt"),
                  type = "AllNonN",
                  outPDF = tf)
      expect_true(file.exists(tf))
    }
  )
})
