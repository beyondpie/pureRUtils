test_that("integrateWithScRNASeq works", {
  snapSeurat <- readRDS(test_path("testdata", "testSnapSeurat.rds"))
  mbSeurat <- readRDS(test_path("testdata", "testMouseBrainOrgSeurat.rds"))
  coEmbed <- integrateWithScRNASeq(snapSeurat = snapSeurat,
                                   rnaSeurat = mbSeurat,
                                   eigDims = 1:20,
                                   snapAssay = "gmat",
                                   preprocessSnap = T,
                                   preprocessRNA = T,
                                   reso = 0.5)
  expect_equal(ncol(coEmbed), ncol(snapSeurat) + ncol(mbSeurat))
})

test_that("getOverlapMatrix works", {
  meta <- readRDS(test_path("testdata", "testOvlp.rds"))
  score <- getOverlapMatrix(meta = meta,
                            ident = "idents",
                            atacCol = "atacLabel",
                            rnaCol = "rnaLabel")
  expect_equal(dim(score),
               c(length(unique(meta$atacLabel)) - 1,
                 length(unique(meta$rnaLabel)) - 1))
})

test_that("getOverlapMatrixOld works", {
  meta <- readRDS(test_path("testdata", "testOvlp.rds"))
  score <- getOverlapMatrixOld(meta = meta,
                            ident = "idents",
                            atacCol = "atacLabel",
                            rnaCol = "rnaLabel")
  expect_equal(dim(score),
               c(length(unique(meta$atacLabel)) - 1,
                 length(unique(meta$rnaLabel)) - 1))
})

