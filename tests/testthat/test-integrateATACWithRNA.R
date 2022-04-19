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
