
test_that("grepPatternList works",{
  patternList <- c("R1_001.fastq.gz", "R1.fastq.gz", "R1.fastq",
                   "F3.csfasta")
  x <- c("nanocage_ACAGAT_carp_embryo_R1.fastq.gz",
         "nanocage_CACGAT_carp_hk_control_R1.fastq.gz",
         "nanocage_CACTGA_carp_hk_Tborreli1_R1.fastq.gz",
         "nanocage_CTGACG_carp_hk_Tborreli2_R1.fastq.gz")
  expect_true(all(grepPatternList(patternList, x, combine="or")))
  expect_false(any(grepPatternList(patternList, x, combine="and")))
  })