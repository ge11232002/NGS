library(NGS)

test_that("getPairedReads returns the right filename",{
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                      "nanocage_ACAGAT_carp_embryo_R1.fastq")
  expect_equal(basename(getPairedReads(reads1)), 
               "nanocage_ACAGAT_carp_embryo_R2.fastq")
})

test_that("trimFastq generated the right fastq", {
  reads1 <- file.path(system.file("extdata", package="NGS"), "fastq",
                      "nanocage_ACAGAT_carp_embryo_R1.fastq")
  ## minTailQuality = 20
  trimmedFastqs <- trimFastq(reads1, outputDir=".", paired=TRUE, nReads=-1L,
                             minTailQuality=20, trimLeft=0L, trimRight=0L,
                             mc.cores=1L)
  library(ShortRead)
  library(IRanges)
  trimmedReads <- readFastq(trimmedFastqs)
  expect_equal(as.matrix(table(width(trimmedReads))),
               matrix(c(2,1,1,1,1,23,21), 
                      dimnames=list(c("105", "107", "109", "119", "124", 
                                      "128", "129")))
  )
  
  file.remove(trimmedFastqs, NGS:::getPairedReads(trimmedFastqs))
})
