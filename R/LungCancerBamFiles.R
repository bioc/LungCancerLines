LungCancerBamFiles <- function() {
  files <- dir(system.file("extdata", package = "LungCancerLines"),
               "\\.concordant_uniq\\.bam$", full.names = TRUE)
  names(files) <- sub("\\..*", "", basename(files))
  BamFileList(files)
}
