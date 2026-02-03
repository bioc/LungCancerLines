LungCancerBamFiles <- function() {
  files <- dir(system.file("extdata", package = "LungCancerLines"),
               "\\.analyzed\\.bam$", full.names = TRUE)
  names(files) <- sub("\\..*", "", basename(files))
  BamFileList(files)
}
