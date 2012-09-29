LungCancerBamFiles <- function() {
  files <- dir(system.file("extdata", package = "LungCancerLines"),
               "\\.bam$", full.names = TRUE)
  names(files) <- sub("\\..*", "", basename(files))
  BamFileList(files)
}
