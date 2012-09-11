LungCancerFastqFiles <- function() {
  files <- dir(system.file("extdata", package = "LungCancerLines"),
               "\\.fastq\\.gz$", full.names = TRUE)
  lines <- sub("\\..*", "", basename(files))
  ends <- sub(".*-([^.]*).*", "\\1", basename(files))
  names(files) <- paste(lines, ends, sep = ".")
  files
}
