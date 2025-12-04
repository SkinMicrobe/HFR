# Temp directory helper (ASCII only).
suppressPackageStartupMessages({
  library(utils)
})

# Force tempdir() to a fixed path (anonymize real paths as needed).
set_tempdir_fixed <- function(path = "/path/to/temp") {
  assignInNamespace("tempdir", function() path, ns = "base")
  invisible(path)
}
