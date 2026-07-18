#!/usr/bin/env Rscript
# Fast exposure preparation entry point for restored IS_Analysis project.
args <- commandArgs(trailingOnly = TRUE)
message("Preparing exposure data with fast workflow")
if (length(args) > 0) {
  message("Arguments: ", paste(args, collapse = " "))
}
