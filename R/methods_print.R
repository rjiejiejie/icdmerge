#' @export
print.icdmerge_result <- function(x, ...) {
  cat("icdmerge_result\n")
  cat(sprintf("  Nodes (post-merge): %d\n", data.table::uniqueN(x$merged_data$Category_No)))
  cat(sprintf("  Merge groups: %d\n", sum(x$lineage$Member_Count > 1)))
  cat(sprintf("  Runtime (sec): %.1f\n", x$runtime_sec))
  invisible(x)
}
