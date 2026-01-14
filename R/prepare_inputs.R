#' Prepare inputs: prevalence filter, mapping, sparse matrix
#'
#' @keywords internal
prepare_inputs <- function(raw_data, codebook,
                           prevalence_cutoff,
                           block_col = "block",
                           verbose = TRUE) {

  # 先用 base 方式保证取列不依赖 data.table NSE
  raw_data <- as.data.frame(raw_data)
  codebook <- as.data.frame(codebook)

  # 清理列名空格
  names(raw_data) <- trimws(names(raw_data))
  names(codebook) <- trimws(names(codebook))

  req_raw <- c("eid", "Category_No")
  req_cb  <- c("No", "Frequency", "Codes", "Names", block_col)
  miss_raw <- setdiff(req_raw, names(raw_data))
  miss_cb  <- setdiff(req_cb,  names(codebook))
  if (length(miss_raw) > 0) stop("raw_data missing columns: ", paste(miss_raw, collapse = ", "))
  if (length(miss_cb)  > 0) stop("codebook missing columns: ", paste(miss_cb,  collapse = ", "))

  # ✅ 取列：base 写法，绝对不会出现 .SD 问题
  dt <- raw_data[, c("eid", "Category_No"), drop = FALSE]
  dt <- data.table::as.data.table(dt)
  dt <- unique(dt)
  dt[, Category_No := as.character(Category_No)]

  cb <- data.table::as.data.table(codebook)
  cb[, No := as.character(No)]
  cb[, (block_col) := as.character(cb[[block_col]])]
  cb[is.na(cb[[block_col]]) | cb[[block_col]] == "", (block_col) := "Unknown"]

  total_patients <- data.table::uniqueN(dt$eid)
  min_freq_threshold <- ceiling(total_patients * prevalence_cutoff)

  if (verbose) {
    cat(sprintf("N=%d, prevalence_cutoff=%.4f, min_freq=%d\n",
                total_patients, prevalence_cutoff, min_freq_threshold))
  }

  valid_cb <- cb[cb[["Frequency"]] >= min_freq_threshold]
  if (nrow(valid_cb) == 0) stop("No nodes remain after prevalence filtering.")

  valid_diseases <- valid_cb$No
  current_data <- dt[Category_No %in% valid_diseases]
  current_data <- unique(current_data)

  fixed_nodes <- sort(unique(current_data$Category_No))
  n_nodes <- length(fixed_nodes)
  node_map_fixed <- data.table::data.table(Category_No = fixed_nodes, col_idx = 1:n_nodes)

  unique_eids <- sort(unique(current_data$eid))
  n_patients <- length(unique_eids)
  eid_map <- data.table::data.table(eid = unique_eids, row_idx = 1:n_patients)

  tmp_blk <- valid_cb[valid_cb$No %in% fixed_nodes, c("No", block_col), with = FALSE]
  data.table::setnames(tmp_blk, c("No", block_col), c("Category_No", "Block"))
  block_map_dt <- unique(tmp_blk, by = "Category_No")
  block_map_dt <- merge(node_map_fixed, block_map_dt, by = "Category_No", all.x = TRUE)
  block_map_dt[is.na(Block) | Block == "", Block := "Unknown"]

  tmp_full <- merge(current_data, eid_map, by = "eid")
  tmp_full <- merge(tmp_full, node_map_fixed, by = "Category_No")

  M_full <- Matrix::sparseMatrix(
    i = tmp_full$row_idx,
    j = tmp_full$col_idx,
    x = 1,
    dims = c(n_patients, n_nodes)
  )

  lineage_dt <- data.table::data.table(
    Current_ID = as.character(valid_cb$No),
    Original_Codes = valid_cb$Codes,
    Original_Names = valid_cb$Names,
    Block = valid_cb[[block_col]],
    Member_Count = 1
  )
  lineage_dt[is.na(Block) | Block == "", Block := "Unknown"]
  lineage_dt <- unique(lineage_dt, by = "Current_ID")

  list(
    current_data = current_data,
    fixed_nodes = fixed_nodes,
    block_map_dt = block_map_dt,
    lineage_dt = lineage_dt,
    M_full = M_full
  )
}
