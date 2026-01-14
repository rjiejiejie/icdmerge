#' Block-wise disease node merging with PAC-selected consensus kmeans
#'
#' @description
#' Build global structural profiles for diseases, perform block-wise consensus k-means clustering
#' under bootstrap resampling, select K using PAC, and optionally allow full-block merge (K=1)
#' when K=2 consensus is nearly uniform.
#'
#' @param raw_data data.table/data.frame with columns eid, Category_No
#' @param codebook data.table/data.frame with columns No, Frequency, Codes, Names, block
#' @param prevalence_cutoff numeric, e.g. 0.001 for 0.1%
#' @param n_resamples integer number of resamples
#' @param sample_frac numeric fraction of patients per resample
#' @param pac_x1,pac_x2 numeric PAC interval bounds
#' @param pac_no_merge numeric if min PAC among K>=2 > this, set K*=1
#' @param homo_q,homo_tau numeric homogeneity gate for allowing full-block merge (K=1)
#' @param k_max_per_block integer maximum K considered per block
#' @param kmeans_nstart_resample,kmeans_nstart_final,kmeans_itmax integers
#' @param min_valid_frac numeric minimal valid fraction for each (block,K)
#' @param seed integer random seed
#' @param write_files logical whether to write CSV outputs
#' @param output_dir character output directory when write_files=TRUE
#' @param verbose logical
#'
#' @return An object of class \code{icdmerge_result}
#' @export
run_block_merge <- function(raw_data, codebook,
                            prevalence_cutoff = 0.001,
                            n_resamples = 1000,
                            sample_frac = 0.8,
                            pac_x1 = 0.1,
                            pac_x2 = 0.9,
                            pac_no_merge = 0.10,
                            homo_q = 0.10,
                            homo_tau = 0.90,
                            k_max_per_block = 10,
                            kmeans_nstart_resample = 5,
                            kmeans_nstart_final = 50,
                            kmeans_itmax = 100,
                            min_valid_frac = 0.5,
                            seed = 123,
                            write_files = FALSE,
                            output_dir = NULL,
                            verbose = TRUE) {

  raw_data <- data.table::as.data.table(raw_data)
  codebook <- data.table::as.data.table(codebook)

  req_raw <- c("eid", "Category_No")
  req_cb  <- c("No", "Frequency", "Codes", "Names", "block")
  missing_raw <- setdiff(req_raw, names(raw_data))
  missing_cb  <- setdiff(req_cb,  names(codebook))
  if (length(missing_raw) > 0) stop("raw_data missing columns: ", paste(missing_raw, collapse=", "))
  if (length(missing_cb)  > 0) stop("codebook missing columns: ", paste(missing_cb,  collapse=", "))

  if (write_files) {
    if (is.null(output_dir) || !nzchar(output_dir)) stop("write_files=TRUE requires output_dir.")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }

  set.seed(seed)
  t0 <- Sys.time()

  if (verbose) {
    cat("========== icdmerge run_block_merge ==========\n")
    cat(sprintf("Start time: %s\n", format(t0)))
    cat(sprintf("prevalence_cutoff=%.4f | n_resamples=%d | sample_frac=%.2f\n",
                prevalence_cutoff, n_resamples, sample_frac))
    cat(sprintf("PAC: x1=%.2f, x2=%.2f, no_merge=%.2f | HOMO: q=%.2f, tau=%.2f\n",
                pac_x1, pac_x2, pac_no_merge, homo_q, homo_tau))
  }

  # -------------------
  # 1) prepare
  # -------------------
  if (verbose) cat("\n[1/4] Preparing inputs...\n")
  prep <- prepare_inputs(
    raw_data = raw_data,
    codebook = codebook,
    prevalence_cutoff = prevalence_cutoff,
    block_col = "block",
    verbose = verbose
  )

  if (verbose) {
    cat(sprintf("Prepared: patients=%d | diseases=%d | blocks=%d\n",
                nrow(prep$M_full), ncol(prep$M_full),
                length(setdiff(unique(prep$block_map_dt$Block), "Unknown"))))
  }

  # -------------------
  # 2) consensus + PAC
  # -------------------
  if (verbose) cat("\n[2/4] Running consensus resampling + PAC selection...\n")
  cons <- consensus_pac_select_k(
    M_full = prep$M_full,
    fixed_nodes = prep$fixed_nodes,
    block_map_dt = prep$block_map_dt,
    lineage_dt = prep$lineage_dt,
    n_resamples = n_resamples,
    sample_frac = sample_frac,
    pac_x1 = pac_x1,
    pac_x2 = pac_x2,
    pac_no_merge = pac_no_merge,
    homo_q = homo_q,
    homo_tau = homo_tau,
    k_max_per_block = k_max_per_block,
    kmeans_nstart_resample = kmeans_nstart_resample,
    kmeans_nstart_final = kmeans_nstart_final,
    kmeans_itmax = kmeans_itmax,
    min_valid_frac = min_valid_frac,
    seed = seed,
    verbose = verbose
  )

  if (verbose) {
    kdist <- cons$kstar[, .N, by = K_star][order(K_star)]
    cat("\nK* distribution:\n")
    print(kdist)
    if ("MergeAll" %in% names(cons$kstar)) {
      cat(sprintf("MergeAll triggered blocks: %d\n", sum(cons$kstar$MergeAll %in% TRUE)))
    }
  }

  # -------------------
  # 3) merge
  # -------------------
  if (verbose) cat("\n[3/4] Applying merges...\n")
  merged <- apply_merges(
    current_data = prep$current_data,
    lineage_dt = prep$lineage_dt,
    kstar_all = cons$kstar,
    final_membership_list = cons$final_membership_list,
    verbose = verbose
  )

  # -------------------
  # 4) finalize + write
  # -------------------
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  res <- list(
    params = list(
      prevalence_cutoff = prevalence_cutoff,
      n_resamples = n_resamples,
      sample_frac = sample_frac,
      pac_x1 = pac_x1, pac_x2 = pac_x2,
      pac_no_merge = pac_no_merge,
      homo_q = homo_q, homo_tau = homo_tau,
      k_max_per_block = k_max_per_block,
      kmeans_nstart_resample = kmeans_nstart_resample,
      kmeans_nstart_final = kmeans_nstart_final,
      kmeans_itmax = kmeans_itmax,
      min_valid_frac = min_valid_frac,
      seed = seed
    ),
    pac_curves = cons$pac_curves,
    kstar = cons$kstar,
    membership = cons$membership,
    update_map = merged$update_map,
    merged_data = merged$current_data,
    lineage = merged$lineage_dt,
    runtime_sec = runtime
  )
  class(res) <- "icdmerge_result"

  if (write_files) {
    if (verbose) cat("\n[4/4] Writing output files...\n")
    data.table::fwrite(res$pac_curves,  file.path(output_dir, "PAC_Curves_ByBlock.csv"))
    data.table::fwrite(res$kstar,       file.path(output_dir, "PAC_Summary_Kstar.csv"))
    data.table::fwrite(res$membership,  file.path(output_dir, "Final_Block_Memberships.csv"))
    data.table::fwrite(res$update_map,  file.path(output_dir, "Merge_Update_Map.csv"))
    data.table::fwrite(res$merged_data, file.path(output_dir, "Final_Merged_Data.csv"))
    data.table::fwrite(res$lineage,     file.path(output_dir, "Final_Merged_Dict.csv"))
    if (verbose) cat(sprintf("Files written to: %s\n", normalizePath(output_dir, winslash = "/")))
  }

  if (verbose) {
    cat("\n========== icdmerge finished ==========\n")
    cat(sprintf("Runtime: %.1f sec\n", runtime))
    cat(sprintf("Post-merge nodes: %d | merge groups: %d\n",
                data.table::uniqueN(res$merged_data$Category_No),
                sum(res$lineage$Member_Count > 1)))
  }

  res
}

