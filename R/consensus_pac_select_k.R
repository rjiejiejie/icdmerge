#' Consensus clustering + PAC selection within blocks
#'
#' @description
#' For each resample, compute global structural profiles for diseases (block rows × all diseases columns),
#' run kmeans for each candidate K within block, accumulate co-clustering counts, compute consensus probability,
#' then select K by PAC. Additionally, allow full-block merge (K=1) when K=2 consensus is nearly uniform.
#'
#' @param M_full Sparse patient × disease matrix.
#' @param fixed_nodes Character vector of disease ids in column order of M_full.
#' @param block_map_dt data.table with columns Category_No, col_idx, Block.
#' @param lineage_dt data.table used only for consistency (not required in current computation).
#' @param n_resamples integer number of resamples.
#' @param sample_frac numeric sample fraction of patients per resample.
#' @param pac_x1,pac_x2 numeric PAC interval bounds.
#' @param pac_no_merge numeric: if min PAC among K>=2 is above this, set K*=1.
#' @param homo_q,homo_tau numeric: full-block merge gate based on K=2 consensus quantile.
#' @param k_max_per_block integer maximum K considered per block.
#' @param kmeans_nstart_resample,kmeans_nstart_final,kmeans_itmax kmeans controls.
#' @param min_valid_frac minimal fraction of resamples that must be valid for a given (block,K).
#' @param seed integer.
#' @param verbose logical.
#'
#' @return list with pac_curves, kstar, membership, final_membership_list
#' @keywords internal
consensus_pac_select_k <- function(M_full,
                                   fixed_nodes,
                                   block_map_dt,
                                   lineage_dt,
                                   n_resamples,
                                   sample_frac,
                                   pac_x1, pac_x2,
                                   pac_no_merge,
                                   homo_q, homo_tau,
                                   k_max_per_block,
                                   kmeans_nstart_resample,
                                   kmeans_nstart_final,
                                   kmeans_itmax,
                                   min_valid_frac,
                                   seed = 123,
                                   verbose = TRUE) {

  set.seed(seed)

  # blocks to process
  block_map_dt <- data.table::as.data.table(block_map_dt)
  blocks <- sort(unique(block_map_dt$Block))
  blocks <- blocks[blocks != "Unknown"]
  n_blocks <- length(blocks)

  if (verbose) {
    cat(sprintf("Consensus/PAC: processing %d blocks\n", n_blocks))
  }

  # precompute per-block indices and candidate K list
  blk_cols_list <- stats::setNames(vector("list", n_blocks), blocks)
  blk_ids_list  <- stats::setNames(vector("list", n_blocks), blocks)
  Kcand_list    <- stats::setNames(vector("list", n_blocks), blocks)
  nblk_list     <- stats::setNames(integer(n_blocks), blocks)

  for (blk in blocks) {
    dt <- block_map_dt[Block == blk, list(Category_No, col_idx)]
    blk_cols_list[[blk]] <- dt$col_idx
    blk_ids_list[[blk]]  <- dt$Category_No
    nblk <- length(dt$col_idx)
    nblk_list[[blk]] <- nblk

    if (nblk < 3) {
      Kcand_list[[blk]] <- integer(0)
    } else {
      kmax <- min(k_max_per_block, nblk - 1)
      Kcand_list[[blk]] <- 2:kmax
    }
  }

  # containers: consensus counts and valid counts
  consensus_counts <- stats::setNames(vector("list", n_blocks), blocks)
  valid_counts     <- stats::setNames(vector("list", n_blocks), blocks)

  for (blk in blocks) {
    nblk <- nblk_list[[blk]]
    Ks <- Kcand_list[[blk]]
    if (length(Ks) == 0) {
      consensus_counts[[blk]] <- NULL
      valid_counts[[blk]] <- NULL
    } else {
      consensus_counts[[blk]] <- lapply(Ks, function(k) matrix(0, nblk, nblk))
      names(consensus_counts[[blk]]) <- as.character(Ks)
      valid_counts[[blk]] <- stats::setNames(as.list(rep(0L, length(Ks))), as.character(Ks))
    }
  }

  n_patients <- nrow(M_full)
  min_valid_n <- ceiling(n_resamples * min_valid_frac)

  if (verbose) {
    cat(sprintf("Resamples=%d, sample_frac=%.2f, min_valid_n=%d\n", n_resamples, sample_frac, min_valid_n))
  }

  # main loop: resample -> block -> K
  pb <- utils::txtProgressBar(min = 0, max = n_resamples, style = 3)

  for (t in 1:n_resamples) {
    sampled_rows <- sample.int(n_patients, size = floor(n_patients * sample_frac), replace = FALSE)
    M_sub <- M_full[sampled_rows, , drop = FALSE]

    for (blk in blocks) {
      Ks <- Kcand_list[[blk]]
      if (length(Ks) == 0) next

      blk_cols <- blk_cols_list[[blk]]

      # block × all-diseases structural profile
      X_blk <- Matrix::crossprod(M_sub[, blk_cols, drop = FALSE], M_sub)
      X_blk <- as.matrix(X_blk)
      X_blk <- row_zscore(X_blk)

      # distinct row check to avoid kmeans errors
      n_distinct <- nrow(unique(X_blk))

      for (k in Ks) {
        if (n_distinct < k) next

        km <- stats::kmeans(X_blk, centers = k, nstart = kmeans_nstart_resample, iter.max = kmeans_itmax)
        coc <- cocluster_mat(km$cluster, k)

        consensus_counts[[blk]][[as.character(k)]] <- consensus_counts[[blk]][[as.character(k)]] + coc
        valid_counts[[blk]][[as.character(k)]] <- valid_counts[[blk]][[as.character(k)]] + 1L
      }
    }

    utils::setTxtProgressBar(pb, t)
  }
  close(pb)

  if (verbose) cat("\nResampling done. Selecting K* by PAC + homogeneity gate...\n")

  # outputs
  pac_curves <- list()
  kstar_summary <- list()
  final_membership_list <- list()

  for (blk in blocks) {
    nblk <- nblk_list[[blk]]
    ids <- blk_ids_list[[blk]]
    Ks <- Kcand_list[[blk]]

    merge_all <- FALSE
    homo_q2 <- NA_real_
    valid_n2 <- NA_integer_

    # small blocks: default K=1, no full-block merge (can be extended later)
    if (nblk < 3 || length(Ks) == 0) {
      final_membership_list[[blk]] <- stats::setNames(rep(1L, nblk), ids)
      pac_curves[[blk]] <- data.table::data.table(
        Block = blk, K = NA_integer_, PAC = NA_real_, valid_n = NA_integer_,
        n_nodes = nblk, K_star = 1L
      )
      kstar_summary[[blk]] <- data.table::data.table(
        Block = blk, n_nodes = nblk, K_star = 1L,
        PAC_min = NA_real_, MergeAll = FALSE, HomoQ2 = homo_q2, ValidN2 = valid_n2
      )
      next
    }

    # compute PAC curve for candidate Ks
    pac_dt <- data.table::rbindlist(lapply(Ks, function(k) {
      vn <- valid_counts[[blk]][[as.character(k)]]
      if (is.null(vn) || vn < min_valid_n) {
        data.table::data.table(K = k, PAC = NA_real_, valid_n = vn)
      } else {
        cnt <- consensus_counts[[blk]][[as.character(k)]]
        prob <- cnt / vn
        data.table::data.table(K = k, PAC = calc_pac(prob, pac_x1, pac_x2), valid_n = vn)
      }
    }))

    pac_dt[, Block := blk]
    pac_dt[, n_nodes := nblk]

    # homogeneity gate: based on K=2 consensus quantile
    if ("2" %in% names(consensus_counts[[blk]])) {
      valid_n2 <- valid_counts[[blk]][["2"]]
      if (!is.null(valid_n2) && valid_n2 >= min_valid_n) {
        prob2 <- consensus_counts[[blk]][["2"]] / valid_n2
        homo_q2 <- calc_homo_q(prob2, q = homo_q)
        if (!is.na(homo_q2) && homo_q2 >= homo_tau) {
          merge_all <- TRUE
        }
      }
    }

    # choose K*
    if (merge_all) {
      K_star <- 1L
      pac_min <- if (all(is.na(pac_dt$PAC))) NA_real_ else min(pac_dt$PAC, na.rm = TRUE)
    } else {
      if (all(is.na(pac_dt$PAC))) {
        K_star <- 1L
        pac_min <- NA_real_
      } else {
        pac_min <- min(pac_dt$PAC, na.rm = TRUE)
        if (is.na(pac_min) || pac_min > pac_no_merge) {
          K_star <- 1L
        } else {
          eps <- 0.01
          cand <- pac_dt[!is.na(PAC) & PAC <= pac_min + eps, K]
          K_star <- min(cand)
        }
      }
    }

    pac_dt[, K_star := K_star]
    pac_curves[[blk]] <- pac_dt

    kstar_summary[[blk]] <- data.table::data.table(
      Block = blk, n_nodes = nblk, K_star = K_star,
      PAC_min = pac_min, MergeAll = merge_all, HomoQ2 = homo_q2, ValidN2 = valid_n2
    )

    # final clustering on full data if K*>1
    if (K_star == 1L) {
      final_membership_list[[blk]] <- stats::setNames(rep(1L, nblk), ids)
    } else {
      blk_cols <- blk_cols_list[[blk]]
      X_full_blk <- Matrix::crossprod(M_full[, blk_cols, drop = FALSE], M_full)
      X_full_blk <- as.matrix(X_full_blk)
      X_full_blk <- row_zscore(X_full_blk)

      if (nrow(unique(X_full_blk)) < K_star) {
        final_membership_list[[blk]] <- stats::setNames(rep(1L, nblk), ids)
        kstar_summary[[blk]][, `:=`(K_star = 1L, MergeAll = FALSE)]
        pac_dt[, K_star := 1L]
      } else {
        km_full <- stats::kmeans(X_full_blk, centers = K_star, nstart = kmeans_nstart_final, iter.max = 200)
        cl <- km_full$cluster
        names(cl) <- ids
        final_membership_list[[blk]] <- cl
      }
    }
  }

  pac_all   <- data.table::rbindlist(pac_curves, fill = TRUE)
  kstar_all <- data.table::rbindlist(kstar_summary, fill = TRUE)

  membership_dt <- data.table::rbindlist(lapply(names(final_membership_list), function(blk) {
    cl <- final_membership_list[[blk]]
    data.table::data.table(Block = blk, Disease = names(cl), Cluster = as.integer(cl))
  }), fill = TRUE)

  list(
    pac_curves = pac_all,
    kstar = kstar_all,
    membership = membership_dt,
    final_membership_list = final_membership_list
  )
}
