test_that("run_block_merge() smoke test on toy data", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("Matrix")

  library(data.table)
  library(icdmerge)

  set.seed(123)

  # ---- Toy data generator with block-wise comorbidity structure ----
  # Designed to be small + fast, but still likely to produce merges.
  simulate_toy <- function(
    n_patients = 300,
    blocks = c("Block_A", "Block_B", "Block_C"),
    diseases_per_block = 6
  ) {
    stopifnot(diseases_per_block %% 2 == 0)

    disease_ids <- unlist(lapply(seq_along(blocks), function(i) {
      sprintf("D%s_%02d", LETTERS[i], 1:diseases_per_block)
    }))
    block_vec <- rep(blocks, each = diseases_per_block)
    subgroup <- rep(rep(c(1, 2), each = diseases_per_block / 2), times = length(blocks))

    # patient latent factors
    g  <- rnorm(n_patients)  # global tendency
    bA <- rnorm(n_patients)
    bB <- rnorm(n_patients)
    bC <- rnorm(n_patients)
    sA1 <- rnorm(n_patients); sA2 <- rnorm(n_patients)
    sB1 <- rnorm(n_patients); sB2 <- rnorm(n_patients)
    sC1 <- rnorm(n_patients); sC2 <- rnorm(n_patients)

    block_factor <- function(blk) {
      if (blk == "Block_A") return(bA)
      if (blk == "Block_B") return(bB)
      bC
    }

    # baseline prevalence
    base <- runif(length(disease_ids), min = -2.6, max = -1.1)

    # weights: keep merges plausible
    w_global <- 0.8
    w_block  <- 1.1
    w_sub    <- 2.0

    # mild role shift
    role_shift <- rep(0, length(disease_ids))
    role_shift[block_vec == "Block_A" & subgroup == 1] <- 0.15
    role_shift[block_vec == "Block_B" & subgroup == 1] <- 0.10
    role_shift[block_vec == "Block_C" & subgroup == 2] <- 0.12

    eid <- sprintf("P%04d", 1:n_patients)
    out_list <- vector("list", length(disease_ids))

    for (j in seq_along(disease_ids)) {
      blk <- block_vec[j]
      bf <- block_factor(blk)

      sf <- switch(
        blk,
        Block_A = if (subgroup[j] == 1) sA1 else sA2,
        Block_B = if (subgroup[j] == 1) sB1 else sB2,
        Block_C = if (subgroup[j] == 1) sC1 else sC2
      )

      eta <- base[j] + w_global * g + w_block * bf + w_sub * sf + role_shift[j]
      p <- 1 / (1 + exp(-eta))
      y <- rbinom(n_patients, size = 1, prob = p)

      if (sum(y) > 0) {
        out_list[[j]] <- data.table(eid = eid[y == 1], Category_No = disease_ids[j])
      } else {
        out_list[[j]] <- NULL
      }
    }

    raw_data <- unique(rbindlist(out_list, use.names = TRUE, fill = TRUE))

    codebook <- data.table(
      No = disease_ids,
      block = block_vec,
      Codes = disease_ids,
      Names = paste("Toy disease", disease_ids)
    )
    freq_dt <- raw_data[, .(Frequency = uniqueN(eid)), by = Category_No]
    codebook <- merge(codebook, freq_dt, by.x = "No", by.y = "Category_No", all.x = TRUE)
    codebook[is.na(Frequency), Frequency := 0L]

    list(raw_data = raw_data, codebook = codebook)
  }

  toy <- simulate_toy()

  # ---- Run quickly ----
  # Keep resamples tiny to ensure test is fast.
  res <- run_block_merge(
    raw_data = toy$raw_data,
    codebook = toy$codebook,
    prevalence_cutoff = 0.02,   # higher cutoff for small toy
    n_resamples = 5,            # very small for speed
    k_max_per_block = 3,        # small K search
    write_files = FALSE,
    verbose = FALSE
  )

  # ---- Assertions: class + key fields ----
  expect_s3_class(res, "icdmerge_result")

  expected_fields <- c(
    "pac_curves", "kstar", "membership", "update_map",
    "merged_data", "lineage", "runtime_sec"
  )
  expect_true(all(expected_fields %in% names(res)))

  # kstar sanity
  expect_true(is.data.table(res$kstar))
  expect_true(all(c("Block", "n_nodes", "K_star") %in% names(res$kstar)))
  expect_true(all(res$kstar$K_star >= 1))

  # update map schema (may be empty, but should have columns)
  expect_true(is.data.table(res$update_map))
  expect_true(all(c("Old_ID", "New_ID") %in% names(res$update_map)))

  # merged_data schema
  expect_true(is.data.table(res$merged_data))
  expect_true(all(c("eid", "Category_No") %in% names(res$merged_data)))

  # lineage schema
  expect_true(is.data.table(res$lineage))
  expect_true(all(c("Current_ID", "Member_Count") %in% names(res$lineage)))
  expect_true(all(res$lineage$Member_Count >= 1))
})
