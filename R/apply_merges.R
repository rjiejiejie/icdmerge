#' Apply merges based on membership and K* table
#'
#' @keywords internal
apply_merges <- function(current_data, lineage_dt, kstar_all, final_membership_list, verbose = TRUE) {

  current_data <- data.table::as.data.table(current_data)
  lineage_dt   <- data.table::as.data.table(lineage_dt)
  kstar_all    <- data.table::as.data.table(kstar_all)

  update_map <- data.table::data.table(Old_ID = character(), New_ID = character())

  for (blk in names(final_membership_list)) {
    row <- kstar_all[Block == blk]
    if (nrow(row) == 0) next

    K_star <- row$K_star[1]
    merge_all <- isTRUE(row$MergeAll[1])

    cl <- final_membership_list[[blk]]
    if (length(cl) < 2) next

    if (K_star > 1L) {
      clusters <- split(names(cl), cl)
      for (cid in names(clusters)) {
        members <- clusters[[cid]]
        if (length(members) <= 1) next

        old_records <- lineage_dt[Current_ID %in% members]
        if (nrow(old_records) == 0) next

        leader_code <- old_records$Original_Codes[1]
        prefix <- strsplit(leader_code, "\\|")[[1]][1]
        blk_safe <- sanitize_id(blk)

        new_id <- paste0(prefix, "_", blk_safe, "_K", K_star, "_C", cid, "_Merged")

        update_map <- rbind(update_map, data.table::data.table(Old_ID = members, New_ID = new_id), fill = TRUE)

        new_entry <- data.table::data.table(
          Current_ID = new_id,
          Original_Codes = paste(old_records$Original_Codes, collapse = "|"),
          Original_Names = paste(old_records$Original_Names, collapse = "|"),
          Block = blk,
          Member_Count = sum(old_records$Member_Count)
        )

        lineage_dt <- lineage_dt[!Current_ID %in% members]
        lineage_dt <- rbind(lineage_dt, new_entry, fill = TRUE)
      }

    } else if (K_star == 1L && merge_all) {
      members <- names(cl)
      old_records <- lineage_dt[Current_ID %in% members]
      if (nrow(old_records) > 0 && length(members) > 1) {
        leader_code <- old_records$Original_Codes[1]
        prefix <- strsplit(leader_code, "\\|")[[1]][1]
        blk_safe <- sanitize_id(blk)

        new_id <- paste0(prefix, "_", blk_safe, "_K1_AllMerged")

        update_map <- rbind(update_map, data.table::data.table(Old_ID = members, New_ID = new_id), fill = TRUE)

        new_entry <- data.table::data.table(
          Current_ID = new_id,
          Original_Codes = paste(old_records$Original_Codes, collapse = "|"),
          Original_Names = paste(old_records$Original_Names, collapse = "|"),
          Block = blk,
          Member_Count = sum(old_records$Member_Count)
        )

        lineage_dt <- lineage_dt[!Current_ID %in% members]
        lineage_dt <- rbind(lineage_dt, new_entry, fill = TRUE)
      }
    }
  }

  # apply mapping to current_data
  if (nrow(update_map) > 0) {
    current_data[, Category_No := as.character(Category_No)]

    to_chg <- current_data[Category_No %in% update_map$Old_ID]
    unchg  <- current_data[!Category_No %in% update_map$Old_ID]

    if (nrow(to_chg) > 0) {
      to_chg <- merge(to_chg, update_map, by.x = "Category_No", by.y = "Old_ID", all.x = TRUE)
      to_chg[, Category_No := New_ID][, New_ID := NULL]
    }

    current_data <- unique(rbind(unchg, to_chg, fill = TRUE))

  } else if (verbose) {
    cat("No merges found.\n")
  }
  if (verbose) {
    n_pairs <- nrow(update_map)
    n_new <- data.table::uniqueN(update_map$New_ID)
    cat(sprintf("Merging summary: mapped %d old nodes into %d new merged nodes.\n", n_pairs, n_new))
  }

  list(current_data = current_data, lineage_dt = lineage_dt, update_map = update_map)
}
