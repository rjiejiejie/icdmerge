# Avoid R CMD check notes for data.table's non-standard evaluation (NSE)
utils::globalVariables(c(
  "Block", "Category_No", "Current_ID", "K", "K_star",
  "New_ID", "No", "PAC", "col_idx", "n_nodes"
))
