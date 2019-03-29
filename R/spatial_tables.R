#####################################
# DEFINITION OF VERIFICATION TABLES #
#####################################
# cspatial_tables. Column names may not have a space or full stop.
# anyway, we hardcode per table
#' return table description for spatial verification score tables

spatial_score_table <- function(tab) {
  standard_fields <- c("model"="CHARACTER", "prm"="CHARACTER",
                       "fcdate"="REAL", "leadtime"="REAL")

  # BASIC:  scores without extra scales, thresholds etc.
  column_names <- switch(tab,
                    "basic" = c("ets", "hk", "f", "bias", "mse", "S", "A", "L"),
                    "fuzzy" = c("scale", "threshold", "fss", "hk", "ets"),
                    stop("Unknown table specification ", tab))
  columns <- structure(rep("REAL", length(column_names)), names=column_names)
  result <- list(name=tab,
                 fields = c(standard_fields, columns),
                 primary = names(standard_fields))
  # tibble template
  result$template <- lapply(result$fields, function(x) switch(x, "CHARACTER"="NA_character_", "REAL"="NA_real_", NA))
  result
}
# TODO: info table, ...

