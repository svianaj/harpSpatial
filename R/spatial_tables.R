#####################################
# DEFINITION OF VERIFICATION TABLES #
#####################################
# spatial_tables. Column names may not have a space or full stop.
# anyway, we hardcode per table

# return table description for spatial verification score tables
# @param tab a table name ("basic" or "fuzzy")
# not exported
spatial_score_table <- function(tab = c("basic", "fuzzy")) {
  tab <- match.arg(tab)

  standard_fields <- c(
    "model"    = "CHARACTER",
    "prm"      = "CHARACTER",
    "fcdate"   = "REAL",
    "fctime"   = "REAL",
    "leadtime" = "REAL"
  )
  extra_fields <- switch(tab,
    "basic" = character(0),
    "fuzzy" = c("scale" = "REAL", "threshold" = "REAL"),
    stop("Unknown table specification ", tab)
  )
  # BASIC:  scores without extra scales, thresholds etc.
  score_names <- switch(tab,
    "basic" = c("ets", "hk", "f", "bias", "mse", "S", "A", "L"),
    "fuzzy" = c("fss", "hk", "ets")
  )
#  score_fields <- structure(rep("REAL", length(score_names)), names=score_names)
  score_fields <- rep("REAL", length(score_names))
  names(score_fields) <- score_names
  result <- list(
    name = tab,
    fields = c(standard_fields, extra_fields, score_fields),
    primary = c(names(standard_fields), names(extra_fields))
  )
  result
}

# TODO: info table, ...

