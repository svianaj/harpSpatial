#' Run spatial verification on a (for now) deterministic forecast
#'
#' @param start_date Date of the first forecast to read.
#' @param end_date Date of the last forecast to read.
#' @param det_model The name of the deterministic model. Maybe expressed as a
#'   vector if more than one model is wanted.
#' @param parameter The parameters to read as a character vector.
#' @param lead_time The lead times to read as a numeric vector.
#' @param by The time between forecasts. Should be a string of a number followed
#'   by a letter, where the letter gives the units - may be d for days, h for
#'   hours or m for minutes.
#' @param fc_file_path The top level path for the forecast files to read.
#' @param fc_file_template The file type to generate the template for. Can be
#'   "harmoneps_grib", "harmeoneps_grib_fp", "harmoneps_grib_sfx", "meps_met",
#'   "harmonie_grib", "harmonie_grib_fp", "harmone_grib_sfx", "vfld", "vobs", or
#'   "fctable". If anything else is passed, it is returned unmodified. In this
#'   case substitutions can be used. Available substitutions are {YYYY} for
#'   year, \{MM\} for 2 digit month with leading zero, \{M\} for month with no
#'   leading zero, and similarly \{DD\} or \{D\} for day, \{HH\} or \{H\} for
#'   hour, \{mm\} or \{m\} for minute. Also \{LDTx\} for lead time and \{MBRx\}
#'   for ensemble member where x is the length of the string including leading
#'   zeros - can be omitted or 2, 3 or 4. Note that the full path to the file
#'   will always be file_path/template.
#' @param fc_file_format The format of the files to read. Can be e.g. "fa" or "grib".
#' @param fc_options A list with format-specific options for the reader function.
#' @param fc_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param fc_accumulation The accumulation type of the forecast. This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h".
#' @param ob_file_path The top level path for the forecast files to read.
#' @param ob_file_template The file type to generate the template for. Can be
#'   "harmoneps_grib", "harmeoneps_grib_fp", "harmoneps_grib_sfx", "meps_met",
#'   "harmonie_grib", "harmonie_grib_fp", "harmone_grib_sfx", "vfld", "vobs", or
#'   "fctable". If anything else is passed, it is returned unmodified. In this
#'   case substitutions can be used. Available substitutions are {YYYY} for
#'   year, \{MM\} for 2 digit month with leading zero, \{M\} for month with no
#'   leading zero, and similarly \{DD\} or \{D\} for day, \{HH\} or \{H\} for
#'   hour, \{mm\} or \{m\} for minute. Also \{LDTx\} for lead time and \{MBRx\}
#'   for ensemble member where x is the length of the string including leading
#'   zeros - can be omitted or 2, 3 or 4. Note that the full path to the file
#'   will always be file_path/template.
#' @param ob_file_format The format of the files to read. Can be e.g. "hdf5" or "grib".
#' @param ob_options A list with format-specific options for the reader function.
#' @param ob_interp_method Interpolation method to be used when transforming a forecast
#'   field to the verification grid.
#' @param ob_accumulation The accumulation type of the observation (or reference). This is only used for
#'   accumulated parameters (e.g. precipitation). NULL signifies that the field is accumulated
#'   from the start of the model run. That is probably rare for observations. 
#'   Otherwise this should be a string containing a numerical value
#'   and a time unit, e.g. "15m" or "1h".
#' @param verif_domain A \code{geodomain} that defines the common verification grid.
#' @param return_data          = TRUE,
#' @param thresholds Thresholds used for FSS, ...
#' @param box_sizes Scales used for fuzzy methods like FSS. A vector of box sizes.
#'   All values must be odd integers (so the central point is really in the center of a box).
#' @param sqlite_path If specified, SQLite files are generated and written to
#'   this directory.
#' @param ... Arguments dependent on \code{file_format}, passed to to read_grid_interpolate.
#'   (More info to be added).
#'
#' @return A list containting two data frames: \code{spatial_scores} and
#'   \code{spatial_threshold_scores}.
#' @export

verify_spatial <- function(
  start_date,
  end_date,
  det_model,
  parameter,
  lead_time            = seq(0, 36, 3),
  lt_unit              = "h",
  by                   = "6h",
  fc_file_path         = "",
  fc_file_template     = "",
  fc_file_format       = "fatar",
  fc_options           = list(),
  fc_interp_method     = "closest",
  fc_accumulation      = NULL,
  ob_file_path         = "",
  ob_file_template     = "",
  ob_file_format       = "hdf5",
  ob_options           = list(),
  ob_interp_method     = "closest",
  ob_accumulation      = "15m",
  verif_domain         = NULL,
  use_mask             = FALSE,
  box_sizes            = c(1, 3, 5, 11, 21), # box sizes must be odd!
  thresholds           = c(0.1, 1, 5, 10),
  sqlite_path          = NULL,
  return_data          = TRUE,
  ...)
{
  prm <- harpIO::parse_harp_parameter(parameter)
  by_secs <- harpIO:::units_multiplier(by) * readr::parse_number(by)

  # For efficiency, we use a slightly counter-intuitive loop order
  # we don't loop over forecast date and then lead time,
  # because that would cause excessive re-reading (or caching) of observations.
  # Rather, we loop over all observation times 
  # and then over all forecasts valid for those times.
  # Close to start_date and end_date you must make sure not to beyond the time window.
  # TODO: to be even more efficient, you should try to
  #       - open (&parse) FC fields only once
  #       - read accumulated fields only once (e.g. "acc3h = 6h - 3h" also re-uses 
  #                                            the 3h field) 
  # Alternative strategy: loop by fcdate, and cache all obs in a list by leadtime
  #     next fcdate -> "ldt -= by" ; drop negative ldt ; read missing obs
  # For an accumulated variable (precip), the minimum lead time is 
  #     the accumulation time. Otherwise zero.

  # some date handling first : create vectors of date_time class
  sdate <- lubridate::ymd_hm(start_date, tz="UTC", truncated=4)
  if (missing(end_date)) {
    edate <- sdate
  } else {
    edate <- lubridate::ymd_hm(end_date, tz="UTC", truncated=4)
    # add the last fc hour if the date is just YYYYMMDD
    if (nchar(end_date) < 10) edate <- edate + (24 * 3600 / by_secs - 1) * by_secs
  }
  # convert lead_time to seconds and remove lead_times smaller than accum
  # we don't have 3h precip at 0h forecast.
  # also, we probably want lead_time in steps of the accumulation
  lt_scale <- harpIO:::units_multiplier(lt_unit)
  lead_time <- lead_time * lt_scale
  if (prm$accum > 0) {
    lead_time <- lead_time[which(lead_time >= prm$accum & lead_time %% prm$accum == 0)]
  }
  all_fc_dates <- seq(as.numeric(sdate), as.numeric(edate), by_secs) %>% 
                      lubridate::as_datetime()
  all_ob_dates <- (rep(all_fc_dates, each=length(lead_time)) + lead_time ) %>%
                      unique() %>%
                      sort()

  message(paste(all_fc_dates, collapse=" "))
  message(paste(all_ob_dates, collapse=" "))
  message(paste(lead_time / lt_scale, collapse=" "))
  # the re-gridding weights will come here:
  init <- list()

  # The read functions
  # FIXME: most read functions can't deal with accumulated parameters like AccPcp1h
  # we will need special "accumulator functions"

  get_ob <- function(obdate) {
    obfile <- get_filenames(file_date = format(obdate, "%Y%m%d%H"),
                            file_path = ob_file_path,
                            file_template=ob_file_template,
                            parameter = parameter)
    message("reading ", obfile)
    do.call(harpIO::read_grid, 
            c(list(filename=obfile, file_format=ob_file_format,
                   parameter = parameter), ob_options))
  }

  get_fc <- function(fcdate, lead_time) {
    fcfile <- get_filenames(file_date = fcdate,
                            lead_time = lead_time,
                            parameter = parameter,
                            det_model = det_model,
                            file_path = fc_file_path,
                            file_template = fc_file_template)
    message("reading ", fcfile, "lead_time ", lead_time)
    do.call(harpIO::read_grid,
            c(list(filename=fcfile, file_format=fc_file_format, 
                   parameter = parameter, lead_time=lead_time),
                       fc_options))
  }
  # We will write to SQL only at the end (more efficient),
  # Define the tables  and vectors at full size, so you don't have to add rows.
  ncases <- length(all_fc_dates) * length(lead_time)
  message("ncases= ", ncases)

  # Most efficient (?): create full table at the start
  # it may be too long if there are missing cases, but that is not so bad

  scores <- list()
  for(sc in c("basic", "fuzzy")) {
    t1 <- spatial_score_table(sc)
    nrows <- switch(sc,
                    "basic" = ncases,
                    "fuzzy" = ncases * length(thresholds) * length(box_sizes),
                    NA)
    template <- lapply(t1$fields, 
                       function(x) switch(x, 
                                          "CHARACTER" = NA_character_,
                                          "INTEGER"   = NA_integer_, 
                                          "REAL"      = NA_real_,
                                          NA_real_))
    scores[[sc]] <- do.call(tibble::tibble, c(template, .rows = nrows))
    scores[[sc]]$model <- det_model
    scores[[sc]]$prm   <- parameter
  }

  #scores_basic <- vector("list", ncases)
  #scores_fuzzy <- vector("list", ncases)
  # MAIN LOOP
  case <- 1
  for (ob in seq_along(all_ob_dates)) {  # (obdate in all_ob_dates) looses POSIXct class
    obdate <- all_ob_dates[ob]
    message("=====\nobdate:", format(obdate, "%Y%m%d-%H%M"))
    obfield <- get_ob(obdate)
    if (inherits(obfield, "try-error")) { # e.g. missing observation
      if (harpenv$verbose) cat("Observation not found. Skipping.\n")
      next
    }
    if (prm$accum > 0) { # an accumulated field like precipitation
      if (is.null(ob_accumulation) || ob_accumulation < 0) {
        # RARE: observation is an accumulated field (e.g. reference run)
        # TODO: This does not look very useful, unless get_ob() can somehow deal with it.
        warning("Accumulated observation fields not yet validated.", immediate.=TRUE)
        obfield <- obfield - get_ob(obdate - prm$accum)
      } else {
        ostep <- readr::parse_number(ob_accumulation) * harpIO:::units_multiplier(ob_accumulation)
        if (ostep == prm$accum) { # this is easy !
          # nothing to do
        } else if  (ostep > prm$accum) { # this is easy !
          stop("The chosen accumulation time is smaller than that of the observations!")
        } else {
          nstep <- prm$accum / ostep
          for (i in 1:(nstep-1)) {
            obfield <- obfield + get_ob(obdate - i*ostep)
          }
        }
      }
    }

    # convert to common verification grid
    if (!is.null(ob_interp_method)) {
      if (is.null(init$regrid_ob)) {
        message("Initialising ob regridding.")
        init$regrid_ob <- meteogrid::regrid.init(olddomain = obfield,
                                                 newdomain = verif_domain,
                                 method = ob_interp_method)
      }
      obfield <- meteogrid::regrid(obfield, weights = init$regrid_ob)
    }

    # find forecasts valid for this date/time
    # intersect drops the POSIXct class
    # valid_fc_dates <- intersect(obdate - lead_time, all_fc_dates)
    valid_fc_dates <- (obdate - lead_time)[which((obdate - lead_time) %in% all_fc_dates)]
    message("valid FC dates: ", paste(valid_fc_dates, collapse=" "))
    # inner loop
    for (fc in seq_along(valid_fc_dates)) {
      fcdate <- valid_fc_dates[fc]
      ldt <- (as.numeric(obdate) - as.numeric(fcdate)) # in seconds !
      message("   +++ fcdate = ", format(fcdate,"%Y%m%d-%H%M"),
              " +++ ldt = ", ldt/lt_scale, lt_unit)

      fcfield <- get_fc(fcdate, ldt/lt_scale)
      if (inherits(fcfield, "try-error")) { # e.g. missing forecast run
        if (harpenv$verbose) message("..... Forecast not found. Skipping.",
                                     .immediate = TRUE)
        next
      }
      # TODO: what if the forecast model needs "accumulating" rather than "decumulating"
      #       e.g. when verifying INCA against radar
      if (prm$accum > 0) {
        if (is.null(fc_accumulation) || fc_accumulation < 0) {
          if (ldt > prm$accum) { # if ldt==accum, you don't need to de-cumulate
            fcfield <- fcfield - get_fc(fcdate, (ldt - prm$accum)/lt_scale)
          }
        } else {
          stop("Can not (yet) create accumulated field ", parameter,
               " from model output fc_accumulation ", fc_accumulation)
          fstep <- readr::parse_number(fc_accumulation) * harpIO:::units_multiplier(fc_accumulation)
          if (fstep == prm$accum) { # this is easy !
          # nothing to do
          } else if  (fstep > prm$accum) { # this is easy !
            stop("The chosen accumulation time is smaller than that of the forecasts!")
          } else {
            nstep <- prm$accum / fstep
            for (i in 1:(nstep-1)) {
              fcfield <- fcfield + get_fc(fcdate, (ldt - i*fstep)/lt_scale)
            }
          }

        }
      }
      # convert forecast to common verification grid
      if (!is.null(fc_interp_method)) {
        if (is.null(init$regrid_fc)) {
          message("Initialising fc regridding.")
          init$regrid_fc <- meteogrid::regrid.init(olddomain = fcfield,
                                                   newdomain = verif_domain,
                                   method = fc_interp_method)
        }
        fcfield <- meteogrid::regrid(fcfield, weights = init$regrid_fc)
      }

      #############################
      ### NOW WE COMPUTE SCORES ###
      #############################

      # Basic (non-threshold)
      sc <- verify_basic(obfield, fcfield)
      scores[["basic"]]$fcdate[case]   <- fcdate
      scores[["basic"]]$leadtime[case] <- ldt/lt_scale
      scores[["basic"]][case, names(sc)] <- sc

      # "fuzzy" (threshold & box_size)
      sc <- verify_fuzzy(obfield, fcfield, thresholds, box_sizes)
      intv <- seq_len(dim(sc)[1]) + (case - 1) * dim(sc)[1]
      scores[["fuzzy"]]$fcdate[intv]   <- fcdate
      scores[["fuzzy"]]$leadtime[intv] <- ldt/lt_scale
      scores[["fuzzy"]][intv, names(sc)] <- sc

      case <- case + 1
    } # fcdate
  } #obdate
  if (case < ncases + 1) {
    message("There were ", ncases + 1 - case, "missing cases out of ", ncases, ".")
    ncases <- case - 1
  }

#  table_basic <- do.call(rbind, scores_basic)
#  table_fuzzy <- do.call(rbind, scores_fuzzy)

  ## write to SQLite
  if (!is.null(sqlite_path)) {
    db_file <- paste0(sqlite_path, "/harp_spatial.sqlite")
    message("Writing to SQLite file ", db_file)
    db <- harpIO:::dbopen(db_file)

    for (sc in c("basic", "fuzzy")) {
      # check for score table and create if necessary
      harpIO:::create_table(db, sc, scores[[sc]],
                            primary=c("model", "prm", "fcdate", "leadtime"))
      harpIO:::dbwrite(db, sc, scores[[sc]][1:ncases,])
    }

    harpIO:::dbclose(db)
  }

  if (return_data) invisible(scores)
  else invisible(NULL)
}


