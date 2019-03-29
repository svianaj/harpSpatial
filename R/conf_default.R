### we define a clean R environment for all HARP settings
### this can be passed and accessed very efficiently
### Make sure that all functions use harpenv as their run environment
### so they know exactly the variables defined in it
### fcdate, obdate etc are NOT in this environment.
### But begindate, enddate, prm are.

###################################################
# PATH statements: hard coded or from environment #
###################################################

if (!exists("HARPDIR")) {
  HARPDIR = Sys.getenv("HARPDIR")
  if (HARPDIR=="") HARPDIR = paste0(Sys.getenv("HOME"),"/Harp")
}
RHOME = paste0(HARPDIR,"/spatial/R/")
CONFDIR = paste0(HARPDIR,"/spatial/conf/")

time_unit = "h"  # FC range & lead times etc expressed in hours. Alternatives: "m" (minutes), "s" (seconds)
                 # internally, HARP uses seconds, time_unit is used in SQLite output

verbose = TRUE   # option -v
prm = "AccPcp1h" # OPTIONAL: you may (re-)set this with HRP_PARAMETER or -p command line option

# read domains BEFORE using them!
dfile = file.path(CONFDIR,"domains_default.R")
if (file.exists(dfile)) source(dfile, local=TRUE)

# SQLite uses file locking to check whether a file is opened by another user before writing
# You *must* turn this off on NFS and Lustre file systems (file locking doesn't work well on network file systems)
# But: this means that two processes are allowed to write at the same time, causing corruption of the db file.
# FUTURE: maybe make this VF_ specific, to have different values for FC_ and OB_
sqlite.lock = TRUE # set FALSE if the file is on a NFS or Lustre file system!


##################
# FORECAST FILES #
##################

FC_path = paste0(Sys.getenv("SCRATCH"),"/FC")
FC_format = "fa" # fa, grib, hdf5, inca_ascii, ...

FC_model = "default"  # could be used as SQL identifier or in FC_file. Can be re-set with -m option

# a function that returns the path/name of a single FC file
FC_file = function(fcdate, ldt=0, ...) { # note that ldt is always in SECONDS
  fullpath = paste0(FC_path, "/", format(fcdate,"%Y/%m/%d"))
  filename = sprintf("FC%s+%04i",format(fcdate,"%Y%m%d%H"), ldt/3600)
  paste0(fullpath, "/", filename)
}

FC_range = 54 # (*3600) forecast range in time_unit
FC_inc = 1    # (*3600) time increment of FC files
FC_cycle_inc = 24 # 1 forecast/day

# ONLY used for accumulated variables like precipitation (ignored in other cases)
FC_interval = 0 # for what time interval (in time_unit) is the value valid?
                # <=0 value means that the value is accumulated from forecast start
                # >0 means the value is valid for a period of X seconds. E.g. a radar scan may have 300
                # but a grib file with 1h accumulated precipitation would have 3600

####################################
# OBSERVATIONS (radar, nowcast...) #
####################################

OB_path = paste(Sys.getenv("SCRATCH"), "OB", sep="/")
OB_format = "hdf5" # FA, grib, grib1, grib2, hdf5...

OB_file = function(obdate) {
  fullpath = paste(OB_path, format(obdate,"%Y/%m/%d"), sep="/")
  filename = sprintf("OB%s.hdf",format(obdate,"%Y%m%d%H%M"))
  paste(fullpath, filename, sep="/")
}

OB_interval = 300/3600 # frequency in *time_unit* (e.g. hours) of the obs files (300 = every 5 minutes), 
                       # <=0 NOT POSSIBLE in OBS
                       # Only used for accumulated variables like precipitation

# we define different decoder formats for different kind of input data. much cleaner!
OB_rawtype = "reflectivity" 
# "cm" (cloudmask) "lwr" (MSG cloud temp), dBZ, rainrate, ...

# data path inside hdf5:
OB_datapath = "dataset1/data1/data" # ODIM default for hdf5

# for non-standard binary formats: write your own decoding function and domain specification
#OB_domain = 
#OB_xxx_decode

#OB_sqlite_path = 
#OB_sqlite_file =  function()

################
# VERIFICATION #
################

VF_sqlite_path = "/tmp"
VF_sqlite_file = function(...) {
  fullpath = VF_sqlite_path
  filename = "HARP_spatial.sqlite"
#  filename = sprintf("VERIF_%s_%s-%s.sqlite",
#                     prm, format(begindate,"%Y%m%d"),format(enddate,"%Y%m%d"))
  paste(fullpath, filename, sep="/")
}

#we assume : all leadtimes, and for precip steps of "accum", not FC_inc
# time step for verification: usually either equal to FC_inc or to the accumulation period of the variable
# But here you can override it. For instance if you have hourly files with 6h accumulated rainfall.
# DEFAULT: FC_inc for non-accumulated variables and accumulation time for e.g. rain.
# Can never be smaller than FC_inc, of course.
#VF_step = FC_inc


VF_thresholds = list(
     "T2m"=c(-10, 0 , 10, 20, 30),
     "AccPcp1h"=c(0.1, 1, 5, 10),
     "AccPcp3h"=c(0.1, 1, 5, 10, 20))

## QUESTION: is it worth making the score lists flexible?
# local methods: aggregated over all the grid points
# some of these correspond simply to a fuzzy score with window size 1
VF_list = c("basic","fuzzy","sal")
VF_window_sizes = c(1, 3, 5, 11, 21) # box sizes must be odd!

# formatting of dates when writing to SQLite
# DON'T CHANGE UNLESS YOU REALLY KNOW WHAT YOU ARE DOING
VF_fcdate_format = function(x) as.numeric(format(x,"%Y%m%d"))
VF_fctime_format = function(x) as.numeric(format(x,"%H"))


############################
# COMMON VERIFICATION GRID #
############################

# common domain for FC and OBS?
# also specify interpolation (upscaling...) method

#common_grid = 

FC_regrid_method = "bilin" # bilinear interpolation
OB_regrid_method = "mean"  # UPSCALING: take mean of all values falling in the new grid box

# how convert to common grid?
# either define FC_to_common() and OB_to_common() explicitely here, or load the defaults:
source(file.path(RHOME,"interpolators.R"), local=TRUE) # this defines FC_to_common() and OB_to_common()

#################################
# other configurations and code #
#################################

source(file.path(RHOME,"db_tables.R"), local=TRUE)
source(file.path(RHOME,"scores.R"), local=TRUE)
source(file.path(RHOME,"decoders.R"), local=TRUE)
# add any extra configuration files (local file formats etc.)

#######
# END #
#######

