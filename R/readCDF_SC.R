#########################################################################
#     GcDuo - R package for GCxGC processing and PARAFAC analysis
#     Copyright (C) 2023 Maria Llambrich & Lluc Semente
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################


##########################
#' readCDF2
#'
#' Read a single CDF file and transform the data into a 3D array.
#'
#' @param pathFolderCDF path to the folder containing the CDF to process.
#' @param modulationTime modulation time used in the data.
#' @param mzRange range values used in the m/z data acquisition.
#' @param ncore number of parallel processes to start. Defaults to number of available cores - 1.
#' @return A 3 dimensional array containing all the data in the files following the structure: |sample, mz, time2D, time1D|
#' @export

readCDF2 <- function(filePath, modulationTime, mzRange, ncore=NULL)
{
  #browser()
  suppressWarnings(library(doParallel))
  # File information Extraction
  nc <- ncdf4::nc_open(filePath)
  #variable for intensity values
  int <- ncdf4::ncvar_get(nc, "intensity_values", start = 1, count = -1)
  #variable acquisition time
  time <- ncdf4::ncvar_get(nc, "scan_acquisition_time")
  #variable mass to charge ratio values
  mz <- ncdf4::ncvar_get(nc, "mass_values")
  #variable total intensity
  tic <- ncdf4::ncvar_get(nc, "total_intensity")
  point <- ncdf4::ncvar_get(nc, "point_count")
  #close connection to cdf file
  ncdf4::nc_close(nc)

  ## int is intensity for m:z; mz is mass values for m:z
  ## tic is total ion count; time is acquisition time; point is point count. Mix of D1 and D2?

  ## Number of time points recorded
  N <- length(time)

  #Convert time to minutes
  time <- time/60
  ## Difference between time steps (mins) == frequency
  f_mostr <- abs(time[1] - time[2]) #frequency
  ## Modulation time as a proportion of gap between time steps
  ## I have removed the rounding step
  cycles <- (modulationTime/60)/f_mostr

  # frequency * modulation = dimension retention time 2
  dim2d <- cycles
  # dimension retention time 1
  dim1d <- N/ cycles

  mz_nodec <- round(mz,0)
  if(tail(sort(unique(mz_nodec)),1) - tail(sort(unique(mz_nodec)),2)[1] == 1){
    mz_seq <- seq(min(unique(mz_nodec)), max(unique(mz_nodec)))
  } else {
    mz_seq <- seq(min(unique(mz_nodec)), tail(sort(unique(mz_nodec)),2)[1] + 1)
  }

  period <- (max(mz_seq) + 1) - min(mz_seq) # mz period
  mz_range <- range(mz_seq) #mz range
  mzRange <- mzRange

  # clusters
  if (is.null(ncore)){
    ncore <- parallel::detectCores(logical = TRUE) - 1
  }
  cl <- parallel::makeCluster(ncore)
  doSNOW::registerDoSNOW(cl)

  # progress bar
  pb <- txtProgressBar(max = length(time), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # processing iteration
  if (foreach::getDoParRegistered())
  {
    int_list <- foreach::foreach(i = iterators::icount(length(time)), .options.snow = opts) %dopar% getSpectrumIntensity(mzs = mz, points =  point, index =  i, ints = int, mz_seq = mz_seq)
  }
  parallel::stopCluster(cl)


  # Shape data conveniently
  gc_data <- array(dim = c(length(int_list), length(int_list[[1]])))
  for (spec in 1:nrow(gc_data))
  {
    gc_data[spec, ] <- int_list[[spec]]
  }
  rm(int_list)


  ## Is the first time point greater than the modulation time? If so...
  if(time[1] > modulationTime/60) {
    ## Divide the first time point by the modulation time in sec
    need_pos <- (time[1]/(modulationTime/60))
    start_pos <- (need_pos - trunc(need_pos))*cycles
    # Filling time and intensity axis with zeros to ensure the folding is correct

    ## Does this add in time points missing if the start is shorter than the modulations??
    needval <- seq(0, (time[1]), length.out = need_pos) #correct time at the endings

    filled.data<-c(rep(min(gc_data), length.out = length(needval)*length(mz_seq)), as.vector(t(gc_data)))
  } else {
    filled.data<-as.vector(t(gc_data))[1:(round(period,0) * round(dim2d,0) * round(dim1d,0))]
  }
  # Filling time and intensity axis with zeros to ensure the folding is correct
  #needtime <- seq(0, (time[1]), by = mean(diff(time))) #correct time at the endings
  # posdelay <- length(needtime)/dim2d
  # misspos <- (trunc(posdelay) - 1) * dim2d
  # needpos <- length(needtime) - misspos
  # time_cor <- c(tail(needtime, n = needpos), time)
  time_pos <- seq(1, length(time), by = dim2d)
  time_values1d <- time[time_pos] #get time values retention time 1
  # dim1d <- length(time_values1d)
  time_values2d <- round(seq(0, (modulationTime - (modulationTime/dim2d)), length.out = round(dim2d, 0)), 3) #get time values retention time 2

  #create the 3D array

  ####
  array3D <- array(data = filled.data,
                   ## period is m-z range, dim2d is acquisition values for second dimension time, dim1d is acquisition values for first dimension time
                   dim = c(round(period,0), round(dim2d,0), round(dim1d,0)),
                   dimnames = list(mz_seq, time_values2d, time_values1d))

  if (length(mz_seq) > length(seq(mzRange[1], mzRange[2]))) {
    mz1 <- which(mz_seq == mzRange[1])
    mz2 <- which(mz_seq == mzRange[2])

    message(paste("mz dimension reduced from ", mz_range[1], "-", mz_range[2],
                  "to ", mzRange[1], "-", mzRange[2]))

    array3D <- array3D[mz1:mz2,,]

  } else {
    mz1 <- 1
    mz2 <- length(mz_seq)
  }

  rm(gc_data)

  results <- list(data = array3D, mz_seq = mz_seq[mz1:mz2], time_values2d = time_values2d, time_values1d = time_values1d)

  return(results)
}

