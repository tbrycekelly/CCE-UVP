library(data.table)
library(RColorBrewer)
library(readxl)
library(openxlsx)

#### init.raw()
# Generates a blank UVP cast object. Only values will be the file name values used to load the data in load.raw()
init.raw = function(Dir, File){
  
  ## Return
  list(data.raw = data.frame(), 
       data = data.frame(), 
       Flux = data.frame(),
       meta = list(), 
       ctd = list(ctd.raw = NA, ctd.inter = NA),
       export = list(sedtrap = NA, thorium = NA),
       params = list(A = NA,
                     b = NA,
                     dz = NA, #depth bin width
                     bins = rep(NA, times = 1), #Vector of left bounds of size bins
                     delta.bins = NA,#Width of size bins
                     n = NA#number of size bins
       ),
       Files = list(Dir = Dir, File = File, Time = date()) ## Information about the raw data files
  )
}

#### load.raw()
# Main loading function for the raw data file. Loads the raw Particle file, matches the relevant metadata, and
# calculates ESD. As with all primary functions, it takes a UVP cast object and returns an appended UVP cast object.
load.raw = function(MasterList, meta, verbose = FALSE) {
  
  if(verbose) {
    cat(paste('Attempting to load:\n', MasterList$Files$File))
  }
  
  ## Read the data
  data = as.data.frame(fread(paste0(MasterList$Files$Dir, MasterList$Files$File)))
  
  metadata = list()
  
  ## Deal with metadata information (copy over the good stuff)
  l = which(meta$'Particle filename' == MasterList$Files$File)
  if (length(l) == 1) {
    metadata$Lat = meta$latitude[l]
    metadata$Lon = meta$longitude[l]
    metadata$Time = as.POSIXct(paste0(substr(MasterList$Files$File, 1, 4), '-',
                                      substr(MasterList$Files$File, 5, 6), '-',
                                      substr(MasterList$Files$File, 7, 8), ' ',
                                      substr(MasterList$Files$File, 9, 10), ':',
                                      substr(MasterList$Files$File, 11, 12), ':',
                                      substr(MasterList$Files$File, 13, 14)), tz = 'UTC')
    
    metadata$Pixel = meta$acq_pixel[l]  ## [mm] / px
    metadata$MetricConvCoef = meta$acq_aa[l]
    metadata$MetricConvExp = meta$acq_exp[l]
    metadata$Stn = meta$stationid[l]
    metadata$Cruise = meta$cruise[l]
    metadata$ImgVol = meta$acq_volimage[l]
  }
  else if (length(l) > 1) {
    metadata$Lat = meta$latitude[l[1]]
    metadata$Lon = meta$longitude[l[1]]
    metadata$Time = as.POSIXct(paste0(substr(MasterList$Files$File, 1, 4), '-',
                                      substr(MasterList$Files$File, 5, 6), '-',
                                      substr(MasterList$Files$File, 7, 8), ' ',
                                      substr(MasterList$Files$File, 9, 10), ':',
                                      substr(MasterList$Files$File, 11, 12), ':',
                                      substr(MasterList$Files$File, 13, 14)), tz = 'UTC')
    
    metadata$Pixel = meta$acq_pixel[l[1]]
    metadata$MetricConvCoef = meta$acq_aa[l[1]]
    metadata$MetricConvExp = meta$acq_exp[l[1]]
    metadata$Stn = meta$stationid[l[1]]
    metadata$Cruise = meta$cruise[l[1]]
    metadata$ImgVol = meta$acq_volimage[l[1]]
    warning(paste0('Multiple metadata matches: ', MasterList$Files$File, ' (', length(l), ')'))
  }
  else {
    metadata$Lat = NA
    metadata$Lon = NA
    metadata$Time = as.POSIXct(paste0(substr(MasterList$Files$File, 1, 4), '-',
                                      substr(MasterList$Files$File, 5, 6), '-',
                                      substr(MasterList$Files$File, 7, 8), ' ',
                                      substr(MasterList$Files$File, 9, 10), ':',
                                      substr(MasterList$Files$File, 11, 12), ':',
                                      substr(MasterList$Files$File, 13, 14)), tz = 'UTC')
    metadata$Pixel = NA
    metadata$MetricConvCoef = NA
    metadata$MetricConvExp = NA
    metadata$Stn = NA
    metadata$Cruise = NA
    metadata$ImgVol = -1
    warning(paste0('No metadata matches: ', file))
  }
  
  ## Fix cruise names (manual fix):
  if (metadata$Cruise == 'lter2008') { metadata$Cruise = 'CCEP0810'}
  if (metadata$Cruise == 'ccelter_2011') { metadata$Cruise = 'CCEP1106'}
  if (metadata$Cruise == 'ccelter_2012') { metadata$Cruise = 'CCEP1208'}
  if (metadata$Cruise == 'ccelter_2014') { metadata$Cruise = 'CCEP1408'}
  if (metadata$Cruise == 'ccelter_2016') { metadata$Cruise = 'CCEP1604'}
  if (metadata$Cruise == 'sn201_ccelter_2017') { metadata$Cruise = 'CCEP1706'}
  
  if (verbose) {
    print(metadata)
  }
  
  metadata$Cycle = NA
  
  ## Basic Calculations
  data$area.converted = metric.conversion(data$area, metadata$MetricConvCoef, metadata$MetricConvExp) ## Pixel [px2] to [mm2]
  data$ESD = sqrt(data$area.converted / 3.141593) * 2
  
  ## Remove extra fields
  data$greylimit1 = NULL
  data$greylimit2 = NULL
  data$greylimit3 = NULL
  
  if(verbose) {
    cat(' --> File loaded. \n\n')
  }
  
  ## Pack up
  MasterList$data.raw = data
  MasterList$meta = metadata
  
  ## Return
  MasterList
}

#### make.bins
# Binning function that takes the raw particle data contained in the UVP data object and bins it accoring
# to the bin parameters in the object.
make.bins = function(MasterList, verbose = FALSE) {
  ## unpack & initialize
  data = MasterList$data.raw
  bins = MasterList$params$bins
  particles = data.frame(Depth = 0)
  
  depths = ceiling(MasterList$data.raw$depth / MasterList$params$dz) ## Rounds up so 1,2,3,4,5 = depth 5 -----  6,7,8,9,10 = 10
  
  ## Initialize New Bin Columns
  for (bin in bins) {
    particles[[paste0('Bin', bin)]] = NA
  }
  
  ## For each depth
  for (i in unique(depths)) {
    
    #n = delta C/delta d
    n = rep(0, length(bins))
    k = which(depths == i) # All entries we want at depth x
    vol = 0
    for (j in unique(depths[k])) {
      kk = which(depths == j)
      vol = vol + data$imgcount[kk[1]] * MasterList$meta$ImgVol
    }
    
    ## For each bin size
    for (j in 2:length(bins)) {
      ## Find which entries are at the right depth and within the bin size
      l = which(depths == i & data$ESD * 1e3 >= bins[j-1] & data$ESD * 1e3 < bins[j]) 
      
      delta_C = sum(data$nbr[l], na.rm = TRUE) / vol
      delta_s = MasterList$params$delta.bins[j-1] / 1000
      n[j-1] = (delta_C/delta_s)
    }
    
    ## Add new row and place proper data
    particles = rbind(particles, particles[1,])
    particles$Depth[nrow(particles)] = i * MasterList$params$dz #Depth
    particles[nrow(particles), c(2:(1+length(bins)))] = n
  }
  
  ## Add particle df to list
  MasterList$data = particles[-1,]
  ## Return
  MasterList
}


#### calc.flux
# Take UVP cast object and calcuclate the Fluxes from the binned data and the stored A and b values
# following Guidi et al. (2008)
calc.flux = function(MasterList, verbose = FALSE) {
  
  Flux = data.frame(Depth = MasterList$data$Depth)
  
  #If no alternative A or b given, sets to Guidi 2008 default
  if (is.na(MasterList$params$A)) {
    MasterList$params$A = 12.5
    MasterList$params$b = 3.81
  }
  
  
  ## Initialize the columns
  for (i in 1:(length(MasterList$params$bins) - 1)) {
    Flux[[paste0('Flux', MasterList$params$bins[i])]] = NA
  }
  
  ## Actually calculate values
  for (row in 1:nrow(Flux)) {
    
    ## For each particle class...
    for (i in 1:(length(MasterList$params$bins) - 1)) {
      diameter = (MasterList$params$bins[i] + (MasterList$params$bins[i+1] - MasterList$params$bins[i]) / 2) #Bins in um. Guidi also used midpoint of the bins so potential improvement is to use true mean diameter.
      
      flx = flux.model(n = MasterList$data[row,i+1],
                       d = diameter/1000,
                       A = MasterList$params$A,
                       b = MasterList$params$b)  * (MasterList$params$delta.bins[i]/1000) #diameter back to mm
      
      Flux[[paste0('Flux', MasterList$params$bins[i])]][row] = flx / MasterList$params$dz
    }
  }
  
  ## Add total
  Flux$Total = apply(Flux[,2:ncol(Flux)], 1, sum)
  
  ## Pack
  MasterList$Flux = Flux
  
  ## Return
  MasterList
}


#### filter.by.size()
# As the name implies, this filters the raw particle data based on ESD cutoffs. 
filter.by.size = function(MasterList, d.min = 0, d.max = 1.5) {
  
  MasterList$data.raw = MasterList$data.raw[MasterList$data.raw$ESD >= d.min & MasterList$data.raw$ESD <= d.max,] #As long as we use Guidi's A and b we should use the size range for which they were calibrated.
  MasterList$params$d.min = d.min
  MasterList$params$d.max = d.max
  
  ## Return
  MasterList
}


#########################
#### Cruise objects  ####
#########################

#### run.all.data
# A wrapper function to automate the UVP cast creation. Compiles a list of raw.file names 
# (i.e. a directory of files) into a cruise object. k is a sizing parameter for the bin resolution.
run.all.data = function (raw.files, A, b, dz, ctd, k = 3, verbose = FALSE){
  cruise = list()
  
  for (i in 1:length(raw.files)){
    
    ##Initialize empty data structure, set Files
    MasterList = init.raw(Dir = input.dir, File = raw.files[i])
    MasterList$params$dz = dz
    MasterList$params$A = A
    MasterList$params$b = b
    
    ## Load our data: set data.raw, fill in meta
    MasterList = load.raw(MasterList, meta, verbose = verbose)
    MasterList = filter.by.size(MasterList, d.min = .102, d.max = 1.500)
    
    k = k #Our current bins (increase for higher resolution, Guidi = 3)
    bins = c(2^(c((5*k):(11*k))/k))
    
    bins = floor(bins)
    delta.bins = diff(bins)
    
    ## Save bin values
    MasterList$params$bins = bins
    MasterList$params$delta.bins = delta.bins
    MasterList$params$n = length(bins)
    
    ## Time to get started:
    MasterList = make.bins(MasterList)
    
    ## Calculate the flux (5m bins)
    MasterList = calc.flux(MasterList)
    
    ## Add export values
    MasterList = add.export(MasterList, thorium = thorium, sedtrap = sedtrap) ## Currently only works for P0810
    
    MasterList = add.ctd(MasterList, ctd)
    
    ## Save
    cruise[[paste0('Cast', substr(MasterList$Files$File, 1, 14))]] = MasterList
  }
  #Return
  cruise
}

#### check.casts()
# An ad-hoc function for checking that each UVP cast object is associated with a lagrangian cycle number.
# Probably not needed for most implementations of this code.
check.casts = function(cruise, verbose = FALSE) {
  
  names = function() {
    out = c()
    
    for (cast in 1:length(cruise)) {
      if (!is.na(cruise[[cast]]$meta$Cycle)) {
        if (!cruise[[cast]]$meta$Cycle %in% out) {
          out = c(out, cruise[[cast]]$meta$Cycle)
        }
      }
    }
    ## Return
    out
  }
  
  if (verbose) {
    print('    Before: ')
    print(names())
  }
  
  k = c()
  for (cast in 1:length(cruise)) {
    if (is.na(cruise[[cast]]$meta$Cycle) | cruise[[cast]]$meta$Cycle > 10 | cruise[[cast]]$meta$Cycle < 1) {
      k = c(k,cast)
    }
  }
  if (length(k) > 0) {
    cruise = cruise[-k]
  }
  
  if (verbose) {
    print('    After: ')
    print(names())
  }
  
  ## Return
  cruise
}


############################
#### Helper Functions  #####
############################

#### metric.conversion()
# Apply metric conversion as in Stemman et al. (2002).
metric.conversion = function(x, coeff, exp) {
  ## Converts UVP area to actual size
  #
  # From Stemmann et al. 2002:
  # "The metric surface (Y) as a function of the pixel 
  #  surface (X)is Y = 0.00139 * X^1.43 (Stemmann, 1998)"
  #  1998 is his inaccessible thesis so this is what we have.
  #  This calibration is also for UVP2. May have changed.
  #
  #  Update: Picheral et al. 2010 recalculated the coefficient and 
  #  and exponent values for the above power law specifically for UVP
  #  5. I've updated the equation thusly.
  #
  # I think this converts x [area in px2] to y [area in mm2]
  #0.003 * x ^ 1.3348
  coeff * x ^ exp
}

#### flux.model()
#Flux calculation per bin
#From meeting with Mike on 1/16/19. According to Guidi 2008 
#n = deltaC/deltaS (we call it delta d). Right now we are treating n 
#as just the particle concentration, so when we go to multiply delta d
#the units for flux become mgC/m2/d/um. The delta d is to remove that unit,
#not add it.
flux.model = function(n, d, A, b) {
  n * A * d ^ b
}

