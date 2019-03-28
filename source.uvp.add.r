
#### add.ctd()
# Function adds CTD downcast values to a UVP cast. Parameters copied are manually specified in the loop
# below. Several interpolation methods are show. Only copies if the downcast datetime agrees with the 
# UVP datetime to within 30 minutes.
add.ctd = function(MasterList, ctd, verbose = FALSE) {
  if (verbose) {
    cat(paste('  Adding CTD data to', MasterList$Files$File, '.\n'))
  }
  
  ## Add the cast time
  ctd.times = unique(ctd$Time)
  dt = abs(as.numeric(difftime(MasterList$meta$Time, ctd.times, units = 'mins'))) ## Calc time between current cast and possible ctd.times
  
  if (min(dt) < 30) {
    l = which(ctd$Time == ctd.times[which.min(dt)])
    
    #Adding cycle to the metadata because this is a more convenient file to do
    #it with than the ecotaxa metadata
    MasterList$meta$Cycle = as.numeric(unique(gsub("[^0-9\\.]", "", ctd$Cycle[l[1]])))
    
    #### Option 1
    ## No interpolation - Just copy all of the data
    temp1 = data.frame(Depth = ctd$Depth[l],
                       S = ctd$Salinity[l],
                       T = ctd$Temperature[l],
                       Fl = ctd$Fluorescence.V[l],
                       O2 = ctd$Oxygen.uM[l],
                       Trans = ctd$Transmission[l])
    
    
    #### Option 3
    ## Interpolated to the same depth as the bin (and averaged)
    depths = MasterList$Flux$Depth
    temp3 = data.frame(Depth = depths, S = NA, T = NA, Fl = NA, O2 = NA, Trans = NA)
    
    for (j in 1:nrow(temp3)) {
      if (j == 1) {
        k = which(ctd$Depth[l] <= depths[j])
      } else {
        k = which(ctd$Depth[l] <= depths[j] & ctd$Depth[l] > depths[j-1])
      }
      temp3$S[j] = mean(ctd$Salinity[l[k]])
      temp3$T[j] = mean(ctd$Temperature[l[k]])
      temp3$Fl[j] = mean(ctd$Fluorescence.V[l[k]])
      temp3$O2[j] = mean(ctd$Oxygen.uM[l[k]])
      temp3$Trans[j] = mean(ctd$Transmission[l[k]])
    }
    
    ## Save Data
    MasterList$ctd$ctd.raw = temp1
    MasterList$ctd$ctd.inter = temp3
    
  } else {
    warning(paste0('Warning: no CTD cast available for UVP cast: ', MasterList[i]))
    MasterList$ctd = NA
    MasterList$meta$Cycle = NA
  }
  if (verbose) {
    cat('!! Done with CTD.\n')
  }
  
  ## Return
  MasterList
}


#### add.export()
# Function to merge UVP cast object with dataframes (usually loaded from an excel file) containing sediment
# trap or thorium export values. Column names are hardcoded below.Values are matched by datetime stamps.
add.export = function(MasterList, thorium, sedtrap) {
  
  l = which(sedtrap$Deployment.Time < MasterList$meta$Time & sedtrap$Recovery.Datetime > MasterList$meta$Time)
  
  if (length(l) > 0) {
    MasterList$export$sedtrap = sedtrap[l,]
  }
  
  ## Add Thorium data based on closest match to cast time.
  dt = abs(as.numeric(difftime(MasterList$meta$Time, thorium$Collection.Datetime, units = 'hours')))
  k = which(thorium$Collection.Datetime == thorium$Collection.Datetime[which.min(dt)])
  
  if (length(k) >  0) {
    MasterList$export$thorium = thorium[k,]
  }
  
  ## Return
  MasterList
}


#### add.export.fp()
# Function to merge UVP cast with a dataframe contianing FP observations. Entries are binned using the same bins
# as the UVP data.
add.export.fp = function(cruise, fp) {
  for (cast in 1:length(cruise)) {
    cruise[[cast]]$export$fp = NA
    
    l = which(fp$studyName == cruise[[cast]]$meta$Cruise & fp$Cycle == cruise[[cast]]$meta$Cycle)
    
    ## there are FP observations
    if (length(l) > 0) {
      
      ## Depths of those observations
      depths = unique(fp$Depth[l])
      
      ## Copy first row of Flux df for template, set depths
      cruise[[cast]]$export$fp = cruise[[cast]]$Flux[1:length(depths),]
      cruise[[cast]]$export$fp$Depth = depths
      
      ## For each depth and bin size, find entries
      for (d in 1:length(depths)) {
        for (bin in 1:(length(cruise[[cast]]$params$bins - 1))) {
          
          ## Set default value
          cruise[[cast]]$export$fp[d, bin + 1] = 0
          
          ## Entries we care about (l referenced)
          k = which(fp$ESD.um[l] >= cruise[[cast]]$params$bins[bin] & fp$ESD.um[l] < cruise[[cast]]$params$bins[bin+1] & fp$Depth[l] == depths[d])
          
          if (length(k) > 0) {
            cruise[[cast]]$export$fp[d, bin + 1] = sum(fp$Mass.pg[l[k]] * fp$Conversion.Factor[l[k]]) / 1e9 ## pg -> mg
          }
        }
      }
      ## Recalcualte total column
      cruise[[cast]]$export$fp$Total = 0
      cruise[[cast]]$export$fp$Total = apply(cruise[[cast]]$export$fp[,-1], 1, sum)
    }
  }
  
  ## Return
  cruise
}


#### add.spectra.model()
# Calculates the size-spectral slope for each bin and saves it in the UVP cast object. This function
# is designed to run on cruise objects (a list of multiple UVP cast objects).
add.spectra.model = function(cruise) { 
  for (cast in 1:length(cruise)){
    delta = cruise[[cast]]$params$delta.bins
    bins = cruise[[cast]]$params$bins[1:length(delta)]
    
    ## Initialize
    spectrum = data.frame(Depth = cruise[[cast]]$data$Depth, m = NA)
    
    for (i in 1:nrow(spectrum)) {
      x = bins + delta/2
      y = as.numeric(cruise[[cast]]$data[i,2:ncol(cruise[[cast]]$data)])[1:length(x)]
      
      y = log(y)
      x = log(x)
      
      x = x[is.finite(y)]
      y = y[is.finite(y)]
      
      if (length(y) > 0) {
        model = lm(y ~ x)
        spectrum$m[i] = as.numeric(model$coefficients[2])
      }
    }
    cruise[[cast]]$params$nslope = spectrum
  }
  
  #Return
  cruise
}
