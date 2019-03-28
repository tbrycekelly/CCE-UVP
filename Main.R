## If these dont work for you then change your working directory
source('source.uvp.r')
source('source.uvp.add.r')
source('source.uvp.plots.r')
source('source.r')

#### Load the Data, C directory
input.dir = '../UVP Data/All PAR files/'
input.dir.meta = '../UVP Data/All Metadata files/'

raw.files = list.files(input.dir)                     #All years
meta.files = list.files(input.dir.meta)

#Reading in only first metadata file.
meta = as.data.frame(fread(paste0(input.dir.meta, meta.files[1])))

#Combining with the rest.
for (i in 2:length(meta.files)) {
  meta = rbind(meta, as.data.frame(fread(paste0(input.dir.meta, meta.files[i]))))
}


#### CTD and Sed Trap File section
#ctd.file = '../Other Data/CTD Downcast Data.xlsx'
#ctd = read.xlsx(ctd.file)
#ctd$Time = as.POSIXct(ctd$Datetime.GMT, tz = 'UTC') ## Get timestamps
#save(ctd, file = '../Other Data/CTD.rdata')

#sedtrap = read.xlsx('../Other Data/Sediment Trap.xlsx')
#sedtrap$Deployment.Time = as.POSIXct(sedtrap$Deployment.Time)
#sedtrap$Recovery.Datetime = as.POSIXct(sedtrap$Recovery.Datetime)
#save(sedtrap, file = '../Other Data/SedTrap.rdata')

# fp = read.xlsx('../Other Data/Sediment Trap Fecal Pellet Enumeration.xlsx')
# fp$ESD.um = 2 * (fp$Vol.um3 * 3 / 4 / 3.14159)^(1/3)
# fp = fp[!is.na(fp$ESD.um) & fp$ESD.um > 10 & fp$studyName > 705,]
# fp$studyName[fp$studyName == 1604] = 'CCEP1604'
# fp$studyName[fp$studyName == '810'] = 'CCEP0810'
# save(fp, file = '../Other Data/FP.rdata')

# gel =  read_excel("C:/Users/Christian/Dropbox/UVP/Other Data/Gel Traps/7.5X All.xls")
# save(gel, file = '../Other Data/Gel.rdata')

# BlendTh = read.xlsx('../Other Data/SedTrapThBlendedFlux_Trans.xlsx')
# BlendTh = BlendTh[BlendTh$Cruise > 705,]
# BlendTh$Cruise[BlendTh$Cruise == 810] = 'CCEP0810'
# BlendTh$Cruise[BlendTh$Cruise == 1106] = 'CCEP1106'
# BlendTh$Cruise[BlendTh$Cruise == 1208] = 'CCEP1208'
# BlendTh$Cruise[BlendTh$Cruise == 1408] = 'CCEP1408'
# BlendTh$Cruise[BlendTh$Cruise == 1604] = 'CCEP1604'
# save(BlendTh, file = '../Other Data/Blended Thorium.rdata')



### Load Ancillaries
load('../Other Data/SedTrap.rdata')
load('../Other Data/Thorium.interp.all.rdata')
load('../Other Data/CTD.rdata')
load('../Other Data/FP.rdata')
load('../Other Data/Gel.rdata')
load('../Other Data/Blended Thorium.rdata')

#sup = read.xlsx('../Other Data/Supplemental.xlsx')
#save(sup, file = '../rstates/Supplemental.rdata')

##############################
#### Run for one to test #####
##############################

##Initialize empty data structure, set Files
MasterList = init.raw(Dir = input.dir, File = raw.files[269]) # Picked a random file.
MasterList$params$dz = 5
MasterList$params$A = 13.7
MasterList$params$b = 1.55

## Load our data: set data.raw, fill in meta
MasterList = load.raw(MasterList, meta, TRUE)
MasterList = filter.by.size(MasterList, d.min = 0.102, d.max = 1.5) # Filtering to 2nd smallest size bin


k = 3 #Our current bins (increase for higher resolution, Guidi = 3)
bins = c(2^(c((5*k):(11*k))/k))

bins = floor(bins)
delta.bins = diff(bins)

## Save bin values
MasterList$params$bins = bins
MasterList$params$delta.bins = delta.bins
MasterList$params$n = length(bins)


##########################
## Time to get started:
##########################

MasterList = make.bins(MasterList)

## Calculate the flux (5m bins)
MasterList = calc.flux(MasterList)

## Add export values
MasterList = add.export(MasterList, thorium = thorium, sedtrap = sedtrap)

##Add ctd metadata
MasterList = add.ctd(MasterList, ctd, verbose = TRUE)

par(mfrow=c(1,1))
plot.raw.profile(MasterList)
plot.stacked(MasterList, ylim = c(0,100))
plot.ctd(MasterList)



######################
### Full Cruise ######
######################

##Looping through all the files for a given year
cruises = run.all.data(raw.files, A = 13.45, b = 1.35, dz = 5, ctd = ctd, verbose = FALSE, k = 3)

cruises = check.casts(cruises, TRUE)
cruises = add.spectra.model(cruises)
cruises = add.export.fp(cruises, fp)

save(cruises, file = '../rstates/cruises.rdata')
load('../rstates/cruises.rdata')
