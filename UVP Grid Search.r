source('source.uvp.r')
source('source.uvp.add.r')
source('source.uvp.plots.r')
source('source.uvp.post.r')
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

### Load Ancillaries
load('../Other Data/SedTrap.rdata')
load('../Other Data/Thorium.interp.all.rdata')
load('../Other Data/CTD.rdata')
load('../Other Data/FP.rdata')
load('../rstates/Supplemental.rdata')


grid = expand.grid(A = seq(10, 18, by = 1), b = seq(1.5, 4.5, by = 0.5), score = 0, r2 = 0)
grid = expand.grid(A = seq(12.2, 13.8, by = 0.15), b = seq(1.4, 1.6, by = 0.025), score = 0, R2 = 0)
verbose = FALSE


for (i in 1:nrow(grid)) {
  print(paste('Grid entry ', i, ' (', Sys.time(), ')  ::   A = ', grid$A[i], '  B = ', grid$b[i]))
  
  k = 3
  ##Looping through all the files for a given year
  cruises = run.all.data(raw.files, A = grid$A[i], b = grid$b[i], dz = 5, ctd = ctd, verbose = FALSE, k = k)
  
  if (verbose) {
    print('Raw data loaded.')
  }
  
  cruises = check.casts(cruises, TRUE)
  cruises = add.spectra.model(cruises)
  cruises = add.export.fp(cruises, fp)
  
  if (verbose) {
    print('Casts checked.')
  }
  
  st.comp = build.st.comparison(cruises, average = TRUE)
  #fp.comp = build.fp.comparison(cruises, average = TRUE)
  #fp.st.comp = build.fp.st.comparison(cruises, average = TRUE)
  
  if (verbose) {
    print('Comparison Generated.')
  }
  
  
  #Plot cycle averaged
  png(filename = paste0("../Images/Grid/New Small (11) -", i,".png"))
  
  plot(st.comp$ST, st.comp$Total, pch=16, xlim=c(0,600), ylim=c(0,600), xaxs='i', yaxs='i',
       ylab = 'UVP Flux', xlab = 'ST Flux', main = paste0('A = ', cruises[[1]]$params$A, '     B = ', cruises[[1]]$params$b))
  add.one.to.one()
  add.error.bars(st.comp$ST, st.comp$ST.sd, st.comp$Total, st.comp$Total.sd)
  
  model = bootstrap(st.comp$ST, st.comp$ST.sd, st.comp$Total, st.comp$Total.sd, n = 1000)
  add.boot.conf(model, trendline = TRUE)
  
  score = sum((log(st.comp$Total) - log(st.comp$ST))^2) ## Guidi's Score
  
  #score = sum(((st.comp$Total-st.comp$ST)/st.comp$ST.sd)^2, na.rm = TRUE)
  r2 = 1 - sum((st.comp$Total - st.comp$ST * median(model$m) + median(model$b))^2) /
    sum((st.comp$Total - mean(st.comp$Total, na.rm = TRUE))^2)
  
  grid$score[i] = score
  grid$r2[i] = r2
  
  text(50, 575, paste0('m = ', floor(median(model$m)*100)/100))
  text(50, 545, paste0('b = ', floor(median(model$b)*100)/100))
  text(50, 515, paste0('S = ', floor(score)))
  text(50, 490, paste0('R2 = ', floor(r2*100)/100))
  
  dev.off()
  
  
  png('../Images/Full Grid Search Status.png')
  
  plot(grid$A, grid$b, col = make.seq.pal(grid$score, n = 16, min = 0, max = 30, rev = TRUE, pal = 'ocean.ice'), pch=15, cex=1.2,
       xlim = c(7.3, 20), ylim = c(1, 4.8), xlab = 'A', ylab = 'b', main = 'Misfit')
  points(grid$A, grid$b, pch = 0, cex = 1.2, col = 'grey')
  text(grid$A, grid$b, c(1:nrow(grid)), cex= 0.3)
  
  ## Colorbar
  add.colorbar(0, 30, 9, 1.7, 9, 3.7, pal = 'ocean.ice', rev = TRUE, n = 16)
  
  ## Reference values
  points(12.1, 3.81, cex=2, pch = 16)
  points(grid$A[which.min(grid$score)], grid$b[which.min(grid$score)], pch=2, cex=3)
  
  
  dev.off()
}

#### Summary Figure

pdf('../Images/Full Grid Search Example.pdf')

plot(grid$A, grid$b, col = make.seq.pal(grid$score, n = 16, min = 10, max = 50, rev = TRUE), pch=15, cex=1.2,
     xlim = c(12, 14), ylim = c(1.4, 1.8), xlab = 'A', ylab = 'b', main = 'Misfit')
points(grid$A, grid$b, pch = 0, cex = 1.2, col = 'grey')
add.colorbar(0, 300, 9, 1.7, 9, 3.7, pal = 'YlOrRd', rev = TRUE, n = 16)
points(12.1, 3.81, cex=2, pch = 16)
points(grid$A[which.min(grid$score)], grid$b[which.min(grid$score)], pch=2, cex=3)
text(grid$A, grid$b, c(1:nrow(grid)), cex= 0.3)


plot(grid$A, grid$b, col = make.seq.pal(grid$r2, n = 16, min = 0, max = 1, rev = TRUE), pch=15, cex=1.2,
     xlim = c(7.3, 20), ylim = c(1, 4.8), xlab = 'A', ylab = 'b', main = 'R2')
points(grid$A, grid$b, pch = 0, cex = 1.2, col = 'grey')
add.colorbar(0, 1, 9, 1.7, 9, 3.7, pal = 'YlOrRd', rev = TRUE, n = 16)
points(12.1, 3.81, cex=2, pch = 16)
points(grid$A[which.min(grid$score)], grid$b[which.min(grid$score)], pch=2, cex=3)
points(grid$A[which.max(grid$r2)], grid$b[which.max(grid$r2)], pch=2, cex=3)
text(grid$A, grid$b, c(1:nrow(grid)), cex= 0.3)

dev.off()
