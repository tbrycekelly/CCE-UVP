## Set of useful R functions for general use in plotting, analyzing and
## converting.
##
## Author: Thomas Bryce Kelly (tbk14 at fsu.edu)
## https://urldefense.proofpoint.com/v2/url?u=http-3A__about.tkelly.org_&d=DwIGAg&c=MNHwOqQ8N1u91SoMLfIblwuGXKgp50OPUXjl8uRAbak&r=OqmwiruNRrBVjwRDz_rOIW4lt52uD5YZ94BCkXBkcSU&m=A8ifeaSsB7OUd2OIZ_RVFtX7Y6rtEsOOC77zcbZ0C6s&s=-_9Yx_KpWVvj1uMgQs3TVGvsru9RvaB7BSu7mj9usJQ&e=
##
## Dept of Earth, Ocean & Atmospherical Sciences
## Florida State University
##
## Center for Ocean & Atmospheric Prediction Studies
## Florida State University
##
## National High Magnetic Field Laboratory
## Florida State University


##############################
## Date Times ################
##############################

## Helper function for converting the date time stamps.
conv.time.excel = function (x, tz = "GMT") {
    as.POSIXct(x*86400, origin = "1899-12-30 00:00:00", tz = tz)
}

conv.time.unix = function(x, tz='GMT') {
    as.POSIXct(x, origin="1970-01-01", tz=tz)
}

conv.time.matlab = function (x, tz = "GMT") {
    as.POSIXct((x - 1)*86400, origin = "0000-01-01", tz = tz)
}

## Find the indicies of the closest times for each entry of x
which.closest.time = function(x, y) {
    if (length(y) > 1) {
        l = c()
        for (i in 1:length(x)) {
            l.new = which.min(as.numeric(difftime(x[i], y, units='mins'))^2)
            l = c(l, l.new)
        }
    } else {
        l = which.min(as.numeric(difftime(x, y, units='mins'))^2)
    }
    l
}


which.unique = function(x) {
    which(!duplicated(x))
}    

##############################
## Plotting ################
##############################

get.pal = function(n=10, pal='ocean.haline', rev = FALSE) {
    pal = do.call(pal, list(n = n))
        
    if (rev) {
        pal = rev(pal)
    }
    pal
}

        
make.pal = function(x, n = 10, min = NA, max = NA, pal='ocean.haline', rev = FALSE) {
    cols = get.pal(n+1, pal = pal, rev = rev)
    
    if (is.na(min)) {  ## set a minimum
        min = base::min(x, na.rm=TRUE)
    }
    if (is.na(max)) { ## Set a maximum
        max = base::max(x, na.rm=TRUE)
    }
    
    ## Force min and max
    x[x < min] = min
    x[x > max] = max
    
    x = (x-min) * n / (max-min) ## Scale so x falls between [0 -> n]  
    cols[floor(x)+1] # return the colors
}

make.pal.custom = function(x, n = 10, pal = c('white', 'blue'),  min = NA, max = NA, rev = FALSE) {
    cols = colorRampPalette(pal)(n)
    
    if (is.na(min)) {  ## set a minimum
        min = base::min(x, na.rm=TRUE)
    }
    if (is.na(max)) { ## Set a maximum
        max = base::max(x, na.rm=TRUE)
    }
    
    ## Force min and max
    x[x < min] = min
    x[x > max] = max
    
    x = (x-min) * n / (max-min) ## Scale so x falls between [0 -> n]  
    cols[floor(x)+1] # return the colors
    
}

## make div and seq for numeric values
make.div.pal = function(x, n = 10, min = NA, max = NA, pal = 'coolwarm', rev = FALSE) {
    make.pal(x, n, min, max, pal, rev)
}

make.seq.pal = function(x, n = 10, min = NA, max = NA, pal = 'ocean.haline', rev = FALSE) {
    make.pal(x, n, min, max, pal, rev)
}


## Make pal for categorical data
make.qual.pal = function(x=100, pal='tol', rev = FALSE) {
    
    ## Determine numeric values
    a = sapply(x, function(xx) {which(xx == unique(x))})
    
    cols = get.pal(n = length(unique(x)), pal = pal, rev = rev)
    cols[a]
}


add.colorbar = function(min, max, X1, Y1, X2, Y2, pch=15, cex=1, cex.text=1, col.bg=NA,  n=10,
                        pal='ocean.algae', text.dist = 0.15, rev = FALSE, log = FALSE, base = 10,
                        n.tick = 5, cex.tick = 1, col.tick = 'black') {
    
    x = seq(X1, X2, length.out = 200)
    y = seq(Y1, Y2, length.out = 200)
    
    if (!is.na(col.bg)) { 
        points(x, y, col = col.bg, cex = cex + 0.2, pch=pch)
    }
    
    col = make.pal(x = c(1:200), n = n, pal = pal, rev = rev)
    
    if (log) {
        k = c(1:9)
        k = c(k/1000, k/100, k/10, k, 10*k, 100 *k , 1000 * k, 1e4*k)
        k = k[which( k >= min & k <= max)]
        k = log(k, base)
        print(k)

        k.min = log(min, base)
        k.max = log(max, base)

        k.x = approx(c(1:length(x)), x, xout = length(x) * (k - k.min) / (k.max - k.min))$y
        k.y = approx(c(1:length(x)), y, xout = length(x) * (k - k.min) / (k.max - k.min))$y

        points(k.x, k.y, pch = 3, cex = cex.tick, col = col.tick)
        
    } else if (!is.na(n.tick)) {
        k = c(0:n.tick)
        
        k.x = approx(c(1:length(x)), x, xout = length(x) * k / n.tick, rule = 2)$y
        k.y = approx(c(1:length(x)), y, xout = length(x) * k / n.tick, rule = 2)$y
        
        points(k.x, k.y, pch = 3, cex = cex.tick, col = col.tick)
    }
    
    
    ## Color bar points
    points(x, y, cex = cex, pch=pch, col = col)
    
    text(X1 - (X2 - X1) * text.dist, Y1 - (Y2 - Y1) * text.dist, min, cex = cex.text)
    text(X2 + (X2 - X1) * text.dist, Y2 + (Y2 - Y1) * text.dist, max, cex = cex.text)
}


add.colorbar.log = function(min, max, X1, Y1, X2, Y2, base, pch=15, cex=1, cex.text=1, col.bg=NA,  n=10,
                        pal='ocean.algae', text.dist = 0.15, rev = FALSE, n.tick = 5, cex.tick = 1, col.tick = 'black') {
  
  x = seq(X1, X2, length.out = 200)
  y = seq(Y1, Y2, length.out = 200)
  
  if (!is.na(col.bg)) { 
    points(x, y, col = col.bg, cex = cex + 0.2, pch=pch)
  }
    
  col = make.pal(x = base^seq(log(min, base), log(max, base), length.out = 200),
                 n = n, pal = pal, rev = rev)
    
    if (!is.na(n.tick)) {
        k = c(0:n.tick)
        
        k.x = approx(c(1:length(x)), x, xout = length(x) * k / n.tick, rule = 2)$y
        k.y = approx(c(1:length(x)), y, xout = length(x) * k / n.tick, rule = 2)$y
        
        points(k.x, k.y, pch = 3, cex = cex.tick, col = col.tick)
    }
    
    
  points(x, y, col = col, cex = cex, pch=pch)

    
  
  text(X1 - (X2 - X1) * text.dist, Y1 - (Y2 - Y1) * text.dist, min, cex = cex.text)
  text(X2 + (X2 - X1) * text.dist, Y2 + (Y2 - Y1) * text.dist, max, cex = cex.text)
}

##############################
## Mapping (cartesian) #######
##############################

plot.map = function(lon, lat, main = '', xlab='Longitude', ylab='Latitude', col='black', pch=20,
                    cex=1, zoom = 0, lat.offset = 0, lon.offset = 0, bounds = NULL) {
    
    if (is.null(bounds)) {
        lat.max = max(lat, na.rm = TRUE) + zoom + lat.offset
        lat.min = min(lat, na.rm = TRUE) - zoom + lat.offset
        lon.max = max(lon, na.rm = TRUE) + zoom + lon.offset
        lon.min = min(lon, na.rm = TRUE) - zoom + lon.offset
    } else{
        lat.min = bounds[3]
        lat.max = bounds[4]
        lon.min = bounds[1]
        lon.max = bounds[2]
    }
    newmap <- getMap(resolution = "high")
    
    plot(newmap, xlim = c(lon.min, lon.max), ylim = c(lat.min, lat.max), asp = 1)
    points(lon, lat, main = main, col = col, xlab = xlab, ylab = ylab, pch = pch, cex = cex)
}

add.map.grid = function(by = 1, axis.by = 2, col='#00000020') {
    
    for (i in seq(-180, 180, by=by)) {
        lines(rep(i,2), c(0,60), lty=2, col=col)
    }

    for (i in seq(-90, 90, by=by)) {
        lines(c(-200,-100), rep(i,2), lty=2, col=col)
    }
    
    axis(side = 2, at = seq(10, 60, by = axis.by))
    axis(side = 1, at = seq(-200, -100, by = axis.by))
    axis(side = 4, at = seq(10, 60, by = axis.by))
    axis(side = 3, at = seq(-200, -100, by = axis.by))
}


##############################
## Statistics ################
##############################

jackknife = function(x, y, n=1000) {
    
    k = c(which(is.na(x)), which(is.na(y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        y = y[-k]
    }
    
    res = data.frame(m=0, b=0)
    
    for (i in 1:n) {
        l = sample(1:length(x), replace = TRUE)
        temp.model = lm(y[l] ~ x[l])
        res = rbind(res, c(coef(temp.model)[2], coef(temp.model)[1]))
    }
    res= res[-1,]
    res
}

bootstrap = function(x, s.x, y, s.y, n = 1000) {
    
    k = c(which(is.na(x)), which(is.na(s.x)), which(is.na(y)), which(is.na(s.y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        s.x = s.x[-k]
        y = y[-k]
        s.y = s.y[-k]
    }
    
    result = data.frame(m = 0, b = 0)
    num.x = length(x)
    
    for (i in 1:n) {
        l = sample(x = c(1:num.x), replace = TRUE)
        
        ## Generate random sampling
        temp.x = rnorm(num.x, x[l], s.x[l])
        temp.y = rnorm(num.x, y[l], s.y[l])
        
        model = lm(temp.y ~ temp.x)
        
        result = rbind(result, rev(coefficients(model)))
    }
    result = result[-1,]
    result
}

bootstrap.2 = function(x, s.x, y, s.y, n = 1000) {
  
    k = c(which(is.na(x)), which(is.na(s.x)), which(is.na(y)), which(is.na(s.y)))
    if (length(k) > 0) {
        k = unique(k)
        x = x[-k]
        s.x = s.x[-k]
        y = y[-k]
        s.y = s.y[-k]
    }
    
    result = data.frame(m = 0, b = 0)
    num.x = length(x)
    
    for (i in 1:n) {
        l = sample(x = c(1:num.x), size = num.x, replace = TRUE)
        
        ## Generate random sampling
        temp.x = rnorm(num.x, x[l], s.x[l])
        temp.y = rnorm(num.x, y[l], s.y[l])
        
        model = lmodel2(temp.y ~ temp.x, nperm = 0)
        
        result = rbind(result, as.numeric(model$regression.results[3,c(3,2)]))
    }
    result = result[-1,]
    result
}

bootstrap.lmodel2.both = function(x, s.x, y, s.y, n=100) {
    mod = lmodel2(y ~ x)$regression.results
    results.OLS = mod[1,c(2:3)]
    results.MA = mod[2,c(2:3)]
    
    for (i in 1:n) {
        new.x = rnorm(length(x), x, s.x)
        new.y = rnorm(length(y), y, s.y)
        
        mod = lmodel2(new.y ~ new.x)$regression.results
        
        results.OLS = rbind(results.OLS, mod[1,c(2:3)])
        results.MA = rbind(results.MA, mod[2,c(2:3)])
    }
    list(OLS = results.OLS, MA = results.MA)
}


#### CONFIDENCE AND TRENDLINES

add.boot.trendline.median = function(model, new.x, col = 'black', lty = 2) {
    m = median(model$m, rm.na = TRUE)
    b = median(model$b, na.rm = TRUE)
    lines(new.x, new.x * m + b, col = col, lty = lty)
}

add.boot.trendline = function(model, new.x, col = 'black', lty = 2, lwd=1) {
    res = c()
    for (i in 1:length(new.x)) {
        res = c(res, median(model$m * new.x[i] + model$b, na.rm=TRUE))
    }
    lines(new.x, res, col = col, lty = lty, lwd = lwd)
}


add.boot.trendline2 = function(res, new.x, col = 'black', lty = 2) {
    m = mean(res[,2], rm.na = TRUE)
    b = mean(res[,1], na.rm = TRUE)
    lines(new.x, new.x * m + b, col = col, lty = lty)
}

add.lm.conf = function(x, name, model, col = '#50505030', level = 0.95, log = FALSE) {
    if(!log) {
        dat = data.frame(a = c(1:length(x)))
        dat[[name]] = x
        pred = predict(model, interval='confidence', newdata = dat, level = level)
        polygon(c(x,rev(x)), c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col)
    } else {
        dat = data.frame(a = c(1:length(x)))
        dat[[name]] = x
        pred = predict(model, interval='confidence', newdata = dat, level = level)
        polygon(c(exp(x),rev(exp(x))), c(pred[,"lwr"], rev(pred[,"upr"])), border = NA, col = col)
    }
}

add.boot.conf = function(model, new.x = c(-1e3:1e3), col = '#55555540', conf = c(0.025, 0.975),
                         border = FALSE, trendline = FALSE) {
    x.upper = new.x
    y.upper = c()
    y.lower = c()
    y.mid = c()

    for (x in x.upper) {
        y = x * model$m + model$b
        y.mid = c(y.mid, median(y))
        
        ## Quantiles
        y = quantile(y, conf)
        y.lower = c(y.lower, y[1])
        y.upper = c(y.upper, y[2])
    }
    
    polygon(x = c(x.upper, rev(x.upper)), y = c(y.upper, rev(y.lower)), col = col, border = border)
    
    if (trendline) {
        lines(x.upper, y.mid, lty=2)
    }
}

get.boot.vals = function(model, x, conf = 0.5) {
    y = rep(0, length(x))
    
    for (i in 1:length(y)) {
        yy = x[i] * model$m + model$b
        y[i] = quantile(yy, probs = conf, na.rm = TRUE)[[1]]
    }
    
    y
}

add.boot.trendline.sig = function(model, new.x, p = 0.01, sig = 0, lty = 2, lwd = 1, col = 'black') {
    if(ecdf(model$m)(sig) < p | ecdf(model$m)(sig) > 1 - p) {
        add.boot.trendline(model, new.x, col = col, lty = lty, lwd = lwd)
    }
}

####################
#### Statistics ####
####################

#### Calc R Squared
calc.boot.r.squared = function(model, x, y) {
    ## rm NAs
    l = which(is.na(x) | is.na(y))
    if (length(l) > 0) {
        x = x[-l]
        y = y[-l]
    }
    
    ## R2 = 1 - SSR / SS
    1 - boot.ssr(model = model, x = x, y = y) / calc.boot.ss(y = y)
}

#### Sum of Squared Residuals
calc.boot.ssr = function(model, x, y) {
    ## rm NAs
    l = which(is.na(x) | is.na(y))
    if (length(l) > 0) {
        x = x[-l]
        y = y[-l]
    }
    
    ## SSR
    y.pred = rep(0, length(x))
    
    for (i in 1:length(x)) {
        y.pred[i] = median(model$m) * x[i] + median(model$b)
    }
    
    sum((y.pred - y)^2)
}

#### Sum of Squares
calc.boot.ss = function(y) {
    y = y[!is.na(y)] ## Remove NAs
    
    ## Calc Total Sum of Squares
    sum((mean(y) - y)^2)
}


calc.r.squared = function(x, y){
    1 - sum((x - y)^2) / sum((x - mean(x))^2)
}

########################
#### Misc Plot Stuff ###
########################

add.error.bars = function(x, s.x, y, s.y, col = 'black') {
    
    if (length(s.x) == 1) {
        s.x = rep(s.x, length(x))
    }
    if (length(s.y) == 1) {
        s.y = rep(s.y, length(y))
    }
    
    for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]),
              y = c(y[i] + s.y[i], y[i] - s.y[i]),
              col = col)
        
        lines(x = c(x[i] + s.x[i], x[i] - s.x[i]),
              y = c(y[i], y[i]),
              col = col)
    }
}

add.error.bars.logx = function(x, s.x, y, s.y, base, col = 'black') {
  
  if (length(s.x) == 1) {
    s.x = rep(s.x, length(x))
  }
  if (length(s.y) == 1) {
    s.y = rep(s.y, length(y))
  }
  
  ## Remove NAs
  l = !is.na(x) & !is.na(y)
  x = x[l]
  s.x = s.x[l]
  y = y[l]
  s.y = s.y[l]
  
  for (i in 1:length(x)) {
    if(!is.na(s.y[i])) {
      lines(x = rep(log(x[i], base), 2),
            y = c(y[i] + s.y[i], y[i] - s.y[i]),
            col = col)
    }
    
    if(!is.na(s.x[i])) {
      lines(x = log(c(x[i] + s.x[i], x[i] - s.x[i]), base),
            y = rep(y[i], 2),
            col = col)
    }
  }
}

add.error.bars.logy = function(x, s.x, y, s.y, base, col = 'black') {
  
  if (length(s.x) == 1) {
    s.x = rep(s.x, length(x))
  }
  if (length(s.y) == 1) {
    s.y = rep(s.y, length(y))
  }
  
  ## Remove NAs
  l = !is.na(x) & !is.na(y)
  x = x[l]
  s.x = s.x[l]
  y = y[l]
  s.y = s.y[l]
  
  for (i in 1:length(x)) {
    if(!is.na(s.y[i])) {
      lines(x = rep(x[i], 2),
            y = log(c(y[i] + s.y[i], y[i] - s.y[i]), base),
            col = col)
    }
    
    if(!is.na(s.x[i])) {
      lines(x = c(x[i] + s.x[i], x[i] - s.x[i]),
            y = rep(log(y[i], base), 2),
            col = col)
    }
  }
}

add.one.to.one = function(scale = 1e7, lty=2, col='#00000080', slope = 1) {
    lines(c(-scale, scale), c(-scale*slope, scale*slope), lty=lty, col=col)
}



## Install the package if needed:
#install.packages('cmocean', repos='https://urldefense.proofpoint.com/v2/url?u=http-3A__cran.us.r-2Dproject.org&d=DwIGAg&c=MNHwOqQ8N1u91SoMLfIblwuGXKgp50OPUXjl8uRAbak&r=OqmwiruNRrBVjwRDz_rOIW4lt52uD5YZ94BCkXBkcSU&m=A8ifeaSsB7OUd2OIZ_RVFtX7Y6rtEsOOC77zcbZ0C6s&s=HV6H_jB_gohewcGJehHluAOL9hXKAyWQt63wTOeePk0&e=')

## Load the relevent packages
library(ncdf4)  # For reading in the NCEP wind fields
#library(R.matlab)  # If you need to read in matlab .mat files
library(openxlsx)  # If you need to read in .xlsx files
library(RColorBrewer)
library(compiler)  # required for JIT (below)
library(rworldmap)
library(rworldxtra)
library(rgdal)
library(lmodel2)
library(pals)

## Enable compilation (speed gain?)
enableJIT(3)
