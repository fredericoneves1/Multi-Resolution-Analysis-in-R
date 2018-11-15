setwd('C:/GURUNG/')

# Load Libraries
# --------------
#library('wavelets')
library('tiff')
library('rgdal')
library('raster')
library('rtiff')
library('waveslim')

##########################################
#      calculate Albedo                  #
##########################################

# Landsat 4 MTL
# -------------
deg2rad <- function(deg) {(deg * pi) / (180)}
landsat4 <- read.csv('DATA/landsat4.csv', header = TRUE)
jd_ls4 <-106 # Julian day
d_ls4 <- 1.0033 
esun_ls4 <- c(1958, 1826, 1554, 1033, 214.7, 0, 80.7)
sun_elev_ls4 <- 27.68842959
zenith_ls4 <- deg2rad(90 - sun_elev_ls4)
bands <- c(1,3,4,5,7)

# Convert to TOA
# --------------
p <- list()
for(i in 1:5){
  k = bands[i]
  str_name <- paste('DATA/SUBSET_B' , as.character(k) , '.tif', sep='')
  t = readTIFF(str_name, as.is=TRUE)
  print(paste("Reading", as.character(str_name)))
  # Convert DN to Radiance
  t <- landsat4$G_RESCALE[k]*t + landsat4$B_RESCALE[k]
  # Convert Radiance to TOA
  t <- (pi*t*d_ls4^2)/(esun_ls4[k]*cos(zenith_ls4))
  # Save t
  p[[i]] <- t
  print(paste("Writing p", as.character(i), sep = ''))
}
rm(t)
#
# Convert to Albedo
# -----------------
alpha <- (0.356*p[[1]] + 0.130*p[[2]] + 0.373*p[[3]] + 0.085*p[[4]] + 0.072*p[[5]] - 0.0018)/1.016

# Save Albedo as georeferenced TIFF
# ---------------------------------
t <- raster('DATA/SUBSET_B1.tif')
study_extent <- t@extent
study_coord  <- t@crs
d <- raster(alpha)
d@extent <- study_extent
d@crs <- study_coord
set.seed(1)
d <- ratify(d)
writeRaster(d, 'DATA/albedo.tif', datatype='INT1U', overwrite=TRUE)
rm(t,d)
#
# Perform MRA on Band 3
# ---------------------
t <- alpha
#rm(alpha)
n_row <-nrow(t)
n_col <-ncol(t)
print(paste("Dimension of image are: ", as.character(n_row), " x ", 
            as.character(n_col), sep =""))
#
# Perform 1-D MRA 
# ---------------
# Choose nrow = 1918 as transect
# Sampling interval = 30 m, same as resolution
#
wf <- "d4"
ltest <- matrix(nrow = 1, ncol= n_col)
tline <- t[1918,1:n_col]                                             # For West-East Transect
# tline <- as.data.frame(tline)                                      # Only needeD for 'wavelets' package
wt <- modwt(tline, n.levels = 8, wf)
wt.bw <- brick.wall(wt.bw, wf)
wt.bw <- phase.shift(wt, wf, inv = FALSE)
wt.var1 <- wave.variance(wt.bw, type="eta3", p=0.025)
plot(wt.var1$upper[1:8], type ='l', ylim = c(0,0.0003))
lines(wt.var1$wavevar[1:8], type ='l', ylim = c(0,0.0003))
lines(wt.var1$lower[1:8], type ='l', ylim = c(0,0.0003))

t_alp <- t
# cumulative sum of Squares
# -------------------------
m <- length(wt.bw) - 1
n <- length(tline)
wt.sss <- wt.bw[1:8]
for(i in 1:m){
  a <- wt.bw[[i]]^2
  asum <- a
  for(j in 2:n){
    test <- asum[j-1] + a[j]
    if(is.na(test) == FALSE){
      asum[j] <- asum[j-1] + a[j]
    }
  }
  wt.sss[[i]] <- asum
}

# Plot
# ----
par(mfrow=c(8,1))
par(mar = rep(2, 4))
for(i in 8:1){
  plot.ts(wt.sss[[i]])
}
par(mfrow=c(1,1))

##
#
#
##########################################
#       Calculate BT                     #
##########################################
#
# Read Thermal band
# -----------------
# Convert to Brightness
# --------------
k <- 6 # Band number
K1 <- 671.62 # K1_CONSTANT_BAND_6
K2 <- 1284.30 #K2_CONSTANT_BAND_6 
str_name <- paste('DATA/SUBSET_B' , as.character(k) , '.tif', sep='')
t = readTIFF(str_name, as.is=TRUE)
print(paste("Reading", as.character(str_name)))
# Convert DN to Radiance
t <- landsat4$G_RESCALE[k]*t + landsat4$B_RESCALE[k]
# Convert Radiance to TOA
t <- K2/log( (K1/t) + 1)
#
# Perform MRA on Thermal band for 30 m spacing
# --------------------------------------------
# Choose nrow = 1918 as transect
# Sampling interval = 30 m, same as resolution
#
wf <- "d4"
ltest <- matrix(nrow = 1, ncol= n_col)
tline <- t[1918,1:n_col]                                             # For West-East Transect
# tline <- as.data.frame(tline)                                      # Only needeD for 'wavelets' package
wt <- modwt(tline, n.levels = 8, wf)
wt.bw <- brick.wall(wt.bw, wf)
wt.bw <- phase.shift(wt, wf, inv = FALSE)
wt.var1 <- wave.variance(wt.bw, type="eta3", p=0.025)

minVar <- min(wt.var1$wavevar[1:8])
diff <- max(wt.var1$wavevar[1:8]) - minVar
wt.var2 <- (wt.var1 - minVar)/diff


plot(wt.var2$upper[1:8], type ='l' )#, ylim = c(0,0.35))
lines(wt.var2$wavevar[1:8], type ='l')#, ylim = c(0,0.35))
lines(wt.var2$lower[1:8], type ='l')#, ylim = c(0,0.35))


# cumulative sum of Squares
# -------------------------
m <- length(wt.bw) - 1
n <- length(tline)
wt.sss <- wt.bw[1:8]
for(i in 1:m){
  a <- wt.bw[[i]]^2
  asum <- a
  for(j in 2:n){
    test <- asum[j-1] + a[j]
    if(is.na(test) == FALSE){
      asum[j] <- asum[j-1] + a[j]
    }
  }
  wt.sss[[i]] <- asum
}

# Plot
# ----
par(mfrow=c(8,1))
par(mar = rep(2, 4))
for(i in 8:1){
  plot.ts(wt.sss[[i]])
}
par(mfrow=c(1,1))

plot(results[1:8,1], type='l', ylim = c(-0.01,0))
plot(results[1:8,2], type='l', ylim = c(-0.01,0))
plot(results[1:8,3], type='l', ylim = c(-0.01,0))




