setwd('C:/GURUNG/GEOG639/PROJECT/')

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















# Perform MRA on image 
# --------------------------
rm(wt, wt.bw, wt.sss, wt.var1)
p <- 1 
q <- n_row                                                             # For west-to-east,  use n_row; 
                                                                       #for north-to-south, use n_col;

l1 <- matrix(nrow = n_row, ncol= n_col)
l2 <- matrix(nrow = n_row, ncol= n_col)
l3 <- matrix(nrow = n_row, ncol= n_col)
l4 <- matrix(nrow = n_row, ncol= n_col)
l5 <- matrix(nrow = n_row, ncol= n_col)
l6 <- matrix(nrow = n_row, ncol= n_col)
l7 <- matrix(nrow = n_row, ncol= n_col)
l8 <- matrix(nrow = n_row, ncol= n_col)
for(i in p:q){
  i_start = 1 
  i_end   = n_col 
  
  tline <- t[i,i_start:i_end] # For West-East
  
  #print(paste("max", max(tline), "min", min(tline)))
  #tline <- as.data.frame(tline)
  wt <- modwt(tline, n.levels = 8, wf)
  wt.bw <- brick.wall(wt, wf)
  wt.align <- phase.shift(wt.bw, wf, inv = FALSE)
  # Level 1
  l1[i,1:n_col] <- wt.align[[1]]  
  l2[i,1:n_col] <- wt.align[[2]]  
  l3[i,1:n_col] <- wt.align[[3]]  
  l4[i,1:n_col] <- wt.align[[4]]  
  l5[i,1:n_col] <- wt.align[[5]]  
  l6[i,1:n_col] <- wt.align[[6]]  
  l7[i,1:n_col] <- wt.align[[7]]  
  l8[i,1:n_col] <- wt.align[[8]]   
  #
  # l1[i,1:n_col] <- wt.align@W[[1]]
  # l2[i,1:n_col] <- wt.align@W[[2]]
  # l3[i,1:n_col] <- wt.align@W[[3]]
  # l4[i,1:n_col] <- wt.align@W[[4]]
  # l5[i,1:n_col] <- wt.align@W[[5]]
  # l6[i,1:n_col] <- wt.align@W[[6]]
  # l7[i,1:n_col] <- wt.align@W[[7]]
  # l8[i,1:n_col] <- wt.align@W[[8]]  
}

str_name <- "IMAGES/final_"

writeTiff(l1, paste(str_name, as.character(1), ".tif", sep = ""))
writeTiff(l2, paste(str_name, as.character(2), ".tif", sep = ""))
writeTiff(l3, paste(str_name, as.character(3), ".tif", sep = ""))
writeTiff(l4, paste(str_name, as.character(4), ".tif", sep = ""))
writeTiff(l5, paste(str_name, as.character(5), ".tif", sep = ""))
writeTiff(l6, paste(str_name, as.character(6), ".tif", sep = ""))
writeTiff(l7, paste(str_name, as.character(7), ".tif", sep = ""))
writeTiff(l8, paste(str_name, as.character(8), ".tif", sep = ""))

rm(l1, l2, l3, l4, l5, l6 ,l7, l8)


