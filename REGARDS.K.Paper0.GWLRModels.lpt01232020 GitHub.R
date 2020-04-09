rm(list = ls(all = T))

#########################################
#########################################  libraries + functions
#########################################
#########################################
library(spdep)
library(maptools)
library(lattice)
library(rgdal)
library(GWmodel)
library(plyr)
library(RColorBrewer)
library(spgwr)
library(tigris)
library(GISTools)

#########################################
#########################################  individual level data
#########################################
#########################################
# importing final analytic dataset regardsmerge.csv
regardsmerge <- read.csv("/Volumes/REGARDS/Tabb/Paper 0/Data/regardsmerge1.0.csv") # 17889 REGARDS black and white participants from 3 regions

impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

regardsmerge <- ddply(regardsmerge, ~ REGION.x, transform, 
                      Tract_Urban_Pct = impute.mean(Tract_Urban_Pct))

#MERGE SHAPE AND REGARDS FILE
#TABLE OF REGARDS STATE VARIABLE
#MERGE SHAPE FILES WITH DF 
#SPCBIND BY STATE ABBREV.

n.final <- dim(regardsmerge)[1]

#####################################
##################################### stratified datasets
#####################################
# black and white only REGARDS dataset
regardsmerge_black <- regardsmerge[regardsmerge$Race == "B",]
regardsmerge_white <- regardsmerge[regardsmerge$Race == "W",]
n.final.black <- dim(regardsmerge_black)[1]
n.final.white <- dim(regardsmerge_white)[1]



#####################################
##################################### spatial data frames 
#####################################
#create dp.locat
regardsmergeALL.coord0 <- subset(regardsmerge, select = c(X,Y))

#removing observations with missing coordinates
regardsmergeALL.coord <- na.omit(regardsmergeALL.coord0) #0 observations removed

#SpatialPointsDataFrame
# no projection set here - ie, used default settings
regardsmergeALL.sp <- SpatialPointsDataFrame(regardsmergeALL.coord,regardsmerge)

# stratified SpatialPointsDataFrame
# blacks
#create dp.locat
regardsmergeALL.coord0.black <- subset(regardsmerge_black, select = c(X,Y))

#removing observations with missing coordinates
regardsmergeALL.coord.black <- na.omit(regardsmergeALL.coord0.black) #0 observations removed

#SpatialPointsDataFrame
# no projection set here - ie, used default settings
regardsmergeALL.sp.black <- SpatialPointsDataFrame(regardsmergeALL.coord.black,regardsmerge_black)

# whites
#create dp.locat
regardsmergeALL.coord0.white <- subset(regardsmerge_white, select = c(X,Y))

#removing observations with missing coordinates
regardsmergeALL.coord.white <- na.omit(regardsmergeALL.coord0.white) #0 observations removed

#SpatialPointsDataFrame
# no projection set here - ie, used default settings
regardsmergeALL.sp.white <- SpatialPointsDataFrame(regardsmergeALL.coord.white,regardsmerge_white)


#####################################
##################################### OLS Models
#####################################
ALL.model1 <- lm(scale(regardsmerge$sum7) ~  relevel(factor(regardsmerge$Race), ref = "W") + scale(regardsmerge$Age) + factor(regardsmerge$Gender))
summary(ALL.model1)

regardsmerge$REGION <- as.factor(regardsmerge$REGION)
ALL.model2 <- lm(scale(regardsmerge$sum7) ~  relevel(factor(regardsmerge$REGION), ref = "0") + relevel(factor(regardsmerge$Race), ref = "W") + scale(regardsmerge$Age) + factor(regardsmerge$Gender))
summary(ALL.model2)

#xy.BELT <- coordinates(regardsmergeBELT.sp)

#lm.morantest(BELT.model1, nb2listw(regardsmergeBELT, style = "W"))

# stratified analyses
ALL.model1.black <- lm(scale(regardsmerge_black$sum7) ~  scale(regardsmerge_black$Age) + factor(regardsmerge_black$Gender))
summary(ALL.model1.black)
ALL.model1.white <- lm(scale(regardsmerge_white$sum7) ~  scale(regardsmerge_white$Age) + factor(regardsmerge_white$Gender))
summary(ALL.model1.white)


regardsmerge_black$REGION <- as.factor(regardsmerge_black$REGION)
regardsmerge_white$REGION <- as.factor(regardsmerge_white$REGION)

ALL.model2.black <- lm(scale(regardsmerge_black$sum7) ~  relevel(factor(regardsmerge_black$REGION), ref = "0") + scale(regardsmerge_black$Age) + factor(regardsmerge_black$Gender))
summary(ALL.model2.black)
ALL.model2.white <- lm(scale(regardsmerge_white$sum7) ~  relevel(factor(regardsmerge_white$REGION), ref = "0") + scale(regardsmerge_white$Age) + factor(regardsmerge_white$Gender))
summary(ALL.model2.white)


#####################################
##################################### GWR: Fixed and Adaptive Models
#####################################

plot(regardsmergeALL.sp, cex=0.3, pch=20) 
plot(regardsmergeALL.sp.black, cex=0.3, pch=20) 
plot(regardsmergeALL.sp.white, cex=0.3, pch=20) 


###############
###############
#---ALL---#
###############
###############
# model 1
#gwr formulas
formula.m1.ALL.gwr <- scale(regardsmerge$sum7) ~  relevel(factor(regardsmerge$Race), ref = "W") + scale(regardsmerge$Age) + factor(regardsmerge$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m1 <- gwr.sel(formula = formula.m1.ALL.gwr, data=regardsmergeALL.sp, adapt = T)
Sys.time()
# 1 hour and 11 mins
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m1 <- gwr(formula = formula.m1.ALL.gwr, adapt = bw.ALL.m1, data = regardsmergeALL.sp, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 1 day and 4 hours (laptop: 9.39am day 1 - 21.02pm day 2) - old stats
# 1 day and 13 hours (day one: 15.51 - day two: 5.30)
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")
# 8 mins to save image

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m1 <- BFC99.gwr.test(gwr.ALL.m1)
Sys.time()
# 3 hours and 15 mins. - old stats
# 3.5 hours (5.38am - 9.04am)
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")
# 8 mins to save image

# model 2
#gwr formulas
formula.m2.ALL.gwr <- scale(regardsmerge$sum7) ~  relevel(factor(regardsmerge$REGION), ref = "0") + relevel(factor(regardsmerge$Race), ref = "W") + scale(regardsmerge$Age) + factor(regardsmerge$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m2 <- gwr.sel(formula = formula.m2.ALL.gwr, data=regardsmergeALL.sp, adapt = T)
Sys.time()
# 1 hour and 30 mins - old stats
# 2 hours (9.12am - 11.11am)
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m2 <- gwr(formula = formula.m2.ALL.gwr, adapt = bw.ALL.m2, data = regardsmergeALL.sp, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 1 day and 13 hours (day one: 18.47 - day two: 7.10)
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")
# 19 mins

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m2 <- BFC99.gwr.test(gwr.ALL.m2)
Sys.time()
# 3 hours and 20 mins
save.image("/Volumes/REGARDS/Tabb/Paper 0/08192019RESULTS.RData")
# 16 mins

# all data stored up to this point, as of 8.23.2019

# Moran's I for GWR residuals
knearneigh.regards <- knearneigh(regardsmergeALL.sp)
knn2nb.regards <- knn2nb(knearneigh.regards)
nb2listw.regards <- nb2listw(knn2nb.regards)
#Sys.time()
#gwr.morantest.m1 <- gwr.morantest(x = gwr.ALL.m1, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
#Sys.time()
#gwr.morantest.m2 <- gwr.morantest(x = gwr.ALL.m2, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
# memory issues are faced when running gwr.morantest(); "Error: vector memory exhausted (limit reached?)"


# check the residuals for autocorrelation
set.seed(1234)
moran.m1 <- moran.test(gwr.ALL.m1$SDF$gwr.e, nb2listw.regards, alternative = "greater", zero.policy = T) 
# Moran I statistic standard deviate = 1.7074, p-value = 0.04387
# Moran I statistic       Expectation          Variance 
# 1.604081e-02            -5.590340e-05        8.887843e-05 

moran.m2 <- moran.test(gwr.ALL.m2$SDF$gwr.e, nb2listw.regards, alternative = "greater", zero.policy = T)
# Moran I statistic standard deviate = 2.4015, p-value = 0.008163
# Moran I statistic       Expectation          Variance 
# 2.258475e-02            -5.590340e-05       8.887844e-05 

set.seed(1234)
moran.m1.mc <- moran.mc(x = gwr.ALL.m1$SDF$gwr.e, listw = nb2listw.regards, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01
moran.m2.mc <- moran.mc(x = gwr.ALL.m2$SDF$gwr.e, listw = nb2listw.regards, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01


save.image("/Volumes/REGARDS/Tabb/Paper 0/09172019RESULTS.RData")
# 6 mins (iMac)


# all data stored up to this point, as of 9.17.2019


#####################################
##################################### GWLR table output
#####################################
# gwr coefficient results
gwr.ALL.m1
gwr.ALL.m2

# anova results
BFC99.gwr.ALL.m1
BFC99.gwr.ALL.m2

# moran's I results
moran.m1
moran.m2

#########################################
#########################################  race maps
#########################################
#########################################
# Model 1
# gwr race coefficients
regardsmergeALL.sp$race.m1 <- gwr.ALL.m1$SDF$"relevel(factor(regardsmerge$Race), ref = \"W\")B"

# gwr race standard errors
regardsmergeALL.sp$race.se.m1 <- gwr.ALL.m1$SDF$"relevel(factor(regardsmerge$Race), ref = \"W\")B_se"

# t-statistic (coefficient/SE) for race
regardsmergeALL.sp$race.t.m1 <- regardsmergeALL.sp$race.m1/regardsmergeALL.sp$race.se.m1

# significant t-statistics for race
regardsmergeALL.sp$race.p.m1 <- as.factor(ifelse(regardsmergeALL.sp$race.t.m1 < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp$race.t.m1 > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp$race.p.m1) # ALL SIGNIFICANT, except 8

# percentage of local estimates that have negative race effects
# 100% of the local estimates at each site have negative race effects for the fully adjusted models (model 3)
sum(regardsmergeALL.sp$race.m1 < 0) / n.final * 100 # 100%


# Model 2
# gwr race coefficients
regardsmergeALL.sp$race.m2 <- gwr.ALL.m2$SDF$"relevel(factor(regardsmerge$Race), ref = \"W\")B"

# gwr race standard errors
regardsmergeALL.sp$race.se.m2 <- gwr.ALL.m2$SDF$"relevel(factor(regardsmerge$Race), ref = \"W\")B_se"

# t-statistic (coefficient/SE) for race
regardsmergeALL.sp$race.t.m2 <- regardsmergeALL.sp$race.m2/regardsmergeALL.sp$race.se.m2

# significant t-statistics for race
regardsmergeALL.sp$race.p.m2 <- as.factor(ifelse(regardsmergeALL.sp$race.t.m2 < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp$race.t.m2 > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp$race.p.m2) # ALL SIGNIFICANT

# percentage of local estimates that have negative race effects
# 100% of the local estimates at each site have negative race effects for the fully adjusted models (model 3)
sum(regardsmergeALL.sp$race.m2 < 0) / n.final * 100 # 100%


# export shapefiles for Steve Melly to create the smoothed maps
#?writeOGR
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")



#########################################
#########################################  maps of REGARDS GWR data 
#########################################  
#########################################
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
#####################
#Set Things Up
#####################
bCRS<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# may need to use this code to transform things
#us=spTransform(usb, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
#plot(us,axes=TRUE)


#################
#State stuff
#statesb<-readShapePoly("/Users/lpp22/Dropbox/Research/County Health Rankings/CHR 2017 Manuscript Share Folder/Shapefiles/cb_2016_us_state_5m.shp", proj4string=CRS("+proj=longlat +init=epsg:4326 +ellps=GRS80"))
#statesb=statesb[statesb$STATEFP != '72' & statesb$STATEFP != '78' & 
#                  statesb$STATEFP != '02' & statesb$STATEFP != '15' & 
#                  statesb$STATEFP != '60' & 
#                  statesb$STATEFP != '66' & statesb$STATEFP != '69',]

#plot(statesb,axes=TRUE)
#statesb = statesb[order(statesb$STATEFP),]
#statesb$ID=1:length(statesb) #I used these elsewhere
#states=spTransform(statesb , CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

data(phenology) # dataset contains a shapefile for US states

plot(us_states)
#str(us_states) # check the projection; which is +proj=longlat +ellps=WGS84

plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp, cex=0.1, pch = 20, col = "blue", add = T)


# subset the states shapefile for non-belt and belt
# belt states: Alabama, Arkansas, Georgia, Louisiana, Mississippi, North Carolina, South Carolina and Tennessee 

spState <- list("sp.polygons", us_states)

#spplot(regardsmergeALL.sp, "race", cuts=c(-0.9, -0.5, -0.4, 0, 0.9), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState)

#contour(regardsmergeALL.sp$race, lwd = 3, add = T)



#########################################
#########################################  initial plots of REGARDS data 
#########################################  NOT FOR PUBLICATION PURPOSES - DUE TO CONFIDENTIALITY (unless we jitter the lat/long values)
#########################################
#Counties: https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
pdf("/Volumes/REGARDS/Tabb/Paper 0/Maps/gwr.race.map.pdf")
spplot(regardsmergeALL.sp, "race", cuts=quantile(regardsmergeALL.sp$race), cex=0.15, main = "", sp.layout = spState)
#spplot(regardsmergeALL.sp, "race", cuts=c(-0.555, -0.5, -0.45, -0.40), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState.ALL)
dev.off()




#########################################
#########################################  stratified analyses
#########################################  
#########################################

###############
###############
#---BLACKS----#
###############
###############
# model 1
#gwr formulas
formula.m1.ALL.gwr.black <- scale(regardsmerge_black$sum7) ~  scale(regardsmerge_black$Age) + factor(regardsmerge_black$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m1.black <- gwr.sel(formula = formula.m1.ALL.gwr.black, data=regardsmergeALL.sp.black, adapt = T)
Sys.time()
# 8 mins
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m1.black <- gwr(formula = formula.m1.ALL.gwr.black, adapt = bw.ALL.m1.black, data = regardsmergeALL.sp.black, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 1 hour
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 12 mins to save image

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m1.black <- BFC99.gwr.test(gwr.ALL.m1.black)
Sys.time()
# 6 mins
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 12 mins to save image

# model 2
#gwr formulas
formula.m2.ALL.gwr.black <- scale(regardsmerge_black$sum7) ~  relevel(factor(regardsmerge_black$REGION), ref = "0") + scale(regardsmerge_black$Age) + factor(regardsmerge_black$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m2.black <- gwr.sel(formula = formula.m2.ALL.gwr.black, data=regardsmergeALL.sp.black, adapt = T)
Sys.time()
# 7 mins 
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 17 mins

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m2.black <- gwr(formula = formula.m2.ALL.gwr.black, adapt = bw.ALL.m2.black, data = regardsmergeALL.sp.black, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 1 hour and 10 mins
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 19 mins

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m2.black <- BFC99.gwr.test(gwr.ALL.m2.black)
Sys.time()
# 6 mins
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 19 mins


# Moran's I for GWR residuals
knearneigh.regards.black <- knearneigh(regardsmergeALL.sp.black)
knn2nb.regards.black <- knn2nb(knearneigh.regards.black)
nb2listw.regards.black <- nb2listw(knn2nb.regards.black)
#Sys.time()
#gwr.morantest.m1 <- gwr.morantest(x = gwr.ALL.m1, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
#Sys.time()
#gwr.morantest.m2 <- gwr.morantest(x = gwr.ALL.m2, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
# memory issues are faced when running gwr.morantest(); "Error: vector memory exhausted (limit reached?)"


# check the residuals for autocorrelation
set.seed(1234)
moran.m1.black <- moran.test(gwr.ALL.m1.black$SDF$gwr.e, nb2listw.regards.black, alternative = "greater", zero.policy = T) 
# Moran I statistic standard deviate = 1.4328, p-value = 0.07596
#Moran I statistic       Expectation          Variance 
#     0.0236263216     -0.0001721170      0.0002758994 

moran.m2.black <- moran.test(gwr.ALL.m2.black$SDF$gwr.e, nb2listw.regards.black, alternative = "greater", zero.policy = T)
# Moran I statistic standard deviate = 1.3741, p-value = 0.08471
# Moran I statistic       Expectation          Variance 
# 0.0226518147            -0.0001721170       0.0002758994

set.seed(1234)
moran.m1.mc.black <- moran.mc(x = gwr.ALL.m1.black$SDF$gwr.e, listw = nb2listw.regards.black, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01
moran.m2.mc.black <- moran.mc(x = gwr.ALL.m2.black$SDF$gwr.e, listw = nb2listw.regards.black, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01
# no sig. remaining spatial variability

save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 11 mins

###############
###############
#---WHITES----#
###############
###############
# model 1
#gwr formulas
formula.m1.ALL.gwr.white <- scale(regardsmerge_white$sum7) ~  scale(regardsmerge_white$Age) + factor(regardsmerge_white$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m1.white <- gwr.sel(formula = formula.m1.ALL.gwr.white, data=regardsmergeALL.sp.white, adapt = T)
Sys.time()
# 24 mins
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 11 mins

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m1.white <- gwr(formula = formula.m1.ALL.gwr.white, adapt = bw.ALL.m1.white, data = regardsmergeALL.sp.white, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 8 hours (2020-01-25 19:17:16 EST to 2020-01-26 04:26:32 EST)
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 12 mins to save image

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m1.white <- BFC99.gwr.test(gwr.ALL.m1.white)
Sys.time()
# 1 hour
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 12 mins to save image

# model 2
#gwr formulas
formula.m2.ALL.gwr.white <- scale(regardsmerge_white$sum7) ~  relevel(factor(regardsmerge_white$REGION), ref = "0") + scale(regardsmerge_white$Age) + factor(regardsmerge_white$Gender)

# gwr bandwidths
Sys.time()
bw.ALL.m2.white <- gwr.sel(formula = formula.m2.ALL.gwr.white, data=regardsmergeALL.sp.white, adapt = T)
Sys.time()
# 1 hour 
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 14 mins

# gwr models
Sys.time()
require(parallel)
cl <- makeCluster(detectCores())
gwr.ALL.m2.white <- gwr(formula = formula.m2.ALL.gwr.white, adapt = bw.ALL.m2.white, data = regardsmergeALL.sp.white, se.fit=T, hatmatrix=T, cl = cl)
Sys.time()
# 10 hours and 40 mins (2020-01-26 07:04:19 EST to 2020-01-26 17:43:33 EST)
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 26 mins

# global test results (OLS vs. GWR)
Sys.time()
BFC99.gwr.ALL.m2.white <- BFC99.gwr.test(gwr.ALL.m2.white)
Sys.time()
# 1 hour
save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
# 19 mins


# Moran's I for GWR residuals
knearneigh.regards.white <- knearneigh(regardsmergeALL.sp.white)
knn2nb.regards.white <- knn2nb(knearneigh.regards.white)
nb2listw.regards.white <- nb2listw(knn2nb.regards.white)
#Sys.time()
#gwr.morantest.m1 <- gwr.morantest(x = gwr.ALL.m1, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
#Sys.time()
#gwr.morantest.m2 <- gwr.morantest(x = gwr.ALL.m2, lw = nb2listw.regards, zero.policy = F)
#Sys.time()
# XX hours XX mins
# memory issues are faced when running gwr.morantest(); "Error: vector memory exhausted (limit reached?)"


# check the residuals for autocorrelation
set.seed(1234)
moran.m1.white <- moran.test(gwr.ALL.m1.white$SDF$gwr.e, nb2listw.regards.white, alternative = "greater", zero.policy = T) 
# Moran I statistic standard deviate = 2.2795, p-value = 0.01132
#Moran I statistic       Expectation          Variance 
#     2.611687e-02     -8.280202e-05      1.321078e-04


moran.m2.white <- moran.test(gwr.ALL.m2.white$SDF$gwr.e, nb2listw.regards.white, alternative = "greater", zero.policy = T)
# Moran I statistic standard deviate = 2.6056, p-value = 0.004586
#Moran I statistic       Expectation          Variance 
#     2.986526e-02     -8.280202e-05      1.321078e-04 

set.seed(1234)
moran.m1.mc.white <- moran.mc(x = gwr.ALL.m1.white$SDF$gwr.e, listw = nb2listw.regards.white, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01
moran.m2.mc.white <- moran.mc(x = gwr.ALL.m2.white$SDF$gwr.e, listw = nb2listw.regards.white, nsim = 99, zero.policy = T) # statistic = 0.67818, observed rank = 100, p-value = 0.01
# sig. spatial variability

save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")


#########################################
#########################################  black maps
#########################################
#########################################
# Model 1
# gwr intercept coefficients
regardsmergeALL.sp.black$int.m1.black <- gwr.ALL.m1.black$SDF$'(Intercept)'

# gwr intercept standard errors
regardsmergeALL.sp.black$int.se.m1.black <- gwr.ALL.m1.black$SDF$"(Intercept)_se"

# t-statistic (coefficient/SE) for intercept
regardsmergeALL.sp.black$int.t.m1.black <- regardsmergeALL.sp.black$int.m1.black/regardsmergeALL.sp.black$int.se.m1.black

# significant t-statistics for intercept
regardsmergeALL.sp.black$int.p.m1.black <- as.factor(ifelse(regardsmergeALL.sp.black$int.t.m1.black < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp.black$int.t.m1.black > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp.black$int.p.m1.black)

# percentage of local estimates that have negative intercept effects
sum(regardsmergeALL.sp.black$int.m1.black < 0) / n.final.black * 100 # 81.9%


# Model 2
# gwr intercept coefficients
regardsmergeALL.sp.black$int.m2.black <- gwr.ALL.m2.black$SDF$'(Intercept)'

# gwr intercept standard errors
regardsmergeALL.sp.black$int.se.m2.black <- gwr.ALL.m2.black$SDF$"(Intercept)_se"

# t-statistic (coefficient/SE) for intercept
regardsmergeALL.sp.black$int.t.m2.black <- regardsmergeALL.sp.black$int.m2.black/regardsmergeALL.sp.black$int.se.m2.black

# significant t-statistics for intercept
regardsmergeALL.sp.black$int.p.m2.black <- as.factor(ifelse(regardsmergeALL.sp.black$int.t.m2.black < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp.black$int.t.m2.black > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp.black$int.p.m2.black)

# percentage of local estimates that have negative intercept effects
sum(regardsmergeALL.sp.black$int.m2.black < 0) / n.final.black * 100 # 68.3%

# export shapefiles for Steve Melly to create the smoothed maps
#?writeOGR
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")



#########################################
#########################################  maps of REGARDS GWR data 
#########################################  
#########################################
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
#####################
#Set Things Up
#####################
bCRS<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# may need to use this code to transform things
#us=spTransform(usb, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
#plot(us,axes=TRUE)


#################
#State stuff
#statesb<-readShapePoly("/Users/lpp22/Dropbox/Research/County Health Rankings/CHR 2017 Manuscript Share Folder/Shapefiles/cb_2016_us_state_5m.shp", proj4string=CRS("+proj=longlat +init=epsg:4326 +ellps=GRS80"))
#statesb=statesb[statesb$STATEFP != '72' & statesb$STATEFP != '78' & 
#                  statesb$STATEFP != '02' & statesb$STATEFP != '15' & 
#                  statesb$STATEFP != '60' & 
#                  statesb$STATEFP != '66' & statesb$STATEFP != '69',]

#plot(statesb,axes=TRUE)
#statesb = statesb[order(statesb$STATEFP),]
#statesb$ID=1:length(statesb) #I used these elsewhere
#states=spTransform(statesb , CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

data(phenology) # dataset contains a shapefile for US states

plot(us_states)
#str(us_states) # check the projection; which is +proj=longlat +ellps=WGS84

plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp, cex=0.1, pch = 20, col = "blue", add = T)
plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp.black, cex=0.1, pch = 20, col = "blue", add = T)
plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp.white, cex=0.1, pch = 20, col = "red", add = T)


# subset the states shapefile for non-belt and belt
# belt states: Alabama, Arkansas, Georgia, Louisiana, Mississippi, North Carolina, South Carolina and Tennessee 

spState <- list("sp.polygons", us_states)

#spplot(regardsmergeALL.sp, "race", cuts=c(-0.9, -0.5, -0.4, 0, 0.9), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState)

#contour(regardsmergeALL.sp$race, lwd = 3, add = T)



#########################################
#########################################  initial plots of REGARDS data 
#########################################  NOT FOR PUBLICATION PURPOSES - DUE TO CONFIDENTIALITY (unless we jitter the lat/long values)
#########################################
#Counties: https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
pdf("/Volumes/REGARDS/Tabb/Paper 0/Maps/gwr.int.black.map.pdf")
spplot(regardsmergeALL.sp.black, "int.m2.black", cuts=quantile(regardsmergeALL.sp.black$int.m2.black), cex=0.15, main = "", sp.layout = spState)
#spplot(regardsmergeALL.sp, "TBD", cuts=c(-0.555, -0.5, -0.45, -0.40), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState.ALL)
dev.off()

#########################################
#########################################  white maps
#########################################
#########################################
# Model 1
# gwr intercept coefficients
regardsmergeALL.sp.white$int.m1.white <- gwr.ALL.m1.white$SDF$'(Intercept)'

# gwr intercept standard errors
regardsmergeALL.sp.white$int.se.m1.white <- gwr.ALL.m1.white$SDF$"(Intercept)_se"

# t-statistic (coefficient/SE) for intercept
regardsmergeALL.sp.white$int.t.m1.white <- regardsmergeALL.sp.white$int.m1.white/regardsmergeALL.sp.white$int.se.m1.white

# significant t-statistics for intercept
regardsmergeALL.sp.white$int.p.m1.white <- as.factor(ifelse(regardsmergeALL.sp.white$int.t.m1.white < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp.white$int.t.m1.white > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp.white$int.p.m1.white)

# percentage of local estimates that have negative intercept effects
sum(regardsmergeALL.sp.white$int.m1.white < 0) / n.final.white * 100 # 52.4%


# Model 2
# gwr intercept coefficients
regardsmergeALL.sp.white$int.m2.white <- gwr.ALL.m2.white$SDF$'(Intercept)'

# gwr intercept standard errors
regardsmergeALL.sp.white$int.se.m2.white <- gwr.ALL.m2.white$SDF$"(Intercept)_se"

# t-statistic (coefficient/SE) for intercept
regardsmergeALL.sp.white$int.t.m2.white <- regardsmergeALL.sp.white$int.m2.white/regardsmergeALL.sp.white$int.se.m2.white

# significant t-statistics for intercept
regardsmergeALL.sp.white$int.p.m2.white <- as.factor(ifelse(regardsmergeALL.sp.white$int.t.m2.white < -1.96, "t < -1.96", ifelse(regardsmergeALL.sp.white$int.t.m2.white > 1.96, "t > 1.96", "NS")))
table(regardsmergeALL.sp.white$int.p.m2.white)

# percentage of local estimates that have negative intercept effects
sum(regardsmergeALL.sp.white$int.m2.white < 0) / n.final.white * 100 # 40.7%

# export shapefiles for Steve Melly to create the smoothed maps
#?writeOGR
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")
#writeOGR(obj = neighmerge3NC.sp, dsn = "/Volumes/Tabb_MESA/Steve/GWR/Update 3/", layer = "neighmerge3NC.sp", driver = "ESRI Shapefile")



#########################################
#########################################  maps of REGARDS GWR data 
#########################################  
#########################################
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
#####################
#Set Things Up
#####################
bCRS<-CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# may need to use this code to transform things
#us=spTransform(usb, CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))
#plot(us,axes=TRUE)


#################
#State stuff
#statesb<-readShapePoly("/Users/lpp22/Dropbox/Research/County Health Rankings/CHR 2017 Manuscript Share Folder/Shapefiles/cb_2016_us_state_5m.shp", proj4string=CRS("+proj=longlat +init=epsg:4326 +ellps=GRS80"))
#statesb=statesb[statesb$STATEFP != '72' & statesb$STATEFP != '78' & 
#                  statesb$STATEFP != '02' & statesb$STATEFP != '15' & 
#                  statesb$STATEFP != '60' & 
#                  statesb$STATEFP != '66' & statesb$STATEFP != '69',]

#plot(statesb,axes=TRUE)
#statesb = statesb[order(statesb$STATEFP),]
#statesb$ID=1:length(statesb) #I used these elsewhere
#states=spTransform(statesb , CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"))

data(phenology) # dataset contains a shapefile for US states

plot(us_states)
#str(us_states) # check the projection; which is +proj=longlat +ellps=WGS84

plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp, cex=0.1, pch = 20, col = "blue", add = T)
plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp.black, cex=0.1, pch = 20, col = "blue", add = T)
plot(us_states, border = "lightgrey")
plot(regardsmergeALL.sp.white, cex=0.1, pch = 20, col = "red", add = T)


# subset the states shapefile for non-belt and belt
# belt states: Alabama, Arkansas, Georgia, Louisiana, Mississippi, North Carolina, South Carolina and Tennessee 

spState <- list("sp.polygons", us_states)

#spplot(regardsmergeALL.sp, "race", cuts=c(-0.9, -0.5, -0.4, 0, 0.9), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState)

#contour(regardsmergeALL.sp$race, lwd = 3, add = T)



#########################################
#########################################  initial plots of REGARDS data 
#########################################  NOT FOR PUBLICATION PURPOSES - DUE TO CONFIDENTIALITY (unless we jitter the lat/long values)
#########################################
#Counties: https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
#States:   https://www.census.gov/geo/maps-data/data/cbf/cbf_state.html
pdf("/Volumes/REGARDS/Tabb/Paper 0/Maps/gwr.int.white.map.pdf")
spplot(regardsmergeALL.sp.white, "int.m2.white", cuts=quantile(regardsmergeALL.sp.white$int.m2.white), cex=0.15, main = "", sp.layout = spState)
#spplot(regardsmergeALL.sp, "TBD", cuts=c(-0.555, -0.5, -0.45, -0.40), cex=0.3, main = "Non-Stroke Belt Region", sp.layout = spState.ALL)
dev.off()




white.est.se <- regardsmergeALL.sp.white[,c("TRACT_KEY", "X", "Y", "int.m1.white", "int.se.m2.white")]
black.est.se <- regardsmergeALL.sp.black[,c("TRACT_KEY", "X", "Y", "int.m1.black", "int.se.m2.black")]

white.est.se@data
black.est.se@data

write.csv(x = white.est.se@data, file = "/Volumes/.sophos_safeguard_REGARDS/Tabb/SteveMelly/GWRResults2/01232020RESULTS.white.csv")
write.csv(x = black.est.se@data, file = "/Volumes/.sophos_safeguard_REGARDS/Tabb/SteveMelly/GWRResults2/01232020RESULTS.black.csv")


save.image("/Volumes/.sophos_safeguard_REGARDS/Tabb/Paper 0/01232020RESULTS.RData")
