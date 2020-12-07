#Libraries
library(spgwr)
library(spatstat)
library(tmap)
library(gstat)
library(sf)
library(raster)
library(rgdal)
library(e1071)
library(spdep)

tmaptools::palette_explorer() #Tool for selecting pallettes

#Set working directory
dir <- "C:/Users/mehnazjahid/Desktop/UVic/GEO518/finalProject/Census"
setwd(dir)

#Reading in particulate matter dataset
#Read in PM2.5 data:
pm2.5 <- readOGR("Pm25Sample.shp") 
pm2.5 <- spTransform(pm2.5, CRS("+init=epsg:26910"))
View(pm2.5@data)
#Reading in dissemination tract and income data
#Read in census income data:
income <- read.csv("Income.csv")  
#Select only ID and Income columns:
colnames(income) <- c("DAUID", "Income") 
#Read in dissemination tract shapefile:
census.tracts <- readOGR("BC_DA.shp") 
census.tracts26910<- spTransform(census.tracts, CRS("+init=epsg:26910"))
#Merge income and dissemination data:
income.tracts <- merge(census.tracts,income, by = "DAUID") 
#Determine the number of columns in the dataframe:
nrow(income.tracts)
#Remove NA values:
income.tracts <- income.tracts[!is.na(income.tracts$Income),]
#Reproject the data:
income.tracts <- spTransform(income.tracts, CRS("+init=epsg:26910"))
View(income.tracts@data)
#Create choropleth map of income:
map_Income <- tm_shape(income.tracts) +
  tm_polygons(col = "Income",
              title = "Median Income",
              style = "jenks",
              palette = "-viridis", n = 6) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"), legend.text.size=.7, legend.title.size= .9)

map_Income

tm_shape(census.tracts26910) + 
  tm_polygons() +
  tm_shape(pm2.5) +
  tm_dots(col="PM25", palette = "-RdBu", 
          title="PM2.5", size=0.3) + 
  tm_legend(legend.outside=TRUE)+
  tm_compass(type="arrow", position=c("right", "TOP"), show.labels = 1, size=0.9, text.size=0.75)+
  tm_grid(x=NA, y= NA, n.x=5, n.y=5,lines = FALSE, labels.size = 0.5)+
  tm_scale_bar(position=c("right", "bottom"))+
  tm_layout(main.title= "Sampled PM2.5 in Metro Vancouver", main.title.size=1,  main.title.position= "center")


#Create a grid called grd to use in your interpolation
# Create an empty grid where n is the total number of cells
grd <- as.data.frame(spsample(pm2.5, "regular", n=5000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
# Create SpatialPixel object:
gridded(grd)     <- TRUE  
# Create SpatialGrid object:
fullgrid(grd)    <- TRUE  
#Reproject the grid:
proj4string(grd) <- proj4string(income.tracts)

### descriptive stats of pm2.5
pm25<- pm2.5$PM25
hist(pm25, xlab="PM2.5", main="Histogram of PM2.5")
min(pm25)
max(pm25)
mean(pm25)
median(pm25)
sd(pm25)
shapiro.test(pm25)$p.value

### descriptive stats of income
Income<- income.tracts$Income
hist(Income, xlab="Income", main="Histogram of income")
min(Income)
max(Income)
mean(Income)
median(Income)
sd(Income)
shapiro.test(Income)$p.value

### point pattern analysis of PM2.5
##QUADRAT ANALYSIS
coordsPM25<- as.data.frame(geom(pm2.5)[,c(2,3)]) # coordinates of PM2.5 dataset
pm25.ext <- as.matrix(extent(pm2.5)) # creating extent
window <- as.owin(list(xrange = pm25.ext[1,], yrange = pm25.ext[2,])) # observation window
pm25.ppp <- ppp(x = coordsPM25$x, y = coordsPM25$y, window = window) # creating ppp object
quads <- 8 # determine the number of qusdrats 
qcount <- quadratcount(pm25.ppp, nx = quads, ny = quads)

plot(pm25.ppp, pch = "+", cex = 0.5)
plot(qcount, add = T, col = "red")

qcount.df <- as.data.frame(qcount)
View(qcount.df)

qcount.df <- plyr::count(qcount.df,'Freq') ##Second, count the number of quadrats with a distinct number of points.
View(qcount.df)
colnames(qcount.df) <- c("x","f") ##Change the column names so that x=number of points and f=frequency of quadrats with x point.

sum.f.x2 <- sum((qcount.df$f)*(qcount.df$x)^2)
M <- quads*quads
N <- sum(qcount.df$f*qcount.df$x)
sum.fx.2 <- (sum(qcount.df$f*qcount.df$x))^2
VAR <- (sum.f.x2 -(sum.fx.2/M))/(M-1)
MEAN <- N/M
VMR <- VAR/ MEAN
chi.square = VMR*(M-1)
p = 1 - pchisq(chi.square, (M - 1))

### obj 1:nspatial segregation of income data
# Global Moran's I
incm<- income.tracts[!is.na(income.tracts@data$Income), ]
incm.nb <- poly2nb(incm) # constructing neighbors list from polygon list using Queen's case
incm.net <- nb2lines(incm.nb, coords=coordinates(income.tracts)) # creating a spatial line data frame
crs(incm.net) <- crs(income.tracts) # setting the CRS as the data

tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(incm.net) + tm_lines(col='red') #mapping the Queen's neighbourhood


incm.nb2 <- poly2nb(income.tracts, queen = FALSE) #Rook's case
incm.net2 <- nb2lines(incm.nb2, coords=coordinates(income.tracts))
crs(incm.net2) <- crs(income.tracts)

tm_shape(income.tracts) + tm_borders(col='lightgrey') + 
  tm_shape(incm.net) + tm_lines(col='blue', lwd = 2) +
  tm_shape(incm.net2) + tm_lines(col='yellow', lwd = 2)

incm.lw2 <- nb2listw(incm.nb2, zero.policy = TRUE, style = "W") # weight matrix of Rook's case
print.listw(incm.lw2, zero.policy = TRUE) # printing weight matrix
incm.lw <- nb2listw(incm.nb, zero.policy = TRUE, style = "W") # weight matrix of Queen's case
print.listw(incm.lw, zero.policy = TRUE) # printing weight matrix

milw <- moran.test(income.tracts$Income, incm.lw, zero.policy = TRUE) # queen's
milw2 <- moran.test(income.tracts$Income, incm.lw2, zero.policy = TRUE) # queen's

mIlw <- milw$estimate[[1]] # extracting the global Moran's I value of queen's case
mIlw2 <- milw2$estimate[[1]] # extracting the global Moran's I value of rook's case
eIlw <- milw$estimate[[2]]# extracting the expected Moran's I value of Queen's case
eIlw2 <- milw2$estimate[[2]]# extracting the expected Moran's I value of Rook's case
varlw <- milw$estimate[[3]] # extracting the variance of Queen's case
varlw2 <- milw2$estimate[[3]] # extracting the variance of Rook's case
zlw <- (mIlw-eIlw)/sqrt(varlw) # calculating the test statistic of Queen's case
zlw2 <- (mIlw2-eIlw2)/sqrt(varlw2) # calculating the test statistic of Rook's case

#local Moran's I
lisa.test <- localmoran(income.tracts$Income, incm.lw, zero.policy = T) # calculating Local Moran's I Queen's
lisa.test2 <- localmoran(income.tracts$Income, incm.lw2, zero.policy = T) # calculating Local Moran's I Rook's

income.tracts$IiQueen <- lisa.test[,1] # storing local Moran's I value for each polygon in the dataset as a variable called Ii for Queens
income.tracts$E.IiQueen<- lisa.test[,2] # storing the expected value of local Moran's I.
income.tracts$Var.IiQueen<- lisa.test[,3] # storing the variance 
income.tracts$Z.IiQueen<- lisa.test[,4] # storing the z-value
income.tracts$P.Queen<- lisa.test[,5] # storing the p-value
income.tracts$IiRook <- lisa.test2[,1] # storing local Moran's I value for each polygon in the dataset as a variable called Ii for Rooks
income.tracts$E.IiRook<- lisa.test2[,2] # storing the expected value of local Moran's I. for Rooks
income.tracts$Var.IiRook<- lisa.test2[,3] # storing the variance Rooks
income.tracts$Z.IiRook<- lisa.test2[,4] # storing the z-value Rooks
income.tracts$P.Rook<- lisa.test2[,5] # storing the p-value
income.tracts$sigQueen<- "not significant"
income.tracts$sigRook<- "not significant"
income.tracts$sigQueen[income.tracts$Z.IiQueen <= -1.96 & income.tracts$IiQueen<income.tracts$E.IiQueen| income.tracts$Z.IiQueen >= 1.96 & income.tracts$IiQueen<income.tracts$E.IiQueen] = "-ve autocorrelation"
income.tracts$sigQueen[income.tracts$Z.IiQueen <= -1.96 & income.tracts$IiQueen>income.tracts$E.IiQueen| income.tracts$Z.IiQueen >= 1.96 & income.tracts$IiQueen>income.tracts$E.IiQueen] = "+ve autocorrelation" 
income.tracts$sigRook[income.tracts$Z.IiRook <= -1.96 & income.tracts$IiRook<income.tracts$E.IiRook| income.tracts$Z.IiRook >= 1.96 & income.tracts$IiRook<income.tracts$E.IiRook] = "+ve autocorrelation"
income.tracts$sigRook[income.tracts$Z.IiRook <= -1.96 & income.tracts$IiRook>income.tracts$E.IiRook| income.tracts$Z.IiRook >= 1.96 & income.tracts$IiRook>income.tracts$E.IiRook] = "+ve autocorrelation" 
View(income.tracts@data)

count(income.tracts$sigQueen=="+ve autocorrelation")/length(income.tracts$sigQueen)
count(income.tracts$sigQueen=="-ve autocorrelation")/length(income.tracts$sigQueen)
count(income.tracts$sigQueen=="not significant")/length(income.tracts$sigQueen)

count(income.tracts$sigRook=="+ve autocorrelation")/length(income.tracts$sigRook)
count(income.tracts$sigRook=="-ve autocorrelation")/length(income.tracts$sigRook)
count(income.tracts$sigRook=="not significant")/length(income.tracts$sigRook)
########################
map_LISA.Queen <- tm_shape(income.tracts) + # mapping the local Moran's I
  tm_polygons(col = "sigQueen", 
              title = "Local Moran's I (Queen's)", 
              style = "cat", 
              palette = "viridis"
              ) +
  tm_legend(legend.position = c("LEFT", "BOTTOM"), legend.text.size=.6, legend.title.size= .9)
map_LISA.Queen

map_LISA.Rook <- tm_shape(income.tracts) + # mapping the local Moran's I
  tm_polygons(col = "sigRook", 
              title = "Local Moran's I (Rook's)", 
              style = "cat", 
              palette = "viridis")+
  tm_legend(legend.position = c("LEFT", "BOTTOM"), legend.text.size=.6, legend.title.size= .9)
map_LISA.Rook

### spatial interpolation by kriging
# Create an empty grid where n is the total number of cells
f.0 <- as.formula(PM25 ~ 1)
#Create variogram
var.pm <- variogram(f.0, pm2.5, cloud = FALSE) 

dat.fit.exp  <- fit.variogram(var.pm, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=0.015, model="Exp", range=8000, nugget=0))
dat.fit.gau  <- fit.variogram(var.pm, fit.ranges = FALSE, fit.sills = FALSE,
                              vgm(psill=0.015, model="Gau", range=8000, nugget=0))
dat.fit.sph  <- fit.variogram(var.pm, fit.ranges = FALSE, fit.sills = FALSE,
                              vgm(psill=0.015, model="Sph", range=8000, nugget=0))

par(mfrow=c(1,3))
plot(var.pm, dat.fit.exp, main= "Exponential")
plot(var.pm, dat.fit.gau, main= "Gaussian")
plot(var.pm, dat.fit.sph, main= "Spherical")
par(mfrow=c(1,1))
dat.krg <- krige(f.0, pm2.5, grd, dat.fit.exp)
# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, income.tracts)
# Plot the map
tm_shape(r.m) + 
  tm_raster(n=10, palette="-RdBu",  
            title="Predicted PM2.5", midpoint = NA) +
  tm_shape(pm2.5) + tm_dots(size=0.06)+
tm_legend(legend.outside=TRUE)

#Extract average pm2.5 for each polygon
income.tracts$Pm2.5 <- round(extract(r, income.tracts, fun = mean)[,1], 5)
View(income.tracts@data)

######Linear Regression##########
#Let's say your dataset with both PM2.5 and Income 
#are stored in a dataset called income.tracts.
#Plot income and PM2.5 from the income.tracts dataset you created
plot(income.tracts$Pm2.5~income.tracts$Income)

#Notice that there are a lot of 0's in this dataset. If you decide to remove them, use the following line:
income.tracts.no0 <-  income.tracts[which(income.tracts$Pm2.5 > 0), ]

#Now plot the data again
plot(income.tracts.no0$Pm2.5~income.tracts.no0$Income)

#Perform a linear regression on the two variables. You should decide which one is dependent.
lm.model <- lm(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
#Add the regression model to the plot you created
plot(income.tracts.no0$Income~income.tracts.no0$Pm2.5)
abline(lm.model, col = "red")
#Get the summary of the results
summary(lm.model)

#add the fitted values to your spatialpolygon dataframe
income.tracts.no0$predictlm <- lm.model$fitted.values

#You want to determine if the model residuals are spatially clustered. 
#add the residuals to your spatialpolygon dataframe
income.tracts.no0$residuals <- residuals.lm(lm.model)

#Observe the result to make sure it looks correct
head(income.tracts.no0)

#Now, create choropleth map of residuals
map_resid <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "residuals",
              title = "Residuals",
              style = "jenks",
              palette = "Oranges", n = 6)+
  tm_legend(legend.outside=TRUE)
map_resid

## Check whether there is any spatial autocorrelation in residuals
# global Moran's I on residuals
resid<- income.tracts.no0[!is.na(income.tracts.no0@data$residuals), ]
resid.nb <- poly2nb(resid) # constructing neighbors list from polygon list using Queen's case
resid.net <- nb2lines(resid.nb, coords=coordinates(income.tracts.no0)) # creating a spatial line data frame
crs(resid.net) <- crs(income.tracts.no0) # setting the CRS as the data

resid.nb2 <- poly2nb(income.tracts.no0, queen = FALSE) #Rook's case
resid.net2 <- nb2lines(resid.nb2, coords=coordinates(income.tracts.no0))
crs(resid.net2) <- crs(income.tracts.no0)

resid.lw <- nb2listw(resid.nb, zero.policy = TRUE, style = "W") # weight matrix of Queen's case
print.listw(resid.lw, zero.policy = TRUE) # printing weight matrix
resid.lw2 <- nb2listw(resid.nb2, zero.policy = TRUE, style = "W") # weight matrix of Rook's case
print.listw(resid.lw2, zero.policy = TRUE) # printing weight matrix

milwresid <- moran.test(income.tracts.no0$residuals, resid.lw, zero.policy = TRUE) # queen's
milw2resid <- moran.test(income.tracts.no0$residuals, resid.lw2, zero.policy = TRUE) # Rook's

zlwresid <- (milwresid$estimate[[1]]-milwresid$estimate[[2]])/sqrt(milwresid$estimate[[3]]) # calculating the test statistic of Queen's case
zlw2resid <- (milw2resid$estimate[[1]]-milw2resid$estimate[[2]])/sqrt(milw2resid$estimate[[3]]) # calculating the test statistic of Rook's case

####Geographically Weighted Regression
#Let's say you are continuing with 
#your data from the regression analysis. 
#The first thing you need to do is to add the 
#polygon coordinates to the spatialpolygondataframe.
#You can obtain the coordinates using the 
#"coordinates" function from the sp library
income.tracts.no0.coords <- sp::coordinates(income.tracts.no0)
#Observe the result:
head(income.tracts.no0.coords)
#Now add the coordinates back to the spatialpolygondataframe
income.tracts.no0$X <- income.tracts.no0.coords[,1]
income.tracts.no0$Y <- income.tracts.no0.coords[,2]

###Determine the bandwidth for GWR: this will take a while
GWRbandwidth <- gwr.sel(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                        data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y),adapt=T) 

###Perform GWR on the two variables with the bandwidth determined above
###This will take a looooooong while
gwr.model = gwr(income.tracts.no0$Income~income.tracts.no0$Pm2.5, 
                data=income.tracts.no0, coords=cbind(income.tracts.no0$X,income.tracts.no0$Y), 
                adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 

#Print the results of the model
gwr.model

#Look at the results in detail
results<-as.data.frame(gwr.model$SDF)
head(results)

#Now for the magic. Let's add our local r-square values to the map
income.tracts.no0$localr <- results$localR2

#Create choropleth map of r-square values
map_r2 <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "localr",
              title = "R2 values",
              style = "jenks",
              palette = "viridis", n = 6)+
  tm_legend(legend.outside=TRUE)
map_r2

#Time for more magic. Let's map the coefficients
income.tracts.no0$coeff <- results$income.tracts.no0.Pm2.5
View(income.tracts.no0@data)
#Create choropleth map of the coefficients
map_coef <- tm_shape(income.tracts.no0) +
  tm_polygons(col = "coeff",
              title = "Coefficients",
              style = "jenks",
              palette = "-Greens", n = 6)+
  tm_legend(legend.outside=TRUE)
map_coef
