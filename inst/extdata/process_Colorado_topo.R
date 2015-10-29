### Data obtained from http://earthexplorer.usgs.gov/
### Download the files gt30w140n40.tif and gt30w140n90.tif

library(rgdal)
library(dplyr)
library(ggplot2)

### Colorado boundary
xlo <- -110
xhi <- -101
ylo <- 36
yhi <- 42

### Load 
X1 <- readGDAL("~/Downloads/gt30w140n40.tif")
X2 <- readGDAL("~/Downloads/gt30w140n90.tif")

meuse.grid$part.a.sel = meuse.grid$part.a
meuse.grid$part.a.sel[meuse.grid$dist >= 0.1] = NA
image(meuse.grid[,,"part.a.sel"],useRasterImage=F)

df1 <- cbind(coordinates(X1),X1@data) %>% 
       filter(x > xlo & x < xhi & y > ylo & y < yhi)
  
df2 <- cbind(coordinates(X2),X2@data) %>% 
  filter(x > xlo & x < xhi & y > ylo & y < yhi)

df <- rbind(df1,df2)
ggplot(df) + geom_tile(aes(x,y,fill=band1)) + scale_fill_gradientn(colours = terrain.colors(10))
head(coordinates(X))

DEM_Colorado <- df
save(DEM_Colorado,file="data/Colorado_topo.rda")

