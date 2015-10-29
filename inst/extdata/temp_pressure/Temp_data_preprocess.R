library(ggplot2)
library(bicon)
data(temps)
X <- read.csv("~/Downloads/555499.csv",header=T) %>%
  filter( Time_Observation2 < 9999 &
            Time_Observation < 9999) %>%
  mutate(STATION = as.character(STATION)) 
X$DATE2 <- ISOdate(substr(X$DATE,1,4),
                   substr(X$DATE,5,6),
                   substr(X$DATE,7,8),
                   tz="GMT")
Sep19_2004 <- filter(X,DATE == 20040919)
# & 
#                        Measurement_Flag == 9999 & 
#                        Measurement_Flag2 == 9999 & 
#                        Quality_Flag == 9999 & 
#                        Quality_Flag2 == 9999 & 
#                        Source_Flag == "0" &
#                        Source_Flag2 == "0" & 


for(i in 1:nrow(temps)) {
  disparities <- sqrt((Sep19_2004$LONGITUDE - temps$lon[i])^2 + (Sep19_2004$LATITUDE - temps$lat[i])^2)
  min_disparity <- min(disparities)
  idx <- which.min(disparities)
  temps$STATION[i] <- Sep19_2004$STATION[idx]
  temps$min_disparity[i] <- min_disparity
}

Y2 <- left_join(temps,X) %>%
  group_by(DATE) %>%
  mutate(TMAX_anom = TMAX/10 - mean(TMAX/10),
         TMIN_anom = TMIN/10 - mean(TMIN/10)) %>% 
  data.frame()

Y <- filter(Y2,DATE == 20040919)
lm_max <- lm(TMAX_anom ~ ELEVATION,Y)
coeff_max <- lm_max$coefficients
lm_min <- lm(TMIN_anom ~ ELEVATION,Y)
coeff_min <- lm_min$coefficients
Y <- mutate(Y, 
            TMAX_anom2 = TMAX_anom - (coeff_max[1] + coeff_max[2]*ELEVATION),
            TMIN_anom2 = TMIN_anom - (coeff_min[1] + coeff_min[2]*ELEVATION))

lm_max <- lm(TMAX ~ ELEVATION,Y)
coeff_max <- lm_max$coefficients
lm_min <- lm(TMIN ~ ELEVATION,Y)
coeff_min <- lm_min$coefficients
Y <- mutate(Y, 
            TMAX_anom3 = (TMAX - (coeff_max[1] + coeff_max[2]*ELEVATION))/10,
            TMIN_anom3 = (TMIN - (coeff_min[1] + coeff_min[2]*ELEVATION))/10)


coords <- dplyr::select(Y,lon,lat) %>% as.matrix()

## Check that all is OK (not quite)
plot(temps$lon,temps$lat)
lines(Y$LONGITUDE,Y$LATITUDE,type='p',col='red')
par(mfrow=c(2,2))
ymin_int = 10.85
ymax_int = 27.2
plot((Y$lon-Y$LONGITUDE) / diff(range(Y$lon)))
plot((Y$lat-Y$LATITUDE) / diff(range(Y$lat)))
plot((Y$maxT - (Y$TMAX/10 -ymax_int)) / diff(range(Y$maxT)))
plot((Y$minT - (Y$TMIN/10 -ymin_int)) / diff(range(Y$minT)))

## Now get all dates

     

## Do some variogram analysis
library(sp)
library(gstat)
Ysp2 <- data.frame(Y2)
Ysp <- data.frame(Y)
coordinates(Ysp) <- ~lon + lat
maxT.vgm = variogram(maxT~ELEVATION, Ysp)
maxT.fit = fit.variogram(maxT.vgm, model = vgm(1, "Exp", 1, 1))
plot(maxT.vgm, maxT.fit)

minT.vgm = variogram(minT~1, Ysp2)
minT.fit = fit.variogram(minT.vgm, model = vgm(1, "Exp", 1, 1))
plot(minT.vgm, minT.fit)

Y2_cleaned <- group_by(Y2,DATE) %>%
             summarise(n = length(TMAX)) %>%
            filter(n == 94) %>%
             left_join(Y2)

PP <- xtabs(data=Y2_cleaned,TMAX_anom ~ STATION + DATE)

### EOF analysis (deprecated)
Y2_stobj <- stConstruct(t(as(PP,"matrix")),
                        space = list(values=1:nrow(PP)),
                        time=unique(Y2_cleaned$DATE2),
                        SpatialObj = SpatialPoints(coords),
                        interval=TRUE)
eof.data = eof(Y2_stobj)
ggplot(filter(Y2,DATE==20040919) %>% mutate(EOF1 = eof.data["EOF1"]$EOF1)) +
  geom_point(aes(lon,lat,colour=EOF1),size=4) + 
  scale_colour_gradient2(low="blue",mid="light yellow",high="red")

### OR
PPt <- t(PP)
decomp <- svd(PPt)
EOFs <- decomp$v
EOF1 <- cbind(select(Y,LONGITUDE,LATITUDE),EOF = EOFs[,1])

### BACK TO REAL LIFE
library(corpcor)

PP <- xtabs(data=Y2_cleaned,TMAX_anom ~ STATION + DATE)
PP2 <- xtabs(data=Y2_cleaned,TMIN_anom ~ STATION + DATE)
distances <- fields::rdist(coords)
stopifnot(all(rownames(PP) == Y$STATION)) # CHECK ordering hasn't changed

#PPX <- cov(t(PP))
library(corpcor)
PPX2 <- cov.shrink(t(PP))
PPX2 <- cov.shrink(t(PP),)
#PPcor <- cov2cor(PPX2)
#PPcor <- cor.shrink(t(PP))
PPX <- spcov::spcov(as(PPX2,"matrix"),as(PPX2,"matrix"),0.01,step.size=0.1)$Sigma

#PPX <- tcrossprod(PP)/94
PPcor <- cov2cor(PPX)
disp <- 2 - 2 * PPcor
plot(distances,disp)

library(EnviroStat)
tx_coords <- coords/100
idx <- 1:94
dispx <- disp
#dispx <- (outer(Y$TMAX_anom,Y$TMAX_anom,FUN = "-"))
#dispx <- dispx / max(dispx) * 2
#diag(dispx) <- runif(n = 94,min = 0,max=0.2)
png(file="~/Desktop/mygraphic.png",width=2000,height=1750)
sg.est2 <- Falternate3(dispx[idx,idx], tx_coords[idx,], max.iter = 100,
                        alter.lim = 100, model = 1,t0=200,a = c(0.85,1.2))
dev.off()
#coords.grid <- Fmgrid(range(tx_coords[,1]),
#                      range(tx_coords[,2]))
coords.grid <- Fmgrid(c(-1.16,-0.9),
                      c(0.33,0.52))
deform <- Ftransdraw(disp = disp[idx,idx], Gcrds = tx_coords[idx,],
                     MDScrds = sg.est2$ncoords,
                     gridstr = coords.grid,lambda = 0)
Tspline <- sinterp(tx_coords[idx,], sg.est2$ncoords, lam = 0.0001 )
par(mfrow = c(1, 1))
Tgrid <- bgrid(start = apply(tx_coords[idx,],2,mean), xmat = tx_coords[idx,],
               coef = Tspline$sol)
tempplot <- setplot(tx_coords[idx,])
text(tx_coords[idx,], labels = 1:length(idx))
draw(Tgrid, fs = TRUE,pts = F,lcolor=c(1,1))
