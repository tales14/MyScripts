########### CONVEX HULL -------------
#### TALES MARTINS - 18/06/19 --

## PACOTES
library(rgbif)
library(sp)
library(dismo)
library(raster)
library(maptools)
library(GeoRange)

### BAIXANDO OS PONTOS
## dismo
#t <- gbif(genus = "salvia", species = "pratensis subsp. pratensis", 
#              geo = TRUE, removeZeros = TRUE, end = 10000)
#t1 <- t[which(t$country == "Germany"),]
#t2 <- t[which(t$country == "France"),]
#t3 <- t[which(t$country == "Italy"),]
#pontos <- rbind(t1[,c("lon", "lat")], t2[,c("lon", "lat")], t3[,c("lon", "lat")])

## rgbif
r <-occ_search(scientificName = "yucca filamentosa", limit = 1000, return = 'data',
               fields = 'minimal', continent = "north_america", hasCoordinate = TRUE)

pontos <- r[,c("decimalLongitude", "decimalLatitude")]
names(pontos) <- c("lon", "lat")

## FAZENDO O MAPA INICIAL
data("wrld_simpl")
plot(wrld_simpl, xlim=c(-80,80), ylim=c(-80,80), axes = TRUE, col = "lightyellow",
     main = "")

points(pontos$lon, pontos$lat, pch = 16, col = "red")

## LIMPANDO DUPLICATAS E NAs
dups <- duplicated(pontos)
dups_row <- which(dups==TRUE)
length(dups_row)
pontos_nodups <- pontos[-dups_row,]
nas <- which(is.na(pontos_nodups == TRUE))
length(nas)
pontos_clean <- pontos_nodups[-nas,]

pontos_clean <- pontos_nodups

## IDENTIFICANDO PONTOS QUE APARECERAM FORA DO CONTINENTE
identify(pontos_clean$lon, pontos_clean$lat)
pontos_clean <- pontos_clean[-c(278, 391),]

## FAZENDO O MAPA SEM DUPLICATAS E SEM NAs
plot(wrld_simpl, xlim=c(-80,80), ylim=c(-80,80), axes = TRUE, col = "lightyellow",
     main = "")

points(pontos_clean$lon, pontos_clean$lat, 
       pch = 16, col = "green")

## FAZENDO O POLIGONO
hull <- convHull(pontos_clean)
plot(hull, col = "red", add = T)

## VERIFICANDO A AREA DO POLIGONO EM GRAUS QUADRADOS
h <- hull@polygons
h@polygons[[1]]@area
plot(pontos_clean, pch = 20)
lines(h@polygons[[1]]@Polygons[[1]]@coords, col = "red")

## VERIFICANDO A AREA EM KM QUADRADOS
CHullAreaEarth(longs = pontos_clean$lon, lats = pontos_clean$lat)


##################################### FIM ################################

