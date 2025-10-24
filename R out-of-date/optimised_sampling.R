if(!require(rgdal)){install.packages("rgdal"); library(rgdal)} 
if(!require(raster)){install.packages("raster"); library(raster)} 
if(!require(rgeos)){install.packages("rgeos"); library(rgeos)} 
if(!require(spsann)){install.packages("spsann"); library(spsann)} 
if(!require(sp)){install.packages("sp"); library(sp)} 
if(!require(ICSNP)){install.packages("ICSNP"); library(ICSNP)} 
if(!require(spatstat)){install.packages("spatstat"); library(spatstat)} 
if(!require(dplyr)){install.packages("dplyr"); library(dplyr)} 
if(!require(tibble)){install.packages("tibble"); library(tibble)} 
if(!require(tidyr)){install.packages("tidyr"); library(tidyr)} 
if(!require(geosphere)){install.packages("geosphere"); library(geosphere)} 
if(!require(Rfast)){install.packages("Rfast"); library(Rfast)} 
library(ggplot2)
library(plot3D)
library(rasterVis)
library(RColorBrewer)

### CARICO I LAYER

# area studio
site <- readOGR(dsn = 'gis zone campionamento', layer = 'campionamento meno esclusione')


# ignoranza

igno_map <- raster("MATERIALE PER LA FUNZIONE/Mappa Ignoranza 2 Km.tif")
#igno_map <- raster::crop(igno_map, site)
#igno_map <- raster::mask(igno_map, extent(site))
plot(igno_map)

# ndvi

ndvi_map <- raster("MATERIALE PER LA FUNZIONE//MAPPA NDVI area studio_28m.tif")

newproj <- '+init=EPSG:3035'
ndvi_map <- projectRaster(ndvi_map, crs=newproj) # riproietto al sistema 3035
plot(ndvi_map)


#####################################
### DEFINISCO LA FUNZIONE ###########

sampleboost <- function(ndvi, ignorance, boundary, samp_strategy, nplot, areaplot, perm, ndvi.weight, igno.weight, dist.weight){
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  } # funzione per normalizzare
  
  
  start_time <- Sys.time() # misuro il tempo: start
  result<-list()
  distanze<-matrix(ncol=1, nrow = perm)
  check <- c()
  boundary <- spTransform(boundary, crs(ignorance))
  
  pb <- txtProgressBar(min = 0, max = perm, style = 3)
  for (i in 1:perm){
    punti_random <- spsample(boundary, n=nplot, type= samp_strategy, iter = 5)
    sampling_points <- as(punti_random, "data.frame")
    xy <- sampling_points[,c(1,2)]
    
    spdf <- SpatialPointsDataFrame(coords = xy, data = sampling_points,
                                   proj4string = crs(ignorance))
    
    
    spdf_buffer <- gBuffer(spdf, width=sqrt(areaplot/pi), byid=TRUE )
    
    # Test self intersection
    combos <- combn(nrow(spdf_buffer@data),2)
    int <- c()
    
    for(k in 1:ncol(combos)){
      ii <- combos[1,k]
      j <- combos[2,k]
      
      int[k] <- gIntersects(spdf_buffer[ii, ], spdf_buffer[j,]) # questa riga salta quando i poligoni non si intersecano
      
    }
    
    if (any(int) == TRUE) {check[[i]] <- TRUE} else {check[[i]] <- FALSE}
    
    spectral_values <- raster::extract(ndvi, spdf) # campiono i valori del raster di NDVI
    igno_values <- raster::extract(ignorance, spdf)  # campiono i valori del rater di ignoranza
    
    ## Calcolare distanze con CRS metrico
    # 1. Ottengo a ppp object dai dati
    m <- ppp(xy$x, xy$y, range(xy$x), range(xy$y))
    # 2. Calcolo una matrice di distanze euclidea
    pairwise_distances <- pairdist.ppp(m)
    #distanze[[i]] <- mean(pairwise_distances) # MEDIA DELLE DISTANZE
    pairwise_distances[pairwise_distances == 0] <- NA # sostituisco gli 0 della diagonale con NA
    distanze[[i]] <- min(pairwise_distances, na.rm =TRUE) # distanza minima tra coppie di plot (idea Giuseppe)
    distance_values <- rep(distanze[[i]], nplot)
    
    
    estratti <- data.frame(coordinates(spdf),spectral_values, igno_values,  distance_values)
    names(estratti) <- c("x", "y", "ndvi", "ignorance", "distances")
    
    estratti$INTERSECTION <- check[[i]]
    result[[i]]<-data.frame(estratti)
    
    setTxtProgressBar(pb, i)
    
  }
  
  new_mat<-plyr::ldply(result, data.frame)
  new_mat$try<-as.factor(rep(1:perm, each= nplot))
  
  agg1<-aggregate(new_mat$ndvi,by=list(new_mat$try),FUN=var)
  agg_igno<-aggregate(new_mat$ignorance,by=list(new_mat$try),FUN=mean)
  
  
  agg2<-data.frame(agg1, distanze, agg_igno[[2]], unlist(check))
  colnames(agg2)<-c('Try','Variance','Mean Dist', 'Mean Ignorance', "INTERSECTION")
  agg2 <- na.omit(agg2)
  
  agg2$ndvi_score <- agg2$Variance * ndvi.weight #applico il valore ponderato
  agg2$igno_score <- agg2$`Mean Ignorance` * igno.weight
  agg2$spatial_score <- agg2$`Mean Dist` * dist.weight
  
  agg2$ndvi_norm <- normalize(agg2$ndvi_score) # normalizzo
  agg2$igno_norm <- normalize(agg2$igno_score)
  agg2$spatial_norm <- normalize(agg2$spatial_score)
  
  agg2$FINAL_SCORE <- agg2$ndvi_norm + agg2$igno_norm + agg2$spatial_norm
  
  agg2 <- agg2[agg2$INTERSECTION=="FALSE",] ## elimino le configurazioni spaziali dove c'è intersezione
  
  
  ordered_solutions <- agg2[order(agg2[,'FINAL_SCORE'], decreasing = TRUE),]
  Index <- as.numeric(ordered_solutions[1,1])
  sol <- subset(new_mat[new_mat$try %in% Index,])
  sol2 <- subset(agg2[agg2$Try %in% Index,])
  
  ## Plot best solution
  
  xy_out1 <- sol[,c(1,2)]
  out1_points <- SpatialPointsDataFrame(coords = xy_out1, data = sol, proj4string = crs(boundary))
  out1_buffers <- gBuffer(out1_points, width=sqrt(areaplot/pi), byid=TRUE )
  
  site2 <- spTransform(site, newproj)
  
  plot(site2)
  plot(out1_buffers, add=TRUE)
  
  
  mytheme <- rasterTheme(region = brewer.pal(9, "YlGn"))
  p <- rasterVis::levelplot(ndvi, layers=1, margin = list(FUN = median), par.settings = mytheme)+
    latticeExtra::layer(sp.points(out1_points, lwd= 1.5, col='black'))
  
  p1 <- rasterVis::levelplot(ignorance, layers=1, margin = list(FUN = median), par.settings = RdBuTheme)+
    latticeExtra::layer(sp.points(out1_points, lwd= 0.8, col='darkgray'))
  
  
  if (nrow(new_mat) > 1000) {new_mat <- sample_n(new_mat, 1000)}
  
  
  p2 <-  ggplot(new_mat, aes(x = ndvi, group = try)) +
    geom_density(colour = "lightgrey")+
    theme(legend.position = "none")+
    geom_density(data = sol, aes(x = ndvi, colour = "red"))
  
  
  index_graph <- match(max(agg2$FINAL_SCORE), agg2$FINAL_SCORE)
  scatter3D(agg2$ndvi_norm, agg2$igno_norm, agg2$spatial_norm, bty = "b2", colvar=agg2$FINAL_SCORE, xlab="NDVI", ylab= "IGNORANCE", zlab = "SPACE" ,clab = c("Multiobjective", "Sampling Optimisation"))
  scatter3D(x = agg2$ndvi_norm[index_graph], y = agg2$igno_norm[index_graph], z = agg2$spatial_norm[index_graph], add = TRUE, colkey = FALSE, 
            pch = 18, cex = 3, col = "black") # riga 149 e 150 non funzionano
  
  end_time <- Sys.time()
  message()
  message(paste0(end_time-start_time))
  
  return(list("Full matrix"=new_mat, "Aggregated matrix"=agg2, "Best"= sol, "Variance of sampling points"=sol2[,'Variance'],
              "Mean Ignorance" = sol2[,'igno_score'],
              "Spatial Median of Distance"= sol2[,'Mean Dist'], "Final score"= sol2[,'FINAL_SCORE'], p, p1, p2))
  
  
}


# Uso la funzione


out1 <- sampleboost(ndvi = ndvi_map, ignorance = igno_map, samp_strategy='random', nplot= 12,  areaplot = 10^6, perm = 100, boundary=site,
                    ndvi.weight = 1, igno.weight=1, dist.weight=1)

out1

write.csv(out1$Best, "punti finali campionamento/12plot-100perm-01.csv") # salvo il csv
