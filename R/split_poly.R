# https://gis.stackexchange.com/questions/375345/dividing-polygon-into-parts-which-have-equal-area-using-r
# Dividing Polygon into parts which have equal area

library(sf)
library(mapview)
library(tidyverse)
library(dismo)
library(osmdata)
library(mapview)

split_poly <- function(sf_poly, n_areas){
  # create random points
  points_rnd <- st_sample(sf_poly, size = 10000)
  #k-means clustering
  points <- do.call(rbind, st_geometry(points_rnd)) %>%
    as_tibble() %>% setNames(c("lon","lat"))
  k_means <- kmeans(points, centers = n_areas)
  # create voronoi polygons
  # voronoi_polys <- dismo::voronoi(k_means$centers, ext = sf_poly)
  voronoi_polys <- dismo::voronoi(k_means$centers)
  # clip to sf_poly
  crs(voronoi_polys) <- crs(sf_poly)
  voronoi_sf <- st_as_sf(voronoi_polys)
  equal_areas <- st_intersection(voronoi_sf, sf_poly)
  equal_areas$area <- st_area(equal_areas)
  return(equal_areas)
}



