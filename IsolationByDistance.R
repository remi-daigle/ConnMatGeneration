# Isolation by distance

require(sf)
require(igraph)
require(tidyverse)


pu <- st_read("pulayer_BC_marine_hx_conservationfeats.shp")

spec <- read.table("Marxan_hexagon_pu/input/spec.dat", sep = "\t", header=T)
spec$maxdispersal <- c(rep(25,14),200,50,200,200,50)

el <- st_buffer(pu,5) %>% 
  st_intersects(pu) %>% 
  as.data.frame() 

p2p <- st_distance(st_centroid(pu[39,]),st_centroid(pu[1,])) %>% 
  as.numeric()


master_el <- data.frame()
for(f in spec$id){
  feat_ids <- pu$PUID[as.data.frame(pu)[paste0("FEAT_",f)]>0]
  g <- graph_from_edgelist(as.matrix(el))
  E(g)$weight <- if_else(as.data.frame(pu)[paste0("FEAT_",f)][el$col.id,]>0,p2p,p2p*3)
  
  d <- distances(g)
  d[d>spec$maxdispersal[spec$id==f]*1000] <- Inf
  d[d==0] <- p2p
  
  p <- 1/d
  p <- p/rowSums(p)
  
  g2 <- graph.adjacency(p,weighted=TRUE)
  el2 <- get.data.frame(g2)
  el2 <- el2[el2$to %in% feat_ids & el2$from %in% feat_ids,] %>% 
    mutate(type=paste0("FEAT_",f))
  master_el <- bind_rows(master_el,el2)
}

master_el <- master_el %>% 
  mutate(habitat = type,
         id1 = from,
         id2 = to,
         value = weight) %>% 
  select(habitat,id1,id2, value)
write.csv(master_el,"IsolationByDistance_edgelist.csv",row.names=FALSE)
