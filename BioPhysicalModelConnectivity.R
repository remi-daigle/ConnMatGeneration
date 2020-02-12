# Bio-physical model connectivity

require(sf)
require(igraph)
require(data.table)
require(tidyverse)

pu <- st_read("pulayer_BC_marine_hx_conservationfeats.shp")

# this data is from: Daigle, R.M. (2015): Cuke Particle Positions - 1998-2007 - Grid sites. Figshare. doi: 10.6084/m9.figshare.1405701
hdd <- "D:/MPA_particles/output"
years <- 1998:2012
days <- 152:182
pld <- 14

filelist <- list.files(file.path(hdd,paste0("G",years)),
                       recursive=TRUE,
                       full.names = TRUE)
#get only correct PLD's
filelist <- filelist[grepl(paste0(pld,".csv",collapse = "|"),filelist)]
#get only correct days
filelist <- filelist[sapply(strsplit(filelist,"/"),'[',5) %in% as.character(days)]

filedf <- filelist %>% 
  strsplit("/") %>% 
  unlist() %>% 
  matrix(ncol=7,byrow=TRUE) %>% 
  as.data.frame() %>% 
  mutate(#V6 = as.numeric(V6),
         file = file.path(V1,V2,V3,V4,V5,V6,V7)) %>% 
  arrange(V4,V5,V7,V6)

filelist <- filedf$file

endpoint <- rbindlist(lapply(filelist, function(x){
  print(x)
  df <- fread(x,select=c(1,2,4))
  setnames(df,c("long","lat","keep"))
  # split <- strsplit(x,"/")[[1]]
  # df[, ':='(
  #   year = as.numeric(gsub("G","",split[4])),
  #   day = as.numeric(split[5]),
  #   subset = as.numeric(split[6]),
  #   pld = as.numeric(gsub("\\.csv","",gsub("para    1","",split[7])))
  # )]
})
)

filelist <- list.files("BC_release_locations/",recursive = TRUE,full.names = TRUE) 
rl <- rbindlist(lapply(filelist, function(x){
  df <- fread(x,select = 1:2)
  setnames(df,c("long","lat"))
  df[, ':='(
    subset = as.numeric(gsub("\\.txt","",gsub("BC_release_locations/rl_","",x)))
  )]
})
) %>% 
  mutate(releaseID=row_number()) 

rl_sf <- rl %>% 
  st_as_sf(coords = c("long","lat"),
           crs = "+proj=longlat +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% 
  st_transform(st_crs(pu))


rl_pu <- as.data.frame(st_intersects(rl_sf,pu)) %>% 
  as.data.table() %>% 
  right_join(rl,by=c("row.id"="releaseID")) %>% 
  select(-long,-lat,-subset,-row.id)

end_rl <- cbind(endpoint,rl_pu) %>% 
  filter(!is.na(col.id)) %>% 
  mutate(id1 = col.id,
         row.id = as.numeric(row.names(.)),
         lat=ifelse(keep==-3,NA,lat)) %>% 
  select(-keep,-col.id) %>% 
  st_as_sf(coords = c("long","lat"),
           crs = "+proj=longlat +no_defs +ellps=WGS84 +towgs84=0,0,0",
           na.fail = FALSE) %>% 
  st_transform(st_crs(pu))

el <- as.data.frame(st_intersects(end_rl,pu)) %>% 
  as.data.table() %>% 
  right_join(select(as.data.table(end_rl),-geometry), by=c("row.id"="row.id")) %>% 
  mutate(id2 = col.id) %>% 
  group_by(id1) %>% 
  mutate(total=n()) %>% 
  ungroup() %>% 
  group_by(id1,id2) %>% 
  summarise(value=n()/min(total)) %>% 
  filter(!is.na(id2))

write.csv(el,"BioPhysical_Flow_EdgeList.csv",row.names = FALSE)  




