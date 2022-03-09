source('../scripts/DH_00_config.R')
## ---- Urbanisation.loadData
load(file = paste0(DATA_PATH, 'primary/EA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/data.EA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/data.MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/OH.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/data.OH.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/IN.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/data.IN.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/urban.sf.RData'))
## ----end

sf_use_s2(FALSE)
## ---- Alter.urbanisation
bb <- st_bbox(c(xmin=130.928,ymax=-12.489,xmax=130.939,ymin=-12.497), crs=st_crs(4326)) %>%
    st_as_sfc()

urban.sf <-
    urban.sf %>% st_difference(bb)
## ----end

## ---- EA.distance
data.EA.sf <-
    data.EA.sf %>% 
    mutate(Distance = data.EA.sf %>%
               st_distance(urban.sf %>%
                           st_union()) %>%
               as.vector())
save(data.EA.sf, file = paste0(DATA_PATH, 'processed/data.EA.sf.RData'))
## ----end
## ---- EA.distance.map
b <- data.EA.sf %>% st_nearest_feature(urban.sf)
a <- data.EA.sf %>% st_nearest_points(urban.sf[b,], pairwise = TRUE)
urban.crop.sf <- urban.sf %>% st_crop(urban.sf[b,])
g <- ggplot() +
    geom_sf(data=EA.sf, fill = scales::hue_pal()(1), alpha=0.3) +
    geom_sf(data=urban.crop.sf) +
    geom_sf(data=data.EA.sf) +
    geom_sf(data=a, color='red') +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/EA.urban.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.crop.sf)))
## ----end

## ---- MWA.distance
data.MWA.sf <-
    data.MWA.sf %>% 
    mutate(Distance = data.MWA.sf %>%
               st_distance(urban.sf %>%
                           st_union()) %>%
               as.vector())
save(data.MWA.sf, file = paste0(DATA_PATH, 'processed/data.MWA.sf.RData'))
## ----end
## ---- MWA.distance.map
b <- data.MWA.sf %>% st_nearest_feature(urban.sf)
a <- data.MWA.sf %>% st_nearest_points(urban.sf[b,], pairwise = TRUE)
urban.crop.sf <- urban.sf %>% st_crop(urban.sf[b,])
g <- ggplot() +
    geom_sf(data=MWA.sf, fill = scales::hue_pal()(1), alpha=0.3) +
    geom_sf(data=urban.crop.sf) +
    geom_sf(data=data.MWA.sf) +
    geom_sf(data=a, color='red') +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/MWA.urban.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.crop.sf)))
## ----end

## ---- OH.distance
data.OH.sf <-
    data.OH.sf %>% 
    mutate(Distance = data.OH.sf %>%
               st_distance(urban.sf %>%
                           st_union()) %>%
               as.vector())
save(data.OH.sf, file = paste0(DATA_PATH, 'processed/data.OH.sf.RData'))
## ----end
## ---- OH.distance.map
b <- data.OH.sf %>% st_nearest_feature(urban.sf)
a <- data.OH.sf %>% st_nearest_points(urban.sf[b,], pairwise = TRUE)
urban.crop.sf <- urban.sf %>% st_crop(urban.sf[b,])
g <- ggplot() +
    geom_sf(data=OH.sf, fill = scales::hue_pal()(1), alpha=0.3) +
    geom_sf(data=urban.crop.sf) +
    geom_sf(data=data.OH.sf) +
    geom_sf(data=a, color='red') +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH.urban.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.crop.sf)))
## ----end


## ---- full.domain.IN
full.coords.IN.sf <- IN.sf %>% st_sample(size=10000, type='regular') %>%
    st_as_sf(crs = st_crs(4326)) #%>% 
## full.coords.IN.sf <- IN.sf %>% st_make_grid(cellsize = 0.003) %>% st_intersection(IN.sf %>% st_union()) %>% st_centroid()
g1 <- full.coords.IN.sf %>%
    ggplot() +
    geom_sf(data=IN.sf) +
    geom_sf(shape='+') +
    coord_sf(expand=FALSE) +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/IN.grid.sf.png'),
       g1,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(IN.sf)))
g.urban.IN <- full.coords.IN.sf %>%
    ggplot() +
    geom_sf(data=IN.sf, fill=NA) +
    geom_sf(shape='+') +
    coord_sf(expand=FALSE) +
    theme_bw()
g.urban.IN 
## ----end
## ---- full.domain.IN.distance
full.coords.IN.sf <- full.coords.IN.sf %>%
    mutate(Distance = full.coords.IN.sf %>%
               st_distance(urban.sf %>% st_union()) %>%
               as.vector())
save(full.coords.IN.sf, file = paste0(DATA_PATH, 'processed/full.coords.IN.sf.RData'))
b <- data.IN.sf %>% st_nearest_feature(urban.sf)
urban.crop.sf <- urban.sf %>% st_crop(urban.sf[b,])
g <- full.coords.IN.sf %>%
    ggplot() +
    geom_sf(data=IN.sf, fill = scales::hue_pal()(1), alpha=0.3) +
    geom_sf(data=urban.crop.sf) +
    ## geom_sf(aes(color=sqrt(Distance))) +
    geom_sf(aes(color=Distance)) +
    scale_color_gradientn('Distance (km)', colors=tim.colors(64), trans = scales::sqrt_trans(),
                          labels = function(x) x/1000) +
    theme_bw() +
    theme(legend.position = c(0.01,0.01),
          legend.justification=c(0,0)) +
    guides(color = guide_colorbar(ticks.colour='white', title.position='left', title.theme=element_text(angle=90)))
g
ggsave(filename = paste0(OUTPUT_PATH, 'figures/IN.urban_distances.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.crop.sf)))
g.urban.distance.IN <- g + coord_sf(expand=FALSE)
## ----end
## ---- full.domain.coordinates
full.coords.IN <- full.coords.IN.sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    st_drop_geometry()
save(full.coords.IN, file = paste0(DATA_PATH, 'processed/full.coords.IN.RData'))
## ----end
## ---- full.domain.distance.raster
IN.distance.raster <- full.coords.IN.sf %>%
    st_rasterize(st_as_stars(st_bbox(full.coords.IN.sf),
                             n = 10000,
                             values = NA_real_))
ggplot() + geom_stars(data=IN.distance.raster)
write_stars(IN.distance.raster, paste0(DATA_PATH, 'processed/IN.distance.raster.tif'))
IN.sf.out <- IN.sf %>% dplyr::select(OBJECTID, Zone_Name)
write_sf(IN.sf.out, paste0(DATA_PATH, 'processed/IN.shp'))
## ----end


## ---- full.domain.OH
full.coords.OH.sf <- OH.sf %>% st_sample(size=10000, type='regular') %>%
    st_as_sf(crs = st_crs(4326)) #%>% 
g1 <- full.coords.OH.sf %>%
    ggplot() +
    geom_sf(data=OH.sf) +
    geom_sf(shape='+') +
    coord_sf(expand=FALSE) +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH.grid.sf.png'),
       g1,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(OH.sf)))
g.urban.OH <- full.coords.OH.sf %>%
    ggplot() +
    geom_sf(data=OH.sf, fill=NA) +
    geom_sf(shape='+') +
    coord_sf(expand=FALSE) +
    theme_bw()
## ----end
## ---- full.domain.OH.distance
full.coords.OH.sf <- full.coords.OH.sf %>%
    mutate(Distance = full.coords.OH.sf %>%
               st_distance(urban.sf %>% st_union()) %>%
               as.vector())
save(full.coords.OH.sf, file = paste0(DATA_PATH, 'processed/full.coords.OH.sf.RData'))
b <- data.OH.sf %>% st_nearest_feature(urban.sf)
urban.crop.sf <- urban.sf %>% st_crop(urban.sf[b,])
g <- full.coords.OH.sf %>%
    ggplot() +
    geom_sf(data=OH.sf, fill = scales::hue_pal()(1), alpha=0.3) +
    geom_sf(data=urban.crop.sf) +
    ## geom_sf(aes(color=sqrt(Distance))) +
    geom_sf(aes(color=Distance)) +
    scale_color_gradientn('Distance (km)', colors=tim.colors(64), trans = scales::sqrt_trans(),
                          labels = function(x) x/1000) +
    theme_bw() +
    theme(legend.position = c(0.01,0.99),
          legend.justification=c(0,1)) +
    guides(color = guide_colorbar(ticks.colour='white', title.position='left', title.theme=element_text(angle=90)))
g
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH.urban_distances.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.crop.sf)))
g.urban.distance.OH <- g + coord_sf(expand=FALSE)
## ----end

## ---- full.domain.coordinates.OH
full.coords.OH <- full.coords.OH.sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    st_drop_geometry()
save(full.coords.OH, file = paste0(DATA_PATH, 'processed/full.coords.OH.RData'))
## ----end
## ---- full.domain.distance.raster
OH.distance.raster <- full.coords.OH.sf %>%
    st_rasterize(st_as_stars(st_bbox(full.coords.OH.sf),
                             n = 10000,
                             values = NA_real_))
ggplot() + geom_stars(data=OH.distance.raster)
write_stars(OH.distance.raster, paste0(DATA_PATH, 'processed/OH.distance.raster.tif'))
OH.sf.out <- OH.sf %>% dplyr::select(OBJECTID, Zone_Name)
write_sf(OH.sf.out, paste0(DATA_PATH, 'processed/OH.shp'))
## ----end

## ---- Grid
g <- wrap_plots(g.urban.OH, g.urban.IN, ncol=1) + 
            plot_annotation(tag_levels = 'a', tag_suffix=')') 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH_grid.png'),
       g,
       width = 5.25,
       height = 8)
## ----end
## ---- UrbanDistances
g <- wrap_plots(g.urban.distance.OH, g.urban.distance.IN, ncol=1) + 
            plot_annotation(tag_levels = 'a', tag_suffix=')') 
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH_urbanDistances.png'),
       g,
       width = 5.25,
       height = 8)
## ----end
