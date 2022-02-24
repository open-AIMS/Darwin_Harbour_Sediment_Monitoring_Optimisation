source('../scripts/DH_00_config.R')

## ---- LoadData
## ---- Shapefiles
## ---- East Arm shapefiles
EA.sf <- read_sf(paste0(DATA_PATH, '/primary/GIS/East_Arm_Sediment_Sampling.shp')) %>% 
    st_transform(crs=st_crs(4326))

g <- EA.sf %>%
    ggplot() +
    geom_sf() +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/EA.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(EA.sf)))
save(EA.sf, file = paste0(DATA_PATH, 'primary/EA.sf.RData'))
## ----end
## ---- Outer harbour shapefiles
## OH.sf <- read_sf(paste0(DATA_PATH, 'primary/GIS/OuterHarbour_EastArm_ShoalBay.shp')) %>%
OH.sf <- read_sf(paste0(DATA_PATH, 'primary/GIS/OuterHarbour_EastArm_ShoalBay1.shp')) %>%
    filter(!is.na(Zone_Name)) %>%
    st_make_valid() %>%
    st_transform(crs=st_crs(4326)) %>%
    st_difference(EA.sf)

g <- OH.sf %>%
    ggplot() +
    geom_sf() +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(OH.sf)))
save(OH.sf, file = paste0(DATA_PATH, 'primary/OH.sf.RData'))
## ----end
## ---- Middle harbour and West Arm shapefiles
MWA.sf <- read_sf(paste0(DATA_PATH, 'primary/GIS/MA_WA_StudyArea.shp')) %>% 
    st_transform(crs=st_crs(4326)) %>% 
    st_difference(EA.sf %>%
                  filter(OBJECTID == 0) %>%
                  mutate(geometry = geometry + c(0.003,0.003)) %>%
                  st_set_crs(4326) %>%
                  st_buffer(0.002))

g <- MWA.sf %>%
    ggplot() +
    geom_sf() +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/MWA.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(MWA.sf)))
save(MWA.sf, file = paste0(DATA_PATH, 'primary/MWA.sf.RData'))
## ----end
## ---- Combine shapefiles
DH.sf <- OH.sf %>% 
    st_union(MWA.sf) %>%
    st_union(EA.sf) 
g <- DH.sf %>%
    ggplot() +
    geom_sf() +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(DH.sf)))

## ----end
## ---- Urbanisation shapefiles
urban.sf <- read_sf(paste0(DATA_PATH, '/primary/GIS/2020_Built_Focal.shp')) %>% 
    st_transform(crs=st_crs(4326))
g <- urban.sf %>%
    ggplot() +
    geom_sf() +
    theme_bw()
ggsave(filename = paste0(OUTPUT_PATH, 'figures/urban.sf.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(urban.sf)))

save(urban.sf, file = paste0(DATA_PATH, 'primary/urban.sf.RData'))
## ----end
## ---- Darwin Harbour background map
DH.background <- get_stamenmap(c(left=130.61, bottom=-12.75, right=131.07, top = -12.17), zoom=11,
                   maptype = 'terrain-background'
                   )
g <- ggmap(DH.background)

ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH.background.png'),
       g,
       width = fig.width,
       height = (fig.width/bb_aspect_ratio(attr(DH.background,'bb')))[[1]])
## ----end
## ---- Shapefiles map
g <- ggmap(DH.background) +
    geom_sf(data = EA.sf, inherit.aes = FALSE, aes(fill='EA'), alpha=0.3) +
    geom_sf(data = MWA.sf, inherit.aes = FALSE, aes(fill='MWA'), alpha=0.3) +
    geom_sf(data = OH.sf, inherit.aes = FALSE, aes(fill='OH'), alpha=0.3) +
    geom_sf(data = urban.sf, inherit.aes = FALSE, alpha=0.3) +
    theme_bw() +
    scale_fill_manual('Region', breaks = c('EA','MWA','OH'),
                      labels = c('Inner Harbour - East', 'Middle Harbour - West', 'Outer Harbour'),
                      values = scales::hue_pal()(3)[1:3])+
    annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
        pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
        style = north_arrow_fancy_orienteering) +
    theme(axis.title = element_blank(),
          legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
          legend.position = c(0.01,0.99),
          legend.justification = c(0,1))
g
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH_map.png'),
       g,
       width = fig.width*1.5,
       height = (fig.width*1.5/bb_aspect_ratio(attr(DH.background,'bb')))[[1]])

## ----end
## ----end

## ---- Sediment data
## ---- East Arm
## List the sheets
## ---- East Arm.sheets
excel_sheets(paste0(DATA_PATH, 'primary/EA_SP1A_CuPbZn_ML.xlsx'))
## ----end
## Read in the data
## ---- East Arm.load
data.EA <- read_excel(path = paste0(DATA_PATH, 'primary/EA_SP1A_CuPbZn_ML.xlsx'),
                          sheet = 'ALL data',
                          trim_ws = TRUE)
glimpse(data.EA)
## ----end
## Tidy up the data
## ---- East Arm.tidy
data.EA <- data.EA %>%
    dplyr::select(Region = `Pilot Region`,
                  Site = `Site number`,
                  Report_Card_Region = `Report card region`,
                  Sample_Name = `Sample`,
                  Site_ID,
                  Site_Type = `Site type`,
                  Longitude,
                  Latitude,
                  Fe_Al = `Fe/Al`,
                  Cu = `Cu (mg/kg)`,
                  Zn = `Zn (mg/kg)`,
                  Pb = `Pb (mg/kg)`,
                  Cu.norm = `nCu(Al,Fe)`,
                  Zn.norm = `nZn(Al,Fe)`,
                  Pb.norm = `nPb(Al,Fe)`
                  ) %>%
    mutate(Normalise_by = ifelse(Fe_Al>1.3, 'Fe', 'Al'))
save(data.EA, file = paste0(DATA_PATH, 'primary/data.EA.RData'))
## ----end
data.EA
## ---- East Arm.summary
load(file = paste0(DATA_PATH, 'primary/data.EA.RData'))
data.EA
## ----end
## ---- East Arm.sf
data.EA.sf <- data.EA %>%
    st_as_sf(coords = c('Longitude', 'Latitude'), crs = st_crs(4326))
save(data.EA.sf, file = paste0(DATA_PATH, 'primary/data.EA.sf.RData'))
## ----end
## ---- East Arm.map
g <- ggplot() +
    geom_sf(data = EA.sf, inherit.aes = FALSE, alpha=0.3) +
    geom_sf(data = data.EA.sf, inherit.aes = FALSE, aes(fill = Site_Type), shape=21) +
    scale_fill_manual('Site type', breaks=c('designated','random'),
                      values = c('red','black'))+
    theme_bw() +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(axis.title = element_blank(),
          legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
          legend.position = c(0.99,0.99),
          legend.justification = c(1,1))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/EA_map.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(EA.sf)))
## ----end
## ---- East Arm.map2
g <-
    data.EA.sf %>%
    pivot_longer(cols=RESPONSES, names_to='Var', values_to = 'Value') %>%
    group_by(Var) %>%
    nest() %>% 
    mutate(g = map2(.x = data, .y = Var, .f = function(x, y) {
        x %>% ggplot() +
            geom_sf(data = EA.sf, inherit.aes = FALSE, alpha=0.3) +
            geom_sf(inherit.aes = FALSE, aes(fill = Value), shape=21) +
            scale_fill_gradientn(y, colors = tim.colors(64), trans = scales::pseudo_log_trans()) +
            ## scale_fill_manual('Site type', breaks=c('designated','random'),
            ##                   values = c('red','black'))+
            theme_bw() +
            annotation_scale(location = "bl", width_hint = 0.25) +
            annotation_north_arrow(location = "bl", which_north = "true", 
                                   pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                                   style = north_arrow_fancy_orienteering) +
            theme(axis.title = element_blank(),
                  legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
                  legend.position = c(0.99,0.99),
                  legend.justification = c(1,1))
    }
    ))
## wrap_plots(g$g, ncol=1)
ggsave(filename = paste0(OUTPUT_PATH, 'figures/EA_map2.png'),
       wrap_plots(g$g, ncol=1),
       width = fig.width,
       height = length(RESPONSES) * fig.width/bb_aspect_ratio(st_bbox(EA.sf)))
## ----end

## ----end
## ---- Middle Harbour and West Arm
## List the sheets
## ---- Middle Harbour and West Arm.sheets
excel_sheets(paste0(DATA_PATH, 'primary/MWA_SP1C_CuPbZn_ML.xlsx'))
## ----end
## Read in the data
## ---- Middle Harbour and West Arm.load
data.MWA <- read_excel(path = paste0(DATA_PATH, 'primary/MWA_SP1C_CuPbZn_ML.xlsx'),
                          sheet = 'ALL_Data',
                          trim_ws = TRUE)
glimpse(data.MWA)
## ----end
## Tidy up the data
## ---- Middle Harbour and West Arm.tidy
data.MWA <- data.MWA %>%
    dplyr::select(Region = `Pilot Region`,
                  Site = `Site number`,
                  Report_Card_Region = `Report card Region`,
                  Sample_Name = `Sample number`,
                  ## Site_ID,
                  Site_Type = `Site type`,
                  Longitude,
                  Latitude,
                  Fe_Al = `Fe/Al`,
                  Cu = `Cu  (mg/kg)`,
                  Zn = `Zn  (mg/kg)`,
                  Pb = `Pb  (mg/kg)`,
                  Cu.norm = `nCu(Al,Fe)`,
                  Zn.norm = `nZn(Al,Fe)`,
                  Pb.norm = `nPb(Al,Fe)`
                  ) %>%
    mutate(Normalise_by = ifelse(Fe_Al>1.3, 'Fe', 'Al'))
save(data.MWA, file = paste0(DATA_PATH, 'primary/data.MWA.RData'))
## ----end
## ---- Middle Harbour and West Arm.summary
load(file = paste0(DATA_PATH, 'primary/data.MWA.RData'))
data.MWA
## ----end
## ---- Middle Harbour and West Arm.sf
data.MWA.sf <- data.MWA %>%
    st_as_sf(coords = c('Longitude', 'Latitude'), crs = st_crs(4326))
save(data.MWA.sf, file = paste0(DATA_PATH, 'primary/data.MWA.sf.RData'))
## ----end
## ---- Middle Harbour and West Arm.map
g <- ggplot() +
    geom_sf(data = MWA.sf, inherit.aes = FALSE, alpha=0.3) +
    geom_sf(data = data.MWA.sf, inherit.aes = FALSE, aes(fill = Site_Type), shape=21) +
    scale_fill_manual('Site type', breaks=c('designated','random'),
                      values = c('red','black'))+
    theme_bw() +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(axis.title = element_blank(),
          legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
          legend.position = c(0.99,0.99),
          legend.justification = c(1,1))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/MWA_map.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(MWA.sf)))
## ----end
## ---- Middle Harbour and West Arm.map2
g <-
    data.MWA.sf %>%
    pivot_longer(cols=RESPONSES, names_to='Var', values_to = 'Value') %>%
    group_by(Var) %>%
    nest() %>% 
    mutate(g = map2(.x = data, .y = Var, .f = function(x, y) {
        x %>% ggplot() +
            geom_sf(data = MWA.sf, inherit.aes = FALSE, alpha=0.3) +
            geom_sf(inherit.aes = FALSE, aes(fill = Value), shape=21) +
            scale_fill_gradientn(y, colors = tim.colors(64), trans = scales::pseudo_log_trans()) +
            ## scale_fill_manual('Site type', breaks=c('designated','random'),
            ##                   values = c('red','black'))+
            theme_bw() +
            annotation_scale(location = "bl", width_hint = 0.25) +
            annotation_north_arrow(location = "bl", which_north = "true", 
                                   pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                                   style = north_arrow_fancy_orienteering) +
            theme(axis.title = element_blank(),
                  legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
                  legend.position = c(0.99,0.99),
                  legend.justification = c(1,1))
    }
    ))
## wrap_plots(g$g, ncol=1)
ggsave(filename = paste0(OUTPUT_PATH, 'figures/MWA_map2.png'),
       wrap_plots(g$g, ncol=1),
       width = fig.width,
       height = length(RESPONSES) * fig.width/bb_aspect_ratio(st_bbox(EA.sf)))
## ----end

## ----end
## ---- Outer Harbour
## List the sheets
## ---- Outer Harbour.sheets
excel_sheets(paste0(DATA_PATH, 'primary/OHR_SP1B_CuPbZn_ML_Designated.xlsx'))
## ----end
## Read in the data
## ---- Outer Harbour.load
data.OH <- read_excel(path = paste0(DATA_PATH, 'primary/OHR_SP1B_CuPbZn_ML_Designated.xlsx'),
                          sheet = 'All_Data',
                          trim_ws = TRUE)
glimpse(data.OH)

data.OH.designated <- read_excel(path = paste0(DATA_PATH, 'primary/OHR_SP1B_CuPbZn_ML_Designated.xlsx'),
                          sheet = 'designated',
                          trim_ws = TRUE)
## ----end
## Tidy up the data
## ---- Outer Harbour.tidy
data.OH.designated <- data.OH.designated %>%
    dplyr::select(
               Site = `Site number`,
               Longitude,
               Latitude,
               Site_Type = `Site type`
           ) 
data.OH <- data.OH %>%
    dplyr::select(Region = `Pilot region`,
                  Site = `Site number`,
                  Report_Card_Region = `Report card region`,
                  Sample_Name = `Sample Name`,
                  Site_ID,
                  Longitude,
                  Latitude,
                  Fe_Al = `Fe/Al`,
                  Cu = `Cu (mg/kg)`,
                  Zn = `Zn (mg/kg)`,
                  Pb = `Pb (mg/kg)`,
                  Cu.norm = `nCu(Al,Fe)`,
                  Zn.norm = `nZn(Al,Fe)`,
                  Pb.norm = `nPb(Al,Fe)`) %>%
    mutate(Normalise_by = ifelse(Fe_Al>1.3, 'Fe', 'Al')) %>%
    left_join(data.OH.designated) %>%
    mutate(Site_Type=ifelse(is.na(Site_Type), 'random', 'designated'))

save(data.OH, file = paste0(DATA_PATH, 'primary/data.OH.RData'))
## ----end
## ---- Outer Harbour.summary
load(file = paste0(DATA_PATH, 'primary/data.OH.RData'))
data.OH
## ----end
## ---- Outer Harbour.sf
data.OH.sf <- data.OH %>%
    st_as_sf(coords = c('Longitude', 'Latitude'), crs = st_crs(4326))
save(data.OH.sf, file = paste0(DATA_PATH, 'primary/data.OH.sf.RData'))
## ----end
## ---- Outer Harbour.map
g <- ggplot() +
    geom_sf(data = OH.sf, inherit.aes = FALSE, alpha=0.3) +
    geom_sf(data = data.OH.sf, inherit.aes = FALSE, aes(fill = Site_Type), shape=21) +
    scale_fill_manual('Site type', breaks=c('designated','random'),
                      values = c('red','black'))+
    theme_bw() +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(axis.title = element_blank(),
          legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
          legend.position = c(0.01,0.99),
          legend.justification = c(0,1))
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH_map.png'),
       g,
       width = fig.width,
       height = fig.width/bb_aspect_ratio(st_bbox(OH.sf)))
## ----end
## ---- Outer Harbour.map2
g <-
    data.OH.sf %>%
    pivot_longer(cols=RESPONSES, names_to='Var', values_to = 'Value') %>%
    group_by(Var) %>%
    nest() %>% 
    mutate(g = map2(.x = data, .y = Var, .f = function(x, y) {
        x %>% ggplot() +
            geom_sf(data = OH.sf, inherit.aes = FALSE, alpha=0.3) +
            geom_sf(inherit.aes = FALSE, aes(fill = Value), shape=21) +
            scale_fill_gradientn(y, colors = tim.colors(64), trans = scales::pseudo_log_trans()) +
            ## scale_fill_manual('Site type', breaks=c('designated','random'),
            ##                   values = c('red','black'))+
            theme_bw() +
            annotation_scale(location = "bl", width_hint = 0.25) +
            annotation_north_arrow(location = "bl", which_north = "true", 
                                   pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                                   style = north_arrow_fancy_orienteering) +
            theme(axis.title = element_blank(),
                  legend.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
                  legend.position = c(0.99,0.99),
                  legend.justification = c(1,1))
    }
    ))
## wrap_plots(g$g, ncol=1)
ggsave(filename = paste0(OUTPUT_PATH, 'figures/OH_map2.png'),
       wrap_plots(g$g, ncol=1),
       width = fig.width,
       height = length(RESPONSES) * fig.width/bb_aspect_ratio(st_bbox(EA.sf)))
## ----end

## ----end

## ---- Inner harbour (EA + MWA)
data.IN.sf <- data.EA.sf %>%
    rbind(data.MWA.sf %>% mutate(Site_ID=NA))
save(data.IN.sf, file = paste0(DATA_PATH, 'primary/data.IN.sf.RData'))

IN.sf <- EA.sf %>%
    st_union(MWA.sf)
IN.sf %>%
    ggplot() +
    geom_sf()
save(IN.sf, file = paste0(DATA_PATH, 'primary/IN.sf.RData'))
## ----end
## ----fullMap
g <- ggmap(DH.background) +
    ## geom_sf(data = IN.sf, inherit.aes = FALSE, aes(fill='IN'), alpha=0.3, show.legend = 'rect') +
    ## geom_sf(data = OH.sf, inherit.aes = FALSE, aes(fill='OH'), alpha=0.3, show.legend = 'rect') +
    geom_sf(data = IN.sf, inherit.aes = FALSE, aes(fill='IN'), alpha=0.9) +
    geom_sf(data = OH.sf, inherit.aes = FALSE, aes(fill='OH'), alpha=0.9) +
    geom_sf(data = urban.sf, aes(color = 'urban'), inherit.aes = FALSE, alpha=0.3, show.legend = 'line') +
    geom_sf(data = data.IN.sf, inherit.aes = FALSE, color="black", size = rel(2.2)) + 
    geom_sf(data = data.IN.sf, inherit.aes = FALSE, aes(color = Site_Type), show.legend = 'point') + 
    geom_sf(data = data.OH.sf, inherit.aes = FALSE, color="black", size = rel(2.2)) + 
    geom_sf(data = data.OH.sf, inherit.aes = FALSE, aes(color = Site_Type), show.legend = 'point') + 
    theme_bw() +
    scale_fill_manual('Region', breaks = c('IN','OH'),
                      labels = c('Inner Harbour', 'Outer Harbour'),
                      values = c('#dfc27d','#80cdc1'),
                      ## values = scales::hue_pal()(2)[1:2],
                      guide = guide_legend(title = NULL, order = 1,
                                           override.aes = list(linetype = 'blank', shape = NA))) +
    scale_color_manual('', breaks = c('designated', 'random', 'urban'),
                       labels = c('Fixed sites', 'Free sites', 'Water infrastructure'),
                       ## values = c('blue','red','grey40'),
                       values = c('#003f5c','#bc5090','grey40'),
                       guide = guide_legend(title = NULL, order = 2,
                                            override.aes = list(#fill = c('white','white','white'),
                                                                fill = c('#003f5c','#bc5090','white'),
                                                                color=c('black','black','grey40'),
                                                                ## color= c('blue','red','grey40'),
                                                                alpha = 1,
                                                                linetype = c('blank', 'blank', 'solid'),
                                                                shape = c(21,21,NA)))) +
    annotation_scale(location = "bl", width_hint = 0.25) +
    annotation_north_arrow(location = "bl", which_north = "true", 
                           pad_x = unit(0.25, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(axis.title = element_blank(),
          legend.box.background = element_rect(fill='white', color='grey20', size = rel(0.5)),
          legend.background = element_rect(fill = 'white', color='white'),
          legend.box.margin = margin(1,1,1,1, unit='pt'),
          legend.spacing.y = unit(0, 'cm'),
          legend.title = element_text(margin = margin(0,0,0.5,0,'line')),
          legend.position = c(0.01,0.99),
          legend.justification = c(0,1))
g
ggsave(filename = paste0(OUTPUT_PATH, 'figures/DH_mapFull.png'),
       g,
       width = fig.width*1.0,
       height = ((fig.width*1.0)*bb_aspect_ratio(attr(DH.background,'bb')))[[1]],
       dpi=300)

## ----end
## ----end

## ----end
