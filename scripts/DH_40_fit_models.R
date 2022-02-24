source('../scripts/DH_00_config.R')

load(file = paste0(DATA_PATH, 'processed/data.EA.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/EA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.IN.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/IN.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.IN.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.IN.RData'))

## Functions
## ---- fitModelFunctions
DH_make_file_attr <- function() {
    FILE_ATTR <<-  paste0('_',METHOD, '__', REGION, '___',
                         ifelse(NORMALISED, 'norm', 'raw'), '____',
                         ifelse(WEIGHTED, 'wtd','unwtd'))
    FILE_ATTR
}

DH_obs_coords <- function(data.sf) {
    data.sf %>%
        st_cast('POINT') %>%
        mutate(Longitude = st_coordinates(.)[,1],
               Latitude = st_coordinates(.)[,2]) %>%
        st_drop_geometry() %>%
        dplyr::select(Longitude, Latitude, Distance) %>%
        as.matrix()
}

DH_make_mesh <- function(obs.coords, bndry) {
    inla.mesh.2d(loc = obs.coords,
                 boundary = bndry,
                 max.edge = c(500, 500),
                 cutoff = 0.01
                 )
}

## ----end


METHOD <<- 'krige'  # gam, spde
REGION <<- 'IN' #OH
NORMALISED <<- FALSE
WEIGHTED <<- FALSE
RESPONSES <<- c('Cu', 'Zn', 'Pb')

REGION.sf <- rlang::parse_expr(paste0(REGION,'.sf')) %>% rlang::eval_bare()
data.REGION.sf <- rlang::parse_expr(paste0('data.',REGION,'.sf')) %>% rlang::eval_bare()
data.REGION.df <- data.REGION.sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    suppressMessages() %>%
    suppressWarnings()
full.coords.REGION <- rlang::parse_expr(paste0('full.coords.',REGION)) %>% rlang::eval_bare()

sf_use_s2(FALSE)
boundary <- REGION.sf %>%
    st_union() %>% 
    as('Spatial') %>%
    suppressMessages()


candi <-  cbind(x=data.REGION.df$Longitude, y=data.REGION.df$Latitude)
covars <- cbind(Cu = data.REGION.sf$Cu,
                Zn = data.REGION.sf$Zn,
                Pb = data.REGION.sf$Pb,
                Distance = data.REGION.sf$Distance)

schedule <- scheduleSPSANN(initial.temperature = 0.001,
                           cellsize = 0,
                           initial.acceptance=c(0.4,1),
                           chains=100)
tN <- data.REGION.df %>% filter(Site_Type=='random') %>% nrow()
sN <- round(seq(2, tN, by=20), -1)
sN[1] <- 2
for (n in sN[3:5]) {
    points <- list(fixed = data.REGION.df %>%
                       filter(Site_Type=='designated') %>%
                       dplyr::select(Longitude, Latitude) %>%
                       as.matrix(),
                   free = n)
    
    WEIGHTED <<- FALSE
    set.seed(123)
    res = optimUSER(
        points = points,
        candi = candi,
        fun = myFUN,
        newdata.full = newdata.full,
        Responses = RESPONSES,
        apply_weights = WEIGHTED,
        schedule = schedule,
        plotit = FALSE,
        progress = "txt",
        boundary = boundry)
    save(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(),'.RData'))

    WEIGHTED <<- TRUE
    set.seed(123)
    res = optimUSER(
        points = points,
        candi = candi,
        fun = myFUN,
        newdata.full = newdata.full,
        Responses = RESPONSES,
        apply_weights = WEIGHTED,
        schedule = schedule,
        plotit = FALSE,
        progress = "txt",
        boundary = boundry)
    save(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(),'.RData'))
}

