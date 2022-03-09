source('../scripts/DH_00_config.R')

cli_h1(paste0('REGION: ',REGION, '\t NORMALISED: ', NORMALISED))

## METHOD <<- 'krige'  # gam, spde
## REGION <<- 'IN' #OH
## NORMALISED <<- FALSE
## WEIGHTED <<- FALSE
## RESPONSES <<- c('Cu', 'Zn', 'Pb')

## ---- fitloadData
load(file = paste0(DATA_PATH, 'processed/data.',REGION,'.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/',REGION,'.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.',REGION,'.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.',REGION,'.RData'))
## ----end

## Functions
## ---- fitModelFunctions

DH_make_file_attr <- function(WEIGHTED=FALSE) {
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
myFUN <- function(x, newdata.full, n_pts, Responses, apply_weights) {
    pts <- x[,'id']
    iter_sites <- cbind(x, covars[pts,]) %>% as.data.frame() %>%
        mutate(Longitude = x, Latitude = y) %>%
        st_as_sf(coords=c('Longitude', 'Latitude'), crs = st_crs(4326))

    METRIC <- vector('list', length = length(Responses))
    names(METRIC) <- Responses
    for ( Response in Responses) {
        newdata.sub <- DH_make_and_fit_model_krige(iter_sites, response=Response)
        METRIC[[Response]] <- DH_model_VE(newdata.full, newdata.sub, Response=Response, apply_weights = apply_weights, denom_type = 1)
        METRIC[[Response]] <- METRIC[[Response]] %>% st_drop_geometry() %>% pull(VE)
    }
    BEST_VE <- METRIC %>% map_dbl(~.x)
    imp <- BEST_VE
    ## print(imp)
    imp.wch <- which(imp==max(imp))
    ## print(imp.wch)
    metric <- BEST_VE[imp.wch]
    metric <- mean(BEST_VE)
    ## print(metric)
    ## ifelse(metric<0, 0, metric)
    -metric
}

## ----end

## ---- fitPrepareData

REGION.sf <- rlang::parse_expr(paste0(REGION,'.sf')) %>% rlang::eval_bare()
data.REGION.sf <- rlang::parse_expr(paste0('data.',REGION,'.sf')) %>% rlang::eval_bare()
data.REGION.df <- data.REGION.sf %>%
    mutate(Longitude = st_coordinates(.)[,1],
           Latitude = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    suppressMessages() %>%
    suppressWarnings()
full.coords.REGION <- rlang::parse_expr(paste0('full.coords.',REGION)) %>% rlang::eval_bare()
full.coords.REGION.sf <- rlang::parse_expr(paste0('full.coords.',REGION, '.sf')) %>% rlang::eval_bare()

## Handle normalisation (if it is active)
if (NORMALISED) {
    data.REGION.sf <- data.REGION.sf %>%
        filter(Normalise_by == NORMALISED_BY)
    data.REGION.df <- data.REGION.df %>%
        filter(Normalise_by == NORMALISED_BY) %>%
        droplevels()
}
## newdata.full <- DH_fit_full(REGION.sf = IN.sf,
##                             data.REGION.sf = data.IN.sf,
##                             full.coords.REGION=full.coords.IN)
newdata.full <- DH_fit_full(REGION.sf = REGION.sf,
                            data.REGION.sf = data.REGION.sf,
                            full.coords.REGION=full.coords.REGION,
                            plotit = FALSE)

save(newdata.full, file=paste0(DATA_PATH,'processed/newdata.full', DH_make_file_attr(),'.RData'))
saveRDS(newdata.full, file=paste0(DATA_PATH,'processed/newdata.full', DH_make_file_attr(),'.rds'))

## ----end

## ---- fitPrepareBoundary

sf_use_s2(FALSE)
boundary <- REGION.sf %>%
    st_union() %>% 
    as('Spatial') %>%
    suppressMessages()

## ----end

## ---- fitPrepareSANN

candi <-  cbind(x=data.REGION.df$Longitude, y=data.REGION.df$Latitude)
## covars <- cbind(Cu = data.REGION.sf$Cu,
##                 Zn = data.REGION.sf$Zn,
##                 Pb = data.REGION.sf$Pb,
##                 Distance = data.REGION.sf$Distance)
covars <- data.REGION.sf %>%
    dplyr::select(RESPONSES, Distance) %>%
    st_drop_geometry() %>%
    as.matrix()

schedule <- scheduleSPSANN(initial.temperature = 0.001,
                           cellsize = 0,
                           ## initial.acceptance=c(0.4,1),
                           initial.acceptance=c(0.5,1),
                           ## initial.acceptance=c(0.9,1),
                           chains=100)
tN <- data.REGION.df %>% filter(Site_Type=='random') %>% nrow()
if (tN <=100) sN <- round(seq(2, tN, by=10), -1)
if (tN >100) sN <- round(seq(2, tN, by=20), -1)
sN[1] <- 2

## ----end

## ---- fitPrepareClusters
## cl <- parallel::makeCluster(3)
## doParallel::registerDoParallel(cl)
## handlers(global = TRUE)
## writeLines(c(""), "log.txt")
## ----end

## Loop through a range of sample sizes All analyses include the
## 'fixed' sites.  Each loop increases the number of 'free' sites
## added to the 'fixed' sites and determines the optimum
## configuration.

## ---- optLoop
    optLoop <- function(points, candi, fun, newdata.full,
                        Responses, apply_weights, schedule,
                        plotit, progress) {
        set.seed(123)
        ll = 0
        tt = 0
        while (ll == 0) {
            cat(paste("Try",tt," from iteration:", n, "\n"))
            res = optimUSER(
                points = points,
                candi = candi,
                fun = myFUN,
                newdata.full = newdata.full,
                Responses = RESPONSES,
                apply_weights = apply_weights,
                schedule = schedule,
                plotit = plotit,
                progress = progress)#,
            tt = tt + 1
            if (res$spsann$chains$used!=1) ll = 1
        }
        res
    }
    ## ----end
## ---- fitSampleSizeLoop
## foreach (n=sN, .packages = PACKAGES) %dopar% {
foreach (n=sN, .packages = PACKAGES) %do% {
    ## sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",n,"\n"))
    print(paste("Starting iteration",n,"\n"))

    ## Define the 'fixed' site locations
    points <- list(fixed = data.REGION.df %>%
                       filter(Site_Type=='designated') %>%
                       dplyr::select(Longitude, Latitude) %>%
                       as.matrix(),
                   free = n)

    ## ---- Unweighted
    if (!DRY_RUN) {
        ## res = optLoop(
        ##         points = points,
        ##         candi = candi,
        ##         fun = myFUN,
        ##         newdata.full = newdata.full,
        ##         Responses = RESPONSES,
        ##         apply_weights = FALSE,
        ##         schedule = schedule,
        ##         plotit = FALSE,
        ##         progress = "txt")#,
        set.seed(123)
        ll = 0
        tt = 0
        while (ll == 0) {
            cat(paste("Try",tt," from iteration:", n, " Unweighted\n"))
            print(paste("Try",tt," from iteration:", n, " Unweighted\n"))
            res = optimUSER(
                points = points,
                candi = candi,
                fun = myFUN,
                newdata.full = newdata.full,
                Responses = RESPONSES,
                apply_weights = FALSE,
                schedule = schedule,
                plotit = FALSE,
                progress = "txt")
            tt = tt + 1
            if (res$spsann$chains$used!=1) ll = 1
        }
       ## boundary = boundry)
        save(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(WEIGHTED=FALSE),'.RData'))
        saveRDS(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(WEIGHTED=FALSE),'.rds'))
    }
    cli_h2(paste0('REGION: ',REGION, '\t NORMALISED: ', NORMALISED, '\t Size: ', n, '\t Weighted: ', FALSE))
    ## ----end
    ## ---- Weighted
    if (!DRY_RUN) {
        set.seed(123)
        ll = 0
        tt = 0
        while (ll == 0) {
            cat(paste("Try",tt," from iteration:", n, " Weighted\n"))
            print(paste("Try",tt," from iteration:", n, " Weighted\n"))
            res = optimUSER(
                points = points,
                candi = candi,
                fun = myFUN,
                newdata.full = newdata.full,
                Responses = RESPONSES,
                apply_weights = TRUE,
                schedule = schedule,
                plotit = FALSE,
                progress = "txt",
                boundary = boundry)
            tt = tt + 1
            if (res$spsann$chains$used!=1) ll = 1
        }
        save(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(WEIGHTED=TRUE),'.RData'))
        saveRDS(res, file=paste0(DATA_PATH, 'processed/res:',n,':', DH_make_file_attr(WEIGHTED=TRUE),'.rds'))
    }
    cli_h2(paste0('REGION: ',REGION, '\t NORMALISED: ', NORMALISED, '\t Size: ', n, '\t Weighted: ', TRUE))
    ## ----end
}
## sink()
## ----end
