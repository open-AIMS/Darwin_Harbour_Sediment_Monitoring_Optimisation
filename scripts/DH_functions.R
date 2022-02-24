GOLDEN_RATIO <- (1 + sqrt(5))/2

bb_aspect_ratio <- function(bb) {
    (bb[3] - bb[1])/(bb[4] - bb[2])
}



DH_fit_full <- function(REGION.sf, data.REGION.sf, full.coords.REGION, plotit = FALSE) {
    if (METHOD == 'INLA') {
        ## Create a boundary
        ## ---- Common setup
        sf_use_s2(FALSE)
        bndry <- REGION.sf %>%
            st_union() %>%
            st_coordinates() %>% 
            `[`(,c('X','Y'))
        save(bndry, file=paste0(DATA_PATH, 'processed/bndry_',REGION,'.RData'))
        load(file=paste0(DATA_PATH, 'processed/bndry_',REGION,'.RData'))
        
        obs.coords <- DH_obs_coords(data.REGION.sf)
        save(obs.coords, file=paste0(DATA_PATH, 'processed/obs.coords_',REGION,'.RData'))
        load(file=paste0(DATA_PATH, 'processed/obs.coords_',REGION,'.RData'))
        ## mesh <- DH_make_mesh(obs.coords, bndry)
        mesh <- inla.mesh.2d(
            loc = obs.coords[,c('Longitude','Latitude')],
            ## boundary = bndry,
            max.edge = c(0.015, 0.05),
            cutoff = 0.001,
            offset = c(0.05, 0.05)
        )
        save(mesh, file=paste0(DATA_PATH, 'processed/mesh_',REGION,'.RData'))
        load(file=paste0(DATA_PATH, 'processed/mesh_',REGION,'.RData'))
        
        g <- ggplot() +
            gg(mesh) +
            geom_sf(data=IN.sf, fill='orange', alpha=0.1) +
            geom_sf(data=data.IN.sf) +
            theme_bw() +
            theme(axis.title = element_blank())
        ggsave(filename=paste0(OUTPUT_PATH,'figures/mesh_', REGION, '.png'), g, width=6, height = 6, dpi = 300)
        ggsave(filename=paste0(OUTPUT_PATH,'figures/mesh_', REGION, '.pdf'), g, width=6, height = 6, dpi = 300)
        ## full.coords <- IN.sf %>% st_sample(size=10000, type='regular') %>%
        ##     st_as_sf() %>%
        ##     mutate(Longitude = st_coordinates(.)[,1],
        ##            Latitude = st_coordinates(.)[,2]) %>%
        ##     st_drop_geometry()
        
        proj.grid <- inla.mesh.projector(mesh, loc=as.matrix(full.coords.REGION[,c('Longitude', 'Latitude')]))
        proj.grid.obs <- inla.mesh.projector(mesh, loc=as.matrix(obs.coords))
        spde <- inla.spde2.matern(mesh, alpha = 2)
        i.spatial <- inla.spde.make.index('spatial.field',
                                          n.spde = spde$n.spde)
        ## ----end
        newdata.full <- vector('list', length = length(RESPONSES))
        names(newdata.full) <- RESPONSES
        newdata.obs.full <- newdata.full
        for (Response in RESPONSES) {
            ## ---- full model
            mod <- DH_make_and_fit_model(mesh, obs.coords, data.REGION.sf, response=Response)
            newdata.full.i <- DH_predict_model(proj.grid, mod, full.coords.REGION, Response)
            newdata.obs.full.i <- DH_predict_model(proj.grid.obs, mod, obs.coords, Response)
            ## ----end
            ## ---- predict
            mod.rf <- inla.mesh.project(proj.grid,
                                        mod$summary.random$spatial.field$mean +
                                        mod$summary.fixed['Intercept', 'mean'])
            mod.obs.rf <- inla.mesh.project(proj.grid.obs,
                                            mod$summary.random$spatial.field$mean +
                                            mod$summary.fixed['Intercept', 'mean'])
            newdata.full[[Response]] <- full.coords.REGION %>%
                as.data.frame() %>%
                mutate({{Response}} := exp(mod.rf)) %>% 
                mutate(Distance = full.coords.REGION$Distance)
            newdata.obs.full[[Response]] <- obs.coords %>%
                as.data.frame() %>%
                mutate({{Response}} := exp(mod.obs.rf)) %>%
                mutate(Distance = obs.coords[,'Distance'])
            ## ----end
        }
        ## ---- collate
        newdata.full <- newdata.full %>%
            reduce(left_join) %>%
            dplyr::select(-matches('.*\\.([xy]|pred|var)$')) %>%
            suppressMessages() %>%
            suppressWarnings()
        ## ----end

        ## ---- plots
        xbreaks <- REGION.sf %>% st_bbox() %>% `[`(c(1,3)) %>% 
            scales::breaks_width(width=0.1)()
        ybreaks <- REGION.sf %>% st_bbox() %>% `[`(c(2,4)) %>% 
            scales::breaks_width(width=0.1)()
        g1 <- vector('list', length(RESPONSES))
        names(g1) <- RESPONSES
        for (resp in RESPONSES) {
            lims <- c(newdata.full %>% pull(resp)) %>% range()
            lims[1] <- floor(lims[1])
            lims[2] <- ceiling(lims[2])
            g1[[resp]] <-
                ggplot() +
                geom_raster(data=newdata.full, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
                geom_sf(data=REGION.sf %>% st_union() %>% suppressMessages() %>% suppressWarnings(), fill=NA) +
                scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
                scale_x_continuous(breaks = xbreaks) +
                scale_y_continuous(breaks = ybreaks) +
                coord_sf(expand=FALSE) +
                theme_bw() +
                theme(axis.title = element_blank()) %>%
                suppressWarnings() %>%
                suppressMessages()
        }
        ## wrap_plots(g1, ncol=2)
        
        g2 <- ggplot() +
            geom_sf(data = data.REGION.sf, aes(color = Site_Type), size=1) +
            geom_sf(data = REGION.sf %>% st_union() %>% suppressMessages(), fill=NA) +
            scale_x_continuous(breaks = xbreaks) +
            scale_y_continuous(breaks = ybreaks) +
            coord_sf(expand=FALSE) +
            theme_bw() +
            scale_color_manual('Site type', breaks = c('designated', 'random'), labels = c('Fixed', 'Free'),
                               values = scales::hue_pal()(2))
        g1[[length(g1)+1]] <- g2
        
        g <- wrap_plots(c(g1)) +
            plot_annotation(tag_levels = 'a', tag_suffix=')') &
            theme(plot.tag.position=c(0,1))
        ggsave(filename=paste0(OUTPUT_PATH, 'figures/',REGION,'.full_', METHOD,'.png'), g,
               width = 10, height = 10, dpi=300)
        ggsave(filename=paste0(OUTPUT_PATH, 'figures/',REGION,'.full_', METHOD,'.pdf'), g,
               width = 10, height = 10)
        ## ----end
    } else if (METHOD =='krige') {
        ## ---- full model
        newdata.full.k <- vector('list', length = length(RESPONSES))
        names(newdata.full.k) <- RESPONSES
        for (Response in RESPONSES) {
            newdata.full.k[[Response]] <- DH_make_and_fit_model_krige(data.REGION.sf, response=Response) %>%
                mutate(Distance = full.coords.REGION$Distance)
        }
        ## ----end
        ## ---- collate
        newdata.full.k <- newdata.full.k %>%
            reduce(st_join) %>%
            dplyr::select(-matches('.*\\.([xy]|pred|var)$')) %>%
            suppressMessages() %>%
            suppressWarnings()
        ## ----end
        ## ---- plots
        if (plotit) {
            xbreaks <- REGION.sf %>% st_bbox() %>% `[`(c(1,3)) %>% 
                scales::breaks_width(width=0.1)()
            ybreaks <- REGION.sf %>% st_bbox() %>% `[`(c(2,4)) %>% 
                scales::breaks_width(width=0.1)()
            g1 <- vector('list', length(RESPONSES))
            names(g1) <- RESPONSES
            for (resp in RESPONSES) {
                lims <- c(newdata.full.k %>% pull(resp)) %>% range()
                lims[1] <- floor(lims[1])
                lims[2] <- ceiling(lims[2])
                g1[[resp]] <-
                    ggplot() +
                    geom_raster(data=newdata.full.k, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
                    geom_sf(data=REGION.sf %>% st_union() %>% suppressMessages() %>% suppressWarnings(), fill=NA) +
                    scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
                    scale_x_continuous(breaks = xbreaks) +
                    scale_y_continuous(breaks = ybreaks) +
                    coord_sf(expand=FALSE) +
                    theme_bw() +
                    theme(axis.title = element_blank()) %>%
                    suppressWarnings() %>%
                    suppressMessages()
            }
            ## wrap_plots(g1, ncol=2)
            
            g2 <- ggplot() +
                geom_sf(data = data.REGION.sf, aes(color = Site_Type), size=1) +
                geom_sf(data = REGION.sf %>% st_union() %>% suppressMessages(), fill=NA) +
                scale_x_continuous(breaks = xbreaks) +
                scale_y_continuous(breaks = ybreaks) +
                coord_sf(expand=FALSE) +
                theme_bw() +
                scale_color_manual('Site type', breaks = c('designated', 'random'), labels = c('Fixed', 'Free'),
                                   values = scales::hue_pal()(2))
            g1[[length(g1)+1]] <- g2
            
            g <- wrap_plots(c(g1)) +
                plot_annotation(tag_levels = 'a', tag_suffix=')') &
                theme(plot.tag.position=c(0,1))
            suppressWarnings(print(g))
            
            ggsave(filename=paste0(OUTPUT_PATH, 'figures/',REGION,'.full_', METHOD,'.png'),
                   g,
                   width = 10, height = 10, dpi=300) %>%
                suppressWarnings()
            ggsave(filename=paste0(OUTPUT_PATH, 'figures/',REGION,'.full_', METHOD,'.pdf'), g,
                   width = 10, height = 10) %>%
                suppressWarnings()
        }
        ## ----end
        newdata.full <- newdata.full.k
    }
    newdata.full
}




DH_make_and_fit_model_gam <- function(data.sf, response) {
    ## form <- y ~ te(Latitude) + te(Longitude) + ti(Longitude, Latitude) 
    form <- y ~ te(Latitude) + te(Longitude) + te(Longitude, Latitude) 
    ## form <- y ~ s(Longitude, Latitude) 
    ## form <- y ~ ti(Longitude, Latitude) 
    data.sf <- data.sf %>% mutate(y := .[[response]],
                                  Longitude = st_coordinates(.)[,1],
                                  Latitude = st_coordinates(.)[,2]) %>%
        st_drop_geometry() %>%
        as.data.frame()
    mod <- gam(form, data=data.sf, family = "Gamma")
    mod
}


## This one
DH_make_and_fit_model_krige <- function(data.sf, response) {
  full.coords.REGION.sf <<- full.coords.REGION.sf %>%
    st_set_crs(st_crs(data.sf)) %>%
    suppressMessages() %>%
    suppressWarnings()
  form <- log(y) ~ 1
  data.sf <- data.sf %>% mutate(y := .[[response]])
  krige(form, locations = data.sf, newdata = full.coords.REGION.sf, debug.level = 0) %>%
        mutate({{response}} := exp(var1.pred),
               Longitude = st_coordinates(.)[,1],
               Latitude = st_coordinates(.)[,2]) %>%
        st_join(full.coords.REGION.sf %>% dplyr::select(Distance)) %>%
        suppressMessages() %>%
        suppressWarnings()
}


DH_predict_model <- function(proj.grid, mod, full.coords.IN, Response) {
    mod.rf <- inla.mesh.project(proj.grid,
                                mod$summary.random$spatial.field$mean +
                                mod$summary.fixed['Intercept', 'mean'])
    full.coords.IN %>%
        as.data.frame() %>%
        mutate({{Response}} := exp(mod.rf))
}
DH_predict_model_gam <- function(mod, full.coords.IN, Response) {
    full.coords.IN %>%
        as.data.frame() %>%
        mutate({{Response}} := as.vector(predict(mod, newdata = full.coords.IN, type='response')))
}
DH_make_and_fit_model <- function(mesh, obs.coords, data.sf, response) {
    spde <- inla.spde2.matern(mesh,
                              alpha = 2)
    i.spatial <- inla.spde.make.index('spatial.field',
                                      n.spde = spde$n.spde)
    A.est <- inla.spde.make.A(mesh = mesh,
                              loc = obs.coords,
                              )
    e.list <- list(c(i.spatial, list(Intercept = 1)))

    stack.est <- inla.stack(data = list(y = data.sf %>% pull(response)),
                            A = list(A.est),
                            effects = e.list,
                            remove.unused = FALSE,
                            tag = 'est')
    form <- y ~  0 + Intercept + f(spatial.field, model = spde)
    mod <- inla(form,
                data = inla.stack.data(stack.est),
                family= 'gamma',
                control.predictor = list(compute = TRUE,
                                         link = 1,
                                         A = inla.stack.A(stack.est)
                                         ),
                control.compute = list(config = TRUE),
                verbose = FALSE,
                silent = 2L)
    mod
}


DH_model_predictions <- function(proj.grid, mod, coords, Response) {
    mod.rf <- inla.mesh.project(proj.grid,
                                    mod$summary.random$spatial.field$mean +
                                    mod$summary.fixed['Intercept', 'mean'])
    coords %>%
        as.data.frame() %>%
        mutate({{Response}} := exp(mod.rf))
    
}

DH_model_MSE <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                        apply_weights = apply_weights, 
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(Diff = .[[Response]] - A) %>%
        summarise(MSE = (sum((Diff^2)*Weight)))   
}
DH_model_RMSE <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(Diff = .[[Response]] - A) %>%
        summarise(RMSE = sqrt(sum((Diff^2)*Weight)))   
}
DH_model_PERC <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                            apply_weights = apply_weights, 
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(PERC = abs(.[[Response]] - A)/.[[Response]],
               Weight = Weight/sum(Weight)) %>%
        summarise(PERC = sum(PERC*Weight))   
}

DH_model_SS <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                            apply_weights = apply_weights, 
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(Weight = Weight/sum(Weight)) %>%
        summarise(
            SS_FULL = sum((Weight*(mean(.[[Response]])-.[[Response]]))^2),
            SS_SUB = sum((Weight*(mean(A)-A))^2)) %>%
        mutate(SS = SS_SUB/SS_FULL)
}
DH_model_VE <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE, denom_type=1) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                            apply_weights = apply_weights,
                            denom_type = denom_type,
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(Weight = ifelse(apply_weights, Weight/sum(Weight), 1)#,
                                                     #B := .[[Response]],
                                                     #A = scales::rescale(A, to=range(B)),
                                                     #B = B
               ) %>%
        summarise(
            SSD = sum(Weight*(.[[Response]] - A)^2),
            ## SSD = sum(Weight*(B-A)^2),
            SST2 = sum(Weight*(.[[Response]] - mean(.[[Response]]))^2),
            SST1 = sum(Weight*(.[[Response]] - mean(A))^2),
            ## SST1 = sum(Weight*(B-mean(A))^2),
            ## SST2 = sum(Weight*(B-mean(B))^2),
            Mean_Weight = mean(Weight)) %>%
        mutate(VE = (1 - (SSD/ifelse(denom_type==1, SST1, SST2)))) %>%
        suppressWarnings() %>%
        suppressMessages()
    
}
DH_model_COR <- function(newdata.full, newdata.sub, Response, apply_weights = FALSE) {
    newdata.full %>% mutate(A := newdata.sub %>% pull(Response),
                            apply_weights = apply_weights, 
                            Weight = 1/ifelse(apply_weights, ifelse(Distance==0, 1, Distance), 1)) %>%
        mutate(Weight = Weight/sum(Weight)) %>%
        summarise(
            COR = cov.wt(cbind(.[[Response]], A), wt = Weight, cor = TRUE)$cor[2,1]
        ) %>%
        suppressWarnings() %>%
        suppressMessages()

}



DH_loop <- function(fixed_sites, free_sites, grid.list,
                    method, method_applies_to, apply_weights, file_attr='') {
    Responses = grid.list$Responses                  # Response name(s) 
    proj.grid = grid.list$proj.grid                  # full projection grid
    full.coords = grid.list$full.coord               # full domain grid
    newdata.full = grid.list$newdata.full            # predictions from full models
    proj.grid.obs = grid.list$proj.grid.obs          # projection grid of full observations
    obs.coords = grid.list$obs.coords                # obsevations grid
    newdata.obs.full = grid.list$newdata.obs.full    # predictions at full observations

    tm <- Sys.time()
    
    
    nK <- nrow(free_sites)                           # for the purpose of the loop size (iterations)
    CONFIGS <- vector('list', length=nK)             # to store the best configs per iteration
    tK <- sum(nrow(free_sites):1)                    # for the purpose of a progress bar
    itK <- 0
    ENERGY <- NULL
    iFree <- NA
    
    ## pb1 <- cli_progress_bar("Cleaning data", total = tK+1)
    
    for (k in 1:nK) {
        print(k)
        if (nrow(free_sites) ==0) return(CONFIGS)
        wch <- 1:nrow(free_sites)

        newdata.sub <- vector('list', length = length(wch))
        names(newdata.sub) <- wch
        IMP <- BEST_VEs <- CONFIG <- METRIC <- newdata.sub
        ENERGY.i <- NULL
        pb2 <- cli_progress_bar("Cleaning data", total = length(wch))
        for (i in wch) {                             # loop over the candidate sites
            iSITE <- free_sites[i,] %>% pull(N) 
            data.sub.sf <- fixed_sites %>% rbind(free_sites[i,])
            newdata.sub[[i]] <- vector('list', length = length(Responses))
            names(newdata.sub[[i]]) <- Responses
            METRIC[[i]] <- newdata.sub[[i]]
            for ( Response in Responses) {
                newdata.sub[[i]][[Response]] <- DH_make_and_fit_model_krige(data.sub.sf, response=Response) %>%
                    mutate(Distance = full.coords$Distance)
                METRIC[[i]][[Response]] <- DH_model_VE(newdata.full, newdata.sub[[i]][[Response]], Response, apply_weights)
                METRIC[[i]][[Response]] <- METRIC[[i]][[Response]] %>%
                    cbind(DH_model_COR(newdata.full, newdata.sub[[i]][[Response]], Response, apply_weights))
            }
            CONFIG[[i]] <- data.sub.sf
            BEST_VE <- METRIC[[i]] %>% map_dbl(~.x$VE)
            BEST_COR <- METRIC[[i]] %>% map_dbl(~.x$COR)
            if (k == 1) {
                imp <- BEST_VE
            } else {
                imp <- BEST_VE - CONFIGS[[k-1]]$BEST_VEs
            }
            imp.wch <- which(imp==max(imp))
            ENERGY.i <- ENERGY.i %>% rbind(data.frame(Free_sites = length(c(na.omit(iFree),i)), iter = i,
                                                      VE = BEST_VE[imp.wch],
                                                      ## VE = mean(BEST_VE),
                                                      COR = BEST_COR[imp.wch],
                                                      ## COR = mean(BEST_COR)))
                                                      Response = Responses[imp.wch]) %>%
                                           cbind(VE=t(BEST_VE)) %>% cbind(COR=t(BEST_COR)))
            BEST_VEs[[i]] <- BEST_VE
            IMP[[i]] <- imp
            itK <- itK + 1
            cli_progress_update(id = pb2)
            ## cli_progress_update(id = pb1)
        }
        iBEST <- which(ENERGY.i$VE == max(ENERGY.i$VE))# %>% filter(VE==max(VE))
        BEST_ENERGY = ENERGY.i[iBEST,]
        BEST_CONFIG = CONFIG[[iBEST]]
        BEST_NEWDATA = newdata.sub[[iBEST]] %>%
            reduce(st_join) %>%
            dplyr::select(-matches('.*\\.([xy]|pred|var)$')) %>%
            suppressMessages() %>%
            suppressWarnings()
        BEST_VE <- BEST_VEs[[iBEST]]
        BEST_IMP <- IMP[[iBEST]]

        fixed_sites <- BEST_CONFIG
        free_sites <- free_sites[-iBEST,]
        iFree = na.omit(c(iFree, iBEST))
        
        CONFIGS[[k]] <- list(
            k = k,
            BEST_ENERGY = BEST_ENERGY,
            BEST_CONFIG = BEST_CONFIG,
            BEST_VEs = BEST_VE,
            BEST_IMP = BEST_IMP,
            BEST_NEWDATA = BEST_NEWDATA,
            BEST_RESP = names(BEST_IMP)[which(BEST_IMP==max(BEST_IMP))]
        )
        ENERGY <- ENERGY %>% rbind(ENERGY.i)
        save(CONFIGS, file=paste0(DATA_PATH, 'processed/CONFIGS',file_attr,'.RData'))
        save(ENERGY, file=paste0(DATA_PATH, 'processed/ENERGY',file_attr,'.RData'))
    }
    tm2 <- Sys.time()
    ## fixed_sites <- BEST_CONFIG
    print(tm2-tm)
}

DH_energy_plot <- function(ENERGY) {
    ENERGY.l <- ENERGY %>% group_by(Free_sites) %>% mutate(wch = max(iter))
    g3 <- ggplot() +
        geom_hline(yintercept=0.8, linetype='dashed') +
        geom_hline(yintercept=0.9, linetype='dashed') +
        geom_line(data=ENERGY.l, aes(y=VE, x=(Free_sites + iter/wch), color='Best')) +
        ## geom_point(data=BEST_ENERGY, aes(y=VE, x=(Free_sites + iter/length(wch)), color=Response), show.legend = TRUE) +
        geom_line(data=ENERGY.l, aes(y=COR, x=(Free_sites + iter/wch)), color='blue') +
        scale_y_continuous(limits=c(0,1)) +
        scale_color_manual('Response', breaks = c(Responses,'Best'), values = scales::hue_pal()(length(Responses)+1)) +
        theme_bw()
    E <- ENERGY.l %>%
        dplyr::select(Free_sites, iter, wch, matches('VE\\.')) %>%
        pivot_longer(cols = matches('VE\\.'), names_prefix='VE.')
    g3 <- g3 + geom_line(data=E, aes(y=value, x=(Free_sites + iter/wch), color=name), alpha=0.5)
    g3
}

DH_energy_plot1 <- function(ENERGY) {
    ## VE - variance explained
    ENERGY.VE <- ENERGY %>%
        dplyr::select(Free_sites, iter, matches('VE\\.')) %>%
        pivot_longer(cols = matches('VE\\.'), names_prefix='VE.')
    ENERGY.VE.resp <- ENERGY.VE %>% group_by(Free_sites, name) %>% summarise(value = mean(value)) %>% suppressMessages()
    ENERGY.VE.mean <- ENERGY.VE.resp %>% group_by(Free_sites) %>% summarise(value = mean(value)) %>% suppressMessages()
    g3a <- ggplot() +
        geom_hline(yintercept=0.8, linetype='dashed') +
        geom_hline(yintercept=0.9, linetype='dashed') +
        geom_line(data=ENERGY.VE.resp, aes(y=value, x=(Free_sites), color=name)) +
        geom_line(data=ENERGY.VE.mean, aes(y=value, x=(Free_sites), color='mean')) +
        scale_x_continuous('Number of free sites in the design') +
        scale_y_continuous('Variance explained (%)', limits=c(0,1), labels = function(x) x*100) +
        scale_color_manual('Response', breaks = c(Responses,'mean'), values = scales::hue_pal()(length(Responses)+1)) +
        theme_bw() +
        theme(legend.position = c(0.99, 0.01), legend.justification = c(1,0))
    ## COR - correlation
    ENERGY.COR <- ENERGY %>%
        dplyr::select(Free_sites, iter, matches('COR\\.')) %>%
        pivot_longer(cols = matches('COR\\.'), names_prefix='COR.')
    ENERGY.COR.resp <- ENERGY.COR %>% group_by(Free_sites, name) %>% summarise(value = mean(value)) %>% suppressMessages()
    ENERGY.COR.mean <- ENERGY.COR.resp %>% group_by(Free_sites) %>% summarise(value = mean(value)) %>% suppressMessages()
    g3b <- ggplot() +
        geom_hline(yintercept=0.8, linetype='dashed') +
        geom_hline(yintercept=0.9, linetype='dashed') +
        geom_line(data=ENERGY.COR.resp, aes(y=value, x=(Free_sites), color=name)) +
        geom_line(data=ENERGY.COR.mean, aes(y=value, x=(Free_sites), color='mean')) +
        scale_x_continuous('Number of free sites in the design') +
        scale_y_continuous('Correlation', limits=c(0,1)) +
        scale_color_manual('Response', breaks = c(Responses,'mean'), values = scales::hue_pal()(length(Responses)+1)) +
        theme_bw() +
        theme(legend.position = c(0.99, 0.01), legend.justification = c(1,0))
    wrap_plots(g3a, g3b)
}

DH_sites_plot <- function(REGION.sf, dat.sf, CONFIGS, iter=1) {
 g4 <- ggplot() +
        geom_sf(data=REGION.sf, color='grey', fill=NA) +
        geom_sf(data=dat.sf, shape = 21, fill=NA, color='black') +
        ## geom_sf(data=dat.sf %>% filter(Site_Type=='designated'), color = 'blue') +
        geom_sf(data=CONFIGS[[iter]]$BEST_CONFIG %>% filter(Site_Type=='designated'), color = 'blue') +
        geom_sf(data=CONFIGS[[iter]]$BEST_CONFIG %>% filter(Site_Type=='random'), color = 'red') +
     theme_bw()
 g4
}

DH_patterns_plot <- function(REGIONS.sf, Responses, newdata.full, CONFIGS, iter = 1) {
    g1 <- vector('list', length(Responses))
    names(g1) <- Responses
    g2 <- g1
    for (resp in Responses) {
        lims <- c(newdata.full %>% pull(resp), CONFIGS[[iter]]$BEST_NEWDATA %>% pull(resp)) %>% range()
        lims[1] <- floor(lims[1])
        lims[2] <- ceiling(lims[2])
        g1[[resp]] <-
            ggplot() +
            geom_raster(data=newdata.full, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
            geom_sf(data=REGIONS.sf %>% st_union() %>% suppressMessages() %>% suppressWarnings(), fill=NA) +
            scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
            theme_bw() %>%
            suppressWarnings() %>%
            suppressMessages()
        g2[[resp]] <- ggplot() +
            ## geom_sf(data=BEST_NEWDATA, aes(color=Cu)) +
            geom_raster(data=CONFIGS[[iter]]$BEST_NEWDATA, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
            geom_sf(data=REGIONS.sf %>% st_union(), fill=NA) +
            scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
            theme_bw() 
    }

    g <- c(g1, g2)
    wrap_plots(g)
    }



DH_loop <- function(fixed_sites, free_sites, grid.list, method, method_applies_to, apply_weights, iFree = NA, ENERGY = NULL, CONFIGS = NULL) {
    ## ensure there are two graphics devices active
    ## while (length(dev.list())<2) {
    ##     dev.new()
    ## }
    ## dev2 <- dev.cur()
    ## dev1 <- dev.prev()

    Responses = grid.list$Responses                  # Response name(s) 
    proj.grid = grid.list$proj.grid                  # full projection grid
    full.coords = grid.list$full.coord               # full domain grid
    newdata.full = grid.list$newdata.full            # predictions from full models
    proj.grid.obs = grid.list$proj.grid.obs          # projection grid of full observations
    obs.coords = grid.list$obs.coords                # obsevations grid
    newdata.obs.full = grid.list$newdata.obs.full    # predictions at full observations
     
    if (is.null(CONFIGS)) { 
        CONFIGS <- list()
    }
    k = length(CONFIGS)+1
    
    if (nrow(free_sites) == 0 ) return(1) 
    wch <- 1:nrow(free_sites)
    ## print(wch)
    ## iFixed = which(dat.sf$Site_Type == 'designated')
    
    newdata.sub <- vector('list', length=length(wch))
    names(newdata.sub) <- wch
    IMP <- CONFIG <- METRIC <- newdata.sub
    ENERGY.l <- NULL
    for (i in wch) {
        iSITE <- free_sites[i,] %>% pull(N) 
        ## cat(paste0('Exploring the addition of the following free sites: ', paste(paste0(SITES_TO_REMOVE, collapse=','), iSITE, sep = ','), ' '))
        data.sub.sf <- fixed_sites %>% rbind(free_sites[i,])#dat.sf[c(iFixed, na.omit(iFree), i),]
        newdata.sub[[i]] <- vector('list', length = length(Responses))
        names(newdata.sub[[i]]) <- Responses
        METRIC[[i]] <- newdata.sub[[i]]
        for ( Response in Responses) {
            newdata.sub[[i]][[Response]] <- DH_make_and_fit_model_krige(data.sub.sf, response=Response) %>%
                mutate(Distance = full.coords$Distance)
            METRIC[[i]][[Response]] <- DH_model_VE(newdata.full, newdata.sub[[i]][[Response]], Response, apply_weights)
            METRIC[[i]][[Response]] <- METRIC[[i]][[Response]] %>%
                cbind(DH_model_COR(newdata.full, newdata.sub[[i]][[Response]], Response, apply_weights))
        }
        CONFIG[[i]] <- data.sub.sf
        BEST_VE <- METRIC[[i]] %>% map_dbl(~.x$VE)
        BEST_COR <- METRIC[[i]] %>% map_dbl(~.x$COR)
        if (k == 1) {
            imp <- BEST_VE
        } else {
            imp <- BEST_VE - CONFIGS[[k-1]]$BEST_IMP
        }
        imp.wch <- which(imp==max(imp))
        ENERGY.l <- ENERGY.l %>% rbind(data.frame(Free_sites = length(c(na.omit(iFree),i)), iter = i,
                                                  VE = BEST_VE[imp.wch],
                                                  ## VE = mean(BEST_VE),
                                                  COR = BEST_COR[imp.wch],
                                                  ## COR = mean(BEST_COR)))
                                                  Response = Responses[imp.wch]) %>%
                                       cbind(VE=t(BEST_VE)) %>% cbind(COR=t(BEST_COR)))
        IMP[[i]] <- BEST_VE #imp
    }
    iBEST <- which(ENERGY.l$VE == max(ENERGY.l$VE))# %>% filter(VE==max(VE))
    BEST_ENERGY = ENERGY.l[iBEST,]
    BEST_CONFIG = CONFIG[[iBEST]]
    BEST_NEWDATA = newdata.sub[[iBEST]] %>%
        reduce(st_join) %>%
        dplyr::select(-matches('.*\\.([xy]|pred|var)$')) %>%
        suppressMessages() %>%
        suppressWarnings()
    ENERGY <- ENERGY %>% rbind(ENERGY.l)
    BEST_IMP <- IMP[[iBEST]]
    CONFIGS[[k]] <- list(
        k = k,
        BEST_ENERGY = BEST_ENERGY,
        BEST_CONFIG = BEST_CONFIG,
        #BEST_VE = BEST_VE
        BEST_IMP = BEST_IMP,
        BEST_RESP = names(BEST_IMP)[which(BEST_IMP==max(BEST_IMP))]
        )
    
    ## dev.set(dev1)
    g3 <- ggplot() +
        geom_hline(yintercept=0.8, linetype='dashed') +
        geom_hline(yintercept=0.9, linetype='dashed') +
        geom_line(data=ENERGY, aes(y=VE, x=(Free_sites + iter/length(wch)), color='Best')) +
        geom_point(data=BEST_ENERGY, aes(y=VE, x=(Free_sites + iter/length(wch)), color=Response), show.legend = TRUE) +
        geom_line(data=ENERGY, aes(y=COR, x=(Free_sites + iter/length(wch))), color='blue') +
        scale_y_continuous(limits=c(0,1)) +
        scale_color_manual('Response', breaks = c(Responses,'Best'), values = scales::hue_pal()(length(Responses)+1)) +
        theme_bw()
    ## coo <- seq(1,length=length(Responses), by = -0.05)
    ## names(coo) <- Responses
    g3 <- g3 + geom_richtext(data=NULL, aes(y=1, x=1,
                                            label = paste0(names(BEST_IMP), '=', round(BEST_IMP,3), collapse='<br>')), vjust=1,hjust = 0) +
        geom_line(data=ENERGY, aes(y=VE, x=(Free_sites + iter/length(wch))), color='black') 
    ## r <- Responses[1]
    ## g3 <- g3 + geom_line(data=ENERGY, aes(y=VE.Cu, x=(Free_sites + iter/length(wch)), color=!!sym(r)))
    ## r <- Responses[2]
    ## g3 <- g3 + geom_line(data=ENERGY, aes(y=VE.Zn, x=(Free_sites + iter/length(wch)), color=!!sym(r)))
    ## r <- Responses[3]
    ## g3 <- g3 + geom_line(data=ENERGY, aes(y=VE.Pb, x=(Free_sites + iter/length(wch)), color=!!sym(r)))
    E <- ENERGY %>%
        dplyr::select(Free_sites, iter, matches('VE\\.')) %>%
        pivot_longer(cols = matches('VE\\.'), names_prefix='VE.')
    g3 <- g3 + geom_line(data=E, aes(y=value, x=(Free_sites + iter/length(wch)), color=name))

                                        #opar <- par(mfrow = c(1,2))
    #plot(VE ~ I(Free_sites + iter/length(wch)), data=ENERGY, type='l', ylim=c(0,1))
    #points(VE ~ I(Free_sites + iter/length(wch)), data = BEST_ENERGY, col='red', pch=16)
    #lines(COR ~ I(Free_sites + iter/length(wch)), data = ENERGY, col='blue', pch=16)

    ## plot(st_geometry(IN.sf), border='gray')
    ## plot(st_geometry(dat.sf), add=TRUE)
    ## plot(st_geometry(data.sub.sf %>% filter(Site_Type=='designated')), add=TRUE, pch=16, col='blue')
    ## plot(st_geometry(data.sub.sf %>% filter(Site_Type=='random')), add=TRUE, pch=16, col='red')
    
    ## plot(st_geometry(IN.sf), border = 'gray')
    ## plot(st_geometry(dat.sf), add=TRUE)
    ## plot(st_geometry(data.sub.sf %>% filter(Site_Type=='designated')), add=TRUE, pch=16, col='blue')
    ## plot(st_geometry(BEST_CONFIG %>% filter(Site_Type=='random')), add=TRUE, pch=16, col='red')
    fixed_sites <- BEST_CONFIG
    free_sites <- free_sites[-iBEST,]
    iFree = na.omit(c(iFree, iBEST))
    
    g4 <- ggplot() +
        geom_sf(data=IN.sf, color='grey', fill=NA) +
        geom_sf(data=dat.sf, shape = 21, fill=NA, color='black') +
        geom_sf(data=data.sub.sf %>% filter(Site_Type=='designated'), color = 'blue') +
        geom_sf(data=BEST_CONFIG %>% filter(Site_Type=='random'), color = 'red') +
        theme_bw()
    
    ## dev.set(dev2)
    ## g1 <- ggplot() +
    ##     geom_raster(data=newdata.full %>%
    ##                     pivot_longer(cols=Responses, names_to='Response', values_to='Value'),
    ##                 aes(y=Latitude, x=Longitude, fill=Value)) +
    ##     geom_sf(data=IN.sf %>% st_union(), fill=NA) +
    ##     scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=c(0,80)) +
    ##     facet_wrap(~Response) +
    ##     theme_bw() 
    g1 <- vector('list', length(Responses))
    names(g1) <- Responses
    g2 <- g1
    for (resp in Responses) {
        lims <- c(newdata.full %>% pull(resp), BEST_NEWDATA %>% pull(resp)) %>% range()
        lims[1] <- floor(lims[1])
        lims[2] <- ceiling(lims[2])
        g1[[resp]] <-
            ggplot() +
            geom_raster(data=newdata.full, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
            geom_sf(data=IN.sf %>% st_union() %>% suppressMessages() %>% suppressWarnings(), fill=NA) +
            scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
            theme_bw() %>%
            suppressWarnings() %>%
            suppressMessages()
        g2[[resp]] <- ggplot() +
            ## geom_sf(data=BEST_NEWDATA, aes(color=Cu)) +
            geom_raster(data=BEST_NEWDATA, aes(y=Latitude, x=Longitude, fill=!!sym(resp))) +
            geom_sf(data=IN.sf %>% st_union(), fill=NA) +
            scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
            theme_bw() 
    }

    g <- c(g1, g2)
    suppressWarnings(print((g3 + g4) /(wrap_plots(g))))
    ## print((g3 + g4) /(g1 + g2))    
    ## print(g4)
    ## print(nrow(ENERGY))
    ## Sys.sleep(10)
    ## print(BEST_CONFIG %>% filter(Site_Type=='random'))
    return(
        DH_loop(
            fixed_sites,
            free_sites,
            ## dat.sf,
            grid.list = list(Responses = Responses,
                             proj.grid = proj.grid,
                             full.coords = full.coords,
                             newdata.full = newdata.full,
                             proj.grid.obs = proj.grid.obs,
                             obs.coords = obs.coords,
                             newdata.obs.full = newdata.obs.full),
            method = method,
            method_applies_to = method_applies_to,
            apply_weights = apply_weights,
            iFree = iFree,
            ENERGY = ENERGY,
            CONFIGS = CONFIGS
        )
    )
}

DH_loss_loop <- function(dat.sf,
                         grid.list,
                         method,                     # information loss method ('MSE', ...)
                         method_applies_to,          # apply method to full domain grid ('full grid') or observations ('obs')
                         apply_weights               # apply spatial weights (Distance) to metric calculations
                         ) { 
 
    Response = grid.list$Response                    # Response name 
    proj.grid = grid.list$proj.grid                  # full projection grid
    full.coords = grid.list$full.coord               # full domain grid
    newdata.full = grid.list$newdata.full            # predictions from full models
    proj.grid.obs = grid.list$proj.grid.obs          # projection grid of full observations
    obs.coords = grid.list$obs.coords                # obsevations grid
    newdata.obs.full = grid.list$newdata.obs.full    # predictions at full observations

    
    if (length(which(dat.sf$Site_Type =='random')) == 0 ) return(1)
    wch <- which(dat.sf$Site_Type == 'random')
    print(wch)
    
    obs.sub.coords <- DH_obs_coords(dat.sf)
    newdata.sub <- vector('list', length=length(wch))
    names(newdata.sub) <- wch
    ## MSE.obs <- MSE <- newdata.obs.sub <- newdata.sub #vector('list', length=length(wch))
    METRIC <- newdata.sub
    ## MSE <- vector('list', length = length(wch))
    ## MSE.obs <- vector('list', length = length(wch))
    for (i in wch) {
        ## i <- wch[i]  
        iSITE <- dat.sf[i,] %>% pull(N) 
        cat(paste0('Exploring the removal of sites: ', paste(paste0(SITES_TO_REMOVE, collapse=','), iSITE, sep = ','), ' '))
        data.sub.sf <- dat.sf[-i,]
        ## mod.sub <- DH_make_and_fit_model(mesh, obs.sub.coords[-i,], data.sub.sf, response=Response)
        ## newdata.sub[[i]] <- DH_model_predictions(proj.grid, mod.sub, full.coords, Response)
        newdata.sub[[i]] <- DH_make_and_fit_model_krige(data.sub.sf, response=Response) %>%
            mutate(Distance = full.coords$Distance)
        if (method_applies_to == 'full grid') {  ## Full grid
            METRIC[[i]] <- DH_model_MSE(newdata.full, newdata.sub[[i]], Response, apply_weights) %>%
                mutate(Obs = i)
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_PERC(newdata.full, newdata.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_SS(newdata.full, newdata.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_VE(newdata.full, newdata.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_COR(newdata.full, newdata.sub[[i]], Response, apply_weights))
        } else {                     ## Observations only
            newdata.obs.sub[[i]] <- DH_model_predictions(proj.grid.obs, mod.sub, obs.coords, Response, apply_weights)
            METRIC[[i]] <- DH_model_MSE(newdata.obs.full, newdata.obs.sub[[i]], Response) %>%
                mutate(Obs = i) 
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_PERC(newdata.obs.full, newdata.obs.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_SS(newdata.obs.full, newdata.obs.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_VE(newdata.obs.full, newdata.obs.sub[[i]], Response, apply_weights))
            METRIC[[i]] <- METRIC[[i]] %>%
                cbind(DH_model_COR(newdata.obs.full, newdata.obs.sub[[i]], Response, apply_weights))
        }
        cat(paste0('Metric:', METRIC[[i]]$VE),'\n')
        
    }
    METRIC1 <- do.call('rbind', METRIC)
    ## MSE1.obs <- do.call('rbind', MSE.obs)

    ## Assess which is obs is associated with the lowest loss of information
    if (method == 'MSE') METRIC <- METRIC1 %>% filter(MSE == min(MSE))
    if (method == 'RMSE') METRIC <- METRIC1 %>% filter(RMSE == min(RMSE))
    if (method == 'PERC') METRIC <- METRIC1 %>% filter(PERC == min(PERC))
    if (method == 'SS') METRIC <- METRIC1 %>% filter(SS == min(SS))
    if (method == 'VE') METRIC <- METRIC1 %>% filter(VE == max(VE))
    if (method == 'COR') METRIC <- METRIC1 %>% filter(COR == max(COR))
    site_with_min_loss <- METRIC %>% pull(Obs)
    iSITE <- dat.sf[site_with_min_loss,] %>% pull(N) 
    ## Record this site number in global vector
    SITES_TO_REMOVE <<- c(SITES_TO_REMOVE, iSITE)
    ## MSE1.obs %>% filter(MSE == min(MSE)) %>% pull(Obs)
    
    RESULTS_LIST[[length(RESULTS_LIST)+1]] <<- list(dat = dat.sf[-site_with_min_loss,],
                                                   iSITE = iSITE,
                                                   newdata.sub = newdata.sub[[site_with_min_loss]],
                                                   metric = METRIC,
                                                   method = method,
                                                   method_applies_to = method_applies_to,
                                                   apply_weights = apply_weights)
    return(
        DH_loss_loop(dat.sf[-site_with_min_loss,],
                     grid.list = list(Response = Response,
                                      proj.grid = proj.grid,
                                      full.coords = full.coords,
                                      newdata.full = newdata.full,
                                      proj.grid.obs = proj.grid.obs,
                                      obs.coords = obs.coords,
                                      newdata.obs.full = newdata.obs.full),
                     method = method,
                     method_applies_to = method_applies_to,
                     apply_weights = apply_weights
        )
    )
    
}












##############################################################################################
## The following functions are taken from:
## source('https://haakonbakka.bitbucket.io/functions-barriers-dt-models-march2017.R')
## They are reproduced here incase the aforementioned dissappears..
##############################################################################################
source('../scripts/functions-barriers-dt-models-march2017.R')

fitINLA.barriermodel = function(bndry, data, var, mesh.type='points',prior.range = c(1, .5), prior.sigma = c(3, 0.01), max.edge=0.02) {

    data = data %>% mutate_(.dots=setNames(var, 'Value'))
    #max.edge=2000 #max.edge=0.02
    bound.outer=0.01
    if (mesh.type=='boundary') {
        mesh = inla.mesh.2d(boundary = bndry,
                            loc=cbind(data$Longitude, data$Latitude),
                            #max.edge = c(1,2)*max.edge,
                            max.edge = c(0.5,1)*max.edge,
                            cutoff = 0.005,
                            offset = c(max.edge, bound.outer))
    } else if (mesh.type=='points') {
        mesh = inla.mesh.2d(loc=cbind(data$Longitude, data$Latitude),
                            max.edge = c(0.5,1)*max.edge)#,
                            #cutoff = 0.0005,
                            #offset = c(max.edge, bound.outer*4))
    }
    ## plot(mesh, main="Our mesh", lwd=0.5)
    ## points(data$Longitude,data$Latitude, col="red")
    ## lines(bndry, col='black',cex=3)

    A.i.s = inla.spde.make.A(mesh, loc=cbind(data$Longitude, data$Latitude))
    stk = inla.stack(data=list(y=data$Value), 
                     effects=list(s=1:mesh$n,
                                  m = rep(1, nrow(data))),
                     A=list(A.i.s, 1),
                     remove.unused = FALSE, tag='est')

    ## The stationary model
    #prior.range = c(1, .5)
    #prior.sigma = c(3, 0.01)
    spde = inla.spde2.pcmatern(mesh, prior.range=prior.range, prior.sigma=prior.sigma)
    #spde = inla.spde2.matern(mesh)
    hyper.iid = list(prec = list(prior='pc.prec', param=prior.sigma))

    mesh = dt.mesh.addon.posTri(mesh)
    ## - compute the triangle positions
    posTri = SpatialPoints(mesh$posTri)
    proj4string(posTri) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84'
    #proj4string(posTri) <- '+proj=utm +zone=52 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
    normal = over(bndry, posTri, returnList=T)
                                        # - checking which mesh triangles are inside the normal area
    normal = unlist(normal)
    Omega = dt.Omega(list(normal, 1:mesh$t), mesh)
    Omega.SP = dt.polygon.omega(mesh, Omega)
    #plot(Omega.SP[[2]], col="grey", main="The barrier region (in grey)")

    Q.barrier = dt.create.Q(mesh, Omega, 
                            fixed.ranges = c(NA, 0.5))
                                        # - We fix the barrier range to a different value than we 
                                        #   used for simulations
                                        # - - Why? It does not matter, as long as it is 'small' 
                                        #     the models are very
                                        #     similar
                                        # - - This shows that you do not need to know the 
                                        #     true 'barrier range'!
                                        # - time: Ca 1 min

    log.prior = dt.create.prior.log.exp(
        prior.param = c(-log(prior.sigma[2])/prior.sigma[1], -log(prior.range[2])/prior.range[1])) 
                                        #c(-log(0.01)/3, -log(0.5)*6))
                                        # - The prior parameters are the lambdas in the exponential 
                                        #   priors for standard deviation and inverse-range
                                        # - the first is log(prob)/exceed, the second log(prob)*exceed
                                        # - the second is exponential for inverse range, therefore multiplication!

    barrier.model = dt.inla.model(
        Q = Q.barrier, log.prior=log.prior)

    mod = list(shortname="barrier-model")
    mod$formula = y~ -1+m + f(s, model=barrier.model)

    stack=stk

    ## Running all the models
    ## Initial values
                                        # - speeds up computations
                                        # - improves accuracy of computations
                                        # - set these to NULL the first time you run a model
    mod$init = c(NULL, NULL, NULL)

    mod$res = inla(mod$formula,
                   data=inla.stack.data(stack),
                   control.predictor=list(A=inla.stack.A(stack), link=1, compute=TRUE),
                   family="gamma", 
                   control.family = list(hyper = hyper.iid, link='log'),
                   control.inla= list(int.strategy = "eb"),
                   control.mode=list(restart=T, theta=mod$init))
    #summary(mod$res)

    #print(paste(round(mod$res$internal.summary.hyperpar$mode, 3), collapse = ','))

    field = mod$res$summary.random$s$mean + mod$res$summary.fixed['m', 'mean']
    xlim = bndry@bbox[1, ] 
    ylim = bndry@bbox[2, ]
    proj = inla.mesh.projector(mesh, xlim = xlim, 
                               ylim = ylim, dims=c(300, 300))
    field.proj = inla.mesh.project(proj, field)
    zlim = range(field.proj, na.rm=TRUE)

    ## image.plot(list(x = proj$x, y=proj$y, z = field.proj), 
    ##            xlim = xlim, ylim = ylim, asp=1) 
    ## contour(x = proj$x, y=proj$y,
    ##         z = field.proj,
    ##         levels=seq(zlim[1], zlim[2],length.out = 10),
    ##         add=TRUE, drawlabels=F, col="white")
    ## plot(Omega.SP[[2]], add=T, border="black", col="white")
    ## points(data$Longitude, data$Latitude)
    coords.grid = as.matrix(expand.grid(Longitude=proj$x, Latitude=proj$y))
    field.proj = inla.mesh.projector(mesh, loc=coords.grid)
    newdata = data.frame(coords.grid, fit = exp(inla.mesh.project(field.proj, field)))
    return(list(fit=newdata, mesh=mesh, Omega.SP=Omega.SP, mod=mod))
}


length_width_ratios <- function(base_width = 1.2, ob.sf, plot_dims, g, ncol=NULL) {
    plot.ratio <- 1/bb_aspect_ratio(st_bbox(ob.sf))
    ## plot_dims <- ggplot2::wrap_dims(length(unique(fit$Date)), ncol = ncol) 
    if (any(class(g) %in% c('gg','ggplot'))) g <- ggplotGrob(g)
    bh <- grid::convertHeight(sum(g$heights), unitTo='in', valueOnly=TRUE)
    bw <- grid::convertWidth(sum(g$widths), unitTo='in', valueOnly=TRUE)

    plot_width = (base_width * plot_dims[1]) + bw
    plot_height = (base_width*plot.ratio*plot_dims[2]) + bh
    list(base_width = base_width, plot.ratio = plot.ratio[[1]], plot_dims = plot_dims, bh = bh, bw = bw,
         fullwidth = plot_width[[1]], fullheight = plot_height[[1]])
}
