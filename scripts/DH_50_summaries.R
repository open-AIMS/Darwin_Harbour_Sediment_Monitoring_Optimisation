source('../scripts/DH_00_config.R')

load(file = paste0(DATA_PATH, 'processed/data.EA.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/EA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.IN.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.OH.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/IN.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/OH.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.IN.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.IN.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.OH.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/full.coords.OH.RData'))

## ---- functions
## https://learnui.design/tools/data-color-picker.html#palette
DH_fit_model_for_all_responses_krige <- function(x, responses) {
    lapply(responses, function(r) DH_make_and_fit_model_krige(x, response=r)) %>%
        reduce(st_join) %>%
        dplyr::select(-matches('.*\\.([xy]|pred|var)$')) %>%
        suppressMessages() %>%
        suppressWarnings()
}

DH_VE_for_all_responses <- function(newdata.full, newdata.sub, responses, apply_weights) {
    #print(apply_weights)
    VE <- lapply(responses, function(r) {
        DH_model_VE(newdata.full %>% st_drop_geometry(),
                 newdata.sub %>% st_drop_geometry(),
                 Response=r, apply_weights = apply_weights, denom_type=1)$VE
    }
    ) %>%
        reduce(cbind) %>%
        as.data.frame() %>%
        bind_cols(Mean = rowMeans(.)) %>%
        suppressMessages() %>%
        suppressWarnings()
    colnames(VE) <- c(responses,'Mean')
    VE
}
DH_singleRasterMap <- function(data.df, response, REGION.sf, lims = NULL) {
    xbreaks <- REGION.sf %>% st_bbox() %>% `[`(c(1,3)) %>% 
        scales::breaks_width(width=0.1)()
    ybreaks <- REGION.sf %>% st_bbox() %>% `[`(c(2,4)) %>% 
        scales::breaks_width(width=0.1)()
    ggplot() +
        geom_raster(data=data.df, aes(y=Latitude, x=Longitude, fill=!!sym(response))) +
        scale_fill_gradientn(colors=fields::tim.colors(64), trans=scales::pseudo_log_trans(), limits=lims) +
        geom_sf(data=REGION.sf, fill=NA) + 
        scale_x_continuous(breaks = xbreaks) +
        scale_y_continuous(breaks = ybreaks) +
        coord_sf(expand=FALSE) +
        theme_bw() +
        theme(
            axis.title = element_blank(),
            axis.text = element_text(margin =margin(l='2', unit='lines')))
}
DH_dualRasterMap <- function(data.df, newdata.full, response, REGION.sf, nrow = 1) {
    lims <- c(newdata.full %>% pull(response), data.df %>% pull(response)) %>% range()
    lims[1] <- floor(lims[1])
    lims[2] <- ceiling(lims[2])
    ## full data
    g1 <- DH_singleRasterMap(newdata.full, response = response, REGION.sf = REGION.sf, lims = lims)
    ## subset data
    g2 <- DH_singleRasterMap(data.df, response = response, REGION.sf = REGION.sf, lims = lims)
    wrap_plots(g1,g2, nrow = nrow)
}
DH_dualRasterMap_list <- function(data.df, newdata.full, response, REGION.sf, nrow = 1) {
    lims <- c(newdata.full %>% pull(response), data.df %>% pull(response)) %>% range()
    lims[1] <- floor(lims[1])
    lims[2] <- ceiling(lims[2])
    ## full data
    g1 <- DH_singleRasterMap(newdata.full, response = response, REGION.sf = REGION.sf, lims = lims)
    ## subset data
    g2 <- DH_singleRasterMap(data.df, response = response, REGION.sf = REGION.sf, lims = lims)
    list(g1, g2)
}
DH_dualRasterMap_for_all_responses <- function(data.df, newdata.full, responses, REGION.sf, nrow = NULL) {
    g <- lapply(responses, function(r)
        DH_dualRasterMap(data.df, newdata.full, response=r, REGION.sf, nrow = nrow)
        )
    names(g) <- responses
    g
}
DH_dualRasterMap_for_all_responses_list <- function(data.df, newdata.full, responses, REGION.sf, nrow = NULL) {
    g <- lapply(responses, function(r)
        DH_dualRasterMap_list(data.df, newdata.full, response=r, REGION.sf, nrow = nrow)
        )
    names(g) <- responses
    g
}

DH_singleConfigMap <- function(dat.sf, REGION.sf, data.REGION.sf) {
    xbreaks <- REGION.sf %>% st_bbox() %>% `[`(c(1,3)) %>% 
        scales::breaks_width(width=0.1)()
    ybreaks <- REGION.sf %>% st_bbox() %>% `[`(c(2,4)) %>% 
        scales::breaks_width(width=0.1)()
    nFree <- dat.sf %>% filter(Site_Type=='random') %>% nrow()
    ggplot() +
        geom_sf(data=REGION.sf, color='grey', fill=NA) +
        ## geom_sf(data=data.REGION.sf, shape = 21, fill=NA, color='black') +
        ## geom_sf(data=dat.sf %>% filter(Site_Type=='designated'), color = 'blue') +
        ## geom_sf(data=dat.sf %>% filter(Site_Type=='random'), color = 'red') +
        geom_sf(data=data.REGION.sf, shape = 21, fill=NA, aes(color='Unselected')) +
        geom_sf(data=dat.sf %>% filter(Site_Type=='designated'), aes(color = 'Fixed')) +
        geom_sf(data=dat.sf %>% filter(Site_Type=='random'), aes(color = 'Free')) +
        scale_x_continuous(breaks = xbreaks) +
        scale_y_continuous(breaks = ybreaks) +
        scale_color_manual('Sites', labels = c('Fixed', 'Free', 'Unselected'),
                           values=c('#003f5c','#bc5090','#ffa600')) +
        coord_sf(expand=FALSE) +
        theme_bw() +
        theme(axis.title=element_blank()) +
        ## annotate(geom='text', x=-Inf, y=-Inf, label = paste0('Number of free sites: ', nFree), hjust = 0, vjust = 0)
        ggtitle(label = paste0('Number of free sites: ', nFree))
}

## ----end

METHOD <<- 'krige'  # gam, spde
REGION <<- 'IN' #OH
RESPONSES <<- c('Cu', 'Zn', 'Pb')
NORMALISED <<- FALSE

for (REGION in c('IN', 'OH')) {
    print(REGION)
    REGION <<- REGION
    for (NORMALISED in c(FALSE, TRUE)) {
        print(NORMALISED)
        NORMALISED <<- NORMALISED
        if (NORMALISED) {
            RESPONSES <<- c('Cu.norm', 'Zn.norm', 'Pb.norm')
            NORMALISED_BY <<- ifelse(REGION=='IN', 'Al', 'Fe')
        } else {
            RESPONSES <<- c('Cu', 'Zn', 'Pb')
        }
        ## ---- prepareDatas
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

        norm_string <- 'raw'
        if (NORMALISED) {
            data.REGION.sf <- data.REGION.sf %>%
                filter(Normalise_by == NORMALISED_BY)
            data.REGION.df <- data.REGION.df %>%
                filter(Normalise_by == NORMALISED_BY) %>%
                droplevels()
            norm_string <- 'norm'
        }
        ## ----end
        ## ---- fit_full_model_for_all_responses
        newdata.full <- DH_fit_model_for_all_responses_krige(data.REGION.sf,
                                                             responses = RESPONSES)
        full.config <- DH_singleConfigMap(data.REGION.sf,
                                          REGION.sf = REGION.sf,
                                          data.REGION.sf = data.REGION.sf)
        ## ----end
        ## ---- load and compile spsann results
        files <- list.files(path = paste0(DATA_PATH, 'processed'),
                            pattern = paste0('^res:.*__',REGION,'___',norm_string,'.*rds'), full.names = TRUE)

        files_data <- tibble(files)
        results <- files_data %>%
            mutate(data = map(files, ~ readRDS(.)),
                   points = map(data, function(x) x$points),
                   nFree = map(files, ~gsub('.*res:([^:]*).*','\\1',., perl=TRUE)),
                   nFree = as.numeric(nFree),
                   Weighted = map(files, ~gsub('.*____(.*)\\.rds', '\\1', ., perl = TRUE)),
                   Norm = map(files, ~gsub('.*___(.*)____.*\\.rds', '\\1', ., perl = TRUE)),
                   Label = paste0(nFree, '_', Weighted, '__', Norm)
                   ) %>%
            unnest(nFree, Weighted, Norm, Label, points) %>%
            dplyr::select(-files, -data) %>%
            group_by(nFree, Weighted, Norm, Label) %>%
            nest()
        ## ----end
        ## ---- perform krige on configs
        results_with_krige <- results %>%
            mutate(Sub = map(data, function(x) data.REGION.sf[x$id,])) %>%
            mutate(Krige = map(.x=Sub, ~DH_fit_model_for_all_responses_krige(., responses = RESPONSES)))
        ## ----end
        ## ---- perform VE
        results <- results_with_krige %>%
            mutate(VE = map2(.x=Krige, .y=Weighted,
                             ~DH_VE_for_all_responses(.,
                                                      newdata.full = newdata.full,
                                                      responses = RESPONSES,
                                                      apply_weights=ifelse(.y=='wtd', TRUE, FALSE)
                                                      )
                             )
                   )
        ## Table data
        tab <- results %>%
            ungroup() %>% 
            unnest(VE) %>%
            dplyr::select(-Label, -data, -Sub, -Krige) %>%
            arrange(nFree, Weighted) %>%
            mutate(Region = REGION) %>%
            dplyr::select(Region, Norm, everything())
        write_csv(tab, file=paste0(OUTPUT_PATH, 'tables/tab_',REGION,'__', NORMALISED, '.csv'))
        ## ----end
        ## ---- map_config
        results <- results %>%
           mutate(Config = map(Sub,  
                                ~DH_singleConfigMap(., REGION.sf = REGION.sf, data.REGION.sf = data.REGION.sf)
                                )
                   )
        fdims <- results %>%
            filter(Weighted == 'wtd') %>%
            ungroup %>%
            group_by(nFree) %>%
            count() %>%
            nrow() %>%
            wrap_dims()
        wch <- which(fdims==max(fdims))
        ## make the maximum dimensions the number of rows
        g <- results %>%
            filter(Weighted == 'wtd') %>%
            arrange(nFree) %>%
            pull(Config) %>%
            wrap_plots(guides = 'collect', nrow = fdims[wch]) +
            plot_annotation(tag_levels = 'a', tag_suffix=')') 

        plot.size <- length_width_ratios(base_width = 4.5, REGION.sf, fdims, g) 

        ## plot.ratio <- bb_aspect_ratio(st_bbox(REGION.sf))
        ## base_width <- 4.5
        ## a = results[1,]$Config[[1]]
        ## (bh <- grid::convertHeight(sum(ggplotGrob(a)$heights), unitTo='in', valueOnly=TRUE))
        ## (bw <- grid::convertWidth(sum(ggplotGrob(a)$widths), unitTo='in', valueOnly=TRUE))
        ggsave(paste0(OUTPUT_PATH, 'figures/config_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___wtd.pdf'),
               g,
               width = plot.size$fullwidth, #(base_width*fdims[-wch]) + bw,
               height = plot.size$fullheight #(base_width*plot.ratio + bh)*fdims[wch]
               )
        ggsave(paste0(OUTPUT_PATH, 'figures/config_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___wtd.png'),
               g,
               width = plot.size$fullwidth,
               height = plot.size$fullheight,
               dpi = 100
               )

        fdims <- results %>%
            filter(Weighted == 'unwtd') %>%
            ungroup %>%
            group_by(nFree) %>%
            count() %>%
            nrow() %>%
            wrap_dims()
        ## wch <- which(fdims==max(fdims))
        g <- results %>%
            filter(Weighted == 'unwtd') %>%
            arrange(nFree) %>%
            pull(Config) %>%
            wrap_plots(guides = 'collect', nrow = fdims[wch]) +
            plot_annotation(tag_levels = 'a', tag_suffix=')') 

        plot.size <- length_width_ratios(base_width = 4.5, REGION.sf, fdims, g, ncol=6) 
        ## plot.ratio <- bb_aspect_ratio(st_bbox(REGION.sf))
        ## base_width <- 4.5
        ## a = results[1,]$Config[[1]]
        ## (bh <- grid::convertHeight(sum(ggplotGrob(a)$heights), unitTo='in', valueOnly=TRUE))
        ## (bw <- grid::convertWidth(sum(ggplotGrob(a)$widths), unitTo='in', valueOnly=TRUE))
        ggsave(paste0(OUTPUT_PATH, 'figures/config_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___unwtd.pdf'),
               g,
               width = plot.size$fullwidth, #(base_width*fdims[-wch]) + bw,
               height = plot.size$fullheight #(base_width*plot.ratio + bh)*fdims[wch]
               )
        ggsave(paste0(OUTPUT_PATH, 'figures/config_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___unwtd.png'),
               g,
               width = plot.size$fullwidth,
               height = plot.size$fullheight,
               dpi = 100
               )
        ## Table data
        walk2(.x=results$Sub,
              .y=results$Label, 
              .f = function(df, lab) {
                  df <- df %>%
                      mutate(Longitude = st_coordinates(.)[,1],
                             Latitude = st_coordinates(.)[,2]) %>%
                      st_drop_geometry() %>%
                      dplyr::select(Region, Report_Card_Region, Site, Site_ID, Site_Type, Longitude, Latitude, Distance)
                  save(df, file=paste0(DATA_PATH, 'processed/config_',REGION,'__',lab,'.RData'))
                  write_csv(df, file=paste0(OUTPUT_PATH, 'tables/config_',REGION,'__', lab, '.csv'))
              }
              )
        ## ----end
        ## ---- map_pattern
        results <- results %>%
            mutate(g = map(Krige,  
                           ~DH_dualRasterMap_for_all_responses(., newdata.full,
                                                               response = RESPONSES,
                                                               REGION.sf = REGION.sf)
                           )
                   )

        results <- results %>% mutate(gs = map(g, ~wrap_plots(., ncol=1) +
                                                      plot_annotation(tag_levels = 'a', tag_suffix=')') &
                                                      theme(plot.tag.position=c(0,1),
                                                            plot.tag = element_text(hjust=1, vjust=0))
                                               ))

        plot.size <- length_width_ratios(base_width = 4.5,
                                         REGION.sf,
                                         plot_dims=c(2,length(results[1,]$g[[1]])),
                                         results[1,]$gs[[1]], ncol=6) 
        
        ## plot.ratio <- bb_aspect_ratio(st_bbox(REGION.sf))
        ## base_width <- 4.5
        ## a = results[1,]$g[[1]][[1]]$patches$plots[[1]]
        ## (bh <- grid::convertHeight(sum(ggplotGrob(a)$heights), unitTo='in', valueOnly=TRUE))
        ## (bw <- grid::convertWidth(sum(ggplotGrob(a)$widths), unitTo='in', valueOnly=TRUE))
        ##                                 #bw = 1.554958*5
        ##                                 #bh = 0.5184851*5
        walk2(.x=results$gs,
              .y=results$Label,
              ~ggsave(paste0(OUTPUT_PATH, 'figures/maps_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___',.y,'.pdf'),
                      .x,
                      width = plot.size$fullwidth,
                      height = plot.size$fullheight
                      ## width = base_width*2,
                      ## width = (base_width*2) + bw,
                      ## height = (base_width*length(RESPONSES))*plot.ratio
                      ## height = (base_width*plot.ratio + bh)*length(RESPONSES)
                      )
              )
        walk2(.x=results$gs,
              .y=results$Label,
              ~ggsave(paste0(OUTPUT_PATH, 'figures/maps_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___',.y,'.png'),
                      .x,
                      width = plot.size$fullwidth,
                      height = plot.size$fullheight,
                      ## width = (base_width*2) + bw,
                      ## height = (base_width*plot.ratio + bh)*length(RESPONSES),
                      dpi = 100
                      )
              )
        ## ----end
        ## ---- coplot
        results <- results %>%
            mutate(Coplot = map2(g, Config,
                                 ~wrap_plots(
    wrap_plots(full.config +
               theme(title = element_text(size=rel(0.75)),
                     plot.margin = margin(t=0.5, r=4, unit = 'lines')),
               .y + theme(title = element_text(size=rel(0.75)),
                          plot.margin = margin(t=0.5, r=4, unit = 'lines')),
               ncol = 2, guides = 'collect') &
    theme(legend.position = 'bottom'),
    wrap_plots(.x, ncol=1),
    ncol=1, heights = c(1,3)) + plot_annotation(tag_levels = 'a', tag_suffix=')') &
        theme(plot.tag.position=c(0,1),
              plot.tag = element_text(hjust=1, vjust=0, size=12)
              )))

        plot.size <- length_width_ratios(base_width = 4.5,
                                         REGION.sf,
                                         plot_dims=c(2,length(RESPONSES)+1),
                                         results[1,]$Config[[1]], ncol=NULL) 
        ## plot.ratio <- bb_aspect_ratio(st_bbox(REGION.sf))
        ## base_width <- 4.5
        ## a = results[1,]$Coplot[[1]][[1]]$patches$plots[[1]]
        ## (bh <- grid::convertHeight(sum(ggplotGrob(a)$heights), unitTo='in', valueOnly=TRUE))
        ## (bw <- grid::convertWidth(sum(ggplotGrob(a)$widths), unitTo='in', valueOnly=TRUE))
        ##                                 #bw = 1.554958*5
        ##                                 #bh = 0.5184851*5
        walk2(.x=results$Coplot,
              .y=results$Label,
              ~ggsave(paste0(OUTPUT_PATH, 'figures/coplot_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___',.y,'.pdf'),
                      .x,
                      width = plot.size$fullwidth,
                      height = plot.size$fullheight,
                      ## width = (base_width*2) + bw,
                      ## height = (base_width*plot.ratio + bh)*length(RESPONSES)
                      )
              )
        walk2(.x=results$Coplot,
              .y=results$Label,
              ~ggsave(paste0(OUTPUT_PATH, 'figures/coplot_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___',.y,'.png'),
                      .x,
                      width = plot.size$fullwidth,
                      height = plot.size$fullheight,
                      ## width = (base_width*2) + bw,
                      ## height = (base_width*plot.ratio + bh)*length(RESPONSES),
                      dpi = 100
                      )
              )
        ## ----end
        ## ---- coplot4report

        r1 <- results %>%
            mutate(g2 = map2(.x = Krige, .y = Config,
                             .f = function(.x, .y) {
                                 s1 <- full.config +
                                     guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                                                                 override.aes = list(size = rel(4),
                                                                                     shape = c(19,19,21))))
                                 s2 <- .y +
                                     guides(color = guide_legend(title.position = 'top', title.hjust = 0.5,
                                                                 override.aes = list(size = rel(4),
                                                                                     shape = c(19,19,21))))
                                 sA <- DH_dualRasterMap_for_all_responses_list(.x, newdata.full,
                                                                               response = RESPONSES,
                                                                               REGION.sf = REGION.sf,
                                                                               nrow = 2)
                                 s3 <- sA[[1]][[1]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 s4 <- sA[[1]][[2]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 s5 <- sA[[2]][[1]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 s6 <- sA[[2]][[2]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 s7 <- sA[[3]][[1]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 s8 <- sA[[3]][[2]] +
                                     theme(legend.key.width = unit(1.5, 'cm')) +
                                     guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0.5))
                                 ## wrap_plots(s1,s2,
                                 ##            s3, s4,
                                 ##            s5, s6,
                                 ##            s7, s8,
                                 ##            ncol = 4,
                                 ##            byrow = FALSE) +
                                 ##     plot_annotation(tag_levels = 'a', tag_suffix = ')') &
                                 ##     theme(plot.tag.position = c(0,1),
                                 ##           plot.tag = element_text(hjust = 1, vjust = 0, size=12))
                                 c1 <- wrap_plots(s1, s2, guides = 'collect', nrow=2) &
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = rel(1.5)),
                                           legend.title = element_text(size = rel(1.5)),
                                           ) #+
                                     #guides(color = guide_legend(title.position = 'top', title.hjust = 0.5))
                                 c2 <- wrap_plots(s3, s4, guides = 'collect', nrow=2) &
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = rel(1.5)),
                                           legend.title = element_text(size = rel(1.5)))
                                 c3 <- wrap_plots(s5, s6, guides = 'collect', nrow=2) &
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = rel(1.5)),
                                           legend.title = element_text(size = rel(1.5)))
                                 c4 <- wrap_plots(s7, s8, guides = 'collect', nrow=2) &
                                     theme(legend.position = 'bottom',
                                           legend.text = element_text(size = rel(1.5)),
                                           legend.title = element_text(size = rel(1.5)))
                                 wrap_plots(c1, c2, c3, c4, ncol=4, byrow=FALSE) +
                                     plot_annotation(tag_levels = 'a', tag_suffix = ')') &
                                     theme(plot.tag.position = c(0,1),
                                           plot.tag = element_text(hjust = 1, vjust = 0, size=16))
                             }
                             )
                   )
        r1[1,]$g2
        
        plot.size <- length_width_ratios(base_width = 4.5,
                                         REGION.sf,
                                         plot_dims=c(length(RESPONSES)+1, 2),
                                         r1[1,]$g2[[1]], ncol=NULL) 
        walk2(.x=r1$g2,
              .y=r1$Label,
              ~ggsave(paste0(OUTPUT_PATH, 'figures/coplot_alt_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'___',.y,'.png'),
                      .x,
                      width = plot.size$fullwidth,
                      height = plot.size$fullheight,
                      ## width = (base_width*2) + bw,
                      ## height = (base_width*plot.ratio + bh)*length(RESPONSES),
                      dpi = 100
                      )
              )
        ## ----end
        ## ---- VE_plot
        g <- results %>% unnest(VE) %>%
            mutate(Weighted = forcats::fct_recode(Weighted, Unweighted = 'unwtd', Weighted = 'wtd'))  %>%
            ungroup() %>%
            pivot_longer(cols = c('Mean',RESPONSES,)) %>%
            mutate(name = factor(name, levels = c('Mean',RESPONSES))) %>%
            ggplot() +
            geom_hline(yintercept = 0.95, linetype='dashed') +
            geom_hline(yintercept = 0.9, linetype='dashed') +
            geom_hline(yintercept = 0.8, linetype='dashed') +
            geom_line(aes(y=value, x=nFree, color=name, alpha=name)) +
            ## scale_color_discrete('Response') +
            ## scale_color_manual('Response', values = tim.colors(n=length(c(RESPONSES, 1)))) +
            viridis::scale_color_viridis('Response', discrete=TRUE) +
            ## scale_color_brewer('Response', palette = "Dark2") +
            scale_color_manual('Response', values = c('#003f5c','#7a5195','#ef5675','#ffa600')) +
            scale_alpha_manual('Response', values = c(1, rep(1, length(RESPONSES)))) +
            scale_y_continuous('Variance explained (%)', labels = function(x) x*100) +
            ## scale_x_continuous('Number of free sites added to fixed sites', breaks = unique(results$nFree)) +
            scale_x_continuous('Number of free sites added to fixed sites', breaks = seq(0, max(results$nFree), by=20)) +
            facet_wrap(~Weighted) +
            theme_bw() +
            theme(legend.position = c(0.99, 0.01), legend.justification = c(1,0),
                  strip.background = element_rect(fill="#0078E6"),
                  strip.text = element_text(color = 'white', face = 'bold', size = rel(1.25)),
                  axis.title = element_text(size = rel(1.25)),
                  axis.title.x = element_text(margin=margin(t=1, unit='lines')),
                  axis.title.y = element_text(margin=margin(r=1, unit='lines')))
        g
        ggsave(filename = paste0(OUTPUT_PATH, 'figures/VE_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'.pdf'),
               g, width=8, height=8/GOLDEN_RATIO)
        ggsave(filename = paste0(OUTPUT_PATH, 'figures/VE_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'.png'),
               g, width=8, height=8/GOLDEN_RATIO,
               dpi = 100)

        ## separate for weigthed and unweighted
        df1 <- results %>% unnest(VE) %>%
            mutate(Weighted = forcats::fct_recode(Weighted, Unweighted = 'unwtd', Weighted = 'wtd'))  %>%
            ungroup() %>%
            pivot_longer(cols = c('Mean',RESPONSES,)) %>%
            mutate(name = factor(name, levels = c('Mean',RESPONSES))) %>%
            dplyr::select(-Label, -data, -Sub, -Krige, - Config, -g, -gs, -Coplot) %>%
            group_by(Weighted, Norm) %>%
            nest()
        walk2(.x = df1$data, .y = df1$Weighted,
              .f = function(.x, .y) {
                  g <- .x %>% ggplot() +
                      geom_hline(yintercept = 0.95, linetype='dashed') +
                      geom_hline(yintercept = 0.9, linetype='dashed') +
                      geom_hline(yintercept = 0.8, linetype='dashed') +
                      geom_line(aes(y=value, x=nFree, color=name, alpha=name)) +
                      ## scale_color_discrete('Response') +
                      ## scale_color_manual('Response', values = tim.colors(n=length(c(RESPONSES, 1)))) +
                      viridis::scale_color_viridis('Response', discrete=TRUE) +
                      ## scale_color_brewer('Response', palette = "Dark2") +
                      scale_color_manual('Response', values = c('#003f5c','#7a5195','#ef5675','#ffa600')) +
                      scale_alpha_manual('Response', values = c(1, rep(1, length(RESPONSES)))) +
                      scale_y_continuous('Variance explained (%)', labels = function(x) x*100) +
                      ## scale_x_continuous('Number of free sites added to fixed sites', breaks = unique(results$nFree)) +
                      scale_x_continuous('Number of free sites added to fixed sites', breaks = seq(0, max(results$nFree), by=20)) +
                      ## facet_wrap(~Weighted) +
                      theme_bw() +
                      theme(legend.position = c(0.99, 0.01), legend.justification = c(1,0),
                            strip.background = element_rect(fill="#0078E6"),
                            strip.text = element_text(color = 'white', face = 'bold', size = rel(1.25)),
                            axis.title = element_text(size = rel(1.25)),
                            axis.title.x = element_text(margin=margin(t=1, unit='lines')),
                            axis.title.y = element_text(margin=margin(r=1, unit='lines')))
                  ggsave(filename = paste0(OUTPUT_PATH, 'figures/VE_',REGION,'__',ifelse(NORMALISED, 'NORMALISED', 'RAW'),'_',.y,'.png'),
                         g, width=8, height=8/GOLDEN_RATIO,
                         dpi = 100)
              }
              )        
        ## ----end
    }
}
## END
