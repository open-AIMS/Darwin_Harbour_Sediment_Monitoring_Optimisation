source('../scripts/DH_00_config.R')

cli_alert('Starting analyses')

METHOD <<- 'krige'  # gam, spde
RESPONSES <<- c('Cu', 'Zn', 'Pb')
WEIGHTED <<- FALSE
DRY_RUN <<- FALSE

## ---- fitLoadData
load(file = paste0(DATA_PATH, 'processed/data.EA.sf.RData'))
load(file = paste0(DATA_PATH, 'processed/data.MWA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/EA.sf.RData'))
load(file = paste0(DATA_PATH, 'primary/MWA.sf.RData'))
## ----end

## ---- REGION_loop
{
    for (REGION in c('IN', 'OH')[1]) {
        REGION <<- REGION
        ## ---- NORMALIZED_loop
        for (NORMALISED in c(FALSE, TRUE)[2]) {
            NORMALISED <<- NORMALISED
            if (NORMALISED) {
                RESPONSES <<- c('Cu.norm', 'Zn.norm', 'Pb.norm')
                NORMALISED_BY <<- ifelse(REGION=='IN', 'Al', 'Fe')
            } else {
                RESPONSES <<- c('Cu', 'Zn', 'Pb')
            }
            source('DH_40_fitModels.R')
        }
        ## ----end
    }
}
## ----end


