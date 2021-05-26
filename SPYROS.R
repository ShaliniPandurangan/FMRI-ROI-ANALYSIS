load("/Users/usi/Desktop/Data_files/residuals/RESI1.Rdata")
load("/Users/usi/Desktop/Data_files/Yfit_diagnosed.Rdata")
load("/Users/usi/Desktop/Data_files/BIC_NETWORK_DIAG_NO_WEIGHTS.Rdata")
load("/Users/usi/Desktop/Data_files/model_results_diag_nw_bic.Rdata")
load("/Users/usi/Desktop/04-Data/fMRI-ROI-time-series/fMRI-ROI-time-series.RData")
source('~/Desktop/Data_files/functions.R')

YFITD= data.frame(Yfit_diagnosed)
YFITD$ID = rep(1:7, each=404)
YFITD = data.frame(do.call(rbind,lapply(split(YFITD, YFITD$ID), function(x){
  #if (nrow(x)==404) return(x[-1,]) else return(rbind(x[-1,],rep(NA,71)))
  return(x[-1,]) 
} )))
rownames(YFITD) <- NULL

RESIDUALS_fixed = do.call(rbind,lapply(split(data.frame(BIC_NETWORK_DIAG_RESI), c(rep(1:4,each=403),rep(5,402),rep(6:7,each=403))), function(y){
  if (nrow(y)==403) return(y) else return(rbind(y,rep(NA,70)))
}))

diagnostics(x=RESIDUALS_fixed, data.fit=YFITD, AR.MAT= t(BIC_NETWORK_DIAG_NO_WEIGHTS), lm.fit=model_results_diag_nw_bic)
