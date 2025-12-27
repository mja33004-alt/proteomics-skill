normalize(iMat)<- function(iMat=iMat){
  source(limma)
  source(qvalue)
  source(MASS)
  source(preprocessCore)
  source(vsn)
  # Normalize so all medians are the same
  meds<-apply(iMat, 2, median, na.rm=TRUE)
  MnMat<-sweep(iMat, 2, meds/mean(meds), FUN="/")
  
  # Quantile Normalization
  QnMat <- 2^normalize.quantiles(log2(iMat))
  # Cyclic Loess Normalization
  CLnMat <- 2^normalizeCyclicLoess(data.frame(log2(iMat)), method = "pairs")
  # Lowess Normalization
  LFnMat <- 2^normalizeCyclicLoess(data.frame(log2(iMat)))
  # Variance Stabilization
  VSNnMat <- 2^justvsn(iMat)
  #meanSdPlot(VSN)
  
  #Sample Loading Normalization
  ## Do SL normalizations then create the reference vector
  dMat <- data.frame(iMat)
  #Separate the TMT data by experiment
  exp1 <- dplyr::select(dMat, ends_with(".1"))
  exp2 <- dplyr::select(dMat, ends_with(".2"))
  exp3 <- dplyr::select(dMat, ends_with(".3"))
  exp4 <- dplyr::select(dMat, ends_with(".4"))
  
  # calculate the global target scaling value
  target <- mean(c(colSums(exp1, na.rm=TRUE), colSums(exp2, na.rm=TRUE), 
                   colSums(exp3, na.rm=TRUE), colSums(exp4, na.rm=TRUE)))
  
  # do the sample loading normalizations (scale to the target value)
  norm_facs <- target / colSums(exp1, na.rm=TRUE)
  exp1 <- sweep(exp1, 2, norm_facs, FUN = "*")
  norm_facs <- target / colSums(exp2, na.rm=TRUE)
  exp2 <- sweep(exp2, 2, norm_facs, FUN = "*")
  norm_facs <- target / colSums(exp3, na.rm=TRUE)
  exp3 <- sweep(exp3, 2, norm_facs, FUN = "*")
  norm_facs <- target / colSums(exp4, na.rm=TRUE)
  exp4 <- sweep(exp4, 2, norm_facs, FUN = "*")
  data <- cbind(exp1, exp2, exp3, exp4)
  
  # Do TMM normalization
  # see exactly what TMM does with SL data
  # data_sl <- na.omit(data) # have to throw out NAs for calcNormFactors which is not ideal for later analysis
  # sl_tmm <- calcNormFactors(data_sl, na.rm=TRUE)
  # data_sl_tmm <- sweep(data_sl, 2, sl_tmm, FUN = "/")
  
  SLnMat <- as.matrix(data)
  
  Return(nMats =
    list(c(MnMat = MnMat, QnMat = QnMat, CLnMat = CLnMat, 
           LFnMat = LFnMat, VSNnMat = VSNnMat, SLnMat = SLnMat)))
}