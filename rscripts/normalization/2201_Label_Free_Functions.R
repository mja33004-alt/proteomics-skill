freq <- function(myMat, myTitle){
  require(reshape2)
  data.frame(myMat) %>%
    melt() %>%
    ggplot(aes(x = value, color = variable)) + 
    geom_freqpoly(na.rm = TRUE, bins = 50) +
    ggtitle(myTitle) + xlim(0,500) + ylim(0,150) + theme(legend.position = "none") + 
    xlab("Intensity")
}

Comp <- function(myMat, myTitle){
  require(dplyr)
  require(limma)
  require(qvalue)
  resTable<-data.frame(Comparison=c("Syria", "Iraq", "Saudi"),
                       q10=rep(0,3), q30=rep(0,3), estimated_changed=rep(0,3))
  
  OKtest<-apply(myMat, 1, function(x) any(!is.na(x)))
  
  ## Reorder sample table
  sampleTable <- left_join(data.frame(columnName = colnames(myMat), stringsAsFactors = FALSE),
                           sampleTable, 
                           by = c("columnName"))
  all(sampleTable$columnName == colnames(myMat))
  
  design <- model.matrix(~ 0 + sampleTable$ID.3 + 
                           sampleTable$group)
  colnames(design) <- gsub("sampleTable[$]","",colnames(design))
  colnames(design) <- gsub("group", "", colnames(design))
  colnames(design) <- gsub("ID.3", "patient", colnames(design))
  
  contrasts <- makeContrasts(contalt - control, case - control, case - contalt, levels = design)
  #rownames(contrasts) <- c("control", "contalt", "case")
  fit <- lmFit(log2(myMat[OKtest, sampleTable$columnName]), design=design)
  # summary(decideTests(fit))
  fit <- contrasts.fit(fit, contrasts)
  ebout <- eBayes(fit, trend = TRUE)
  
  resTable<-data.frame(Comparison=colnames(contrasts),
                       q10=rep(0,3), q30=rep(0,3), estimated_changed=rep(0,3))
  
  pValues <- list()
  qValues <- list()
  tStatistics <- list()
  qObj <- list()
  for(i in 1:ncol(contrasts)){
    pValues[[i]] <- rep(NA, nrow(myMat))
    pValues[[i]][OKtest] <- ebout$p.value[,i]
    qObj[[i]] <- qvalue(ebout$p.value[,i])
    qValues[[i]] <- rep(NA, nrow(myMat))
    qValues[[i]][OKtest] <- qObj[[i]]$qvalues
    tStatistics[[i]] <- rep(NA, nrow(myMat))
    tStatistics[[i]][OKtest] <- ebout$t[,i]
    
    resTable$q10[i] <- sum(!is.na(qValues[[i]]) & qValues[[i]] < 0.1)
    resTable$q30[i] <- sum(!is.na(qValues[[i]]) & qValues[[i]] < 0.3)
    resTable$estimated_changed[i] <- round(sum(OKtest) * (1 - qObj[[i]]$pi0))
  }
  
  names(pValues) <- colnames(contrasts)
  names(qValues) <- colnames(contrasts)
  names(tStatistics) <- colnames(contrasts)
  
  qqFrame <- data.frame(A = apply(myMat[,sampleTable$group =="control"], 1, mean, na.rm = TRUE),
                        B = apply(myMat[,sampleTable$group == "contalt"], 1, mean, na.rm = TRUE),
                        C = apply(myMat[,sampleTable$group == "case"], 1, mean, na.rm = TRUE)) 
  
  volcanoPlotList <- list()
  MAplotList <- list()
  plotFrameList <- list()
  titleVec <- colnames(contrasts)
  for(i in 1:ncol(contrasts)){
    myGroups <- row.names(contrasts)[contrasts[,i] != 0]
    myFC <- rep(NA, nrow(myMat))
    myFC[OKtest] <- 2^ebout$coefficients[,i]
    plotFrameA <- data.frame(accession = proteins$Accession,
                             description = sapply(strsplit(proteins$Description, " OS"), '[', 1),
                             FC = myFC,
                             Intensity = rowMeans(myMat[,sampleTable$group %in% myGroups], na.rm = TRUE),
                             Abundance = rowMeans(myMat[,sampleTable$group %in% myGroups], na.rm = TRUE),
                             logpValue = -log10(pValues[[i]]),
                             pValue = pValues[[i]],
                             affected = !is.na(qValues[[i]]) & qValues[[i]] < 0.3,
                             qValue = qValues[[i]],
                             logqValue = -log10(qValues[[i]]),
                             tStatistics = tStatistics[[i]])
    plotFrameList[[i]] <- plotFrameA
    volcanoPlotList[[i]] <- ggplot(plotFrameA, aes(x = FC, y = logqValue, text = description, color = affected)) +
      geom_point(na.rm = TRUE) +
      scale_color_manual(values = c("black", "red")) +
      scale_x_continuous(trans = "log2", breaks = 2^seq(-2,2,1), labels = 2^seq(-2,2,1)) +
      ylim(c(0,2)) +
      theme(legend.position = "none") +
      labs(title = myTitle, subtitle = titleVec[i])
    
    MAplotList[[i]] <- ggplot(plotFrameA, aes(x = Abundance, y = FC, text = description, color = affected)) +
      geom_point(na.rm = TRUE) +
      geom_smooth(data = subset(plotFrameA, affected == FALSE), aes(group = affected), formula = y~x, method=stats::loess, fullrange = TRUE) + 
      scale_color_manual(values = c("black", "red")) +
      scale_x_continuous(trans = "log10") +
      #scale_y_continuous(trans = "log2", breaks = 2^seq(-2,2,1), labels = 2^seq(-2,2,1)) +
      ylim(c(0,4)) +
      theme(legend.position = "none")+
      labs(title = myTitle, subtitle = titleVec[i])
  }
  
  return(list(r = resTable, p = pValues, q = qValues, t = tStatistics, 
              v = volcanoPlotList, m = MAplotList, qq = qqFrame, qObj = qObj,
              plotFrameList = plotFrameList))
}


merge <- function(myMat, myComp, myFile){
my <- cbind(myMat, proteins$Accession)
my <- cbind(my, sapply(strsplit(proteins$Description, " OS"), '[', 1))
my <- as.data.frame(my, stringsAsFactors = FALSE)
my <- my %>%
  rename(accession = V52, description = V53)
df <- data.frame(myComp$plotFrameList)

df <- left_join(df, my, by = c("accession" = "accession", "description" = "description"))
write.csv(df, "Excel/Test.csv")
write.csv(my, file = myFile)
}

plotcheck <- function(access, desc){
  my %>% filter(accession == access) %>% 
    select(-c(52,53)) %>%
    pivot_longer(cols = everything()) %>%
    replace_na(list(value = 0)) %>%
    mutate(value = as.numeric(value)) %>%
    left_join(as.data.frame(sampleTable), by = c("name" = "columnName")) %>%
    ggplot(aes(x = ID.3, y = value, fill = as.factor(group))) + 
    geom_histogram(position = position_dodge(), stat = "identity") + 
    scale_fill_discrete(name = "Treatment Group") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    facet_grid(rows = vars(ID.3)) + 
    xlab("Patient ID") + ylab("Reporter Ion Intensity") + ggtitle(desc, access)
}

tableFun <- function(myComp, x){
  contaltcontrol <- as.data.frame(myComp$plotFrameList[1]) %>% 
    select(description, FC, qValue) %>% 
    arrange(qValue) %>% 
    head(x)
  Table1 <- contaltcontrol %>% kable(col.names = c("Protein", "Fold Change", "q-Value"))
  
  casecontrol <- as.data.frame(myComp$plotFrameList[2]) %>% 
    select(description, FC, qValue) %>% 
    arrange(qValue) %>% 
    head(x)
  Table2<- casecontrol %>% kable(col.names = c("Protein", "Fold Change", "q-Value"))
  
  casecontalt <- as.data.frame(myComp$plotFrameList[3]) %>% select(description, FC, qValue) %>% 
    arrange(qValue) %>% 
    head(x)
  Table3 <- casecontalt %>% kable(col.names = c("Protein", "Fold Change", "q-Value"))
  
  dataframes <- list(contaltcontrol, casecontrol, casecontalt)
  tableList <- list(Table1, Table2, Table3)
  return(list(dataframes = dataframes, tableList = tableList))
}

#1 Contalt - Control
#2 Case - Control
#3 Case - Contalt

comptables <- function(mytab1, mytab2){  
  compTable <- data.frame(Comparison = c("Contalt - Control", "Case - Control", "Case - Contalt"), Percent  = rep(0,3))
  for(i in 1:3){
    compTable$Percent[i] <- sum(mytab1$dataframes[[i]][,1] %in% mytab2$dataframes[[i]][,1])
  }
  return(compTable)
}

dfformat <- function(myComp, myNorm){
  ResComp <- myComp$r
  ReturnComp <- ResComp %>% mutate(Normalization = myNorm) %>%
    melt(id.vars = c("Normalization", "Comparison"))
  return(ReturnComp)
}

getProteinTable <- function(x){
  return(dplyr::filter(proteins, geneName %in% strsplit(x, ",")[[1]]))
}
