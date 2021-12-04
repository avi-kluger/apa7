# FOR TESTING
# data = x
# alphaVector <- c(1, 1, 1, 1, .87, .91, .77, .96, .91, .59, .66, .82, .92, .64, .48)
# disattenuated = FALSE
# listwiseDeletion = TRUE
# reliabilities = c(rep(1, ncol(data)))

apa7corr <- function(data = x,
                     useLabels = FALSE,
                     listwiseDeletion = FALSE,
                     disattenuated = FALSE, 
                     reliabilities = c(rep(1, ncol(data))),
                     threeStars = FALSE
                     ) {
  library(tidyverse)
  library(sjlabelled)
  library(sjPlot)
  library(CTT) # for disattenuated correlations
  library(weights)  # for removing leading zeors with rd
  library(Hmisc) # for correlation p values
  
  if(disattenuated == TRUE & 
     sum(reliabilities) == ncol(data)) {
  warning('You requested disattenuated correlations but did not supply your own 
 vector of reliabilities. Disattenuated correlations are NOT printed.')
    disattenuated <-  FALSE
  }
  
  # Prepare a header & footnote
  if(listwiseDeletion == FALSE){
     Col.header  <-  c("Measure", "N", "Mean",	"SD", 1:ncol(data))
  }else{
     Col.header  <-  c("Measure", "Mean",	"SD", 1:ncol(data))
  }
  
  ft.note  <- "<i>Note</i>. "
  ft.N     <- paste0(" <i>N</i> = ", nrow(data), ".")
  
  if(sum(reliabilities) < ncol(data)) {
    ft.star  <- " Values in the diagonal are reliabilities.
                <br>*<i>p</i>< .05. **<i>p</i>< .01."
  }else{
    ft.star  <- "<br>*<i>p</i>< .05. **<i>p</i>< .01."
  }
  if(threeStars == TRUE) 
  ft.star <- paste(ft.star, " ***<i>p</i>< .001.")
  ft.miss  <- " Correlations for variables with missing data are based on
                pairwise deletion."
  ft.above <- " Values above the diagonal are disattenuated correlations."
  ft.diag  <- " Values in the diagonal are reliabilities."

  ft <- paste0(ft.note, ft.star)

  if(listwiseDeletion == TRUE & disattenuated == FALSE){
    ft <- paste0(ft.note, ft.N, ft.miss, ft.star)
  }
  if(listwiseDeletion == TRUE & disattenuated == TRUE){
    ft <- paste0(ft.note, ft.N, ft.miss, ft.above, ft.star)
  }
  if(listwiseDeletion == FALSE & disattenuated == TRUE){
    ft <- paste0(ft.note, ft.N, ft.above, ft.star)
  }   
  
  if(listwiseDeletion == TRUE) data <- all_of(data) %>% drop_na

  # Build a data frame with the descriptive statistics
  set_label(data) <-
    ifelse(get_label(data) == "", names(data), get_label(data))
  
  if(useLabels == TRUE)  {
    corTable          <- data.frame(Measure = get_label(data))
  } else {
    corTable          <- data.frame(Measure = names(data))
  }
  
  corTable$Measure  <- paste(1:nrow(corTable), corTable$Measure, sep = ". ")
  corTable$N        <- sapply(data, function(x) sum(!is.na(x)))
  corTable$Mean     <- sapply(data, mean, na.rm = TRUE)
  corTable$SD       <- sapply(data, sd, na.rm = TRUE)
  
  # Compute correlations in a data that contains missing values
  corMatrix       <- cor(data, use = "pairwise.complete.obs")
  
  # Create a table of significance stars below the diagonal: * and **
  p <- rcorr(as.matrix(data))
  p <- as.matrix(p$P)
  pCorr <- matrix("", nrow = nrow(p), ncol = nrow(p))
  pCorr[p < .05]   <- "*"
  pCorr[p < .01]   <- "**"
  if(threeStars == TRUE) pCorr[p < .001]  <- "***"
  
  #Remove stars from the upper triangle above the diagonal
  pCorr[upper.tri(pCorr)] <- ""
  
  corMatrix         <- disattenuated.cor(corMatrix, unlist(reliabilities))
  corMatrix         <- as.data.frame(corMatrix)
  corMatrix         <- sapply(corMatrix, rd, 2)
  
  if (disattenuated == FALSE) corMatrix[upper.tri(corMatrix)] <- ""
  
  # Fix: drop reliabilities of 1.00 and change correlations of .000000 to .00
  for (i in 1:ncol(corMatrix)) {
    corMatrix[i, i] <- ifelse(corMatrix[i, i] == "1.00", "-", corMatrix[i, i])
  }
  
  # Drop disattenated correlations that are identical to raw correlations
  for (i in 1:ncol(corMatrix)) {
    for (j in 1:ncol(corMatrix))
    corMatrix[i, j] <- ifelse(corMatrix[i, i] == "-" & 
                              corMatrix[j, j] == "-" &  j > i,
                              "", corMatrix[i, j])
  }
  corMatrix <- ifelse(corMatrix == ".000000", ".00", corMatrix)
  
  # Wrap reliabilities with parentheses using a loop
  for (i in 1:ncol(corMatrix))
    corMatrix[i, i] <- paste0("(", corMatrix[i, i], ")")
  
  corMatrix [pCorr == "*"]   <- paste0(corMatrix [pCorr == "*"], "*")
  corMatrix [pCorr == "**"]  <-
    paste0(corMatrix [pCorr == "**"], "**")
  corMatrix [pCorr == "***"] <-
    paste0(corMatrix [pCorr == "***"], "***")
  
  corTable            <- cbind(corTable, corMatrix)
  row.names(corTable) <- NULL
  colnames(corTable)  <- c("Measure", "N", "Mean",	"SD", 1:nrow(corTable))
  if(listwiseDeletion == TRUE)  corTable <- corTable %>% select(-N)
  corTable[, c("Mean", "SD")] <- round(corTable[, c("Mean", "SD")], 2)
  
  tab_df(
    corTable,
    title = "<b>Table X.</b> <br>
       Descriptive Statistics and Correlations for Study Variables",
    col.header = Col.header,
    # use.viewer = FALSE,
    show.footnote = TRUE,
    footnote = ft,
    CSS = list(css.footnote = 'text-align: left;',
               css.caption  = 'text-align: left;')
  )
}