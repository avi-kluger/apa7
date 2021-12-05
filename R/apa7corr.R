# FOR  TESTING
# data("attitude")
# attitude[5, 2:3] <- NA
# data <- attitude
# useLabels = NULL
# listwiseDeletion = FALSE
# disattenuated = TRUE
# reliabilities = c(rep(.9, ncol(data)))
# threeStars = FALSE
# myLabel <- c(
#   "Overall rating (the Dependent Variable)",
#   "Handling of employee complaints",
#   "Does not allow special privileges",
#   "Opportunity to learn",
#   "Raises based on performance",
#   "Too critical",
#   "Advancement"
# )
# alphaVector <- c(.7, .83, 1, 1, .9, .89, .63)
# useLabels = myLabel
# listwiseDeletion = TRUE
# disattenuated = TRUE
# reliabilities = alphaVector


apa7corr <- function (data, useLabels = NULL, listwiseDeletion = FALSE,
          disattenuated = FALSE, reliabilities = c(rep(1, ncol(data))),
          threeStars = FALSE)
{
  N <- x <- NULL

  if (!requireNamespace("sjPlot", quietly = TRUE))
    install.packages("sjPlot")
  if (!requireNamespace("CTT", quietly = TRUE))
    install.packages("CTT")
  if (!requireNamespace("weights", quietly = TRUE))
    install.packages("weights")
  if (!requireNamespace("Hmisc", quietly = TRUE))
    install.packages("Hmisc")
  if (!requireNamespace("psychometric", quietly = TRUE))
    install.packages("psychometric")

  if (disattenuated == TRUE & sum(reliabilities) == ncol(data)) {
    warning("You requested disattenuated correlations but did not supply your own\n vector of reliabilities. Disattenuated correlations are NOT printed.")
    disattenuated <- FALSE
  }
  if (listwiseDeletion == TRUE) {
    if (nrow(data) == nrow(na.omit(data))){
      warning("You asked for listwise deletion, but there are no missing data. Processing countinues with listwiseDeletion == FALSE.")
      listwiseDeletion <-  FALSE
          } else {
      data <- na.omit(data)
   }
  }

  corMatrix <- cor(data, use = "pairwise.complete.obs")
  corMatrix <- CTT::disattenuated.cor(corMatrix, unlist(reliabilities))
  corMatrix <- as.data.frame(corMatrix)

  ciMatrix <- matrix(nrow = nrow(corMatrix), ncol = ncol(corMatrix))

    for (i in 1:ncol(corMatrix)) {
    for (j in 1:ncol(corMatrix)) {
      if(j > i){
      r <- corMatrix[i, j]
      if(r >= 1) {ciMatrix[i, j] <- corMatrix[i, j]
        # print(c(i, j))
        # print(ciMatrix[i, j])
      } else {
        N <- nrow(na.omit(data[, c(i,j)]))
        ciMatrix[i, j] <- psychometric::CIr(r, N, level = 0.95)[2]
        # print(c(i, j))
        # print(ciMatrix[i, j])
      }
    }
  }
}

  divergence.moderate <- max(corMatrix) > 1 | max(ciMatrix, na.rm = TRUE) > .9
  divergence.minor    <- max(corMatrix) < 1 | max(ciMatrix, na.rm = TRUE) > .8

  if (listwiseDeletion == TRUE | sum(is.na(x)) == 0) {
    Col.header <- c("Measure", "Mean", "SD", 1:ncol(data))
  } else {
    Col.header <- c("Measure", "N", "Mean", "SD", 1:ncol(data))
  }
  ft.note  <- "<i>Note</i>. "
  ft.N     <- paste0(" <i>N</i> = ", nrow(data), ".")
  ft.miss  <- " Correlations for variables with missing data are based on pairwise deletion."
  ft.diag  <- " Values in the diagonal are reliabilities."
  ft.above <- " Values above the diagonal are disattenuated correlations."
  ft.bold  <- " Disattenuated correlations for which the upper limit of their confidence interval > .80 are printed in <b>bold</b>. This suggests marginal problem with divergent validity."
  ft.score <- " Disattenuated correlations > 1 and those for which the upper limit of their confidence interval > .90 are printed in <b><u>underscore and bold</u></b>.This suggests moderate to severe problem with divergent validity."
  ft.star  <- "<br>*<i>p</i>< .05. **<i>p</i>< .01."
  ft.3star <- " ***<i>p</i>< .001."


  if (listwiseDeletion == TRUE)          ft.N     <- ""
  if (listwiseDeletion == FALSE)         ft.miss  <- ""
  if (sum(reliabilities) == ncol(data))  ft.diag  <- ""
  if (disattenuated == FALSE)            ft.above <- ""
  if (disattenuated == FALSE)            ft.bold  <- ""
  if (disattenuated == FALSE)            ft.score <- ""
  if (disattenuated == TRUE & divergence.minor == FALSE)    ft.bold  <- ""
  if (disattenuated == TRUE & divergence.moderate == FALSE) ft.score <- ""
  if (ft.3star == FALSE)                 ft.3star <- ""

  ft <- paste0(ft.note, ft.N, ft.miss, ft.diag,
               ft.above, ft.bold, ft.score, ft.star, ft.3star)

  if (is.null(useLabels)) {
    corTable <- data.frame(Measure = names(data))
  } else {
    corTable <- data.frame(Measure = useLabels)
  }
  corTable$Measure <- paste(1:nrow(corTable), corTable$Measure,
                            sep = ". ")
  corTable$N <- sapply(data, function(x) sum(!is.na(x)))
  corTable$Mean <- sapply(data, mean, na.rm = TRUE)
  corTable$SD <- sapply(data, sd, na.rm = TRUE)

  p <- Hmisc::rcorr(as.matrix(data))
  p <- as.matrix(p$P)
  pCorr <- matrix("", nrow = nrow(p), ncol = nrow(p))
  pCorr[p < 0.05] <- "*"
  pCorr[p < 0.01] <- "**"
  if (threeStars == TRUE)
    pCorr[p < 0.001] <- "***"
  pCorr[upper.tri(pCorr)] <- ""

  corMatrix <- sapply(corMatrix, weights::rd, 2)
  if (disattenuated == FALSE)
    corMatrix[upper.tri(corMatrix)] <- ""
  for (i in 1:ncol(corMatrix)) {
    corMatrix[i, i] <- ifelse(corMatrix[i, i] == "1.00",
                              "-", corMatrix[i, i])
  }
  for (i in 1:ncol(corMatrix)) {
    for (j in 1:ncol(corMatrix)) {
        if(corMatrix[i, i] == "-" & corMatrix[j, j] == "-" & j > i){
          corMatrix[i, j] <- ""
      }
    }
  }

  for (i in 1:ncol(corMatrix)) {
    for (j in 1:ncol(corMatrix)) {
      if(corMatrix[i, j] >= 1) {
        corMatrix[i, j] <- paste0("<b><u>", corMatrix[i, j],  "</u></b>")
      } else {
            if(!is.na(ciMatrix[i, j] ) & ciMatrix[i, j]> .90){
            corMatrix[i, j] <- paste0("<b><u>", corMatrix[i, j], "</u></b>")
            } else {
              if(!is.na(ciMatrix[i, j] ) & ciMatrix[i, j]> .80){
                corMatrix[i, j] <- paste0("<b>", corMatrix[i, j],  "</b>")
              }
      }
    }
   }
  }

  corMatrix <- ifelse(corMatrix == ".000000", ".00",
                      corMatrix)
  for (i in 1:ncol(corMatrix)) corMatrix[i, i] <- paste0("(",
                                                         corMatrix[i, i], ")")
  corMatrix[pCorr == "*"] <- paste0(corMatrix[pCorr == "*"], "*")
  corMatrix[pCorr == "**"] <- paste0(corMatrix[pCorr == "**"], "**")
  corMatrix[pCorr == "***"] <- paste0(corMatrix[pCorr ==  "***"], "***")


  corTable <- cbind(corTable, corMatrix)
  row.names(corTable) <- NULL
  colnames(corTable) <- c("Measure", "N", "Mean", "SD", 1:nrow(corTable))
  if (listwiseDeletion == TRUE | sum(is.na(data)) == 0)
    corTable <- corTable[, -c(grep("N", names(corTable)))]
  corTable[, c("Mean", "SD")] <- round(corTable[,c("Mean", "SD")], 2)

  tab <- sjPlot::tab_df(corTable, title = "<b>Table X.</b> <br>\n
              Descriptive Statistics and Correlations for Study Variables",
              col.header = Col.header, show.footnote = TRUE, footnote = ft,
              CSS = list(css.footnote = "text-align: left;",
                         css.caption = "text-align: left;"))
  tab$page.complete <- gsub("double", "1px solid",
                            tab$page.complete)
  tab
}
