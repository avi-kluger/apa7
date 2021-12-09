# FOR  TESTING
# myLabel <- c(
#   "Overall rating (the Dependent Variable)",
#   "Handling of employee complaints",
#   "Does not allow special privileges",
#   "Opportunity to learn",
#   "Raises based on performance",
#   "Too critical",
#   "Advancement"
# )
#
# data[5, 2:3] <- NA
# data$complaints <- -1 * data$complaints
# reliabilities = c(rep(.9, ncol(data)))

# data("attitude")
# data <- attitude
# useLabels = NULL
# listwiseDeletion = FALSE
# disattenuated = TRUE
# threeStars = FALSE
# dagger = TRUE
# alphaVector <- c(.7, .83, 1, 1, .9, .89, .63)
# reliabilities = alphaVector

# useLabels = myLabel

apa7corr <- function (data, useLabels = NULL, listwiseDeletion = FALSE,
          disattenuated = FALSE, reliabilities = c(rep(1, ncol(data))),
          dagger = FALSE, threeStars = FALSE)
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
      warning("You asked for listwise deletion, but there are no missing data. \nProcessing countinues with listwiseDeletion == FALSE.")
      listwiseDeletion <-  FALSE
          } else {
      data <- na.omit(data)
   }
  }

  corMatrix <- cor(data, use = "pairwise.complete.obs")
  ciMatrix  <- matrix(nrow = nrow(corMatrix), ncol = ncol(corMatrix))

  for (i in 1:ncol(corMatrix)) {
    for (j in 1:ncol(corMatrix)) {
      if(j != i){
        r <- abs(corMatrix[i, j])
        N <- nrow(na.omit(data[, c(i,j)]))
        ciMatrix[i, j] <- psychometric::CIr(r, N, level = 0.95)[2]
        # print(c(i, j))
        # print(ciMatrix[i, j])
      }  else {
        ciMatrix[i, j] <- as.numeric(reliabilities) [i]
        }
      }
    }


  corMatrix <- CTT::disattenuated.cor(corMatrix, unlist(reliabilities))
  corMatrix <- as.data.frame(corMatrix)
  corMatrixNumeric <- as.data.frame(corMatrix)
  corMatrixNumeric <- as.data.frame(corMatrixNumeric)

  ciMatrix <- CTT::disattenuated.cor(ciMatrix, unlist(reliabilities))

  divergence.severe   <- ciMatrix[which.max(abs(ciMatrix))] > 1
  divergence.moderate <- ciMatrix[which.max(abs(ciMatrix))] > .9
  divergence.minor    <- ciMatrix[which.max(abs(ciMatrix))] > .8

  if (listwiseDeletion == TRUE | sum(is.na(data)) == 0) {
    Col.header <- c("Measure", "Mean", "SD", 1:ncol(data))
  } else {
    Col.header <- c("Measure", "N", "Mean", "SD", 1:ncol(data))
  }

  if (is.null(useLabels)) {
    corTable <- data.frame(Measure = names(data))
  } else {
    corTable <- data.frame(Measure = useLabels)
  }
  corTable$Measure <- paste(1:nrow(corTable), corTable$Measure, sep = ". ")
  corTable$N       <- format(sapply(data, function(x) sum(!is.na(x))), big.mark = ",")
  corTable$Mean    <- sapply(data, mean, na.rm = TRUE)
  corTable$SD      <- sapply(data, sd, na.rm = TRUE)

  p <- Hmisc::rcorr(as.matrix(data))
  p <- as.matrix(p$P)
  pCorr <- matrix("", nrow = nrow(p), ncol = nrow(p))
  pCorr[p < 0.1]  <- "\u2020"
  pCorr[p < 0.05] <- "*"
  pCorr[p < 0.01] <- "**"
  if (threeStars == TRUE) pCorr[p < 0.001] <- "***"
  pCorr[upper.tri(pCorr)] <- ""

  corMatrix <- sapply(corMatrixNumeric, weights::rd, 2)
  if (disattenuated == FALSE)     corMatrix[upper.tri(corMatrix)] <- ""
  for (i in 1:ncol(corMatrix)) {
    corMatrix[i, i] <-
      ifelse(corMatrix[i, i] == "1.00", "-", corMatrix[i, i])
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
      if (!is.na(ciMatrix[i, j] ) & j > i) {
      if(ciMatrix[i, j] > 1) {
        corMatrix[i, j] <- paste0("<i><b><u>", corMatrix[i, j], "</i></u></b>")
      } else {
        if(ciMatrix[i, j] > .90){
          corMatrix[i, j] <- paste0("<b><u>", corMatrix[i, j], "</u></b>")
        } else {
          if(ciMatrix[i, j] > .80){
            corMatrix[i, j] <- paste0("<b>", corMatrix[i, j], "</b>")
               }
             }
          }
       }
    }
  }


    corMatrix <- ifelse(corMatrix == ".000000", ".00", corMatrix)
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

  ft.note  <- "<i>Note</i>. "
  ft.N     <- paste0(" <i>N</i> = ", nrow(data), ".")
  ft.miss  <- " Correlations for variables with missing data are based on pairwise deletion."
  ft.diag  <- " Values in the diagonal are reliabilities."
  ft.above <- " Values above the diagonal are disattenuated correlations."
  ft.empty <- " For pairs of variables for which reliability is unavailable for both, the cell for the disattenuated correlations is left empty."
  ft.dis1  <- " Disattenuated correlations for which the upper limit of their confidence interval > "
  ft.dis2  <- " problem with divergent validity."
  ft.bold  <- paste0(ft.dis1, ".80 are printed in <b>bold</b> and suggests a marginal ", ft.dis2)
  ft.score <- paste0(ft.dis1, ".90 are printed in <b><u>underscore and bold</u></b> and suggests a moderate ", ft.dis2)
  ft.itali <- paste0(ft.dis1, "1 are printed in <i><b><u>italics, underscore and bold</u></b></i> and suggests a severe ", ft.dis2)
  ft.sig   <- "<br>"
  ft.dag   <- "†<i>p</i> < .10."
  ft.star  <- "*<i>p</i> < .05. **<i>p</i> < .01."
  ft.3star <- " ***<i>p</i> < .001."


  if (listwiseDeletion == FALSE & sum(is.na(data)) > 0) ft.N    <- ""
  if (listwiseDeletion == FALSE)         ft.miss  <- ""
  if (sum(reliabilities) == ncol(data))  ft.diag  <- ""
  if (disattenuated == FALSE)            ft.above <- ""
  ft.empty <- ifelse(sum(corMatrix == "(-)") > 1 & disattenuated == TRUE,
                     ft.empty, "")
  if (disattenuated == FALSE)            ft.bold  <- ""
  if (disattenuated == FALSE)            ft.score <- ""
  if (disattenuated == FALSE |
      (disattenuated == TRUE & divergence.minor == FALSE))    ft.bold  <- ""
  if (disattenuated == FALSE |
      (disattenuated == TRUE & divergence.moderate == FALSE)) ft.score <- ""
  if (disattenuated == FALSE |
      (disattenuated == TRUE & divergence.severe == FALSE))   ft.itali <- ""
  if (dagger == FALSE)                   ft.dag   <- ""
  if (threeStars == FALSE)               ft.3star <- ""

  ft <- paste0(ft.note, ft.N, ft.miss, ft.diag,
               ft.above, ft.empty, ft.bold, ft.score, ft.itali,
               ft.sig, ft.dag, ft.star, ft.3star)

  tab <- sjPlot::tab_df(corTable, title = "<b>Table X.</b> <br>\n
              Descriptive Statistics and Correlations for Study Variables",
              col.header = Col.header, show.footnote = TRUE, footnote = ft,
              CSS = list(css.footnote = "text-align: left;",
                         css.caption = "text-align: left;"),
              encoding = "Windows-1252")
  tab$page.complete <- gsub("double", "1px solid",
                            tab$page.complete)
  tab
}
