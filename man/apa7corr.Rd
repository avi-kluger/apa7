\name{apa7corr}
\alias{apa7corr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an APA 7th correlation table
%%  ~~function to create  ~~
}
\description{Creates a correlation table following the 7th Edition of the Publication Manual of the APA (Table 7.10 on p. 215) and variants including reliabilities in the diagonal and disattenuated correlations above the diagonal.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
apa7corr(data,
         useLabels = NULL,
         listwiseDeletion = FALSE,
         disattenuated = FALSE,
         reliabilities = c(rep(1, ncol(data))),
         threeStars = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric data frame.}
  \item{useLabels}{NULL (default). If labels are supplied (e.g., useLabels = myLabels), the table will present the labels instead of the column names.}
  \item{listwiseDeletion}{FALSE (default). If TRUE, observations with missing data will be deleted, the N columns will be deleted, and the total N will be presented in the Note}
  \item{disattenuated}{FALSE (default). If TRUE, correlations above the diagonal will be disattenuated based on the reliability vector.}
  \item{reliabilities}{The default uses 1.00 for reliabilities and does not print them. You can add your vector of reliabilities (a numerical vector). Supplying your vector will result in showing them in the diagonal. If you supply your own vector and choose disattenuated = TRUE, disattenuated correlations will be printed above the diagonal. If you have some variables with reliabilities and some without (e.g., demographics such as age), use 1.00 for the variables with no reliability. In this case, disattenuated correlations will be printed only for pair of variables for at least one of the reliabilities is < 1.}
  \item{threeStars}{FALSE (default). If TRUE, three starts will be added for p < .001.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
https://apastyle.apa.org/products/publication-manual-7th-edition

Rönkkö, M., & Cho, E. (2020). An updated guideline for assessing discriminant validity. Organizational Research Methods. https://doi.org/10.1177/1094428120968614

%% ~put references to the literature/web site here ~
}
\author{Avraham (Avi) N. Kluger <avraham.kluger@huji.ac.il>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
apa7corr(attitude)

# Prepare labels
myLabel <- c(
  "Overall rating (the Dependent Variable)",
  "Handling of employee complaints",
  "Does not allow special privileges",
  "Opportunity to learn",
  "Raises based on performance",
  "Too critical",
  "Advancement"
)
apa7corr(attitude,
         useLabels = myLabel)

# Make-up a reliability vector (e.g., the reliability for Overall rating = .70).
# Let's assume that two variables have no information about reliability using 1.
alphaVector <- c(.7, .83, 1, 1, .9, .89, .63)
apa7corr(attitude,
         useLabels = myLabel,
         reliabilities = alphaVector)

# Let's assume missing data
attitude[5, 2:3] <- NA
apa7corr(attitude,
         useLabels = myLabel,
         listwiseDeletion = TRUE)

# Ask for disattenuated correlations, BUT note the warning.
apa7corr(attitude,
         useLabels = myLabel,
         listwiseDeletion = FALSE,
         disattenuated = TRUE)

# Fix the above problem and add a reliability vector.
# Note that for pairs of variables without reliability,
# the place for the disattenuated correlation is empty.
apa7corr(attitude,
         useLabels = myLabel,
         listwiseDeletion = FALSE,
         disattenuated = TRUE,
         reliabilities = alphaVector)

# The same as above with listwise deletion
apa7corr(attitude,
         useLabels = myLabel,
         listwiseDeletion = TRUE,
         disattenuated = TRUE,
         reliabilities = alphaVector)

# Add three stars for p < .001
apa7corr(attitude,
         useLabels = myLabel,
         listwiseDeletion = TRUE,
         disattenuated = TRUE,
         reliabilities = alphaVector,
         threeStars = TRUE)

## The function is currently defined as
 function(data,
          useLabels = NULL,
          listwiseDeletion = FALSE,
          disattenuated = FALSE,
          reliabilities = c(rep(1, ncol(data))),
          threeStars = FALSE)
{}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

