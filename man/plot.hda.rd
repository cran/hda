\name{plot.hda}
\alias{plot.hda}
\title{Plot transformed data}
\description{
Visualizes the scores on selected components of the  
discriminant space of reduced dimension.
}
\usage{
\method{plot}{hda}(x, comps = 1:x$reduced.dimension, scores = TRUE, col = x$grouping, ...)
}
\arguments{
  \item{x}{An object of class \code{hda}.}
  \item{comps}{A vector of component ids for which the data should be displayed.}
  \item{scores}{Logical indicating whether the scores in the projected space should be plotted. 
                If FALSE estimated densities are plotted.}
  \item{col}{Color vector for the data to be displayed. Per default, different colors represent the classes.}
  \item{\dots}{Further arguments to be passed to the plot function.}
}
\details{
Scatterplots of the scores or estimated densities.        
}
\value{No value is returned.
}
\references{
Kumar, N. and Andreou, A. (1998): \emph{Heteroscedastic discriminant analysis and reduced rank HMMs 
for improved speech recognition.} Speech Communication 25, pp.283-297.

Szepannek G., Harczos, T., Klefenz, F. and Weihs, C. (2009): \emph{Extending features for automatic speech recognition 
by means of auditory modelling.} In: Proceedings of European Signal Processing Conference (EUSIPCO) 2009, Glasgow, pp.1235-1239. 
}

\author{Gero Szepannek}
\seealso{\code{\link{hda}}, \code{\link{predict.hda}}, \code{\link{showloadings}}}

\examples{
library("mvtnorm")
library("MASS")

# simulate data for two classes
n           <- 50
meana       <- meanb <- c(0,0,0,0,0)
cova        <- diag(5)
cova[1,1]   <- 0.2
for(i in 3:4){
  for(j in (i+1):5){
    cova[i,j] <- cova[j,i] <- 0.75^(j-i)}
  }
covb       <- cova
diag(covb)[1:2]  <- c(1,0.2)

xa      <- rmvnorm(n, meana, cova)
xb      <- rmvnorm(n, meanb, covb)
x       <- rbind(xa,xb)
classes <- as.factor(c(rep(1,n), rep(2,n)))
## rotate simulated data
symmat <- matrix(runif(5^2),5)
symmat <- symmat + t(symmat)
even   <- eigen(symmat)$vectors
rotatedspace <- x \%*\% even
plot(as.data.frame(rotatedspace), col = classes)

# apply heteroscedastic discriminant analysis and plot data in discriminant space
hda.res <- hda(rotatedspace, classes)

# plot scores
plot(hda.res)
}

\keyword{classif}
\keyword{multivariate}
