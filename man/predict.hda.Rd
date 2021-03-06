\name{predict.hda}
\alias{predict.hda}
\title{Heteroscedastic discriminant analysis}
\description{Computes linear transformation of new data into lower dimensional discriminative space using some model produced by \code{\link{hda}}.}
\usage{
\method{predict}{hda}(object, newdata, alldims = FALSE, task = c("dr", "c"), ...)
}
\arguments{
  \item{object}{Model resulting from a call of \code{hda}.}
  \item{newdata}{A matrix or data frame to be transformed into lower dimensional space of the same dimension as the data used for building the model.}
  \item{alldims}{Logical flag specifying whether the result should contain only the reduced space (default) or should also 
                include the redundant dimensions and thus be of the same dimension as the input data. In this case the reduced
                space is given by the first \code{newdim} columns.}
  \item{task}{\code{"dr"} for standard application of the \code{hda} model to \code{newdata}. Choose \code{"c"} for classification of new data. 
              This is an interface to predict function of \code{\link[e1071]{naiveBayes}}. 
              The option can be chosen if \code{crule = TRUE} has been specified in the \code{hda()} call.}
  \item{\dots}{Further arguments to be passed to the \code{\link[e1071]{naiveBayes}} predict function.}        
}

\value{
  If option \code{type = "dr"} the transformed data are returned. For \code{type = "c"} both the transformed data as well as
  the resulting object of the naive Bayes prediction are returned.
}
\references{
Kumar, N. and Andreou, A. (1998): \emph{Heteroscedastic discriminant analysis and reduced rank HMMs 
for improved speech recognition.} Speech Communication 25, pp. 283-297.

Szepannek G., Harczos, T., Klefenz, F. and Weihs, C. (2009): \emph{Extending features for automatic speech recognition 
by means of auditory modelling.} In: Proceedings of European Signal Processing Conference (EUSIPCO) 2009, Glasgow, pp. 1235-1239. 
}

\author{Gero Szepannek}
\seealso{\code{\link{hda}}, \code{\link{showloadings}}, \code{\link{plot.hda}}}

\examples{
library(mvtnorm)
library(MASS)

# simulate data for two classes
n           <- 50
meana       <- meanb <- c(0,0,0,0,0)
cova        <- diag(5)
cova[1,1]   <- 0.2
for(i in 3:4){
  for(j in (i+1):5){cova[i,j] <- cova[j,i] <- 0.75^(j-i)}
  }
covb       <- cova
diag(covb)[1:2]  <- c(1,0.2)

xa      <- rmvnorm(n,meana,cova)
xb      <- rmvnorm(n,meanb,covb)
x       <- rbind(xa,xb)
classes <- as.factor(c(rep(1,n),rep(2,n)))
# rotate simulated data
symmat <- matrix(runif(5^2),5)
symmat <- symmat + t(symmat)
even   <- eigen(symmat)$vectors
rotatedspace <- x \%*\% even

# apply heteroscedastic discriminant analysis and plot data in discriminant space
hda.res <- hda(rotatedspace, classes)

# simulate new data
xanew      <- rmvnorm(n,meana,cova)
xbnew      <- rmvnorm(n,meanb,covb)
xnew       <- rbind(xanew,xbnew)
classes <- as.factor(c(rep(1,n),rep(2,n)))
newrotateddata <- x \%*\% even
plot(as.data.frame(newrotateddata), col = classes)

# transform new data 
prediction <- predict(hda.res, newrotateddata)
plot(as.data.frame(prediction), col = classes)

# predict classes for new data on automatically computed naive Bayes classification rule 
# this requires package e1071
hda.res2 <- hda(rotatedspace, classes, crule = TRUE)
prediction2 <- predict(hda.res2, newrotateddata, task = "c")
prediction2
}
\keyword{classif}
\keyword{multivariate}
