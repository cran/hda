\name{showloadings}
\alias{showloadings}
\title{Loadings plot for heteroscedastic discriminant analysis}
\description{
Visualizes the loadings of the original variables on the components of the transformed 
discriminant space of reduced dimension.
}
\usage{
showloadings(object, comps = 1:object$reduced.dimension)
}
\arguments{
  \item{object}{An object of class \code{hda}.}
  \item{comps}{A vector of component ids for which the loadings should be displayed.}
}
\details{
Scatterplots where the loading of any variable on any hda component is represented by one point of different 
colour and shape. Gives an idea of what variables do mainly contribute to the different discriminant components.      
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
\seealso{\code{\link{hda}}, \code{\link{predict.hda}}, \code{\link{plot.hda}}}

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
plot(as.data.frame(rotatedspace), col = classes)

# apply heteroscedastic discriminant analysis and plot data in discriminant space
hda.res <- hda(rotatedspace, classes)

# visualize loadings
showloadings(hda.res)
}

\keyword{classif}
\keyword{multivariate}
