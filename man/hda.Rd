\name{hda}
\alias{hda}
\alias{hda.default}
\alias{hda.formula}
\alias{print.hda}
\title{Heteroscedastic discriminant analysis}
\description{Computes a linear transformation loadings matrix for discrimination of classes with unequal covariance matrices.}
\usage{
hda(x, ...)
\method{hda}{default}(x, grouping, newdim = 1:(ncol(x)-1), crule = FALSE, 
             reg.lamb = NULL, reg.gamm = NULL, initial.loadings = NULL, 
             sig.levs = c(0.05,0.05), noutit = 7, ninit = 10, verbose = TRUE, ...)
\method{hda}{formula}(formula, data = NULL, ...)
}
\arguments{
  \item{x}{A matrix or data frame containing the explanatory variables. The method is restricted to numerical data.}
  \item{grouping}{A factor specifying the class for each observation.}
  \item{formula}{A formula of the form \code{grouping ~ x1 + x2 + ...} That is, the response is the grouping factor and the right hand side specifies the (non-factor) discriminators.}
  \item{data}{Data frame from which variables specified in formula are to be taken.}
  \item{newdim}{Dimension of the discriminative subspace. The class distributions are assumed to be equal in the remaining dimensions. 
                Alternatively, a vector of integers can be specified which is then computed until for the first time both tests on equal means as well 
                as homoscedasticity do not reject. This option is to be be applied with care and the resulting dimension should be 
                checked manually.}  
  \item{crule}{Logical specifying whether a \code{\link[e1071]{naiveBayes}} classification rule should be computed. Requires package \code{e1071}.}
  \item{reg.lamb}{Parameter in [0,1] for regularization towards equal covariance matrix estimations of the classes (in the original space): 
                  0 means equal covariances, 1 (default) means complete heteroscedasticity.}
  \item{reg.gamm}{Similar to \code{reg.lambd}: parameter for shrinkage towards diagonal covariance matrices of equal variance in all variables where 0 means diagonality. Default is no shrinkage.}
  \item{initial.loadings}{Initial guess of the matrix of loadings. Must be quadratic of size \code{ncol(x)} Default is the identity matrix. By specification of \code{initial.loadings = "random"} a random orthonormal matrix will be generated using \code{qr.Q(qr())} of a random matrix with uniformly distributed elements.}
  \item{sig.levs}{Vector of significance levels for eqmean.test (position 1) and homog.test (pos. 2) to stop search for an appropriate dimension of the reduced space.}
  \item{noutit}{Number iterations of the outer loop, i.e. iterations of the likelihood. Default is 7.}
  \item{ninit}{Number of iterations of the inner loop, i.e. reiterations of the loadings matrix within one iteration step of the likelihood.}
  \item{verbose}{Logical indicating whether iteration process should be displayed.}
  \item{\dots}{For \code{hda.formula}: Further arguments passed to function \code{hda.default} such as \code{newdim}. For \code{hda.default}: currently not used.}
}
\details{
The function returns the transformation that maximizes the likelihood if the classes are normally distributed 
but differ only in a \code{newdim} dimensional subspace and have equal distributions in the remaining dimensions  
(see Kumar and Andreou, 1998). The scores are uncorrelated for all classes. The algorithm is implemented as it is proposed by 
Burget (2006). Regularization is computed as proposed by Friedman et al. (1989) and Szepannek et al. (2009).

}
\value{
  Returns an object of class \emph{hda}. 
  \item{hda.loadings}{Transformation matrix to be post-multiplied to new data.}
  \item{hda.scores}{Input data after hda transformation. Reduced discriminative space are the first \code{newdim} dimensions.}
  \item{grouping}{Corresponding class labels for \code{hda.scores} data. Identical to input grouping.}
  \item{class.dist}{Estimated class means and covariance matrices in the transformed space.}
  \item{reduced.dimension}{Input parameter: dimension of the reduced space.}
  \item{naivebayes}{Object of class \code{naiveBayes} trained on input data in the reduced space for classification 
        of new (transformed) data. Its computation must be specified by input the parameter \code{crule}.}
  \item{comp.acc}{Matrix of accuracies per component and class: reports up to which degree each class k can be classified (\eqn{P(f_k>f_{\ne k})}{P(f[k]>f[!k])}) correctly according to the estimated (normal) distribution in any single component in the identified subspace. Meaningful for reasons of interpretability as HDA is invariant to reordering of the components.}        
  \item{vlift}{Returns the variable importance in terms of ratio between the accuracy \code{comp.acc} and the resulting accuracy that results if single variable loadings are set to 0. The first element describes overall accuracy lift where  the second element is an array of dimension (number of classes, number of components in reduced space, number of variables) specifying the lifts for recognition each class separately.}
  \item{reg.lambd}{Input regularization parameter.}
  \item{reg.gamm}{Input regularization parameter.}
  \item{eqmean.test}{Test on equal means of the classes in the remaining dimensions like in \code{manova} based on Wilk's lambda.}
  \item{homog.test}{Test on homoscedasticity of the classes in the remaining dimensions (see e.g. Fahrmeir et al., 1984, p.75.)}
  \item{hda.call}{(Matched) function call.}
  \item{initial.loadings}{Initialization of the loadings matrix.}
  \item{trace.dimensions}{Matrix of p values for different subspace dimensions (as specified in \code{newdim}).}
}
\references{
Burget, L. (2006): \emph{Combination of speech features using smoothed heteroscedastic discriminant analysis.} 
Proceedings of Interspeech 2004, pp. 2549-2552.

Fahrmeir, L. and Hamerle, A. (1984): \emph{Multivariate statistische Verfahren.} de Gruyter, Berlin.

Friedman, J. (1989): \emph{Regularized discriminant analysis.} JASA 84, 165-175.

Kumar, N. and Andreou, A. (1998): \emph{Heteroscedastic discriminant analysis and reduced rank HMMs 
for improved speech recognition.} Speech Communication 25, pp.283-297.

Szepannek G., Harczos, T., Klefenz, F. and Weihs, C. (2009): \emph{Extending features for automatic speech recognition 
by means of auditory modelling.} In: Proceedings of European Signal Processing Conference (EUSIPCO) 2009, Glasgow, pp.1235-1239. 
}

\author{Gero Szepannek}
\seealso{\code{\link{predict.hda}}, \code{\link{showloadings}}, \code{\link{plot.hda}}}

\examples{
library(mvtnorm)
library(MASS)

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
x       <- rbind(xa, xb)
classes <- as.factor(c(rep(1,n), rep(2,n)))
# rotate simulated data
symmat <- matrix(runif(5^2),5)
symmat <- symmat + t(symmat)
even   <- eigen(symmat)$vectors
rotatedspace <- x \%*\% even
plot(as.data.frame(rotatedspace), col = classes)

# apply linear discriminant analysis and plot data on (single) discriminant axis
lda.res <- lda(rotatedspace, classes)
plot(rotatedspace \%*\% lda.res$scaling, col = classes, 
     ylab = "discriminant axis", xlab = "Observation index")

# apply heteroscedastic discriminant analysis and plot data in discriminant space
hda.res <- hda(rotatedspace, classes)
plot(hda.res$hda.scores, col = classes)

# compare with principal component analysis
pca.res  <- prcomp(as.data.frame(rotatedspace), retx = TRUE)
plot(as.data.frame(pca.res$x), col=classes)

# Automatically build classification rule
# this requires package e1071
hda.res2 <- hda(rotatedspace, classes, crule = TRUE)
}
\keyword{classif}
\keyword{multivariate}
