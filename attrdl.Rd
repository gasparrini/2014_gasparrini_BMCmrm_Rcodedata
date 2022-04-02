\name{attrdl}
\alias{attrdl}

\title{ Attributable Risk from a DLNM }

\description{
This function computes estimates of attributable risk (fractions and numbers) from a distributed lag non-linear model (DLNM).
}

\usage{
attrdl(x, basis, cases, model=NULL, coef=NULL, vcov=NULL, model.link=NULL,
  type="af", dir="back", tot=TRUE, cen, range=NULL, sim=FALSE, nsim=5000)
}

\arguments{
  \item{x }{ exposure vector or matrix of exposure histories for which the attributable risk needs to be computed.}
  \item{basis }{ an object of class \code{"crossbasis"} used for fitting the model.}
  \item{cases }{ the vector of cases or matrix of future cases corresponding to \code{x}.}
  \item{model }{ the fitted model, to be defined with a log link function.}
  \item{coef, vcov, model.link }{ user-provided coefficients, (co)variance matrix and model link for the prediction.}
  \item{type }{ measure of attributable risk, either \code{"af"} (fraction, default) or \code{"an"} (number).}
  \item{dir }{ direction used for computing the attributable risk, either \code{"back"} (backward, default) or \code{"forw"} (forward).}
  \item{tot }{ logical. If \code{TRUE} (default), the total estimate is returned.}
  \item{cen }{ numeric scalar. The reference value used as counterfactual scenario.}
  \item{range }{ the range of exposure. If \code{NULL}, the whole range is used.}
  \item{sim }{ logical. If \code{TRUE}, simulation samples are returned.}
  \item{nsim }{ number of simulation samples (only for \code{nsim=TRUE}).}
}

\details{
This function computes the attributable fraction or number for a specific exposure scenario and associated cases, given an estimated exposure-lag-response association defined by a DLNM. Either \emph{forward} or \emph{backward} versions of attributable risk measures are available in this setting. The method is described by Gasparrini and Leone (2014), see references below. The function works in combination with other functions in the package \pkg{dlnm}, which is assumed to be available.

The exposure and cases are provided by the arguments \code{x} and \code{cases}, respectively. The original cross-basis and fitted model containg it used for estimation are provided by the arguments \code{basis} and \code{model}, respectively. Alternatively, the user can provide estimated coefficients and (co)variance matrix with \code{coef} and \code{vcov}.

The function works both with time series and non-time series data. In a time series setting, both \code{x} and \code{cases} represent a complete series of ordered observations. More generally, the user can apply this function for any kind of data: in this case \code{x} must be a matrix of lagged exposures when \code{dir="back"}, and \code{cases} must be a matrix of future cases when \code{dir="forw"}. The function can compute the total attributable risk (if \code{tot=TRUE}, the default) or the contribution for each observation. The argument \code{cen} defines the value used as counterfactual scenario.

If \code{sim=TRUE}, the function computes samples of the attributable risk measures by simulating from the assumed normal distribution of the estimated coefficients (only implemented for total estimates). These samples can be used to defined empirical confidence intervals.
}

\value{
By default, a numeric scalar representing the total attributable fraction or number. If \code{sim=TRUE}, a vector of the simulated samples with length \code{nsim}. If \code{tot=FALSE}, a vector with contributions for all the observations (see Note below). These quantities are defined versus a counterfactual scenario defined through the argument \code{cen}.
}

\references{
Gasparrini A, Leone M. Attributable risk from distributed lag models. \emph{BMC Medical Research Methodology}. 2014; [freely available \href{http://www.ag-myresearch.com/2014_gasparrini_bmcmrm.html}{here}]
}

\author{Antonio Gasparrini, \email{antonio.gasparrini@lshtm.ac.uk}}

\note{
The function handles missing values in both the \code{x} and \code{cases} objects, excluding incomplete observations (also due to lagging) accordingly. However, the total attributable number is rescaled to match the fraction using as denominator the total observed number in \code{cases}. This approach uses the all the available information even in the presence of missing values in \code{x}. All of this under the assumption that the missing mechanism is unrelated with both exposure and cases values.

The functions can be also used with estimates from DLNMs reduced to the overall cumulative exposure-response through the function \code{crossreduce} in the package \pkg{dlnm}. In this case, the modified coefficients and (co)variance matrix of the reduced cross-basis in \code{basis} must be passed using the arguments \code{coef} and \code{vcov}, together with the information on the link function in \code{model.link}. This option can be useful when the original estimates from the full cross-basis are not available any more, for example following a meta-analysis. Given the lag-specific estimates are not available in this case, only the forward version of attributable risk (\code{dir="forw"}) can be computed. See Gasparrini and Leone (2014) for further info.
}

\section{Warnings}{
The code composing this function has not been systematically tested. The presence of bugs cannot be ruled out. Also, although written generically for working in different scenarios and data, the function has not been tested in contexts different than the example included in the article from Gasparrini and Leone (2014). It is responsibility of the user to check the reliability of the results in different applications.
}

\seealso{
Other functions in the package \pkg{dlnm}. See \href{https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata}{here} for an updated version of this function.
}

\examples{
# load the package dlnm and the function attrdl
library(dlnm)
source("attrdl.R")

# define the cross-basis and fit the model
cb <- crossbasis(chicagoNMMAPS$temp, lag=30, argvar=list(fun="bs",
  knots=c(-10,3,18)), arglag=list(knots=c(1,3,10)))
library(splines)
model <- glm(death ~  cb + ns(time, 7*14) + dow, 
  family=quasipoisson(), chicagoNMMAPS)

# global backward attributable risk of temperature (number and fraction)
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,type="an",cen=21)
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,cen=21)

# global forward attributable fraction
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,dir="forw",cen=21)

# empirical confidence intervals
afsim <- attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,cen=21,
  sim=TRUE,nsim=1000)
quantile(afsim,c(2.5,97.5)/100)

# attributable fraction component due to heat and cold
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,cen=21,range=c(21,100))
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,cen=21,range=c(-100,21))

# daily attributable deaths in the second month
attrdl(chicagoNMMAPS$temp,cb,chicagoNMMAPS$death,model,type="an",
  tot=FALSE,cen=21)[31:60]
}

