\name{profileExpr}
\alias{profileExpr}
\title{ Read Rprof Profile Data }
\description{
  Return profile data collected while evaluating an expression.
}
\usage{
profileExpr(expr, GC = TRUE, srclines = TRUE, memory = TRUE)
}
\arguments{
  \item{expr}{Valid \R expression to be profiled.}
  \item{srclines}{logical; include source information, if available, or not.}
  \item{GC}{logical; include GC information or not.}
  \item{memory}{logical; include memory use information or not.}
}
\value{
  R representation of Rprof data.
}
\details{
   \code{profileExpr} uses \code{\link{Rprof}} to profile execution of
   \code{expr} and returns the profile data read in using
   \code{\link{readProfileData}}.  By default GC and source information
   are included in the profile data. If memory profiling is supported,
   memory use information is also included by default.
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{funSummary}},
  \code{\link{srcSummary}},
  \code{\link{plotProfileCallGraph}},
}
\examples{
## This defines the function rw()
source(system.file("samples", "rw.R", package="proftools"))

## Execute the function and collect profile data
pd <- profileExpr(rw(200000))

## Examine the profiel data
funSummary(pd)
callSummary(pd)
hotPaths(pd)
srcSummary(pd)
plotProfileCallGraph(pd)
}
\keyword{programming}
\keyword{utilities}
