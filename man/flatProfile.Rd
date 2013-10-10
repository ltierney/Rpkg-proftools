\name{flatProfile}
\alias{flatProfile}
\title{Flat Profile for Rprof Profile Data}
\description{
Computes a flat profile reflecting time spent in functions themselves
(self) and functions plus callees (total).
}
\usage{
flatProfile(pd, byTotal = TRUE, GC = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{byTotal}{logical; sort by total time if true, self time if not. }
  \item{GC}{logical; include GC information or not.}
}
\value{
  A matrix with one row per function recorded in the profile data.
}
\details{
  If \code{byTotal} is true then the result is analogous to the
  \code{by.total} component of the result returned by \code{summaryRprof}.
  Otherwise, the result is analogous to the \code{by.self} component
  returned by \code{summaryRprof} but with an additional cumulative 
  self times column.  The result returned when \code{byTotal} is not true
  is analogous to the flat profile produced by \code{gprof}.
}
\references{ User manual for \code{gprof}, the GNU profiler. }
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{readProfileData}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{printProfileCallGraph}},
  \code{\link{profileCallGraph2Dot}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
flatProfile(pd)
flatProfile(pd, FALSE)
}
\keyword{programming}
\keyword{utilities}
