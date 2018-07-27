\name{readProfileData}
\alias{readProfileData}
\alias{readRStudioProfileCacheData}
\title{ Read Rprof Profile Data }
\description{
  Reads in Rprof profile data for further processing.
}
\usage{
readProfileData(filename = "Rprof.out")
readRStudioProfileCacheData()
}
\arguments{
  \item{filename}{Name of a file produced by \code{Rprof()} }
}
\value{
  R representation of Rprof data,
}
\details{
   \code{readProfileData} reads the data in the file produced by
   \code{Rprof} into a data structure for processing by other functions.
   The details of the structure are subject to change.

   \code{readRStudioProfileCacheData} returns the data for the most
   recent profiling run in the RStudio profile cache, or \code{NULL}
   if no data is available.
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{flatProfile}},
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
