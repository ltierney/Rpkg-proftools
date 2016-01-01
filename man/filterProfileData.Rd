\name{filterProfileData}
\alias{filterProfileData}
\title{ Filter Profile Data }
\description{
  Allow profile data to be filtered on several criteria.
}
\usage{
filterProfileData(pd, ..., normalize = FALSE, regex = FALSE)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{...}{filter specifications is \code{filter = value} form as
    described below.}
  \item{normalize}{logical; if true the total hit count is set to the
    total number of hits in the reduced profile data; otherwise the
    original value is retained.}
  \item{regex}{logical; if true the specifications in \code{select}, and
    \code{omit}, and \code{focus} filters are treated as regular
    expressions; otherwise exact matches are required.}
}
\value{A reduced profile data structure.}
\details{
  This function can be used to make plots and summaries more readable or
  relevant by removing functions that are not of direct interest or have
  low hit counts.

  Filters are specified in \code{filter = value} form with \code{value}
  typically specifying a filter level or argument. Possible filters and
  their argument values are:
  \describe{
    \item{\code{select}}{character vector specifying names of functions;
      call stacks not containing functions matching any of these names
      are dropped.}
    \item{\code{omit}}{character vector specifying names of functions; call
      stacks containing functions matching any of these names are dropped.}
    \item{\code{focus}}{character vector specifying names of functions; call
      stacks not containing functions matching any of these names are
      dropped, and functions at the bottom of the stack not matching the
      \code{focus} specification are dropped.}
    \item{\code{skip}}{integer; the number of elements to trim from the
      bottom of the stacks.}
    \item{\code{maxdepth}}{integer; stacks are truncated to have at most
      \code{maxdepth} elements.}
    \item{\code{self.pct}}{numeric; functions at the bottom of the
      stacks with self percentages below this value are removed.}
    \item{\code{total.pct}}{numeric; functions at the top of the stacks with
      total percentages below this value are removed.}
    \item{\code{interval}}{inter vector of length 2 specifying first and last
      sample to use.}
    \item{\code{merge.pct}}{numeric; functions at the top of the stacks
      are removed and stack traces merged until each retained trace accounts
      for at least this percentage of run time.}
  }
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{flatProfile}},
  \code{\link{readProfileData}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{profileCallGraph2Dot}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
plotProfileCallGraph(pd)
plotProfileCallGraph(filterProfileData(pd, self.pct = 1))
plotProfileCallGraph(filterProfileData(pd, self.pct = 1, total.pct = 10))
plotProfileCallGraph(filterProfileData(pd, select = "glm", self.pct=1,
                                       total.pct=10))
}
\keyword{programming}
\keyword{utilities}
