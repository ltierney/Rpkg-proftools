\name{filterProfileData}
\alias{filterProfileData}
\title{ Filter Profile Data }
\description{
  Allow profile data to be filtered on several criteria.
}
\usage{
filterProfileData(pd, select, omit, maxdepth = NA,
                  self.pct = 0, total.pct = 0, focus = FALSE,
                  normalize = FALSE, regex = FALSE)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{select}{character vector specifying names of functions; call
                stacks not containing functions matching any of these
                names are dropped.}
  \item{omit}{character vector specifying names of functions; call
              stacks containing functions matching any of these names
              are dropped.}
  \item{maxdepth}{integer; stacks are truncated to have at most
                  \code{maxdepth} elements.}
  \item{self.pct}{logical; functions at the bottom of the stacks with
                  self percentages below this value are removed.}
  \item{total.pct}{logical; functions at the top of the stacks with
                   total percentages below this value are removed.}
  \item{focus}{logical; for use with the \code{select} argument. If
               true, functions at the bottom of the stack not matching the
               \code{select} specification are dropped.}
  \item{normalize}{logical; if true the total hit count is set to the
                   total number of hits in the reduced profile data;
                   otherwise the original value is retained.}
  \item{regex}{logical; if true the specifications in \code{select} and
               \code{omit} are trated as regular expressions; otherwise
               exact matches are required.}
}
\value{A reduced profile data structure.}
\details{
  This function can be used to make plots an summaries more readable or
  relavant by removing functions that are not of direct interest or have
  low hit counts.
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
plotProfileCallGraph(pd, style=google.style, score = "total")
plotProfileCallGraph(filterProfileData(pd, self.pct=1),
                     style=google.style, score = "total")
plotProfileCallGraph(filterProfileData(pd, self.pct=1, total.pct=10),
                     style=google.style, score = "total")
plotProfileCallGraph(filterProfileData(pd, select = "glm", self.pct=1,
                                       total.pct=10),
                      style=google.style, score = "total")
}
\keyword{programming}
\keyword{utilities}
