\name{hotPaths}
\alias{hotPaths}
\title{ Hot Paths in Profile Data }
\description{
  Computes and displays hot paths in profiling data.
}
\usage{
hotPaths(pd, value = c("pct", "time", "hits"), self = TRUE,
         srclines = TRUE, GC = FALSE, memory = FALSE, maxdepth = 10,
         self.pct = 0, total.pct = 0, short = ". ", nlines = NA)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{value}{character; show result as percentage, time, or hits.}
  \item{self}{logical; include self time for each stack in the result.}
  \item{srclines}{logical; include source information, if available, or not.}
  \item{GC}{logical; include GC information or not.}
  \item{memory}{logical; include memory use information or not.}
  \item{maxdepth}{integer; stacks are truncated to have at most
                  \code{maxdepth} elements.}
  \item{self.pct}{numeric; stacks with self percent values below this
                  level are dropped.}
  \item{total.pct}{numeric; stacks with total percent values below this
                  level are dropped.}
  \item{short}{character; abbreviation to be used for functions lower on
               the stack.}
  \item{nlines}{integer; number of lines to show. The lines shown are
                the ones with the highest total percentage.}
}
\value{A data frame designed to produce a useful printed result.}
\details{
The hot path ordering sorts stacks in the profile data first by the
frequency with which the bottom functions on the stack are called, with
highest frequency first, then within each bottom function by the
frequency of the bottom two, and so on. Examining the result of
\code{hotPaths} starting with low values of \code{maxdepth} and then
moving to higher levels is a useful way to expore where the
computational effort is concentrated.
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{flatProfile}},
  \code{\link{filterProfileData}},
  \code{\link{readProfileData}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{profileCallGraph2Dot}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
hotPaths(pd)
hotPaths(pd, maxdepth = 8)
}
\keyword{programming}
\keyword{utilities}
