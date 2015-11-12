\name{rgl-package}
\title{Tools for examining R profile output}
\alias{proftools-package}
\alias{proftools}
\description{
  Tools for examining and displaying ouptut from the \code{Rprof} R
  profiling tool.
}
\details{
  \pkg{proftools} provides a set of tools for summarizing and displaying
  time profile outpus produced by R's \code{\link{Rprof}}.

  The starting point for a profiling analysis using \pkg{proftools} is
  to profile code using \code{\link{Rprof}} and then use
  \code{\link{readProfileData}} to read in the profile data into a
  sutable format for furhter processing. An alternative is to use the
  \code{\link{profileExpr}} function to handle profiling and reading in
  one step. The function \code{\link{filterProfileData}} can be used to
  narrow the profile data to particular regions of interest.

  The summary functions \code{\link{funSummary}} and
  \code{\link{callSummary}} produce summaries at the function and call
  level. \code{\link{pathSummary}} produces a summary for each unique
  call stack, or path; and \code{\link{hotPaths}} identifies produces
  path data ordered to show the hottest paths first.

  If source information is recorded when profiling then
  \code{\link{srcSummary}} to show profiling by source lines, and
  \code{\link{annotateSource}} produces an annotated version of the
  source files.

  The \code{\link{plot}} method for profile data objects can produce
  call graphs, tree maps, flame graphs, and time graphs; the
  \code{type} argument choses the particular visualization to
  produce. These graphs can also be produced by the functions
  \code{\link{plotProfileCallGraph}}, \code{\link{calleeTreeMap}}, and
  \code{\link{flameGraph}}.

  The function \code{\link{writeCallgrindFile}} writes a file for use
  by the code{kcachegrind} program available on some operating
  systems.

  \code{\link{flatProfile}}

}
%%\seealso{****}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
funSummary(pd)
callSummary(pd)
pathSummary(pd)
hotPaths(pd)
plot(pd)
plot(filterProfileData(pd, focus = "glm", self.pct=1, total.pct=10))
}
\keyword{programming}
\keyword{utilities}
