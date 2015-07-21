\name{rgl-package}
\title{Tools for examining R profile output}
\alias{proftools-package}
\alias{proftools}
\description{
  Tools for examining and displaying ouptut from the \code{Rprof} R
  profiling tool.
}
\details{
  \pkg{proftools} is

  \code{profileExpr}
  \code{readProfileData}
  \code{filterProfileData}
  
  Summaries:
  \code{callSummary} \code{funSummary} \code{pathSummary} \code{hotPaths}
  \code{flatProfile}

  Source-level summaries:
  \code{srcSummary} \code{annotateSource}
  
  Visualizations: **** need plot method?
  \code{plotProfileCallGraph} \code{google.style} \code{plain.style}
  \code{calleeTreeMap} \code{flameGraph}
  \code{printProfileCallGraph} \code{profileCallGraph2Dot}  
  
  Cachegrind support:
  \code{writeCallgrindFile}   

  %%See the first example below to display the ChangeLog.
}
%%\seealso{****}
%%\examples{
%%file.show(system.file("NEWS", package = "proftools"))
%%example(hotPaths)
%%}
\keyword{programming}
\keyword{utilities}
