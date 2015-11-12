\name{plot.proftools_profData}
\alias{plot.proftools_profData}
\title{Plot Profile Data}
\description{
  Plot method for objects of class \code{"proftools_profData"}.
}
\usage{
\method{plot}{proftools_profData}(x, type = c("call", "tree", "flame", "time"),
                                  \dots)
}
\arguments{
  \item{x}{an object of class \code{"proftools_profData"}.}
  \item{type}{the type of plot to be drawn; default is a call graph.}
  \item{\dots}{additional arguments for specific plot types.}
}
\details{
  Depending on the \code{type} argument the profile data \code{plot}
  method calls \code{link{plotProfileCallGraph}},
  \code{\link{calleeTreeMap}}, or \code{\link{flameGraph}} with
  \code{order = "hot"} or \code{orter = "time"}.
}
\seealso{\code{\link{plotProfileCallGraph}},
  \code{\link{calleeTreeMap}}, \code{\link{flameGraph}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
plot(pd)
plot(pd, style = plain.style)
plot(pd, type = "call")
plot(pd, type = "flame")
}
\keyword{hplot}
\keyword{ts}
