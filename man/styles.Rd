\name{styles}
\alias{plain.style}
\alias{google.style}
\title{Style Specifications for Call Graphs}
\description{
  Styles providing coordinated settings of display parameters for the
  call graph display functions \code{plotProfileCallGraph} and
  \code{profileCallGraph2Dot}.
}
\usage{
plain.style
google.style
}
\value{
  A list containing the following components:

    \item{layout}{The layout method to use: One of \code{"dot"},
                \code{"neato"}, and \code{"twopi"}.}
  \item{score}{ character string specifying whether to use total time
                or self time for coloring nodes/edges; no color used if
                missing. }
  \item{transfer}{ function; maps score values in unit interval to unit
                   interval }
  \item{nodeColorMap,edgeColorMap}{ character vectors of color
    specifications as produced by \code{rainbow}; \code{transfer} of
    score is mapped to color }
  \item{mergeCycles}{ logical; whether to merge each cycle of recursion into
                      a single node }
  \item{edgesColored}{ logical; whether to color edges }
  \item{rankDir}{ The direction that the plot is laid out in, one of
                  either \code{"TB"} for Top-to-Bottom or \code{"LR"} for
                  Left-to-Right.  The default value is \code{"LR"}.  This
                  argument is only useful for \code{dot} layouts.}
  \item{nodeDetails,edgeDetails}{logical; whether count information
    should be shown.}
  \item{nodeSizeScore}{character; value to encode in the size of the nodes.}
  \item{edgeSizeScore}{character; value to encode in the width of the edges.}
  \item{shape}{ character; node shape.}
}
\details{
  The \code{plain.style} style corresponds to the default parameter setitngs
  in the display functions. It can be used as the basis for creating a
  new custom style.

  The \code{google.style} style is based on the display style used in the
  \code{pprof} tool from the Google Performance Tools suite.
}
\references{ \url{http://google-perftools.googlecode.com/svn/trunk/doc/cpuprofile.html}.}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{flatProfile}},
  \code{\link{summaryRprof}},
  \code{\link{readProfileData}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{printProfileCallGraph}},
  \code{\link{profileCallGraph2Dot}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
plotProfileCallGraph(pd, style = plain.style)
plotProfileCallGraph(pd, style = google.style)
}
\keyword{programming}
\keyword{utilities}
