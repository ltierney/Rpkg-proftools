\name{plotProfileCallGraph}
\alias{plotProfileCallGraph}
\title{ Plot Call Graph for Rprof Profile Data }
\description{
  Uses the graph and Rgraphviz packages to plot a call graph for
  profile data produced by \code{Rprof}.
}
\usage{
plotProfileCallGraph(pd, layout = "dot",
                     score = c("none", "total", "self"),
                     transfer = function(x) x, nodeColorMap = NULL,
                     edgeColorMap = NULL, mergeCycles = FALSE,
                     edgesColored = FALSE, rankDir = c("TB", "LR"),
                     nodeDetails = FALSE, edgeDetails = FALSE,
                     nodeSizeScore = c("none", "total", "self"),
                     edgeSizeScore = c("none", "total"),
                     shape = "ellipse", style, ...)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
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
  \item{style}{named list of values for arguments \code{score} through
               \code{layout} to use if not explicitly supplied.}
  \item{\dots}{ additional arguments for the \code{graphNEL} \code{plot}
                method.}
}
\value{Used for side effect.}
\details{
  Color is used to encode the fraction of total or self time spent in
  each function or call.  The scores used correspond to the values in
  the printed representation produced by \code{printProfileCallGraph}.
  For now, see the \code{gprof} manual for further details.  The color
  encoding for a score \code{s} and a color map \code{m} is
  \code{m[ceiling(length(m) * transfer(s))]}

  A style can be specified to set options to a set of cvalues that work
  well together.
}
\note{
  Because of lazy evaluation, nested calls like \code{f(g(x))}
  appear in the profile graph as \code{f} or one of its callees
  calling \code{g}.
}

\references{
  User manual for \code{gprof}, the GNU profiler.

  Graphviz: \url{http://www.research.att.com/sw/tools/graphviz/}
}

\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{readProfileData}},
  \code{\link{flatProfile}},
  \code{\link{profileCallGraph2Dot}}
  \code{\link{printProfileCallGraph}}
  \code{\link{plain.style}}
  \code{\link{google.style}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
plotProfileCallGraph(pd)
plotProfileCallGraph(pd, score = "total")
plotProfileCallGraph(pd, style = google.style, score = "total")
}
\keyword{programming}
\keyword{utilities}
