\name{callTreeGraphs}
\alias{flameGraph}
\alias{calleeTreeMap}
\alias{identify.proftools_calleeTreeMap}
\alias{identify.proftools_flameGraph}
\title{ Call Tree Graphs }
\description{
  Produce a flame graph or a callee tree map for the call tree in a profile
  stack trace.
}
\usage{
flameGraph(pd, svgfile, order = c("hot", "alpha", "time"),
           colormap = NULL, srclines = FALSE, cex = 0.75,
           main = "Call Graph", tooltip = FALSE)

calleeTreeMap(pd, srclines = FALSE, cex = 0.75, colormap = NULL,
              main = "Callee Tree Map", squarify = FALSE, border = NULL)

\method{identify}{proftools_calleeTreeMap}(x, n = 1, print = FALSE, ...)

\method{identify}{proftools_flameGraph}(x, n = 1, print = FALSE,
                                        outline = FALSE, ...)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{svgfile}{character; name for SVG output file.}
  \item{order}{character; see details below.}
  \item{colormap}{a function or \code{NULL}; see details below.}
  \item{srclines}{logical; include source information, if available, or not.}
  \item{cex}{numeric character expansion value.}
  \item{main}{character; plot title.}
  \item{tooltip}{logical; whether SVG should show details in tooltips.}
  \item{border}{character or \code{NULL}; border color for rectangles.}
  \item{squarify}{logical; whether the squarified tiling algorithm should
    be used.}
  \item{x}{\code{flameGraph} of \code{calleeTreeMap} object to use.}
  \item{n}{integer; number of items to identify.}
  \item{print}{logical; whether to print result on each click.}
  \item{outline}{logical; whether to outline rectangles corresponding to
    identified call.}
  \item{...}{Further arguments for the identify methods; currently ignored.}
}
\value{Objects that can be used with \code{identify}.}
\details{
\code{calleeTreeMap} shows a tree map of the calls in a stack trace's
call tree. The tiling algorithm used depends on the \code{squarify}
argument.  If \code{squarify} is \code{TRUE} then the
\emph{squarified} algorithm is used; otherwise, the longer side is
partitioned.

\code{flameGraph} produces a \emph{flame graph} of the call tree. The
vertical positions of rectangles represent call depth on the
stack. The widths of the rectangles represent the amount of time spent
in a call at a particular call or set of calls at a particular depth.
The \code{order} argument determines the ordering of call rectangles
at a particular level within a call at the lower level. The
\code{"alpha"} ordering orders the calls alphabetically; \code{"hot"}
uses the hot path ordering with the call with the largest amount of
time first. The \code{"time"} ordering preserves the original time
ordering within the stack trace file.

Default colors are based on the \code{rainbow} palette. Alternative
colors can be specified by a \code{colormap} function. This function
is called with three argument vectors, with one element for each
rectangle. The arguments are the function name, call depth, and number
of hits.

The results returned by \code{flameGraph} and \code{calleeTreeMap} can
be passed to \code{identify}. The \code{identify} method for
\code{calleeTreeMap} returns a list of the call stacks for the
identified rectangles. The \code{identify} method for \code{flameGraph}
returns a character vector of the labels of the identified rectangles.
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{flatProfile}},
  \code{\link{filterProfileData}},
  \code{\link{readProfileData}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{profileCallGraph2Dot}},
  \code{\link{hotPaths}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
flameGraph(pd)
calleeTreeMap(pd)
}
\keyword{programming}
\keyword{utilities}
