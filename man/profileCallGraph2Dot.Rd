\name{profileCallGraph2Dot}
\alias{profileCallGraph2Dot}
\title{ Write Call Graph for Rprof Profile Data to Graphviz Dot File }
\description{
   Prints a Graphviz .dot file representation of the call graph for profile
   data produced by \code{Rprof}.
}
\usage{
profileCallGraph2Dot(pd, score = c("none", "total", "self"),
                     transfer = function(x) x, nodeColorMap = NULL,
                     edgeColorMap = NULL, filename = "Rprof.dot",
                     landscape = FALSE, mergeCycles = FALSE,
                     edgesColored = FALSE,
                     rankDir = c("TB", "LR"),
                     nodeDetails = FALSE, edgeDetails = FALSE,
                     nodeSizeScore = c("none", "total", "self"),
                     edgeSizeScore = c("none", "total"),
                     center = FALSE, size, shape = "ellipse",
                     layout = "dot", style, GC = TRUE,
                     maxnodes = NA, total.pct = 0)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{score}{ character string specifying whether to use total time
                or self time for coloring nodes/edges; no color used if
                missing. }
  \item{transfer}{ function; maps score values in unit interval to unit
                   interval }
  \item{nodeColorMap,edgeColorMap}{ character vectors of color
    specifications as produced by \code{rainbow}; \code{transfer} of
    score is mapped to color }
  \item{filename}{ name of \code{.dot} file }
  \item{landscape}{ logical; whether to add the \code{rotate=90} option to
                    the \code{.dot} file }
  \item{mergeCycles}{ logical; whether to merge each cycle of recursion into
                      a single node }
  \item{edgesColored}{ logical; whether to color edges }
  \item{rankDir}{ character; value to use for the \code{rankdir=} option to
                  specify the direction that the plot is laid out using
                  the \code{dot} layout; must be either \code{"TB"} for
                  Top-to-Bottom or \code{"LR"} for Left-to-Right.  The
                  default value is \code{"LR"}.}
  \item{nodeDetails,edgeDetails}{logical; whether count information
    should be shown.}
  \item{nodeSizeScore}{character; value to encode in the size of the nodes.}
  \item{edgeSizeScore}{character; value to encode in the width of the edges.}
  \item{center}{ logical; whether to add the \code{center=1} option to
                 the \code{.dot} file. }
  \item{size}{ character; string to add as \code{size=} option in the
               \code{.dot} file. }
  \item{shape}{ character; node shape.}
  \item{layout}{ character; layout method to use.}
  \item{style}{named list of values for arguments \code{score} through
               \code{layout} to use if not explicitly supplied.}
  \item{GC}{logical; include GC information or not.}
  \item{maxnodes}{integer; maximal number of nodes to use; nodes with
                  lower total hit counts are dropped.}
  \item{total.pct}{numeric; if positive, nodes with hit percentages
    below this level are dropped.}	      
}
\value{Used for side effect.}
\details{
  Writes the call graph as a Graphviz \code{.dot} file. Color is used
  to encode the fraction total or self time spent in each function or
  call.  The scores used correspond to the values in the printed
  representation produced by \code{printProfileCallGraph}.  For now,
  see the \code{gprof} manual for further details.  The color encoding
  for a score \code{s} and a color map \code{m} is
  \code{ceiling(length(m) * transfer(s))}

  A style can be specified to set options to a set of cvalues that work
  well together.
}
\note{
  Because of lazy evaluation, nested calls like \code{f(g(x))}
  appear in the profile graph as \code{f} or one of its callees
  calling \code{g}.

  Calls to functions with names containing a \code{|} character are
  dropped as they cause problems in the \pkg{graph} package.
}

\references{
  User manual for \code{gprof}, the GNU profiler.

  Graphviz: \url{https://graphviz.gitlab.io/download/}
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{readProfileData}},
  \code{\link{flatProfile}},
  \code{\link{plotProfileCallGraph}},
  \code{\link{printProfileCallGraph}}
  \code{\link{plain.style}}
  \code{\link{google.style}}
}
\examples{
pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
tmp <- tempfile()
profileCallGraph2Dot(pd, filename = tmp)
file.show(tmp)
unlink(tmp)

\dontrun{
## If you have graphviz installed on a UNIX-like system
## then in R do:

tmp.dot <- tempfile()
tmp.pdf <- tempfile()

profileCallGraph2Dot(pd, filename = tmp.dot)
system(sprintf("dot -Tpdf -o%s %s", tmp.pdf, tmp.dot))
browseURL(sprintf("file://%s", tmp.pdf))

profileCallGraph2Dot(pd, filename = tmp.dot)
system(sprintf("dot -Tpdf -o%s %s", tmp.pdf, tmp.dot))
browseURL(sprintf("file://%s", tmp.pdf))

unlink(tmp.dot)
unlink(tmp.pdf)
  }
}
\keyword{programming}
\keyword{utilities}
