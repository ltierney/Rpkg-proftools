\name{printProfileCallGraph}
\alias{printProfileCallGraph}
\title{ Print Call Graph for Rprof Profile Data }
\description{
  Prints a representation of the call graph for profile data produced
  by \code{Rprof}.  Output can be directed to a connection or a file.
}
\usage{
printProfileCallGraph(pd, file = stdout(), percent = TRUE)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{file}{ a connection or the name of the file where the profile
               graph will be written. }
  \item{percent}{ logical; if true use percent of total time; otherwise
                  use time in seconds }
}
\value{Used for side effect.}
\details{

  \code{printProfileCallGraph} produces a printed representation of
  the call graph for profile data produced by \code{Rprof}.  The
  representation is analogous to the call graph produced by
  \code{gprof} with a few minor changes.  eventually more complete
  documentation of the format will be provided here; for now, reading
  the \code{gprof} manual section on the call graph should help
  understanding this output.  The output is similar enough to
  gprof output for the \code{cgprof} script to be able to produce a
  visual representation of the call graph via Graphviz.
}
\note{
  Because of lazy evaluation, nested calls like \code{f(g(x))}
  appear in the profile graph as \code{f} or one of its callees
  calling \code{g}.
}
\references{
  User manual for \code{gprof}, the GNU profiler.

  \code{cgprof}: \url{http://mvertes.free.fr/}

  Graphviz: \url{http://www.research.att.com/sw/tools/graphviz/}
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
printProfileCallGraph(pd)
\dontrun{
## If you have graphviz and cgprof installed on a UNIX-like system
## then in R do:

pd <- readProfileData(system.file("samples", "glmEx.out", package="proftools"))
printProfileCallGraph(pd, "foo.graph")

## and then in a shell do (to use the interactive dotty):

cgprof -TX foo.graph

## or (to create a postscript version and view with gv):

cgprof -Tps foo.graph > foo.ps
gv foo.ps
  }
}
\keyword{programming}
\keyword{utilities}
