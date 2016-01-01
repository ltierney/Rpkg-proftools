\name{writeCallgrindFile}
\alias{writeCallgrindFile}
\title{ Write Out Profile Data in Callgrind Format }
\description{
  Writes the profile data in callgrind format suitable for use with
  \code{kcachegrind} or \code{qcachegrind}.
}
\usage{
writeCallgrindFile(pd, file = "Rprof.cg", GC = TRUE, dropSource = TRUE)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{file}{ a connection or the name of the file where the callgrind
               output will be written. }
  \item{GC}{ logical; if true include GC information. }
  \item{dropSource}{ logical; if true drop initial stack entried from a
    \code{source} call and add a top level marker. }
}
\value{Used for side effect.}
\details{

  The callgrind format is used by Valgrind's callgrind tool. The KDE
  tool \code{kcachegrind} can be used to display the file;
  \code{kcachegrind} displays the summary statistics, a call graph, and
  annotated source code if source information is
  available. \code{kcachegrind} is available in Linux and Windows; on
  Mac OSX \code{qcachegrind} is available.
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
tmp <- tempfile()
writeCallgrindFile(pd, file = tmp)
file.show(tmp)
unlink(tmp)
\dontrun{
## If you have kcachegrind installed on a UNIX-like system then do:
tmp <- tempfile()
writeCallgrindFile(pd, file = tmp)
system(sprintf("kcachegrind %s", tmp))
unlink(tmp)
  }
}
\keyword{programming}
\keyword{utilities}
