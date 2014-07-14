\name{annotateSource}
\alias{annotateSource}
\title{ Annotate Soruce Files }
\description{
  Annotates source files with profile information.
}
\usage{
annotateSource(pd, value = c("pct", "time", "hits"), GC = TRUE,
               sep = ":  ", show = TRUE, ...)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{value}{character; show result as percentage, time, or hits.}
  \item{GC}{logical; include GC information or not.}
  \item{sep}{character; separator between profile info and source
    lines.}
  \item{show}{logical; if true, show files with \code{file.show}.}
  \item{...}{additional arguments for \code{file.show}.}	   
}
\value{A list of character vectors of the annotated file lines.}
\details{
For lines that appear in the stack trace the percent time and,
optionally, GC time are shown before each line.
}
\author{ Luke Tierney }

\seealso{
  \code{\link{Rprof}},
  \code{\link{summaryRprof}},
  \code{\link{flatProfile}},
  \code{\link{filterProfileData}},
  \code{\link{readProfileData}},
  \code{\link{srcSummary}}
}
\examples{
## This defines the function rw()
source(system.file("samples", "rw.R", package="proftools"))

## Execute the function and collect profile data
Rprof(tmp <- tempfile(), gc.profiling = TRUE, line.profiling = TRUE)
w <- rw(200000)
Rprof(NULL)
pd <- readProfileData(tmp)
unlink(tmp)

## Annotate the sources
annotateSource(pd)
}
\keyword{programming}
\keyword{utilities}
