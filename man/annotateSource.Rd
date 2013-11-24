\name{annotateSource}
\alias{annotateSource}
\title{ Annotate Soruce Files }
\description{
  Annotates source files with profile information.
}
\usage{
annotateSource(pd, GC = TRUE, sep = ":  ", show = TRUE, ...)
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{GC}{logical; include GC information or not.}
  \item{sep}{character; separator between s profile info and source
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
\keyword{programming}
\keyword{utilities}
