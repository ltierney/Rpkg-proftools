\name{summaries}
\alias{funSummary}
\alias{callSummary}
\alias{pathSummary}
\alias{srcSummary}
\title{Basic Profile Data Summaries}
\description{
  Functions to summarize profile data.
}
\usage{
funSummary(pd, byTotal = TRUE, value = c("pct", "time", "hits"),
           srclines = TRUE, GC = TRUE, self.pct = 0, total.pct = 0)
callSummary(pd, byTotal = TRUE, value = c("pct", "time", "hits"),
            srclines = TRUE, GC = TRUE, self.pct = 0, total.pct = 0)
pathSummary(pd, value = c("pct", "time", "hits"),
            srclines = FALSE, GC = TRUE, total.pct = 0, ...)
srcSummary(pd, byTotal = TRUE, value = c("pct", "time", "hits"),
           GC = TRUE, total.pct = 0, source = TRUE, width = getOption("width"))
}
\arguments{
  \item{pd}{profile data as returned by \code{readProfileData}.}
  \item{byTotal}{logical; sort by total or by self time.}
  \item{value}{character; show result as percentage, time, or hits.}
  \item{srclines}{logical; include source information, if available, or not.}
  \item{GC}{logical; include GC information or not.}
  \item{self.pct}{numeric; functions at the bottom of the stacks with
                  self percentages below this value are removed.}
  \item{total.pct}{numeric; functions at the top of the stacks with
                   total percentages below this value are removed.}
  \item{source}{logical; if true, include source lines if available.}
  \item{width}{maximal line length to use; source lines are abbreviated
    to fit if necessary.}
  \item{...}{arguments to control path formatting; not useful at this
    time.}
}
\value{
  A data frame of summary information.
}
\details{
  \code{funSummary} returns a summary of the time spent in each
  function, or each call site if source information is available and
  requested. It is similar to \code{flatProfile}.

  \code{callSummary} provides a breakdown by calls, again with an option
  of distinguishing call and callee sites if source information is
  available.

  \code{pathSummary} returns a summary of time spent in each unique call
  path contained in the profile data.

  For profile data containing source information \code{srcSummary}
  returns a summary of time spent in each file line recorded in the
  profile data.
}
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
funSummary(pd)
callSummary(pd)
pathSummary(pd)
}
\keyword{programming}
\keyword{utilities}
