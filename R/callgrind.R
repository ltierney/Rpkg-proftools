###
### Writing callgrind file
###

## For a set of source references determine if they have a common file
## index and return that index. If they do not have a common index
## then return NA.
commonFile <- function(refs) {
    fn <- unique(refFN(refs))
    if (length(fn) == 1 && ! is.na(2))
        fn
    else
        NA
}

## For each caller check whether all calls have a common file index.
## If they do, then assume this is the file in which the caller is
## defined.  Otherwise treat the caller's home file as unkown. The
## result returned by this function is a named vector with one element
## per function for which the home file is assumed known.  The names
## are the names of the callers, and the values are the indices of the
## files in whicn the callers are defined.  For leaf calls NA sites
## are ignored. Possibly a disagreement of the leaf call site with
## other calls should be ignored as well.
homeFileMap <- function(pd, cc) {
    lsites <- leafCallRefs(pd)
    site <- c(cc$callee.site, lsites$site)
    caller <- c(cc$caller, lsites$fun)
    map <- tapply(site, caller, commonFile)
    map[! is.na(map)]
}

leafCallRefs <- function(pd) {
    ln <- sapply(pd$stacks, length)
    stacks <- pd$stacks[ln > 0]
    refs <- pd$refs[ln > 0]
    lfuns <- sapply(stacks, function(x) x[length(x)])
    lrefs <- sapply(refs, function(x) x[length(x)])
    goodrefs <- ! is.na(lrefs)
    list(fun = lfuns[goodrefs], site = lrefs[goodrefs])
}

## Collect the data for the callgrind output. The basic data is
##
##     fc = function counts and leaf call references
##     cc = call counts and call site references
##
## To fc we add the index of the home file for each function (NA if
## not known) in fl.
##
## To cc we add in cfl the index of the home file of the function
## called (the callee), and in cln the line number of the call (in the
## caller's file). If the caller's file is considered unknown, then
## the line number is NA.
##
## If we do not want GC information in the output then we set the
## gcself entries in fc to zero, since the output functions only
## generate GC output for positive gcself counts.
getCGdata <- function(pd, GC) {
    fc <- getCGselfData(pd)
    cc <- callCounts(pd, TRUE, FALSE)

    hfm <- homeFileMap(pd, cc)

    fc$fl <- hfm[match(fc$fun, names(hfm))]
    cc$cfl <- hfm[match(cc$callee, names(hfm))]
    cc$cln <- ifelse(is.na(match(cc$caller, names(hfm))),
                     NA, refLN(cc$callee.site))

    if (! GC)
        fc$gcself <- 0

    list(fc = fc, cc = cc, gcself = sum(fc$gcself),
         funs = sort(unique(fc$fun)),
         files = pd$files)
}

getCGselfData <- function(pd) {
    leafFun <- function(x) {
        n <- length(x)
        if (n > 0)
            x[n]
        else
            "<TOP>"
    }
    fc <- funCounts(pd, FALSE)
    fc <- fc[, c("fun", "site", "self", "gcself")]
    f <- sapply(pd$stacks, leafFun)
    r <- sapply(pd$refs, function(x) x[length(x)])
    fs <- aggregateCounts(data.frame(fun = f, site = r),
                          data.frame(self = pd$counts, gcself = pd$gccounts))
    rbind(fs, fc[fc$self == 0, ])
}

writeSelfEntry <- function(con, fun, fc, files) {
    fn <- fc$fl[fc$fun == fun][1]
    file <- if (is.na(fn)) "??" else files[fn]
    self <- fc$self[fc$fun == fun]
    gcself <- fc$gcself[fc$fun == fun]
    site <- fc$site[fc$fun == fun]
    line <- ifelse(is.na(site), 0, refLN(site))

    cat(sprintf("\nfl=%s\nfn=%s\n", file, fun), file = con)
    cat(sprintf("%d %d\n", line, self - gcself), sep = "", file = con)

    gcself <- sum(gcself)
    if (gcself > 0)
        cat(sprintf("cfl=??\ncfn=<GC>\ncalls=%d 0\n0 %d\n", gcself, gcself),
            sep = "", file = con)
}
    
writeCallEntries <- function(con, fun, cc, files) {
    fcc <- cc[cc$caller == fun, ]
    cfun <- fcc$callee
    tot <- fcc$total
    file <- ifelse(is.na(fcc$cfl), "??", files[fcc$cfl])
    line <- ifelse(is.na(fcc$cln), 0, fcc$cln)

    cat(sprintf("cfl=%s\ncfn=%s\ncalls=%d 0\n%d %d\n",
                file, cfun, tot, line, tot),
        sep = "", file = con)
}

writeFunEntries <- function(con, fun, data) {
    fc <- data$fc
    cc <- data$cc
    files <- data$files
    writeSelfEntry(con, fun, fc, files)
    writeCallEntries(con, fun, cc, files)
}

writeGCEntry <- function(con, data)  {
    gcself <- data$gcself
    if (gcself > 0)
        cat(sprintf("\nfl=??\nfn=<GC>\n0 %d\n", gcself), file = con)
}

writeCG <- function(con, pd, GC = TRUE) {
    if (is.character(con)) {
        con <- file(con, "w")
        on.exit(close(con))
    }

    data <- getCGdata(pd, GC)
    
    cat("events: Hits\n", file = con)

    for (fun in data$funs)
        writeFunEntries(con, fun, data)

    writeGCEntry(con, data)
}

writeCallgrindFile <- function(pd, file = "Rprof.cg", GC = TRUE,
                               dropSource = TRUE) {
    if (dropSource)
        pd <- dropSourceFrames(pd)
    writeCG(file, pd, GC)
}
