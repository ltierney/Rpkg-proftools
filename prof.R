###
### Read profile data
###

readPDheader <- function(con) {
    header <- readLines(con, 1)

    pat <- ".*sample.interval=([[:digit:]]+).*"
    if (grepl(pat, header))
        interval <- as.integer(sub(pat, "\\1", header))
    else
        stop("not a valid Rprof file")

    haveMem <- grepl("memory profiling:", header)
    haveGC <- grepl("GC profiling:", header)
    haveRefs <- grepl("line profiling:", header)

    if (haveMem)
        stop("memory profiling is currently not supported")

    list(interval = interval, haveGC = haveGC, haveRefs = haveRefs)
}

readPDlines <- function(con, hdr) {
    stacks <-  readLines(con)
    if (hdr$haveRefs) {
        fstacks <- grepl("#File ", stacks)
        files <- sub("#File [[:digit:]]: (.+)", "\\1", stacks[fstacks])
        stacks <- stacks[! fstacks]
    }
    else files <- NULL

    ## record and strip out GC info
    inGC <- grepl("<GC>", stacks)
    stacks <- sub("\"<GC>\" ", "", stacks)
    
    ustacks <- unique(stacks)
    trace <- match(stacks, ustacks)

    list(files = files, stacks = ustacks, trace = trace, inGC = inGC)
}

splitStacks <- function(data) {
    stacks <- data$stacks
    rsstacks <- lapply(strsplit(stacks, " +"),
                 function(x) sub("\"(.+)\"", "\\1", rev(x)))
    refpat <- "[[:digit:]]+#[[:digit:]]+$"
    refs <- lapply(rsstacks, function(x) grepl(refpat, x))
    stacks <- lapply(rsstacks, function(x) x[! grepl(refpat, x)])
    stackrefs <- lapply(rsstacks, function(x) {
        isref <- grepl(refpat, x)
        stack <- x[! isref]
        refs <- rep(NA_character_, length(stack) + 1)
        k <- 1
        i <- 1
        n <- length(x)
        while (i <= n) {
            if (isref[i]) {
                refs[k] <- x[i]
                i <- i + 2
            }
            else
                i <- i + 1
            k <- k + 1
        }
        refs
    })
    data$stacks <- stacks
    data$refs <- stackrefs
    data
}

countStacks <- function(trace, inGC)
    list(counts = as.numeric(table(trace)),
         gccounts = as.numeric(tapply(inGC, trace, sum)))

readPD <- function(file) {
    con <- file(file, "r")
    on.exit(close(con))
    
    hdr <- readPDheader(con)
    data <- readPDlines(con, hdr)
    sdata <- splitStacks(data)
    counts <- countStacks(sdata$trace, sdata$inGC)
    c(hdr, sdata, counts)
}


###
### Flame graph and time graph
###

flameGraph <- function(stacks, counts, reorder = TRUE) {
    mx <- max(sapply(stacks, length))

    ## For 'standard' flame graph order the stacks so they are
    ## alphabetical within lines within calls, with missing entires
    ## frst. This does a lexicographic sort by sorint on the top entry
    ## first, then the next, and do on; since the sorts are stable
    ## this keeps the top levels sorted within the lower ones.
    if (reorder) {
        ord <- seq_along(stacks)
        for (i in (mx : 1))
            ord <- ord[order(sapply(stacks[ord], `[`, 1), na.last = FALSE)]
        stacks <- stacks[ord]
        counts <- counts[ord]
    }

    plot(c(0, sum(counts)), c(0, mx), type = "n",
         axes = FALSE, xlab = "", ylab = "")

    ## half-em for half-character offset used by text() when pos
    ## argument is used.
    hm <- 0.5 * strwidth("m")

    for (k in 1 : mx) {
        runs <- rle(sapply(stacks, `[`, k))
        lens <- runs$lengths
        n <- length(lens)
        csums <- cumsum(tapply(counts, rep(1 : n, lens), sum))

        top <- k
        bottom <- k - 1
        left <- 0
        right <- 0

        ## should be able to vectorize this loop
        for (i in 1 : n) {
            left <- right
            right <- csums[i]
            label <- runs$values[i]
            if (! is.na(label))
                rect(left, bottom, right, top,
                     col = rgb(runif(1), runif(1), runif(1)))
            if (strheight(label) <= 1 &&
                strwidth(label) + 2 * hm <= 0.8 * (right - left))
                text(left, bottom + 0.4, label, pos = 4)
        }
    }
}

flameGraph <- function(stacks, counts, reorder = TRUE) {
    mx <- max(sapply(stacks, length))

    ## For 'standard' flame graph order the stacks so they are
    ## alphabetical within lines within calls, with missing entires
    ## frst. This does a lexicographic sort by sorint on the top entry
    ## first, then the next, and do on; since the sorts are stable
    ## this keeps the top levels sorted within the lower ones.
    if (reorder) {
        ord <- seq_along(stacks)
        for (i in (mx : 1))
            ord <- ord[order(sapply(stacks[ord], `[`, 1), na.last = FALSE)]
        stacks <- stacks[ord]
        counts <- counts[ord]
    }

    plot(c(0, sum(counts)), c(0, mx), type = "n",
         axes = FALSE, xlab = "", ylab = "")

    ## half-em for half-character offset used by text() when pos
    ## argument is used.
    hm <- 0.5 * strwidth("m")

    for (k in 1 : mx) {
        runs <- rle(sapply(stacks, `[`, k))
        lens <- runs$lengths
        n <- length(lens)
        csums <- cumsum(tapply(counts, rep(1 : n, lens), sum))

        top <- k
        bottom <- k - 1
        left <- c(0, csums[-n])
        right <- csums
        cols <- rgb(runif(n), runif(n), runif(n))
        labels <- runs$values

        show <- ! is.na(labels)
        if (any(show))
            rect(left[show], bottom, right[show], top, col = cols[show])

        show <- (! is.na(labels) &
                strheight(labels) <= 1 &
                strwidth(labels) + 2 * hm <= 0.8 * (right - left))
        if (any(show))
            text(left[show], bottom + 0.4, labels[show], pos = 4)
    }
}

## produce a flame graph from an Rprof file
fg <- function(file) {
    d <- readPD(file)
    flameGraph(d$stacks, d$counts)
}

## produce a time graph (like profr) from an Rprof file
tg <- function(file) {
    d <- readPD(file)
    r <- rle(d$trace)
    stacks <- d$stacks[r$values]
    counts <- r$lengths
    flameGraph(stacks, counts, FALSE)
}    


###
### Function and call summaries
###

fact2char <- function(d) {
    for (i in seq_along(d))
        if (is.factor(d[[i]]))
            d[[i]] <- as.character(d[[i]])
    d
}

aggregateCounts <- function(entries, ...) {
    counts <- data.frame(...)
    clean <- function(x)
        if (any(is.na(x)))
            factor(as.character(x), exclude = "")
        else
            x
    fact2char(aggregate(counts, lapply(entries, clean), sum))
}

rbindEntries <- function(entries)
    as.data.frame(do.call(rbind, entries), stringsAsFactors = FALSE)
    
mergeCounts <- function(data, leafdata) {
    val <- merge(data, leafdata, all = TRUE)
    val$self[is.na(val$self)] <- 0
    val$gcself[is.na(val$gcself)] <- 0
    val
}

entryCounts <- function(pd, lineFun, leafFun, control) {
    stacks <- pd$stacks
    refs <- pd$refs
    counts <- pd$counts
    gccounts <- pd$gccounts

    which <- seq_along(stacks)
    
    doLine <- function(i) lineFun(stacks[[i]], refs[[i]], control)
    entries <- lapply(which, doLine)
    edf <- rbindEntries(entries)

    reps <- unlist(lapply(entries, nrow))
    tot <- rep(counts, reps)
    gctot <- rep(gccounts, reps)
    aedf <- aggregateCounts(edf, total = tot, gctotal = gctot)

    doLeaf <- function(i) leafFun(stacks[[i]], refs[[i]], control)
    ledf <- rbindEntries(lapply(which, doLeaf))
    aledf <- aggregateCounts(ledf, self = counts, gcself = gccounts)

    mergeCounts(aedf, aledf)
}

lineFuns <- function(line, refs, useSite) {
    if (useSite) {
        n <- length(line)
        site <- refs[-(n + 1)]
    }
    else site <- NA_character_
    unique(cbind(fun = line, site))
}

leafFun <- function(line, refs, useSite) {
    n <- length(line)
    fun <- line[n]
    site <- if (useSite) refs[n] else NA_character_
    cbind(fun, site)
}

funCounts <- function(pd, useSite = TRUE)
    entryCounts(pd, lineFuns, leafFun, useSite)
    
lineCalls <- function(line, refs, cntrl) {
    n <- length(line)
    if (n > 0) {
        caller <- line[-n]
        callee <- line[-1]
        if (cntrl$useCalleeSite)
            callee.site <- refs[-c(1, n + 1)]
        else
            callee.site <- NA_character_
        if (cntrl$useCallerSite)
            caller.site <- refs[-c(n, n + 1)]
        else
            caller.site <- NA_character_
    }
    else
        caller <- callee <- callee.site <- caller.site <- character()
    unique(cbind(caller, callee, caller.site, callee.site))
}

leafCall <- function(line, refs, cntrl) {
    n <- length(line)
    if (n > 0) {
        caller <- line[n - 1]
        callee <- line[n]
        if (cntrl$useCalleeSite)
            callee.site <- refs[n]
        else
            callee.site <- NA_character_
        if (cntrl$useCallerSite)
            caller.site <- refs[n - 1]
        else
            caller.site <- NA_character_
    }
    else
        caller <- callee <- callee.site <- caller.site <- character()
    cbind(caller, callee, caller.site, callee.site)
}

callCounts <- function(pd, useCalleeSite = TRUE, useCallerSite = FALSE) {
    cntrl <- list(useCalleeSite = useCalleeSite, useCallerSite = useCallerSite)
    entryCounts(pd, lineCalls, leafCall, cntrl)
}


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
## defined.  Otherwise tread the caller's home file as unkown. The
## result returned by this function is a named vector with one element
## per funciton for which the home file is assumed know.  The names
## are the names of the callers, and the values are the indices of the
## files in whicn the callers are defined.
homeFileMap <- function(cc) {
    map <- tapply(cc$callee.site, cc$caller, commonFile)
    map <- map[! is.na(map)]
    map
}

## Extract the file indices and line numbers from source references of
## the form FN#LN.
refFN <- function(refs)
    as.integer(sub("([[:digit:]]+)#[[:digit:]]+", "\\1", refs))

refLN <- function(refs)
    as.integer(sub("[[:digit:]]+#([[:digit:]]+)", "\\1", refs))

## Collect the data for the callgrind output. The basic data is
##
##     fc = function counts
##     cc = call counts and call site references
##
## To fc we add the indes of the xome file for each cunction (NA if
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
    fc <- funCounts(pd, FALSE)
    cc <- callCounts(pd, TRUE, FALSE)

    hfm <- homeFileMap(cc)

    fc$fl <- hfm[match(fc$fun, names(hfm))]
    cc$cfl <- hfm[match(cc$callee, names(hfm))]
    cc$cln <- ifelse(is.na(match(cc$caller, names(hfm))),
                     NA, refLN(cc$callee.site))

    if (! GC)
        fc$gcself <- 0

    list(fc = fc, cc = cc, gcself = sum(fc$gcself), files = pd$files)
}

writeSelfEntry <- function(con, fun, fc, files) {
    fn <- fc$fl[fc$fun == fun]
    file <- if (is.na(fn)) "??" else files[fn]
    self <- fc$self[fc$fun == fun]
    gcself <- fc$gcself[fc$fun == fun]

    cat(sprintf("\nfl=%s\nfn=%s\n0 %d\n",  file, fun, self - gcself),
        file = con)

    if (gcself > 0)
        cat(sprintf("cfn=<GC>\ncalls=%d 0\n0 %d\n", gcself, gcself),
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
        cat(sprintf("\nfn=<GC>\n0 %d\n", gcself), file = con)
}

writeCG <- function(con, pd, GC = TRUE) {
    if (is.character(con)) {
        con <- file(con, "w")
        on.exit(close(con))
    }

    data <- getCGdata(pd, GC)
    
    cat("events: Hits\n", file = con)

    for (fun in data$fc$fun)
        writeFunEntries(con, fun, data)

    writeGCEntry(con, data)
}


###
### Source reference summaries
###

lineSites <- function(line, refs, useSite)
    unique(cbind(site = refs))

leafSite <- function(line, refs, useSite) {
    n <- length(refs)
    cbind(site = refs[n])
}

siteCounts <- function(pd)
    entryCounts(pd, lineSites, leafSite, TRUE)


###
### Experimental stuff
###


countHits <- function(stacks, counts) {
    stacks <- lapply(stacks, function(x) x[! is.na(x)])
    uitems <- unique(unlist(stacks))
    ln <- sapply(uitems,
                 function(y) sapply(stacks,
                                    function(x) as.integer(y %in% x)))
    t(ln) %*% do.call(cbind, counts)
}

countSelfHits <- function(stacks, counts) {
    uitems <- unique(unlist(lapply(stacks, function(x) unique(x[! is.na(x)]))))
    ln <- sapply(uitems,
                 function(y) sapply(stacks,
                                    function(x) {
                                        n <- length(x)
                                        if (n > 0 && identical(y, x[n]))
                                            1
                                        else
                                            0
                                    }))
    t(ln) %*% do.call(cbind, counts)
}

recodeRefs <- function(refs, files, na.value = NA) {
    fn <- refFN(refs)
    ln <-refLN(refs)
    ifelse(is.na(refs), na.value, paste(basename(files)[fn], ln, sep = ":"))
}

recodeRefsList <- function(refs, files)
    sapply(refs, recodeRefs, files)

formatTrace <- function(trace, maxlen = 50, skip = 0) {
    if (skip > 0)
        trace <- trace[-(1 : skip)]
    out <- paste(trace, collapse = " -> ")
    while (nchar(out) > maxlen && length(trace) > 1) {
        trace <- trace[-1]
        out <- paste("... ", paste(trace, collapse = " -> "))
    }
    out
}

d <- readPD("Rprof-lmfit-new.out")
cts <- list(counts = d$counts, gccounts = d$gccounts)
countHits(d$stacks, cts)

countSelfHits(d$stacks, cts)
countSelfHits(d$refs, cts)

## this drops the final ref
v <- mapply(function (x, y)
            ifelse(is.na(y), x, paste0(x, " (", y,")")),
            d$stacks, recodeRefsList(d$refs, d$files))
countHits(v, cts)

## this attributes the final ref to "<Unknown>"
v <- mapply(function (x, y) {
            x <- c(x, "<Unknown>")
            ifelse(is.na(y), x, paste0(x, " (", y,")"))
            },
            d$stacks, recodeRefsList(d$refs, d$files))
countHits(v, cts)

printPaths <- function(pd, n, ...) {
    ord = rev(order(pd$counts))
    if (! missing(n) && length(ord) > n)
        ord <- ord[1 : n]
    tot <- sum(pd$counts)
    pct <- round(100 * pd$counts[ord] / tot, 1)
    gcpct <- round(100 * pd$gccounts[ord] / tot, 1)
    paths <- sapply(pd$stacks[ord], formatTrace, ...)
    mapply(function(x, y, z) cat(sprintf("%5.1f %5.1f   %s\n", x, y, z)),
           pct, gcpct, paths)
    invisible(NULL)
}

pathSummary <- function(pd, ...) {
    paths <- sapply(pd$stacks, formatTrace, ...)
    apd <- aggregateCounts(list(paths = paths),
                           counts = pd$counts, gccounts = pd$gccounts)
    apd <- apd[rev(order(apd$counts)),]
    tot <- sum(pd$counts)
    pct <- round(100 * apd$counts / tot, 1)
    gcpct <- round(100 * apd$gccounts / tot, 1)
    val <- data.frame(total.pct = pct, gc.pct = gcpct)
    rownames(val) <- apd$paths
    val
}

## **** option of pct, hits, time

## **** subset -- only stacks containing X or starting with Y
## **** merge -- cycles or subtrees

## **** Hot paths -- print nicely, but also allow examination of
## **** source refs and such.

## **** Maybe print with source refs?
## **** Reorder paths, revise time index, so hottest one is first?

## **** flexible ways of examining/displaying fc, cc?

## **** show hot files, hot lines, calls within line
## **** show hot paths -- tree view with collapsing of some kind
## **** show call tree

## **** figure out how to generate call graphs as in proftools
## **** figure out how to generate gprof output
## **** allow pct, counts, or time in final output
## **** use [...] instead of <...> for GC and Anonymous?

## **** would be useful if checkUsage could warn for non-namespace-globals

