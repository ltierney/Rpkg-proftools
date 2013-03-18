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

aggregateCounts <- function(entries, counts) {
    dcounts <- as.data.frame(counts)
    clean <- function(x)
        if (any(is.na(x)))
            factor(as.character(x), exclude = "")
        else
            x
    fact2char(aggregate(dcounts, lapply(entries, clean), sum))
}

rbindEntries <- function(entries)
    as.data.frame(do.call(rbind, entries), stringsAsFactors = FALSE)
    
mergeCounts <- function(data, leafdata) {
    val <- merge(data, leafdata, all = TRUE)
    val$self[is.na(val$self)] <- 0
    val$gcself[is.na(val$gcself)] <- 0
    val
}

entryCounts0 <- function(pd, fun, control, names) {
    stacks <- pd$stacks
    refs <- pd$refs
    counts <- pd$counts
    gccounts <- pd$gccounts

    which <- seq_along(stacks)
    
    doLine <- function(i) fun(stacks[[i]], refs[[i]], control)
    entries <- lapply(which, doLine)
    edf <- rbindEntries(entries)

    reps <- unlist(lapply(entries, nrow))
    tot <- rep(counts, reps)
    gctot <- rep(gccounts, reps)
    ct <- cbind(tot, gctot, deparse.level = 0)
    colnames(ct) <- names
    aggregateCounts(edf, ct)
}
    
entryCounts <- function(pd, lineFun, leafFun, control) {
    aedf <- entryCounts0(pd, lineFun, control, c("total", "gctotal"))
    aledf <- entryCounts0(pd, leafFun, control, c("self", "gcself"))
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
    if (n > 1) {
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
    if (n > 1) {
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

pathSummary <- function(pd, value = c("pct", "time", "hits"), ...) {
    value = match.arg(value)
    paths <- sapply(pd$stacks, formatTrace, ...)

    ## need to aggregate in case some collapsed paths are identical
    ## or some paths differ only in source references.
    apd <- aggregateCounts(list(paths = paths),
                           cbind(counts = pd$counts, gccounts = pd$gccounts))
    apd <- apd[rev(order(apd$counts)),]

    counts <- apd$counts
    gccounts <- apd$gccounts
    paths <- apd$paths

    if (value == "pct") {
        tot <- sum(counts)
        pct <- round(100 * counts / tot, 1)
        gcpct <- round(100 * gccounts / tot, 1)
        data.frame(total.pct = pct, gc.pct = gcpct, row.names = paths)
    }
    else if (value == "time") {
        delta <- pd$interval / 1.0e6
        tm <- counts * delta
        gctm <- gccounts * delta
        data.frame(total.time = tm, gc.time = gctm, row.names = paths)
    }
    else
        data.frame(total.hits = counts, gc.hits = gccounts, row.names = paths)
}

keepIDX <- function(idx) {
    if (is.character(idx))
        which(sapply(d$stacks, function(s) any(idx %in% s)))
    else if (is.logical(idx))
        which(idx)
    else 
         idx
}

subsetPD <- function(pd, which) {
    keep <- keepIDX(which)

    pd$stacks <- pd$stacks[keep]
    pd$refs <- pd$refs[keep]
    pd$counts <- pd$counts[keep]
    pd$gccounts <- pd$gccounts[keep]

    tkeep <- which(pd$trace %in% keep)
    pd$inGC <- pd$inGC[tkeep]
    pd$trace <- pd$trace[tkeep]

    map <- match(seq_along(keep), keep)
    pd$trace <- map[pd$trace]

    pd
}

focusPD  <- function(pd, which) {
    keep <- keepIDX(which)
    nkeep <- setdiff(seq_along(pd$stacks), keep)
    
    pd$stacks[nkeep] <- list("<Other>")
    pd$refs[nkeep] <- list(c(NA_character_, NA_character_))

    compactPD(pd)
}

compactPD <- function(pd) {
    key <- mapply(c, pd$stacks, pd$refs)
    map <- match(key, unique(key))
    ct <- aggregateCounts(data.frame(key = map),
                          cbind(counts = pd$counts, gccounts = pd$gccounts))
    ct <- ct[order(ct$key),] ## may not be needed
    invmap <- match(unique(key), key)
    pd$stacks <- pd$stacks[invmap]
    pd$refs <- pd$refs[invmap]
    pd$counts <- ct$counts
    pd$gccounts <- ct$gccounts
    pd$trace <- map[pd$trace]
    pd
}

checkStackRefs <- function(val, nf) {
    s <- val$stack
    r <- val$refs
    if (length(s) == 0)
        stop("stacks must have at least one entry")
    if (length(r) != length(s) + 1)
        stop("stack and source references do not match")
    fn <- refFN(r)
    fn <- fn[! is.na(fn)]
    if (length(fn) > 0 && (min(fn) < 1 || max(fn) > nf))
        stop("invalid source references produced.")
    val
}

transformPD <- function(pd, fun) {
    stacks <- pd$stacks
    refs <- pd$refs
    nf <- length(pd$files)
    for (i in seq_along(pd$stacks)) {
        val <- checkStackRefs(fun(stacks[[i]], refs[[i]], i), nf)
        stacks[[i]] <- val$stack
        refs[[i]] <- val$refs
    }
    pd$stacks <- stacks
    pd$refs <- refs

    compactPD(pd)
}

skipIDX <- function(pd, what) {
    if (is.character(what)) {
        findFirst <- function(s) {
            idx <- match(what, s)
            if (any(! is.na(idx)))
                min(idx, na.rm = TRUE) - 1
            else
                0
        }
        idx <- sapply(pd$stacks, findFirst)
    }
    else
        idx <- ifelse(is.na(what), 0, what)
    if (length(idx) != length(pd$stacks))
        idx <- rep(idx, length = length(pd$stacks))
    idx
}
    
skipPD <- function(pd, what, merge = FALSE) {
    idx <- skipIDX(pd, what)

    skip <- function(stack, refs, i) {
        n <- idx[i]
        if (n > 0) {
            if (n < length(stack)) {
                skip <- 1 : n
                stack <- stack[-skip]
                refs <- refs[-skip]
            }
            else {
                stack <- "<Other>"
                refs <- c(NA_character_, NA_character_)
            }
        }
        else if (merge) {
            stack <- "<Other>"
            refs <- c(NA_character_, NA_character_)
        }
        list(stack = stack, refs = refs)
    }

    transformPD(pd, skip)
}

pruneIDX <- function(pd, what) {
    if (is.character(what)) {
        findPrune <- function(s) {
            idx <- match(what, s)
            if (any(! is.na(idx)))
                length(s) - min(idx, na.rm = TRUE)
            else
                0
        }
        idx <- sapply(pd$stacks, findPrune)
    }
    else
        idx <- ifelse(is.na(what), 0, what)
    if (length(idx) != length(pd$stacks))
        idx <- rep(idx, length = length(pd$stacks))
    idx
}

prunePD <- function(pd, what, merge = FALSE) {
    idx <- pruneIDX(pd, what)

    prune <- function(stack, refs, i) {
        n <- idx[i]
        if (n > 0) {
            slen <- length(stack)
            if (n < slen) {
                drop <- (slen - n + 1) : slen
                stack <- stack[-drop]
                refs <- refs[-drop]
            }
            else {
                stack <- "<Other>"
                refs <- c(NA_character_, NA_character_)
            }
        }
        list(stack = stack, refs = refs)
    }

    transformPD(pd, prune)
}
    
dl <- subsetPD(d, sapply(d$stacks, function(s) "lm.fit" %in% s))
    
d0 <- d
d0$stacks[[12]] <- d0$stacks[[13]] <- "<Other>"
d0$refs[[12]] <- d0$refs[[13]] <- c(NA_character_, NA_character_)
d0 <- compactPD(d0)

f1 <- function(stack, refs, i) {
    if (".completeToken" %in% stack)
        list(stack = "<Other>", refs = c(NA_character_, NA_character_))
    else
        list(stack = stack, refs = refs)
}

d1 <- transformPD(d, f1)
stopifnot(identical(compactPD(d0), d1))

f2 <- function(stack, refs, i) {
    path <- c("source", "withVisible", "eval", "eval")
    n <- length(path)
    idx <- 1 : n
    if (length(stack) >= n  && identical(stack[idx], path)) {
        if (length(stack) > n)
            list(stack = stack[-idx], refs = refs[-idx])
        else
            list(stack = "<Other>", refs = c(NA_character_, NA_character_))
    }
    else
        list(stack = stack, refs = refs)
}

d2 <- transformPD(d1, f2)

## **** have pathSummary optionally trim the top instead of the bottom

## **** pull path control to pathSummary

## **** need control argument in transformPD function?

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

## **** pathSummary is breaking ties alphabetically -- change with skip, maxlen

## **** use compactPD for initial pd merge?
## **** use compactPD for chunk-wise reading?

## **** would be useful if checkUsage could warn for non-namespace-globals
