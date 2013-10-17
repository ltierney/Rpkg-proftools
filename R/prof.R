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

readProfileData <- function(filename = "Rprof.out")
    readPD(filename)

readPDlines <- function(con, hdr) {
    stacks <-  readLines(con)
    if (hdr$haveRefs) {
        fstacks <- grepl("#File ", stacks)
        files <- sub("#File [[:digit:]]: (.+)", "\\1", stacks[fstacks])
        stacks <- stacks[! fstacks]
    }
    else files <- NULL

    ## remove any lines with only a <GC> entry
    onlyGC <- grepl("^\"<GC>\" $", stacks)
    if (any(onlyGC))
        stacks <- stacks[! onlyGC]

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

countStacks <- function(trace, inGC) {
    counts <- as.numeric(table(trace))
    gccounts <- as.numeric(tapply(inGC, trace, sum))
    list(counts = counts, gccounts = gccounts, total = sum(counts))
}

readPD <- function(file) {
    con <- file(file, "r")
    on.exit(close(con))
    
    hdr <- readPDheader(con)
    data <- readPDlines(con, hdr)
    sdata <- splitStacks(data)
    counts <- countStacks(sdata$trace, sdata$inGC)
    structure(c(hdr, sdata, counts), class =  "proftools_profData")
}

print.proftools_profData <- function(x, n = 6, ...) {
    cat(sprintf("Samples: %d\n", x$total))
    cat(sprintf("Time:    %.2f\n", x$total * (x$interval / 1000000)))
    fs <- funSummary(x)
    print(head(fs, n), ...)
    if (n < nrow(fs)) cat("...\n")
}


###
### Operations on profile data
###

filterProfileData <- function(pd, select, omit,
                              maxdepth = NA,
                              self.pct = 0, total.pct = 0,
                              focus = FALSE,
                              normalize = FALSE,
                              regex = FALSE) {
    if (! is.na(maxdepth))
        pd <- prunePD(pd, to = maxdepth)
    pd <- focusPD(pd, self.pct = self.pct, total.pct = total.pct)
    if (! missing(select))
        if (focus)
            pd <- focusPD(pd, select, regex = regex)
        else
            pd <- subsetPD(pd, select, regex = regex)
    if (! missing(omit))
        pd <- subsetPD(pd, omit = omit, regex = regex)
    if (normalize)
        pd$total <- sum(pd$counts)
    pd
}

subsetIDX <- function(idx, pd, regex) {
    isIn <- if (regex) patMatchAny else function(x, s) x %in% s
    if (is.character(idx))
        which(sapply(pd$stacks, function(s) any(isIn(idx, s))))
    else if (is.logical(idx))
        which(idx)
    else 
        idx
}

patMatchAny <- function(p, x) {
    v <- grepl(p[1], x)
    for (q in p[-1])
        v <- v | grepl(q, x)
    v
}

subsetPD <- function(pd, select, omit, regex = TRUE) {
    stackIDX <- seq_along(pd$stacks)
    if (missing(select))
        select <- stackIDX
    keep <- subsetIDX(select, pd, regex)
    if (! missing(omit))
        keep <- setdiff(keep, subsetIDX(omit, pd, regex))
    
    pd$stacks <- pd$stacks[keep]
    pd$refs <- pd$refs[keep]
    pd$counts <- pd$counts[keep]
    pd$gccounts <- pd$gccounts[keep]

    traceKeep <- which(pd$trace %in% keep)
    pd$inGC <- pd$inGC[traceKeep]
    pd$trace <- pd$trace[traceKeep]

    map <- match(stackIDX, keep)
    pd$trace <- map[pd$trace]

    pd
}

focusPD <- function(pd, which, self.pct = 0, total.pct = 0, regex = TRUE) {
    if (self.pct > 0 || total.pct > 0) {
        fs <- funSummary(pd, srclines = FALSE)
        funs <- rownames(fs)
        if (self.pct > 0)
            pd <- focusPD(pd, funs[fs$self.pct >= self.pct], regex = FALSE)
        if (total.pct > 0) {
            high <- funs[fs$total.pct >= total.pct]
            pd <- focusPD(pd, high, regex = FALSE)
            pd <- trimTop(pd, setdiff(funs, high))
        }
    }
    if (missing(which))
        pd
    else
        skipPD(subsetPD(pd, which, regex = regex), which, regex = regex)
}
trimTop <- function(pd, funs) {
    to <- sapply(pd$stacks, function(s) {
        notLow <- ! s %in% funs
        if (any(notLow)) max(which(notLow)) else 0
    })
    pd <- prunePD(pd, to = to)
    drop <- sapply(pd$stacks, function(s) all(s %in% funs))
    subsetPD(pd, omit = drop)
}

compactPD <- function(pd) {
    ## **** need SIMPLIFY=FALSE here since mapply creates a matrix if all
    ## **** elements of the result happen to be the same length. Might
    ## **** be more robust to use paste() rather than c() (probably
    ## **** what I intended originally).
    key <- mapply(c, pd$stacks, pd$refs, SIMPLIFY = FALSE)
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

skipIDX <- function(pd, what, regex) {
    if (regex)
        isIn <- function(s, x) sapply(s, function(s) any(patMatchAny(x, s)))
    else
        isIn <- function(s, x) s %in% x
    if (is.character(what)) {
        findFirst <- function(s) {
            idx <- isIn(s, what)
            if (any(idx))
                min(which(idx)) - 1
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

OtherFunsToken <- "<Other>"

skipPD <- function(pd, what, merge = FALSE, regex = TRUE) {
    idx <- skipIDX(pd, what, regex)

    skip <- function(stack, refs, i) {
        n <- idx[i]
        if (n > 0) {
            if (n < length(stack)) {
                skip <- 1 : n
                stack <- stack[-skip]
                refs <- refs[-skip]
            }
            else {
                stack <- OtherFunsToken
                refs <- c(NA_character_, NA_character_)
            }
        }
        else if (merge) {
            stack <- OtherFunsToken
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

prunePD <- function(pd, to, by, merge = FALSE) {
    if (missing(by))
        idx <- pruneIDX(pd, if (is.character(to)) to else -to)
    else
        idx <- pruneIDX(pd, by)

    prune <- function(stack, refs, i) {
        n <- idx[i]
        slen <- length(stack)

        if (n < 0)
            if (-n < slen)
                n <- slen + n
            else
                n <- 0

        if (n > 0) {
            if (n < slen) {
                drop <- (slen - n + 1) : slen
                stack <- stack[-drop]
                refs <- refs[-(drop + 1)]
            }
            else {
                stack <- OtherFunsToken
                refs <- c(NA_character_, NA_character_)
            }
        }
        list(stack = stack, refs = refs)
    }

    transformPD(pd, prune)
}

mergeGC <- function(pd) {
    ostacks <- pd$stacks
    hasGC <- which(pd$gccounts > 0)

    pd$stacks <- c(pd$stacks, lapply(pd$stacks[hasGC], c, "<GC>"))
    pd$refs <- c(pd$refs, lapply(pd$refs[hasGC], c, NA_character_))
    pd$counts <- c(pd$counts - pd$gccounts, pd$gccounts[hasGC])
    pd$gccounts <- rep(0, length(pd$counts))

    map <- match(seq_along(ostacks), hasGC) + length(ostacks)
    pd$trace[pd$inGC] <- map[pd$trace[pd$inGC]]

    pd$inGC[] <- FALSE
    pd$haveGC <- FALSE

    pd
}


###
### Hot path summaries
###

pathAbbrev <- function(paths, short) {
    pad <- function(n) paste(rep(short, n), collapse = "")
    sapply(strsplit(paths, " -> "),
           function(path) {
               n <- length(path)
               label <- if (n > 1) paste0(pad(n - 1), path[n]) else path
           })
}

stripRefs <- function(pd) {
    pd$refs <- lapply(pd$refs, function(r) rep(NA_character_, length(r)))
    pd <- compactPD(pd)
    pd$haveRefs <- FALSE
    pd
}

percent <- function(x, y) round(100 * x / y, 2)

## The Hot Path order orders each level according to the number of
## hits within the call chain upt to that point.
hotPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1)) {
        key <- sapply(stacks[ord], `[`, i)
        tbl <- aggregate(list(val = -counts[ord]), list(key = key), sum)
        val <- tbl$val[match(key, tbl$key)]
        ord <- ord[order(val, na.last = TRUE)]
    }
    ord
}

UnknownFunToken <- "??"

hotPathData <- function(pd) {
    files <- pd$files
    pathLabels <- function(s, t) {
        n <- length(s)
        if (n == 0)
            funLabels(UnknownFunToken, t[1], files)
        else if (is.na(t[n + 1]))
            funLabels(s, t[1:n], files)
        else
            funLabels(c(s, UnknownFunToken), t, files)
    }

    pl <- mapply(pathLabels, pd$stacks, pd$refs, SIMPLIFY = FALSE)
    ord <- hotPathOrd(pl, pd$counts)
    stacks <- pl[ord]
    counts <- pd$counts[ord]
    gccounts <- pd$gccounts[ord]

    pathData <- function(k) {
        s <- stacks[[k]]
        keys <- sapply(1 : length(s),
                       function(i) paste(s[1 : i], collapse = " -> "))
        data.frame(key = keys,
                   count = counts[k],
                   gccount = gccounts[k])
    }
    tbl <- do.call(rbind, lapply(1:length(stacks), pathData))

    ## **** turn off stringsAsFactore in aggregate?
    data <- aggregate(list(count = tbl$count, gccount = tbl$gccount),
                      list(key = tbl$key),
                      sum)
    data$key <- as.character(data$key)

    selfidx <- match(data$key, sapply(stacks, paste, collapse = " -> "))
    data$self <- ifelse(is.na(selfidx), 0, counts[selfidx])
    data$gcself <- ifelse(is.na(selfidx), 0, gccounts[selfidx])
    data
}

hotPathsPct <- function(data, self, gc, grandTotal) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.pct = percent(data$count, grandTotal),
                      stringsAsFactors = FALSE)
    if (gc) val$gc.pct = percent(data$gccount, grandTotal)
    if (self) {
        val$self.pct = percent(data$self, grandTotal)
        if (gc) val$gcself.pct = percent(data$gcself, grandTotal)
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPathsHits <- function(data, self, gc) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.hits = data$count,
                      stringsAsFactors = FALSE)
    if (gc) val$gc.hits = data$gccount
    if (self) {
        val$self.hits = data$self
        if (gc) val$gcself.hits = data$gcself
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPathsTime <- function(data, self, gc, delta) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.time = data$count * delta,
                      stringsAsFactors = FALSE)
    if (gc) val$gc.time = data$gccount * delta
    if (self) {
        val$self.time = data$self * delta
        if (gc) val$gcself.time = data$gcself * delta
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPaths <- function(pd, value = c("pct", "time", "hits"),
                     self = TRUE, srclines = FALSE, GC = FALSE,
                     maxdepth = NA, self.pct = 0, total.pct = 0,
                     short = ". ", nlines = 20) {
    value <- match.arg(value)
    if (! is.na(maxdepth))
        pd <- prunePD(pd, maxdepth)

    if (! srclines && pd$haveRefs) pd <- stripRefs(pd)

    data <- hotPathData(pd)

    if (! is.na(nlines) && length(data$count) > nlines) {
        top <- head(order(data$count, decreasing = TRUE), nlines)
        data <- data[1:nrow(data) %in% top,]
    }

    if (self.pct > 0)
        data <- data[data$self >= pd$total * (self.pct / 100),]
    if (total.pct > 0)
        data <- data[data$count >= pd$total * (total.pct / 100),]

    data$pa <- pathAbbrev(data$key, short)

    if (value == "pct")
        hotPathsPct(data, self, GC && pd$haveGC, pd$total)
    else if (value == "time")
        hotPathsTime(data, self, GC && pd$haveGC, pd$interval / 1.0e6)
    else
        hotPathsHits(data, self, GC && pd$haveGC)
}

print.proftools_hotPaths <- function(x, ..., right = FALSE, row.names = FALSE)
    print.data.frame(x, ..., right = right, row.names = row.names)


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

lineRefs <- function(line, refs, useSite)
    unique(cbind(fun = "", refs = refs))

## leafRef <- function(line, refs, useSite) {
##     n <- length(line)
##     cbind(fun = "", refs = refs[n + 1])
## }
leafRef <- function(line, refs, useSite)
    cbind(fun = "", refs = NA_character_)

refCounts <- function(pd) {
    val <- entryCounts(pd, lineRefs, leafRef, TRUE)
    val$fun <- val$self <- val$gcself <- NULL
    val[! is.na(val$refs), ]
}

funSummaryPct <- function(fc, label, gc, grandTotal) {
    pct <- percent(fc$total, grandTotal)
    spct <- percent(fc$self, grandTotal)
    if (gc) {
        gcpct <- percent(fc$gctotal, grandTotal)
        sgcpct <- percent(fc$gcself, grandTotal)
        data.frame(total.pct = pct, gc.pct = gcpct,
                   self.pct = spct, gcself.pct = sgcpct,
                   row.names = label)
    }
    else
        data.frame(total.pct = pct, self.pct = spct, row.names = label)
}

funSummaryTime <- function(fc, label, gc, delta) {
    tm <- fc$total * delta
    stm <- fc$self * delta
    if (gc) {
        gctm <- fc$gctotal * delta
        sgctm <- fc$gcself * delta
        data.frame(total.time = tm, gc.time = gctm,
                   self.time = stm, gcself.time = sgctm,
                   row.names = label)
    }
    else
        data.frame(total.time = tm, self.time = stm, row.names = label)
}

funSummaryHits <- function(fc, label, gc) {
    if (gc)
        data.frame(total.hits = fc$total, gc.hits = fc$gctotal,
                   self.hits = fc$self, gcself.hits = fc$gcself,
                       row.names = label)
    else
        data.frame(total.hits = fc$total, self.hits = fc$self,
                   row.names = label)
}

## Extract the file indices and line numbers from source references of
## the form FN#LN.
refFN <- function(refs)
    as.integer(sub("([[:digit:]]+)#[[:digit:]]+", "\\1", refs))

refLN <- function(refs)
    as.integer(sub("[[:digit:]]+#([[:digit:]]+)", "\\1", refs))

funLabels <- function(fun, site, files) {
    if (all(is.na(site)))
        fun
    else {
        file <- basename(files[refFN(site)])
        line <- refLN(site)
        funsite <- sprintf("%s (%s:%d)", fun, file, line)
        ifelse(is.na(site), fun, funsite)
    }
}

funSummary <- function(pd, byTotal = TRUE,
                       value = c("pct", "time", "hits"),
                       srclines = TRUE,
                       GC = TRUE, self.pct = 0, total.pct = 0) {
    value <- match.arg(value)

    fc <- funCounts(pd, srclines)
    if (byTotal)
        fc <- fc[rev(order(fc$total)), ]
    else
        fc <- fc[rev(order(fc$self)), ]

    if (self.pct > 0)
        fc <- fc[fc$self >= pd$total * (self.pct / 100),]
    if (total.pct > 0)
        fc <- fc[fc$total >= pd$total * (total.pct / 100),]

    label <- funLabels(fc$fun, fc$site, pd$files)

    if (value == "pct")
        funSummaryPct(fc, label, GC && pd$haveGC, pd$total)
    else if (value == "time")
        funSummaryTime(fc, label, GC && pd$haveGC, pd$interval / 1.0e6)
    else
        funSummaryHits(fc, label, GC && pd$haveGC)
}

callSummary <- function(pd, byTotal = TRUE,
                        value = c("pct", "time", "hits"),
                        srclines = TRUE,
                        GC = TRUE, self.pct = 0, total.pct = 0) {
    value <- match.arg(value)

    cc <- callCounts(pd, srclines, srclines)
    if (byTotal)
        cc <- cc[rev(order(cc$total)), ]
    else
        cc <- cc[rev(order(cc$self)), ]

    if (self.pct > 0)
        cc <- cc[cc$self >= pd$total * (self.pct / 100),]
    if (total.pct > 0)
        cc <- cc[cc$total >= pd$total * (total.pct / 100),]

    caller.label <- funLabels(cc$caller, cc$caller.site, pd$files)
    callee.label <- funLabels(cc$callee, cc$callee.site, pd$files)
    label <- paste(caller.label, callee.label, sep = " -> ")

    if (value == "pct")
        funSummaryPct(cc, label, GC && pd$haveGC, pd$total)
    else if (value == "time")
        funSummaryTime(cc, label, GC && pd$haveGC, pd$interval / 1.0e6)
    else
        funSummaryHits(cc, label, GC && pd$haveGC)
}


###
### Flame graph and time graph
###

## For 'standard' flame graph order the stacks so they are
## alphabetical within lines within calls, with missing entires
## first. This does a lexicographic sort by sorting on the top entry
## first, then the next, and do on; since the sorts are stable this
## keeps the top levels sorted within the lower ones.
alphaPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1))
        ord <- ord[order(sapply(stacks[ord], `[`, i), na.last = FALSE)]
    ord
}

## The next two functions compute the data to be used for drawing as a
## data frame with columns left, bottom, right, top, col, and label.
fgDataLine <- function(k, stacks, counts) {
    runs <- rle(sapply(stacks, `[`, k))
    lens <- runs$lengths
    n <- length(lens)
    csums <- cumsum(tapply(counts, rep(1 : n, lens), sum))

    top <- k
    bottom <- k - 1
    left <- c(0, csums[-n])
    right <- csums
    cols <- rgb(runif(n), runif(n), runif(n))
    label <- runs$values

    show <- ! is.na(label)
    nshow <- sum(show)
    data.frame(left = left[show], bottom = rep(bottom, nshow),
               right = right[show], top = rep(top, nshow),
               col = cols[show], label = label[show],
               stringsAsFactors = FALSE)
}

fgData <- function(stacks, counts, reorder = c("alpha", "hot", "no"),
                   colormap) {
    mx <- max(sapply(stacks, length))

    reorder <- match.arg(reorder)
    if (reorder != "no") {
        ## ord <- seq_along(stacks)
        ## for (i in (mx : 1))
        ##     ord <- ord[order(sapply(stacks[ord], `[`, i), na.last = FALSE)]
        if (reorder == "hot")
            ord <- hotPathOrd(stacks, counts)
        else
            ord <- alphaPathOrd(stacks, counts)
        stacks <- stacks[ord]
        counts <- counts[ord]
    }

    val <- do.call(rbind, lapply(1 : mx, fgDataLine, stacks, counts))
    if (! is.null(colormap))
        val$col <- colormap(val$label)
    val
}

## This is computes the data with fdData and draws the graph with base
## graphics. This is the bit we would need to change for grid, maybe
## ggplot2, and for svg output.
flameGraph <- function(stacks, counts, reorder, colormap) {
    fdg <- fgData(stacks, counts, reorder, colormap)
    left <- fdg$left
    bottom <- fdg$bottom
    right <- fdg$right
    top <- fdg$top
    col <- fdg$col
    label <- fdg$label

    plot(c(min(left), max(right)), c(min(bottom), max(top)),
         type = "n", axes = FALSE, xlab = "", ylab = "")

    rect(left, bottom, right, top, col = col)
    
    ## half-em for half-character offset used by text() when pos
    ## argument is used.
    hm <- 0.5 * strwidth("m")

    show <- (strheight(label) <= 0.9 &
             strwidth(label) + 2 * hm <= 0.8 * (right - left))
    if (any(show))
        text(left[show], bottom[show] + 0.4, label[show], pos = 4)
}

htmlencode <- function(x)
    sub(">", "&gt;", sub("<", "&lt;", x))

svgFlameGraph <- function(file, stacks, counts, reorder, colormap) {
    fdg <- fgData(stacks, counts, reorder, colormap)
    mx <- max(fdg$top)
    totalCount <- max(fdg$right)
    counts <- fdg$right-fdg$left
    percents <- round(counts*100/totalCount, 2)
    widths <- round(percents*1180/100, 2)
    y <- 33 + (mx-fdg$top)*16
    x <- 10+round(fdg$left*1180/totalCount, 2)
    col <- fdg$col
    labels <- htmlencode(fdg$label)
    
    svgCode = paste("<rect x=\"", x, "\" y=\"", y, 
    "\" width=\"", widths, "\" height=\"15.0\" fill=\"", col, 
    "\" rx=\"2\" ry=\"2\" onmouseover=\"s('", labels, " (",
    counts, " samples, ", percents, "%)')\" onmouseout=\"c()\" />", 
    sep="")

    show <- (! is.na(labels) & 10*nchar(labels)<widths)
    if (any(show))
        svgCode = append(svgCode, paste("<text text-anchor=\"\" x=\"", 
        x[show]+3, "\" y=\"", y[show]+10.5, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\" onmouseover=\"s('", 
        labels[show], " (", counts[show], " samples, ", percents[show], 
        "%)')\" onmouseout=\"c()\" >", labels[show],"</text>", sep=""))
        
    writeFile(file, svgCode, mx)
}
## This writes the header of the svg file
writeFile <- function(file, svgCode, mx){
    write(c(paste("<?xml version=\"1.0\" standalone=\"no\"?>
    <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">
    <svg version=\"1.1\" width=\"1200\" height=\"", 16*mx+66, "\" onload=\"init(evt)\" viewBox=\"0 0 1200 ", 16*mx+66, "\" xmlns=\"http://www.w3.org/2000/svg\" >
    <defs >
    <linearGradient id=\"background\" y1=\"0\" y2=\"1\" x1=\"0\" x2=\"0\" >
        <stop stop-color=\"#eeeeee\" offset=\"5%\" />
        <stop stop-color=\"#eeeeb0\" offset=\"95%\" />
    </linearGradient>
    </defs>
    <style type=\"text/css\">
    rect[rx]:hover { stroke:black; stroke-width:1; }
    text:hover { stroke:black; stroke-width:1; stroke-opacity:0.35; }
    </style>
    <script type=\"text/ecmascript\">
    <![CDATA[
    var details;
    function init(evt) { details = document.getElementById(\"details\").firstChild; }
    function s(info) { details.nodeValue = info; }
    function c() { details.nodeValue = ' '; }
    ]]>
    </script>
    <rect x=\"0.0\" y=\"0\" width=\"1200.0\" height=\"", 16*mx+66, "\" fill=\"url(#background)\"  />
    <text text-anchor=\"middle\" x=\"600\" y=\"24\" font-size=\"17\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\"  >Call Graph</text>
    <text text-anchor=\"left\" x=\"10\" y=\"", 16*mx+50, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\"  >Function:</text>
    <text text-anchor=\"\" x=\"70\" y=\"", 16*mx+50, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\" id=\"details\" > </text>", sep=""), svgCode, "</svg>"), file = file)
}

## Merge non-NA refs into stacks. Add <Internal> for lef refs.
refStacks <- function(d) {
    rs <-function(s, r, d) {
        if (is.na(r[length(r)]))
            r <- r[-length(r)]
        else
            s <- c(s, "<Internal>")
        fl <- d$files[refFN(r)]
        ln <- refLN(r)
        paste0(s, ifelse(is.na(r), "", sprintf(" (%s:%d)", fl, ln)))
    }
    lapply(1:length(d$stacks), function(i) rs(d$stacks[[i]], d$refs[[i]], d))
}

## produce a flame graph from an Rprof file
fg <- function(file, svgfile, reorder = c("hot", "alpha", "no"),
               colormap = NULL, srclines = FALSE) {
    reorder <- match.arg(reorder)
    if (is.character(file))
        d <- readPD(file)
    else
        d <- file
    counts <- d$counts
    stacks <- if (srclines) refStacks(d) else d$stacks
    if (! missing(svgfile))
        svgFlameGraph(svgfile, stacks, counts, reorder, colormap)
    else
        flameGraph(stacks, counts, reorder, colormap)
}

## produce a time graph (like profr) from an Rprof file
tg <- function(file, svgfile, colormap = NULL, srclines = FALSE) {
    if (is.character(file))
        d <- readPD(file)
    else
        d <- file
    r <- rle(d$trace)
    stacks <- if (srclines) refStacks(d) else d$stacks
    tstacks <- stacks[r$values]
    counts <- r$lengths
    if (! missing(svgfile))
        svgFlameGraph(svgfile, tstacks, counts, "no", colormap)
    else
        flameGraph(tstacks, counts, "no", colormap)
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
## defined.  Otherwise treat the caller's home file as unkown. The
## result returned by this function is a named vector with one element
## per function for which the home file is assumed know.  The names
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
## To fc we add the indes of the home file for each cunction (NA if
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
    fc <- funCounts(pd, FALSE)
    fc$total <- fc$gctotal <- NULL
    f <- sapply(pd$stacks, function(x) x[length(x)])
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

writeCallgrindFile <- function(pd, file = "Rprof.cg", GC = TRUE)
    writeCG(file, pd, GC)


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

formatTrace <- function(trace, maxlen = 50, skip = 0, trimtop = FALSE) {
    if (skip > 0)
        trace <- trace[-(1 : skip)]
    out <- paste(trace, collapse = " -> ")
    if (trimtop)
        while (nchar(out) > maxlen && length(trace) > 1) {
            trace <- trace[-length(trace)]
            out <- paste(paste(trace, collapse = " -> "), "... ")
        }
    else
        while (nchar(out) > maxlen && length(trace) > 1) {
            trace <- trace[-1]
            out <- paste("... ", paste(trace, collapse = " -> "))
        }
    out
}

printPaths <- function(pd, n, ...) {
    ord = rev(order(pd$counts))
    if (! missing(n) && length(ord) > n)
        ord <- ord[1 : n]
    tot <- pd$total
    pct <- percent(pd$counts[ord], tot)
    gcpct <- percent(pd$gccounts[ord], tot)
    paths <- sapply(pd$stacks[ord], formatTrace, ...)
    mapply(function(x, y, z) cat(sprintf("%5.1f %5.1f   %s\n", x, y, z)),
           pct, gcpct, paths, SIMPLIFY = FALSE)
    invisible(NULL)
}

pathSummaryPct <- function(apd, gc, tot) {
    counts <- apd$counts
    gccounts <- apd$gccounts
    paths <- apd$paths

    pct <- percent(counts, tot)
    if (gc) {
        gcpct <- percent(gccounts, tot)
        data.frame(total.pct = pct, gc.pct = gcpct, row.names = paths)
    }
    else
        data.frame(total.pct = pct, row.names = paths)
}
                           
pathSummaryTime <- function(apd, gc, delta) {
    counts <- apd$counts
    gccounts <- apd$gccounts
    paths <- apd$paths
    tm <- counts * delta
    if (gc) {
        gctm <- gccounts * delta
        data.frame(total.time = tm, gc.time = gctm, row.names = paths)
    }
    else
        data.frame(total.time = tm, row.names = paths)
}

pathSummaryHits <- function(apd, gc) {
    hits <- apd$counts
    gchits <- apd$gccounts
    paths <- apd$paths
    if (gc)
        data.frame(total.hits = hits, gc.hits = gchits, row.names = paths)
    else
        data.frame(total.hits = hits, row.names = paths)
}

pathSummary <- function(pd, value = c("pct", "time", "hits"),
                        srclines = FALSE, GC = TRUE, total.pct = 0, ...) {
    value <- match.arg(value)

    if (srclines && pd$haveRefs) {
        files <- pd$files ## shorter: as.character(seq_along(pd$files)
        rstacks <- mapply(function(a, b) funLabels(a, b, files),
                          pd$stacks,
                          sapply(pd$refs, function(x) x[-length(x)]),
                          SIMPLIFY = FALSE)
        paths <- sapply(rstacks, formatTrace, ...)
    }
    else
        paths <- sapply(pd$stacks, formatTrace, ...)

    ## need to aggregate in case some collapsed paths are identical
    ## or some paths differ only in source references.
    apd <- aggregateCounts(list(paths = paths),
                           cbind(counts = pd$counts, gccounts = pd$gccounts))
    apd <- apd[rev(order(apd$counts)),]

    if (total.pct > 0)
        apd <- apd[apd$counts >= pd$total * (total.pct / 100),]

    if (value == "pct")
        pathSummaryPct(apd, GC && pd$haveGC, pd$total)
    else if (value == "time")
        pathSummaryTime(apd, GC && pd$haveGC, pd$interval / 1.0e6)
    else
        pathSummaryHits(apd, GC && pd$haveGC)
}

srcSummary <- function(pd, byTotal = TRUE,
                       value = c("pct", "time", "hits"),
                       GC = TRUE, total.pct = 0) {
    if (! pd$haveRefs)
        stop("profile data does not contain source information")

    value <- match.arg(value)

    rc <- refCounts(pd)

    file <- basename(pd$files[refFN(rc$refs)])
    line <- refLN(rc$refs)
    label <- paste(file, line, sep = ":")

    if (GC && pd$haveGC)
        val <- cbind(total = rc$total, gctotal = rc$gctotal)
    else
        val <- cbind(total = rc$total)
    rownames(val) <- label

    if (total.pct > 0)
        val <- val[val[,1] >= pd$total * (total.pct / 100),]

    if (value == "pct") {
        tot <- pd$total
        colnames(val) <- paste(colnames(val), "pct", sep = ".")
        as.data.frame(percent(val, tot))
    }
    else if (value == "time") {
        delta <- pd$interval / 1.0e6
        colnames(val) <- paste(colnames(val), "time", sep = ".")
        as.data.frame(val * delta)
    }
    else {
        colnames(val) <- paste(colnames(val), "hits", sep = ".")
        as.data.frame(val)
    }

}
