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

    list(interval = interval,
         haveGC = haveGC, haveRefs = haveRefs, haveMem  = haveMem)
}

readProfileData <- function(filename = "Rprof.out")
    readPD(filename)

readRStudioProfileCacheData <- function() {
    profDir <- options("profvis.prof_output")$profvis.prof_output
    if (is.null(profDir)) return(NULL)
    stackFiles <- list.files(profDir, pattern = "*.Rprof", full.names = TRUE)
    if (length(stackFiles) == 0) return(NULL)
    details <- file.info(stackFiles)
    newest <- stackFiles[order(details$mtime, decreasing = TRUE)[1]]
    readProfileData(newest)
}

readPDlines <- function(con, hdr) {
    stacks <-  readLines(con)
    if (hdr$haveRefs) {
        fstacks <- grepl("#File ", stacks)
        files <- sub("#File [[:digit:]]: (.+)", "\\1", stacks[fstacks])
        stacks <- stacks[! fstacks]
    }
    else files <- NULL

    if (hdr$haveMem && length(stacks) > 0) {
        memstuff <- sub(":([0-9]+:[0-9]+:[0-9]+:[0-9]+):.*", "\\1", stacks)
        memnums <- sapply(strsplit(memstuff, ":"), as.numeric)
        mem <- as.vector(t(memnums) %*%  c(8, 8, 1, 0)) / 1048576
        stacks <- substr(stacks, nchar(memstuff) + 3, nchar(stacks))
    }
    else mem <- NULL

    ## remove any lines with only a <GC> entry
    onlyGC <- grepl("^\"<GC>\" $", stacks)
    if (any(onlyGC)) {
        stacks <- stacks[! onlyGC]
        if (! is.null(mem))
            mem <- mem[! onlyGC]
    }

    ## record and strip out GC info
    inGC <- grepl("<GC>", stacks)
    stacks <- sub("\"<GC>\" ", "", stacks)
    
    ustacks <- unique(stacks)
    trace <- match(stacks, ustacks)

    list(files = files, stacks = ustacks, trace = trace, inGC = inGC, mem = mem)
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

countStacks <- function(trace, inGC, mem) {
    counts <- as.numeric(table(trace))
    gccounts <- as.numeric(tapply(inGC, trace, sum))
    total <- sum(counts)
    if (is.null(mem))
        alloc <- rep(NA_real_, length(counts))
    else
        alloc <- as.numeric(tapply(c(0, pmax(diff(mem), 0)), trace, sum))
    list(counts = counts, gccounts = gccounts, total = total, alloc = alloc)
}

readPD <- function(file) {
    con <- file(file, "r")
    on.exit(close(con))
    
    hdr <- readPDheader(con)
    data <- readPDlines(con, hdr)
    sdata <- splitStacks(data)
    counts <- countStacks(sdata$trace, sdata$inGC, sdata$mem)
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

filterProfileData <- function(pd, select, omit, focus,
                              skip = NA,
                              maxdepth = NA,
                              self.pct = 0, total.pct = 0,
                              interval,
                              normalize = FALSE,
                              regex = FALSE) {
    if (!is.na(skip))
        pd <- skipPD(pd, skip)
    if (! is.na(maxdepth))
        pd <- prunePD(pd, to = maxdepth)
    pd <- focusPD(pd, self.pct = self.pct, total.pct = total.pct)
    if (! missing(select))
        pd <- subsetPD(pd, select, regex = regex)
    if (! missing(focus))
        pd <- focusPD(pd, focus, regex = regex)
        else
    if (! missing(omit))
        pd <- subsetPD(pd, omit = omit, regex = regex)
    if (! missing(interval))
        pd <- timeSubsetPD(pd, interval)
    if (normalize)
        pd$total <- sum(pd$counts)
    pd
}
filterProfileData <- function(pd, ..., normalize = FALSE, regex = FALSE)
{
    fargs <- list(...)
    fnames <- names(fargs)
    if (is.null(fnames) || any(nchar(fnames) == 0))
        stop("all filters must be named")
    for (i in seq_along(fargs)) {
        pd <- switch(fnames[i],
                     skip = skipPD(pd, fargs[[i]]),
                     maxdepth = prunePD(pd, to = fargs[[i]]),
                     self.pct = focusPD(pd, self.pct = fargs[[i]]),
                     total.pct = focusPD(pd, total.pct = fargs[[i]]),
                     focus = focusPD(pd, fargs[[i]], regex = regex),
                     select = subsetPD(pd, fargs[[i]], regex = regex),
                     omit = subsetPD(pd, omit = fargs[[i]], regex = regex),
                     interval = timeSubsetPD(pd, fargs[[i]]),
                     merge.pct = mergeStacks(pd, fargs[[i]] / 100),
                     stop("unknown filter: ", fnames[i]))
    }
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
    pd$alloc <- pd$alloc[keep]

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
        funs <- fs$fun
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
    if (length(pd$stacks) == 0) return(pd)
    ## **** need SIMPLIFY=FALSE here since mapply creates a matrix if all
    ## **** elements of the result happen to be the same length. Might
    ## **** be more robust to use paste() rather than c() (probably
    ## **** what I intended originally).
    key <- mapply(c, pd$stacks, pd$refs, SIMPLIFY = FALSE)
    map <- match(key, unique(key))
    ct <- aggregateCounts(data.frame(key = map),
                          cbind(counts = pd$counts,
                          gccounts = pd$gccounts,
                          alloc = pd$alloc))
    ct <- ct[order(ct$key),] ## may not be needed
    invmap <- match(unique(key), key)
    pd$stacks <- pd$stacks[invmap]
    pd$refs <- pd$refs[invmap]
    pd$counts <- ct$counts
    pd$gccounts <- ct$gccounts
    pd$alloc <- ct$alloc
    pd$trace <- map[pd$trace]
    pd
}

checkStackRefs <- function(val, nf) {
    s <- val$stack
    r <- val$refs
    #if (length(s) == 0)
    #    stop("stacks must have at least one entry")
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

timeSubsetPD <- function(pd, interval) {
    lo <- interval[1]
    hi <- interval[2]
    idx <- lo : hi
    pd$trace <- pd$trace[idx]
    pd$inGC <- pd$inGC[idx]
    sidx <- unique(pd$trace)
    pd$trace <- match(seq_along(pd$stacks), sidx)[pd$trace]
    pd$stacks <- pd$stacks[sidx]
    pd$refs <- pd$refs[sidx]
    pd$mem <- pd$mem[idx]
    cs <- countStacks(pd$trace,pd$inGC,pd$mem)
    pd$counts <- cs$counts
    pd$gccounts <- cs$gccounts
    pd$alloc<-cs$alloc
    pd
}

mergeStacks <- function(pd, p){
    trim <- function(x) {
        n <- length(x)
        if (n > 0) x[-n] else x
    }
    repeat {
        lns <- sapply(pd$stacks, length)
        low <- pd$counts < p * pd$total
        if (sum(low) == 0 || max(lns[low]) <=1)
            return(pd)
        to <- ifelse(low & lns == max(lns[low]), lns - 1, lns)
        pd <- prunePD(pd, to = to)
    }
}

dropSourceFrames <- function(pd, addTOP = TRUE) {
    if (length(pd$stacks) > 0 &&
        all(sapply(pd$stacks, length) > 4) &&
        pd$haveRefs &&
        isTRUE(all(diff(refFN(sapply(pd$refs,`[`, 5))) == 0)) &&
        all(sapply(pd$stacks, `[`, 1) == "source") &&
        all(sapply(pd$stacks, `[`, 2) == "withVisible") &&
        all(sapply(pd$stacks, `[`, 3) == "eval") &&
        all(sapply(pd$stacks, `[`, 4) == "eval")) {
        pd <- filterProfileData(pd, skip = 4)
        if (addTOP) {
            tref <- sprintf("%d#1", refFN(pd$refs[[1]][1]))
            pd$refs <- lapply(pd$refs, function(x) c(tref, x))
            pd$stacks <- lapply(pd$stacks, function(x) c("<TOP>", x))
        }
    }
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
    alloc <- pd$alloc[ord]
    if (is.null(alloc)) ## an old pd maybe?
        alloc <- rep(NA_real_, length(counts))

    pathData <- function(k) {
        s <- stacks[[k]]
        keys <- sapply(1 : length(s),
                       function(i) paste(s[1 : i], collapse = " -> "))
        data.frame(key = keys,
                   count = counts[k],
                   gccount = gccounts[k],
                   alloc = alloc[k])
    }
    tbl <- do.call(rbind, lapply(seq_along(stacks), pathData))

    ## **** turn off stringsAsFactore in aggregate?
    data <- aggregate(list(count = tbl$count, gccount = tbl$gccount, alloc = tbl$alloc),
                      list(key = tbl$key),
                      sum)
    data$key <- as.character(data$key)

    selfidx <- match(data$key, sapply(stacks, paste, collapse = " -> "))
    data$self <- ifelse(is.na(selfidx), 0, counts[selfidx])
    data$gcself <- ifelse(is.na(selfidx), 0, gccounts[selfidx])
    data$allocself <- ifelse(is.na(selfidx), 0, alloc[selfidx])
    data
}

hotPathsPct <- function(data, self, gc, memory, grandTotal) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.pct = percent(data$count, grandTotal),
                      stringsAsFactors = FALSE)
    if (gc) val$gc.pct = percent(data$gccount, grandTotal)
    if (memory) val$alloc <- data$alloc
    if (self) {
        val$self.pct = percent(data$self, grandTotal)
        if (gc) val$gcself.pct = percent(data$gcself, grandTotal)
        if (memory) val$allocself <- data$allocself
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPathsHits <- function(data, self, gc, memory) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.hits = data$count,
                      stringsAsFactors = FALSE)
    if (gc) val$gc.hits = data$gccount
    if (memory) val$alloc <- data$alloc
    if (self) {
        val$self.hits = data$self
        if (gc) val$gcself.hits = data$gcself
        if (memory) val$allocself <- data$allocself
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPathsTime <- function(data, self, gc, memory, delta) {
    pa <- data$pa
    val <- data.frame(path = sprintf(sprintf("%%-%ds", max(nchar(pa))), pa),
                      total.time = data$count * delta,
                      stringsAsFactors = FALSE)
    if (gc) val$gc.time = data$gccount * delta
    if (memory) val$alloc <- data$alloc
    if (self) {
        val$self.time = data$self * delta
        if (gc) val$gcself.time = data$gcself * delta
        if (memory) val$allocself <- data$allocself
    }
    class(val) <- c("proftools_hotPaths", "data.frame")
    val
}

hotPaths <- function(pd, value = c("pct", "time", "hits"),
                     self = TRUE, srclines = TRUE, GC = FALSE, memory = FALSE,
                     maxdepth = 10, self.pct = 0, total.pct = 0,
                     short = ". ", nlines = NA) {
    if (length(pd$stacks) == 0) return(NULL)
    value <- match.arg(value)
    if (! is.na(maxdepth))
        pd <- prunePD(pd, maxdepth)
    GC <- GC && pd$haveGC
    memory <- memory && pd$haveMem

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
        hotPathsPct(data, self, GC, memory, pd$total)
    else if (value == "time")
        hotPathsTime(data, self, GC, memory, pd$interval / 1.0e6)
    else
        hotPathsHits(data, self, GC, memory)
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
    if (nrow(dcounts) == 0) return(dcounts)     
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
    val$allocself[is.na(val$allocself)] <- 0
    val
}

entryCounts0 <- function(pd, fun, control, names) {
    stacks <- pd$stacks
    refs <- pd$refs
    counts <- pd$counts
    gccounts <- pd$gccounts
    alloc <- pd$alloc
    if (is.null(alloc)) ## an old pd maybe?
        alloc <- rep(NA_real_, length(counts))

    which <- seq_along(stacks)
    
    doLine <- function(i) fun(stacks[[i]], refs[[i]], control)
    entries <- lapply(which, doLine)
    edf <- rbindEntries(entries)

    reps <- unlist(lapply(entries, nrow))
    tot <- rep(counts, reps)
    gctot <- rep(gccounts, reps)
    altot <- rep(alloc, reps)
    ct <- cbind(tot, gctot, altot, deparse.level = 0)
    colnames(ct) <- names
    aggregateCounts(edf, ct)
}

entryCounts <- function(pd, lineFun, leafFun, control) {
    aedf <- entryCounts0(pd, lineFun, control,
                         c("total", "gctotal", "alloc"))
    aledf <- entryCounts0(pd, leafFun, control,
                          c("self", "gcself", "allocself"))
    mergeCounts(aedf, aledf)
}

lineFuns <- function(line, refs, useSite) {
    n <- length(line)
    if (useSite)
        site <- refs[-(n + 1)]
    else site <- rep(NA_character_, n)
    unique(cbind(fun = line, site))
}

leafFun <- function(line, refs, useSite) {
    n <- length(line)
    fun <- line[n]
    site <- if (useSite) refs[n] else rep(NA_character_, n)
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
    val$fun <- val$self <- val$gcself <- val$allocself <- NULL
    val[! is.na(val$refs), ]
}

funSummaryPct <- function(fc, locs, gc, memory, grandTotal, call = FALSE) {
    total <- data.frame(total.pct = percent(fc$total, grandTotal))
    self <- data.frame(self.pct = percent(fc$self, grandTotal))
    if (gc) {
        total$gc.pct = percent(fc$gctotal, grandTotal)
        self$gcself.pct = percent(fc$gcself, grandTotal)
    }
    if (memory) {
        total$alloc <- fc$alloc
        self$allocself <- fc$allocself
    }
    cls <- if (call) "proftools_callSummary" else "proftools_funSummary"
    structure(data.frame(locs, total, self, stringsAsFactors = FALSE),
              class = c(cls, "data.frame"))
}

funSummaryTime <- function(fc, locs, gc, memory, delta, call = FALSE) {
    total <- data.frame(total.time = fc$total * delta)
    self <- data.frame(self.time = fc$self * delta)
    if (gc) {
        total$gc.time <- fc$gctotal * delta
        self$gcself.time <- fc$gcself * delta
    }
    if (memory) {
        total$alloc <- fc$alloc
        self$allocself <- fc$allocself
    }
    cls <- if (call) "proftools_callSummary" else "proftools_funSummary"
    structure(data.frame(locs, total, self, stringsAsFactors = FALSE),
              class = c(cls, "data.frame"))
}

funSummaryHits <- function(fc, locs, gc, memory, call = FALSE) {
    total <- data.frame(total.hits = fc$total)
    self <- data.frame(self.hits = fc$self)
    if (gc) {
        total$gc.hits <- fc$gctotal
        self$gcself.hits <- fc$gcself
    }
    if (memory) {
        total$alloc <- fc$alloc
        self$allocself <- fc$allocself
    }
    cls <- if (call) "proftools_callSummary" else "proftools_funSummary"
    structure(data.frame(locs, total, self, stringsAsFactors = FALSE),
              class = c(cls, "data.frame"))
}

format.proftools_funSummary <- function (x, ...) {
    file <- x$file
    line <- x$line
    x$file <- x$line <- NULL
    v <- format(as.data.frame(x), digits = 2, justify = "left")
    if (! is.null(file) && any(! is.na(file))) {
        fun <- x$fun
        funsite <- sprintf("%s (%s:%d)", fun, file, line)
        v$fun <- ifelse(is.na(file), fun, funsite)
        names(v)[1] <- "fun (file:line)"
    }
    v
}

print.proftools_funSummary <- function(x, ...) {
    x <- format(x)
    print.data.frame(x[-1], ..., row.names = x[[1]])
}

format.proftools_callSummary <- function(x, ...) {
    fun1 <- x$caller
    fun2 <- x$callee
    file1 <- x$caller.file
    line1 <- x$caller.line
    file2 <- x$callee.file
    line2 <- x$callee.line
    x$caller <- NULL
    x$caller.file <- x$caller.line <- NULL
    x$callee.file <- x$callee.line <- NULL
    v <- format(as.data.frame(x), digits = 2, justify = "left")
    if (nrow(v) == 0) return(v)
    if ((! is.null(file1) && any(! is.na(file1))) ||
        (! is.null(file2) && any(! is.na(file2)))) {
        funsite1 <- sprintf("%s (%s:%d)", fun1, file1, line1)
        funsite1 <- ifelse(is.na(file1), fun1, funsite1)
        funsite2 <- sprintf("%s (%s:%d)", fun2, file2, line2)
        funsite2 <- ifelse(is.na(file2), fun2, funsite2)
        call <- "caller (file:line) -> callee (file:line)"
    }
    else {
        funsite1 <- fun1
        funsite2 <- fun2
        call <- "caller -> callee"
    }
    v[[1]] <- paste(funsite1, "->", funsite2)
    names(v)[1] <- call
    v
}

print.proftools_callSummary <- function(x, ...) {
    x <- format(x)
    print.data.frame(x[-1], row.names = x[[1]])
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
funLocs <- function(fun, site, files) {
    if (length(site) == 0 || all(is.na(site)))
        data.frame(fun, stringsAsFactors = FALSE)
    else {
        file <- basename(files[refFN(site)])
        line <- refLN(site)
        data.frame(fun, file, line, stringsAsFactors = FALSE)
    }
}
callLocs <- function(caller, caller.site, callee, callee.site, files) {
    if (is.null(files))
        data.frame(caller, callee, stringsAsFactors = FALSE)
    else {
        caller.file <- basename(files[refFN(caller.site)])
        caller.line <- refLN(caller.site)
        callee.file <- basename(files[refFN(callee.site)])
        callee.line <- refLN(callee.site)
        data.frame(caller, caller.file, caller.line,
                   callee, callee.file, callee.line, stringsAsFactors = FALSE)
    }
}

funSummary <- function(pd, byTotal = TRUE,
                       value = c("pct", "time", "hits"),
                       srclines = TRUE, GC = TRUE, memory = FALSE,
                       self.pct = 0, total.pct = 0) {
    value <- match.arg(value)
    GC <- GC && pd$haveGC
    memory <- memory && pd$haveMem

    fc <- funCounts(pd, srclines)
    if (byTotal)
        fc <- fc[rev(order(fc$total)), ]
    else
        fc <- fc[rev(order(fc$self)), ]

    if (self.pct > 0)
        fc <- fc[fc$self >= pd$total * (self.pct / 100),]
    if (total.pct > 0)
        fc <- fc[fc$total >= pd$total * (total.pct / 100),]

    locs <- funLocs(fc$fun, fc$site, pd$files)

    if (value == "pct")
        funSummaryPct(fc, locs, GC, memory, pd$total)
    else if (value == "time")
        funSummaryTime(fc, locs, GC, memory, pd$interval / 1.0e6)
    else
        funSummaryHits(fc, locs, GC, memory)
}

callSummary <- function(pd, byTotal = TRUE,
                        value = c("pct", "time", "hits"),
                        srclines = TRUE, GC = TRUE, memory = FALSE,
                        self.pct = 0, total.pct = 0) {
    value <- match.arg(value)
    GC <- GC && pd$haveGC
    memory <- memory && pd$haveMem

    cc <- callCounts(pd, srclines, srclines)
    if (byTotal)
        cc <- cc[rev(order(cc$total)), ]
    else
        cc <- cc[rev(order(cc$self)), ]

    if (self.pct > 0)
        cc <- cc[cc$self >= pd$total * (self.pct / 100),]
    if (total.pct > 0)
        cc <- cc[cc$total >= pd$total * (total.pct / 100),]

    locs <- callLocs(cc$caller, cc$caller.site, cc$callee, cc$callee.site,
                     pd$files)

    if (value == "pct")
        funSummaryPct(cc, locs, GC, memory, pd$total, TRUE)
    else if (value == "time")
        funSummaryTime(cc, locs, GC, memory, pd$interval / 1.0e6, TRUE)
    else
        funSummaryHits(cc, locs, GC, memory, TRUE)
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
### Collecting profile data
###

profileExpr <- function(expr, GC = TRUE, srclines = TRUE, memory = TRUE) {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    Rprof(tmp, gc.profiling = GC, line.profiling = srclines,
          memory.profiling = memory)
    expr
    Rprof(NULL)
    pd <- readProfileData(tmp)
    mydepth <- length(sys.calls())
    skipPD(pd, mydepth)
}

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

pathSummaryPct <- function(apd, gc, memory, tot) {
    total <- data.frame(total.pct = percent(apd$counts, tot))
    if (gc)
        total$gc.pct <- percent(apd$gccounts, tot)
    if (memory)
        total$alloc <- apd$alloc
    data.frame(total, row.names = apd$paths)
}
                           
pathSummaryTime <- function(apd, gc, memory, delta) {
    total <- data.frame(total.time = apd$counts * delta)
    if (gc)
        total$gc.time <- apd$gccounts * delta
    if (memory)
        total$alloc <- apd$alloc
    data.frame(total, row.names = apd$paths)
}

pathSummaryHits <- function(apd, gc, memory) {
    total <- data.frame(total.hits = apd$counts)
    if (gc)
        total$gc.hits = apd$gccounts
    if (memory)
        total$alloc <- apd$alloc
    data.frame(total, row.names = apd$paths)
}

pathSummary <- function(pd, value = c("pct", "time", "hits"),
                        srclines = FALSE, GC = TRUE, memory = FALSE,
                        total.pct = 0, ...) {
    value <- match.arg(value)
    GC <- GC && pd$haveGC
    memory <- memory && pd$haveMem

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
                           cbind(counts = pd$counts,
                                 gccounts = pd$gccounts,
                                 alloc = pd$alloc))
    apd <- apd[rev(order(apd$counts)),]

    if (total.pct > 0)
        apd <- apd[apd$counts >= pd$total * (total.pct / 100),]

    if (value == "pct")
        pathSummaryPct(apd, GC, memory, pd$total)
    else if (value == "time")
        pathSummaryTime(apd, GC, memory, pd$interval / 1.0e6)
    else
        pathSummaryHits(apd, GC, memory)
}

trimOrPad <- function(str, width) {
    ifelse(nchar(str) > width,
           paste(substr(str, 1, width - 4), "..."),
           sprintf("%-*s", width, str))
}

srcSummary <- function(pd, byTotal = TRUE,
                       value = c("pct", "time", "hits"),
                       GC = TRUE, memory = FALSE, total.pct = 0,
                       source = TRUE, width = getOption("width")) {
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

    ## order rows of val and rc by row number within file and alphabetically
    ## by file.
    ord <- order(line)
    val <- val[ord, , drop = FALSE]
    rc <- rc[ord, , drop = FALSE]
    file <- file[ord]
    val <- val[order(file), , drop = FALSE]
    rc <- rc[order(file), , drop = FALSE]

    if (total.pct > 0)
        val <- val[val[,1] >= pd$total * (total.pct / 100), , drop = FALSE]

    if (value == "pct") {
        tot <- pd$total
        colnames(val) <- paste(colnames(val), "pct", sep = ".")
        val <- as.data.frame(percent(val, tot))
    }
    else if (value == "time") {
        delta <- pd$interval / 1.0e6
        colnames(val) <- paste(colnames(val), "time", sep = ".")
        val <- as.data.frame(val * delta)
    }
    else {
        colnames(val) <- paste(colnames(val), "hits", sep = ".")
        val <- as.data.frame(val)
    }

    if (memory && pd$haveMem)
       val$alloc <- rc$alloc

    if (source) {
        pad <- 5 ## allow for padding print.data.frame might insert
        width <- width - max(nchar(capture.output(print(val)))) - pad
        src <- unlist(
            lapply(seq_along(pd$files),
                   function(fn) {
                       fname <- pd$files[fn]
                       which <- refFN(rc$refs) == fn
                       if (file.exists(fname))
                           readLines(fname)[refLN(rc$refs[which])]
                       else
                           rep("<source not available>", sum(which))
               }))
        cbind(val, data.frame(source = trimOrPad(src, width),
                              stringsAsFactors = FALSE))
    }
    else val
}

annotateSource <- function(pd, value = c("pct", "time", "hits"),
                           GC = TRUE, sep = ":  ", show = TRUE, ...) {
    if (! is.null(pd$files)) {
        r <- refCounts(pd)
        file <- refFN(r$ref)
        line <- refLN(r$ref)
        value <- match.arg(value)
        delta <- pd$interval / 1.0e6
        if (GC && pd$haveGC)
            if(value == "pct")
                s <- sprintf(" %5.2f%% %5.2f%%  ",
                             round(100 * (r$total / pd$total), 2),
                             round(100 * (r$gctotal / pd$total), 2))
            else if(value == "time")
                s <- sprintf(" %5.2f %5.2f  ", round(r$total*delta, 2), 
                             round(r$gctotal*delta, 2))
            else
                s <- sprintf(" %5d %5d  ", r$total, r$gctotal)
        else
            if(value == "pct")
                s <- sprintf(" %5.2f%%  ", round(100 * (r$total / pd$total), 2))
            else if(value == "time") 
                s <- sprintf(" %5.2f  ", round(r$total*delta, 2))
            else
                s <- sprintf(" %5d  ", r$total)
        ann <- lapply(seq_along(pd$files),
                      function(fn) {
                          fname <- pd$files[fn]
                          if (file.exists(fname)) {
                              flines <- readLines(fname)
                              b <- rep(paste0(rep(" ", nchar(s[1])),
                                              collapse = ""),
                                       length(flines))
                              b[line[file == fn]] <- s[file == fn]
                              paste(b, flines, sep = sep)
                          }
                          else "<not available>"
                      })
        names(ann) <- pd$files
        if (show) {
            showAnnotation(ann, ...)
            invisible(ann)
        }
        else ann
    }
}

## **** option to only show lines with annotation, plus or minus a few?
## **** show line numbers?
showAnnotation <- function(ann, width = getOption("width"), ...) {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    for (i in seq_along(ann)) {
        a <- ann[[i]]
        if (is.na(width))
            ta <- a
        else
            ta <- ifelse(nchar(a) > width,
                         paste(substr(a, 1, width - 4), "..."),
                         a)
        writeLines(ta, tmp)
        file.show(tmp, ...)
    }
}
