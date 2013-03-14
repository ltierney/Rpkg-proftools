readPD <- function(file) {
    ## read, check, and record header
    header <- readLines(file, 1)

    pat <- ".*sample.interval=([[:digit:]]+).*"
    if (grepl(pat, header))
        interval <- as.integer(sub(pat, "\\1", header))
    else
        stop("not a valid Rprof file")

    haveMem <- grepl("memory profiling:", header)
    haveGC <- grepl("GC profiling:", header)
    haveLines <- grepl("line profiling:", header)

    if (haveMem)
        stop("memory profiling is currently not supported")

    ## read lines and filter out file information
    stacks <- readLines(file)[-1]
    if (haveLines) {
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

    list(gc = haveGC, lines = haveLines, files = files,
         stacks = ustacks, trace = trace, inGC = inGC)
}

cleanStacks <- function(stacks) {
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
    ## Maybe split out file numbers, line numbers
    list(stacks = stacks, refs = stackrefs)
}

countStacks <- function(trace, inGC)
    list(counts = as.numeric(table(trace)),
         gccounts = as.numeric(tapply(inGC, trace, sum)))

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

recodeRefs <- function(refs, files) {
    files <- basename(files)
    recode <- function(refs) {
        fi <- as.integer(sub("([[:digit:]]+)#[[:digit:]]+", "\\1", refs))
        li <- as.integer(sub("[[:digit:]]+#([[:digit:]]+)", "\\1", refs))
        ifelse(is.na(refs), NA, paste(files[fi],  li, sep = ":"))
    }
    sapply(refs, recode)
}

formatTrace <- function(trace, maxlen = 50) {
    out <- paste(trace, collapse = " -> ")
    while (nchar(out) > maxlen && length(trace) > 1) {
        trace <- trace[-1]
        out <- paste("... ", paste(trace, collapse = " -> "))
    }
    out
}

d <- readPD("Rprof-lmfit-new.out")
s <- cleanStacks(d$stacks)
ct <- countStacks(d$trace, d$inGC)
countHits(s$stacks, ct)

countSelfHits(s$stacks, ct)
countSelfHits(s$refs, ct)

## this drops the final ref
v <- mapply(function (x, y) ifelse(is.na(y), x, paste(x, y)), s$stacks, s$refs)
v <- mapply(function (x, y)
            ifelse(is.na(y), x, paste0(x, " (", y,")")),
            s$stacks, recodeRefs(s$refs, d$files))
countHits(v, ct)

## this attributes the final ref to "<Unknown>"
v <- mapply(function (x, y) {
            x <- c(x, "<Unknown>")
            ifelse(is.na(y), x, paste0(x, " (", y,")"))
            },
            s$stacks, recodeRefs(s$refs, d$files))
countHits(v, ct)

printPaths <- function(stacks, counts, n) {
    ord = rev(order(counts$counts))
    if (! missing(n) && length(ord) > n)
        ord <- ord[1 : n]
    tot <- sum(counts$counts)
    pct <- round(100 * counts$counts[ord] / tot, 1)
    gcpct <- round(100 * counts$gccounts[ord] / tot, 1)
    paths <- sapply(stacks[ord], formatTrace)
    mapply(function(x, y, z) cat(sprintf("%5.2f %5.2f   %s\n", x, y, z)),
           pct, gcpct, paths)
    invisible(NULL)
}

## Hot paths -- print ncely, but also allow examination of source refs and such
## Maybe print with source refs?
## Reorder paths, revise time index, so hottest one is first?

## need call counts, fun -> fun and site->fun
## need both since multiple sites, same fun can happen in same stack trace
## collapse cycles within stack trace?

## show hot files, hot lines, calls within line
## show hot paths -- tree view with collapsing of some kind
## show call graph
## show call tree
## flame graph

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
    s <- cleanStacks(d$stacks)
    ct <- countStacks(d$trace, d$inGC)
    stacks <- s$stacks
    counts <- ct$counts
    flameGraph(stacks, counts)
}

## produce a time graph (like profr) from an Rprof file
tg <- function(file) {
    d <- readPD(file)
    s <- cleanStacks(d$stacks)
    r <- rle(d$trace)
    stacks <- s$stacks[r$values]
    counts <- r$lengths
    flameGraph(stacks, counts, FALSE)
}    

lineCalls <- function(i) {
    line <- s$stacks[[i]]
    refs <- s$refs[[i]]
    n <- length(line)
    if (n > 0) {
        caller <- line[-n]
        callee <- line[-1]
        site <- refs[-c(1, n + 1)]
    }
    else
        caller <- callee <- site <- character()
    unique(cbind(caller, callee, site))
}

leafCall <- function(i) {
    line <- s$stacks[[i]]
    refs <- s$refs[[i]]
    n <- length(line)
    if (n > 0) {
        caller <- line[n - 1]
        callee <- line[n]
        site <- refs[n]
    }
    else
        caller <- callee <- site <- character()
    cbind(caller, callee, site)
}

calls <- lapply(lapply(seq_along(s$stacks), lineCalls), as.data.frame)
reps <- unlist(lapply(calls, nrow))

cdf <- do.call(rbind, calls)
cdf$counts <- rep(ct$counts, reps)
cdf$gccounts <- rep(ct$gccounts, reps)

f1 <- cdf$caller
f2 <- cdf$callee
fs <- factor(as.character(cdf$site), exclude = "")

aggregate(cdf[-(1:3)], list(caller = f1, callee = f2, site = fs), sum)

## **** need to optionally include where  of caller

do.call(rbind, lapply(seq_along(s$stacks), leafCall))

## **** compute total. self, gctotal, gcself for calls w and w/o sites
## **** figure out how to write out callgrind from this
## **** figure out how to generate call graphs as in proftools

## **** wrap into function
## **** make result be char not factor
useSite <- FALSE

stacks <- s$stacks
refs <- s$refs
counts <- ct$counts
gccounts <- ct$gccounts

lineFuns <- function(i) {
    line <- stacks[[i]]
    linerefs <- refs[[i]]
    if (useSite) {
        n <- length(line)
        site <- linerefs[-(n + 1)]
    }
    else site <- NA_character_
    unique(cbind(fun = line, site))
}

leafFun <- function(i) {
    line <- stacks[[i]]
    linerefs <- refs[[i]]
    n <- length(line)
    fun <- line[n]
    site <- if (useSite) linerefs[n] else NA_character_
    cbind(fun, site)
}

funs <- lapply(lapply(seq_along(stacks), lineFuns), as.data.frame)
reps <- unlist(lapply(funs, nrow))

fdf <- do.call(rbind, funs)
fdf$total <- rep(counts, reps)
fdf$gctotal <- rep(gccounts, reps)

ff <- fdf$fun
fs <- factor(as.character(fdf$site), exclude = "")

afdf <- aggregate(fdf[-(1:2)], list(fun = ff, site = fs), sum)

sfdf <- as.data.frame(do.call(rbind, lapply(seq_along(stacks), leafFun)))
sfdf$self <- counts
sfdf$gcself <- gccounts

sff <- sfdf$fun
sfs <- factor(as.character(sfdf$site), exclude = "")

asfdf <- aggregate(sfdf[-(1:2)], list(fun = sff, site = sfs), sum)

mfdf <- merge(afdf, asfdf, all = TRUE)
mfdf$self[is.na(mfdf$self)] <- 0
mfdf$gcself[is.na(mfdf$gcself)] <- 0
