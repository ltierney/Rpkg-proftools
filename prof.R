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
## Reorderp paths, revise time index, so hottest one is first?

## need call counts, fun -> fun and site->fun
## need both since multiple sites, same fun can happen in same stack trace
## collapse cycles within stack trace?

## show hot files, hot lines, calls within line
## show hot paths -- tree view with collapsing of some kind
## show call graph
## show call tree
## flame graph


