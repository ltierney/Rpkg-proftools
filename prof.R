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

