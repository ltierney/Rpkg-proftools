alphaPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1))
        ord <- ord[order(sapply(stacks[ord], `[`, i), na.last = FALSE)]
    ord
}

hotPathOrdOne <- function(k, stacks, counts, ord) {
    s <- sapply(stacks[ord], `[`, k)
    v <- aggregate(list(total = counts[ord]), list(fun = s), sum)
    order(v$total[match(s, v$fun)], decreasing = TRUE, na.last=FALSE)
}

hotPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1)) {
        ord1 <- hotPathOrdOne(i, stacks, counts, ord)
        ord <- ord[order(ord1, na.last = FALSE)]
    }
    ord
}

hotPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1)) {
        s <- sapply(stacks[ord], `[`, i)
        v <- aggregate(list(total = -counts[ord]), list(fun = s), sum)
        ss <- v$total[match(s, v$fun)]
        ord <- ord[order(ss, na.last = TRUE)]
    }
    ord
}

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

fgData <- function(stacks, counts, reorder = TRUE) {
    mx <- max(sapply(stacks, length))

    if (is.logical(reorder))
        reorder <- if (reorder) "alpha" else "no"
    if (reorder != "no") {
        if (reorder == "hot")
            ord <- hotPathOrd(stacks, counts)
        else
            ord <- alphaPathOrd(stacks, counts)
        stacks <- stacks[ord]
        counts <- counts[ord]
    }
    
    do.call(rbind, lapply(1 : mx, fgDataLine, stacks, counts))
}

fg <- function(file, reorder = "alpha") {
    if (is.character(file))
        d <- readPD(file)
    else
        d <- file
    flameGraph(d$stacks, d$counts, reorder)
}

## **** color according to fun?
## remap stacks, refs, counts, gccounts
## remap trace
## stacks <- stacks[ord]
## counts <- counts[ord]
## refs <- refs[ord]
## ... trace ...
