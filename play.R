alphaPathOrd <- function(stacks, counts) {
    mx <- max(sapply(stacks, length))
    ord <- seq_along(stacks)
    for (i in (mx : 1))
        ord <- ord[order(sapply(stacks[ord], `[`, i), na.last = FALSE)]
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
    
    val <- do.call(rbind, lapply(1 : mx, fgDataLine, stacks, counts))
    if (FALSE) {
        nm <- unique(val$label)
        nc <- length(nm)
        cm <- rainbow(nc)
        cm <- sample(heat.colors(nc), nc)
        val$col <- cm[match(val$label, nm)]
        nr <- nrow(val)
        val$col <- sample(heat.colors(nr), nr)
    }
    val
}

fg <- function(file, reorder = "alpha") {
    if (is.character(file))
        d <- readPD(file)
    else
        d <- file
    flameGraph(d$stacks, d$counts, reorder)
}

## **** reduce strheight limit more?
## **** subsetting variations; record total counts
## **** heat.colors;
## **** color according to fun?
## remap stacks, refs, counts, gccounts
## remap trace
## stacks <- stacks[ord]
## counts <- counts[ord]
## refs <- refs[ord]
## ... trace ...
