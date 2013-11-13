splitcalls <- function(pd) {
    stacks <- pd$stacks
    counts <- pd$counts
    first <- sapply(stacks, `[`, 1)
    spstacks <- split(stacks, first)
    spcounts <- split(counts, first)
    totals <- lapply(spcounts, sum)  ## maybe sapply?
    merge <- function(total, stacks, counts) {
        keep <- sapply(stacks, function(s) length(s) > 1)
        stacks <- lapply(stacks[keep], `[`, -1)
        counts <- counts[keep]
        list(total = total, stacks = stacks, counts = counts)
    }
    mapply(merge, totals, spstacks, spcounts, SIMPLIFY=FALSE)
}

makeTree <- function(fun, pd) {
    if (length(pd$stacks) > 0) {
        spcalls <- splitcalls(pd)
        funs <- names(spcalls)
        calls <- lapply(seq_along(spcalls),
                        function(i) makeTree(funs[i], spcalls[[i]]))
    }
    else calls <- NULL
    list(fun = fun, hits = pd$total, calls = calls)
}

drawNode <- function(n, left, bottom, right, top, cex = 0.75, rotated = FALSE,
                     showLabel = TRUE) {
    pad <- strwidth("M", cex = cex)* 0.3 ## 0.01
    hits <- n$hits
    ##**** fill with color
    ## rect(left, bottom, right, top)
    col <- rgb(runif(1), runif(1), runif(1))
    rect(left, bottom, right, top, col = col, border = NA)
    ##**** do not write text if rect is too small
    if (showLabel) {
        label <- n$fun
        wd <- strwidth(label, cex = cex) + 2 * pad
        ht <- strheight(label, cex = cex) + 2 * pad
        if (rotated) {
            if (wd <= top - bottom && ht <= right - left)
                text(right - pad, top - pad, label, adj = c(0, 1), cex = cex, srt = -90)
        }
        else {
            if (wd <= right - left && ht <= top - bottom)
                text(left + pad, top - pad, label, adj = c(0, 1), cex = cex)
        }
    }
}

drawTree <- function(n, left, bottom, right, top, cex = 0.75, rotated = TRUE) {
    drawNode(n, left, bottom, right, top, cex, rotated)
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    if (nc > 0) {
        callhits <- sapply(calls, function(x) x$hits) / hits
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        cumhits <- cumsum(callhits)
        pad <- min(right - left, top - bottom) * 0.01
        left <- left + pad
        bottom <- bottom + pad
        right <- right - pad
        top <- top - pad
        if (rotated) {
            bottoms <- bottom + (top - bottom) * c(0, cumhits[-nc])
            tops <- bottom + (top - bottom) * cumhits
            lefts <- rep(left, nc)
            rights <- rep(right, nc)
        }
        else {
            bottoms <-rep(bottom, nc)
            tops <- rep(top, nc)
            lefts <- left + (right - left) * c(0, cumhits[-nc])
            rights <- left + (right - left) * cumhits
        }
        mapply(function(n, left, bottom, right, top)
                   drawTree(n, left, bottom, right, top, cex, ! rotated),
               calls, lefts, bottoms, rights, tops)
    }
    NULL
}

## split on longer edge; use longer edge for label of leaf calls
drawTree <- function(n, left, bottom, right, top, cex = 0.75) {
    rotated <- (right - left) < (top - bottom)
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    drawNode(n, left, bottom, right, top, cex,
             if (nc > 0) ! rotated else rotated)
    if (nc > 0) {
        callhits <- sapply(calls, function(x) x$hits) / hits
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        cumhits <- cumsum(callhits)
        pad <- min(right - left, top - bottom) * 0.01
        left <- left + pad
        bottom <- bottom + pad
        right <- right - pad
        top <- top - pad
        if (rotated) {
            bottoms <- bottom + (top - bottom) * c(0, cumhits[-nc])
            tops <- bottom + (top - bottom) * cumhits
            lefts <- rep(left, nc)
            rights <- rep(right, nc)
        }
        else {
            bottoms <-rep(bottom, nc)
            tops <- rep(top, nc)
            lefts <- left + (right - left) * c(0, cumhits[-nc])
            rights <- left + (right - left) * cumhits
        }
        mapply(function(n, left, bottom, right, top)
                   drawTree(n, left, bottom, right, top, cex),
               calls, lefts, bottoms, rights, tops)
    }
    NULL
}

## Keep label at top if it fits; only show if enough self time
drawTree <- function(n, left, bottom, right, top, cex = 0.75) {
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    callhits <- if (nc > 0) sapply(calls, function(x) x$hits) / hits else 0
    selffrac <- 1 - sum(callhits)
    drawNode(n, left, bottom, right, top, cex, FALSE, FALSE)
    lpad <- strwidth("M", cex = cex)* 0.3 ## 0.01
    label <- n$fun
    wd <- strwidth(label, cex = cex) + 2 * lpad
    ht <- strheight(label, cex = cex) + 2 * lpad
    if (wd <= right - left && ht <= selffrac * (top - bottom)) {
        text(left + lpad, top - lpad, label, adj = c(0, 1), cex = cex)
        top <- top - selffrac * (top - bottom)
    }
    else if (wd <= top - bottom && ht <= selffrac * (right - left)) {
        text(right - lpad, top - lpad, label, adj = c(0, 1), cex = cex, srt = -90)
right <- right - selffrac * (right - left)
    }
    else top <- top - selffrac * (top - bottom)
    if (nc > 0) {
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        pad <- min(right - left, top - bottom) * 0.01
        left <- left + pad
        bottom <- bottom + pad
        right <- right - pad
        top <- top - pad
        width <- right - left
        height <- top - bottom
        cumhits <- cumsum(callhits) / sum(callhits)
        if (width < height) {
            bottoms <- bottom + height * c(0, cumhits[-nc])
            tops <- bottom + height * cumhits
            lefts <- rep(left, nc)
            rights <- rep(right, nc)
        }
        else {
            bottoms <-rep(bottom, nc)
            tops <- rep(top, nc)
            lefts <- left + width * c(0, cumhits[-nc])
            rights <- left + width * cumhits
        }
        mapply(function(n, left, bottom, right, top)
                   drawTree(n, left, bottom, right, top, cex),
               calls, lefts, bottoms, rights, tops)
    }
    NULL
}

drawTree <- function(n, left, bottom, right, top, cex = 0.75) {
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    callhits <- if (nc > 0) sapply(calls, function(x) x$hits) / hits else 0
    selffrac <- 1 - sum(callhits)
    drawNode(n, left, bottom, right, top, cex, FALSE, FALSE)
    lpad <- strwidth("M", cex = cex)* 0.3 ## 0.01
    label <- n$fun
    wd <- strwidth(label, cex = cex) + 2 * lpad
    ht <- strheight(label, cex = cex) + 2 * lpad
    if (wd <= right - left && ht <= selffrac * (top - bottom)) {
        text(left + lpad, top - lpad, label, adj = c(0, 1), cex = cex)
        top <- top - selffrac * (top - bottom)
    }
    else if (wd <= top - bottom && ht <= selffrac * (right - left)) {
        text(right - lpad, top - lpad, label, adj = c(0, 1),
             cex = cex, srt = -90)
        right <- right - selffrac * (right - left)
    }
    else top <- top - selffrac * (top - bottom)
    if (nc > 0) {
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        pad <- min(right - left, top - bottom) * 0.01
        s <- splitRect(callhits,
                       left + pad, bottom + pad, right - pad, top - pad)
        mapply(function(n, left, bottom, right, top)
                   drawTree(n, left, bottom, right, top, cex),
               calls, s$left, s$bottom, s$right, s$top)
    }
    NULL
}

splitRect <- function(score, left, bottom, right, top) {
    stopifnot(all(score > 0))
    cumscore <- cumsum(score) / sum(score)
    nc <- length(cumscore)
    width <- right - left
    height <- top - bottom
    if (width < height) {
        bottoms <- bottom + height * c(0, cumscore[-nc])
        tops <- bottom + height * cumscore
        lefts <- rep(left, nc)
        rights <- rep(right, nc)
    }
    else {
        bottoms <-rep(bottom, nc)
        tops <- rep(top, nc)
        lefts <- left + width * c(0, cumscore[-nc])
        rights <- left + width * cumscore
    }
    list(left = lefts, bottom = bottoms, right = rights, top = tops)
}

testRect <- function(score, left, bottom, right, top) {
    plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "", axes = FALSE)
    rect(left, bottom, right, top)
    s <- splitRect(score, left, bottom, right, top)
    mapply(rect, s$left, s$bottom, s$right, s$top)
    NULL
}

calleeTreeMap <- function(pd, cex = 0.6) {
    plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "", axes = FALSE)
    drawTree(makeTree("", pd), 0, 0, 1, 1, cex = cex)
}

squarifiedTiles <- function(v, left, bottom, right, top) {
    squarify <- function(children, row,  w) {
        child <- children[1] ## **** check for length 0?
        newRow <- c(row, child)
        if (length(row) == 0 || worst(row, w) >= worst(newRow, w)) {
            remainingChildren <- children[-1]
            if (length(remainingChildren) == 0)
                layoutrow(newRow)
            else
                squarify(remainingChildren, newRow, w)
        }
        else {
            layoutrow(row)
            squarify(children, numeric(0), width())
        }
    }

    width <- function() min(right - left, top - bottom)

    layoutrow <- function(r) {
        if (right - left >= top - bottom) {
            a <- sum(r) / (top - bottom)
            b <- r / a
            rects <<- rbind(rects,
                            cbind(left,
                                  c(bottom, bottom + cumsum(b)[-length(r)]),
                                  left + a,
                                  bottom + cumsum(b)))
            left <<- left + a
        }
        else {
            b <- sum(r) / (right - left)
            a <- r / b
            rects <<- rbind(rects,
                            cbind(c(left, left + cumsum(a)[-length(r)]),
                                  bottom,
                                  left + cumsum(a),
                                  bottom + b))
            bottom <<- bottom + b
        }
    }

    worst <- function(R, w) {
        s <- sum(R)
        max(w^2 * R / s^2, s^2 / (w^2 * R))
    }

    stopifnot(all(v > 0) && all(diff(v) <= 0))

    v <- (v / sum(v)) * (right - left) * (top - bottom)

    rects <- NULL
    squarify(v, numeric(0), min(right - left, top - bottom))
    list(left = rects[,1], bottom = rects[,2],
         right = rects[,3], top = rects[,4])
}

makeTreeMapData <- function(n, left, bottom, right, top,
                            cex = 0.75, depth = 1, tile = splitRect) {
    hits <- n$hits
    calls <- n$calls
    nc <- length(calls)
    callhits <- if (nc > 0) sapply(calls, function(x) x$hits) / hits else 0
    selffrac <- 1 - sum(callhits)
    lpad <- strwidth("M", cex = cex)* 0.3 ## 0.01
    label <- n$fun
    wd <- strwidth(label, cex = cex) + 2 * lpad
    ht <- strheight(label, cex = cex) + 2 * lpad
    if (wd <= right - left && ht <= selffrac * (top - bottom)) {
        showLabel <- TRUE
        rotate <- FALSE
        labX <- left + lpad
        labY <- top - lpad
        sTop <- top - selffrac * (top - bottom)
        sRight <- right
    }
    else if (wd <= top - bottom && ht <= selffrac * (right - left)) {
        showLabel <- TRUE
        rotate <- TRUE
        labX <- right - lpad
        labY <- top - lpad
        sTop <- top
        sRight <- right - selffrac * (right - left)
    }
    else {
        showLabel <- FALSE
        rotate <- FALSE
        labX <- 0
        labY <- 0
        sTop <- top - selffrac * (top - bottom)
        sRight <- right
    }
    if (nc > 0) {
        calls <- calls[order(callhits, decreasing = TRUE)]
        callhits <- sort(callhits, decreasing = TRUE)
        pad <- min(right - left, top - bottom) * 0.01
        s <- tile(callhits, left + pad, bottom + pad, sRight - pad, sTop - pad)
        rest <- mapply(function(n, left, bottom, right, top)
                       makeTreeMapData(n, left, bottom, right, top,
                                       cex, depth + 1, tile),
                       calls, s$left, s$bottom, s$right, s$top,
                       SIMPLIFY = FALSE)
    }
    else rest <- NULL
    d <- data.frame(left = left, bottom = bottom, right = right, top = top,
                    hits = hits, depth = depth,
                    label = label, labX = labX, labY = labY,
                    showLabel = showLabel, rotate = rotate,
                    stringsAsFactors = FALSE)
    if (is.null(rest)) d
    else rbind(d, do.call(rbind, rest))
}

calleeTreeMap <- function(pd, srclines = FALSE, cex = 0.75, colormap = NULL,
                          main = "Callee Tree Map", squarify = FALSE) {
    plot(c(0,1), c(0,1), type = "n", xlab = "", ylab = "",
         axes = FALSE, main = main)
    if (srclines)
        pd$stacks <- refStacks(pd)
    tile <- if (squarify) squarifiedTiles else splitRect
    v <-makeTreeMapData(makeTree("", pd), 0, 0, 1, 1, tile = tile)
    nc <- nrow(v)
    cmap <- if (! is.null(colormap)) colormap else default.cmap
    fun <- sub(" .*", "", v$label)
    col <- cmap(fun, val$depth, val$hits)
    rect(v$left, v$bottom, v$right, v$top, col = col)
    vlr <- subset(v, showLabel & rotate)
    if (nrow(vlr) > 0)
        with(vlr, text(labX, labY, label, srt = -90, adj = c(0, 1), cex = cex))
    vlnr <- subset(v, showLabel & ! rotate)
    if (nrow(vlnr) > 0)
        with(vlnr, text(labX, labY, label, adj = c(0, 1), cex = cex))
    invisible(structure(v, class = c("proftools_calleeTreeMap", class(v))))
}

ctmIdentify <- function(p, n = 1, print = FALSE) {
    val <- list()
    while (n > 0) {
        n <- n - 1
        loc <- locator(1)
        if (! is.null(loc)) {
            idx <- which(loc$x >= p$left & loc$x <= p$right &
                         loc$y >= p$bottom & loc$y <= p$top)
            if (length(idx) > 0) {
                stack <- p$label[idx][-1]
                if (print)
                    cat(stack, "\n")
                val <- c(val, stack)
            }
            else break
        }
        else break
    }
    val
}

identify.proftools_calleeTreeMap <- ctmIdentify
