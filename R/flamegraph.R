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

default.cmap <- function(x, ...) {
    y <- unique(x)
    clr <- rainbow(length(y))
    clr[match(x, y)]
}

fgData <- function(data, reorder = c("alpha", "hot", "time"),
                   colormap) {
    stacks <- data$stacks
    mx <- max(sapply(stacks, length))

    reorder <- match.arg(reorder)
    if (reorder == "time") {
        runs <- rle(data$trace)
        stacks <- stacks[runs$values]
        counts <- runs$lengths
    }
    else {
        ## ord <- seq_along(stacks)
        ## for (i in (mx : 1))
        ##     ord <- ord[order(sapply(stacks[ord], `[`, i), na.last = FALSE)]
        counts <- data$counts
        if (reorder == "hot")
            ord <- hotPathOrd(stacks, counts)
        else
            ord <- alphaPathOrd(stacks, counts)
        stacks <- stacks[ord]
        counts <- counts[ord]
    }

    val <- do.call(rbind, lapply(1 : mx, fgDataLine, stacks, counts))
    cmap <- if (! is.null(colormap)) colormap else default.cmap
    fun <- sub(" .*", "", val$label)
    val$col <- cmap(fun, val$top, val$right - val$left)
    val
}

## This is computes the data with fdData and draws the graph with base
## graphics. This is the bit we would need to change for grid, maybe
## ggplot2, and for svg output.
stdFlameGraph <- function(data, reorder, colormap, cex, main) {
    fdg <- fgData(data, reorder, colormap)
    left <- fdg$left
    bottom <- fdg$bottom
    right <- fdg$right
    top <- fdg$top
    col <- fdg$col
    label <- fdg$label

    plot(c(min(left), max(right)), c(min(bottom), max(top)),
         type = "n", axes = FALSE, xlab = "", ylab = "", main = main)

    rect(left, bottom, right, top, col = col)
    
    ## half-em for half-character offset used by text() when pos
    ## argument is used.
    hm <- 0.5 * strwidth("m", cex = cex)

    show <- (strheight(label, cex = cex) <= 0.9 &
             strwidth(label, cex = cex) + 2 * hm <= 0.8 * (right - left))
    if (any(show))
        text(left[show], bottom[show] + 0.4, label[show], pos = 4, cex = cex)

    invisible(structure(list(left = left, right = right,
                             bottom = bottom, top = top,
                             label = label),
                        class = "proftools_flameGraph"))
}

htmlencode <- function(x)
    sub(">", "&gt;", sub("<", "&lt;", x))

svgFlameGraph <- function(file, data, reorder, colormap, main,
                          tooltip) {
    fdg <- fgData(data, reorder, colormap)
    fdg$col <- substr(fdg$col, 1, 7) ## remove 'alpha' component from colors
    mx <- max(fdg$top)
    totalCount <- max(fdg$right)
    counts <- fdg$right-fdg$left
    percents <- round(counts*100/totalCount, 2)
    widths <- round(percents*1180/100, 2)
    y <- 33 + (mx-fdg$top)*16
    x <- 10+round(fdg$left*1180/totalCount, 2)
    col <- fdg$col
    labels <- htmlencode(fdg$label)
    
    svgCode = paste0("<rect x=\"", x, "\" y=\"", y, 
    "\" width=\"", widths, "\" height=\"15.0\" fill=\"", col, 
    "\" rx=\"2\" ry=\"2\" onmouseover=\"s('", labels, " (",
    counts, " samples, ", percents, "%)')\" onmouseout=\"c()\" />")

    show <- (! is.na(labels) & 10*nchar(labels)<widths)
    if (any(show))
        svgCode = append(svgCode, paste0("<text text-anchor=\"\" x=\"", 
        x[show]+3, "\" y=\"", y[show]+10.5, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\" onmouseover=\"s('", 
        labels[show], " (", counts[show], " samples, ", percents[show], 
        "%)')\" onmouseout=\"c()\" >", labels[show],"</text>"))
        
    writeFile(file, svgCode, mx, main, tooltip)
}
## This writes the header of the svg file
writeFile <- function(file, svgCode, mx, main, tooltip){
    write(c(paste0("<?xml version=\"1.0\" standalone=\"no\"?>
    <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">
    <svg version=\"1.1\" width=\"1200\" height=\"", 16*mx+66, "\" onload=\"init(evt)\"", if (tooltip) " onmousemove=\"mouseMove(evt)\" " else " ", "viewBox=\"0 0 1200 ", 16*mx+66, "\" xmlns=\"http://www.w3.org/2000/svg\" >
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
    <![CDATA[",
                  if (tooltip) "
    var bMouseDragging = false;
    var nMouseOffsetX = 0;
    var nMouseOffsetY = 0;
    var detailss;
    bMouseDragging = true;

    var details = document.getElementById(\"details\");
    var detailsBG = document.getElementById(\"detailsBG\");
    if(details) {
        var p = document.documentElement.createSVGPoint();
        p.x = evt.clientX;
        p.y = evt.clientY;
    
        var m = details.getScreenCTM();
        p = p.matrixTransform(m.inverse());
        nMouseOffsetX = p.x - parseInt(details.getAttribute(\"x\"));
        nMouseOffsetY = p.y - parseInt(details.getAttribute(\"y\"));
    }
        
    function mouseMove(evt) { 
        var p = document.documentElement.createSVGPoint();
        p.x = evt.clientX;
        p.y = evt.clientY;
        var bClient = true;
        
        if(bMouseDragging) {
            var details = document.getElementById(\"details\");
            var detailsBG = document.getElementById(\"detailsBG\");
            if(details) {
        
                var m = details.getScreenCTM();
                p = p.matrixTransform(m.inverse());
                var totalOffset = evt.pageY;
                var totalXOffset = evt.pageX;

                details.setAttribute(\"y\", totalOffset - nMouseOffsetY + 18);
                detailsBG.setAttribute(\"y\", totalOffset - nMouseOffsetY);
                
                var bbox = details.getBBox();
                var width = bbox.width;
                var height = bbox.height;
                if(totalXOffset+width >= 1150){
                    details.setAttribute(\"x\", totalXOffset - nMouseOffsetX - 15 - width);
                    detailsBG.setAttribute(\"x\", totalXOffset - nMouseOffsetX - 20 - width);
                }
                else{
                    details.setAttribute(\"x\", totalXOffset - nMouseOffsetX + 20);              
                    detailsBG.setAttribute(\"x\", totalXOffset - nMouseOffsetX + 15);              
                }
                
                detailsBG.setAttribute(\"width\", width+10);
                detailsBG.setAttribute(\"height\", height+10);
                bClient = false;
            }
        }
        
    }
    
    function init() {
        detailss = document.getElementById(\"details\").firstChild; 
        var details = document.getElementById(\"details\");
        var detailsBG = document.getElementById(\"detailsBG\");
        if(details) {
            details.addEventListener(\"mousemove\", mouseMove, false);
        }
    }
    
    function s(info) { detailss.nodeValue = info; }
        function c() { detailss.nodeValue = ' '; }
    ]]>" else "
    var details;
    function init(evt) { details = document.getElementById(\"details\").firstChild; }
    function s(info) { details.nodeValue = info; }
    function c() { details.nodeValue = ' '; }
    ]]>", "
    </script>
    <rect x=\"0.0\" y=\"0\" width=\"1200.0\" height=\"", 16*mx+66, "\" fill=\"url(#background)\"  />
    <text text-anchor=\"middle\" x=\"600\" y=\"24\" font-size=\"17\" font-family=\"Verdana\" font-weight=\"bold\" fill=\"rgb(0,0,0)\"  >", main, "</text>", if (tooltip) "" else paste0("
    <text text-anchor=\"left\" x=\"10\" y=\"", 16*mx+50, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\"  >Function:</text>
    <text text-anchor=\"\" x=\"70\" y=\"", 16*mx+50, "\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\" id=\"details\" > </text>")), svgCode, if (tooltip) "
    <rect x=\"10\" y=\"129\" width=\"1\" height=\"1.0\" fill=\"rgb(250,250,250)\" rx=\"2\" ry=\"2\" id=\"detailsBG\" />
    <text x=\"70\" y=\"162\" font-size=\"12\" font-family=\"Verdana\" fill=\"rgb(0,0,0)\" id=\"details\" > </text>" else "", "
    </svg>"), file = file)
}

## Merge non-NA refs into stacks. Add UnknownFunToken for leaf refs.
refStacks <- function(d) {
    rs <-function(s, r, d) {
        if (is.na(r[length(r)]))
            r <- r[-length(r)]
        else
            s <- c(s, UnknownFunToken)
        fl <- basename(d$files[refFN(r)])
        ln <- refLN(r)
        paste0(s, ifelse(is.na(r), "", sprintf(" (%s:%d)", fl, ln)))
    }
    lapply(1:length(d$stacks), function(i) rs(d$stacks[[i]], d$refs[[i]], d))
}

## produce a flame graph from an Rprof file.
## order = "time" produces a graph like profr.
flameGraph <- function(pd, svgfile, order = c("hot", "alpha", "time"),
                       colormap = NULL, srclines = FALSE, cex = 0.75,
                       main = "Call Graph", tooltip = FALSE) {
    order <- match.arg(order)
    if (is.character(pd))
        pd <- readPD(pd)
    data <- list(counts = pd$counts,
                 stacks = if (srclines) refStacks(pd) else pd$stacks,
                 trace = pd$trace)
    if (! missing(svgfile))
        svgFlameGraph(svgfile, data, order, colormap, main, tooltip)
    else
        stdFlameGraph(data, order, colormap, cex, main)
}

fgIdentify <- function(p) {
    loc <- locator(1)
    if (! is.null(loc)) {
        idx <- which(loc$x >= p$left & loc$x <= p$right &
                     loc$y >= p$bottom & loc$y <= p$top)
        if (length(idx) > 0)
            p$label[idx]
    }
}

fgIdentify <- function(p, n = 1, print = FALSE) {
    val <- NULL
    while (n > 0) {
        n <- n - 1
        loc <- locator(1)
        if (! is.null(loc)) {
            idx <- which(loc$x >= p$left & loc$x <= p$right &
                         loc$y >= p$bottom & loc$y <= p$top)
            if (length(idx) > 0) {
                if (print)
                    cat(p$label[idx], "\n")
                val <- c(val, p$label[idx])
            }
            else break
        }
        else break
    }
    val
}

## Outline clicked things (could redraw before outlining -- would need colors):
fgIdentify <- function(x, n = 1, print = FALSE, outline = FALSE, ...) {
    p <- x
    val <- NULL
    while (n > 0) {
        n <- n - 1
        loc <- locator(1)
        if (! is.null(loc)) {
            idx <- which(loc$x >= p$left & loc$x <= p$right &
                         loc$y >= p$bottom & loc$y <= p$top)
            if (length(idx) > 0) {
                if (outline) {
                    allIDX <- p$label %in% p$label[idx]
                    pp <- lapply(p, `[`, allIDX)
                    rect(pp$left, pp$bottom, pp$right, pp$top, lwd = 3)
                }
                if (print)
                    cat(p$label[idx], "\n")
                val <- c(val, p$label[idx])
            }
            else break
        }
        else break
    }
    val
}

identify.proftools_flameGraph <- fgIdentify
