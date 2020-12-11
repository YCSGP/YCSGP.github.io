##------- generate GraphViz code, return unique nodes in the network -----------

genGraphViz <- function(edge,node,title="",filename=""){
  first <- "graph ER {
	start=1058656;
	epsilon=0.0001;
	splines=true;
	
	node [shape=box,fontsize=6, height=0.25,width=0.25, fixedsize=true];
	fontsize=9;"
  cat(first,"\n", file=filename)
  allsymbol <- c(paste(edge[,"genesymbol1"],"\\n",edge[,"p1"],sep=""),
                 paste(edge[,"genesymbol2"],"\\n",edge[,"p2"],sep=""))
  uniquesymbol <- names(table(allsymbol))
  uniqueNode <- node[node[,"genesymbol"] %in%
                     names(table(c(edge[,"genesymbol1"],edge[,"genesymbol2"]))),]

  for(sym in uniquesymbol)
    cat('"',sym,'";\n', file=filename,sep="",append=T);
  
  for(i in 1:nrow(edge)) {
    line <- paste('"',edge[i,"genesymbol1"],"\\n",edge[i,"p1"],'" -- "',
                  edge[i,"genesymbol2"],"\\n",edge[i,"p2"],'"',sep="")
    cat(line,"\n", file=filename,append=T)
  }
  if(title=="")
    last <- paste("      label = ", '"', edge[1,"nwid"], '";\n}', sep="")
  else
    last <- paste("      label = ", '"', title, '";\n}', sep="")
  cat(last, "\n", file=filename,append=T)

  return(uniqueNode)
}

##------------return a vector of number of neighbors -----------
numNeighbour <- function(nNodes,trueEdge){
  ##browser()
  netNeighbour <- rep(0, nNodes)
  nu <- apply(trueEdge, 1, function(edge) {
    netNeighbour[edge[1]] <<- netNeighbour[edge[1]]+1
    netNeighbour[edge[2]] <<- netNeighbour[edge[2]]+1
  } )
  netNeighbour
}

##------------return a list of vectors of  neighbors for each node-----------
allNeighbour <- function(nNodes,trueEdge){
  rst <- lapply(1:nNodes, function(ni) sort(union(trueEdge[trueEdge[,1]==ni,2],trueEdge[trueEdge[,2]==ni,1])))
  names(rst) <- as.character(1:nNodes)
  rst
}
##------------transform from char to factor for nodes and edges ------
## adapted from 04_allpathway_graphviz.r to get edges and nodes
procGraph <- function(nodes,edges) {
  newNodes <- nodes
  newEdges <- aggregate(as.data.frame(edges),edges[,c("genesymbol1","genesymbol2")],FUN=I,simplify=T)[,c(1,2)]
  ##edges <<- read.csv(edgefile,header=T)
  ##nodes <<- read.csv(nodefile,header=T)
  all.levels <- levels(as.factor(nodes$genesymbol))
  ## make common factor levels for genes
  newEdges$genesymbol1 <-  factor(edges$genesymbol1, levels=all.levels)
  newEdges$genesymbol2 <-  factor(edges$genesymbol2, levels=all.levels)
  newNodes$genesymbol <-  factor(nodes$genesymbol, levels=all.levels)

  trueEdges <- cbind(as.numeric(newEdges$genesymbol1),as.numeric(newEdges$genesymbol2))
  ##-----------------ad hoc , transform pvalue to normal(0,1)
  ##ndata <<- abs(qnorm(1-nodes$wtc/2,0,1))
  ##-----------------use p value directly

  newNodes$nNei <- numNeighbour(nrow(newNodes), trueEdges)
  list(nodes=newNodes, edges=newEdges, trueEdges=trueEdges)
}

plot.network.default <- function (x, attrname = NULL, label = network.vertex.names(x), 
    coord = NULL, jitter = TRUE, thresh = 0, usearrows = TRUE, 
    mode = "fruchtermanreingold", displayisolates = TRUE, interactive = FALSE, 
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, pad = 0.2, 
    label.pad = 0.5, displaylabels = !missing(label), boxed.labels = TRUE, 
    label.pos = 0, label.bg = "white", vertex.sides = 8, vertex.rot = 0, 
    arrowhead.cex = 1, label.cex = 1, loop.cex = 1, vertex.cex = 1, 
    edge.col = 1, label.col = 1, vertex.col = 2, label.border = 1, 
    vertex.border = 1, edge.lty = 1, label.lty = NULL, vertex.lty = 1, 
    edge.lwd = 0, label.lwd = par("lwd"), edge.len = 0.5, edge.curve = 0.1, 
    edge.steps = 50, loop.steps = 20, object.scale = 0.01, uselen = FALSE, 
    usecurve = FALSE, suppress.axes = TRUE, vertices.last = TRUE, 
    new = TRUE, layout.par = NULL, ...) 
{
    bellstate <- options()$locatorBell
    expstate <- options()$expression
    on.exit(options(locatorBell = bellstate, expression = expstate))
    options(locatorBell = FALSE, expression = Inf)
    "%iin%" <- function(x, int) (x >= int[1]) & (x <= int[2])
    if (is.hyper(x)) {
        d <- as.matrix.network(x, matrix.type = "incidence", 
            attrname = attrname)
        n <- sum(dim(d))
        temp <- matrix(0, nrow = n, ncol = n)
        if (is.directed(x)) {
            temp[1:dim(d)[1], (dim(d)[1] + 1):n] <- abs(pmin(d, 
                0))
            temp[(dim(d)[1] + 1):n, 1:dim(d)[1]] <- t(abs(pmax(d, 
                0)))
            d <- temp
        }
        else {
            temp[1:dim(d)[1], (dim(d)[1] + 1):n] <- d
            temp[lower.tri(temp)] <- t(temp)[lower.tri(temp)]
            d <- temp
            usearrows <- FALSE
        }
        if (length(label) == network.size(x)) 
            label <- c(label, paste("e", 1:(n - network.size(x)), 
                sep = ""))
    }
    else if (is.bipartite(x)) {
        n <- network.size(x)
        temp <- as.matrix.network(x, matrix.type = "adjacency", 
            attrname = attrname)
        d <- matrix(0, n, n)
        d[1:NROW(temp), (NROW(temp) + 1):NCOL(d)] <- temp
        d[(NROW(temp) + 1):NCOL(d), 1:NROW(temp)] <- t(temp)
        colnames(d) <- c(rownames(temp), colnames(temp))
        rownames(d) <- c(rownames(temp), colnames(temp))
        usearrows <- FALSE
    }
    else {
        n <- network.size(x)
        d <- as.matrix.network(x, matrix.type = "adjacency", 
            attrname = attrname)
        if (!is.directed(x)) 
            usearrows <- FALSE
    }
    diag <- has.loops(x)
    d[is.na(d)] <- 0
    d.raw <- d
    d <- matrix(as.numeric(d > thresh), n, n)
    if (!is.null(coord)) {
        cx <- coord[, 1]
        cy <- coord[, 2]
    }
    else {
        layout.fun <- try(match.fun(paste("network.layout.", 
            mode, sep = "")), silent = TRUE)
        if (class(layout.fun) == "try-error") 
            stop("Error in plot.network.default: no layout function for mode ", 
                mode)
        temp <- layout.fun(d, layout.par)
        cx <- temp[, 1]
        cy <- temp[, 2]
    }
    if (jitter) {
        cx <- jitter(cx)
        cy <- jitter(cy)
    }
    use <- displayisolates | (((apply(d, 1, sum) + apply(d, 2, 
        sum)) > 0))
    if (is.null(xlab)) 
        xlab = ""
    if (is.null(ylab)) 
        ylab = ""
    if (is.null(xlim)) 
        xlim <- c(min(cx[use]) - pad, max(cx[use]) + pad)
    if (is.null(ylim)) 
        ylim <- c(min(cy[use]) - pad, max(cy[use]) + pad)
    xrng <- diff(xlim)
    yrng <- diff(ylim)
    xctr <- (xlim[2] + xlim[1])/2
    yctr <- (ylim[2] + ylim[1])/2
    if (xrng < yrng) 
        xlim <- c(xctr - yrng/2, xctr + yrng/2)
    else ylim <- c(yctr - xrng/2, yctr + xrng/2)
    baserad <- min(diff(xlim), diff(ylim)) * object.scale
    if (new) {
        plot(0, 0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, 
            ylab = ylab, asp = 1, axes = !suppress.axes, ...)
    }
    if (is.character(vertex.cex) && (length(vertex.cex == 1))) {
        temp <- vertex.cex
        vertex.cex <- rep(get.vertex.attribute(x, vertex.cex), 
            length = n)
        if (all(is.na(vertex.cex))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    else vertex.cex <- rep(vertex.cex, length = n)
    vertex.radius <- rep(baserad * vertex.cex, length = n)
    if (is.character(vertex.sides) && (length(vertex.sides == 
        1))) {
        temp <- vertex.sides
        vertex.sides <- rep(get.vertex.attribute(x, vertex.sides), 
            length = n)
        if (all(is.na(vertex.sides))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    else vertex.sides <- rep(vertex.sides, length = n)
    if (is.character(vertex.border) && (length(vertex.border) == 
        1)) {
        temp <- vertex.border
        vertex.border <- rep(get.vertex.attribute(x, vertex.border), 
            length = n)
        if (all(is.na(vertex.border))) 
            vertex.border <- rep(temp, length = n)
        else {
            if (!all(is.color(vertex.border), na.rm = TRUE)) 
                vertex.border <- as.color(vertex.border)
        }
    }
    else vertex.border <- rep(vertex.border, length = n)
    if (is.character(vertex.col) && (length(vertex.col) == 1)) {
        temp <- vertex.col
        vertex.col <- rep(get.vertex.attribute(x, vertex.col), 
            length = n)
        if (all(is.na(vertex.col))) 
            vertex.col <- rep(temp, length = n)
        else {
            if (!all(is.color(vertex.col), na.rm = TRUE)) 
                vertex.col <- as.color(vertex.col)
        }
    }
    else vertex.col <- rep(vertex.col, length = n)
    if (is.character(vertex.lty) && (length(vertex.lty) == 1)) {
        temp <- vertex.lty
        vertex.lty <- rep(get.vertex.attribute(x, vertex.lty), 
            length = n)
        if (all(is.na(vertex.lty))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    else vertex.lty <- rep(vertex.lty, length = n)
    if (is.character(vertex.rot) && (length(vertex.rot) == 1)) {
        temp <- vertex.rot
        vertex.rot <- rep(get.vertex.attribute(x, vertex.rot), 
            length = n)
        if (all(is.na(vertex.rot))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    else vertex.rot <- rep(vertex.rot, length = n)
    if (is.character(loop.cex) && (length(loop.cex) == 1)) {
        temp <- loop.cex
        loop.cex <- rep(get.vertex.attribute(x, loop.cex), length = n)
        if (all(is.na(loop.cex))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    else loop.cex <- rep(loop.cex, length = n)
    if (is.character(label.col) && (length(label.col) == 1)) {
        temp <- label.col
        label.col <- rep(get.vertex.attribute(x, label.col), 
            length = n)
        if (all(is.na(label.col))) 
            label.col <- rep(temp, length = n)
        else {
            if (!all(is.color(label.col), na.rm = TRUE)) 
                label.col <- as.color(label.col)
        }
    }
    else label.col <- rep(label.col, length = n)
    if (is.character(label.border) && (length(label.border) == 
        1)) {
        temp <- label.border
        label.border <- rep(get.vertex.attribute(x, label.border), 
            length = n)
        if (all(is.na(label.border))) 
            label.border <- rep(temp, length = n)
        else {
            if (!all(is.color(label.border), na.rm = TRUE)) 
                label.border <- as.color(label.border)
        }
    }
    else label.border <- rep(label.border, length = n)
    if (is.character(label.bg) && (length(label.bg) == 1)) {
        temp <- label.bg
        label.bg <- rep(get.vertex.attribute(x, label.bg), length = n)
        if (all(is.na(label.bg))) 
            label.bg <- rep(temp, length = n)
        else {
            if (!all(is.color(label.bg), na.rm = TRUE)) 
                label.bg <- as.color(label.bg)
        }
    }
    else label.bg <- rep(label.bg, length = n)
    if (!vertices.last) 
        network.vertex(cx[use], cy[use], radius = vertex.radius[use], 
            sides = vertex.sides[use], col = vertex.col[use], 
            border = vertex.border[use], lty = vertex.lty[use], 
            rot = vertex.rot[use])
    px0 <- vector()
    py0 <- vector()
    px1 <- vector()
    py1 <- vector()
    e.lwd <- vector()
    e.curv <- vector()
    e.type <- vector()
    e.col <- vector()
    e.hoff <- vector()
    e.toff <- vector()
    e.diag <- vector()
    e.rad <- vector()
    if (is.character(edge.col) && (length(edge.col) == 1)) {
        if (edge.col %in% list.edge.attributes(x)) {
            edge.col <- as.matrix.network.adjacency(x, attrname = edge.col)
            if (!all(is.color(edge.col), na.rm = TRUE)) 
                edge.col <- matrix(as.color(edge.col), NROW(edge.col), 
                  NCOL(edge.col))
        }
    }
    if (is.character(edge.lty) && (length(edge.lty) == 1)) {
        temp <- edge.lty
        edge.lty <- as.matrix.network.adjacency(x, attrname = edge.lty)
        if (all(is.na(edge.lty))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    if (is.character(edge.lwd) && (length(edge.lwd) == 1)) {
        temp <- edge.lwd
        edge.lwd <- as.matrix.network.adjacency(x, attrname = edge.lwd)
        if (all(is.na(edge.lwd))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    if (is.character(edge.curve) && (length(edge.curve) == 1)) {
        temp <- edge.curve
        edge.curve <- as.matrix.network.adjacency(x, attrname = edge.curve)
        if (all(is.na(edge.curve))) 
            stop("Attribute", temp, "had illegal missing values or was not present in plot.graph.default.")
    }
    if (!is.array(edge.col)) 
        edge.col <- array(edge.col, dim = dim(d))
    if (!is.array(edge.lty)) 
        edge.lty <- array(edge.lty, dim = dim(d))
    if (!is.array(edge.lwd)) {
        if (edge.lwd > 0) 
            edge.lwd <- array(edge.lwd * d.raw, dim = dim(d))
        else edge.lwd <- array(1, dim = dim(d))
    }
    if (!is.array(edge.curve)) {
        if (!is.null(edge.curve)) 
            edge.curve <- array(edge.curve * d.raw, dim = dim(d))
        else edge.curve <- array(0, dim = dim(d))
    }
    dist <- as.matrix(dist(cbind(cx, cy)))
    tl <- d.raw * dist
    tl.max <- max(tl)
    for (i in (1:n)[use]) for (j in (1:n)[use]) if (d[i, j]) {
        px0 <- c(px0, as.real(cx[i]))
        py0 <- c(py0, as.real(cy[i]))
        px1 <- c(px1, as.real(cx[j]))
        py1 <- c(py1, as.real(cy[j]))
        e.toff <- c(e.toff, vertex.radius[i])
        e.hoff <- c(e.hoff, vertex.radius[j])
        e.col <- c(e.col, edge.col[i, j])
        e.type <- c(e.type, edge.lty[i, j])
        e.lwd <- c(e.lwd, edge.lwd[i, j])
        e.diag <- c(e.diag, i == j)
        e.rad <- c(e.rad, vertex.radius[i] * loop.cex[i])
        if (uselen) {
            if (tl[i, j] > 0) {
                e.len <- dist[i, j] * tl.max/tl[i, j]
                e.curv <- c(e.curv, edge.len * sqrt((e.len/2)^2 - 
                  (dist[i, j]/2)^2))
            }
            else {
                e.curv <- c(e.curv, 0)
            }
        }
        else {
            e.curv <- c(e.curv, edge.curve[i, j])
        }
    }
    if (diag && (length(px0) > 0) && sum(e.diag > 0)) {
        network.loop(as.vector(px0)[e.diag], as.vector(py0)[e.diag], 
            length = 1.5 * baserad * arrowhead.cex, angle = 25, 
            width = e.lwd[e.diag] * baserad/10, col = e.col[e.diag], 
            border = e.col[e.diag], lty = e.type[e.diag], offset = e.hoff[e.diag], 
            edge.steps = loop.steps, radius = e.rad[e.diag], 
            arrowhead = usearrows, xctr = mean(cx[use]), yctr = mean(cy[use]))
    }
    if (length(px0) > 0) {
        px0 <- px0[!e.diag]
        py0 <- py0[!e.diag]
        px1 <- px1[!e.diag]
        py1 <- py1[!e.diag]
        e.curv <- e.curv[!e.diag]
        e.lwd <- e.lwd[!e.diag]
        e.type <- e.type[!e.diag]
        e.col <- e.col[!e.diag]
        e.hoff <- e.hoff[!e.diag]
        e.toff <- e.toff[!e.diag]
        e.rad <- e.rad[!e.diag]
    }
    if (!usecurve & !uselen) {
        if (length(px0) > 0) 
            network.arrow(as.vector(px0), as.vector(py0), as.vector(px1), 
                as.vector(py1), length = 2 * baserad * arrowhead.cex, 
                angle = 20, col = e.col, border = e.col, lty = e.type, 
                width = e.lwd * baserad/10, offset.head = e.hoff, 
                offset.tail = e.toff, arrowhead = usearrows)
    }
    else {
        if (length(px0) > 0) {
            network.arrow(as.vector(px0), as.vector(py0), as.vector(px1), 
                as.vector(py1), length = 2 * baserad * arrowhead.cex, 
                angle = 20, col = e.col, border = e.col, lty = e.type, 
                width = e.lwd * baserad/10, offset.head = e.hoff, 
                offset.tail = e.toff, arrowhead = usearrows, 
                curve = e.curv, edge.steps = edge.steps)
        }
    }
    if (vertices.last) 
        network.vertex(cx[use], cy[use], radius = vertex.radius[use], 
            sides = vertex.sides[use], col = vertex.col[use], 
            border = vertex.border[use], lty = vertex.lty[use], 
            rot = vertex.rot[use])
    if (displaylabels & (!all(label == "")) & (!all(use == FALSE))) {
        if (label.pos == 0) {
            xoff <- cx[use] - mean(cx[use])
            yoff <- cy[use] - mean(cy[use])
            roff <- sqrt(xoff^2 + yoff^2)
            xhat <- xoff/roff
            yhat <- yoff/roff
        }
        else if (label.pos < 5) {
            xhat <- switch(label.pos, 0, -1, 0, 1)
            yhat <- switch(label.pos, -1, 0, 1, 0)
        }
        else {
            xhat <- 0
            yhat <- 0
        }
        os <- par()$cxy * label.cex
        lw <- strwidth(label[use], cex = label.cex)/2
        lh <- strheight(label[use], cex = label.cex)/2
        if (boxed.labels) {
            rect(cx[use] - lw * (1 + label.pad) + xhat * (lw * 
                (1 + label.pad + 0.2) + vertex.radius[use]), 
                cy[use] - lh * (1 + label.pad) + yhat * (lh * 
                  (1 + label.pad + 0.2) + vertex.radius[use]), 
                cx[use] + lw * (1 + label.pad) + xhat * (lw * 
                  (1 + label.pad + 0.2) + vertex.radius[use]), 
                cy[use] + lh * (1 + label.pad) + yhat * (lh * 
                  (1 + label.pad + 0.2) + vertex.radius[use]), 
                col = label.bg, border = label.border, lty = label.lty, 
                lwd = label.lwd)
        }
        text(cx[use] + xhat * (lw * (1 + label.pad + 0.2) + vertex.radius[use]), 
            cy[use] + yhat * (lh * (1 + label.pad + 0.2) + vertex.radius[use]), 
            label[use], cex = label.cex, col = label.col, offset = 0)
    }
    if (interactive && ((length(cx) > 0) && (!all(use == FALSE)))) {
        os <- c(0.2, 0.4) * par()$cxy
        textloc <- c(min(cx[use]) - pad, max(cy[use]) + pad)
        tm <- "Select a vertex to move, or click \"Finished\" to end."
        tmh <- strheight(tm)
        tmw <- strwidth(tm)
        text(textloc[1], textloc[2], tm, adj = c(0, 0.5))
        fm <- "Finished"
        finx <- c(textloc[1], textloc[1] + strwidth(fm))
        finy <- c(textloc[2] - 3 * tmh - strheight(fm)/2, textloc[2] - 
            3 * tmh + strheight(fm)/2)
        finbx <- finx + c(-os[1], os[1])
        finby <- finy + c(-os[2], os[2])
        rect(finbx[1], finby[1], finbx[2], finby[2], col = "white")
        text(finx[1], mean(finy), fm, adj = c(0, 0.5))
        clickpos <- unlist(locator(1))
        if ((clickpos[1] %iin% finbx) && (clickpos[2] %iin% finby)) {
            cl <- match.call()
            ## Min:  add the following line
            for ( nm in setdiff(names(cl),c("","x", "coord"))) cl[[nm]] <- get(nm)
            cl$interactive <- FALSE
            cl$coord <- cbind(cx, cy)
            cl$x <- x
            return(eval(cl))
        }
        else {
            clickdis <- sqrt((clickpos[1] - cx[use])^2 + (clickpos[2] - 
                cy[use])^2)
            selvert <- match(min(clickdis), clickdis)
            if (all(label == "")) 
                label <- 1:n
            rect(textloc[1], textloc[2] - tmh/2, textloc[1] + 
                tmw, textloc[2] + tmh/2, border = "white", col = "white")
            tm <- "Where should I move this vertex?"
            tmh <- strheight(tm)
            tmw <- strwidth(tm)
            text(textloc[1], textloc[2], tm, adj = c(0, 0.5))
            fm <- paste("Vertex", label[use][selvert], "selected")
            finx <- c(textloc[1], textloc[1] + strwidth(fm))
            finy <- c(textloc[2] - 3 * tmh - strheight(fm)/2, 
                textloc[2] - 3 * tmh + strheight(fm)/2)
            finbx <- finx + c(-os[1], os[1])
            finby <- finy + c(-os[2], os[2])
            rect(finbx[1], finby[1], finbx[2], finby[2], col = "white")
            text(finx[1], mean(finy), fm, adj = c(0, 0.5))
            clickpos <- unlist(locator(1))
            cx[use][selvert] <- clickpos[1]
            cy[use][selvert] <- clickpos[2]
            cl <- match.call()
            ## Min:  add the following line
            for ( nm in setdiff(names(cl),c("","x", "coord"))) cl[[nm]] <- get(nm)
            cl$coord <- cbind(cx, cy)
            cl$x <- x
            ##browser()
            return(eval(cl))
        }
    }
    invisible(cbind(cx, cy))
}
