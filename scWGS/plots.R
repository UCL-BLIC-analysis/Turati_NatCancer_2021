#' Title
#'
#' @param data 
#'
#' @return
#' @export chr.lengths data.frame with "length", "from" and "to" for each chromosome
#'
#' @examples
get.chr.lengths <- function(data) {
  chr.lengths <- aggregate(data$END, by=list(data$CHR), max, simplify = T)
  chr.lengths
  chr.lengths[, 1] <- factor(chr.lengths[, 1])
  chr.lengths
  chr.lengths <- chr.lengths[unique(as.factor(data$CHR)), ]
  chr.lengths
  offset <- chr.lengths
  offset[, 3] <- cumsum(as.numeric(offset[, 2])) - offset[, 2]
  offset[, 4] <- offset[, 3] + offset[, 2]
  chr.lengths <- offset[, 2:4]
  colnames(chr.lengths) <- c("length", "from", "to")
  rownames(chr.lengths) <- offset[, 1]
  rm(offset)

  return(chr.lengths)  
}

plotFrequency <- function(data, showBins = FALSE, showTics = TRUE, showLevels = TRUE) {
  chr.lengths <- get.chr.lengths(data)
  
  x.left <- data$START
  x.right <- data$END
  
  x.left <- x.left + chr.lengths[as.character(data$CHR), "from"]
  x.right <- x.right + chr.lengths[as.character(data$CHR), "from"]
  
  n.samples <- ncol(data) - 3
  par("mar" = c(0, 3, 3, 0) + 0.1)
  plot(NA, xlim = c(min(x.left), max(x.right)), ylim = c(0, n.samples),
       axes = F, xlab = NA, ylab = NA)
  x.margin <-  (max(x.right) - min(x.left)) / 500
  y.margin <-  strheight("legend") * 1.1
  usr <- c(min(x.left), max(x.right) + x.margin, -y.margin, n.samples + y.margin * 1.5)
  par("usr" = usr)
  # rect(usr[1], usr[3], usr[2], usr[4], border = "red", col = "pink")
  rect(chr.lengths[, "from"], 0, chr.lengths[, "to"], n.samples, border = NA, col = c("white", rgb(0.9, 0.9, 0.9))[2 - (1:nrow(chr.lengths)) %% 2])
  if (showBins) {
    segments(x.right, 0, x.right, n.samples, col = rgb(0.7, 0.7, 0.7), lwd = min(25 / length(x.right), 0.1))
  }
  if (showTics) {
    chr <- unique(as.character(data$CHR))
    if (length(chr) == 1) {
      if (chr.lengths[chr, "length"] > 10000000) {
        tics <- seq(10000000, chr.lengths[chr, "length"], by = 10000000)
        tics <- tics + chr.lengths[chr, "from"]
        segments(tics, -y.margin / 4, tics, 0, col = "black", lwd = 1)
      }
    }
  }
  if (showLevels) {
    lapply(1:19, function(x) {lines(c(min(x.left), max(x.right)), c(n.samples * x / 20, n.samples * x / 20), col = "grey", lty = 3)})
  }
  
  s.sc <- apply(data[, 4:ncol(data)],  1, function(x) {return(c(sum(x == 0), sum(x == 1), sum(x == 2), sum(x == 3), sum(x >= 4)))})
  rect(x.left, 0, x.right, s.sc[1, ], col = rgb(0, 0, 0.8), border = NA)
  rect(x.left, s.sc[1, ], x.right, colSums(s.sc[1:2, ]), col = rgb(0.6, 0.6, 1), border = NA)
  rect(x.left, colSums(s.sc[1:2, ]), x.right, colSums(s.sc[1:3, ]), col = rgb(0, 1, 1, 0), border = NA)
  rect(x.left, colSums(s.sc[1:3, ]), x.right, colSums(s.sc[1:4, ]), col = rgb(1, 0.6, 0.6), border = NA)
  rect(x.left, colSums(s.sc[1:4, ]), x.right, colSums(s.sc[1:5, ]), col = rgb(0.8, 0, 0), border = NA)
  if (nrow(chr.lengths) > 1) {
    lines(as.vector(t(cbind(chr.lengths[-c(1), "from"], chr.lengths[-c(1), "from"], NA))),
          rep(c(0, n.samples, NA), times = nrow(chr.lengths) - 1),
          col = "grey", lwd = 1.5)
  }
  rect(min(chr.lengths[, "from"]), 0, max(chr.lengths[, "to"]), n.samples, border = "black", col = NA)
  text((chr.lengths[, "from"] + chr.lengths[, "to"]) / 2, -y.margin / 2, labels = sub("chr", "", row.names(chr.lengths)), cex = 0.7)
  axis(2, at = c(0, 0.5, 1) * n.samples, labels = c("0%", "50%", "100%"), las = 2, tick = TRUE)

  x.usr.length = usr[2] - usr[1]
  box.height <- y.margin
  box.width <- box.height * par("pin")[2] / par("pin")[1] * (usr[2] - usr[1]) / (usr[4] - usr[3])
  unit.width <- box.width * 1.5 + strwidth("n=4+")
  separation.width <- box.width * 1.5
  total.width <- 4 * unit.width + 3 * separation.width
  x.usr.start <- usr[1] + (x.usr.length - total.width) / 2
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0, 0, 0.8))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=0", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.6, 0.6, 1))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=1", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(1, 0.6, 0.6))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=3", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.8, 0, 0))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=4+", adj = c(0, 1))
}

get_annotation_palette <- function() {
  # library(RColorBrewer)
  return(c("blue", "red", "darkorange", "green3", "darkcyan", "magenta", "darkgrey", "purple", "brown",
           "steelblue1", "tomato", "cyan", "yellow", "thistle", "salmon2", "tan4"));
}

plotAnnotationLegend <- function(cex.titles = 1.0, cex.items = 0.8, ...) {
  annotations <- list(...)
  annotations <- annotations[which(lapply(annotations, function(x) {sum(!is.na(x))}) > 0)]
  if (length(annotations) > 0) {
    current.palette <- palette()
    palette(get_annotation_palette())
    par("mar" = c(0, 0, 0, 0))
    plot(NA, xlim = c(0, 1), ylim = c(0, 10), axes = F)

    num.rows = length(annotations);
    padding <- 0.005
    box_width <- 0.01
    annotation_width <- c()
    max.per.row <- c()
    for (a in 1:length(annotations)) {
      these.levels <- levels(annotations[[a]])[unique(annotations[[a]])]
      annotation_width[a] <- max(0.1,
                                 padding + box_width + padding + max(strwidth(these.levels, cex = cex.items)) + padding)
      max.per.row[a] <- max(1, floor(1 / annotation_width[a]))
      num.rows <- num.rows + 0.5 * ceiling(length(these.levels) / max.per.row[a])
    }
    par("usr" = c(0, 1, 0, num.rows + 1))
    par("usr" = c(0, 1,
                  0.8 - strheight("Hello", cex = cex.items) / 2,
                  num.rows + strheight("Hello", font = 2, cex = cex.titles) / 2))
    row <- num.rows
    # print(num.rows)
    for (a in 1:length(annotations)) {
      these.levels <- levels(annotations[[a]])[unique(annotations[[a]])]
      text(padding, row, adj = c(0, 0.5), names(annotations)[a], cex = cex.titles, font = 2)
      row <- row - 0.7
      for (b in 1:(ceiling(length(these.levels) / max.per.row[a]))) {
        # print(paste(a, b, row))
        from <- 1 + (b - 1) * max.per.row[a]
        to <- min(b * max.per.row, length(levels(annotations[[a]])))
        for (c in from:to) {
          level.value <- sort(unique(annotations[[a]]))[c]
          rect(padding + (c -from) * annotation_width[a],
               row - strheight("Legend", cex = cex.items) / 2,
               padding + box_width + (c -from) * annotation_width[a],
               row + strheight("Legend", cex = cex.items) / 2,
               col = get_annotation_palette()[level.value],
               border = NA)
          text(padding + box_width + padding + (c -from) * annotation_width[a], row, adj = c(0, 0.5),
               levels(annotations[[a]])[level.value], cex = cex.items)
        }
        row <- row - 0.5
      }
      row <- row - 0.3
    }
    palette(current.palette)
  }
}

process_tree <- function(tree, data, tree.width = 100) {
  tree.sorted.cells <- c()
  if (length(setdiff(tree$tip.label, colnames(data)))) {
    stop("The tree contains samples that do not exists in the data set")
  }
  tree.arcs <- cbind(tree$edge, tree$edge.length)
  colnames(tree.arcs) <- c("parent", "child", "dist")
  tree.arcs <- cbind(tree.arcs, is.leaf = FALSE, x1 = NA, x2 = 0, y1 = 0, y2 = 0)
  set.coords <- function(tree.arcs, node, depth, loc) {
    arcs <- which(tree.arcs[, 1] == node)
    if (length(arcs) == 0) {
      tree.arcs[which(tree.arcs[, 2] == node), "y2"] <- max(tree.arcs[, "y2"]) + 1
      tree.arcs[which(tree.arcs[, 2] == node), "is.leaf"] <- TRUE
      tree.sorted.cells <<- c(tree.sorted.cells, node)
    }
    for (arc in arcs) {
      new.depth <- depth + tree.arcs[arc, 3]
      tree.arcs[arc, "x1"] <- depth
      tree.arcs[arc, "x2"] <- new.depth
      tree.arcs <- set.coords(tree.arcs, tree.arcs[arc, 2], new.depth)
    }
    if (length(arcs) > 0) {
      y2 <- (min(tree.arcs[arcs, "y2"]) + max(tree.arcs[arcs, "y2"])) / 2
      tree.arcs[which(tree.arcs[, 2] == node), "y2"] <- y2
      tree.arcs[arcs, "y1"] <- y2
    }
    return(tree.arcs)
  }
  root <- setdiff(tree.arcs[, 1], intersect(tree.arcs[, 1], tree.arcs[, 2]))
  tree.arcs <- set.coords(tree.arcs, root, 0)
  
  max.x <- max(tree.arcs[, "x2"], na.rm = TRUE)
  max.y <- max(tree.arcs[, "y2"] + 1)
  # tree.width <- (max(x.right) - min(x.left)) / 10
  tree.arcs[, "x1"] <- tree.arcs[, "x1"] / max.x * tree.width
  tree.arcs[, "x2"] <- tree.arcs[, "x2"] / max.x * tree.width
  # segments(tree.arcs[, "x1"], tree.arcs[, "y1"], tree.arcs[, "x1"], tree.arcs[, "y2"])
  # segments(tree.arcs[, "x1"], tree.arcs[, "y2"], tree.arcs[, "x2"], tree.arcs[, "y2"], col = "red", lwd = 5)

  return(list(tree.arcs, tree.sorted.cells))  
}

plotGenome <- function(data, showBins = FALSE, showTics = TRUE, tree, genes, ...) {
  chr.lengths <- get.chr.lengths(data)
  
  ## x.left and x.right are the start and end coordinates of each bin on the plot
  x.left <- data$START
  x.right <- data$END
  x.left <- x.left + chr.lengths[as.character(data$CHR), "from"]
  x.right <- x.right + chr.lengths[as.character(data$CHR), "from"]
  
  ## ================================================================
  ## Draw a tree for the samples if provided
  ## ================================================================
  tree.arcs <- NA
  tree.sorted.cells <- c()
  tree.width <- 0
  if (!missing(tree)) {
    tree.width <- (max(x.right) - min(x.left)) / 10
    tree.list <- process_tree(tree = tree, data = data, tree.width = tree.width)
    tree.arcs <- tree.list[[1]]
    tree.sorted.cells <- tree.list[[2]]
    data <- data[, c(1:3, tree.sorted.cells + 3)]
  }
  ## ================================================================
  
  samples <- colnames(data)[-c(1:3)]
  n.samples <- length(samples)
  
  par("mar" = c(0, 5, 3, 0) + 0.1)
  plot(NA, xlim = c(min(x.left), max(x.right)), ylim = c(0, n.samples),
       axes = F, xlab = NA, ylab = NA)
  x.margin <-  (max(x.right) - min(x.left)) / 500
  y.margin <-  strheight("legend") * 1.1
  annotations <- list(...)
  annotation.width <- strwidth("m")
  
  usr <- c(min(x.left) - x.margin, max(x.right) + x.margin, -y.margin, n.samples + y.margin * 1.5)
  if (length(annotations) > 0) {
    usr <- c(min(x.left) - x.margin, max(x.right) + x.margin + annotation.width * length(annotations), -y.margin, n.samples + y.margin * 1.5)
  }
  if (!is.null(dim(tree.arcs))) {
    tree.end <- usr[2] + x.margin + tree.width
    usr[2] <- tree.end
  }
  par("usr" = usr)
  
  if (length(annotations) > 0) {
    if (missing(tree)) {
      new.order <- order(...)
      data <- data[, c(1:3, new.order+3)]
    } else {
      new.order <- tree.sorted.cells
    }
    # usr <- c(min(x.left) - x.margin, max(x.right) + x.margin + annotation.width * length(annotations), -y.margin, n.samples + y.margin * 1.5)
    # par("usr" = usr)
    current.palette <- palette()
    palette(get_annotation_palette())
    for (a in 1:length(annotations)) {
      annotations[[a]] <- annotations[[a]][new.order]
      rect(max(x.right) + annotation.width * (a - 1) + annotation.width * 0.2, 0:(n.samples - 1) + 0.1,
           max(x.right) + annotation.width * a, 0:(n.samples - 1) + 0.9,
           col = annotations[[a]], border = NA)
      # text(max(x.right) + annotation.width * (a - 1), ncol(data) + 1.2,
      #      labels = names(annotations)[a], cex = 0.7, srt = 90, pos =4)
      mtext(names(annotations)[a], at = max(x.right) + annotation.width * (a - 0.1),
            line = -0.9,
            adj = 0, padj = 0,
            cex = 0.8, las = 2)
    }
    palette(current.palette)
    # } else {
    #   usr <- c(min(x.left) - x.margin, max(x.right) + x.margin, -y.margin, n.samples + y.margin * 1.5)
    #   par("usr" = usr)
  }
  if (!is.null(dim(tree.arcs))) {
    # usr <- par("usr")
    # tree.end <- usr[2] + x.margin + tree.width
    # usr[2] <- tree.end
    # par("usr" = usr)
    segments(tree.end - tree.arcs[, "x1"], tree.arcs[, "y1"] - 0.5, tree.end - tree.arcs[, "x1"], tree.arcs[, "y2"] - 0.5)
    segments(tree.end - tree.arcs[, "x1"], tree.arcs[, "y2"] - 0.5, tree.end - tree.arcs[, "x2"], tree.arcs[, "y2"] - 0.5)
  }
  
  rect(chr.lengths[, "from"], 0, chr.lengths[, "to"], n.samples, border = NA, col = c("white", rgb(0.9, 0.9, 0.9))[2 - (1:nrow(chr.lengths)) %% 2])
  
  if (showBins) {
    segments(x.right, 0, x.right, n.samples, col = rgb(0.7, 0.7, 0.7), lwd = min(25 / length(x.right), 0.1))
  }
  if (showTics) {
    chr <- unique(as.character(data$CHR))
    if (length(chr) == 1) {
      if (chr.lengths[chr, "length"] > 10000000) {
        tics <- seq(10000000, chr.lengths[chr, "length"], by = 10000000)
        tics <- tics + chr.lengths[chr, "from"]
        segments(tics, -y.margin / 4, tics, 0, col = "black", lwd = 1)
      }
    }
  }
  
  for (cell in 0:(n.samples - 1)) {
    d <- data[, cell + 4]
    draw.rect <- function(y, col) {
      r <- range(y); rect(x.left[r[1]], cell + 0.05, x.right[r[2]], cell + 0.95, col = col, border = NA)
    }
    x <- which(d == 0)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0, 0, 0.8))
    }
    x <- which(d == 1)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0.6, 0.6, 1))
    }
    x <- which(d == 3)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(1, 0.6, 0.6))
    }
    x <- which(d == 4)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0.8, 0, 0))
    }
    x <- which(d == 5)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0.8, 0.8, 0))
    }
    x <- which(d == 6)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0, 0.8, 0.8))
    }
    x <- which(d > 6)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), draw.rect, col = rgb(0.8, 0, 0.8))
    }
  }
  if (nrow(chr.lengths) > 1) {
    lines(as.vector(t(cbind(chr.lengths[-c(1), "from"], chr.lengths[-c(1), "from"], NA))),
          rep(c(0, n.samples, NA), times = nrow(chr.lengths) - 1),
          col = "grey", lwd = 1.5)
  }
  rect(min(chr.lengths[, "from"]), 0, max(chr.lengths[, "to"]), n.samples, border = "black", col = NA)
  text((chr.lengths[, "from"] + chr.lengths[, "to"]) / 2, -y.margin / 2, labels = sub("chr", "", row.names(chr.lengths)), cex = 0.7)
  axis(2, at = 0:(n.samples - 1) + 0.5, labels = colnames(data[, -c(1:3), drop = F]), las = 2, tick = FALSE, cex.axis = 0.5, gap.axis = -1)
  
  x.usr.length = usr[2] - usr[1]
  box.height <- y.margin
  box.width <- box.height * par("pin")[2] / par("pin")[1] * (usr[2] - usr[1]) / (usr[4] - usr[3])
  unit.width <- box.width * 1.5 + strwidth("n=4+")
  separation.width <- box.width * 1.5
  total.width <- 4 * unit.width + 3 * separation.width
  x.usr.start <- usr[1] + (x.usr.length - total.width) / 2
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0, 0, 0.8))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=0", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.6, 0.6, 1))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=1", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(1, 0.6, 0.6))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=3", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.8, 0, 0))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=4+", adj = c(0, 1))
  # legend("top", ncol = 4, bty = "n",
  #        legend = c("n=0", "n=1", "n=3", "n=4+"),
  #        fill = c(rgb(0, 0, 0.8), rgb(0.6, 0.6, 1), rgb(1, 0.6, 0.6), rgb(0.8, 0, 0)))

  if (!missing(genes)) {
    addGeneLine <- function(gene.symbol, chr, location, col = rgb(0, 0, 0, 0.5), to) {
      location <- location + chr.lengths[chr, "from"]
      if (missing(to)) {
        segments(location, 0, location, n.samples, col = col)
      } else {
        to <- to + chr.lengths[chr, "from"]
        rect(location, 0, to, n.samples, border = col, col = adjustcolor(col, alpha.f = 0.2))
      }
      text(location, n.samples, labels = gene.symbol, srt = 90, adj = c(1, 0), col = col)
    }
    
    for (this_gene in names(genes)) {
      if (length(genes[[this_gene]]) == 4) {
        addGeneLine(this_gene, chr = genes[[this_gene]][1], location = as.numeric(genes[[this_gene]][2]),
                    to = as.numeric(genes[[this_gene]][3]), col = genes[[this_gene]][4])
      } else if (length(genes[[this_gene]]) == 3) {
          addGeneLine(this_gene, genes[[this_gene]][1], as.numeric(genes[[this_gene]][2]), genes[[this_gene]][3])
      } else {
        addGeneLine(this_gene, genes[[this_gene]][1], as.numeric(genes[[this_gene]][2]))
      }
    }    
  }

  invisible(data)
}


plotNormGenome <- function(data, showBins = FALSE, showTics = TRUE, tree, genes,
                           cell_order,...) {
  chr.lengths <- get.chr.lengths(data)
  
  ## x.left and x.right are the start and end coordinates of each bin on the plot
  x.left <- data$START
  x.right <- data$END
  x.left <- x.left + chr.lengths[as.character(data$CHR), "from"]
  x.right <- x.right + chr.lengths[as.character(data$CHR), "from"]
  
  ## ================================================================
  ## Draw a tree for the samples if provided
  ## ================================================================
  tree.arcs <- NA
  tree.sorted.cells <- c()
  tree.width <- 0
  if (!missing(tree)) {
    tree.width <- (max(x.right) - min(x.left)) / 10
    tree.list <- process_tree(tree = tree, data = data, tree.width = tree.width)
    tree.arcs <- tree.list[[1]]
    tree.sorted.cells <- tree.list[[2]]
    data <- data[, c(1:3, tree.sorted.cells + 3)]
  } else if (!missing(cell_order)) {
    data <- data[, c(1:3, cell_order + 3)]
  }
  ## ================================================================
  
  samples <- colnames(data)[-c(1:3)]
  n.samples <- length(samples)
  
  par("mar" = c(0, 5, 2, 0) + 0.1)
  plot(NA, xlim = c(min(x.left), max(x.right)), ylim = c(0, n.samples),
       axes = F, xlab = NA, ylab = NA)
  x.margin <-  (max(x.right) - min(x.left)) / 500
  y.margin <-  strheight("legend") * 1.1
  annotations <- list(...)
  annotation.width <- strwidth("m")
  
  usr <- c(min(x.left) - x.margin, max(x.right) + x.margin, -y.margin, n.samples + y.margin * 1.5)
  if (length(annotations) > 0) {
    usr <- c(min(x.left) - x.margin, max(x.right) + x.margin + annotation.width * length(annotations), -y.margin, n.samples + y.margin * 1.5)
  }
  if (!is.null(dim(tree.arcs))) {
    tree.end <- usr[2] + x.margin + tree.width
    usr[2] <- tree.end
  }
  par("usr" = usr)
  
  if (length(annotations) > 0) {
    if (!missing(tree)) {
      new.order <- tree.sorted.cells
    } else if (!missing(cell_order)) {
      new.order <- cell_order
    } else {
      new.order <- order(...)
      data <- data[, c(1:3, new.order+3)]
    }
    # usr <- c(min(x.left) - x.margin, max(x.right) + x.margin + annotation.width * length(annotations), -y.margin, n.samples + y.margin * 1.5)
    # par("usr" = usr)
    current.palette <- palette()
    palette(get_annotation_palette())
    for (a in 1:length(annotations)) {
      annotations[[a]] <- annotations[[a]][new.order]
      rect(max(x.right) + annotation.width * (a - 1) + annotation.width * 0.2, 0:(n.samples - 1) + 0.1,
           max(x.right) + annotation.width * a, 0:(n.samples - 1) + 0.9,
           col = annotations[[a]], border = NA)
      # text(max(x.right) + annotation.width * (a - 1), ncol(data) + 1.2,
      #      labels = names(annotations)[a], cex = 0.7, srt = 90, pos =4)
      mtext(names(annotations)[a], at = max(x.right) + annotation.width * (a - 0.1),
            line = -0.9,
            adj = 0, padj = 0,
            cex = 0.8, las = 2)
    }
    palette(current.palette)
    # } else {
    #   usr <- c(min(x.left) - x.margin, max(x.right) + x.margin, -y.margin, n.samples + y.margin * 1.5)
    #   par("usr" = usr)
  }
  if (!is.null(dim(tree.arcs))) {
    # usr <- par("usr")
    # tree.end <- usr[2] + x.margin + tree.width
    # usr[2] <- tree.end
    # par("usr" = usr)
    segments(tree.end - tree.arcs[, "x1"], tree.arcs[, "y1"] - 0.5, tree.end - tree.arcs[, "x1"], tree.arcs[, "y2"] - 0.5)
    segments(tree.end - tree.arcs[, "x1"], tree.arcs[, "y2"] - 0.5, tree.end - tree.arcs[, "x2"], tree.arcs[, "y2"] - 0.5)
  }
  
  rect(chr.lengths[, "from"], 0, chr.lengths[, "to"], n.samples, border = NA, col = c("white", rgb(0.9, 0.9, 0.9))[2 - (1:nrow(chr.lengths)) %% 2])
  
  if (showBins) {
    segments(x.right, 0, x.right, n.samples, col = rgb(0.7, 0.7, 0.7), lwd = min(25 / length(x.right), 0.1))
  }
  if (showTics) {
    chr <- unique(as.character(data$CHR))
    if (length(chr) == 1) {
      if (chr.lengths[chr, "length"] > 10000000) {
        tics <- seq(10000000, chr.lengths[chr, "length"], by = 10000000)
        tics <- tics + chr.lengths[chr, "from"]
        segments(tics, -y.margin / 4, tics, 0, col = "black", lwd = 1)
      }
    }
  }

  cr <- colorRamp(c(rgb(0, 0, 0.8, 1), rgb(0.6, 0.6, 1, 1), rgb(1, 1, 1, 0), rgb(1, 0.6, 0.6, 1), rgb(0.8, 0, 0, 1)), alpha = T)
  for (cell in 0:(n.samples - 1)) {
    d <- data[, cell + 4]
    d <- d / 2
    d[d > 1] <- 1
    d[d < 0] <- 0
    col <- cr(d)
    col <- rgb(col[, 1], col[, 2], col[, 3], col[, 4], maxColorValue = 255)
    rect(x.left, cell + 0.15, x.right, cell + 0.95, col = col, border = col)
  }
  if (nrow(chr.lengths) > 1) {
    lines(as.vector(t(cbind(chr.lengths[-c(1), "from"], chr.lengths[-c(1), "from"], NA))),
          rep(c(0, n.samples, NA), times = nrow(chr.lengths) - 1),
          col = "grey", lwd = 1.5)
  }
  rect(min(chr.lengths[, "from"]), 0, max(chr.lengths[, "to"]), n.samples, border = "black", col = NA)
  text((chr.lengths[, "from"] + chr.lengths[, "to"]) / 2, -y.margin / 2, labels = sub("chr", "", row.names(chr.lengths)), cex = 0.7)
  axis(2, at = 0:(n.samples - 1) + 0.5, labels = colnames(data[, -c(1:3), drop = F]), las = 2, tick = FALSE, cex.axis = 0.5)
  
  x.usr.length = usr[2] - usr[1]
  box.height <- y.margin
  box.width <- box.height * par("pin")[2] / par("pin")[1] * (usr[2] - usr[1]) / (usr[4] - usr[3])
  unit.width <- box.width * 1.5 + strwidth("n=4+")
  separation.width <- box.width * 1.5
  total.width <- 4 * unit.width + 3 * separation.width
  x.usr.start <- usr[1] + (x.usr.length - total.width) / 2
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0, 0, 0.8))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=0", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.6, 0.6, 1))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=1", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(1, 0.6, 0.6))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=3", adj = c(0, 1))
  x.usr.start <- x.usr.start + unit.width + separation.width
  rect(x.usr.start, usr[4] - y.margin, x.usr.start + box.width, usr[4], border = "black", col = rgb(0.8, 0, 0))
  text(x.usr.start + box.width * 1.5, usr[4], labels = "n=4+", adj = c(0, 1))
  # legend("top", ncol = 4, bty = "n",
  #        legend = c("n=0", "n=1", "n=3", "n=4+"),
  #        fill = c(rgb(0, 0, 0.8), rgb(0.6, 0.6, 1), rgb(1, 0.6, 0.6), rgb(0.8, 0, 0)))
  
  if (!missing(genes)) {
    addGeneLine <- function(gene.symbol, chr, location, col = rgb(0, 0, 0, 0.5)) {
      location <- location + chr.lengths[chr, "from"]
      segments(location, 0, location, n.samples, col = col)
      text(location, n.samples, labels = gene.symbol, srt = 90, adj = c(1, 0), col = col)
    }
    
    for (this_gene in names(genes)) {
      if (length(genes[[this_gene]]) == 3) {
        addGeneLine(this_gene, genes[[this_gene]][1], as.numeric(genes[[this_gene]][2]), genes[[this_gene]][3])
      } else {
        addGeneLine(this_gene, genes[[this_gene]][1], as.numeric(genes[[this_gene]][2]))
      }
    }    
  }
  
  invisible(data)
}


plotChr3 <- function(data, chr, showBins = T, ...) {
  data <- subset(data, CHR == chr)
  if (data[1, "START"] > 1) {
    data <- rbind(data[1, , drop = FALSE], data)
    data[1, "START"] <- 1
    data[1, "END"] <- data[2, "START"] - 1
    data[1, -c(1:3)] <- 2
  }
  data$CHR <- factor(data$CHR)
  chr.lengths <- get.chr.lengths(data)

  data <- plotGenome(data, showBins, ...)
  n.samples <- ncol(data) - 3
  
  addGeneLine <- function(gene.symbol, chr, location, col = rgb(0, 0, 0, 0.5)) {
    location <- location + chr.lengths[chr, "from"]
    segments(location, 0, location, n.samples, col = col)
    text(location, n.samples, labels = gene.symbol, srt = 90, adj = c(1, 0), col = col)
  }
  
  if (chr == "chr9") {
    addGeneLine("CDKN2A", "chr9", (21967751+21995300)/2)
    addGeneLine("PAX5", "chr9", (36833272+37034103)/2)
  }
  if (chr == "chr12") {
    addGeneLine("ETV6", "chr12", (11802788+12048336)/2)
  }
  if (chr == "chr21") {
    addGeneLine("RUNX1", "chr21", (36160098+37376965)/2)
  }
  
}


plotNormChr <- function(data, chr, showBins = T, ...) {
  data <- subset(data, CHR == chr)
  data$CHR <- factor(data$CHR)
  chr.lengths <- get.chr.lengths(data)
  
  data <- plotNormGenome(data, showBins, ...)
  n.samples <- ncol(data) - 3
  
  addGeneLine <- function(gene.symbol, chr, location, col = rgb(0, 0, 0, 0.5)) {
    location <- location + chr.lengths[chr, "from"]
    segments(location, 0, location, n.samples, col = col)
    text(location, n.samples, labels = gene.symbol, srt = 90, adj = c(1, 0), col = col)
  }
  
  if (chr == "chr9") {
    addGeneLine("CDKN2A", "chr9", (21967751+21995300)/2)
    addGeneLine("PAX5", "chr9", (36833272+37034103)/2)
  }
  
}


plotChr2 <- plotChr3

plotChr <- plotChr3

# plotChr <- function(data, chr) {
#   chr.lengths <- get.chr.lengths(data)
#   if (is.na(chr.lengths[chr, 1])) {
#     stop(paste("Chromosome", chr, "is not found in the data frame"))
#   }
# 
#   plot(NA,
#        xlim = c(0, chr.lengths[chr, "length"]),
#        ylim = c(0, 5),
#        axes = FALSE, xlab = NA, ylab = NA)
#   mtext(text = "Copy number", side = 2, line = 1, cex = 1)
#   mtext(text = paste("Chromosome", chr), side = 3, line = 0.6, cex = 1)
#   for (l in 0:5) {
#     mtext(text = l, side = 2, line = -1, at = l, cex = 1, las = 1)
#   }
#   
#   n.samples = ncol(data) - 3
#   rect(0,
#        0,
#        chr.lengths[chr, "length"],
#        5,
#        col = rgb(0.85,0.85,0.85),
#        border = NA)
#   
#   for (l in 0:5) {
#     lines(c(0, chr.lengths[chr, "length"]),
#           c(l, l),
#           col = "black")
#   }
# 
#   ## Add annotation for key genes (a dotted line with the gene name at the top)
#   draw.gene("P16", (21967751+21995300)/2)
#   draw.gene("PAX5", (36833272+37034103)/2)
#   # draw.gene("ETV6", (11802788+12048336)/2 + max.chrpos[12]) 
#   # draw.gene("RUNX1", (36160098+37376965)/2 + max.chrpos[21]) 
# 
#   # Change NAs to -10 (just because diff() would carry NAs forever)
#   data.copy <- data
#   data.copy[is.na(data.copy)] <- -10
# 
#   ## Plot the CNV for each cell
#   for (cell in 1:n.samples) {
#     # Split the data in lists of consecutive identical values
#     # * diff(data.result[x.chr == x.chr.level, cell]) != 0 is true ( = 1) if value changes
#     # * c(1, diff(...)) is to initialise the vector with a 1
#     # * cumsum() changes every time diff is 1, i.e. every time the value changes
#     # * split return a list of values as defined by the categories in cumsum, i.e. a list
#     #   of consecutive identical values.
#     # Each element of the list is a vector of consecutive identical values.
#     list.data <- split(data.copy[data.copy$CHR == chr, cell + 3],
#                        cumsum(c(1, diff(data.copy[data.copy$CHR == chr, cell + 3]) != 0)))
#     
#     bin_num <- 1 # bin counter
#     for (list in list.data) {
#       length <- length(list) # Number of bins with consecutive value
#       value <- list[1] # Actual value of CN
#       if (value != -10) {
#         lines(c(subset(data.copy, CHR == chr)[bin_num, "START"],
#                 subset(data.copy, CHR == chr)[bin_num + length - 1, "END"]),
#               c(value - (ncol(data.result) / 2 - cell) / 35,
#                 value - (ncol(data.result) / 2 - cell) / 35),
#               col = rainbow(n.samples)[cell], lwd = 3, lend = 1)
#       }
#       bin_num <- bin_num + length # Increase bin counter
#     }
#   }
# }
# 
# 
# 
# plotChr2 <- function(data, chr = "chr9") {
#   chr.lengths <- get.chr.lengths(data)
#   data <- subset(data, CHR == chr)
#   n.samples = ncol(data) - 3
#   
#   x.left <- data$START
#   x.right <- data$END
#   
#   x.left <- x.left + chr.lengths[as.character(data$CHR), "from"]
#   x.right <- x.right + chr.lengths[as.character(data$CHR), "from"]
#   
#   plot(NA, xlim = c(chr.lengths[chr, "from"], chr.lengths[chr, "to"]), ylim = c(0, n.samples + 10),
#        axes = F, xlab = NA, ylab = NA)
#   x.margin <-  10
#   y.margin <-  10
#   x.margin <-  (max(x.right) - min(x.left)) / 500
#   y.margin <-  strheight("hello") * 3
#   par("usr" = c(chr.lengths[chr, "from"] - x.margin, chr.lengths[chr, "to"] + x.margin,
#                 - y.margin / 10, n.samples + y.margin))
#   palette(c(rgb(0.9, 0.9, 0.9), palette()[2:length(palette())]))
#   rect(chr.lengths[chr, "from"], 0, chr.lengths[chr, "to"], n.samples, border = NA, col = 1)
#   segments(x.right, 0, x.right, n.samples, col = "black", lwd = 0.1)
#   
#   for (cell in 1:n.samples) {
#     cell.col <- cell + 3
#     hom.deletions <- data.frame(NA, ncol = 2, nrow = 0)
#     het.deletions <- data.frame(NA, ncol = 2, nrow = 0)
#     gains <- data.frame(NA, ncol = 2, nrow = 0)
#     amplifications <- data.frame(NA, ncol = 2, nrow = 0)
#     bin = 1;
#     while (bin <= nrow(data)) {
#       from_bin <- bin
#       while (bin < nrow(data) & data[bin, cell.col] == data[bin + 1, cell.col]) {
#         bin <- bin + 1
#       }
#       to_bin <- bin
#       if (data[bin, cell.col] < 0.5) {
#         hom.deletions <- rbind(hom.deletions, c(x.left[from_bin], x.right[to_bin]))
#         # rect(x.left[from_bin], cell + 0.1, x.right[to_bin], cell + 0.9,
#         #      col = rgb(0, 0, 0.8), border = NA)
#       } else if (data[bin, cell.col] < 1.5) {
#         het.deletions <- rbind(het.deletions, c(x.left[from_bin], x.right[to_bin]))
#         # rect(x.left[from_bin], cell + 0.1, x.right[to_bin], cell + 0.9,
#         #      col = rgb(0.6, 0.6, 1), border = NA)
#       } else if (data[bin, cell.col] > 3.5) {
#         amplifications <- rbind(amplifications, c(x.left[from_bin], x.right[to_bin]))
#         # rect(x.left[from_bin], cell + 0.1, x.right[to_bin], cell + 0.9,
#         #      col = rgb(0.8, 0, 0), border = NA)
#       } else if (data[bin, cell.col] > 2.5) {
#         gains <- rbind(gains, c(x.left[from_bin], x.right[to_bin]))
#         # rect(x.left[from_bin], cell + 0.1, x.right[to_bin], cell + 0.9,
#         #      col = rgb(1, 0.6, 0.6), border = NA)
#       }
#       bin <- bin + 1
#     }
#     rect(gains[, 1], cell - 0.9, gains[, 2], cell - 0.1,
#          col = rgb(1, 0.6, 0.6), border = NA)
#     rect(het.deletions[, 1], cell - 0.9, het.deletions[, 2], cell - 0.1,
#          col = rgb(0.6, 0.6, 1), border = NA)
#     rect(amplifications[, 1], cell - 0.9, amplifications[, 2], cell - 0.1,
#          col = rgb(0.8, 0, 0), border = NA)
#     rect(hom.deletions[, 1], cell - 0.9, hom.deletions[, 2], cell - 0.1,
#          col = rgb(0, 0, 0.8), border = NA)
#     # rect(chr.lengths[, "from"], cell + 0.1, chr.lengths[, "to"], cell + 0.9, border = "black", col = "NA")
#   }
#   # Chromosome separation lines
#   segments(chr.lengths[-c(1), "from"], 0, chr.lengths[-c(1), "from"], n.samples,
#         col = "grey", lwd = 1.5)
#   # Black outline of the plot
#   rect(min(chr.lengths[chr, "from"]), 0, max(chr.lengths[chr, "to"]), n.samples, border = "black", col = NA)
#   text((chr.lengths[chr, "from"] + chr.lengths[chr, "to"]) / 2, -y.margin / 2, labels = sub("chr", "", chr), cex = 0.7)
#   axis(1, at = (chr.lengths[, "from"] + chr.lengths[, "to"]) / 2,
#        labels = sub("chr", "", rownames(chr.lengths)), las = 1, tick = FALSE)
#   axis(2, at = 1:n.samples - 0.5, labels = colnames(data[, -c(1:3)]), las = 2, tick = FALSE)
#   drawGene <- function(gene.symbol, chr, location, col = rgb(0, 0, 0, 0.5)) {
#     location <- location + chr.lengths[chr, "from"]
#     segments(location, 0, location, n.samples, col = col)
#     text(location, n.samples, labels = gene.symbol, srt = 90, adj = c(1, 0), col = col)
#   }
#   # draw.gene("P16", (21967751+21995300)/2)
#   # draw.gene("PAX5", (36833272+37034103)/2)
# 
#   legend("top", ncol = 4, bty = "n",
#          legend = c("n=0", "n=1", "n=3", "n=4+"),
#          fill = c(rgb(0, 0, 0.8), rgb(0.6, 0.6, 1), rgb(1, 0.6, 0.6), rgb(0.8, 0, 0)))
# }


normalise_sizes <- function(these_sizes, original_sizes = these_sizes, norm = "perc") {
  norm.sizes <- these_sizes
  if (norm == "sample") {
    for (i in 1:num.sizes) {
      max.size <- min(max(these_sizes[, i]), 1)
      norm.sizes[, i] <- 0.4 * ((max.size * 0.1 + (these_sizes[, i])) ** 0.5 / (1.1 * max.size) ** 0.5)
    }
  } else if (norm == "perc") {
    max.perc <- max(apply(original_sizes, 2, function(x) { if (sum(x) > 0) { return(x / sum(x)) } else { return(rep(0.1, length(x))) } }))
    for (i in 1:ncol(these_sizes)) {
      max.size <- sum(these_sizes[, i]) * max.perc
      norm.sizes[, i] <- 0.4 * ((max.size * 0.1 + (these_sizes[, i])) ** 0.5 / (1.1 * max.size) ** 0.5)
    }
  } else {
    norm.sizes <- 0.4 * ((max(these_sizes) * 0.1 + (these_sizes)) ** 0.5 / (1.1 * max(these_sizes)) ** 0.5)
    rownames(norm.sizes) <- row.names(these_sizes)
  }
  return(norm.sizes)
}

#' plot.node
#'
#' @param x X-coordinate of the center of the node 
#' @param y Y-coordinate of the center of the node
#' @param sizes a vector of radiuses for the node. If displaying multiple data, each radius will be used for each segment of the node
#' @param counts a vector of actual numbers. This is used to determine empty sets (shown with a white background) and to display them
#' with the plot.counts option.
#' @param plot.counts logical to indicate whether counts should be plotted on the segments
#' @param plot.internal.lines logical to indicate whether to draw or not the internal separation between the segments.
#'
#' @return Nothing
#' @export
#'
#' @examples
plot.node <- function(x, y, sizes, counts, plot.internal.lines = T,
                      plot.counts = F,
                      plot.counts.colour = "grey50",
                      plot.counts.cex = 0.8,
                      node.color = NA) {
  num.sizes <- length(sizes)
  
  theta.margins <- seq(0, 2 * pi, length.out = num.sizes + 1)
  theta <- matrix(NA, nrow = num.sizes, ncol = floor(200 / num.sizes))
  for (i in 1:num.sizes) {
    theta[i, ] <- seq(theta.margins[i], theta.margins[i + 1], length.out = floor(200 / num.sizes))
  }
  
  for (i in 1:num.sizes) {
    size <- as.numeric(sizes[i])
    # tree.nodes[node, "size"] <<- max(tree.nodes[node, "size"], size, na.rm = T)
    x.lines <- x - size * sin(theta[i, ])
    y.lines <- y - size * cos(theta[i, ])
    if (num.sizes > 2) {
      x.polygon <- c(x, x.lines, x)
      y.polygon <- c(y, y.lines, y)
    } else {
      x.polygon <- x.lines
      y.polygon <- y.lines
    }
    if (is.na(node.color)) {
      col <- i
    } else {
      col <- node.color
    }
    if (plot.internal.lines) {
      polygon(x = x.polygon, y = y.polygon, col = ifelse(counts[i] == 0, "white", col), lwd = 1)
    } else {
      polygon(x = x.polygon, y = y.polygon, col = ifelse(counts[i] == 0, "white", col), lty = "blank")
      lines(x = x.lines, y = y.lines, col = "black", lwd = 1)
    }
    if (plot.counts) {
      text(x - size / 2 * sin(theta[i, ncol(theta) / 2]),
           y - size / 2  * cos(theta[i, ncol(theta) / 2]),
           counts[i], cex = plot.counts.cex,
           col = plot.counts.colour)
    }
  }
}

plot.nodes <- function(tree.nodes, counts, plot.counts, node.colors) {
  num.counts <- ncol(counts)
  # palette(c("lightblue", "pink", "yellow3"))
  
  norm.sizes <- normalise_sizes(counts)
  # print(cbind(counts, norm.sizes))

  if (!missing(node.colors)) {
    for (node in rownames(tree.nodes)) {
      plot.node(tree.nodes[node, "x"], tree.nodes[node, "y"],
                as.vector(norm.sizes[node, ]), counts[node, ],
                plot.counts = plot.counts,
                node.color = node.colors[as.numeric(node)])
    }
  } else {
    for (node in rownames(tree.nodes)) {
      plot.node(tree.nodes[node, "x"], tree.nodes[node, "y"],
                as.vector(norm.sizes[node, ]), counts[node, ],
                plot.counts = plot.counts)
    }
  } 
  
} 


#' plot.tree
#'
#' @param tree a data.frame with first column being the parent node ID and the second column the
#' child node ID, one row per connection/arc/relationship. An optional third column can contain
#' the distance between the nodes.
#' @param sizes a data.frame with rownames corresponding to the node IDs and the number of
#' cells in a column. You can have more than one column if you want to represent several samples,
#' like counts of cells before and after treatment for instance.
#' @param root
#' @param square.layout 
#' @param direction 
#' @param main 
#'
#' @return
#' @export
#'
#' @examples
plot.tree <- function(tree, counts, root = 1, direction = "down", main, norm = "perc", lab.cex = 0.5, plot.counts = F,
                      palette = c("lightblue", "pink", "yellow3", "grey"),
                      legend, legend.perc.cex = 0.8, legend.samples.cex = 1) {

  use.colors.for.nodes <- F
  if (!is.null(names(palette)) & ncol(tree) >= 4) {
    use.colors.for.nodes <- T
  }
  palette(palette)
  node.names <- sort(unique(c(tree[, 1], tree[, 2])))

  if (missing(counts)) {
    counts = data.frame(S1 = rep(1, length(node.names)))
    rownames(counts) = node.names
  }

  # Set the label to an empty string if not defined in the data.frame (3rd column)
  if (ncol(tree) == 2) {
    tree <- cbind(tree, "")
  }
  colnames(tree) <- c("parent", "child", "label")

  tree.arcs <- cbind(tree, is.leaf = FALSE, x1 = 0, x2 = 0, y1 = 0, y2 = -1)
  
  # print(tree.arcs)
  # print(sizes)
  set.coords <- function(node, depth = 0) {
    arcs <- which(tree[, 1] == node)
    if (length(arcs) == 0) {
      y2 <- max(tree.arcs[, "y2"]) + 1
      tree.arcs[which(tree.arcs[, 2] == node), "y2"] <<- y2
      tree.arcs[which(tree.arcs[, 2] == node), "is.leaf"] <<- TRUE
    }
    for (arc in arcs) {
      new.depth <- depth + 1
      tree.arcs[arc, "x1"] <<- depth
      tree.arcs[arc, "x2"] <<- new.depth
      # tree.arcs <<- set.coords(tree[arc, 2], new.depth)
      set.coords(tree[arc, 2], new.depth)
    }
    if (length(arcs) > 0) {
      y2 <- (min(tree.arcs[arcs, "y2"]) + max(tree.arcs[arcs, "y2"])) / 2
      tree.arcs[which(tree.arcs[, 2] == node), "y2"] <<- y2
      tree.arcs[arcs, "y1"] <<- y2
    }
    return(tree.arcs)
  }
  # tree.arcs <- set.coords(root, 0)
  set.coords(root)
  # print(tree.arcs)

  x1 <- tree.arcs[, "x1"]
  x2 <- tree.arcs[, "x2"]
  y1 <- tree.arcs[, "y1"]
  y2 <- tree.arcs[, "y2"]
  if (direction == "down") {
    tree.arcs[, "x1"] <- y1
    tree.arcs[, "x2"] <- y2
    tree.arcs[, "y1"] <- max(x2) - x1
    tree.arcs[, "y2"] <- max(x2) - x2
  } else if (direction == "up") {
    tree.arcs[, "x1"] <- y1
    tree.arcs[, "x2"] <- y2
    tree.arcs[, "y1"] <- x1
    tree.arcs[, "y2"] <- x2
  } else if (direction == "left") {
    tree.arcs[, "x1"] <- max(x2) - x1
    tree.arcs[, "x2"] <- max(x2) - x2
  }

  tree.nodes <- tree.arcs[, c("x2", "y2")]
  tree.nodes <- rbind(c(tree.arcs[which(tree.arcs[, 1] == root)[1], "x1"],
                        tree.arcs[which(tree.arcs[, 1] == root)[1], "y1"]),
                      tree.nodes)
  rownames(tree.nodes) = c(root, tree.arcs[, "child"])
  colnames(tree.nodes) = c("x", "y")
  
  # print(tree.nodes)

  # palette(rainbow(max(tree.arcs[, 2])))
  prev.mar <- par("mar")
  par("mar" = c(0, 0, 0, 0))
  max.x <- max(tree.nodes[, "x"], na.rm = TRUE)
  if (!missing(legend) & max.x < 3) {
    max.x <- max.x + 1
  }
  max.y <- max(tree.nodes[, "y"])
  # print(tree.nodes)
  plot(NA, xlim = c(-0.5, max.x + 0.5), ylim = c(-0.5, max.y + 0.5), axes = FALSE, xlab = NA, ylab = NA, asp = T)
  if (!missing(main)) {
    title(main = main, line = -1)
  }
  segments(tree.arcs[, "x1"], tree.arcs[, "y1"], tree.arcs[, "x2"], tree.arcs[, "y2"])

  if (use.colors.for.nodes) {
    plot.nodes(tree.nodes, counts, plot.counts, node.colors = palette)
  } else {
    plot.nodes(tree.nodes, counts, plot.counts)
  }

  ## ----------------------------------------------------------------
  ## Print labels on the branches
  ## ----------------------------------------------------------------
  apply(X = tree.arcs, MARGIN = 1, function(row) {
    # print(row)
    if (row[3] != "") {
      x.diff <- as.numeric(row["x2"]) - as.numeric(row["x1"]) 
      y.diff <- as.numeric(row["y2"]) - as.numeric(row["y1"])
      x <- (as.numeric(row["x1"]) + as.numeric(row["x2"])) / 2
      y <- (as.numeric(row["y1"]) + as.numeric(row["y2"])) / 2
      # x.width <- strwidth(paste0(" ", row[3], " "), cex = lab.cex) / 2
      # y.height <- strheight(row[3], cex = lab.cex) / 2 * 2 
      # rect(x - x.width, y - y.height, x + x.width, y + y.height, col = "white", border = "black")
      if (abs(x.diff) > 0) {
        text(x, y, row[3], cex = lab.cex, srt = atan(y.diff / x.diff) / pi * 180,
             adj = c(0.5, -0.3), font = 2)
      } else {
        text(x, y, row[3], cex = lab.cex, font = 2)
      }
    }
    })
  ## ----------------------------------------------------------------


  ## ----------------------------------------------------------------
  ## Plot legend
  ## ----------------------------------------------------------------
  if (!missing(legend)) {
    if (norm == "perc") {
      max.perc <- max(apply(counts, 2, function(x) { x / sum(x)}), na.rm = TRUE)
      # print(paste("max.perc", max.perc))
      legend.perc <- seq(0, max.perc, by = 0.1)
      # print(legend.perc)
      max.size <- sum(counts[, 1])
      # print("here")
      # print(legend.perc)
      legend.sizes <- 0.4 * ((max.size * 0.1 + (max.size * legend.perc)) ** 0.5 / (1.1 * max.size) ** 0.5)
      # print(legend.sizes)    
      num.counts <- ncol(counts)
      # palette(c("lightblue", "pink", "yellow3"))
      for (i in length(legend.sizes):1) {
        plot.node(x = max.x + 0.2,
                  y = max.y + legend.sizes[1] - sum(legend.sizes[1:i]) / 2,
                  size = matrix(rep(legend.sizes[i], num.counts), nrow = 1),
                  counts = rep(legend.perc[i], num.counts),
                  plot.internal.lines = F,
                  node.color = ifelse(use.colors.for.nodes, "grey", NA))
      }
      for (i in length(legend.sizes):1) {
        text(max.x + 0.2,
             max.y + legend.sizes[1] - sum(legend.sizes[1:i]) / 2 - legend.sizes[i], adj = c(0.5, -0.5),
             cex = legend.perc.cex,
             paste0(" ", legend.perc[i] * 100, "%"))
      }
      if (length(legend) > 1) {
        plot.node(x = max.x - 0.4, y = max.y,
                  size = matrix(seq(from = 0.2, to = 0.3, length.out = num.counts), nrow = 1),
                  counts = legend,
                  plot.counts = T,
                  plot.counts.colour = "black",
                  plot.counts.cex = legend.samples.cex,
                  plot.internal.lines = T,
                  node.color = ifelse(use.colors.for.nodes, "grey", NA))
      }

    }
  }
  ## ----------------------------------------------------------------

  par("mar" = prev.mar)
}

plot_cnv_boxplot <- function(this_cnv, segnorm_data, max = 100) {
  cnv_raw <- segnorm_data %>%
    filter(CHR == this_cnv$chr & START >= this_cnv$start & END <= this_cnv$end) %>%
    select(-c(CHR, START, END))
  par(mar = c(5,2,3,1))
  top <- ncol(cnv_raw)
  these_cells <- strsplit(as.character(this_cnv$cells), ",")[[1]]
  for (s in 1:ceiling(top / max)) {
    this_cnv_raw <- cnv_raw[, ((s - 1) * max + 1):(min((s) * max, top))]
    boxplot(this_cnv_raw * 2, ylim = c(0,5), las = 2,
            main = paste0("Distribution of CN estimates across bins\n",
                          this_cnv$chr, ":", this_cnv$start, "-", this_cnv$end,
                          " (", ifelse(this_cnv$n < 2, "loss", "gain"), ")"),
            cex.axis = 0.5)
    pos <- which(colnames(this_cnv_raw) %in% these_cells)
    if (length(pos) > 0) {
      if (this_cnv$n > 2) {
        rect(pos - 0.5, 0, pos + 0.5, 6, border = "red", col = rgb(1, 0, 0, 0.2))
      } else {
        rect(pos - 0.5, 0, pos + 0.5, 6, border = "blue", col = rgb(0, 0, 1, 0.2))
      }
    }
    abline(h = 0:5)
  }
}

# tree.arcs <- data.frame(
#   parent = c(1,2,3,3,3,6),
#   child =  c(2,3,4,5,6,7),
#   label =  c("CDKN2A(+/-)\nPAX5(+/-)",
#              "CDKN2A(-/-)\n",
#              "FGFR3(+/-)",
#              "\nARHGAP26(+/-)",
#              "chr10:46-49M-",
#              "CDKN1B+"))
# 
# plot.tree(tree.arcs,
#           root = 1,
#           direction = "down",
#           counts = data.frame(s1 = c(0,5,5,5,5,5,15),
#                               s3 = c(0,1,1,1,1,1,10)),
#           legend = c("Day 0", "Day 28"), legend.samples.cex = 0.8,
#           palette = c(rgb(0.2,0.7,0.4), rgb(0.5,1,0.5)),
#           plot.counts = T,
#           norm = "perc",
#           lab.cex = 0.8)
# plot.tree(tree.arcs, root = 1, direction = "down", counts = data.frame(s1 = c(0,5,5,5,5,5,15),
#                                                                       s3 = c(0,1,1,1,1,1,10)), norm = "perc", lab.cex = 0.8)
# plot.tree(tree.arcs, root = 1, direction = "down", counts = data.frame(s1 = c(0,5,5,5,5,5,15),
#                                                                       s3 = c(0,1,1,1,1,1,10)), norm = "total", lab.cex = 0.8)
# plot.tree(tree.arcs, root = 1, direction = "down", counts = data.frame(s1 = c(0,1,1,3,8,1,10),
#                                                                       s3 = c(0,2,4,6,8,10,20)), norm = "total", lab.cex = 0.8)
# plot.tree(tree.arcs, root = 1, direction = "down", counts = data.frame(s1 = c(0,1,1,10,0,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,0,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,0,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1)), norm = "sample")
# 
# plot.tree(tree.arcs, root = 1, direction = "down", counts = data.frame(s1 = c(0,11,11,14,10,11,11),
#                                                                       s3 = c(0,2,4,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,4,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,0,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1)), norm = "total")
# 
# plot.tree(tree.arcs, root = 1, direction = "left", counts = data.frame(s1 = c(0,1,1,10,0,1,1),
#                                                                       s3 = c(0,2,1,1,3,1,1)), norm = "total")
# plot.tree(tree.arcs, root = 1, direction = "right", counts = data.frame(s1 = c(0,1,1,1,0,1,1)), norm = "total")
# 


########################################################################################################################
# TEST plotGenome()
########################################################################################################################
# test_segcopy_data <- matrix(2, nrow = 400, ncol = 20)
# test_segcopy_data[sample.int(length(test_segcopy_data), 50)] <- 0
# test_segcopy_data[sample.int(length(test_segcopy_data), 100)] <- 1
# test_segcopy_data[sample.int(length(test_segcopy_data), 100)] <- 3
# test_segcopy_data[sample.int(length(test_segcopy_data), 50)] <- 4
# set.seed(100)
# c1 <- sample(20, 15)
# c2 <- setdiff(1:20, c1)
# c1.1 <- sample(c1, 10)
# c1.1.1 <- sample(c1.1, 5)
# test_segcopy_data[c(31:40), c1] <- 1
# test_segcopy_data[c(381:400), c2] <- 0
# test_segcopy_data[c(181:200), c1.1] <- 3
# test_segcopy_data[c(141:160), c1.1.1] <- 4
# colnames(test_segcopy_data) <- paste0("cell", 1:20)
# test_segcopy_data <- data.frame(
#   CHR = rep(paste0("chr", c(1:20)), each = 20),
#   START = rep(0:19, 20) * 1000000 + 1,
#   END = rep(1:20, 20) * 1000000,
#   test_segcopy_data
# )
# library(ape)
# test_tree <- as.phylo(hclust(dist(t(test_segcopy_data[, -c(1:3)])), method = "complete"))
# plotGenome(test_segcopy_data, tree = test_tree)

########################################################################################################################
# TEST plotNormGenome()
########################################################################################################################
# test_segnorm_data <- matrix(0.9 + range(rbeta(8000, 1, 1)) * 0.2, nrow = 400, ncol = 20)
# test_segnorm_data[sample.int(length(test_segnorm_data), 50)] <- range(rbeta(50, 1, 1)) * 0.2
# test_segnorm_data[sample.int(length(test_segnorm_data), 100)] <- 0.5 + range(rbeta(100, 1, 1)) * 0.2 - 0.1
# test_segnorm_data[sample.int(length(test_segnorm_data), 100)] <- 1.5 + range(rbeta(100, 1, 1)) * 0.2 - 0.1
# test_segnorm_data[sample.int(length(test_segnorm_data), 50)] <- 2 + range(rbeta(50, 1, 1)) * 0.2 - 0.1
# set.seed(100)
# c1 <- sample(20, 15)
# c2 <- setdiff(1:20, c1)
# c1.1 <- sample(c1, 10)
# c1.1.1 <- sample(c1.1, 5)
# test_segnorm_data[c(31:40), c1] <- 0.5 + rbeta(150, 2, 2) * 0.2 - 0.1
# test_segnorm_data[c(381:400), c2] <- 0 + rbeta(100, 2, 2) * 0.2
# test_segnorm_data[c(181:200), c1.1] <- 1.5 + rbeta(200, 2, 2) * 0.2 - 0.1
# test_segnorm_data[c(141:160), c1.1.1] <- 2 + rbeta(100, 2, 2) * 0.2 - 0.1
# colnames(test_segnorm_data) <- paste0("cell", 1:20)
# test_segnorm_data <- data.frame(
#   CHR = rep(paste0("chr", c(1:20)), each = 20),
#   START = rep(0:19, 20) * 1000000 + 1,
#   END = rep(1:20, 20) * 1000000,
#   test_segnorm_data
# )
# plotNormGenome(test_segnorm_data, tree = test_tree)
# # heatmap.2(as.matrix(t(test_segnorm_data[, -c(1:3)])), trace = 'none', Colv = NA, dendrogram = "row", labCol = NA)

########################################################################################################################
# TEST plotNormGenome()
########################################################################################################################
# test_segnorm_data <- matrix(0.9 + range(rbeta(8000, 1, 1)) * 0.2, nrow = 400, ncol = 20)
# test_segnorm_data[sample.int(length(test_segnorm_data), 50)] <- range(rbeta(50, 1, 1)) * 0.2
# test_segnorm_data[sample.int(length(test_segnorm_data), 100)] <- 0.5 + range(rbeta(100, 1, 1)) * 0.2 - 0.1
# test_segnorm_data[sample.int(length(test_segnorm_data), 100)] <- 1.5 + range(rbeta(100, 1, 1)) * 0.2 - 0.1
# test_segnorm_data[sample.int(length(test_segnorm_data), 50)] <- 2 + range(rbeta(50, 1, 1)) * 0.2 - 0.1
# set.seed(100)
# c1 <- sample(20, 15)
# c2 <- setdiff(1:20, c1)
# c1.1 <- sample(c1, 10)
# c1.1.1 <- sample(c1.1, 5)
# test_segnorm_data[c(31:40), c1] <- 0.5 + rbeta(150, 2, 2) * 0.2 - 0.1
# test_segnorm_data[c(381:400), c2] <- 0 + rbeta(100, 2, 2) * 0.2
# test_segnorm_data[c(181:200), c1.1] <- 1.5 + rbeta(200, 2, 2) * 0.2 - 0.1
# test_segnorm_data[c(141:160), c1.1.1] <- 2 + rbeta(100, 2, 2) * 0.2 - 0.1
# colnames(test_segnorm_data) <- paste0("cell", 1:20)
# test_segnorm_data <- data.frame(
#   CHR = rep(paste0("chr", c(1:20)), each = 20),
#   START = rep(0:19, 20) * 1000000 + 1,
#   END = rep(1:20, 20) * 1000000,
#   test_segnorm_data
# )
# library(ape)
# test_tree <- as.phylo(hclust(dist(t(test_segnorm_data[, -c(1:3)])), method = "complete"))
# plotNormGenome(test_segnorm_data, tree = test_tree)
# plotNormChr(test_segnorm_data, "chr2", tree = test_tree)
# # heatmap.2(as.matrix(t(test_segnorm_data[, -c(1:3)])), trace = 'none', Colv = NA, dendrogram = "row", labCol = NA)
