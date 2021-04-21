library(dplyr)
library(data.table)

get_annotation <- function(segcopy, meta_data, annotation) {
  if (!(annotation %in% colnames(meta_data))) {
    stop("Cannot find annotation ", filter.name, " in the meta_data")
  }
  annotation.values <- as.character(unique(meta_data[colnames(segcopy)[-c(1:3)], annotation]))
  return(annotation.values)  
}

get_table_of_cells_per_library <- function(segcopy, meta_data) {
  libraries <- get_annotation(segcopy = segcopy, meta_data = meta_data, annotation = "Library")
  cells <- sapply(libraries, grep, colnames(segcopy)[4:ncol(segcopy)], value = T)
  if (length(libraries) == 1) {
    cells.df <- data.frame("Library" = libraries,
                           "num.cells" = length(cells),
                           "Cells" = paste(cells, collapse = ", "))
  } else {
    cells.df <- data.frame("Library" = names(cells),
                           "num.cells" = sapply(cells, length),
                           "Cells" = sapply(cells, paste, collapse = ", "))
  }
  return(cells.df)
}


#' read_ginkgo_signal_files
#'
#' @param signal_files A list of file names (full paths) with the SegCopy or similar info
#'
#' @return Returns a data.frame with the data, one bin per row. The first 3 columns should be 'CHR', 
#'   'START', 'END' for the bin coordinates and then one column per cell with copy number estimates
#' @export
#'
#' @examples
#' segcopy_data <- read_ginkgo_signal_files(c("dir1/SegCopy", "dir2/SegCopy"))
#' segnorm_data <- read_ginkgo_signal_files(c("dir1/SegCopy", "dir2/SegCopy"))

read_ginkgo_signal_files <- function(signal_files) {
  data <- fread(signal_files[1], header = T, data.table = F)
  for (signal_file in signal_files[-1]) {
    df <- fread(signal_file, header = T, data.table = F)
    data <- cbind(data, df[, -c(1:3)])
  }
  return(data)
}

read_ginkgo_rawdata_files <- function(rawdata_files) {
  data <- fread(rawdata_files[1], header = T, data.table = F)
  for (rawdata_file in rawdata_files[-1]) {
    df <- fread(rawdata_file, header = T, data.table = F)
    data <- cbind(data, df)
  }
  return(data)
}

read_ginkgo_segstats_files <- function(segstats_files) {
  segstats_data <- read.table(segstats_files[1], header =T)
  for (segstats_file in segstats_files[-c(1)]) {
    segstats_data <- rbind(segstats_data, read.table(segstats_file, header =T))
  }
  return(segstats_data)
}

#' subset_signal_data
#'
#' @param 
#'
#' @return
#' @export
#'
#' @examples
subset_signal_data <- function(signal_data, meta_data, ...) {
  cells <- which_cells(signal_data, meta_data, ...)
  signal_data <- signal_data[, c("CHR", "START", "END", cells)]
  
  return(signal_data)  
}

#' subset_segstats_data
#'
#' @param 
#'
#' @return
#' @export
#'
#' @examples
subset_segstats_data <- function(segstats_data, meta_data, ...) {
  filters <- list(...)
  for (filter.name in names(filters)) {
    if (!(filter.name %in% colnames(meta_data))) {
      stop("Cannot subset signal_data on ", filter.name, " as this field does not exist in the meta_data")
    }
    segstats_data <- segstats_data[rownames(segstats_data) %in% as.character(meta_data[meta_data[, filter.name] == filters[[filter.name]], "Cell"]), ]
  }

  return(segstats_data)  
}

which_cells <- function(signal_data, meta_data, ...) {
  filters <- list(...)
  for (filter.name in names(filters)) {
    if (!(filter.name %in% colnames(meta_data))) {
      stop("Cannot subset signal_data on ", filter.name, " as this field does not exist in the meta_data")
    }
    signal_data <- signal_data[, colnames(signal_data) %in% as.character(meta_data[meta_data[, filter.name] == filters[[filter.name]], "Cell"])]
  }
  
  return(colnames(signal_data))
}

get_logR_from_norm_data <- function(norm.data) {
  logR.data <- norm.data
  logR.data[, -c(1:3)] <- log2(sweep(logR.data[, -c(1:3)],
                                     2,
                                     apply(logR.data[, -c(1:3)], 2, median),
                                     "/"))

  return(logR.data)
}

segcopy.to.cnv.df <- function(segcopy_data, genome, ignore_filters = FALSE) {
  cnv.df <- data.frame(matrix(NA, nrow = 0, ncol = 4), stringsAsFactors = FALSE);
  colnames(cnv.df) <- c("cell", "from", "to", "n")
  smallest.bin.size <- min(segcopy_data[, "END"] - segcopy_data[, "START"])
  add.cnv <- function(y, cell.name, copy) {
    r <- range(y);
    if (segcopy_data[r[1], "CHR"] != segcopy_data[r[2], "CHR"]) {
      last_bin <- max(which(segcopy_data[, "CHR"] == segcopy_data[r[1], "CHR"]))
      add.cnv(r[1]:last_bin, cell.name, copy)
      add.cnv((last_bin + 1):r[2], cell.name, copy)
      return()
    }
    if (ignore_filters) {
      cnv.df <<- rbind(cnv.df, c(cell.name, r[1], r[2], copy), stringsAsFactors = FALSE)
    } else {
      if (segcopy_data[r[1], "CHR"] == "chrY") {
        return()
      }
      if (segcopy_data[r[1], "CHR"] == segcopy_data[r[2], "CHR"]) {
        if ((r[2] - r[1] > 10) | (segcopy_data[r[2], "END"] - segcopy_data[r[1], "START"]) / length(y) < 2 * smallest.bin.size) {
          cnv.df <<- rbind(cnv.df, c(cell.name, r[1], r[2], copy), stringsAsFactors = FALSE)
        } else {
          print(paste0("Skipping CNV on ", cell.name, ", that spans centromere ", segcopy_data[r[1], "CHR"], ":",
                       segcopy_data[r[1], "START"], "-", segcopy_data[r[2], "END"]))
        }
      } else {
        print(paste0("Skipping CNV on ", cell.name, ", that spans chromosomes ends: ",
                     segcopy_data[r[1], "CHR"], " and ", segcopy_data[r[2], "CHR"]))
      }
    }
  }
  
  for (cell in 4:ncol(segcopy_data)) {
    d <- segcopy_data[, cell]
    cell.name <- colnames(segcopy_data)[cell]
    x <- which(d == 0)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 0)
    }
    x <- which(d == 1)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 1)
    }
    x <- which(d == 3)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 3)
    }
    x <- which(d == 4)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 4)
    }
    x <- which(d == 5)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 5)
    }
    x <- which(d >= 6)
    if (length(x) > 0) {
      tapply(x, c(0, cumsum(diff(x) != 1)), add.cnv, cell.name, 6)
    }
  }
  
  cnv.df <- data.frame(cell = as.factor(cnv.df[, 1]), from = as.integer(cnv.df[, 2]),
                       to = as.integer(cnv.df[, 3]), n = as.integer(cnv.df[, 4]))
  
}

## fix.segcopy is to create a new segcopy object based on the cnv.df (after removing offending CNVs)
fix.segcopy_data <- function(data, cnv.df_or_final_cnvs) {
  # This method can take eitehr cnv.df or final_cnvs:
  ## cnv.df structure: a DF with cell, from (#bin), to (#bin) and n (CN state)
  ##              cell  from    to n
  ## 1 VT15b_D701_D501 11074 11162 0
  ## 2 VT15b_D701_D501  1404  1432 1
  ## 3 VT15b_D701_D501  4956  4975 1
  ## 4 VT15b_D701_D501  5909  5913 1
  ## 5 VT15b_D701_D501  7951  7958 1
  ##
  ## final_cnvs structure: a DF with from (#bin), to (#bin), n (CN state), num.cells, cells, chr, start, end
  ##   from   to n num.cells                          cells  chr     start       end
  ##   5914 5921 3         2 VT15b_D701_D501,VT22_D712_D504 chr9  22231215  24275135
  ##   6269 6296 3         1                VT15b_D701_D501 chr9 133929996 141213431
  ##
  fixed.segcopy_data <- data
  fixed.segcopy_data[, -c(1:3)] <- 2
  if ("cell" %in% colnames(cnv.df_or_final_cnvs)) {
    for (i in 1:nrow(cnv.df_or_final_cnvs)) {
      fixed.segcopy_data[cnv.df_or_final_cnvs[i, "from", drop = T]:cnv.df_or_final_cnvs[i, "to", drop = T],
                         as.character(cnv.df_or_final_cnvs[i, "cell"])] <-
        cnv.df_or_final_cnvs[i, "n"]
    }
  } else if ("cells" %in% colnames(cnv.df_or_final_cnvs)) {
    for (i in 1:nrow(cnv.df_or_final_cnvs)) {
      fixed.segcopy_data[cnv.df_or_final_cnvs[i, "from", drop = T]:cnv.df_or_final_cnvs[i, "to", drop = T],
                         unique(strsplit(as.character(cnv.df_or_final_cnvs[i, "cells"]), ",")[[1]])] <-
        cnv.df_or_final_cnvs[i, "n"]
    }
  } else {
    stop("Can't find the right column. Check your call")
  }
  return(fixed.segcopy_data)
}

## Cluster the CNA events as edges may not be very reliable
cluster_cnvs <- function(cnv.df, segcopy_data, matches_threshold = 5) {
  ## cnv.df structure: a DF with cell, from (#bin), to (#bin) and n (CN state)
  ##              cell  from    to n
  ## 1 VT15b_D701_D501 11074 11162 0
  ## 2 VT15b_D701_D501  1404  1432 1
  ## 3 VT15b_D701_D501  4956  4975 1
  ## 4 VT15b_D701_D501  5909  5913 1
  ## 5 VT15b_D701_D501  7951  7958 1
  matches <- vector("list", nrow(cnv.df))
  for (this_cnv in 1:nrow(cnv.df)) {
    if (segcopy_data[cnv.df[this_cnv, "from"], "CHR"] != segcopy_data[cnv.df[this_cnv, "to"], "CHR"]) {
      stop("Wrong start and end for a CNV: it overlaps a chromosome boundary")
    }
    matches[[this_cnv]] <- which(abs(cnv.df[, "from"] - cnv.df[this_cnv, "from"]) <= matches_threshold &
                                   abs(cnv.df[, "to"] - cnv.df[this_cnv, "to"]) <= matches_threshold &
                                   segcopy_data[cnv.df[, "from"], "CHR"] == segcopy_data[cnv.df[this_cnv, "from"], "CHR"] &
                                   cnv.df[, "n"] == cnv.df[this_cnv, "n"])
  }

  final_cnvs = data.frame(matrix(NA, nrow = 0, ncol = 5))
  for (this_cnv in 1:nrow(cnv.df)) {
    ## candidates is the list of all unique matches for the matches of this CNA.
    candidates <- unique(unlist(matches[matches[[this_cnv]]]))
    ## If there is any left (would be removed if this CNA has been accounted for already)
    if (length(candidates) > 0) {
      candidates <- candidates[!is.na(candidates)]
    }
    while(length(candidates) > 0) {
      ## Find the mode "from" and "to" bins by looking at the max in the summary table
      t.from <- table(cnv.df[candidates, "from"])
      bin.from <- as.integer(names(which(t.from == max(t.from))[1]))
      t.to <- table(cnv.df[candidates, "to"])
      bin.to <- as.integer(names(which(t.to == max(t.to))[1]))
      
      # cnv.df[candidates, ]
      
      ## Only consider the CNA that still match the threshold rule w.r.t. the mode "from" and "to" bins
      final <- which(abs(cnv.df[, "from"] - bin.from) <= matches_threshold &
                       abs(cnv.df[, "to"] - bin.to) <= matches_threshold &
                       cnv.df[, "n"] == cnv.df[this_cnv, "n"])
      
      ## Alternatively, just take them all
      # final <- candidates
      
      ## Make sure we only include the ones that were candidate initially!
      final <- final[final %in% candidates]
      if (length(final) == 0) {
        final <- candidates[1]
        bin.from <- cnv.df[final, "from"]
        bin.to <- cnv.df[final, "to"]
      }

      ## Include the new CNA
      final_cnvs <- rbind(final_cnvs,
                          c(bin.from, bin.to, cnv.df[this_cnv, "n"],
                            length(unique(cnv.df[final, "cell"])),
                            paste(unique(cnv.df[final, "cell"]), collapse = ",")),
                          stringsAsFactors = FALSE)
      
      ## Remove the instances used in the new CNA from the list of matches
      for (c in candidates) {
        if (length(matches[[c]]) > 0) {
          matches[[c]] <- matches[[c]][!(matches[[c]] %in% final)]
        }
      }
      for (f in final) {
        matches[[f]] <- NA
      }
      ## Get a new list of candidates with the remaining matches (if any)
      candidates <- unique(unlist(matches[matches[[this_cnv]]]))
    }
  }
  
  colnames(final_cnvs) <- c("from", "to", "n", "num.cells", "cells")
  final_cnvs <- data.frame(from = as.integer(final_cnvs$from),
                           to = as.integer(final_cnvs$to),
                           n = as.integer(final_cnvs$n),
                           num.cells = as.integer(final_cnvs$num.cells),
                           cells = final_cnvs$cells)
  ## TODO: This assumes that the "from" and "to" bins are in the same chromosome. It so happens to be the
  ## case for this dataset.
  if(any(segcopy_data[final_cnvs$from, "CHR"] != segcopy_data[final_cnvs$to, "CHR"])) {
    stop(paste("Not all CNVs start and end in the same chromosome. Assumptions are being violated: ",
               final_cnvs[which(segcopy_data[final_cnvs$from, "CHR"] != segcopy_data[final_cnvs$to, "CHR"]), "from"],
               final_cnvs[which(segcopy_data[final_cnvs$from, "CHR"] != segcopy_data[final_cnvs$to, "CHR"]), "to"],
               final_cnvs[which(segcopy_data[final_cnvs$from, "CHR"] != segcopy_data[final_cnvs$to, "CHR"]), "cells"], collapse = "\n"))
  }
  final_cnvs$chr <- segcopy_data[final_cnvs$from, "CHR"]
  final_cnvs$start <- segcopy_data[final_cnvs$from, "START"]
  final_cnvs$end <- segcopy_data[final_cnvs$to, "END"]
  
  return(final_cnvs)
}

get_clones_from_cnvs <- function(final_cnvs, cells, min.num.cells = 5, plot = F) {
  main_cnvs <- subset(final_cnvs, num.cells >= min.num.cells)

  if (missing(cells)) {
    cells <- sort(unique(strsplit(paste(as.character(final_cnvs$cells), collapse = ","), split = ",")[[1]]))
  }

  ## ==========================================================================
  ## Consolidated CNVs
  ## ==========================================================================
  ## consolidated_cnvs is a DF with a CNV per col. Each row is a cell with 0
  ## and 1 depending on whether that CNV is found on that cell or not. Eg:
  ##        chr1:10M-20M  chr3:45M-50M  chr9:0M-23M
  ## cell1             0             1            1
  ## cell2             1             1            0
  ## cell3             1             0            1
  ## cell4             0             1            1
  ## ==========================================================================
  consolidated_cnvs <- data.frame(matrix(0, nrow = length(cells), ncol = nrow(main_cnvs)), row.names = cells)
  # Set the name of the columns to the location of the CNV
  colnames(consolidated_cnvs) <- paste0(main_cnvs$chr,
                                        ":",
                                        round(main_cnvs$start / 1000000),
                                        "-",
                                        round(main_cnvs$end / 1000000),
                                        "M",
                                        case_when(main_cnvs$n == 0 ~ "(-/-)",
                                                  main_cnvs$n == 1 ~ "(+/-)",
                                                  main_cnvs$n == 2 ~ "(+/+)",
                                                  main_cnvs$n == 3 ~ "(+/+/+)",
                                                  main_cnvs$n > 3 ~ "(amp)"))
  colnames(consolidated_cnvs) <- sub("WGD.+", "WGD", colnames(consolidated_cnvs))
  colnames(consolidated_cnvs) <- sub("HyDp.+", "HyDp", colnames(consolidated_cnvs))
  # Fill up the consolidated CNVs DF using the data in main_cnvs
  null <- apply(main_cnvs, 1, function(x) {
    these.cells <- strsplit(as.character(x[5]), split = ",")[[1]]
    this_cnv <- paste0(as.character(x[6]), ":", round(as.numeric(x[7]) / 1000000), "-", round(as.numeric(x[8]) / 1000000), "M", case_when(x[3] == 0 ~ "(-/-)",
                                                                                                                                          x[3] == 1 ~ "(+/-)",
                                                                                                                                          x[3] == 2 ~ "(+/+)",
                                                                                                                                          x[3] == 3 ~ "(+/+/+)",
                                                                                                                                          x[3] > 3 ~ "(amp)"))
    if (grepl("WGD", this_cnv)) {this_cnv <- "WGD"}
    if (grepl("HyDp", this_cnv)) {this_cnv <- "HyDp"}
    this_cnv <- gsub(pattern = " ", replacement = "", x = this_cnv)
    consolidated_cnvs[these.cells, this_cnv] <<- 1
  })

  # Plot the distribution of consolidated CN events
  if (plot) {
    cur.par <- par(no.readonly = TRUE)
    par("mar" = cur.par$mar + c(2, 4, 0, 0))
    barplot(sort(colSums(consolidated_cnvs)), horiz = T, las = 2, main = "Distribution of consolidated CN events")
    par(cur.par)
  }
  ## ==========================================================================
  
  # Transform <consolidated_cnvs> into a string of 0/1 for each CNA, a sort of barcode
  binary_consolidated_cnvs <- apply(consolidated_cnvs, 1, paste, collapse = "")
  
  # Focus on the clones that are represented by at least n cells.
  clones <- table(binary_consolidated_cnvs)[table(binary_consolidated_cnvs) >= min.num.cells]
  
  ## ==========================================================================
  ## Simplified Consolidated CNVs
  ## ==========================================================================
  ## The simplifed version contains only the CNVs that are considered in any of
  ## the previously defined clones
  ## ==========================================================================
  # Find the CN events that are not used in the definition of any of these clones
  unused_cnvs <- which(colSums(matrix(as.numeric(unlist(strsplit(names(clones), split = ""))), nrow = length(clones), byrow = T)) == 0)
  if (length(unused_cnvs) == ncol(consolidated_cnvs) - 1) {
    simplified_consolidated_cnvs <- data.frame(matrix(consolidated_cnvs[, -c(unused_cnvs)], nrow = 1))
  } else if (length(unused_cnvs) > 0) {
      simplified_consolidated_cnvs <- consolidated_cnvs[, -c(unused_cnvs)]
    } else {
    simplified_consolidated_cnvs <- consolidated_cnvs
  }
  
  # Transform <simplified_consolidated_cnvs> into a string of 0/1 for each CNA, a sort of barcode
  simplified_binary_consolidated_cnvs <- apply(simplified_consolidated_cnvs, 1, paste, collapse = "")

  # Plot the distribution of simplified CN events (as the ones found on common clones)
  if (plot) {
    cur.par <- par(no.readonly = TRUE)
    par("mar" = cur.par$mar + c(0, 4, 0, 0))
    barplot(sort(colSums(simplified_consolidated_cnvs)), horiz = T, las = 2, main = "Distribution of simplified CN events")
    par(cur.par)
  }
  ## ==========================================================================
  
  # Focus on the clones that are represented by at least n cells.
  simplified_clones <- table(simplified_binary_consolidated_cnvs)[table(simplified_binary_consolidated_cnvs) >= min.num.cells]
  
  if (length(simplified_clones) == 0) {
    simplified_clones <- data.frame(matrix(NA, nrow = 0, ncol = 2))
  } else {
    simplified_clones <- data.frame(simplified_clones)
  }
  colnames(simplified_clones) <- c("karyotype", "num.cells")
  simplified_clones$cells <- sapply(simplified_clones$karyotype, function(x) {
    paste(names(which(simplified_binary_consolidated_cnvs == x)), collapse = ",")
  })

  # Change the clone names (karyotype) to the combination of the locations of the CNA in <clones>.
  if (nrow(simplified_clones) > 0) {
    simplified_clones$karyotype <- factor(sapply(as.character(simplified_clones$karyotype), function(x) {
      if (x == "") {
        return("other")
      }
      elements <- strsplit(x, split = "")[[1]]
      name <- paste(colnames(simplified_consolidated_cnvs)[elements == "1"], collapse = " + ")
      # name <- gsub("(\\d+)(\\d)\\d{5}", "\\1.\\2M", name, perl = T)
      if (name == "") {
        name <- "normal"
      }
      return(name)
    }))
  }

  return(simplified_clones)
}

annotate_cells_by_clones <- function(meta_data, clones) {
  # Start building the data DF from the meta_data to have the Library, mouse, treatment...
  data <- meta_data
  data$Clone <- as.character("other")

  if (nrow(clones) > 0) {
    apply(clones, 1, function(x) {
      data[strsplit(x[3], ",")[[1]], "Clone"] <<- x[1]
    })
  }

  return(data)
}