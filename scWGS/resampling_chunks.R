## ---- initial_density_plots
for (x in 1:ncol(old_signal_matrix)) {
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## CN #", x, "\n\n"), quiet = T))
  }
  
  this_cnv <- cnv_table[x, ]
  title <- paste0(this_cnv$chr,
                  ":",
                  round(this_cnv$start/1000000), "M",
                  "-",
                  round(this_cnv$end/1000000), "M",
                  " [", this_cnv$n, "] - ",
                  this_cnv$n.bins,
                  " bins")
  # print(title)
  p <- old_signal_matrix %>%
    as_tibble() %>%
    ggplot() +
    geom_vline(xintercept = 0.5:4.5, col = "black", linetype = "dashed") +
    geom_density(aes(x = get(paste0("V", x)))) +
    geom_point(aes(x = get(paste0("V", x)),
                   y = rnorm(nrow(old_signal_matrix), 0, 0.4)^2)) +
    geom_rug(aes(x = get(paste0("V", x)))) +
    geom_vline(xintercept = starting_thresholds$mean[x], col = "pink") +
    stat_function(fun = dnorm, n = 201,
                  args = list(mean = starting_thresholds$mean[x],
                              sd = starting_thresholds$sd[x]),
                  col = "pink") +
    labs(title = title, x = "mean signal") +
    coord_cartesian(xlim = c(0.1, 4.9)) +
    theme_bw()
  print(p)
  
  plotNormGenome(segrenorm_data %>%
                   filter(CHR == this_cnv$chr & START >= this_cnv$start - 5000000 & END <= this_cnv$end + 5000000),
                 genes = list("from" = c(this_cnv$chr, this_cnv$start, "black"),
                              "to" = c(this_cnv$chr, this_cnv$end, "black")))
  
  plotGenome(single_seg_renorm[["10"]] %>%
               filter(CHR == this_cnv$chr & START >= this_cnv$start - 5000000 & END <= this_cnv$end + 5000000),
             genes = list("from" = c(this_cnv$chr, this_cnv$start, "black"),
                          "to" = c(this_cnv$chr, this_cnv$end, "black")))
}

## ---- pairwise_scatterplots ----
this_signal_matrix <- old_signal_matrix
for (y in 1:ncol(this_signal_matrix)) {
  # for (y in 1:1) {
  # x <- 5
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## CN #", y, " {.tabset}\n\n"), quiet = T))
  }
  title_y <- paste0(cnv_table[y, "chr"],
                    ":",
                    round(cnv_table[y, "start"]/1000000), "M",
                    "-",
                    round(cnv_table[y, "end"]/1000000), "M",
                    " [", cnv_table[y, "n"], "] - ",
                    cnv_table[y, "n.bins"],
                    " bins")
  
  this_cnv <- cnv_table[y, , drop = T]
  normal_cells <- apply(single_seg_renorm[["10"]] %>%
                          filter(CHR == this_cnv$chr & START >= this_cnv$start & END <= this_cnv$end) %>%
                          select(-c(1:3)), 2, function(m) {all(m == 2)})
  affected_cells <- apply(single_seg_renorm[["10"]] %>%
                            filter(CHR == this_cnv$chr & START >= this_cnv$start & END <= this_cnv$end) %>%
                            select(-c(1:3)), 2, function(m) {all(m == this_cnv$n)})
  class = ifelse(normal_cells, "normal", ifelse(affected_cells, "affected", "other"))
  
  for (x in 1:ncol(this_signal_matrix)) {
    if (!is.null(opts_knit$get("output.dir"))) {
      cat(knit_child(text = paste0("\n### CN #", x, "\n\n"), quiet = T))
    }
    title_x <- paste0(cnv_table[x, "chr"],
                      ":",
                      round(cnv_table[x, "start"]/1000000), "M",
                      "-",
                      round(cnv_table[x, "end"]/1000000), "M",
                      " [", cnv_table[x, "n"], "] - ",
                      cnv_table[x, "n.bins"],
                      " bins")
    p <- this_signal_matrix %>%
      as_tibble() %>%
      ggplot() +
      geom_density(aes(x = get(paste0("V", x)), fill = class, col = class), alpha = 0.2) +
      geom_point(aes(x = get(paste0("V", x)), y = get(paste0("V", y)), col = class)) +
      geom_rug(aes(x = get(paste0("V", x)), y = get(paste0("V", y)), col = class)) +
      scale_color_manual(limits = c("normal", "affected", "other"), values = c("black", "red", "blue")) +
      scale_fill_manual(limits = c("normal", "affected", "other"), values = c("black", "red", "blue")) +
      guides(alpha = FALSE) +
      labs(x = title_x, y = title_y) +
      geom_hline(yintercept = 0.5:4.5) +
      geom_vline(xintercept = 0.5:4.5) +
      coord_equal(xlim = c(0.1, 4.9), ylim = c(0.1, 4.9)) +
      theme_bw()
    print(p)
}
}


## ---- compare_original_segcopy_data_to_final_data_matrix ----
compare_original_segcopy_data_to_final_data_matrix <- function(original_segcopy_data, cnv_table, data_matrix, min_overlap = 30) {
  all_cells <- colnames(original_segcopy_data)[-c(1:3)]
  original_cnvs <- get_cnv_table(original_segcopy_data)
  overlapping_cnvs <- original_cnvs
  missing_cnvs <- original_cnvs
  rownames(data_matrix) <- all_cells
  for (i in 1:nrow(original_cnvs)) {
    this_original_cnv <- original_cnvs[i, ]
    this_chr <- this_original_cnv$chr
    this_start <- this_original_cnv$start
    this_end <- this_original_cnv$end
    this_n <- this_original_cnv$n
    final_cnv_events <- cnv_table %>% mutate(idx = row_number()) %>% filter(chr == this_chr) %>%
      filter(n == this_n) %>%
      filter(start < this_end & end > this_start) %>%
      mutate(overlap = (min(this_end, end) - max(this_start, start)) / (this_end - this_start) * 100) %>%
      filter(overlap > min_overlap)
    if (nrow(final_cnv_events) > 0) {
      original_cells <- strsplit(this_original_cnv$cells, split = ",")[[1]]
      final_cells <- unlist(lapply(final_cnv_events$idx, function(idx) {
        names(which(data_matrix[, idx] == this_n))}))
      missing_cells <- setdiff(original_cells, final_cells)
      overlapping_cells <- intersect(original_cells, final_cells)
      
      overlapping_cnvs[i, "cells"] <- paste(overlapping_cells, collapse = ",")
      overlapping_cnvs[i, "num.cells"] <- length(overlapping_cells)
      missing_cnvs[i, "cells"] <- paste(missing_cells, collapse = ",")
      missing_cnvs[i, "num.cells"] <- length(missing_cells)
    } else {
      overlapping_cnvs[i, "cells"] <- ""
      overlapping_cnvs[i, "num.cells"] <- 0
    }
  }
  attr(missing_cnvs, "overlapping_cnvs") <- overlapping_cnvs
  return(missing_cnvs)
}

original_segcopy_data <- single_seg_renorm[["10"]]
missing_cnvs <- compare_original_segcopy_data_to_final_data_matrix(
  original_segcopy_data, cnv_table, resampled_segcopy_data)
overlapping_cnvs <- attr(missing_cnvs, "overlapping_cnvs")
ggplot() +
  geom_density(data = overlapping_cnvs,
               aes(x = (end - start) / 1000000, weight = num.cells/sum(num.cells)),
               col = "blue", bw = 0.05) +
  geom_density(data = missing_cnvs,
               aes(x = (end - start) / 1000000, weight = num.cells / sum(num.cells)),
               col = "red", bw = 0.05) +
  coord_trans(x = "log10") +
  labs(x = "length (Mbs)")
plotFrequency(fix.segcopy_data(original_segcopy_data, overlapping_cnvs))
title("Overlappping")

plotFrequency(original_segcopy_data)
title("Original (single segmentation gamma = 10)")
plotFrequency(multi_seg_renorm[["05"]])
title("Original (multi segmentation gamma = 05)")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], overlapping_cnvs))
title("All overlapping")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], missing_cnvs))
title("All missing")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], missing_cnvs %>% filter((end - start) > 1000000)))
title("1Mb+ missing")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], missing_cnvs %>% filter((end - start) > 3000000)))
title("3Mb+ missing")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], missing_cnvs %>% filter((end - start) > 5000000)))
title("5Mb+ missing")
plotFrequency(fix.segcopy_data(single_seg_renorm[["10"]], missing_cnvs %>% filter((end - start) > 10000000)))
title("10Mb+ missing")


## ---- genome_plots_with_clones
cosmic_list <- function() {
  cosmic_list <- list()
  cosmic <- read.csv("../data/cosmic/grch37/cancer_gene_census.grch37.v90.csv", as.is = T)
  cosmic.loc <- strsplit(cosmic$Genome.Location, ":")
  cosmic$chr <- paste0("chr", sapply(cosmic.loc, function(x) {x[[1]]}))
  cosmic$start <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][1]}))
  cosmic$end <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][2]}))
  for (i in 1:nrow(cosmic)) {
    cosmic_list[[cosmic$Gene.Symbol[i]]] <- c(cosmic$chr[i],
                                              (cosmic$start[i] + cosmic$end[i]) / 2,
                                              scales::muted("red"))
  }
  return(cosmic_list)
}
cosmic_list <- cosmic_list()

for (x in 1:ncol(old_signal_matrix)) {
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## CN #", x, "\n\n"), quiet = T))
  }
  
  this_cnv <- cnv_table[x, ]
  title <- paste0(this_cnv$chr,
                  ":",
                  round(this_cnv$start/1000000), "M",
                  "-",
                  round(this_cnv$end/1000000), "M",
                  " [", this_cnv$n, "] - ",
                  this_cnv$n.bins,
                  " bins")
  cat(title)

  plotNormGenome(segrenorm_data %>%
                   filter(CHR == this_cnv$chr & START >= this_cnv$start - 5000000 & END <= this_cnv$end + 5000000),
                 clones = as.factor(apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})),
                 genes = c(list("from" = c(this_cnv$chr, this_cnv$start, "black"),
                                "to" = c(this_cnv$chr, this_cnv$end, "black")),
                           cosmic_list))
  title("Raw data")
  
  cat("\n<p>Single-sample segmentation (Gamma = 10)</p>\n")
  plotGenome(single_seg_renorm[["10"]] %>%
               filter(CHR == this_cnv$chr & START >= this_cnv$start - 5000000 & END <= this_cnv$end + 5000000),
             clones = as.factor(apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})),
             genes = c(list("from" = c(this_cnv$chr, this_cnv$start, "black"),
                            "to" = c(this_cnv$chr, this_cnv$end, "black")),
                       cosmic_list))
  title("Single-sample segmentation (Gamma = 10)")

  plotGenome(multi_seg_renorm[["05"]] %>%
               filter(CHR == this_cnv$chr & START >= this_cnv$start - 5000000 & END <= this_cnv$end + 5000000),
             clones = as.factor(apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})),
             genes = c(list("from" = c(this_cnv$chr, this_cnv$start, "black"),
                            "to" = c(this_cnv$chr, this_cnv$end, "black")),
                       cosmic_list))
  title("Multi-sample segmentation (Gamma = 5)")
}


## ---- density_plots_post_tree ----
for (x in 1:ncol(old_signal_matrix)) {
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## CN #", x, "\n\n"), quiet = T))
  }
  
  title <- paste0(cnv_table[x, "chr"],
                  ":",
                  round(cnv_table[x, "start"]/1000000), "M",
                  "-",
                  round(cnv_table[x, "end"]/1000000), "M",
                  " [", cnv_table[x, "n"], "] - ",
                  cnv_table[x, "n.bins"],
                  " bins")
  print(title)
  
  p_top <- old_signal_matrix %>%
    as_tibble() %>%
    ggplot() +
    geom_vline(xintercept = 0.5:4.5, col = "black", linetype = "dashed") +
    geom_point(aes(x = get(paste0("V", x)),
                   y = apply(best_sample_data, 1:2, mean)[, x],
                   col = apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))}),
                   pch = as.factor(resampled_segcopy_data[, x])
    )) +
    scale_color_discrete(name = "CN clone") +
    scale_shape_discrete(name = "Final CN") +
    labs(x = x,
         y = "mean CN",
         title = title) +
    geom_vline(xintercept = mean(best_sample_thresholds[x, ]), col = "red") +
    coord_cartesian(ylim = c(min(2, cnv_table[x, "n", drop = T]),
                             max(2, cnv_table[x, "n", drop = T])),
                    xlim = c(0.1, 4.9)) +
    theme_bw()
  # print(p_top)
  
  max_density = max(density(old_signal_matrix[, paste0("V", x)])$y) * 0.9
  
  p <- old_signal_matrix %>%
    as_tibble() %>%
    ggplot() +
    geom_density(data = tibble(t = best_sample_thresholds[x, ]),
                 aes(x = t), col = "pink", fill = "pink") +
    geom_vline(xintercept = mean(best_sample_thresholds[x, ]), col = "red") +
    geom_vline(xintercept = 0.5:4.5, col = "black", linetype = "dashed") +
    geom_density(aes(x = get(paste0("V", x)))) +
    geom_point(aes(x = get(paste0("V", x)),
                   y = abs(vipor::offsetX(get(paste0("V", x)),
                                          width = max_density)),
                   col = apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))}),
                   pch = as.factor(resampled_segcopy_data[, x])
    )) +
    geom_rug(aes(x = get(paste0("V", x)),
                 col = apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})
    )) +
    scale_color_discrete(name = "CN clone") +
    scale_shape_discrete(name = "Final CN") +
    labs(x = "mean signal") +
    coord_cartesian(xlim = c(0.1, 4.9)) +
    theme_bw()
  pg <- plot_grid(
    plot_grid(p_top + theme(legend.position = "none", axis.title.x = element_blank()),
              p + theme(legend.position = "none"),
              nrow = 2, align = "v", rel_heights = c(1, 2)),
    get_legend(p),
    rel_widths = c(4, 1)
  )
  print(pg)
}


## ---- density_plots_per_clone ----
# saveRDS(resampled_segcopy_data, "resampled_segcopy_data.rds")
# resampled_segcopy_data <- readRDS("resampled_segcopy_data.rds")
clones <- apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})
clones_list <- sort(unique(clones))

for (this_clone.id in 1:length(clones_list)) {
  this_clone <- clones_list[this_clone.id]
  this_clone.cns <- strsplit(this_clone, "")[[1]]
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## Clone ", this_clone.id, " {.tabset}\n\n"), quiet = T))
  }
  cat(paste0("Clone ", paste(which(strsplit(this_clone, "")[[1]] == "|"), collapse = "+"), ": ", sum(clones == this_clone), " cells.\n\n"))
  for (x in 1:ncol(old_signal_matrix)) {
    if (!is.null(opts_knit$get("output.dir"))) {
      cat(knit_child(text = paste0("\n### CN #", x, " {.tabset}\n\n"), quiet = T))
    }
    title <- paste0("Clone ", paste(which(strsplit(this_clone, "")[[1]] == "|"), collapse = "+"),
                    " -- ",
                    cnv_table[x, "chr"],
                    ":",
                    round(cnv_table[x, "start"]/1000000), "M",
                    "-",
                    round(cnv_table[x, "end"]/1000000), "M",
                    " [", cnv_table[x, "n"], "] - ",
                    cnv_table[x, "n.bins"],
                    " bins")
    p_top <- old_signal_matrix %>%
      as_tibble() %>%
      ggplot() +
      geom_vline(xintercept = 0.5:4.5, col = "black", linetype = "dashed") +
      geom_point(aes(x = get(paste0("V", x)),
                     y = apply(best_sample_data, 1:2, mean)[, x],
                     col = clones,
                     alpha = ifelse(clones == this_clone, 1, 0),
                     size = ifelse(clones == this_clone, 1, 0),
                     pch = as.factor(resampled_segcopy_data[, x])
      )) +
      scale_color_discrete(name = "CN clone") +
      scale_shape_discrete(name = "Final CN") +
      labs(x = x,
           y = "mean CN",
           title = title) +
      geom_vline(xintercept = mean(best_sample_thresholds[x, ])) +
      scale_alpha_continuous(range = c(0.3, 1)) +
      scale_size_continuous(range = c(1, 3)) +
      coord_cartesian(ylim = c(min(2, cnv_table[x, "n", drop = T]),
                               max(2, cnv_table[x, "n", drop = T])),
                      xlim = c(0.1, 4.9)) +
      theme_bw()
    # print(p_top)
    
    max_density = max(density(old_signal_matrix[, paste0("V", x)])$y) * 0.9
    
    p <- old_signal_matrix %>%
      as_tibble() %>%
      ggplot() +
      geom_vline(xintercept = 0.5:4.5, col = "black", linetype = "dashed") +
      # geom_density(data = as.data.frame(old_signal_matrix)[clones == this_clone, , drop = F],
      #              aes(x = get(paste0("V", x)), col = this_clone),
      #              show.legend = F) +
      geom_density(aes(x = get(paste0("V", x)))) +
      # geom_density(aes(x = get(paste0("V", x)),
      #                 fill = as.factor(resampled_segcopy_data[, x])
      #                 ), alpha = 0.3) +
      geom_point(aes(x = get(paste0("V", x)),
                     y = abs(vipor::offsetX(get(paste0("V", x)),
                                            width = max_density)),
                     col = clones,
                     alpha = ifelse(clones == this_clone, 1, 0),
                     size = ifelse(clones == this_clone, 1, 0),
                     pch = as.factor(resampled_segcopy_data[, x])
      )) +
      geom_rug(aes(x = get(paste0("V", x)),
                   col = apply(resampled_segcopy_data, 1, function(x) {gsub("2", ".", paste(ifelse(x==2, ".", "|"), collapse = ""))})
      )) +
      scale_alpha_continuous(range = c(0.3, 1)) +
      scale_size_continuous(range = c(1, 3)) +
      scale_color_discrete(name = "CN clone") +
      scale_shape_discrete(name = "Final CN") +
      guides(size = F, alpha = F) +
      labs(x = "mean signal") +
      geom_vline(xintercept = mean(best_sample_thresholds[x, ])) +
      coord_cartesian(xlim = c(0.1, 4.9)) +
      theme_bw()
    pg <- plot_grid(
      plot_grid(p_top + theme(legend.position = "none", axis.title.x = element_blank()),
                p + theme(legend.position = "none"),
                nrow = 2, align = "v", rel_heights = c(1, 2)),
      get_legend(p),
      rel_widths = c(4, 2)
    )
    print(pg)
    # break
  }
  # break
}


## ---- cosmic_annotations ----
datatable(NULL)
cosmic_annotations <- tibble()
cosmic <- read.csv("../data/cosmic/grch37/cancer_gene_census.csv", as.is = T)
cosmic.loc <- strsplit(cosmic$Genome.Location, ":")
cosmic$chr <- paste0("chr", sapply(cosmic.loc, function(x) {x[[1]]}))
cosmic$start <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][1]}))
cosmic$end <- as.integer(sapply(cosmic.loc, function(x) {strsplit(x[2], "-")[[1]][2]}))
for (i in 1:nrow(cnv_table)) {
  this_cnv <- cnv_table[i, , drop = FALSE]
  tag <- paste0("cosmic_", this_cnv$chr, "_", this_cnv$start, "_", this_cnv$end)
  location <- paste0(sub("chr", "", this_cnv$chr), ":", this_cnv$start, "-", this_cnv$end)
  label <- paste0(this_cnv$chr, ":", round(this_cnv$start / 1000000), "-", round(this_cnv$end / 1000000), "M")
  if (!is.null(opts_knit$get("output.dir"))) {
    cat(knit_child(text = paste0("\n## ", label, "\n\n"), quiet = T))
  }
  # cat(paste0('<div id="cosmic-', tag, '" class="section level2"><h2>', label, '</h2>'))
  cosmic.strict_matches <- cosmic %>%
    filter(chr == this_cnv$chr) %>%
    filter(start < this_cnv$end) %>%
    filter(end > this_cnv$start)
  cosmic.wider_matches <- cosmic %>%
    filter(chr == this_cnv$chr) %>%
    filter(start < this_cnv$end + 100000) %>%
    filter(end > this_cnv$start - 100000)
  print(htmltools::tagList(
    htmltools::h4("Cosmic"),
    datatable(elementId = label,
              cosmic.wider_matches,
              class = 'display',
              height = "auto",
              width = 800,
              options = list(columnDefs = list(
                list(width = '80px', targets = list(0,1,3,5,6,7,19,20,21)),
                list(width = '150px', targets = list(2,4,8,9,10,11,12,13,14,15,16,17,18))
              ),
              scrollX = TRUE),
              escape = FALSE),
    htmltools::h4("Ensembl"),
    htmltools::a(href = paste0("http://grch37.ensembl.org/Homo_sapiens/Location/Overview?r=", location),
                 target = "_blank",
                 paste("Link to Ensembl GRCh37 for", location))
  ))
  # cat('</div>')
  cosmic_annotations <- rbind(
    cosmic_annotations,
    tibble(lesion = label, genes = paste(cosmic.wider_matches %>% pull(Gene.Symbol), collapse = "/"))
  )
}

## ---- annotated_tree ----
event_names <- paste0(cnv_table$chr,
                      ":",
                      round(cnv_table$start / 1000000),
                      "-",
                      round(cnv_table$end / 1000000),
                      "M",
                      case_when(cnv_table$n == 0 ~ "(-/-)",
                                cnv_table$n == 1 ~ "(+/-)",
                                cnv_table$n == 2 ~ "(+/+)",
                                cnv_table$n == 3 ~ "(+/+/+)",
                                cnv_table$n == 4 ~ "(4)",
                                cnv_table$n == 5 ~ "(5)",
                                cnv_table$n > 5 ~ "(6)"))
event_names <- sub("WGD.+", "WGD", event_names)
event_names <- sub("HyDp.+", "HyDp", event_names)

tree.arcs$events <- sapply(tree.arcs$label, function(l) {
  paste(event_names[as.integer(strsplit(l, ",")[[1]])], collapse = " + ")})

tree.nodes <- tibble(node_id = c(1), depth = 0, x0 = 0)[0, ]
x0 <- 0

tree.labels <- tibble()[0, ]

update_tree.nodes <- function(this_node = 1, depth = 0) {
  children <- tree.arcs %>% filter(parent == this_node) %>% pull(child)
  for (this_child in children) {
    num.events <- tree.arcs %>% filter(child == this_child) %>% pull(label) %>%
      gregexpr(",", .) %>% .[[1]] %>% length()
    update_tree.nodes(this_child, depth + num.events)
    parent_events <- tree.arcs %>% filter(child == this_node) %>%
      pull(events) %>% strsplit(" + ", fixed = T)
    if (length(parent_events) > 0) {
      parent_events <- parent_events[[1]]
    }
    this_child_events <- tree.arcs %>% filter(child == this_child) %>%
      pull(events) %>% strsplit(" + ", fixed = T) %>% .[[1]]
    these_new_events <- setdiff(this_child_events, parent_events)
    for (i in 1:length(these_new_events)) {
      event <- these_new_events[i]
      if (grepl("(+/-)", event, fixed = T) | grepl("(-/-)", event, fixed = T)) {
        tree.labels <<- rbind(tree.labels,
                              tibble(node_id = this_child,
                                     label = event,
                                     num = i,
                                     tot = length(these_new_events),
                                     color = "darkblue")
        )
      } else {
        tree.labels <<- rbind(tree.labels,
                              tibble(node_id = this_child,
                                     label = event,
                                     num = i,
                                     tot = length(these_new_events),
                                     color = "red")
        )
      }
    }
  }
  if (length(children) == 0) {
    tree.nodes <<- tree.nodes %>%
      rbind(tibble(node_id = this_node,
                   depth = depth,
                   x0 = x0,
                   clone = tree.arcs %>% filter(child == this_node) %>%
                     pull("label")))
    x0 <<- x0 + 2
  } else {
    r <- tree.nodes %>% filter(node_id %in% children) %>% pull(x0) %>% range
    this_x0 <- (r[1] + r[2]) / 2
    this_clone <- tree.arcs %>% filter(child == this_node) %>% pull("label")
    if (length(this_clone) == 0) {
      this_clone <- "normal"
    }
    tree.nodes <<- tree.nodes %>%
      rbind(tibble(node_id = this_node,
                   depth = depth,
                   x0 = this_x0,
                   clone = this_clone))
  }
}

update_tree.nodes()

# Include the percentage of cells within the nodes

tree.nodes <- tree.nodes %>% mutate(y0 = max(depth) - depth)
tree.nodes <- tree.nodes %>%
  mutate(counts = attr(tree.arcs, "counts")[tree.nodes$node_id, 1, drop = T]) %>%
  mutate(perc = paste0(round(100 * counts / sum(counts)), "%")) %>%
  mutate(perc = ifelse(perc == "0%", "", perc))

arcs <- tree.arcs %>%
  left_join(tree.nodes, by = c("parent" = "node_id")) %>%
  select(parent, child, x = x0, y = y0, label) %>%
  left_join(tree.nodes, by = c("child" = "node_id")) %>%
  select(parent, child, x, y, xend = x0, yend = y0, label) %>%
  mutate(nudge_x = ifelse(x > xend, -0.2, ifelse(x < xend, 0.2, 0))) %>%
  mutate(nudge_y = ifelse(x > xend, -0.1, ifelse(x < xend, 0.1, 0)))

tree.labels <- arcs %>%
  left_join(tree.labels, by = c("child" = "node_id")) %>%
  mutate(x = x + (xend - x) * num / (tot + 1),
         y = y + (yend - y) * num / (tot + 1)) %>%
  select(x, y, label = label.y, color) %>%
  mutate(fill = "white") %>%
  mutate(label = sub("M.+", "M", label))

library(ggrepel)
ggplot() +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 5, col = "darkgrey") +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 4, col = "lightblue") +
  geom_label_repel(data = tree.labels,
                   aes(label = label, x = x, y = y),
                   hjust = 0.5,
                   size = 2,
                   color = tree.labels$color,
                   fill = tree.labels$fill,
                   min.segment.length = 0,
                   point.padding = NA, seed = 1) +
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 10, col = "black",
             fill = ifelse(tree.nodes$counts > 0, "darkgrey", "white"),
             shape = 21) +
  geom_text(data = tree.nodes,
            aes(x = x0, y = y0, label = perc),
            size = 3) +
  scale_x_continuous(expand = c(.1, .1)) +
  theme_void()

tree.labels <- tree.labels %>%
  left_join(cosmic_annotations, by = c("label" = "lesion")) %>%
  mutate(label = ifelse(genes != "",
                        gsub("([^/]+/[^/]+/[^/]+/[^/]+/)", "\\1\n", genes),
                        label))

ggplot() +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 5, col = "darkgrey") +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = 4, col = "lightblue") +
  geom_label_repel(data = tree.labels,
                   aes(label = label, x = x, y = y),
                   hjust = 0.5,
                   size = 2,
                   color = tree.labels$color,
                   fill = tree.labels$fill,
                   min.segment.length = 0,
                   point.padding = NA, seed = 1) +
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 10, col = "black",
             fill = ifelse(tree.nodes$counts > 0, "darkgrey", "white"),
             shape = 21) +
  geom_text(data = tree.nodes,
            aes(x = x0, y = y0, label = perc),
            size = 3) +
  scale_x_continuous(expand = c(.1, .1)) +
  theme_void()
