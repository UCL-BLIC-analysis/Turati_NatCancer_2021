## Only used for testing. The code used in the analysis is in the chunks below.
rm(list = ls())
load(Sys.glob("resampling_Pt11_100k_cache/html/define_cnvs_*.RData"))
load(Sys.glob("resampling_Pt11_100k_cache/html/get_signal_matrix_*.RData"))
load(Sys.glob("resampling_Pt11_100k_cache/html/resampling_tree_*.RData"))
load(Sys.glob("resampling_Pt11_100k_cache/html/sort_overlapping_events_*.RData"))
tree_title <- "Pt11 diagnosis"
# load(Sys.glob("resampling_Pt4_dia_100k_cache/html/sort_hyperdiploidy_*.RData"))
source("resampling_functions.R")
source("../R/plots.R")



cosmic_annotations <- tibble()
cosmic <- read.csv("../data/cosmic/grch37/cancer_gene_census.grch37.v90.csv", as.is = T)
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
  # cat('</div>')
  cosmic_annotations <- rbind(
    cosmic_annotations,
    tibble(lesion = label, genes = paste(cosmic.wider_matches %>% pull(Gene.Symbol), collapse = "/"))
  )
}


## ---- plot_final_trees ----
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
event_names <- sub(":0-0M", "", event_names, fixed = T)
event_names <- sub("WGD.+", "WGD", event_names)
event_names <- sub("HyDp.+", "HyDp", event_names)

tree.arcs$events <- sapply(tree.arcs$label, function(l) {
  paste(event_names[as.integer(strsplit(l, ",")[[1]])], collapse = " + ")})

update_tree.nodes <- function(this_node = 1, depth = 0, lineage = "") {
  children <- tree.arcs %>% filter(parent == this_node) %>% pull(child)
  this_clone <- tree.arcs %>% filter(child == this_node) %>% pull("label")
  this_lineage = paste(this_clone, lineage, sep = ",") %>%
    strsplit(",") %>% unlist() %>% as.numeric() %>% sort() %>%
    paste(collapse = ",")
  for (this_child in children) {
    num.events <- tree.arcs %>% filter(child == this_child) %>% pull(label) %>%
      strsplit(",") %>% .[[1]] %>% length()
    update_tree.nodes(this_child, depth + num.events, this_lineage)
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
                   x0_min = x0,
                   x0_max = x0,
                   clone = tree.arcs %>% filter(child == this_node) %>%
                     pull("label"),
                   lineage = this_lineage,
                   leaf = T))
    x0 <<- x0 + 1.5
  } else {
    r <- tree.nodes %>% filter(node_id %in% children) %>% pull(x0) %>% range
    this_x0 <- (r[1] + r[2]) / 2
    r <- tree.nodes %>% filter(node_id %in% children) %>% select(x0_min, x0_max) %>% range
    this_x0_min <- r[1]
    this_x0_max <- r[2]
    if (length(this_clone) == 0) {
      this_clone <- "normal"
    }
    tree.nodes <<- tree.nodes %>%
      rbind(tibble(node_id = this_node,
                   depth = depth,
                   x0 = this_x0,
                   x0_min = this_x0_min,
                   x0_max = this_x0_max,
                   clone = this_clone,
                   lineage = this_lineage,
                   leaf = F))
  }
}

tree.nodes <- tibble(node_id = c(1), depth = 0, x0 = 0)[0, ]
x0 <- 0

tree.labels <- tibble()[0, ]

update_tree.nodes()

# Include the percentage of cells within the nodes

tree.nodes <- tree.nodes %>% mutate(y0 = max(depth) - depth)
tree.nodes <- tree.nodes %>%
  mutate(counts = attr(tree.arcs, "counts")[tree.nodes$node_id, 1, drop = T]) %>%
  mutate(perc = paste0(round(100 * counts / sum(counts)), "%")) %>%
  mutate(perc.num = round(100 * counts / sum(counts))) %>%
  mutate(perc = ifelse(perc == "0%", "", perc))

arcs <- tree.arcs %>%
  left_join(tree.nodes, by = c("parent" = "node_id")) %>%
  select(parent, child, x = x0, y = y0, label) %>%
  left_join(tree.nodes, by = c("child" = "node_id")) %>%
  select(parent, child, x, y, xend = x0, yend = y0, label) %>%
  mutate(nudge_x = ifelse(x > xend, -0.2, ifelse(x < xend, 0.2, 0))) %>%
  mutate(nudge_y = ifelse(x > xend, -0.1, ifelse(x < xend, 0.1, 0)))

tree.labels.Orig <- tree.labels
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
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 10, col = "black",
             fill = ifelse(tree.nodes$counts > 0, "darkgrey", "white"),
             shape = 21) +
  geom_text(data = tree.nodes,
            aes(x = x0, y = y0, label = perc),
            size = 3) +
  geom_label_repel(data = tree.labels,
                   aes(label = label, x = x, y = y),
                   hjust = 0.5,
                   size = 2,
                   color = tree.labels$color,
                   fill = tree.labels$fill,
                   min.segment.length = 0,
                   point.padding = NA, seed = 1) +
  scale_x_continuous(expand = c(.1, .1)) +
  theme_void()

tree.labels <- tree.labels.Orig
tree.labels <- arcs %>%
  left_join(tree.labels, by = c("child" = "node_id")) %>%
  mutate(x = xend,
         y = y + (yend - y) * num / (tot + 1)) %>%
  select(x, y, label = label.y, color) %>%
  mutate(fill = "white") %>%
  mutate(label = sub("M.+", "M", label))

tree.labels <- tree.labels %>%
  left_join(cosmic_annotations, by = c("label" = "lesion")) %>%
  mutate(label = ifelse(genes != "",
                        gsub("([^/]+/[^/]+/[^/]+/[^/]+/)", "\\1\n", genes),
                        label))


ggplot() +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = y - 0.3),
               size = 5, col = "darkgrey", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = xend, y = y - 0.3, xend = xend, yend = yend),
               size = 5, col = "darkgrey", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = y - 0.3),
               size = 4, col = "lightblue", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = xend, y = y - 0.3, xend = xend, yend = yend),
               size = 4, col = "lightblue", lineend = "round") +
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

ggplot() +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = y - 0.3),
               size = 5, col = "darkgrey", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = xend, y = y - 0.3, xend = xend, yend = yend),
               size = 5, col = "darkgrey", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = x, y = y, xend = xend, yend = y - 0.3),
               size = 4, col = "lightblue", lineend = "round") +
  geom_segment(data = arcs,
               aes(x = xend, y = y - 0.3, xend = xend, yend = yend),
               size = 4, col = "lightblue", lineend = "round") +
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 5 + tree.nodes$perc.num / max(tree.nodes$perc.num) * 10,
             col = "black",
             fill = ifelse(tree.nodes$counts > 0, "darkgrey", "white"),
             shape = 21) +
  scale_x_continuous(expand = c(.1, .1)) +
  theme_void()


## ---- plot_final_trees_with_shared_events ----

y.start_of_tree <- 13
tree.nodes <- tree.nodes %>% mutate(y0 = 19 - depth * 1.5)

arcs <- tree.arcs %>%
  left_join(tree.nodes, by = c("parent" = "node_id")) %>%
  select(parent, child, x = x0, y = y0, label) %>%
  left_join(tree.nodes, by = c("child" = "node_id")) %>%
  select(parent, child, x, y, xend = x0, yend = y0, label) %>%
  mutate(nudge_x = ifelse(x > xend, -0.2, ifelse(x < xend, 0.2, 0))) %>%
  mutate(nudge_y = ifelse(x > xend, -0.1, ifelse(x < xend, 0.1, 0)))

tree.labels <- tree.labels.Orig

colored.bent_arcs <- arcs %>%
  left_join(tree.labels %>% filter(num == 1), by = c("child" = "node_id")) %>%
  expand_grid(segment = 1:2) %>%
  arrange(y, desc(abs(x - xend)), segment) %>%
  mutate(arc.x = ifelse(segment == 1, x, xend),
         arc.xend = xend,
         arc.y = ifelse(segment == 1, y, y - 0.3),
         arc.yend = ifelse(segment == 1, y - 0.3, y - 1.5)) %>%
  select(arc.x, arc.xend, arc.y, arc.yend, label = label.y, color)

colored.straight_arcs <- arcs %>%
  inner_join(tree.labels %>% filter(num > 1), by = c("child" = "node_id")) %>%
  mutate(arc.x = xend,
         arc.xend = xend,
         arc.y = y + (yend - y) * (num - 1) / tot,
         arc.yend = y + (yend - y) * num / tot) %>%
  select(arc.x, arc.xend, arc.y, arc.yend, label = label.y, color)

tree.labels <- tree.labels %>% group_by(label) %>% mutate(rep = n() > 1) %>% ungroup()

repeated_labels <- arcs %>%
  inner_join(tree.labels %>% filter(rep), by = c("child" = "node_id")) %>%
  mutate(lab.x = xend + 0,
         lab.y = y + (yend - y) * (num - 0.5) / tot) %>%
  select(lab.x, lab.y, label = label.y) %>%
  mutate(mark = letters[as.numeric(factor(label))])

ids <- expand_grid(repeated_labels %>% mutate(id = row_number()), tree.nodes) %>%
  filter(abs(lab.x - x0) < 0.5 & lab.y > y0 & lab.y - y0 < 1) %>% pull(id)
repeated_labels[ids, "lab.y"] <- repeated_labels[ids, "lab.y"] + 0.1
ids <- expand_grid(repeated_labels %>% mutate(id = row_number()), tree.nodes) %>%
  filter(abs(lab.x - x0) < 0.5 & lab.y > y0 & lab.y - y0 < 1 & perc.num > 30) %>% pull(id)
repeated_labels[ids, "lab.y"] <- repeated_labels[ids, "lab.y"] + 0.1

ids <- expand_grid(repeated_labels %>% mutate(id = row_number()), tree.nodes) %>%
  filter(abs(lab.x - x0) < 0.5 & y0 > lab.y & y0 - lab.y < 1) %>% pull(id)
repeated_labels[ids, "lab.y"] <- repeated_labels[ids, "lab.y"] - 0.1
ids <- expand_grid(repeated_labels %>% mutate(id = row_number()), tree.nodes) %>%
  filter(abs(lab.x - x0) < 0.5 & y0 > lab.y & y0 - lab.y < 1 & perc.num > 30) %>% pull(id)
repeated_labels[ids, "lab.y"] <- repeated_labels[ids, "lab.y"] - 0.1


events_order <- enframe(event_names, name = "event_id", "event_tag") %>%
  inner_join(arcs %>%
               inner_join(tree.labels,
                          by = c("child" = "node_id")) %>%
               group_by(label.y) %>%
               summarise(order1 = max(y), order2 = min(xend)) %>%
               arrange(desc(order1), order2) %>%
               mutate(y = -row_number(), event_order = row_number()),
             by = c("event_tag" = "label.y"))

events_order <- enframe(event_names, name = "event_id", "event_name") %>%
  inner_join(tree.nodes %>% select(x0_min, x0_max, clone) %>%
               separate_rows(clone, sep = ",") %>%
               filter(clone != "normal") %>%
               mutate(event_id = as.numeric(clone)) %>%
               arrange(x0_min, desc(x0_max)) %>%
               mutate(y = -row_number(), event_order = row_number()),
             by = c("event_id"))

events_order <- enframe(event_names, name = "event_id", "event_name") %>%
  inner_join(tree.nodes %>% select(x0_min, x0_max, clone) %>%
    separate_rows(clone, sep = ",") %>%
    filter(clone != "normal") %>%
    mutate(event_id = as.numeric(clone)) %>%
    arrange(x0_min, desc(x0_max)) %>%
    mutate(y = -row_number(), event_order = row_number()) %>%
    group_by(event_id) %>% summarise(events_order = min(event_order)) %>%
    arrange(events_order) %>%
    mutate(events_order = row_number(), y = -row_number()),
  by = c("event_id"))

dots <- tree.nodes %>% filter(leaf) %>%
  select(x0, lineage) %>% separate_rows(lineage, sep = ",") %>%
  mutate(lineage = as.numeric(lineage)) %>%
  left_join(events_order, by = c("lineage" = "event_id"))

connected_dots <- tree.nodes %>% select(x0_min, x0_max, clone) %>%
  separate_rows(clone, sep = ",") %>%
  filter(clone != "normal") %>%
  mutate(clone = as.numeric(clone)) %>%
  left_join(events_order, by = c("clone" = "event_id"))


colourCount = length(unique(tree.labels$label))

x <- col2rgb(RColorBrewer::brewer.pal(8, "Dark2"), alpha = FALSE)/255
x <- x + (0.5 - x) * 0
x <- rgb(x[1, ], x[2, ], x[3, ])
getPalette = colorRampPalette(x)

ggplot() +
  # geom_segment(data = arcs,
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              size = 5, col = "darkgrey", lineend = "round") +
  # geom_segment(data = arcs,
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              size = 4, col = "lightblue", lineend = "round") +
  geom_segment(data = colored.bent_arcs,
               aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = label),
               size = 4, lineend = "round") +
  geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
               aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = label),
               size = 4) +
  geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
               aes(x = arc.x, y = arc.y + 0.02, xend = arc.xend, yend = arc.y - 0.02),
               col = "black",
               size = 4) +

  ## Add nodes (grey or white if no cells)
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 4 + tree.nodes$perc.num / max(tree.nodes$perc.num) * 6,
             col = "black",
             fill = ifelse(tree.nodes$counts > 0, "lightgrey", "white"),
             shape = 21) +

  ## Add letters to mark parallel events
  geom_text(data = repeated_labels,
            aes(x = lab.x, y = lab.y, label = mark), col = "black", size = 4,
            fontface = "bold", vjust = 0.4, hjust = "middle") +
  geom_text(data = repeated_labels,
            aes(x = lab.x, y = lab.y, label = mark), col = "white", size = 3,
            fontface = "bold", vjust = 0.4, hjust = "middle") +
  
  ## Add empty matrix underneath the tree
  geom_point(data = tree.nodes %>% filter(leaf) %>% select(x0, lineage) %>% separate_rows(lineage, sep = ",") %>% expand(x0, lineage),
             aes(x = x0, y = -as.numeric(lineage)), col = "grey",
             size = 6) +
  geom_point(data = tree.nodes %>% filter(leaf) %>% select(x0, lineage) %>% separate_rows(lineage, sep = ",") %>% expand(x0, lineage),
             aes(x = x0, y = -as.numeric(lineage)), col = "white",
             size = 5) +

  scale_color_manual("CN envent", values = getPalette(colourCount), limits = event_names, guide = F) +
  
  ## Add colored dots and connecting lines for presence of events
  ## (connection lines mark shared lineage event)
  geom_segment(data = connected_dots,
               aes(x = x0_min, xend = x0_max, y = y, yend = y,
                   col = event_names[clone]),
               size = 2) +
  geom_point(data = dots,
             aes(x = x0, y = y,
                 col = event_names[as.numeric(lineage)]),
             size = 6) +

  geom_label(data = events_order %>%
               cbind(genes = cosmic_annotations$genes) %>%
               mutate(genes = as.character(genes)) %>%
               mutate(label = ifelse(str_length(genes) > 1 & str_length(genes) < 50,
                                     ifelse(grepl("/-", event_name, fixed = T),
                                            paste(genes, "(loss)"),
                                            paste(genes, "(gain)")),
                                     event_name)),
             aes(y = y,
                 col = event_name, label = label),
             x = max(tree.nodes$x0) + 1, hjust = 0, fill = "#EEEEEE",
             fontface = "bold",
             size = 3) +
  geom_segment(data = tree.nodes %>% filter(leaf),
               aes(x = x0, xend = x0, y = y0 - 0.5, yend = -0.5),
               size = 0.5, col = "black", linetype = "dotted") +

  coord_fixed(xlim = c(0, 32), ylim = c(-29, 19)) +
  labs(title = tree_title) +
  theme_void()


## ---- plot_final_trees_small ----

ggplot() +
  # geom_segment(data = arcs,
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              size = 5, col = "darkgrey", lineend = "round") +
  # geom_segment(data = arcs,
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              size = 4, col = "lightblue", lineend = "round") +
  geom_segment(data = colored.bent_arcs,
               aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = label),
               size = 4, lineend = "round") +
  geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
               aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = label),
               size = 4) +
  geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
               aes(x = arc.x, y = arc.y + 0.02, xend = arc.xend, yend = arc.y - 0.02),
               col = "black",
               size = 4) +
  
  ## Add nodes (grey or white if no cells)
  geom_point(data = tree.nodes,
             aes(x = x0, y = y0),
             size = 4 + tree.nodes$perc.num / max(tree.nodes$perc.num) * 6,
             col = "black",
             fill = ifelse(tree.nodes$counts > 0, "lightgrey", "white"),
             shape = 21) +
  
  ## Add letters to mark parallel events
  geom_text(data = repeated_labels,
            aes(x = lab.x, y = lab.y, label = mark), col = "black", size = 4,
            fontface = "bold", vjust = "center", hjust = "middle") +
  geom_text(data = repeated_labels,
            aes(x = lab.x, y = lab.y, label = mark), col = "white", size = 3,
            fontface = "bold", vjust = "center", hjust = "middle") +
  
  scale_color_manual("CN envent", values = getPalette(colourCount), limits = event_names, guide = F) +
  
  coord_fixed(xlim = c(0, 20), ylim = c(1, 19)) +
  labs(title = tree_title) +
  theme_void()
