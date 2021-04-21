# tree is a tibble of (node_id, parent_id, events, size)
plot_cna_tree <- function(tree, tree.arcs, event_labels, title,
                          show_parallel_events = T, show_cna_matrix = T, show_sizes = F,
                          show_event_labels = F, color_by_clone = T,
                          highlights,
                          label_width = 20, x_left = 0, x_right, y_top, y_bottom,
                          x_step = 1.5, y_step = 1.5,
                          min_node_size = 4, max_node_size = 10, branch_size = 4,
                          label_size = 2) {
  if (missing(tree) & !missing(tree.arcs)) {
    tree <- tree.arcs %>%
      rbind(tibble(parent = as.numeric(NA), child = 1, label = "")) %>%
      transmute(node_id = child, parent_id = parent, events = label) %>%
      arrange(node_id) %>%
      mutate(size = attr(tree.arcs, "counts")[, 1, drop = T])
  }
  
  add_tree_coordinates <- function(root_id = 1, x0 = 0, depth = 0) {
    children_node_ids <- tree.nodes %>% filter(parent_id == root_id) %>% pull(node_id)
    events <- tree.nodes %>% filter(parent_id == root_id) %>% pull("events")
    for (this_child_node_id in children_node_ids) {
      num.events <- tree.nodes %>% filter(node_id == this_child_node_id) %>%
        pull(events) %>% strsplit(",") %>% .[[1]] %>% length()
      x0 <- add_tree_coordinates(this_child_node_id, x0, depth + num.events)
    }
    row_id <- which(tree.nodes$node_id == root_id)
    if (length(children_node_ids) == 0) {
      tree.nodes[row_id, c("x0", "x0_min", "x0_max", "depth", "is.leaf")] <<-
        tibble(x0 = x0, x0_min = x0, x0_max = x0, depth = depth, is.leaf = TRUE)
      x0 <- x0 + 1.5
    } else {
      r <- tree.nodes %>% filter(node_id %in% children_node_ids) %>% pull(x0) %>% range
      this_x0 <- (r[1] + r[2]) / 2
      r <- tree.nodes %>% filter(node_id %in% children_node_ids) %>% select(x0_min, x0_max) %>% range
      this_x0_min <- r[1]
      this_x0_max <- r[2]
      tree.nodes[row_id, c("x0", "x0_min", "x0_max", "depth", "is.leaf")] <<-
        tibble(x0 = this_x0, x0_min = this_x0_min, x0_max = this_x0_max, depth = depth, is.leaf = FALSE)
    }
    return(x0)
  }

  root_id <- tree %>% filter(is.na(parent_id)) %>% pull(node_id)
  if (length(root_id) < 1) {
    stop("No root found on the tree. One of the nodes should have parent_id == NA.")
  } else if (length(root_id) > 1) {
    stop("More than one root found on the tree. This method can only plot one tree at a time.")
  }

  tree.nodes <- tree %>%
    mutate(x0 = 0, x0_min = 0, x0_max = 0, depth = 0, is.leaf = F)

  tree.labels <- tibble()[0, ]
  
  # print(tree.nodes)
  add_tree_coordinates(root_id = root_id)
  
  if (missing(y_top)) {
    y_top <- max(tree.nodes$depth) * y_step + 1
  }
  tree.nodes <- tree.nodes %>% mutate(y0 = y_top - depth * y_step)
  # print(tree.nodes)

  tree.plot <- ggplot()

  ## =====================================================================================
  ## BRANCHES and NODES
  ## =====================================================================================
  ## tree arcs is a tibble with node_id, parent_id, x, xend, y, yend, event, event_num
  ## where the branches of the tree are split into events; event_num is the order of that event
  ## in that branch
  tree.arcs <- tree.nodes %>%
    select(node_id, parent_id, events, xend = x0) %>%
    inner_join(tree.nodes %>% select(node_id, x = x0, y = y0),
               by = c("parent_id" = "node_id")) %>%
    select(node_id, parent_id, events, x, xend, y) %>%
    ## Split the events by row and assing an event number to tell first from others
    separate_rows(events) %>% rename(event_id = events) %>%
    mutate(event_id = as.numeric(event_id)) %>%
    group_by(node_id) %>% mutate(event_num = row_number()) %>% ungroup() %>%
    ## Adjust y and yend for each arc
    mutate(y = y - (event_num - 1) * y_step, yend = y - y_step)

  if (missing(event_labels)) {
    event_labels <- paste0("CNA-", unique(tree.arcs$event_id))
    label_width <- 3
  }
  
  ## Add highlighting if requested
  if (!missing(highlights)) {
    highlighted_arcs <- tibble(.rows = 0)
    for (i in 1:length(highlights)) {
      full_arc <- tree.arcs %>% filter((node_id %in% highlights[[i]]) & ((parent_id %in% highlights[[i]]) | (event_num != 1)))
      partial_arc <- tree.arcs %>% filter(node_id %in% highlights[[i]] & !(parent_id %in% highlights[[i]]) & event_num == 1)
      if (nrow(full_arc) > 0) {
        highlighted_arcs <- rbind(
          highlighted_arcs,
          cbind(color = names(highlights)[i],
                mode = "full",
                full_arc))
      }
      if (nrow(partial_arc) > 0) {
        highlighted_arcs <- rbind(
          highlighted_arcs,
          cbind(color = names(highlights)[i],
                mode = "partial",
                partial_arc))
      }
    }
    # print(highlighted_arcs)
    
    # Split the first arc after a branch to create a bent and then...
    highlighted.bent_arcs <- highlighted_arcs %>%
      filter(event_num == 1) %>%
      filter(mode == "full") %>%
      expand_grid(segment = 1:2) %>%
      arrange(y, desc(abs(x - xend)), segment) %>%
      mutate(arc.x = ifelse(segment == 1, x, xend),
             arc.xend = xend,
             arc.y = ifelse(segment == 1, y, y - 0.3),
             arc.yend = ifelse(segment == 1, y - 0.3, y - 1.5)) %>%
      select(arc.x, arc.xend, arc.y, arc.yend, color)
    
    highlighted.bent_arcs <- rbind(
      highlighted.bent_arcs,
      highlighted_arcs %>%
        filter(event_num == 1) %>%
        filter(mode == "partial") %>%
        expand_grid(segment = 1:2) %>%
        arrange(y, desc(abs(x - xend)), segment) %>%
        mutate(arc.x = ifelse(segment == 1, x + (xend - x) * 0.8, xend),
               arc.xend = xend,
               arc.y = ifelse(segment == 1, y - 0.3 * 0.8, y - 0.3),
               arc.yend = ifelse(segment == 1, y - 0.3, y - 1.5)) %>%
        select(arc.x, arc.xend, arc.y, arc.yend, color)
    )
    
    # Continue with straight vertical lines
    highlighted.straight_arcs <- highlighted_arcs %>%
      filter(event_num > 1) %>%
      mutate(arc.x = xend,
             arc.xend = xend,
             arc.y = y,
             arc.yend = yend) %>%
      select(arc.x, arc.xend, arc.y, arc.yend, color)

    tree.plot <- tree.plot +
      geom_segment(data = highlighted.bent_arcs,
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = highlighted.bent_arcs$color,
                   size = branch_size + 8, lineend = "round") +
      geom_segment(data = highlighted.straight_arcs,
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = highlighted.straight_arcs$color,
                   size = branch_size + 8, lineend = "round") #+
      # geom_segment(data = highlighted.bent_arcs,
      #              aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
      #              col = "white",
      #              size = branch_size + 4, lineend = "round") +
      # geom_segment(data = highlighted.straight_arcs,
      #              aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
      #              col = "white",
      #              size = branch_size + 4, lineend = "round")
  }


  
  # Split the first arc after a branch to create a bent and then...
  colored.bent_arcs <- tree.arcs %>%
    filter(event_num == 1) %>%
    expand_grid(segment = 1:2) %>%
    arrange(y, desc(abs(x - xend)), segment) %>%
    mutate(arc.x = ifelse(segment == 1, x, xend),
           arc.xend = xend,
           arc.y = ifelse(segment == 1, y, y - 0.3),
           arc.yend = ifelse(segment == 1, y - 0.3, y - 1.5)) %>%
    select(arc.x, arc.xend, arc.y, arc.yend, event_id)

  # Continue with straight vertical lines
  colored.straight_arcs <- tree.arcs %>%
    filter(event_num > 1) %>%
    mutate(arc.x = xend,
           arc.xend = xend,
           arc.y = y,
           arc.yend = yend) %>%
    select(arc.x, arc.xend, arc.y, arc.yend, event_id)
  
  if (color_by_clone) {
    tree.plot <- tree.plot +
      geom_segment(data = colored.bent_arcs,
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = event_labels[event_id]),
                   size = branch_size, lineend = "round") +
      geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend, col = event_labels[event_id]),
                   size = branch_size) +
      geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
                   aes(x = arc.x, y = arc.y + 0.02, xend = arc.xend, yend = arc.y - 0.02),
                   col = "black",
                   size = branch_size)
  } else {
    tree.plot <- tree.plot +
      geom_segment(data = colored.bent_arcs,
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = "darkgrey",
                   size = branch_size, lineend = "round") +
      geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = "darkgrey",
                   size = branch_size) +
      geom_segment(data = colored.bent_arcs,
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = "lightblue",
                   size = branch_size - 1, lineend = "round") +
      geom_segment(data = colored.straight_arcs %>% arrange(arc.yend),
                   aes(x = arc.x, y = arc.y, xend = arc.xend, yend = arc.yend),
                   col = "lightblue",
                   size = branch_size - 1)
  }
  ## Add nodes (grey or white if no cells)
  tree.plot <- tree.plot +
    geom_point(data = tree.nodes,
               aes(x = x0, y = y0),
               size = min_node_size + tree.nodes$size / max(tree.nodes$size) * (max_node_size - min_node_size),
               col = "black",
               fill = ifelse(tree.nodes$size > 0, "lightgrey", "white"),
               shape = 21)
  if (show_sizes) {
    tree.plot <- tree.plot +
      geom_text(data = tree.nodes %>% filter(size > 0),
                aes(x = x0, y = y0, label = size),
                size = 4)
  }
  ## =====================================================================================


  ## =====================================================================================
  ## Plot labels on the tree unless they will be printed on the CNA matrix
  ## =====================================================================================
  if (!show_cna_matrix & show_event_labels) {
    labels.color <- ifelse(grepl("(loss)", event_labels[tree.arcs$event_id], fixed = T), "blue", "red")
    labels.text <- sub(" \\((loss|gain)\\)", "", event_labels[tree.arcs$event_id])
    labels.text <- gsub("([^/]+/[^/]+/[^/]+/[^/]+/)", "\\1\n", labels.text)
    tree.plot <- tree.plot +
      ggrepel::geom_label_repel(data = tree.arcs,
                                aes(x = xend, y = yend + as.numeric(y_step) / 2),
                                label = labels.text,
                                hjust = 0.5,
                                size = label_size,
                                col = labels.color,
                                min.segment.length = 0,
                                point.padding = NA, seed = 1)
  }
  ## =====================================================================================
  
  
  ## =====================================================================================
  ## COLOR PALETTE FOR THE EVENTS
  ## =====================================================================================
  ## Create a palette based on a toned down version of the Dark2 brewer palette for the events
  colourCount = length(unique(tree.arcs$event_id))
  
  x <- col2rgb(RColorBrewer::brewer.pal(8, "Dark2"), alpha = FALSE)/255
  x <- x + (0.5 - x) * 0
  x <- rgb(x[1, ], x[2, ], x[3, ])
  getPalette = colorRampPalette(x)
  
  tree.plot <- tree.plot +
    scale_color_manual("CN envent", values = getPalette(colourCount), guide = F)
  ## =====================================================================================
  
  
  ## =====================================================================================
  ## LABEL REPEATED (PARALLEL) EVENTS
  ## =====================================================================================
  ## Get a list of repeated events from the tree.arcs tibble and "mark" them as a, b, c...
  if (show_parallel_events) {
    
    labels_for_parallel_events <- tree.arcs %>% group_by(event_id) %>%
      filter(n() > 1) %>% ungroup() %>%
      transmute(event_id, lab.x = xend, lab.y = (y + yend) / 2,
                mark = letters[as.numeric(factor(event_id))])
    # print(labels_for_parallel_events)
    
    ## Nudge these labels up or down depending on the presence of a node or large node above or below the label
    table_of_repeated_labels_and_nodes <- expand_grid(labels_for_parallel_events %>% mutate(id = row_number()),
                                                      tree.nodes)
    ## Nudge labels up
    ids <- table_of_repeated_labels_and_nodes %>%
      filter(abs(lab.x - x0) < 0.5 & lab.y > y0 & lab.y - y0 < 1) %>%
      pull(id)
    labels_for_parallel_events[ids, "lab.y"] <- labels_for_parallel_events[ids, "lab.y"] + 0.1
    
    ## Nudge labels up further if large node
    ids <- table_of_repeated_labels_and_nodes %>%
      filter(abs(lab.x - x0) < 0.5 & lab.y > y0 & lab.y - y0 < 1) %>%
      filter(size > 0.5 * max(tree.nodes$size)) %>%
      pull(id)
    labels_for_parallel_events[ids, "lab.y"] <- labels_for_parallel_events[ids, "lab.y"] + 0.1
    
    ## Nudge labels down
    ids <- table_of_repeated_labels_and_nodes %>%
      filter(abs(lab.x - x0) < 0.5 & y0 > lab.y & y0 - lab.y < 1) %>%
      pull(id)
    labels_for_parallel_events[ids, "lab.y"] <- labels_for_parallel_events[ids, "lab.y"] - 0.1
    
    ## Nudge labels down further if large node
    ids <- table_of_repeated_labels_and_nodes %>%
      filter(abs(lab.x - x0) < 0.5 & y0 > lab.y & y0 - lab.y < 1) %>%
      filter(size > 0.5 * max(tree.nodes$size)) %>%
      pull(id)
    labels_for_parallel_events[ids, "lab.y"] <- labels_for_parallel_events[ids, "lab.y"] - 0.1
    
    tree.plot <- tree.plot +
      geom_text(data = labels_for_parallel_events,
                aes(x = lab.x, y = lab.y, label = mark), col = "black", size = 4,
                fontface = "bold", vjust = 0.4, hjust = "middle") +
      geom_text(data = labels_for_parallel_events,
                aes(x = lab.x, y = lab.y, label = mark), col = "white", size = 3,
                fontface = "bold", vjust = 0.4, hjust = "middle")
  }
  ## =====================================================================================
  

  ## =====================================================================================
  ## CN EVENT MATRIX
  ## =====================================================================================
  if (show_cna_matrix) {
    events_order <- tree.arcs %>%
      left_join(tree.nodes, by = c("node_id" = "node_id")) %>%
      arrange(x0_min, desc(x0_max), desc(y)) %>%
      mutate(event_order = row_number()) %>%
      group_by(event_id) %>%
      summarise(events_order = min(event_order)) %>%
      arrange(events_order) %>%
      mutate(events_order = row_number())
    # print(events_order)
    
    dots <- expand_grid(tree.nodes %>% select(events, x0_min, x0_max),
                        tree.nodes %>% filter(is.leaf) %>% select(x0)) %>%
      filter(x0_min <= x0 & x0 <= x0_max) %>%
      separate_rows(events) %>%
      transmute(x0, event_id = as.numeric(events)) %>%
      left_join(events_order, by = c("event_id")) %>%
      mutate(y = -events_order) %>%
      filter(!is.na(y)) %>% select(x0, y, event_id)
    
    connected_dots <- tree.nodes %>% select(x0_min, x0_max, events) %>%
      separate_rows(events, sep = ",") %>%
      filter(events != "normal") %>%
      mutate(event_id = as.numeric(events)) %>%
      left_join(events_order, by = c("event_id")) %>%
      mutate(y = -events_order) %>%
      filter(!is.na(y)) %>%
      as_tibble()
    
    ## Add empty matrix underneath the tree
    tree.plot <- tree.plot +
      geom_point(data = expand_grid(y = -unique(events_order$events_order),
                                    x = tree.nodes %>% filter(is.leaf) %>% pull(x0) %>% unique()),
                 aes(x = x, y = y), col = "grey",
                 size = 6) +
      geom_point(data = expand_grid(y = -unique(events_order$events_order),
                                    x = tree.nodes %>% filter(is.leaf) %>% pull(x0) %>% unique()),
                 aes(x = x, y = y), col = "white",
                 size = 5) +
      geom_segment(data = connected_dots,
                   aes(x = x0_min, xend = x0_max, y = y, yend = y, col = event_labels[event_id]),
                   size = 2) +
      geom_point(data = dots,
                 aes(x = x0, y = y, col = event_labels[event_id]),
                 size = 6) +
      geom_label(data = events_order %>%
                   mutate(y = -events_order),
                 aes(y = y,
                     col = event_labels[event_id],
                     label = event_labels[event_id]),
                 x = max(tree.nodes$x0) + 1, hjust = 0, fill = "#EEEEEE",
                 fontface = "bold",
                 size = 3) +
      
      geom_segment(data = tree.nodes %>% filter(is.leaf),
                   aes(x = x0, xend = x0, y = y0 - 0.5, yend = -0.5),
                   size = 0.5, col = "black", linetype = "dotted")
    
    if (missing(y_bottom)) {
      y_bottom <- -length(event_labels)
    }
    if (missing(x_right)) {
      x_right <- max(tree.nodes$x0) + 1 + label_width
    }
    
  } else {
    if (missing(y_bottom)) {
      y_bottom <- 0
    }
    if (missing(x_right)) {
      x_right <- max(tree.nodes$x0)
    }
  }
  ## =====================================================================================

  tree.plot <- tree.plot +
    coord_fixed(xlim = c(x_left, x_right),
                ylim = c(y_bottom, y_top),
                clip = "off")

  if (!missing(title)) {
    tree.plot <- tree.plot +
      labs(title = title)
  }

  tree.plot <- tree.plot + theme_void()

  return(tree.plot)
}


# source("../resampling/resampling_functions.R")
# matrix <- matrix(sample(1:2, size = 100, replace = T), nrow = 10)
# tree.arcs <- get_tree_from_data_matrix(matrix, plot = F)
# p <- plot_cna_tree(tree.arcs = tree.arcs)
# # print(p + theme_classic())
