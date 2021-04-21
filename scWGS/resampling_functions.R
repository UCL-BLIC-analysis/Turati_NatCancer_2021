get_sample_segcopy_data_from_signal_matrix <- function(signal_matrix, cns, thresholds) {
  if (ncol(signal_matrix) != length(cns)) {
    stop("The number of CN states does not match the number of events")
  }
  sample_segcopy_data <- matrix(2, nrow = nrow(signal_matrix), ncol = ncol(signal_matrix))
  sample_thresholds <- c()
  for (i in 1:length(cns)) {
    this_cn <- cns[i]
    for (x in 1:10) {
      cells <- c()
      if (this_cn < 2) {
        # Set the nominal threshold w.r.t cn unless otherwise specified
        nominal_threshold <- thresholds$mean[i]
        threshold_sd <- thresholds$sd[i]
        # Add some noise to the threshold but keep it whithin resonable limits
        nominal_threshold <- rnorm(1, nominal_threshold, threshold_sd)
        nominal_threshold <- min(max(nominal_threshold - 0.5, nominal_threshold), nominal_threshold + 0.5)
        # Save this value in the sample_threshold array
        sample_thresholds[i] <- nominal_threshold
        # Assign cells to this CN event based on the nominal threshold (with rnorm(sd = 0.1))
        cells <- signal_matrix[, i] < rnorm(nrow(signal_matrix), mean = nominal_threshold,
                                            sd = threshold_sd)
      }
      if (this_cn > 2) {
        # Set the nominal threshold w.r.t cn unless otherwise specified
        nominal_threshold <- thresholds$mean[i]
        threshold_sd <- thresholds$sd[i]
        # Add some noise to the threshold but keep it whithin resonable limits
        nominal_threshold <- rnorm(1, nominal_threshold, threshold_sd)
        nominal_threshold <- min(max(nominal_threshold - 0.5, nominal_threshold), nominal_threshold + 0.5)
        # Save this value in the sample_threshold array
        sample_thresholds[i] <- nominal_threshold
        # Assign cells to this CN event based on the nominal threshold (with rnorm(sd = 0.1))
        cells <- signal_matrix[, i] > rnorm(nrow(signal_matrix), mean = nominal_threshold,
                                            sd = threshold_sd)
      }
      # Modify the sample_segcopy_data matrix to reflect the cell harbouring this event
      if (sum(cells) > 0) {
        # print(paste(x, sum(cells)))
        sample_segcopy_data[cells, i] <- this_cn
        break
      }
    }
  }
  
  attr(sample_segcopy_data, "sample_thresholds") <- sample_thresholds
  return(sample_segcopy_data)
}

get_tree_from_data_matrix <- function(segcopy_data, plot = T, event_sizes) {
  unique_clones <- as_tibble(rbind(2, segcopy_data)) %>%
    mutate_all(function(x) {ifelse(x == 2, 0, 1)}) %>%
    group_by_all() %>% summarise(num.cells = n())
  num_cells <- unique_clones %>% pull(num.cells)
  num_cells[1] <- num_cells[1] - 1
  unique_clones <- unique_clones %>% select(-num.cells)
  arcs <- matrix(NA, nrow = 0, ncol = 3)

  if (missing(event_sizes)) {
    event_sizes <- rep(1, ncol(segcopy_data))
  }
  event_sizes <- tibble(event_size = event_sizes) %>% mutate(label = as.character(row_number()))
  for (s1 in 1:(nrow(unique_clones) - 1)) {
    for (s2 in (s1 + 1):nrow(unique_clones)) {
      if (all(unique_clones[s1, ] >= unique_clones[s2, ])) {
        dist <- sum(unique_clones[s1, ] - unique_clones[s2, ]) - 1 + 1 / (num_cells[s1] + 1)
        arcs <- rbind(arcs, c(s1, s2, dist))
      } else if (all(unique_clones[s2, ] >= unique_clones[s1, ])) {
        dist <- sum(unique_clones[s2, ] - unique_clones[s1, ]) - 1 + 1 / (num_cells[s1] + 1)
        arcs <- rbind(arcs, c(s1, s2, dist))
      }
    }
  }

  msa_tree <- getMinimumArborescence(1:nrow(unique_clones), arcs = arcs,
                                     show.graph = F, show.data = F)
  
  tree.arcs <- as.data.frame(msa_tree$tree.arcs)[, 1:2]
  colnames(tree.arcs) <- c("parent", "child")
  tree.arcs[, "label"] <- sapply(1:nrow(msa_tree$tree.arcs), function(row) {
    paste(which(unique_clones[tree.arcs[row, 2], ] - unique_clones[tree.arcs[row, 1], ] != 0),
          collapse = ",")
  })
    
  repeat{
    missing_links <- tree.arcs %>%
      group_by(parent) %>% mutate(n.children = n()) %>%
      separate_rows(label, sep = ",") %>%
      group_by(parent,label,n.children) %>% summarise(n = n(), children = paste(child, collapse = ",")) %>%
      filter(n > 1) %>%
      left_join(event_sizes, by = "label") %>%
      group_by(parent, n, n.children, children) %>%
      summarise(label = paste(label, collapse = ","), n.events = n(), total_size = sum(event_size)) %>%
      mutate(missing = n.children - n) %>%
      arrange(missing, desc(n.events), desc(total_size))
    if (nrow(missing_links) > 0) {
      this_subclone <- missing_links[1, ]
      new_node_id <- max(tree.arcs$child) + 1
      events <- strsplit(this_subclone$label, ",")[[1]]
      demoted_children <- tree.arcs %>%
        separate_rows(label, sep = ",") %>%
        filter(parent == this_subclone$parent & label %in% events) %>%
        pull(child) %>% unique()
      tree.arcs <- tree.arcs %>%
        separate_rows(label, sep = ",") %>%
        filter(!((parent == this_subclone$parent) & label %in% events)) %>%
        mutate(parent = ifelse(child %in% demoted_children, new_node_id, parent)) %>%
        rbind(tibble(parent = this_subclone$parent, child = new_node_id, label = events)) %>%
        group_by(parent, child) %>% summarise(label = paste(label, collapse = ",")) %>% ungroup
      num_cells[new_node_id] <- 0
      next
    } else {
      break
    }
  }
  tree.arcs <- as.data.frame(tree.arcs)
  attr(tree.arcs, "counts") <- tibble(n = num_cells)
  
  if (plot) {
    plot.tree(tree.arcs, counts = attr(tree.arcs, "counts"), plot.counts = T)
  }
  
  return(tree.arcs)
}

get_num_of_independent_events_from_tree <- function(tree) {
  if (!("label" %in% colnames(tree))) {
    stop("Tree needs to have labels of event")
  }
  sample_num_events <- tree %>%
    group_by(parent) %>%
    summarise(num_events = length(strsplit(paste(label, collapse = ","), ",")[[1]])) %>%
    pull(num_events) %>%
    sum()
  return(sample_num_events)
}


resample_matrix <- function(signal_matrix, cns,
                            num_loops = 50, num_samples = 1000, num_best = 100,
                            plot_intermediate_trees = F, use_num_events = T,
                            starting_thresholds) {
  c1 <- 0
  c2 <- 0
  best_sample_data <- array(0, dim = c(nrow(signal_matrix), ncol(signal_matrix), 0))
  best_sample_thresholds <- array(0, dim = c(ncol(signal_matrix), 0))
  best_sample_num_clones <- c()
  best_sample_num_events <- c()
  best_sample_rss <- c()
  
  # Current thresholds set to 0.5 closer to 2 than cns unless otherwise specified. This
  # would be:
  # * thr <- 0.5 if cns == 0
  # * thr <- 1.5 if cns == 1
  # * thr <- 2.5 if cns == 3
  # * thr <- 3.5 if cns == 4
  # * etc
  if (!missing(starting_thresholds)) {
    current_thresholds <- starting_thresholds
  } else {
    current_thresholds <- list(
      mean = ifelse(cns < 2, cns + 0.5, cns - 0.5),
      sd = 0.05)
  }
  
  # new_signal_matrix starts being the original ones but gets updated based on the
  # best solutions
  new_signal_matrix <- signal_matrix
  
  for (loop in 1:num_loops) {
    num_data_entries <- dim(best_sample_data)[3]
    if (num_data_entries > 0) {
      new_data <- array(0, dim = c(nrow(signal_matrix),
                                   ncol(signal_matrix),
                                   num_data_entries + num_samples))
      new_data[, , 1:num_data_entries] <- best_sample_data
      best_sample_data <- new_data
      rm(new_data)
      new_best_sample_thresholds <- array(0, dim = c(ncol(signal_matrix),
                                                     num_data_entries + num_samples))
      new_best_sample_thresholds[, 1:num_data_entries] <- best_sample_thresholds
      best_sample_thresholds <- new_best_sample_thresholds
      rm(new_best_sample_thresholds)
    } else {
      best_sample_data <- array(0, dim = c(nrow(signal_matrix),
                                           ncol(signal_matrix),
                                           num_samples))
      best_sample_thresholds <- array(0, dim = c(ncol(signal_matrix),
                                                 num_samples))
    }
    
    for (a in (num_data_entries + 1):(num_data_entries + num_samples)) {
      sample_segcopy_data <- get_sample_segcopy_data_from_signal_matrix(
        new_signal_matrix, cns, current_thresholds)
      sample_thresholds <- attr(sample_segcopy_data, "sample_thresholds")
      unique_clones <- as_tibble(rbind(2, sample_segcopy_data)) %>%
        mutate_all(function(x) {ifelse(x == 2, 0, 1)}) %>%
        group_by_all() %>% summarise(num.cells = n())
      num_cells <- unique_clones %>% pull(num.cells)
      num_cells[1] <- num_cells[1] - 1
      # print(nrow(unique_clones))
      dim(best_sample_data)
      best_sample_data[, , a] <- sample_segcopy_data
      best_sample_thresholds[, a] <- sample_thresholds
      sample_num_clones = nrow(unique_clones)
      sample_rss <- sum((as.vector(signal_matrix) - as.vector(sample_segcopy_data)) ^ 2)
      # print(colSums(unique_clones))
      # print(any(colSums(unique_clones) == 0))
      if (any(colSums(unique_clones) == 0)) {
        best_sample_num_clones <- c(best_sample_num_clones, NA)
        best_sample_rss <- c(best_sample_rss, NA)
        best_sample_num_events <- c(best_sample_num_events, NA)
      } else {
        best_sample_num_clones <- c(best_sample_num_clones, sample_num_clones)
        best_sample_rss <- c(best_sample_rss, sample_rss)
        best_ones <- order(best_sample_num_clones, best_sample_rss)[1:num_best]
        # best_sample_num_events <- c(best_sample_num_events, NA)
        if (use_num_events &
            (sample_num_clones < 50) &
            (sample_num_clones <= max(best_sample_num_clones[best_ones], na.rm = T))) {
          sample_tree <- get_tree_from_data_matrix(sample_segcopy_data, plot = F)
          sample_num_events <- get_num_of_independent_events_from_tree(sample_tree)
          best_sample_num_events <- c(best_sample_num_events, sample_num_events)
          c1 <- c1 + 1
        } else {
          best_sample_num_events <- c(best_sample_num_events, 1)
          c2 <- c2 + 1
        }
      }
    }
    
    # table(l)
    best_ones <- order(best_sample_num_clones, best_sample_num_events, best_sample_rss)[1:num_best]
    best_sample_data       <- best_sample_data[, , best_ones]
    best_sample_thresholds <- best_sample_thresholds[, best_ones]
    best_sample_num_events <- best_sample_num_events[best_ones]
    best_sample_num_clones <- best_sample_num_clones[best_ones]
    best_sample_rss        <- best_sample_rss[best_ones]
    
    resampled_segcopy_data <- round(apply(best_sample_data, 1:2, mean))
    num_wobbly <- sum(!apply(best_sample_data, 1:2, function(x){all(x == x[1])}))
    events_table <- as_tibble(rbind(2, resampled_segcopy_data)) %>%
      mutate_all(function(x) {ifelse(x == 2, 0, 1)})
    unique_clones <- events_table %>%
      group_by_all() %>% summarise(num.cells = n())
    
    print(paste(loop, num_wobbly, nrow(unique_clones), best_sample_num_clones[1], sum(events_table)))
    if (use_num_events) {
      print(paste(c1, c2))
    }
    if (plot_intermediate_trees) {
      get_tree_from_data_matrix(resampled_segcopy_data, plot = T)
    }
    
    new_signal_matrix <- signal_matrix
    for (x in 1:nrow(best_sample_data)) {
      for (y in 1:ncol(best_sample_data)) {
        if (all(best_sample_data[x, y, ] == best_sample_data[x, y, 1])) {
          new_signal_matrix[x, y] <- best_sample_data[x, y, 1]
        }
      }
    }
    if (loop > 10 & num_wobbly == 0) {
      break
    }
    # current_thresholds <- rowMeans(best_sample_thresholds)
  }
  
  attr(best_sample_data, "thresholds") <- best_sample_thresholds
  
  return(best_sample_data)
}


