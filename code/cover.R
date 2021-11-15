suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
library(graphics)
library(RColorBrewer)
suppressWarnings(suppressPackageStartupMessages(library(scales)))

set.seed(19740827)

a <- 12.5
b <- 2.5
l = 100
x <- c((0:l)^2/(l^2) * a - a, a - ((l-1):0)^2/(l^2) * a)
y <- sqrt(b^2 * (1 - x^2 / a^2))

ellipse_template <- tibble(x = c(x,rev(x)[-1]), y = c(y, -y[-1])) %>%
  mutate(order = row_number(), front = y < 0, back = y > 0)
ellipse_template <- rbind(
  cbind(ellipse_template[1:(nrow(ellipse_template) -1), ],
        ellipse_template[-1, ] %>% transmute(xend = x, yend = y)),
  cbind(ellipse_template[nrow(ellipse_template) - 1, ],
        ellipse_template[1, ] %>% transmute(xend = x, yend = y))
)

tornado_bottom <- 0
tornado_top <- 300
x_bias <- tibble(h = seq(tornado_bottom, tornado_top, length.out = 11),
                 x0 = c(0, 0, 30, 45, 20, 36, 45, 80, 120, 60, 75))
x0.loess <- loess(x0 ~ h, x_bias)

grey_ellipses <- tibble()
for (i in 1:300) {
  height <- runif(1, min = tornado_bottom, max = tornado_top)
  # height <- tornado_bottom + rbeta(1, 1.5, 0.8) * (tornado_top - tornado_bottom)
  delta_y <- rnorm(1) - height * 0.05
  delta_x <- rnorm(1, sd = height / 20) + predict(x0.loess, height) + height * 0.1
  grey_ellipses <- grey_ellipses %>%
    rbind(ellipse_template %>%
            mutate(ellipse = i,
                   height = height,
                   size = y * (0.7 + 0.005 * height ^ 1.3),
                   y = y * (0.7 + 0.005 * height ^ 1.3) + height + delta_y,
                   x = x * (0.7 + 0.005 * height ^ 1.3) + delta_x,
                   col = "black"))
}

## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
## Color functions ========================================
## XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

lighter_col <- function(col, factor = 0, alpha = 1) {
  l <- farver::decode_colour(col, to ="hcl")[, "l"]
  col2hcl(col, l = l + min(10, (100 - l)) * factor, alpha = alpha)
}

new_col <- function(base_col, col_factor) {
  return(base_col)
  col2hcl(base_col,
          c = farver::decode_colour(base_col, to = "hcl")[, "c"] *
            (0.3 + col_factor * (1 - 0.3)),
          l = farver::decode_colour(base_col, to = "hcl")[, "l"] *
            (0.5 + col_factor * (1 - 0.5)))
}

## YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY


get_ellipses_arcs <- function(ellipse_template, num, base_colors = c("red", "blue", "green4"), tornado_bottom = 0, tornado_top = 100, xbias = x_bias) {
  x0.loess <- loess(x0 ~ h, x_bias)
  ellipses_arcs <- tibble()
  for (i in 1:num) {
    # Pick a random heigth from tornado bottom to tornado_top
    height <- runif(1, min = tornado_bottom, max = tornado_top)
    height <- tornado_bottom + rbeta(1, 1.1, 0.8) * (tornado_top - tornado_bottom)
    
    # Modify the size of the ellipse
    size_factor <- rnorm(1, 1, 0.01)
    
    # Move the ellipse right/left and top/bottom
    delta_y <- rnorm(1)
    delta_x <- rnorm(1, sd = height / 20) + predict(x0.loess, height)
    
    # Pick segment of the ellipse (l ~ 1/4 of the full ellipse)
    start <- sample(ellipse_template$order, 1)
    length <- rnorm(1, mean = l , sd = l / 10)
    
    # Pick a color depending on the height (favour first color at the bottom)
    height_ratio <- (height - tornado_bottom) / (tornado_top - tornado_bottom)
    base_col <- sample(base_colors, 1, prob = dbeta(height_ratio, 1:length(base_colors), 1))
    
    # Color factor (NOT USED ANYMORE)
    col_factor <- rbeta(1, 1.5, 0.9) * (height - tornado_bottom) / (tornado_top - tornado_bottom)
    col_factor <- runif(1, 0, 1) * (height - tornado_bottom) / (tornado_top - tornado_bottom)
    ellipse_col <- new_col(base_col, col_factor)
    # print(paste(col_factor, ellipse_col))
    this_ellipse_arc <- ellipse_template %>%
      mutate(ellipse = i,
             height = height + delta_y,
             size = y * (0.7 + 0.005 * (size_factor * (height - delta_y)) ^ 1.3),
             y = y * (0.7 + 0.005 * (size_factor * (height - delta_y)) ^ 1.3) + height + delta_y,
             x = x * (0.7 + 0.005 * (size_factor * (height - delta_y)) ^ 1.3) + delta_x,
             yend = yend * (0.7 + 0.005 * (size_factor * (height - delta_y)) ^ 1.3) + height + delta_y,
             xend = xend * (0.7 + 0.005 * (size_factor * (height - delta_y)) ^ 1.3) + delta_x,
             col_factor = col_factor,
             col = ellipse_col)
    
    if (start + length > max(ellipse_template$order)) {
      this_ellipse_arc_data <-
        rbind(
          this_ellipse_arc %>% filter(order >= start),
          this_ellipse_arc %>% filter(order < start + length - max(ellipse_template$order))
        )
    } else {
      this_ellipse_arc_data <- this_ellipse_arc %>%
        filter(order >= start & order < start + length)
    }
    
    ellipses_arcs <- ellipses_arcs %>%
      rbind(this_ellipse_arc_data)
  }
  return(ellipses_arcs)
}

# Color scheme
colors <- c("red", "blue", "green4", "darkorange", "yellow2", "purple2","brown3")
colors <- brewer.pal(9, "Set1")
colors <- c("red", "blue", "green4", "darkorange", "yellow2", "purple2","brown3")
colors <- c("red2", "blue3", "darkgreen", "yellow1")
colors <- brewer.pal(9, "PuBuGn")
colors <- brewer.pal(11, "RdBu")[-c(6)]
colors <- c("#25413C", "#8E5936", "#A69384", "#859866", "#53978F")
colors <- c("#344F59", "#648C8C", "#D2E4D8", "#F3EDDD", "#CB7E54")
colors <- c("#BCBC45", "#7F7F7F", "#58BCCC", "#8E6AB8",
            #"#D57EBF",
            "#529D3E", "#EF8535", "#C53932", "#3B76AF", "#84584E")
colors <- c("#86281C", "#D98977", "#297373", "#D9580D", "#1E91D9")
colors <- c("#2A4759", "#658CBF", "#F2AF5C", "#F2D5CE", "#D98F8F")
colors <- c("#592525", "#03658C", "#595332", "#A65656", "#8C845D")
colors <- c("#382F73", "#C1D0D9", "#3A4001", "#403001", "#732C02")
colors <- c("#253C59", "#5A798C", "#95B3BF", "#BFA004", "#BF8D30")
colors <- c("#557570", "#D6DCD8", "#D39A89", "#A2543D", "#CE7E5D")
colors <- c("#F2F2F2", "#0F2610", "#90A686", "#567339", "#A4A66A")
colors <- c("#47594A", "#F2F2F0", "#A6A381", "#F2E0D0", "#A67B77")
colors <- c("#590202", "#344F59", "#648C8C", "#D2E4D8", "#CB7E54", "#F2E0D0")
hires = T

ellipses_arcs <- get_ellipses_arcs(ellipse_template, num = 2000, base_colors = colors[-1],
                                   tornado_bottom = tornado_bottom,
                                   tornado_top = tornado_top)
# ellipses_arcs <- ellipses_arcs %>% arrange(height, ellipse, col)
# saveRDS(ellipses_arcs, file = "selected_ellipses.20210410.rds")

ellipses_arcs3D <- ellipses_arcs
ellipses_arcs3D <- ellipses_arcs3D %>%
  mutate(line_width = 0.5 + (size - max(size)) / (min(size) - max(size)) * 2.5)
if (hires) {
  ellipses_arcs3D <- rbind(ellipses_arcs3D %>%
                             mutate(line_width = line_width * 1, col = lighter_col(col, -0.1)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.95, col = lighter_col(col, 0)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.9, col = lighter_col(col, 0.13)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.85, col = lighter_col(col, 0.16)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.8, col = lighter_col(col, 0.2)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.75, col = lighter_col(col, 0.3)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.7, col = lighter_col(col, 0.4)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.65, col = lighter_col(col, 0.5)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.6, col = lighter_col(col, 0.6)))
} else {
  ellipses_arcs3D <- rbind(ellipses_arcs3D %>%
                             mutate(line_width = line_width * 1, col = lighter_col(col, -0.1)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.95, col = lighter_col(col, 0)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.9, col = lighter_col(col, 0.13)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.85, col = lighter_col(col, 0.16)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.8, col = lighter_col(col, 0.2)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.75, col = lighter_col(col, 0.3)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.7, col = lighter_col(col, 0.4)),
                           # ellipses_arcs3D %>%
                           #   mutate(line_width = line_width * 0.65, col = lighter_col(col, 0.5)),
                           ellipses_arcs3D %>%
                             mutate(line_width = line_width * 0.6, col = lighter_col(col, 0.6)))
}

ellipses_arcs3D <- ellipses_arcs3D %>% arrange(height, ellipse)


# g <- ggplot() +
#   geom_path(data = all_ellipses,
#             aes(x = x, y = y, group = ellipse, size = size),
#             color = "grey",
#             alpha = 0.3, linejoin = "mitre", show.legend = F) +
#   geom_path(data = selected_ellipses %>% filter(!front),
#             aes(x = x, y = y, group = ellipse, size = size),
#             col = selected_ellipses %>% filter(!front) %>% pull(col),
#             alpha = 1, linejoin = "bevel", lineend = "round", show.legend = F) +
#   geom_path(data = selected_ellipses %>% filter(!back),
#             aes(x = x, y = y, group = ellipse, size = size),
#             col = selected_ellipses %>% filter(!back) %>% pull(col),
#             alpha = 1, linejoin = "bevel", lineend = "round", show.legend = F)
# g + coord_fixed() + theme_void() + scale_size(range = c(3, 0.3)) +
#   theme(plot.background = element_rect(fill = rgb(0.01, 0.05, 0.15)))



#
png("tmp.png")
pmat <- persp(1:50, 1:30, matrix(0, nrow = 50, ncol = 30), zlim = c(0, 0.1), theta = 5, phi = 10)
dev.off()
unlink("tmp.png")
point.size <- function (x, y, z, pmat, delta = 0.1) {
  low <- trans3d(x-delta, y-delta, z-delta, pmat)
  high <- trans3d(x+delta, y+delta, z+delta, pmat)
  avg = 50000*(abs(high$x-low$x) + abs(high$y-low$y))^3
  return(avg)
}


cells <- tibble(orig_x = rep(1:50, times = 30), orig_y = rep(1:30, each = 50))
cells$orig_x <- cells$orig_x + 0.5 * (cells$orig_y %% 2)
cells$orig_x <- cells$orig_x + rnorm(nrow(cells), 0, 0.15)
cells$orig_y <- cells$orig_y + rnorm(nrow(cells), 0, 0.15)
t <- trans3d(cells$orig_x, cells$orig_y, 0, pmat)
cells$x <- t$x
cells$y <- t$y
cells$size = point.size(cells$x, cells$y, 0, pmat)
cells$size <- 2 + 2 * (cells$size - range(cells$size)[1]) / diff(range(cells$size))
cells$x <- cells$x * 300*50/40 + 35
cells$y <- cells$y * 250*30/30 + 70
cells$col_factor <- runif(nrow(cells), 0, 1)
cells$col <- new_col(sample(colors, nrow(cells), replace = T, prob = c(0.5, rep(2, length(colors) - 1))),
                     col_factor = cells$col_factor)


background <- tibble()
t <- trans3d(x = c(0, 0, 51, 51), y = c(0, 31, 31, 0), 0, pmat)
this_area <- tibble(x = t$x, y = t$y,
                    group = 0,
                    color = "lightgrey")
this_area$x <- this_area$x * 300*50/40 + 35
this_area$y <- this_area$y * 250*30/30 + 70
background <- rbind(background, this_area)

for (i in 1:100) {
  d <- car::ellipse(c(5 + 0.1 * i, 0.1 * i), matrix(c(2,1,1,2), nrow = 2), 12 - 0.1 * i, draw = F)
  d[d[, 1] < 0, 1] <- 0
  d[d[, 2] < 0, 2] <- 0
  d[d[, 1] > 41, 1] <- 41
  d[d[, 2] > 31, 2] <- 31
  t <- trans3d(x = d[, 1], y = d[,2 ], 0, pmat)
  this_area <- tibble(x = t$x, y = t$y,
                      group = i,
                      color = muted("orange", c = i / 2, l = 85 - 0.25 * i ))
  this_area$x <- this_area$x * 300*50/40 + 35
  this_area$y <- this_area$y * 250*30/30 + 70
  background <- rbind(background, this_area)
}


add_3D_effect <- function(cells, num_steps = 15, alpha = 1) {
  # Outer color is a tad darker
  cells3D <- cells %>% mutate(color = lighter_col(col, -0.30 * 3, alpha = alpha))
  # cells3D <- cells %>% mutate(color = lighter_col(col, -0.30 * (1 + (size - 2) * 0.3), alpha = alpha))
  for (i in 1:num_steps) {
    step <- 0.75 / num_steps
    size_factor <- 1 - step * i # From 1 (skipped) to 0.25
    lighter_factor <- -0.25 + step * i # From -0.25 (skipped) to 0.50
    cells3D <- rbind(cells3D,
                     cells %>% mutate(color = lighter_col(col, lighter_factor * 3, alpha = alpha),
                                      size = size * size_factor))
    # cells %>% mutate(color = lighter_col(col, lighter_factor * (1 + (size - 2) * 0.3), alpha = alpha),
    #                  size = size * size_factor))
  }
  return(cells3D)
}

all_cells <- cells
if (hires) {
  cells <- add_3D_effect(all_cells, 40)
} else {
  cells <- add_3D_effect(all_cells, 5)
}

cells <- cells %>% arrange(desc(y), desc(x), desc(size))
# cells <- cells %>% filter(sqrt((orig_x - 16)^2 + (orig_y - 11)^2) + sqrt((orig_x - 5)^2 + (orig_y - 3)^2) > 24 | col_factor < 0.1)
cells <- cells %>% filter(sqrt((orig_x - 16)^2 + (orig_y - 11)^2) + sqrt((orig_x - 5)^2 + (orig_y - 3)^2) > 24 | (col %in% colors[1] & col_factor < 0.4))
ggplot() +
  geom_polygon(data = background, aes(x = x, y = y, group = group), fill = background$color) +
  geom_point(data = cells, aes(x = x, y = y),
             size = cells$size,
             col = cells$color) +
  coord_fixed()

create_flying_cell <- function(x, y, size, color, alpha = 1) {
  if (hires) {
    ball <- add_3D_effect(tibble(x = x, y = y, size = size, col = color), size * 10, alpha = alpha)
  } else {
    ball <- add_3D_effect(tibble(x = x, y = y, size = size, col = color), size * 0.5, alpha = alpha)
  }
}

create_flying_cell_with_woosh <- function(x, y, x0, y0, size, color, steps = 10) {
  coord_x <- seq(x0, x, length.out = steps)
  coord_y <- seq(y0, y, length.out = steps)
  sizes <- seq(size * 0.7, size, length.out = steps)
  woosh <- tibble()
  for (i in 1:(steps - 1)) {
    woosh <- rbind(woosh,
                   create_flying_cell(coord_x[i], coord_y[i], sizes[i], color)
    )
  }
  woosh <- woosh %>% arrange(color)
  # for (i in 1:(steps - 1)) {
  #   woosh <- rbind(woosh,
  #                  create_flying_cell(coord_x[i], coord_y[i], sizes[i], color, alpha = 0.01)
  #   )
  # }
  # woosh <- rbind(woosh,
  #                create_flying_cell(x, y, size, color)
  # )
  return(woosh)
}


if (hires) {
  background_dots <- tibble(x = -125 + runif(20000, 0, 1) * 400,
                            y = -125 + runif(20000, 0, 1) * 550,
                            size = runif(20000, 0.5, 10))
} else {
  background_dots <- tibble(x = -125 + runif(200, 0, 1) * 400,
                            y = -125 + runif(200, 0, 1) * 550,
                            size = runif(200, 0.5, 10))
}
background_color <- rgb(0.01, 0.05, 0.15)
background_contrast_color <- rgb(0.1, 0.3, 0.4)
# show_col(c(background_color, background_contrast_color))

flying_cells <- rbind(
  create_flying_cell(200, 310, 4.3, colors[6]),
  create_flying_cell(211, 284, 6, colors[6]),
  create_flying_cell(-60, 280, 7, colors[6]),
  create_flying_cell(120, 287, 8.2, colors[6]),
  
  create_flying_cell(185, 252, 8.5, colors[5]),
  create_flying_cell(85, 261, 8, colors[5]),
  create_flying_cell(110, 155, 4.7, colors[5]),
  create_flying_cell(-8, 146, 5, colors[5]),
  
  create_flying_cell(143, 306, 4.1, colors[4]),
  create_flying_cell(130, 238, 8.1, colors[4]),
  create_flying_cell(35, 220, 7, colors[4]),
  create_flying_cell(172, 219, 6.2, colors[4]),
  
  create_flying_cell(-8, 239, 7.2, colors[3]),
  create_flying_cell(40, 201, 7.2, colors[3]),
  create_flying_cell(52, 140, 7, colors[3]),
  create_flying_cell(-31, 138, 5, colors[3]),
  
  create_flying_cell(-8, 101, 5, colors[2]),
  create_flying_cell(-9, 70, 4, colors[2]),
  create_flying_cell(15, 17, 3.5, colors[2]),
  create_flying_cell(19, 48, 4.5, colors[2]),
  create_flying_cell(52, 55, 4, colors[2]),
  
  # create_flying_cell_with_woosh(-90, 259, -87, 263, 15, colors[5]),
  # create_flying_cell_with_woosh(-50, 210, -48, 211, 9, colors[4]),
  # create_flying_cell_with_woosh(150, 110, 145, 111, 20, colors[3], steps = 6),
  # create_flying_cell_with_woosh(210, 170, 206, 169, 12, colors[2], steps = 7),
  # create_flying_cell_with_woosh(90, 45, 87, 46, 11, colors[1], steps = 15)
  
  create_flying_cell(-90, 259, 12, colors[6]),
  create_flying_cell(-50, 210, 9, colors[5]),
  create_flying_cell(150, 110, 13, colors[4]),
  create_flying_cell(210, 170, 10, colors[3]),
  create_flying_cell(90, 45, 7, colors[2])
)

g <- ggplot() +
  geom_point(data = background_dots,
             aes(x = x, y = y),
             size = background_dots$size,
             col = background_contrast_color, alpha = 0.03) +
  geom_polygon(data = background, aes(x = x, y = y, group = group), fill = background$color) +
  geom_point(data = cells, aes(x = x, y = y),
             size = cells$size,
             col = cells$color) +
  geom_path(data = grey_ellipses,
            aes(x = x, y = y, group = ellipse, size = size),
            color = "grey",
            alpha = 0.3, linejoin = "mitre", show.legend = F) +
  geom_segment(data = ellipses_arcs3D,
               aes(x = x, y = y, xend = xend, yend = yend),
               size = ellipses_arcs3D %>% pull(line_width),
               col = ellipses_arcs3D %>% pull(col),
               alpha = 1, linejoin = "bevel", lineend = "round", show.legend = F) +
  geom_point(data = flying_cells, aes(x = x, y = y),
             size = flying_cells$size,
             col = flying_cells$color) +
  coord_fixed(xlim = c(-100, 250), ylim = c(-100, 400)) +
  scale_size(range = c(3, 0.5)) +
  theme_void() +
  # geom_vline(xintercept = (-10:25) * 10, col = "white", size = 0.1) +
  # geom_vline(xintercept = (-2:5) * 50, col = "white", size = 0.3) +
  # geom_hline(yintercept = (0:50) * 10, col = "white", size = 0.1) +
  # geom_hline(yintercept = (0:10) * 50, col = "white", size = 0.3) +
  theme(plot.background = element_rect(fill = background_color))
# ggsave("cover.pdf", plot = g, device = "pdf", width = 8, height = 11)
# if(hires) {
#   # ggsave("cover.tiff", plot = g, device = "tiff", width = 8, height = 11, dpi = 1200)
#   ggsave("cover.png", plot = g, device = "png", width = 8, height = 11, dpi = 1200)
# }
print(g)
