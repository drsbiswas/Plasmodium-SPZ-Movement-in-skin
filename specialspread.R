library(readxl)
library(dplyr)

# =========================================================
# 1. TRAJECTORY PLOT FUNCTION
# =========================================================

plot_trajectory_base_viewer <- function(data,
                                        track_col = "Clean TrackIDs",
                                        x_col = "Position X (um)",
                                        y_col = "Position Y (um)",
                                        track_id = NULL,
                                        legend_pos = "topright") {
  
  # Pick the track
  ids <- data[[track_col]] %>% na.omit() %>% unique()
  
  if (length(ids) == 0) {
    stop("No valid track IDs found.")
  }
  
  if (is.null(track_id)) {
    track_id <- sample(ids, 1)
  } else {
    if (!track_id %in% ids) {
      stop(paste0("Track ID not found: ", track_id))
    }
  }
  
  message("Using track: ", track_id)
  
  df <- data %>%
    dplyr::filter(.data[[track_col]] == track_id)
  
  if (nrow(df) < 2) {
    stop("Selected track has fewer than 2 points.")
  }
  
  x <- df[[x_col]]
  y <- df[[y_col]]
  
  # Bounds and f
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  
  dx <- xmax - xmin
  dy <- ymax - ymin
  f_val <- sqrt(dx^2 + dy^2)
  
  # Rounded values for display
  dx_r <- round(dx, 2)
  dy_r <- round(dy, 2)
  f_r  <- round(f_val, 2)
  
  f_expr <- substitute(
    sqrt(dx^2 + dy^2) * " =" * f * " µm",
    list(dx = dx_r, dy = dy_r, f = f_r)
  )
  
  # Start/end
  x_start <- x[1]
  y_start <- y[1]
  
  x_end <- x[length(x)]
  y_end <- y[length(y)]
  
  # Plot
  plot(
    x, y,
    type = "l",
    lwd = 2,
    col = "black",
    main = "",
    ylim= c(156,170),
    xlab = "X position (µm)",
    ylab = "Y position (µm)",
    cex.lab = 2.5,
    cex.axis = 2.0
  )
  
  points(x, y, pch = 16, col = "black", cex = 1.0)
  
  # Diagonal spread line
  segments(xmin, ymin, xmax, ymax, col = "red", lwd = 2.5)
  
  # Min/max reference lines
  abline(v = c(xmin, xmax), col = "darkgreen", lwd = 2, lty = 3)
  abline(h = c(ymin, ymax), col = "blue", lwd = 2, lty = 3)
  
  # Start/end markers
  points(x_start, y_start, pch = 21, bg = "white", cex = 2.2)
  points(x_end, y_end, pch = 23, bg = "white", cex = 2.2)
  
  text(x_start, y_start, "Start", pos = 3, cex = 2.0)
  text(x_end, y_end, "End", pos = 3, cex = 2.0)
  
  grid()
  
  # Solid non-transparent legend
  legend(
    legend_pos,
    legend = list(
      paste("Track ID:", track_id),
      "Trajectory",
      "X min/max (vertical)",
      "Y min/max (horizontal)",
      "Overall spatial spread (f) = ",
      f_expr
    ),
    lty = c(NA, 1, 3, 3, 1, NA),
    lwd = c(NA, 2, 2.5, 2.5, 2.5, NA),
    col = c(NA, "black", "darkgreen", "blue", "red", NA),
    pch = c(NA, NA, NA, NA, NA, NA),
    bty = "o",
    bg = "white",
    box.col = "black",
    cex = 1.8
  )
}


# =========================================================
# 2. READ DATA
# =========================================================

dat <- read_excel(
  "~/Desktop/Movement Project/movHMM/Sporozoites_Master_File_Hopp_e21_ver3_EH.xlsx"
)


# =========================================================
# 3. COMPUTE TRACK SPREAD f
# =========================================================

track_spread <- dat %>%
  group_by(Species, `Time Post-Inoculation (min)`, `Clean TrackIDs`) %>%
  summarise(
    f = sqrt(
      (
        max(`Position X (um)`, na.rm = TRUE) -
          min(`Position X (um)`, na.rm = TRUE)
      )^2 +
        (
          max(`Position Y (um)`, na.rm = TRUE) -
            min(`Position Y (um)`, na.rm = TRUE)
        )^2
    ),
    .groups = "drop"
  ) %>%
  mutate(Species_upper = toupper(Species))


# =========================================================
# 4. SPLIT BY SPECIES
# =========================================================

pf_data <- track_spread %>%
  filter(Species_upper == "PF")

py_data <- track_spread %>%
  filter(Species_upper == "PY")

time_pf <- sort(unique(pf_data$`Time Post-Inoculation (min)`))
time_py <- sort(unique(py_data$`Time Post-Inoculation (min)`))


# =========================================================
# 5. COLORS AND LINE TYPES
# =========================================================

cols_pf <- c("black", "blue", "red", "darkgreen", "brown", "purple")[seq_along(time_pf)]
cols_py <- c("black", "blue", "red", "darkgreen", "brown", "purple")[seq_along(time_py)]

lty_pf <- rep(1:6, length.out = length(time_pf))
lty_py <- rep(1:6, length.out = length(time_py))


# =========================================================
# 6. COMMON X RANGE AND HISTOGRAM BREAKS
# =========================================================

x_min <- 10^-2
x_max <- 10^3

common_breaks <- 10^seq(-2, 3, length.out = 41)


# =========================================================
# 7. GLOBAL Y RANGE FOR PANELS B AND C
# =========================================================

max_count_global <- 0

for (sp_data in list(pf_data, py_data)) {
  
  for (tp in unique(sp_data$`Time Post-Inoculation (min)`)) {
    
    fvals <- sp_data$f[
      sp_data$`Time Post-Inoculation (min)` == tp &
        sp_data$f > 0
    ]
    
    if (length(fvals) > 0) {
      h <- hist(fvals, breaks = common_breaks, plot = FALSE)
      max_count_global <- max(max_count_global, h$counts, na.rm = TRUE)
    }
  }
}


# =========================================================
# 8. LOG-X HISTOGRAM FUNCTION
# =========================================================

plot_logx_hist_math <- function(data, time_points, cols, ltys, species_label,
                                x_min, x_max, y_max, breaks) {
  
  plot(
    NA, NA,
    xlim = c(x_min, x_max),
    ylim = c(0, y_max * 1.12),
    log = "x",
    axes = FALSE,
    xlab = expression(paste("Track spread ", italic(f), " (", mu, "m)")),
    ylab = expression("Number of trajectories"),
    cex.lab = 2.5,
    cex.axis = 2.0
  )
  
  # Custom X-axis ticks
  xticks <- 10^seq(-2, 3, by = 1)
  
  axis(
    1,
    at = xticks,
    labels = parse(text = paste0("10^", seq(-2, 3, by = 1))),
    cex.axis = 2.0,
    lwd = 1.5,
    lwd.ticks = 2.5
  )
  
  axis(
    2,
    cex.axis = 2.0,
    lwd = 1.5,
    lwd.ticks = 2.5
  )
  
  box(lwd = 2.0)
  
  # Draw histogram curves
  for (i in seq_along(time_points)) {
    
    fvals <- data$f[
      data$`Time Post-Inoculation (min)` == time_points[i] &
        data$f > 0
    ]
    
    if (length(fvals) > 0) {
      h <- hist(fvals, breaks = breaks, plot = FALSE)
      lines(
        h$mids,
        h$counts,
        col = cols[i],
        lty = ltys[i],
        lwd = 3.0
      )
    }
  }
  
  # Vertical dashed threshold line at f = 5 µm
  abline(v = 5, col = "black", lty = 2, lwd = 2.5)
  
  # Label for threshold line
  text(
    x = 5,
    y = y_max * 1.04,
    labels = expression(italic(f) == 5),
    pos = 2,
    cex = 2
  )
  
  # Species label inside plot
  text(
    x = 0.02,
    y = y_max * 1.04,
    labels = bquote(italic(.(species_label))),
    cex =2.5,
    adj = c(0, 1)
  )
  
  # Legend
  legend(
    "topright",
    legend = paste(time_points, "min"),
    col = cols,
    lty = ltys,
    lwd = 2.5,
    bty = "o",
    cex = 2.0
  )
}


# =========================================================
# 9. PRINT MAX-f SUMMARY
# =========================================================

max_f_values <- track_spread %>%
  group_by(Species_upper) %>%
  summarise(
    Max_f = max(f, na.rm = TRUE),
    Track_with_max_f = `Clean TrackIDs`[which.max(f)],
    Time_point = `Time Post-Inoculation (min)`[which.max(f)],
    .groups = "drop"
  )

print(max_f_values)


# =========================================================
# 10. COMBINED FIGURE: PANELS A, B, C IN ONE ROW
# =========================================================

pdf(
  file = "~/Desktop/Figure_ABC_single_row.pdf",
  width = 16.5,
  height = 6.8,
  family = "Times",
  useDingbats = FALSE
)

# One row, three columns
# oma gives extra outer margin so panel letters are not cut
par(
  mfrow = c(1, 3),
  oma = c(0, 0, 3.5, 0)
)


# -------------------------
# Panel A: trajectory example
# -------------------------

par(mar = c(6, 6, 5.5, 2) + 0.1)

plot_trajectory_base_viewer(
  dat,
  track_id = "PF.160323.30min.C0.1001036601.3",
  legend_pos = "topleft"
)

mtext(
  "A",
  side = 3,
  adj = -0.14,
  line = 2.2,
  font = 2,
  cex = 2.4
)


# -------------------------
# Panel B: P. falciparum
# -------------------------

par(mar = c(6, 6, 5.5, 2) + 0.1)

plot_logx_hist_math(
  pf_data,
  time_pf,
  cols_pf,
  lty_pf,
  "P. falciparum",
  x_min,
  x_max,
  max_count_global,
  common_breaks
)

mtext(
  "B",
  side = 3,
  adj = -0.14,
  line = 2.2,
  font = 2,
  cex = 2.4
)


# -------------------------
# Panel C: P. yoelii
# -------------------------

par(mar = c(6, 6, 5.5, 2) + 0.1)

plot_logx_hist_math(
  py_data,
  time_py,
  cols_py,
  lty_py,
  "P. yoelii",
  x_min,
  x_max,
  max_count_global,
  common_breaks
)

mtext(
  "C",
  side = 3,
  adj = -0.14,
  line = 2.2,
  font = 2,
  cex = 2.4
)

dev.off()

