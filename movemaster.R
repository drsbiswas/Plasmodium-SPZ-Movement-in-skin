
## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  message = FALSE, error = FALSE, warning = FALSE,
  comment = NA
)

## ----load-package, echo = FALSE-----------------------------------------------
library(moveHMM)
library(dplyr)
library(ggplot2)
library(readxl)
set.seed(342)
theme_set(theme_minimal())
getwd()
setwd("~/Desktop")
## ----show-data----------------------------------------------------------------
filtered_data_dplyr<- read_excel("~/Downloads/Sporozoites_Master_File_Hopp_e21_ver3_EH.xlsx")

# Convert time to numeric once
filtered_data_dplyr$time_min_num <- suppressWarnings(
  as.numeric(filtered_data_dplyr$`Time Post-Inoculation (min)`)
)

library(dplyr)
# timewise pf/py -----------------------------------------
raw <- filtered_data_dplyr %>%
  filter(
    tolower(Species) == "py",
    !is.na(time_min_num),
    time_min_num == 60
  ) %>%
  select(
    `Clean TrackIDs`,
    `Time (s)`,
    `Position X (um)`,
    `Position Y (um)`
  )


#------ whole pf/py--------------------------------------
#raw <- filtered_data_dplyr %>%
#  filter(
#    tolower(Species) == "py",
#    !is.na(time_min_num)   # keep all non-NA times
#  ) %>%
#  select(
#    `Clean TrackIDs`,
#    `Time (s)`,
#    `Position X (um)`,
#    `Position Y (um)`
#  )
#
#raw<- filtered_data_dplyr[ keep_rows , 
                                 #c("Clean TrackIDs","Time (s)",
                                  #"Position X (um)","Position Y (um)") ]

# Rename columns to the clean set you wanted
names(raw) <- c("ID","time","x","y")

# optional: save
# write.csv(filtered, "filtered_pf_time5.csv", row.names = FALSE)
  # Expecting columns: X, Y, Time, TrackID
  #rename(x = X, y = Y, time = Time, ID = TrackID)


# Basic checks
stopifnot(all(c("x","y","time","ID") %in% names(raw)))
raw <- raw %>% arrange(ID, time)
# If your time is in seconds (typical), make minutes & centered/scaled versions

total_tracks <- n_distinct(raw$ID, na.rm = FALSE)
total_tracks




raw <- raw %>%
  group_by(ID) %>%
  mutate(
    time_sec = time - min(time, na.rm = TRUE),
    time_min = time_sec/60
  ) %>%
  ungroup()


#total_tracks <- n_distinct(raw$ID, na.rm = TRUE)
#total_tracks


library(dplyr)
#extra<- raw %>%
#  group_by(ID) %>%
#  summarise(n_points = n()) %>%
#  filter(n_points < 3)

# assuming your dataset is in `raw` and has a column ID for track IDs
raw <- raw %>%
  group_by(ID) %>%        # group by track
  filter(n() >2) %>%    # keep only tracks with > 2 rows
  ungroup()

#total_tracks <- n_distinct(raw$ID, na.rm = TRUE)
#total_tracks

raw$x <- as.numeric(raw$x)
raw$y <- as.numeric(raw$y)
raw$time <- as.numeric(raw$time)
raw$ID <- as.factor(raw$ID)

library(dplyr)

#data <- data %>%
#  arrange(ID, time) %>%                 # ensure stable order
#  distinct(ID, time, .keep_all = TRUE)  # drops later duplicates
#total_tracks <- n_distinct(data$ID, na.rm = TRUE)
#total_tracks

#data <- raw%>% filter(time != 0)
#data<- na.omit()
#total_tracks <- n_distinct(data$ID, na.rm = TRUE)
#total_tracks

# If coordinates are in microns/pixels, that's fine—treat as planar Cartesian.
# moveHMM's "UTM" works for any planar coords (units are whatever you provide).
data <- prepData(raw, type = "UTM")
head(data)
#data<- na.omit(data)
total_tracks <- n_distinct(data$ID, na.rm = TRUE)
total_tracks


################################################################################
## ZERO STEP REMOVAL############################################################
################################################################################

# Count zeros among finite steps
n_zero   <- sum(is.finite(data$step) & data$step == 0, na.rm = TRUE)
n_finite <- sum(is.finite(data$step))
pct_zero <- if (n_finite > 0) 100 * n_zero / n_finite else 0

message(sprintf(
  "Zero-step removal: %d of %d finite steps are zero (%.2f%%). Removing zeros.",
  n_zero, n_finite, pct_zero
))

# Remove only exact zeros; keep NAs (first obs per ID often has NA step/angle)
data <- dplyr::filter(data, !(is.finite(step) & step == 0))




## DOUBLE CHECK ################################################################

# Detect zero steps among finite values
is_zero_step <- is.finite(data$step) & data$step == 0
n_zero       <- sum(is_zero_step, na.rm = TRUE)
n_finite     <- sum(is.finite(data$step))
pct_zero     <- if (n_finite > 0) 100 * n_zero / n_finite else 0

message(sprintf(
  "Zero-step check: %d of %d finite steps are zero (%.2f%%). Removing zero steps.",
  n_zero, n_finite, pct_zero
))

# Remove only zero steps; keep NA rows (first obs in each track often has NA step/angle)
data <- dplyr::filter(data, !(is.finite(step) & step == 0))

# (Optional) quick sanity check
stopifnot(sum(is.finite(data$step) & data$step == 0) == 0)

################################################################################
## ELEMINATION OF BUZZING/NOT ACTUALY MOVING SPZ################################
################################################################################

threshold_um <- 5  # keep tracks whose spread is at least 3 µm

track_spread <- data %>%
  dplyr::mutate(ID_chr = as.character(ID)) %>%
  dplyr::group_by(ID_chr) %>%
  dplyr::summarise(
    n_points  = dplyr::n(),
    dx_um     = max(x, na.rm = TRUE) - min(x, na.rm = TRUE),
    dy_um     = max(y, na.rm = TRUE) - min(y, na.rm = TRUE),
    spread_um = sqrt(dx_um^2 + dy_um^2),
    .groups = "drop"
  )

keep_ids <- track_spread %>%
  dplyr::filter(is.finite(spread_um) & spread_um >= threshold_um) %>%
  dplyr::pull(ID_chr)

drop_ids <- setdiff(unique(as.character(data$ID)), keep_ids)

message(sprintf(
  "Track spread filter (>= %.2f µm): kept %d / %d IDs (removed %d = %.1f%%).",
  threshold_um, length(keep_ids), length(keep_ids) + length(drop_ids),
  length(drop_ids),
  100 * length(drop_ids) / (length(keep_ids) + length(drop_ids))
))
if (length(drop_ids) > 0) {
  message(sprintf("Removed IDs (examples): %s", paste(head(drop_ids, 10), collapse = ", ")))
}

data <- data %>%
  dplyr::mutate(ID_chr = as.character(ID)) %>%
  dplyr::filter(ID_chr %in% keep_ids) %>%
  dplyr::select(-ID_chr)

# Safety check
stopifnot(length(unique(as.character(data$ID))) == length(keep_ids))

################################################################################
## 2 STATE FITTING #############################################################
################################################################################

# Step means at 25th/75th percentiles; SDs = half the means (stable starts)
# === helpers ================================================================
pos_steps <- data$step[is.finite(data$step) & data$step > 0]
stopifnot(length(pos_steps) >= 20)

# rough circular concentration from pooled angles (for scaling kappas if desired)
ang <- data$angle[is.finite(data$angle)]
Rbar <- sqrt(mean(cos(ang))^2 + mean(sin(ang))^2)
kappa_pooled <-
  if (Rbar < 0.53) { 2*Rbar + Rbar^3 + 5*Rbar^5/6 } else
    if (Rbar < 0.85) { -0.4 + 1.39*Rbar + 0.43/(1 - Rbar) } else
    { 1/(Rbar^3 - 4*Rbar^2 + 3*Rbar) }

# Keep kappas biologically plausible regardless of Rbar
scale_k <- function(x, lo = 0.2, hi = 6) pmin(pmax(x, lo), hi)




# Step means at 25th/75th percentiles; SDs = half the means (stable starts)
q25 <- as.numeric(quantile(pos_steps, 0.25))
q75 <- as.numeric(quantile(pos_steps, 0.75))
m1  <- q25; m2 <- q75
s1  <- max(m1/2, .Machine$double.eps)
s2  <- max(m2/2, .Machine$double.eps)

stepPar0_2 <- c(m1, m2, s1, s2)

# Angles: means ~ 0 (forward bias), kappa_run > kappa_search
k1 <- scale_k(0.6 * kappa_pooled, lo = 0.3, hi = 2.0)  # search (diffuse)
k2 <- scale_k(1.8 * kappa_pooled, lo = 1.2, hi = 5.0)  # run (persistent)
anglePar0_2 <- c(3.14, 0, k1, k2)

cat("2-state stepPar0:", paste(signif(stepPar0_2,3), collapse=", "), "\n")
cat("2-state anglePar0:", paste(signif(anglePar0_2,3), collapse=", "), "\n")

mod2 <- fitHMM(data=data, nbStates=2, stepPar0=stepPar0_2, anglePar0=anglePar0_2)
print(mod2)
AIC(mod2)
plot(mod2, ask = FALSE, animals = 12)  # overview plots
#plot(mod2, plotCI=T)

################################################################################
## 3 STATE FITTING #############################################################
################################################################################

q <- as.numeric(quantile(pos_steps, c(0.20, 0.50, 0.80)))
m1 <- q[1]; m2 <- q[2]; m3 <- q[3]
s1 <- max(m1/2, .Machine$double.eps)
s2 <- max(m2/2, .Machine$double.eps)
s3 <- max(m3/2, .Machine$double.eps)

#stepPar0_3 <- c(m1, m2, m3, s1, s2, s3)
stepPar0_3 <- c(0.22, 0.6, 1.8, 0.15, 0.39, 1.61)

# Angles: all forward-biased; κ grows from pause/search → run
k1 <- scale_k(0.3 * kappa_pooled, lo = 0.2, hi = 1.0)  # pause/jitter or tight search
k2 <- scale_k(1.0 * kappa_pooled, lo = 0.8, hi = 2.0)  # meander
k3 <- scale_k(2.5 * kappa_pooled, lo = 2.0, hi = 6.0)  # persistent glide
anglePar0_3 <- c(3.14, 0, 0, k1, k2, k3)

cat("3-state stepPar0:", paste(signif(stepPar0_3,3), collapse=", "), "\n")
cat("3-state anglePar0:", paste(signif(anglePar0_3,3), collapse=", "), "\n")

mod3 <- fitHMM(data=data, nbStates=3, stepPar0=stepPar0_3, anglePar0=anglePar0_3)
print(mod3)
AIC(mod3)
AIC(mod2, mod3)
plot(mod3, ask = FALSE, animals =6)  # overview plots
#plot(mod3, plotCI=T)
plotPR(mod3)
#----------------- for py ---------------------------------------------------
# 1) Reset any weird settings from earlier base plots
#par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0,0,0,0))

# 2) If in RStudio, click "Zoom" on the Plots pane, or open a larger device:
#if (.Platform$OS.type == "windows") windows(width = 9, height = 7)
#if (Sys.info()[["sysname"]] == "Darwin") quartz(width = 9, height = 7)  # macOS
#if (Sys.info()[["sysname"]] == "Linux") x11(width = 9, height = 7)

# 3) Now plot
#plotPR(mod3)
#------------------------------------------------------------------------------
# Add most likely state sequence to data
data$state <- factor(viterbi(mod3))
# Plot tracks coloured by state
library(ggplot2)
ggplot(data, aes(x, y, col = state, group = ID)) +
  geom_path() +
  coord_equal()
# Plot step lengths coloured by state
ggplot(data, aes(x = 1:nrow(data), y = step, col = state, group = ID)) +
  geom_point(size = 0.3)
################################################################################
## TRAJECTORIES PLOT ###########################################################
################################################################################
##### multiple##################################################################
library(dplyr)
library(ggplot2)
library(moveHMM)

# --- Build plotting data from your fit ----------------------------------------
real_plot <- data
real_plot$state <- viterbi(mod3)  # numeric 1..K
real_plot <- real_plot %>% mutate(x = as.numeric(x), y = as.numeric(y))

# Custom legend labels for states
state_to_label <- function(s) {
  sn <- suppressWarnings(as.integer(as.character(s)))
  if (any(is.na(sn))) sn <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(s))))
  factor(
    sn,
    levels = c(1, 2, 3),
    labels = c("1. Reverse and Searching",
               "2. Reverse and Meandering",
               "3. Long Travel")
  )
}
real_plot <- real_plot %>% mutate(state_lbl = state_to_label(state))

# --- Choose 4 IDs and lock facet order ---------------------------------------
ids4 <- unique(real_plot$ID)[1:4]             # or sample(unique(real_plot$ID), 4)
df4  <- real_plot %>% filter(ID %in% ids4)

df4_clean <- df4 %>%
  filter(is.finite(x), is.finite(y)) %>%
  group_by(ID) %>%
  mutate(.ord = row_number(), n_pts = n()) %>%
  arrange(ID, .ord) %>%
  ungroup() %>%
  filter(n_pts >= 2) %>%
  mutate(ID = factor(ID, levels = ids4))      # keep facet order = ids4

# Start/end markers AFTER cleaning
starts <- df4_clean %>% group_by(ID) %>% slice_head(n = 1) %>% ungroup() %>% mutate(point = "Start")
ends   <- df4_clean %>% group_by(ID) %>% slice_tail(n = 1) %>% ungroup() %>% mutate(point = "End")
pts    <- bind_rows(starts, ends)

# --- Panel tags A–D (one per facet) ------------------------------------------
panel_tags <- data.frame(
  ID    = factor(ids4, levels = ids4),
  label = c("A", "B", "C", "D"),
  x     = -Inf,
  y     =  Inf
)

# --- Plot ---------------------------------------------------------------------
ggplot(df4_clean, aes(x, y, group = ID)) +
  geom_path(aes(colour = state_lbl), linewidth = 0.6) +
  geom_point(data = pts, aes(x, y, shape = point),
             size = 2.2, colour = "black", fill = "white", stroke = 0.8) +
  # Panel labels in the top-left corner of each facet
  geom_text(data = panel_tags, aes(x, y, label = label),
            inherit.aes = FALSE, hjust = -0.2, vjust = 1.2,
            size = 5, fontface = "bold") +
  facet_wrap(~ ID, ncol = 2, scales = "free") +
  scale_color_discrete(name = "State") +
  scale_shape_manual(values = c(Start = 21, End = 23), name = "Point") +
  labs(x = "Axial Travel (µm)", y = "Transverse Travel (µm)") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_blank(),     # remove facet strip titles
    strip.background = element_blank()
  )
######single####################################################################
library(dplyr)
library(ggplot2)
library(moveHMM)

# --- Prep real data + states ---------------------------------------------------
real_plot <- data
real_plot$state <- viterbi(mod3)
real_plot <- real_plot %>%
  mutate(x = as.numeric(x), y = as.numeric(y)) %>%
  group_by(ID) %>% mutate(idx = row_number()) %>% ungroup()

# State labels as requested
state_to_label <- function(s) {
  sn <- suppressWarnings(as.integer(as.character(s)))
  if (any(is.na(sn))) sn <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(s))))
  factor(sn, levels = c(1,2,3),
         labels = c("State 1", "State 2", "State 3"))
}
real_plot <- real_plot %>% mutate(state_lbl = state_to_label(state))

# --- Pick one ID ---------------------------------------------------------------
id_to_plot <- as.character(unique(real_plot$ID)[1])   # change as needed

df1 <- real_plot %>%
  filter(ID == id_to_plot, is.finite(x), is.finite(y)) %>%
  arrange(idx)

# Per-step segments to ensure continuity
segdf <- df1 %>%
  transmute(x, y, xend = lead(x), yend = lead(y), state_lbl) %>%
  filter(is.finite(x), is.finite(y), is.finite(xend), is.finite(yend))

# Start/End markers
starts <- df1 %>% slice_head(n = 1) %>% mutate(point = "Start")
ends   <- df1 %>% slice_tail(n = 1) %>% mutate(point = "End")
pts    <- bind_rows(starts, ends)

# In-panel ID label position
xr <- diff(range(df1$x, na.rm = TRUE)); yr <- diff(range(df1$y, na.rm = TRUE))
lab_x <- min(df1$x, na.rm = TRUE) + 0.02 * xr
lab_y <- max(df1$y, na.rm = TRUE) - 0.02 * yr

library(ggnewscale)   # install.packages("ggnewscale") if needed

# --- Per-vertex step lengths (for dot colours) ---
points_step <- df1 %>%
  mutate(step_len = if ("step" %in% names(.)) step else
    sqrt((x - dplyr::lag(x))^2 + (y - dplyr::lag(y))^2)) %>%
  filter(is.finite(step_len)) %>%
  transmute(x = x, y = y, step_len)

# Labels/limits for step-length legend
min_step <- min(points_step$step_len, na.rm = TRUE)
max_step <- max(points_step$step_len, na.rm = TRUE)
fmt_um   <- function(v) paste0(formatC(v, format = "f", digits = 2), " µm")

# Dummy data to inject a legend note (drawn invisible but shown in legend)
note_df <- data.frame(
  x = min(df1$x, na.rm = TRUE),
  y = min(df1$y, na.rm = TRUE),
  xend = min(df1$x, na.rm = TRUE),
  yend = min(df1$y, na.rm = TRUE),
  note = "Time after inoculation = 5 min"
)

p_single <- ggplot() +
  # 1) Continuous path via per-step segments (colour = STATE)
  geom_segment(
    data = segdf,
    aes(x = x, y = y, xend = xend, yend = yend, colour = state_lbl),
    linewidth = 1.0, lineend = "round"
  ) +
  scale_color_manual(
    name   = "State",
    values = c("State 1" = "#0072B2", "State 2" = "#E69F00", "State 3" = "#009E73"),
    guide  = guide_legend(order = 1)
  ) +
  
  # 2) New colour scale for step-length dots
  ggnewscale::new_scale_color() +
  geom_point(
    data = points_step,
    aes(x = x, y = y, colour = step_len),
    size = 1.2, alpha = 0.9
  ) +
  scale_color_viridis_c(
    name   = "Step length",
    limits = c(min_step, max_step),
    breaks = c(min_step, max_step),
    labels = c(fmt_um(min_step), fmt_um(max_step)),
    guide  = guide_colourbar(order = 2)
  ) +
  
  # 3) Start/End markers
  geom_point(
    data = pts, aes(x, y, shape = point),
    size = 2.6, colour = "black", fill = "white", stroke = 0.9
  ) +
  scale_shape_manual(
    name  = "Point",
    values = c(Start = 21, End = 23),
    guide  = guide_legend(order = 3)
  ) +
  
  # 4) Legend note: "Time after inoculation = 5 min"
  ggnewscale::new_scale("linetype") +
  geom_segment(
    data = note_df,
    aes(x = x, y = y, xend = xend, yend = yend, linetype = note),
    alpha = 0
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = c("Time after inoculation = 5 min" = "solid"),
    guide  = guide_legend(order = 4, override.aes = list(alpha = 1, colour = "black"))
  ) +
  
  # 5) ID label inside panel
  annotate("text", x = lab_x, y = lab_y,
           label = paste0("ID: ", id_to_plot),
           hjust = 0, vjust = 1.1, size = 5, fontface = "bold") +
  
  labs(x = "X (µm)", y = "Y (µm)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.border        = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor    = element_blank(),
    legend.position     = c(0.98, 0.02),            # bottom-right inside
    legend.justification= c(1, 0),
    legend.background   = element_blank(),          # no legend box
    legend.box.background = element_blank(),
    legend.key          = element_blank()
  )

print(p_single)

## =========================
## BASE R: Single trajectory with state-colored step dots
## =========================

## ====== BASE R: single trajectory with grid + true top-left legend ======

## choose the ID you want (you already have data, mod3)
id_to_plot <- as.character(unique(data$ID)[6])  # change as needed

## viterbi states & subset one ID
state_num <- as.integer(viterbi(mod3))
dat_id <- transform(data, state = state_num)
dat_id <- subset(dat_id, ID == id_to_plot)

## clean & order
x <- as.numeric(dat_id$x); y <- as.numeric(dat_id$y)
ok <- is.finite(x) & is.finite(y)
x <- x[ok]; y <- y[ok]
state <- dat_id$state[ok]
n <- length(x); stopifnot(n >= 2)

## state palette + labels
pal <- c("State 1" = "#0072B2", "State 2" = "#E69F00", "State 3" = "#009E73")
state_lbl <- factor(state, levels = c(1,2,3), labels = names(pal))

## per-step colors (state at the start of each step)
seg_col <- unname(pal[state_lbl[1:(n-1)]])

## set up plot
xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
dx <- diff(xr); dy <- diff(yr)

## set up plot window & custom ticks  ------------------------------------------
ylim <- c(106,114)
xlim <- c(55,61)
xlim <- xr
ylim <- yr
## choose tick positions (edit 'by' to taste)
xticks <- seq(xlim[1], xlim[2], by = 2)   # e.g., 80,82,...,92
yticks <- seq(ylim[1], ylim[2], by = 2)

getwd()
setwd("/Users/sayantanbiswas/Desktop/Movement Project/final plot")
pdf("E.pdf", width = 5, height = 5, onefile = FALSE)

par(mar = c(4.2, 4.2, 1.2, 1.2), mgp = c(2.2, 0.6, 0), xaxs = "i", yaxs = "i")
plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "X (µm)", ylab = "Y (µm)",
     asp = 1, xaxt = "n", yaxt = "n")

## draw axes with your labels
axis(1, at = xticks, labels = sprintf("%.0f", xticks), las = 1,cex.axis = 0.9)
axis(2, at = yticks, labels = sprintf("%.0f", yticks), las = 1, cex.axis = 0.9)
## optional minor ticks (no labels)
axis(1, at = seq(xlim[1], xlim[2], by = 0.5), labels = FALSE, tcl = -0.2)
axis(2, at = seq(ylim[1], ylim[2], by = 0.5), labels = FALSE, tcl = -0.2)

## background grid aligned to major ticks (draw before data)
abline(v = xticks, col = "darkgray",  lty = "dotted", lwd = 1)
abline(h = yticks, col = "darkgrey",  lty = "dotted", lwd = 1)

## continuous path (colored by state)
segments(x[1:(n-1)], y[1:(n-1)], x[2:n], y[2:n], col = seg_col, lwd = 2)

## small dots at end of each step (same state color; no legend)
points(x[2:n], y[2:n], pch = 16, col = seg_col, cex = 0.6)

## start/end markers
points(x[1], y[1], pch = 21, bg = "white", col = "black", cex = 1.2, lwd = 1.2)  # Start

points(x[n], y[n], pch = 24, bg = "white", col = "black", cex = 1.2, lwd = 1.2)  # End

## separate legends
legend("top", inset = 0.02,
       legend = paste0("Track ID: ", id_to_plot),
       bty = "n", text.font = 2, cex = 0.8)

legend("bottomright", inset = 0.00,
       legend = "Time post inoculation = 30 mins",
       bty = "n", text.font = 2, cex = 0.7)
legend("bottomleft", inset = 0.00,
       legend = "P. yoelii",
       bty = "n", text.font = 2, cex = 1)
legend("bottomright", bty = "n",
       legend = c(names(pal), "Start", "End"),
       col    = c(unname(pal), "black", "black"),
       lty    = c(1,1,1, NA, NA),
       lwd    = c(2,2,2, NA, NA),
       pch    = c(NA,NA,NA, 21, 24),
       pt.bg  = c(NA,NA,NA, "white", "white"),
       pt.cex = c(1,1,1, 1.2, 1.2),
       seg.len = 2)

box()  # neat border
dev.off()
################################################################################
## TRAVEL SHARE CALCULATION ####################################################
################################################################################

# pick the fitted model you want to report
mod <- if (exists("mod3")) mod3 else if (exists("mod2")) mod2 else stop("Fit a model first (mod2/mod3).")
mod<- mod3
# Robust way to get K and state labels
sp <- stateProbs(mod)                  # n x K matrix of posteriors
K  <- as.integer(ncol(sp))             # ensure length-1 integer
stopifnot(is.finite(K), K >= 1)
states <- paste0("S", seq_len(K))
df <- mod$data
df$state_hat <- factor(viterbi(mod), levels = seq_len(K), labels = states)

# Viterbi states
df$state_hat <- factor(viterbi(mod), levels = seq_len(K), labels = states)

# Use only finite, positive steps as distance (µm)
mask    <- is.finite(df$step) & df$step > 0
df_step <- df[mask, , drop = FALSE]

# Basic counts
n_ids   <- length(unique(df_step$ID))
n_steps <- nrow(df_step)
tot_um  <- sum(df_step$step)

# ---------- DATASET-LEVEL (Viterbi) ----------
library(dplyr)
dataset_tbl <- df_step |>
  group_by(state_hat) |>
  summarise(distance_um = sum(step), .groups="drop") |>
  right_join(data.frame(state_hat = states), by="state_hat") |>
  mutate(distance_um = ifelse(is.na(distance_um), 0, distance_um),
         percentage  = 100 * distance_um / sum(distance_um))

# ---------- PER-ID (equal-weight) ----------
per_id_tbl <- df_step |>
  group_by(ID, state_hat) |>
  summarise(distance_um = sum(step), .groups="drop") |>
  right_join(expand.grid(ID = unique(df_step$ID), state_hat = states), by=c("ID","state_hat")) |>
  mutate(distance_um = ifelse(is.na(distance_um), 0, distance_um)) |>
  group_by(ID) |>
  mutate(pct = 100 * distance_um / sum(distance_um)) |>
  ungroup()

eq_weight_tbl <- per_id_tbl |>
  group_by(state_hat) |>
  summarise(pct_mean = mean(pct, na.rm=TRUE),
            pct_sd   = sd(pct,   na.rm=TRUE),
            .groups="drop")

# ---------- Main-text strings ----------
fmt1 <- function(x) sprintf("%.1f", x)
dataset_str <- paste(sprintf("%s %s%%", dataset_tbl$state_hat, fmt1(dataset_tbl$percentage)), collapse = "; ")
eq_str      <- paste(sprintf("%s %s ± %s%%", eq_weight_tbl$state_hat,
                             fmt1(eq_weight_tbl$pct_mean), fmt1(eq_weight_tbl$pct_sd)),
                     collapse = "; ")

cat(sprintf("\nTOTAL traveled path: %.1f µm across %d IDs (%d steps).\n", tot_um, n_ids, n_steps))
cat("Dataset-level Viterbi % by state: ", dataset_str, ".\n", sep = "")
cat(sprintf("Per-ID equal-weight mean ± SD: %s (n = %d IDs).\n", eq_str, n_ids))


################################################################################
# MODEL CHECKS #################################################################
#### REVARSAL MOVEMENT #########################################################
################################################################################

# pick your fitted model
mod <- if (exists("mod3")) mod3 else mod2

# 1) Which state is closest to μ ≈ π ?
wrap <- function(a) ((a + pi) %% (2*pi)) - pi
K <- ncol(stateProbs(mod))
ap <- mod$mle$anglePar
mu    <- wrap(ap[seq_len(K)])                 # first K entries are means
kappa <- ap[K + seq_len(K)]                   # next K are concentrations
rev_state <- which.min(cos(mu))               # cos(π) = -1 → smallest = most “reversal”
mu; kappa; rev_state

# 2) Dataset-level evidence of back-tracking irrespective of states
df <- mod$data
mask <- is.finite(df$angle) & is.finite(df$step) & df$step > 0
pct_back_steps  <- mean(cos(df$angle[mask]) < 0) * 100     # % of steps with >90° turn
pct_back_dist   <- 100 * sum(df$step[mask & cos(df$angle) < 0]) / sum(df$step[mask])
c(pct_back_steps = pct_back_steps, pct_back_distance = pct_back_dist)

# 3) If a reversal state exists, what share of travel is in it? (Viterbi)
df$state_hat <- viterbi(mod)
rev_share <- 100 * sum(df$step[mask & df$state_hat == rev_state]) / sum(df$step[mask])
rev_share


# ---------- (Optional) compact tables ----------
# knitr::kable(dataset_tbl, digits = 2, caption = "Dataset-level Viterbi distance share by state")
# knitr::kable(eq_weight_tbl, digits = 2, caption = "Per-ID equal-weight % (mean ± SD) by state")


################################################################################
# --- A) SIMULATE FROM FITTED HMM ---------------------------------------------
################################################################################
# =========================
# OVERLAY: REAL vs SIMULATED (MATCH REAL LENGTHS + ANCHORED STARTS)
# =========================
# Needs: your fitted model `mod3` and your prepped dataset `data` from moveHMM::prepData()

library(moveHMM)
library(dplyr)
library(ggplot2)

set.seed(1234)

# --- Helper: wrap angles to (-pi, pi] ---
wrap <- function(a) ((a + pi) %% (2*pi)) - pi

# --- Helper: rebuild XY from step/turning-angle, starting at (x0,y0) ---
rebuild_xy <- function(df, x0, y0) {
  # df must contain 'step' and 'angle' and optionally 'state'
  n <- nrow(df)
  if (n == 0) return(df)
  bearing <- numeric(n)
  if (n > 1) for (i in 2:n) bearing[i] <- ((bearing[i-1] + df$angle[i] + pi) %% (2*pi)) - pi
  
  x <- numeric(n); y <- numeric(n)
  x[1] <- x0; y[1] <- y0
  if (n > 1) {
    dx <- df$step[-1] * cos(bearing[-1])
    dy <- df$step[-1] * sin(bearing[-1])
    x[-1] <- x0 + cumsum(dx)
    y[-1] <- y0 + cumsum(dy)
  }
  
  # keep all original columns (including state), add new x,y
  df$x <- x
  df$y <- y
  df
}

# 1) Per-ID lengths from the real (prepped) data + real start points
obs_per <- data %>%
  group_by(ID) %>%
  summarise(n = n(), .groups = "drop")

starts <- data %>%
  group_by(ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(ID, x0 = x, y0 = y)

# 2) Simulate one trajectory per real ID, with matching length
sim_list <- lapply(seq_len(nrow(obs_per)), function(i) {
  out <- simData(
    m = mod3,
    nbAnimals    = 1,
    obsPerAnimal = obs_per$n[i],  # scalar OK
    states = TRUE
  )
  out$ID <- as.character(obs_per$ID[i])  # label with the real ID
  # ensure temporal order within each ID
  out <- out %>% group_by(ID) %>% mutate(.ord = row_number()) %>% arrange(.ord) %>% ungroup()
  out
})
sim_raw <- bind_rows(sim_list)

# 3) Anchor each simulated track to the real start (x0,y0) and rebuild XY
sim_anchored <- sim_raw %>%
  group_by(ID) %>%
  group_split() %>%
  lapply(function(df) {
    s <- starts %>% filter(ID == unique(df$ID)) %>% slice(1)
    rebuild_xy(df, x0 = s$x0, y0 = s$y0)
  }) %>%
  bind_rows()

library(dplyr)
library(ggplot2)

# Ensure 'state' column exists (rename if needed) and make it label-ready
sim_plot <- sim_anchored %>%
  { if ("states" %in% names(.)) dplyr::rename(., state = states) else . } %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  )

# Pick 4 IDs (or sample(unique(sim_plot$ID), 4))
ids4 <- unique(sim_plot$ID)[5:8]
df4  <- sim_plot %>% filter(ID %in% ids4)

# Clean & keep groups with ≥ 2 points
df4_clean <- df4 %>%
  filter(is.finite(x), is.finite(y)) %>%
  group_by(ID) %>%
  mutate(n_pts = n()) %>%
  ungroup() %>%
  filter(n_pts >= 2)

# Map state -> requested legend labels (handles 1/2/3 or S1/S2/S3)
state_to_label <- function(s) {
  sn <- suppressWarnings(as.integer(as.character(s)))
  if (any(is.na(sn))) {
    sn <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(s))))
  }
  factor(
    sn,
    levels = c(1, 2, 3),
    labels = c("1. Reverse and Searching",
               "2. Reverse and Meandering",
               "3. Long Travel")
  )
}

df4_clean <- df4_clean %>%
  mutate(state_lbl = state_to_label(state))

# Recompute start/end AFTER cleaning
starts <- df4_clean %>% group_by(ID) %>% slice_head(n = 1) %>% ungroup() %>% mutate(point = "Start")
ends   <- df4_clean %>% group_by(ID) %>% slice_tail(n = 1) %>% ungroup() %>% mutate(point = "End")
pts    <- bind_rows(starts, ends)

# Plot: simulated only, per-panel free scales, no strip titles
ggplot(df4_clean, aes(x, y, group = ID)) +
  geom_path(aes(colour = state_lbl), linewidth = 0.6) +
  geom_point(data = pts, aes(x, y, shape = point),
             size = 2.2, colour = "black", fill = "white", stroke = 0.8) +
  facet_wrap(~ ID, ncol = 2, scales = "free") +
  scale_color_discrete(name = "State") +
  scale_shape_manual(values = c(Start = 21, End = 23)) +
  labs(x = "Axial Travel (µm)", y = "Transverse Travel (µm)") +  # no overall title if you don't want one
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_blank(),        # <-- remove per-panel titles
    strip.background = element_blank()
  )
############single##############################################################
library(dplyr)
library(ggplot2)

# ---------- Build plotting data from your simulated set ----------
# If your simulated object has 'states' (plural), normalize to 'state'
sim_plot <- sim_anchored %>%
  { if ("states" %in% names(.)) dplyr::rename(., state = states) else . } %>%
  mutate(
    x = as.numeric(x),
    y = as.numeric(y)
  ) %>%
  group_by(ID) %>% mutate(idx = dplyr::row_number()) %>% ungroup()

# State labels exactly as requested
state_to_label <- function(s) {
  sn <- suppressWarnings(as.integer(as.character(s)))
  if (any(is.na(sn))) sn <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(s))))
  factor(sn, levels = c(1,2,3), labels = c("State 1","State 2","State 3"))
}
sim_plot <- sim_plot %>% mutate(state_lbl = state_to_label(state))

# ---------- Pick ONE simulated ID to plot ----------
id_to_plot_sim <- as.character(unique(sim_plot$ID)[1])   # change or sample if you like
df1s <- sim_plot %>%
  filter(ID == id_to_plot_sim, is.finite(x), is.finite(y)) %>%
  arrange(idx)

# Build per-step segments to ensure continuity
segdf_s <- df1s %>%
  transmute(x, y, xend = dplyr::lead(x), yend = dplyr::lead(y), state_lbl) %>%
  filter(is.finite(x), is.finite(y), is.finite(xend), is.finite(yend))

# Start/End markers
starts_s <- df1s %>% slice_head(n = 1) %>% mutate(point = "Start")
ends_s   <- df1s %>% slice_tail(n = 1) %>% mutate(point = "End")
pts_s    <- dplyr::bind_rows(starts_s, ends_s)

# In-panel ID label position
xr <- diff(range(df1s$x, na.rm = TRUE)); yr <- diff(range(df1s$y, na.rm = TRUE))
lab_x <- min(df1s$x, na.rm = TRUE) + 0.02 * xr
lab_y <- max(df1s$y, na.rm = TRUE) - 0.02 * yr

# ---------- Plot (legends inside bottom-right) ----------
p_single_sim <- ggplot() +
  geom_segment(
    data = segdf_s,
    aes(x = x, y = y, xend = xend, yend = yend, colour = state_lbl),
    linewidth = 1.0, lineend = "round"
  ) +
  geom_point(
    data = pts_s, aes(x, y, shape = point),
    size = 2.6, colour = "black", fill = "white", stroke = 0.9
  ) +
  annotate("text", x = lab_x, y = lab_y,
           label = paste0("ID: ", id_to_plot_sim),
           hjust = 0, vjust = 1.1, size = 5, fontface = "bold") +
  scale_color_manual(
    name = "State",
    values = c("State 1" = "#0072B2", "State 2" = "#E69F00", "State 3" = "#009E73")
  ) +
  scale_shape_manual(name = "Point", values = c(Start = 21, End = 23)) +
  labs(x = "X (µm)", y = "Y (µm)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = c(0.98, 0.02),   # bottom-right inside
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.box.margin = margin(4, 6, 4, 6)
  )

print(p_single_sim)

# Optional: export as vector PDF
# ggsave("single_simulated_trajectory_state123_bottomright.pdf", p_single_sim,
#        device = cairo_pdf, width = 3.54, height = 2.8, units = "in",
#        dpi = 300, bg = "white", useDingbats = FALSE)


## =========================
## BASE R: Single simulated trajectory with state-colored step dots + grid
## =========================
getwd()
setwd("/Users/sayantanbiswas/Desktop/Movement Project/final plot")
pdf("spz_sim2.pdf", width = 5, height = 5, onefile = FALSE)
## pick simulated ID (change as needed)
simdat <- sim_anchored
if ("states" %in% names(simdat)) simdat$state <- simdat$states  # normalize name

id_to_plot_sim <- as.character(unique(simdat$ID)[6])

## subset & clean
dat_sim <- subset(simdat, ID == id_to_plot_sim)
x <- as.numeric(dat_sim$x); y <- as.numeric(dat_sim$y)
ok <- is.finite(x) & is.finite(y)
x <- x[ok]; y <- y[ok]
state <- as.integer(dat_sim$state[ok])
n <- length(x); stopifnot(n >= 2)

## palette + labels
pal <- c("State 1" = "#0072B2", "State 2" = "#E69F00", "State 3" = "#009E73")
state_lbl <- factor(state, levels = c(1,2,3), labels = names(pal))
seg_col <- unname(pal[state_lbl[1:(n-1)]])  # color of step i → i+1

## axes ranges
xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
dx <- diff(xr); dy <- diff(yr)

par(xaxs = "i", yaxs = "i")
# (A) auto limits:
plot(NA, NA, xlim = xr, ylim = yr, xlab = "X (µm)", ylab = "Y (µm)", asp = 1, cex=1 )
# (B) or force specific window like your real plot:
# plot(NA, NA, xlim = c(80,92), ylim = c(84,96), xlab = "X (µm)", ylab = "Y (µm)", asp = 1)

## background grid
abline(v = axTicks(1), col = "darkgray",  lty = "dotted", lwd = 0.8)
abline(h = axTicks(2), col = "darkgrey", lty = "dotted", lwd = 0.8)

## continuous path (segments by state)
segments(x[1:(n-1)], y[1:(n-1)], x[2:n], y[2:n], col = seg_col, lwd = 2)

## tiny step-end dots (same state color; no legend)
points(x[2:n], y[2:n], pch = 16, col = seg_col, cex = 0.5)

## start / end markers
points(x[1], y[1], pch = 21, bg = "white", col = "black", cex = 1.2, lwd = 1.2)
points(x[n], y[n], pch = 24, bg = "white", col = "black", cex = 1.2, lwd = 1.2)

## separate legends
#legend("top", inset = 0.02,
#       legend = paste0("ID: ", id_to_plot_sim),
#       bty = "n", text.font = 2, cex = 1.1)
#
#legend("left", inset = 0.00,
#       legend = "Time post inoculation = 5 mins",
#       bty = "n", text.font = 2, cex = 0.8)

legend("bottomright", bty = "n",
       legend = c(names(pal), "Start", "End"),
       col    = c(unname(pal), "black", "black"),
       lty    = c(1,1,1, NA, NA),
       lwd    = c(2,2,2, NA, NA),
       pch    = c(NA,NA,NA, 21, 24),
       pt.bg  = c(NA,NA,NA, "white", "white"),
       pt.cex = c(1,1,1, 1.2, 1.2),
       seg.len = 2)

box()  # border
dev.off()
########### origin centred plot#################################################
## recenter trajectory at (0,0)
x <- x - x[1]
y <- y - y[1]

## symmetric limits around 0
xr <- range(x, na.rm = TRUE)
yr <- range(y, na.rm = TRUE)
max_range <- max(abs(c(xr, yr)))   # largest absolute extent
lim <- c(-max_range, max_range)   # symmetric about 0

par(xaxs = "i", yaxs = "i")
plot(NA, NA, xlim = lim, ylim = lim,
     xlab = "X displacement (µm)", ylab = "Y displacement (µm)",
     asp = 1)

## add grid
abline(v = axTicks(1), col = "darkgray", lty = "dotted", lwd = 0.8)
abline(h = axTicks(2), col = "darkgrey", lty = "dotted", lwd = 0.8)

## add central axes (through 0,0)
abline(h = 0, col = "black", lwd = 0.8)
abline(v = 0, col = "black", lwd = 0.8)

## trajectory segments
segments(x[1:(n-1)], y[1:(n-1)], x[2:n], y[2:n],
         col = seg_col, lwd = 1.5)

## dots
points(x[2:n], y[2:n], pch = 16, col = seg_col, cex = 0.4)

## start & end
points(0, 0, pch = 21, bg = "white", col = "black", cex = 1.2, lwd = 1.2) # start at (0,0)
points(x[n], y[n], pch = 24, bg = "white", col = "black", cex = 1.2, lwd = 1.2)

## legends (unchanged)
legend("top", inset = 0.02,
       legend = paste0("ID: ", id_to_plot_sim),
       bty = "n", text.font = 2, cex = 1.1)
legend("left", inset = 0.00,
       legend = "Time post inoculation = 5 mins",
       bty = "n", text.font = 2, cex = 0.8)
legend("bottomright", bty = "n",
       legend = c(names(pal), "Start", "End"),
       col    = c(unname(pal), "black", "black"),
       lty    = c(1,1,1, NA, NA),
       lwd    = c(2,2,2, NA, NA),
       pch    = c(NA,NA,NA, 21, 24),
       pt.bg  = c(NA,NA,NA, "white", "white"),
       pt.cex = c(1,1,1, 1.2, 1.2),
       seg.len = 2)

box()

