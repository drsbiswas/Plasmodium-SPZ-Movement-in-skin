

# Load Libraries

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggExtra)
library(openxlsx)
library(effsize)


# Read file

data_original <- read_excel("~/Downloads/Sporozoites_Master_File_Hopp_e21_ver3_EH.xlsx")
getwd()
setwd("~/Desktop/Movement Project/main text")
getwd()
data_original$time_min_num <- suppressWarnings(
  as.numeric(filtered_data_dplyr$`Time Post-Inoculation (min)`)
)

# timewise pf/py -----------------------------------------
#data_original <- data_original %>%
#  filter(
#    tolower(Species) == "py" #,
#!is.na(time_min_num),
#time_min_num == 120
#  )

# Calculate Final Displacement per track

# -------- 1) Helper: drop "buzzing"/non-moving tracks by spatial spread --------
# Keeps tracks whose overall spread sqrt((maxX-minX)^2 + (maxY-minY)^2) >= threshold_um
# track_col: "Custom TrackID" or "Clean TrackIDs"
# x_col/y_col: "Position X (um)"/"Position Y (um)" in your files
n=5
filter_tracks_by_spread <- function(data,
                                    track_col,
                                    x_col = "Position X (um)",
                                    y_col = "Position Y (um)",
                                    threshold_um = n,
                                    min_points = 1 ,
                                    verbose = TRUE) {
  stopifnot(all(c(track_col, x_col, y_col) %in% names(data)))
  
  track_stats <- data %>%
    dplyr::mutate(track_chr = as.character(.data[[track_col]])) %>%
    dplyr::group_by(track_chr) %>%
    dplyr::summarise(
      n_points  = dplyr::n(),
      dx_um     = max(.data[[x_col]], na.rm = TRUE) - min(.data[[x_col]], na.rm = TRUE),
      dy_um     = max(.data[[y_col]], na.rm = TRUE) - min(.data[[y_col]], na.rm = TRUE),
      spread_um = sqrt(dx_um^2 + dy_um^2),
      .groups   = "drop"
    )
  
  keep_ids <- track_stats %>%
    dplyr::filter(is.finite(spread_um),
                  spread_um >= threshold_um,
                  n_points >= min_points) %>%
    dplyr::pull(track_chr)
  
  drop_ids <- setdiff(unique(as.character(data[[track_col]])), keep_ids)
  
  if (isTRUE(verbose)) {
    total_ids <- length(keep_ids) + length(drop_ids)
    message(sprintf(
      "Track spread filter (>= %.2f um): kept %d / %d IDs (removed %d = %.1f%%).",
      threshold_um, length(keep_ids), total_ids, length(drop_ids),
      if (total_ids > 0) 100 * length(drop_ids) / total_ids else 0
    ))
    if (length(drop_ids) > 0)
      message(sprintf("Removed IDs (examples): %s",
                      paste(utils::head(drop_ids, 10), collapse = ", ")))
  }
  
  out <- data %>%
    dplyr::mutate(track_chr = as.character(.data[[track_col]])) %>%
    dplyr::filter(track_chr %in% keep_ids) %>%
    dplyr::select(-track_chr)
  
  # Safety check
  stopifnot(length(unique(as.character(out[[track_col]]))) == length(keep_ids))
  out
}


threshold_um <- n  # or 3, per your note

# Filter the raw long table once per tracking scheme
data_custom_raw <- filter_tracks_by_spread(
  data_original, track_col = "Custom TrackID",
  threshold_um = threshold_um
)
data_clean_raw  <- filter_tracks_by_spread(
  data_original, track_col = "Clean TrackIDs",
  threshold_um = threshold_um
)
###############################################################################
# Compute duration (max time - min time) for each Clean TrackID
track_duration <- data_clean_raw %>%
  group_by(`Clean TrackIDs`) %>%
  summarise(duration_s = max(`Time (s)`, na.rm = TRUE) - 
              min(`Time (s)`, na.rm = TRUE),
            .groups = "drop")

# Plot histogram of track duration vs frequency
ggplot(track_duration, aes(x = duration_s)) +
  geom_histogram(binwidth = 10,   # adjust binwidth as needed
                 fill = "steelblue", color = "black") +
  labs(title = "Histogram of Track Duration vs Frequency",
       x = "Track duration (s)",
       y = "Frequency") +
  theme_minimal()
# Calculate mean track duration
mean_duration <- mean(track_duration$duration_s, na.rm = TRUE)
print(paste("Mean track duration:", round(mean_duration, 2), "seconds"))

# (If you also plot instantaneous/mean speeds or turning angles,
#  feed their compute_* functions with data_custom_raw / data_clean_raw the same way.)


compute_final_displacement_strict <- function(data, track_col) {
  data %>%
    group_by(Species, `Time Post-Inoculation (min)`, !!sym(track_col)) %>%
    arrange(`Time (s)`, .by_group = TRUE) %>%
    mutate(Time_Diff = `Time (s)` - lag(`Time (s)`)) %>%
    filter(all(Time_Diff > 0, na.rm=TRUE)) %>%
    summarise(
      Final_Displacement = sqrt(
        (last(`Position X (um)`) - first(`Position X (um)`))^2 +
          (last(`Position Y (um)`) - first(`Position Y (um)`))^2
      ),
      .groups = "drop"
    ) %>%
    filter(Final_Displacement > 0)
}

data_trackid    <- compute_final_displacement_strict(data_original, "Custom TrackID")
data_newtrackid <- compute_final_displacement_strict(data_original, "Clean TrackIDs")
# Then recompute your summaries on the FILTERED data:
data_trackid_fd    <- compute_final_displacement_strict(data_custom_raw, "Custom TrackID")
data_newtrackid_fd <- compute_final_displacement_strict(data_clean_raw,  "Clean TrackIDs")


# Merge all data

add_tracking_method <- function(df, label) df %>% mutate(TrackingMethod = label)
combined_data <- bind_rows(
  add_tracking_method(data_trackid_fd, "Custom"),
  add_tracking_method(data_newtrackid_fd, "Clean")
) %>% rename(Time = `Time Post-Inoculation (min)`)


# T-Tests & Cohen's d

all_timepoints <- c(5,10,20,30,60,120)
geometric_mean <- function(x) if(length(x)) exp(mean(log(x), na.rm=TRUE)) else NA_real_

perform_ttest <- function(df, f1, f2, lab1, lab2, tp) {
  g1 <- df %>% filter(!!f1, Time==tp) %>% pull(Final_Displacement)
  g2 <- df %>% filter(!!f2, Time==tp) %>% pull(Final_Displacement)
  if(length(g1)<2 || length(g2)<2) {
    return(tibble(
      Time=tp, Comparison=paste(lab1,"vs",lab2),
      N_group1=length(g1), Mean_group1=ifelse(length(g1),geometric_mean(g1),NA),
      N_group2=length(g2), Mean_group2=ifelse(length(g2),geometric_mean(g2),NA),
      t_stat=NA, p_value=NA, cohen_d=NA
    ))
  }
  tt <- t.test(g1,g2)
  d  <- cohen.d(g1,g2,paired=FALSE)$estimate
  tibble(
    Time=tp, Comparison=paste(lab1,"vs",lab2),
    N_group1=length(g1), Mean_group1=geometric_mean(g1),
    N_group2=length(g2), Mean_group2=geometric_mean(g2),
    t_stat=tt$statistic[[1]], p_value=tt$p.value, cohen_d=d
  )
}

ttests <- list()
for(tp in all_timepoints) {
  ttests[[length(ttests)+1]] <- bind_rows(
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  "PF clean","PY clean",tp),
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PF custom","PY custom",tp),
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  "PF clean","PF custom",tp),
    perform_ttest(combined_data,
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PY clean","PY custom",tp)
  )
}
ttest_df <- bind_rows(ttests)


# Linear models & effect sizes

species_list <- c("PF","PY"); track_list <- c("Clean","Custom")
lm_results   <- list()
for(sp in species_list) for(tm in track_list) {
  df_sub <- combined_data %>% filter(Species==sp, TrackingMethod==tm, Final_Displacement>0)
  if(nrow(df_sub)<2 || length(unique(df_sub$Time))<2) {
    lm_results[[length(lm_results)+1]] <- tibble(
      Species=sp, TrackingMethod=tm,
      slope=NA, p_value=NA, r_value=NA, r_sq=NA, cohen_f2=NA, n=nrow(df_sub)
    )
  } else {
    m   <- lm(log(Final_Displacement) ~ Time, data=df_sub)
    sm  <- summary(m)
    r2  <- sm$r.squared
    f2  <- if(!is.na(r2)&&r2<1) r2/(1-r2) else NA
    lm_results[[length(lm_results)+1]] <- tibble(
      Species=sp, TrackingMethod=tm,
      slope=sm$coefficients[2,1],
      p_value=sm$coefficients[2,4],
      r_value=cor(log(df_sub$Final_Displacement), df_sub$Time),
      r_sq=r2, cohen_f2=f2, n=nrow(df_sub)
    )
  }
}
lm_df <- bind_rows(lm_results)


# Generate Excel file with stats

wb <- createWorkbook()
addWorksheet(wb,"Pairwise_Tests"); addWorksheet(wb,"Linear_Models")
writeData(wb,"Pairwise_Tests",ttest_df); writeData(wb,"Linear_Models",lm_df)
saveWorkbook(wb,"Stat_tests_final_displacement.xlsx",overwrite=TRUE)
cat("Stats written to Excel.\n")


# Plot with dashed regression line and annotation

generate_plot <- function(data_summarized, species_label, track_label, ann_label) {
  df <- data_summarized %>%
    filter(Species==species_label) %>%
    rename(Time=`Time Post-Inoculation (min)`) %>%
    filter(Final_Displacement>0) %>%
    mutate(
      Time_str = factor(paste0(Time," min"),
                        levels=c("5 min","10 min","20 min","30 min","60 min","120 min")),
      Time_num = as.numeric(Time_str)
    )
  
  stats <- df %>%
    group_by(Time_str, Time_num) %>%
    summarise(
      geom_mean = exp(mean(log(Final_Displacement))),
      Count     = n(),
      .groups   = "drop"
    )
  
  # Extended regression line
  fit_means <- lm(log(geom_mean) ~ Time_num, data=stats)
  x0 <- min(stats$Time_num)-0.5; x1 <- max(stats$Time_num)+0.6
  y0 <- exp(predict(fit_means,newdata=data.frame(Time_num=min(stats$Time_num))))
  y1 <- exp(predict(fit_means,newdata=data.frame(Time_num=max(stats$Time_num))))
  line_df <- tibble(x=x0,xend=x1,y=y0,yend=y1)
  
  # Compute raw model stats for Excel
  raw_model <- lm(log(Final_Displacement) ~ Time, data=df)
  rm_sum    <- summary(raw_model)
  slope_raw <- rm_sum$coefficients[2,1]
  p_raw     <- rm_sum$coefficients[2,4]
  r2_raw    <- rm_sum$r.squared
  f2_raw    <- if(!is.na(r2_raw)&&r2_raw<1) r2_raw/(1-r2_raw) else NA
  
  slope_txt <- ifelse(abs(slope_raw)<0.001,
                      formatC(slope_raw,format="e",digits=1),
                      sprintf("%.3f",slope_raw))
  p_txt     <- ifelse(p_raw<0.001,
                      formatC(p_raw,format="e",digits=1),
                      sprintf("%.3f",p_raw))
  f2_txt    <- ifelse(f2_raw<0.001,
                      formatC(f2_raw,format="e",digits=1),
                      sprintf("%.3f",f2_raw))
  
  slope_num <- as.numeric(slope_txt)
  p_num     <- as.numeric(p_txt)
  f2_num    <- as.numeric(f2_txt)
  
  stats_expr <- substitute(
    "(         " ~ italic(m)==M ~ "," ~ italic(p)==P ~ ", "~
      "Cohen's " ~ italic(f)^2==F2 ~ ")",
    list(M=slope_num, P=p_num, F2=f2_num)
  )
  
  color_map <- c(
    "5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
    "30 min"="darkgreen","60 min"="blue","120 min"="#E75480"
  )
  
  ggplot(df, aes(x=Time_str,y=Final_Displacement,fill=Time_str))+
    geom_violin(alpha=0.3,scale="width",width=0.9,color="black")+
    scale_y_log10(
      limits=c(0.007,4000),
      breaks=c(0.01,0.1,1,10,100,1000),
      labels=scales::trans_format("log10", scales::math_format(10^.x)),
      expand=c(0,0)
    )+
    scale_fill_manual(values=color_map)+
    geom_segment(
      data=stats,
      aes(x=Time_num-0.4,xend=Time_num+0.4,
          y=geom_mean,yend=geom_mean),
      inherit.aes=FALSE,color="black",linewidth=1.5
    )+
    geom_segment(
      data=line_df,
      aes(x=x,xend=xend,y=y,yend=yend),
      inherit.aes=FALSE,linetype="dashed",size=2,color="blue"
    )+
    annotate(
      "text",
      x=stats$Time_num,y=400,
      label=paste0("µ = ", round(stats$geom_mean,1)," µm"),
      size=8,hjust=0.5
    )+
    annotate(
      "text",
      x=stats$Time_num,y=200,
      label=paste0("n = ",stats$Count),
      size=8,hjust=0.5
    )+
    
    annotate(
      "text",
      x=0.5,y=2200,
      label=ann_label,
      size=14,hjust=0, fontface="italic"
    )+
    
    annotate("segment", x=0.56, xend=0.73, y=900, yend=900,
             colour="blue", linewidth=1.6, linetype="longdash")+
    
    annotate(
      "text",
      x=0.5+0.08,y=900,
      label=stats_expr,
      parse=FALSE, size=10, hjust=0
    )+
    labs(x="Time Post-Inoculation (min)",y="Final Displacement (µm)")+
    theme_minimal()+
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      axis.line=element_line(colour="black",size=1.25),
      axis.text.x=element_text(size=26,color="black"),
      axis.text.y=element_text(size=26,color="black"),
      axis.title=element_text(size=26),
      axis.ticks.x.bottom=element_line(color="black", size = 1.25),
      axis.ticks.length=unit(7,"pt"),
      axis.ticks.y.left=element_line(color="black", size = 1.25),
      legend.position="none"
    )
}

# Make and save
pf_c <- generate_plot(data_trackid_fd, "PF", "Custom", expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_c <- generate_plot(data_trackid_fd, "PY", "Custom", expression("P. yoelii (f" [threshold] * "= 5µm)"))
pf_c0 <- generate_plot(data_trackid, "PF", "Custom", "P. falciparum")
py_c0 <- generate_plot(data_trackid, "PY", "Custom", "P. yoelii")
panel <- plot_grid(pf_c0, py_c0,pf_c,py_c, ncol=2, nrow= 2, labels=c("A","B","C","D"), label_size=40)
ggsave("plot_CustomTrackID_Pf_Py.pdf",panel,width=28,height=20)

pf_cl <- generate_plot(data_newtrackid_fd, "PF", "Clean", expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_cl <- generate_plot(data_newtrackid_fd, "PY", "Clean", expression("P. yoelii (f" [threshold] * "= 5µm)"))
pf_cl0 <- generate_plot(data_newtrackid, "PF", "Clean", "P. falciparum")
py_cl0 <- generate_plot(data_newtrackid, "PY", "Clean", "P. yoelii")
panel2 <- plot_grid(pf_cl0,py_cl0,pf_cl,py_cl,ncol=2, nrow= 2, labels=c("A","B","C","D"),label_size=40)
ggsave("plot_CleanTrackIDs_Pf_Py.pdf",panel2,width=28,height=20)

cat("All done—plots annotations now match Excel output in the same style.\n")










################################################################
################################################################
################################################################








# Load Libraries

#library(readxl)
#library(dplyr)
#library(ggplot2)
#library(cowplot)
#library(ggExtra)
#library(openxlsx)
#library(effsize)


# Read file

#data_original <- read_excel("Hopp_PF-PY_Master_file.xlsx")


# Calculate Instantaneous Speed

compute_instant_speeds <- function(data, track_col) {
  data %>%
    group_by(Species, `Time Post-Inoculation (min)`, !!sym(track_col)) %>%
    arrange(`Time (s)`, .by_group = TRUE) %>%
    mutate(
      Distance  = sqrt((`Position X (um)` - lag(`Position X (um)`))^2 +
                         (`Position Y (um)` - lag(`Position Y (um)`))^2),
      Time_Diff = `Time (s)` - lag(`Time (s)`)
    ) %>%
    filter(!is.na(Distance), !is.na(Time_Diff)) %>%
    filter(all(Time_Diff > 0, na.rm=TRUE)) %>%
    mutate(Instantaneous_Speed = Distance / Time_Diff) %>%
    ungroup() %>%
    filter(Instantaneous_Speed > 0)
}

data_trackid_is    <- compute_instant_speeds(data_custom_raw, "Custom TrackID")
data_newtrackid_is <- compute_instant_speeds(data_clean_raw, "Clean TrackIDs")



# Combine both data sets for stat tests


add_tracking_method <- function(df, method_label) {
  df %>% mutate(TrackingMethod = method_label)
}

data_custom_is <- add_tracking_method(data_trackid_is,    "Custom")
data_clean_is  <- add_tracking_method(data_newtrackid_is, "Clean")

combined_data_is <- bind_rows(data_custom_is, data_clean_is) %>%
  rename(Time = `Time Post-Inoculation (min)`)


# T-tests and Cohen's d

all_timepoints <- c(5,10,20,30,60,120)
geometric_mean <- function(x) if(length(x)) exp(mean(log(x), na.rm=TRUE)) else NA_real_

perform_ttest <- function(df, f1, f2, lab1, lab2, tp) {
  g1 <- df %>% filter(!!f1, Time==tp) %>% pull(Instantaneous_Speed)
  g2 <- df %>% filter(!!f2, Time==tp) %>% pull(Instantaneous_Speed)
  if(length(g1)<2||length(g2)<2) {
    return(tibble(
      Time=tp,
      Comparison=paste(lab1,"vs",lab2),
      N_group1=length(g1), Mean_group1=ifelse(length(g1), geometric_mean(g1), NA),
      N_group2=length(g2), Mean_group2=ifelse(length(g2), geometric_mean(g2), NA),
      t_stat=NA, p_value=NA, cohen_d=NA
    ))
  }
  tt <- t.test(g1,g2)
  d  <- cohen.d(g1,g2,paired=FALSE)$estimate
  tibble(
    Time=tp,
    Comparison=paste(lab1,"vs",lab2),
    N_group1=length(g1), Mean_group1=geometric_mean(g1),
    N_group2=length(g2), Mean_group2=geometric_mean(g2),
    t_stat=tt$statistic[[1]], p_value=tt$p.value, cohen_d=d
  )
}

ttest_results <- list()
for(tp in all_timepoints) {
  ttest_results[[length(ttest_results)+1]] <- bind_rows(
    perform_ttest(combined_data_is,
                  quote(Species=="PF" & TrackingMethod=="Clean"),
                  quote(Species=="PY" & TrackingMethod=="Clean"),
                  "PF clean","PY clean",tp),
    perform_ttest(combined_data_is,
                  quote(Species=="PF" & TrackingMethod=="Custom"),
                  quote(Species=="PY" & TrackingMethod=="Custom"),
                  "PF custom","PY custom",tp),
    perform_ttest(combined_data_is,
                  quote(Species=="PF" & TrackingMethod=="Clean"),
                  quote(Species=="PF" & TrackingMethod=="Custom"),
                  "PF clean","PF custom",tp),
    perform_ttest(combined_data_is,
                  quote(Species=="PY" & TrackingMethod=="Clean"),
                  quote(Species=="PY" & TrackingMethod=="Custom"),
                  "PY clean","PY custom",tp)
  )
}
ttest_df <- bind_rows(ttest_results)


# Linear Regression per species & tracking method

species_list       <- c("PF", "PY")
trackingmethod_lst <- c("Clean", "Custom")
lm_results <- list()

for(sp in species_list) {
  for(tm in trackingmethod_lst) {
    subdat <- combined_data_is %>% filter(Species==sp, TrackingMethod==tm)
    if(nrow(subdat) < 2) {
      lm_results[[length(lm_results)+1]] <- tibble(
        Species=sp, TrackingMethod=tm,
        slope=NA, p_value=NA, r_value=NA, r_sq=NA,
        cohen_f2=NA, n=nrow(subdat)
      )
      next
    }
    group_means <- subdat %>%
      group_by(Time) %>%
      summarise(geom=exp(mean(log(Instantaneous_Speed), na.rm=TRUE)), .groups="drop") %>%
      filter(!is.na(geom), geom>0)
    if(nrow(group_means) < 2) {
      lm_results[[length(lm_results)+1]] <- tibble(
        Species=sp, TrackingMethod=tm,
        slope=NA, p_value=NA, r_value=NA, r_sq=NA,
        cohen_f2=NA, n=nrow(subdat)
      )
      next
    }
    model <- lm(log(geom) ~ Time, data=group_means)
    sm    <- summary(model)
    r2    <- sm$r.squared
    f2    <- if(!is.na(r2)&&r2<1) r2/(1-r2) else NA
    lm_results[[length(lm_results)+1]] <- tibble(
      Species       = sp,
      TrackingMethod= tm,
      slope         = sm$coefficients[2,1],
      p_value       = sm$coefficients[2,4],
      r_value       = cor(log(group_means$geom), group_means$Time),
      r_sq          = r2,
      cohen_f2      = f2,
      n             = nrow(subdat)
    )
  }
}
lm_df <- bind_rows(lm_results)


# Export all results to a single excel file

wb_out <- createWorkbook()
addWorksheet(wb_out, "Pairwise_Tests")
addWorksheet(wb_out, "Linear_Models")

writeData(wb_out, "Pairwise_Tests", ttest_df)
writeData(wb_out, "Linear_Models", lm_df)

saveWorkbook(wb_out, "Stat_tests_instantaneous_speeds.xlsx", overwrite=TRUE)
cat("\nAll T-tests, Cohen's d, linear regressions (with R^2, Cohen's f^2), etc.\n",
    "written to 'Stat_tests_instantaneous_speeds.xlsx'.\n")


# Plot

generate_plot_instant <- function(
    data_instant,
    species_label="PF",
    track_id_type="Clean",
    annotation_label=""
) {
  data_for_plot <- data_instant %>%
    filter(Species == species_label, TrackingMethod == track_id_type) %>%
    mutate(
      Time_str     = factor(paste0(Time, " min"),
                            levels = c("5 min","10 min","20 min",
                                       "30 min","60 min","120 min")),
      Time_numeric = as.numeric(Time_str)
    ) %>%
    filter(!is.na(Instantaneous_Speed), Instantaneous_Speed > 0)
  
  group_stats <- data_for_plot %>%
    group_by(Time_str, Time_numeric) %>%
    summarise(
      geom_mean = exp(mean(log(Instantaneous_Speed), na.rm=TRUE)),
      Count     = n(),
      .groups   = "drop"
    )
  
  # Extended dashed‐line
  fit_means <- lm(log(geom_mean) ~ Time_numeric, data=group_stats)
  x0 <- min(group_stats$Time_numeric) - 0.5
  x1 <- max(group_stats$Time_numeric) + 0.6
  y0 <- exp(predict(fit_means, newdata=data.frame(Time_numeric=min(group_stats$Time_numeric))))
  y1 <- exp(predict(fit_means, newdata=data.frame(Time_numeric=max(group_stats$Time_numeric))))
  line_df <- tibble(x=x0, xend=x1, y=y0, yend=y1)
  
  
  match_row <- lm_df %>%
    filter(Species == species_label, TrackingMethod == track_id_type)
  slope_reg <- match_row$slope
  p_reg     <- match_row$p_value
  f2_reg    <- match_row$cohen_f2
  
  
  slope_txt <- ifelse(abs(slope_reg) < 0.001,
                      formatC(slope_reg, format="e", digits=1),
                      sprintf("%.3f", slope_reg))
  p_txt     <- ifelse(p_reg < 0.001,
                      formatC(p_reg, format="e", digits=1),
                      sprintf("%.3f", p_reg))
  f2_txt    <- ifelse(f2_reg < 0.001,
                      formatC(f2_reg, format="e", digits=1),
                      sprintf("%.3f", f2_reg))
  
  
  slope_num <- as.numeric(slope_txt)
  p_num     <- as.numeric(p_txt)
  f2_num    <- as.numeric(f2_txt)
  
  
  stats_expr <- substitute(
    "(         " ~ italic(m)==M ~ "," ~ italic(p)==P ~ ", "~
      "Cohen's " ~ italic(f)^2==F2 ~ ")",
    list(M = slope_num, P = p_num, F2 = f2_num)
  )
  
  color_map <- c(
    "5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
    "30 min"="darkgreen","60 min"="blue","120 min"="#E75480"
  )
  
  p <- ggplot(data_for_plot, aes(x=Time_str, y=Instantaneous_Speed)) +
    geom_violin(aes(fill=Time_str), alpha=0.3, scale="width", width=0.9, color="black") +
    scale_x_discrete(labels = function(x) gsub(" min","", x)) +
    scale_y_log10(
      limits=c(0.0007,1000),
      breaks = c(0.001,0.01,0.1,1,10,100),
      expand=c(0,0),
      labels=scales::trans_format("log10", scales::math_format(10^.x)),
    ) +
    scale_fill_manual(values=color_map) +
    geom_segment(
      data=group_stats,
      aes(x=Time_numeric-0.4, xend=Time_numeric+0.4,
          y=geom_mean,       yend=geom_mean),
      color="black", linewidth=1.5
    ) +
    geom_segment(
      data=line_df,
      aes(x=x, xend=xend, y=y, yend=yend),
      linetype="dashed", size=2, color="blue"
    ) +
    annotate(
      "text",
      x=group_stats$Time_numeric, y=50,
      label=paste0("µ = ", round(group_stats$geom_mean,2), " µm/s"),
      size=7.5, hjust=0.5
    ) +
    annotate(
      "text",
      x=group_stats$Time_numeric, y=20,
      label=paste0("n = ", group_stats$Count),
      size=7.5, hjust=0.5
    ) +
    
    annotate(
      "text",
      x=0.5, y=500,
      label=annotation_label,
      size=14, hjust=0, fontface="italic"
    ) +
    
    annotate("segment", x=0.56, xend=0.73, y=150, yend=150,
             colour="blue", linewidth=1.8, linetype="longdash")+
    
    annotate(
      "text",
      x=0.58, y=150,
      label=stats_expr,
      parse=FALSE, size=12, hjust=0
    ) +
    labs(
      x="Time Post-Inoculation (min)",
      y="Instantaneous Speed (µm/s)"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.line           = element_line(colour="black", size=1.25),
      axis.text.x         = element_text(size=26, color="black"),
      axis.text.y         = element_text(size=26, color="black"),
      axis.ticks          = element_line(size=1.25),
      axis.ticks.length   = unit(7, "pt"),
      axis.title          = element_text(size=26),
      legend.position     = "none"
    )
  
  return(p)
}


# Make plots for custom and clean

pf_custom_plot <- generate_plot_instant(
  data_instant     = combined_data_is,
  species_label    = "PF",
  track_id_type    = "Custom",
  annotation_label = "P. falciparum"
)

py_custom_plot <- generate_plot_instant(
  data_instant     = combined_data_is,
  species_label    = "PY",
  track_id_type    = "Custom",
  annotation_label = "P. yoelii"
)

plot_custom_panel <- plot_grid(
  pf_custom_plot, py_custom_plot,
  ncol=2, labels=c("A","B"), label_size=40, align="hv"
)
ggsave("InstantSpeed_Custom_Pf_Py.pdf", plot_custom_panel, width=28, height=9)
cat("Saved instantaneous speed (custom) side-by-side plot.\n")

pf_clean_plot <- generate_plot_instant(
  data_instant     = combined_data_is,
  species_label    = "PF",
  track_id_type    = "Clean",
  annotation_label = "P. falciparum"
)

py_clean_plot <- generate_plot_instant(
  data_instant     = combined_data_is,
  species_label    = "PY",
  track_id_type    = "Clean",
  annotation_label = "P. yoelii"
)

plot_clean_panel <- plot_grid(
  pf_clean_plot, py_clean_plot,
  ncol=2, labels=c("A","B"), label_size=40, align="hv"
)
ggsave("InstantSpeed_Clean_Pf_Py.pdf", plot_clean_panel, width=28, height=9)
cat("Saved instantaneous speed (clean) side-by-side plot.\n")

cat("\nAll done with Instantaneous Speed Analysis!\n")







################################################################
################################################################
################################################################










# Load Libraries

#library(readxl)
#library(dplyr)
#library(ggplot2)
#library(cowplot)
#library(ggExtra)
#library(openxlsx)
#library(effsize)


# Read file

#data_original <- read_excel("Hopp_PF-PY_Master_file.xlsx")


# Calculate mean speed

compute_mean_speeds_strict <- function(data, track_col) {
  data %>%
    group_by(Species, `Time Post-Inoculation (min)`, !!sym(track_col)) %>%
    arrange(`Time (s)`, .by_group = TRUE) %>%
    mutate(
      Distance   = sqrt((`Position X (um)` - lag(`Position X (um)`))^2 +
                          (`Position Y (um)` - lag(`Position Y (um)`))^2),
      Time_Diff  = `Time (s)` - lag(`Time (s)`)
    ) %>%
    filter(!is.na(Distance), !is.na(Time_Diff)) %>%
    filter(all(Time_Diff > 0)) %>%
    mutate(Instantaneous_Speed = Distance / Time_Diff) %>%
    summarise(
      Mean_Speed = mean(Instantaneous_Speed, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    filter(Mean_Speed > 0)
}

data_trackid_ms    <- compute_mean_speeds_strict(data_custom_raw, "Custom TrackID")
data_newtrackid_ms <- compute_mean_speeds_strict(data_clean_raw, "Clean TrackIDs")
data_trackid_ms0    <- compute_mean_speeds_strict(data_original, "Custom TrackID")
data_newtrackid_ms0 <- compute_mean_speeds_strict(data_original, "Clean TrackIDs")


# Combine both data sets for stat tests

add_tracking_method <- function(df, method_label) {
  df %>% mutate(TrackingMethod = method_label)
}

data_custom <- add_tracking_method(data_trackid_ms,    "Custom")
data_clean  <- add_tracking_method(data_newtrackid_ms, "Clean")
data_custom0 <- add_tracking_method(data_trackid_ms0,    "Custom")
data_clean0  <- add_tracking_method(data_newtrackid_ms0, "Clean")

combined_data <- bind_rows(data_custom, data_clean) %>%
  rename(Time = `Time Post-Inoculation (min)`)
combined_data0 <- bind_rows(data_custom0, data_clean0) %>%
  rename(Time = `Time Post-Inoculation (min)`)

# T-tests and Cohen's d

all_timepoints <- c(5,10,20,30,60,120)
geometric_mean <- function(x) if(length(x)==0) NA_real_ else exp(mean(log(x), na.rm=TRUE))

perform_ttest <- function(df, f1, f2, lab1, lab2, tp) {
  g1 <- df %>% filter(!!f1, Time==tp) %>% pull(Mean_Speed)
  g2 <- df %>% filter(!!f2, Time==tp) %>% pull(Mean_Speed)
  if(length(g1)<2 || length(g2)<2) {
    return(tibble(
      Time        = tp,
      Comparison  = paste(lab1, "vs", lab2),
      N_group1    = length(g1),
      Mean_group1 = ifelse(length(g1)==0, NA, geometric_mean(g1)),
      N_group2    = length(g2),
      Mean_group2 = ifelse(length(g2)==0, NA, geometric_mean(g2)),
      t_stat      = NA, p_value = NA, cohen_d = NA
    ))
  }
  tt <- t.test(g1, g2)
  d  <- cohen.d(g1, g2, paired=FALSE)$estimate
  tibble(
    Time        = tp,
    Comparison  = paste(lab1, "vs", lab2),
    N_group1    = length(g1),
    Mean_group1 = geometric_mean(g1),
    N_group2    = length(g2),
    Mean_group2 = geometric_mean(g2),
    t_stat      = tt$statistic[[1]],
    p_value     = tt$p.value,
    cohen_d     = d
  )
}

ttest_results <- list()
for(tp in all_timepoints) {
  ttest_results[[length(ttest_results)+1]] <- bind_rows(
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  "PF clean","PY clean",tp),
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PF custom","PY custom",tp),
    perform_ttest(combined_data,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  "PF clean","PF custom",tp),
    perform_ttest(combined_data,
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PY clean","PY custom",tp)
  )
}
ttest_df <- bind_rows(ttest_results)

ttest_results0 <- list()
for(tp in all_timepoints) {
  ttest_results[[length(ttest_results)+1]] <- bind_rows(
    perform_ttest(combined_data0,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  "PF clean","PY clean",tp),
    perform_ttest(combined_data0,
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PF custom","PY custom",tp),
    perform_ttest(combined_data0,
                  quote(Species=="PF"&TrackingMethod=="Clean"),
                  quote(Species=="PF"&TrackingMethod=="Custom"),
                  "PF clean","PF custom",tp),
    perform_ttest(combined_data0,
                  quote(Species=="PY"&TrackingMethod=="Clean"),
                  quote(Species=="PY"&TrackingMethod=="Custom"),
                  "PY clean","PY custom",tp)
  )
}
ttest_df0 <- bind_rows(ttest_results0)

# Linear Regression per species & tracking method

species_list       <- c("PF", "PY")
trackingmethod_lst <- c("Clean", "Custom")
lm_results         <- list()

for(sp in species_list) {
  for(tm in trackingmethod_lst) {
    subdat <- combined_data %>% filter(Species==sp, TrackingMethod==tm)
    if(length(unique(subdat$Time))<2 || nrow(subdat)<2) {
      lm_results[[length(lm_results)+1]] <- tibble(
        Species=sp, TrackingMethod=tm,
        slope=NA, p_value=NA, r_value=NA, r_sq=NA,
        cohen_f2=NA, n=nrow(subdat)
      )
      next
    }
    model <- lm(log(Mean_Speed) ~ Time, data=subdat)
    sm    <- summary(model)
    r2    <- sm$r.squared
    f2    <- if(!is.na(r2)&&r2<1) r2/(1-r2) else NA
    lm_results[[length(lm_results)+1]] <- tibble(
      Species=sp,
      TrackingMethod=tm,
      slope=sm$coefficients[2,1],
      p_value=sm$coefficients[2,4],
      r_value=cor(log(subdat$Mean_Speed), subdat$Time),
      r_sq=r2,
      cohen_f2=f2,
      n=nrow(subdat)
    )
  }
}
lm_df <- bind_rows(lm_results)

species_list0       <- c("PF", "PY")
trackingmethod_lst0 <- c("Clean", "Custom")
lm_results0         <- list()

for(sp in species_list0) {
  for(tm in trackingmethod_lst0) {
    subdat <- combined_data0 %>% filter(Species==sp, TrackingMethod==tm)
    if(length(unique(subdat$Time))<2 || nrow(subdat)<2) {
      lm_results[[length(lm_results)+1]] <- tibble(
        Species=sp, TrackingMethod=tm,
        slope=NA, p_value=NA, r_value=NA, r_sq=NA,
        cohen_f2=NA, n=nrow(subdat)
      )
      next
    }
    model <- lm(log(Mean_Speed) ~ Time, data=subdat)
    sm    <- summary(model)
    r2    <- sm$r.squared
    f2    <- if(!is.na(r2)&&r2<1) r2/(1-r2) else NA
    lm_results0[[length(lm_results0)+1]] <- tibble(
      Species=sp,
      TrackingMethod=tm,
      slope=sm$coefficients[2,1],
      p_value=sm$coefficients[2,4],
      r_value=cor(log(subdat$Mean_Speed), subdat$Time),
      r_sq=r2,
      cohen_f2=f2,
      n=nrow(subdat)
    )
  }
}
lm_df0 <- bind_rows(lm_results0)
# Export all results to a single Excel file

wb_out <- createWorkbook()
addWorksheet(wb_out, "Pairwise_Tests")
addWorksheet(wb_out, "Linear_Models")
writeData(wb_out, "Pairwise_Tests", ttest_df)
writeData(wb_out, "Linear_Models", lm_df)
saveWorkbook(wb_out, "CStat_tests_mean_speeds.xlsx", overwrite=TRUE)
cat("\nAll T-tests & linear models saved to 'CStat_tests_mean_speeds.xlsx'.\n")


# Plot with dashed regression line and annotation

generate_plot <- function(data_summarized, species_label, track_id_type, annotation_label) {
  df <- data_summarized %>%
    filter(Species==species_label) %>%
    rename(Time=`Time Post-Inoculation (min)`) %>%
    filter(!is.na(Mean_Speed), Mean_Speed>0) %>%
    mutate(
      Time_str     = factor(paste0(Time," min"),
                            levels=c("5 min","10 min","20 min","30 min","60 min","120 min")),
      Time_numeric = as.numeric(Time_str)
    )
  
  group_stats <- df %>%
    group_by(Time_str, Time_numeric) %>%
    summarise(
      geom_mean = exp(mean(log(Mean_Speed), na.rm=TRUE)),
      Count     = n(),
      .groups   = "drop"
    )
  
  # Extended dashed‐line
  fit_means <- lm(log(geom_mean) ~ Time_numeric, data=group_stats)
  x0 <- min(group_stats$Time_numeric) - 0.5
  x1 <- max(group_stats$Time_numeric) + 0.6
  y0 <- exp(predict(fit_means, newdata=data.frame(Time_numeric=min(group_stats$Time_numeric))))
  y1 <- exp(predict(fit_means, newdata=data.frame(Time_numeric=max(group_stats$Time_numeric))))
  line_df <- tibble(x=x0, xend=x1, y=y0, yend=y1)
  
  # Compute stats for annotation
  raw_model <- lm(log(Mean_Speed) ~ Time, data=df)
  rm_sum    <- summary(raw_model)
  slope_raw <- rm_sum$coefficients[2,1]
  p_raw     <- rm_sum$coefficients[2,4]
  r2_raw    <- rm_sum$r.squared
  f2_raw    <- if(!is.na(r2_raw)&&r2_raw<1) r2_raw/(1-r2_raw) else NA
  
  slope_txt <- ifelse(abs(slope_raw)<0.001,
                      formatC(slope_raw,format="e",digits=1),
                      sprintf("%.3f",slope_raw))
  p_txt     <- ifelse(p_raw<0.001,
                      formatC(p_raw,format="e",digits=1),
                      sprintf("%.3f",p_raw))
  f2_txt    <- ifelse(f2_raw<0.001,
                      formatC(f2_raw,format="e",digits=1),
                      sprintf("%.3f",f2_raw))
  
  slope_num <- as.numeric(slope_txt)
  p_num     <- as.numeric(p_txt)
  f2_num    <- as.numeric(f2_txt)
  
  stats_expr <- substitute(
    "(         " ~ italic(m)==M ~ "," ~ italic(p)==P ~ ", "~
      "Cohen's " ~ italic(f)^2==F2 ~ ")",
    list(M = slope_num, P = p_num, F2 = f2_num)
  )
  
  color_map <- c(
    "5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
    "30 min"="darkgreen","60 min"="blue","120 min"="#E75480"
  )
  
  p <- ggplot(df, aes(x=Time_str, y=Mean_Speed, fill=Time_str)) +
    geom_violin(alpha=0.3, scale="width", width=0.9, color="black") +
    scale_x_discrete(labels = function(x) gsub(" min","",x)) +
    scale_y_log10(
      limits = c(0.007,300),
      breaks = c(0.01,0.1,1,10,100),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = c(0,0)
    ) +
    scale_fill_manual(values=color_map) +
    # mean segments
    geom_segment(
      data=group_stats,
      aes(x=Time_numeric-0.4, xend=Time_numeric+0.4,
          y=geom_mean, yend=geom_mean),
      inherit.aes=FALSE, color="black", linewidth=1.5
    ) +
    # extended dashed regression line
    geom_segment(
      data=line_df,
      aes(x=x, xend=xend, y=y, yend=yend),
      inherit.aes=FALSE, linetype="dashed", size=2, color="blue"
    ) +
    # µ and n labels
    annotate("text",
             x=group_stats$Time_numeric, y=40,
             label=paste0("µ = ", round(group_stats$geom_mean,2)," µm/s"),
             size=7.5, hjust=0.5) +
    annotate("text",
             x=group_stats$Time_numeric, y=20,
             label=paste0("n = ",group_stats$Count),
             size=7.5, hjust=0.5) +
    # species name
    annotate("text",
             x=0.5, y=200,
             label=annotation_label,
             size=14, hjust=0, fontface="italic") +
    # blue dashes
    annotate("segment", x=0.56, xend=0.73, y=80, yend=80,
             colour="blue", linewidth=1.6, linetype="longdash")+
    # stats expression
    annotate("text",
             x=0.58, y=80,
             label=stats_expr,
             parse=FALSE, size=10, hjust=0) +
    labs(x="Time Post-Inoculation (min)", y="Mean Speed (µm/s)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour="black", size=1.25),
      axis.text.x      = element_text(size=26, color="black"),
      axis.text.y      = element_text(size=26, color="black"),
      axis.title       = element_text(size=26),
      axis.ticks.x.bottom=element_line(color="black", size = 1.25),
      axis.ticks.length=unit(7,"pt"),
      axis.ticks.y.left=element_line(color="black", size = 1.25),
      legend.position  = "none"
    )
  
  return(p)
}

# 7. MAKE COMBINED PLOTS (PF vs PY), FOR BOTH CUSTOM & CLEAN
pf_c  <- generate_plot(data_trackid_ms,    "PF",    "Custom", expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_c  <- generate_plot(data_trackid_ms,    "PY",    "Custom",  expression("P. yoelii (f" [threshold] * "= 5µm)"))
pf_c0  <- generate_plot(data_trackid_ms0,    "PF",    "Custom", "P. falciparum")
py_c0  <- generate_plot(data_trackid_ms0,    "PY",    "Custom", "P. yoelii")

panel <- plot_grid(pf_c0, py_c0, pf_c, py_c,ncol=2, nrow=2, labels=c("A","B", "C","D"), label_size=40)
ggsave("plot_CustomTrackID_Pf_Py_side_by_side_mean_speed.pdf", panel, width=28, height=20)

pf_cl0 <- generate_plot(data_newtrackid_ms0, "PF",    "Clean",  "P. falciparum")
py_cl0 <- generate_plot(data_newtrackid_ms0, "PY",    "Clean",  "P. yoelii")
pf_cl <- generate_plot(data_newtrackid_ms,"PF","Clean",expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_cl <- generate_plot(data_newtrackid_ms, "PY",    "Clean",  expression("P. yoelii (f" [threshold] * "= 5µm)"))
panel2 <- plot_grid(pf_cl0, py_cl0, pf_cl, py_cl,ncol=2, nrow=2,labels=c("A","B", "C","D"), label_size=40)
ggsave("plot_CleanTrackIDs_Pf_Py_side_by_side_mean_speed.pdf", panel2, width=28, height=20)

cat("\nAll done with mean‐speed analysis!\n")









################################################################
################################################################
################################################################











# Load Libraries

#library(readxl)
library(dplyr)
library(MASS)
#library(ggplot2)
#library(cowplot)
#library(openxlsx)
#library(effsize)
#library(purrr)
#library(tidyr)
#library(ggtext)


# Read Files

#data_original <- read_excel("Hopp_PF-PY_Master_file.xlsx")


# Calculate turning angles

compute_turning_angles <- function(data, track_col) {
  data %>%
    dplyr::group_by(Species, !!rlang::sym(track_col)) %>%
    dplyr::arrange(`Time (s)`, .by_group = TRUE) %>%
    dplyr::mutate(
      Time_Diff = `Time (s)` - dplyr::lag(`Time (s)`),
      dX        = `Position X (um)` - dplyr::lag(`Position X (um)`),
      dY        = `Position Y (um)` - dplyr::lag(`Position Y (um)`)
    ) %>%
    dplyr::filter(!is.na(Time_Diff), Time_Diff > 0) %>%
    dplyr::mutate(
      Angle         = atan2(dY, dX) * (180 / pi),
      Turning_Angle = Angle - dplyr::lag(Angle),
      Turning_Angle = dplyr::if_else(
        Turning_Angle > 180,  Turning_Angle - 360,
        dplyr::if_else(Turning_Angle < -180, Turning_Angle + 360, Turning_Angle)
      ),
      Turning_Angle = abs(Turning_Angle)
    ) %>%
    dplyr::filter(!is.na(Turning_Angle)) %>%
    dplyr::select(-dplyr::any_of(c("Time_Diff","dX","dY","Angle"))) %>%
    dplyr::ungroup()
}


data_custom_angles <- compute_turning_angles(data_original, "Custom TrackID")
data_clean_angles  <- compute_turning_angles(data_original, "Clean TrackIDs")
data_custom_angles_ta <- compute_turning_angles(data_custom_raw, "Custom TrackID")
data_clean_angles_ta  <- compute_turning_angles(data_clean_raw, "Clean TrackIDs")


# Combine Data

combined_data <- bind_rows(
  data_custom_angles %>% mutate(TrackingMethod = "Custom"),
  data_clean_angles  %>% mutate(TrackingMethod = "Clean")
) %>%
  rename(Time = `Time Post-Inoculation (min)`)
combined_data_ta <- bind_rows(
  data_custom_angles_ta %>% mutate(TrackingMethod = "Custom"),
  data_clean_angles_ta  %>% mutate(TrackingMethod = "Clean")
) %>%
  rename(Time = `Time Post-Inoculation (min)`)


# T-tests and Cohen's d

all_timepoints <- c(5,10,20,30,60,120)
test_angle <- function(df, f1, f2, lbl1, lbl2, tp) {
  g1 <- df %>% filter(!!f1, Time==tp) %>% pull(Turning_Angle)
  g2 <- df %>% filter(!!f2, Time==tp) %>% pull(Turning_Angle)
  if(length(g1)<2 || length(g2)<2) return(
    tibble(Time=tp, Comparison=paste(lbl1,"vs",lbl2),
           N1=length(g1), N2=length(g2),
           t_stat=NA, p_value=NA, cohen_d=NA)
  )
  tt <- t.test(g1, g2)
  d  <- effsize::cohen.d(g1, g2)$estimate
  tibble(Time=tp, Comparison=paste(lbl1,"vs",lbl2),
         N1=length(g1), N2=length(g2),
         t_stat=tt$statistic[[1]], p_value=tt$p.value, cohen_d=d)
}

comparisons <- list(
  list(quote(Species=="PF" & TrackingMethod=="Clean"),
       quote(Species=="PY" & TrackingMethod=="Clean"), "PF clean", "PY clean"),
  list(quote(Species=="PF" & TrackingMethod=="Custom"),
       quote(Species=="PY" & TrackingMethod=="Custom"),"PF custom","PY custom"),
  list(quote(Species=="PF" & TrackingMethod=="Clean"),
       quote(Species=="PF" & TrackingMethod=="Custom"),"PF clean","PF custom"),
  list(quote(Species=="PY" & TrackingMethod=="Clean"),
       quote(Species=="PY" & TrackingMethod=="Custom"),"PY clean","PY custom")
)

ttest_df <- map_dfr(all_timepoints, function(tp) {
  map_dfr(comparisons, ~test_angle(combined_data,
                                   .x[[1]], .x[[2]],
                                   .x[[3]], .x[[4]], tp))
})
ttest_df_ta <- map_dfr(all_timepoints, function(tp) {
  map_dfr(comparisons, ~test_angle(combined_data_ta,
                                   .x[[1]], .x[[2]],
                                   .x[[3]], .x[[4]], tp))
})


# Linear Regression

grid <- expand_grid(sp = c("PF","PY"), method = c("Clean","Custom"))
lm_df <- pmap_dfr(grid, function(sp, method) {
  sub <- combined_data %>% filter(Species==sp, TrackingMethod==method)
  grp <- sub %>%
    group_by(Time) %>%
    summarise(mean_angle = mean(Turning_Angle), .groups = "drop")
  if(nrow(grp) < 2) {
    return(tibble(Species=sp, TrackingMethod=method,
                  slope=NA, p_value=NA, r_value=NA,
                  r_sq=NA, cohen_f2=NA, N=nrow(sub)))
  }
  mod <- lm(mean_angle ~ Time, data=grp)
  sm  <- summary(mod)
  r2  <- sm$r.squared
  f2  <- if(r2 < 1) r2/(1-r2) else NA
  tibble(Species=sp, TrackingMethod=method,
         slope=sm$coefficients[2,1],
         p_value=sm$coefficients[2,4],
         r_value=cor(grp$mean_angle, grp$Time),
         r_sq=r2, cohen_f2=f2, N=nrow(sub))
})

lm_df_ta <- pmap_dfr(grid, function(sp, method) {
  sub <- combined_data_ta %>% filter(Species==sp, TrackingMethod==method)
  grp <- sub %>%
    group_by(Time) %>%
    summarise(mean_angle = mean(Turning_Angle), .groups = "drop")
  if(nrow(grp) < 2) {
    return(tibble(Species=sp, TrackingMethod=method,
                  slope=NA, p_value=NA, r_value=NA,
                  r_sq=NA, cohen_f2=NA, N=nrow(sub)))
  }
  mod <- lm(mean_angle ~ Time, data=grp)
  sm  <- summary(mod)
  r2  <- sm$r.squared
  f2  <- if(r2 < 1) r2/(1-r2) else NA
  tibble(Species=sp, TrackingMethod=method,
         slope=sm$coefficients[2,1],
         p_value=sm$coefficients[2,4],
         r_value=cor(grp$mean_angle, grp$Time),
         r_sq=r2, cohen_f2=f2, N=nrow(sub))
})
# Export stats to excel files

wb <- createWorkbook()
addWorksheet(wb, "TTests")
addWorksheet(wb, "LinModels")
writeData(wb, "TTests", ttest_df)
writeData(wb, "LinModels", lm_df)
saveWorkbook(wb, "Stat_tests_turning_angles_Abs.xlsx", overwrite = TRUE)

# Plot turning angles

plot_turning0<- function(df, species, method, label) {
  sub <- df %>%
    filter(Species==species, TrackingMethod==method) %>%
    mutate(
      Time_str = factor(paste0(Time, " min"),
                        levels=c("5 min","10 min","20 min","30 min","60 min","120 min")),
      Time_num = as.numeric(Time_str)
    )
  
  stats <- sub %>%
    group_by(Time_str, Time_num) %>%
    summarise(
      mean_angle = mean(Turning_Angle),
      Count      = n(),
      .groups    = "drop"
    )
  
  # Extended dashed‐line
  
  fit_means <- lm(mean_angle ~ Time_num, data=stats)
  x0        <- min(stats$Time_num) - 0.5
  x1        <- max(stats$Time_num) + 0.6
  y0        <- predict(fit_means, newdata = data.frame(Time_num=min(stats$Time_num)))
  y1        <- predict(fit_means, newdata = data.frame(Time_num=max(stats$Time_num)))
  line_df   <- tibble(x=x0, xend=x1, y=y0, yend=y1)
  
  match_row <- lm_df %>% filter(Species==species, TrackingMethod==method)
  slope_reg <- match_row$slope
  p_reg     <- match_row$p_value
  f2_reg    <- match_row$cohen_f2
  
  slope_txt <- ifelse(abs(slope_reg)<0.001, formatC(slope_reg, format="e", digits=1), sprintf("%.3f", slope_reg))
  p_txt     <- ifelse(p_reg<0.001,         formatC(p_reg,     format="e", digits=1), sprintf("%.3f", p_reg))
  f2_txt    <- ifelse(f2_reg<0.001,        formatC(f2_reg,    format="e", digits=1), sprintf("%.3f", f2_reg))
  
  
  
  slope_num <- as.numeric(slope_txt)
  p_num     <- as.numeric(p_txt)
  f2_num    <- as.numeric(f2_txt)
  
  label_expr <- substitute(
    "(           " ~ italic(m)==M ~ "," ~ italic(p)==P ~ ", Cohen’s " ~ italic(f)^2==F2 ~ ")",
    list(M = slope_num, P = p_num, F2 = f2_num)
  )
  
  col_map   <- c("5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
                 "30 min"="darkgreen","60 min"="blue","120 min"="#E75480")
  shape_map <- c("5 min"=21,"10 min"=3,"20 min"=23,"30 min"=24,"60 min"=22,"120 min"=4)
  
  ggplot(sub, aes(x=Time_str, y=Turning_Angle, color=Time_str, shape=Time_str)) +
    geom_violin(aes(fill=Time_str), alpha=.3, width=.9, color="black") +
    geom_segment(
      data=line_df,
      aes(x=x, xend=xend, y=y, yend=yend),
      inherit.aes=FALSE, linetype="dashed", color="blue", size=2.5
    ) +
    scale_x_discrete(labels=function(x) gsub(" min","",x)) +
    scale_y_continuous(limits=c(0,280), breaks=c(0,45,90,135,180)) +
    scale_color_manual(values=col_map) +
    scale_fill_manual(values=col_map) +
    scale_shape_manual(values=shape_map) +
    labs(x="Time Post-Inoculation (min)", y="Turning Angle (°)") +
    theme_minimal() +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(colour="black", size=1.25),
      axis.text         = element_text(size=26, color="black"),
      axis.ticks.length = unit(7,"pt"),
      axis.ticks        = element_line(size=1.25),
      axis.title        = element_text(size=26),
      legend.position   = "none"
    ) +
    geom_segment(
      data=stats,
      aes(x=Time_num-0.4, xend=Time_num+0.4, y=mean_angle, yend=mean_angle),
      inherit.aes=FALSE, size=1.5, color="black"
    ) +
    
    annotate(
      "text",
      x     = stats$Time_num,
      y     = 230,
      label = paste0("µ = ", round(stats$mean_angle,2),"°"),
      size  = 9
    ) +
    annotate(
      "text",
      x     = stats$Time_num,
      y     = 210,
      label = paste0("n = ", stats$Count),
      size  = 9
    ) +
    annotate(
      "text",
      x     = 0.5,
      y     = 275,
      label = label,
      size  = 14,
      hjust = 0,
      fontface="italic"
    ) +
    
    annotate("text",
             x=0.5+0.08, y=255,
             label="  — — —",
             color="blue",
             fontface = "bold",
             size=11, hjust=0) +
    annotate(
      "text",
      x     = 0.5 + 0.07,
      y     = 255,
      label = label_expr,
      size  = 10,
      hjust = 0
    )
}

plot_turning<- function(df, species, method, label) {
  sub <- df %>%
    filter(Species==species, TrackingMethod==method) %>%
    mutate(
      Time_str = factor(paste0(Time, " min"),
                        levels=c("5 min","10 min","20 min","30 min","60 min","120 min")),
      Time_num = as.numeric(Time_str)
    )
  
  stats <- sub %>%
    group_by(Time_str, Time_num) %>%
    summarise(
      mean_angle = mean(Turning_Angle),
      Count      = n(),
      .groups    = "drop"
    )
  
  # Extended dashed‐line
  
  fit_means <- lm(mean_angle ~ Time_num, data=stats)
  x0        <- min(stats$Time_num) - 0.5
  x1        <- max(stats$Time_num) + 0.6
  y0        <- predict(fit_means, newdata = data.frame(Time_num=min(stats$Time_num)))
  y1        <- predict(fit_means, newdata = data.frame(Time_num=max(stats$Time_num)))
  line_df   <- tibble(x=x0, xend=x1, y=y0, yend=y1)
  
  match_row <- lm_df_ta %>% filter(Species==species, TrackingMethod==method)
  slope_reg <- match_row$slope
  p_reg     <- match_row$p_value
  f2_reg    <- match_row$cohen_f2
  
  slope_txt <- ifelse(abs(slope_reg)<0.001, formatC(slope_reg, format="e", digits=1), sprintf("%.3f", slope_reg))
  p_txt     <- ifelse(p_reg<0.001,         formatC(p_reg,     format="e", digits=1), sprintf("%.3f", p_reg))
  f2_txt    <- ifelse(f2_reg<0.001,        formatC(f2_reg,    format="e", digits=1), sprintf("%.3f", f2_reg))
  
  
  
  slope_num <- as.numeric(slope_txt)
  p_num     <- as.numeric(p_txt)
  f2_num    <- as.numeric(f2_txt)
  
  label_expr <- substitute(
    "(           " ~ italic(m)==M ~ "," ~ italic(p)==P ~ ", Cohen’s " ~ italic(f)^2==F2 ~ ")",
    list(M = slope_num, P = p_num, F2 = f2_num)
  )
  
  col_map   <- c("5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
                 "30 min"="darkgreen","60 min"="blue","120 min"="#E75480")
  shape_map <- c("5 min"=21,"10 min"=3,"20 min"=23,"30 min"=24,"60 min"=22,"120 min"=4)
  
  ggplot(sub, aes(x=Time_str, y=Turning_Angle, color=Time_str, shape=Time_str)) +
    geom_violin(aes(fill=Time_str), alpha=.3, width=.9, color="black") +
    geom_segment(
      data=line_df,
      aes(x=x, xend=xend, y=y, yend=yend),
      inherit.aes=FALSE, linetype="dashed", color="blue", size=2.5
    ) +
    scale_x_discrete(labels=function(x) gsub(" min","",x)) +
    scale_y_continuous(limits=c(0,280), breaks=c(0,45,90,135,180)) +
    scale_color_manual(values=col_map) +
    scale_fill_manual(values=col_map) +
    scale_shape_manual(values=shape_map) +
    labs(x="Time Post-Inoculation (min)", y="Turning Angle (°)") +
    theme_minimal() +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(colour="black", size=1.25),
      axis.text         = element_text(size=26, color="black"),
      axis.ticks.length = unit(7,"pt"),
      axis.ticks        = element_line(size=1.25),
      axis.title        = element_text(size=26),
      legend.position   = "none"
    ) +
    geom_segment(
      data=stats,
      aes(x=Time_num-0.4, xend=Time_num+0.4, y=mean_angle, yend=mean_angle),
      inherit.aes=FALSE, size=1.5, color="black"
    ) +
    
    annotate(
      "text",
      x     = stats$Time_num,
      y     = 230,
      label = paste0("µ = ", round(stats$mean_angle,2),"°"),
      size  = 9
    ) +
    annotate(
      "text",
      x     = stats$Time_num,
      y     = 210,
      label = paste0("n = ", stats$Count),
      size  = 9
    ) +
    annotate(
      "text",
      x     = 0.5,
      y     = 275,
      label = label,
      size  = 14,
      hjust = 0,
      fontface="italic"
    ) +
    
    annotate("text",
             x=0.5+0.08, y=255,
             label="  — — —",
             color="blue",
             fontface = "bold",
             size=11, hjust=0) +
    annotate(
      "text",
      x     = 0.5 + 0.07,
      y     = 255,
      label = label_expr,
      size  = 10,
      hjust = 0
    )
}

# Generate & Save
pf_custom_plot0 <- plot_turning0(combined_data, "PF", "Custom", "P. falciparum")
py_custom_plot0 <- plot_turning0(combined_data, "PY", "Custom", "P. yoelii")
pf_custom_plot <- plot_turning(combined_data_ta, "PF", "Custom", expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_custom_plot <- plot_turning(combined_data_ta, "PY", "Custom", expression("P. yoelii (f" [threshold] * "= 5µm)"))
ggsave("TurningAngle_Custom_Pf_Py_Abs.pdf",
       plot_grid(pf_custom_plot0, py_custom_plot0, pf_custom_plot, py_custom_plot,ncol=2, nrow=2,labels=c("A","B","C","D"), label_size=40),
       width=28, height=20) 

pf_clean_plot0 <- plot_turning0(combined_data, "PF", "Clean", "P. falciparum")
py_clean_plot0 <- plot_turning0(combined_data, "PY", "Clean", "P. yoelii")
pf_clean_plot <- plot_turning(combined_data_ta, "PF", "Clean", expression("P. falciparum (f" [threshold] * "= 5µm)"))
py_clean_plot <- plot_turning(combined_data_ta, "PY", "Clean", expression("P. yoelii (f" [threshold] * "= 5µm)"))
ggsave("TurningAngle_Clean_Pf_Py_Abs.pdf",
       plot_grid(pf_clean_plot0, py_clean_plot0, pf_clean_plot, py_clean_plot,ncol=2, nrow=2, labels=c("A","B","C","D"), label_size=40),
       width=28, height=20)

cat("Pipeline complete!\n")




################################################################
################################################################
################################################################





# Load libraries

library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(car)
library(scales)
library(writexl)


# Read files

#data_original <- read_excel("Hopp_PF-PY_Master_file.xlsx")

data_clean_raw  <- filter_tracks_by_spread(
  data_original, track_col = "Clean TrackIDs",
  threshold_um = threshold_um
)


# Pre-processing
data_prepped0 <- data_original %>%
  rename(
    Time = `Time (s)`,
    X = `Position X (um)`,
    Y = `Position Y (um)`,
    TrackID = `Clean TrackIDs`,
    TimeGroup = `Time Post-Inoculation (min)`
  ) %>%
  filter(!is.na(X), !is.na(Y), !is.na(Time)) %>%
  group_by(Species, TimeGroup, TrackID) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    Time_Adjusted = Time - min(Time),
    Displacement_Squared = (X - first(X))^2 + (Y - first(Y))^2
  ) %>%
  
  filter(all(diff(Time_Adjusted) > 0)) %>%
  ungroup()

data_prepped <- data_clean_raw %>%
  rename(
    Time = `Time (s)`,
    X = `Position X (um)`,
    Y = `Position Y (um)`,
    TrackID = `Clean TrackIDs`,
    TimeGroup = `Time Post-Inoculation (min)`
  ) %>%
  filter(!is.na(X), !is.na(Y), !is.na(Time)) %>%
  group_by(Species, TimeGroup, TrackID) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(
    Time_Adjusted = Time - min(Time),
    Displacement_Squared = (X - first(X))^2 + (Y - first(Y))^2
  ) %>%
  
  filter(all(diff(Time_Adjusted) > 0)) %>%
  ungroup()


# MSD Curves

calculate_msd_curves <- function(df, species_filter) {
  df %>%
    filter(Species == species_filter, Time_Adjusted > 0) %>%
    group_by(TimeGroup, TrackID, Time_Adjusted) %>%
    summarise(Squared_Disp = mean(Displacement_Squared), .groups = "drop") %>%
    group_by(TimeGroup, Time_Adjusted) %>%
    summarise(MSD = mean(Squared_Disp), .groups = "drop") %>%
    mutate(TimeGroup = factor(
      TimeGroup,
      levels = c(5, 10, 20, 30, 60, 120),
      labels = c("5 min", "10 min", "20 min", "30 min", "60 min", "120 min")
    ))
}


# Meta Info

calculate_meta_info <- function(df, species_filter) {
  df %>%
    filter(Species == species_filter, Time_Adjusted > 0) %>%
    group_by(TimeGroup) %>%
    summarise(
      n_tracks = n_distinct(TrackID),
      mean_msd = mean(Displacement_Squared, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(TimeGroup = factor(
      TimeGroup,
      levels = c(5, 10, 20, 30, 60, 120),
      labels = c("5 min", "10 min", "20 min", "30 min", "60 min", "120 min")
    ))
}

library(dplyr)
library(tidyr)
library(writexl)   # if you want the Excel export

# --- Helper to build EAMSD table once (Species × TimeGroup × Time_Adjusted) ---
eammsd_tbl <- function(df) {
  df %>%
    dplyr::filter(Time_Adjusted > 0) %>%
    dplyr::group_by(Species, TimeGroup, Time_Adjusted) %>%
    dplyr::summarise(MSD = mean(Displacement_Squared, na.rm = TRUE), .groups = "drop")
}

# --- Mean MSD slope per species (unweighted and weighted by #tracks) ---
compute_mean_msd_slope0 <- function(data_prepped0, tmax = 30) {
  # Ensemble-averaged MSD surface
  eamsd <- eammsd_tbl(data_prepped0)
  
  # Slopes α (log10 MSD ~ log10 Time) per Species × TimeGroup in the window (0, tmax]
  slopes_by_group <- eamsd %>%
    dplyr::filter(Time_Adjusted <= tmax, MSD > 0, is.finite(MSD)) %>%
    dplyr::group_by(Species, TimeGroup) %>%
    dplyr::group_modify(~{
      d <- .x
      if (nrow(d) < 3) return(tibble::tibble(Slope = NA_real_))
      fit <- stats::lm(log10(MSD) ~ log10(Time_Adjusted), data = d)
      tibble::tibble(Slope = unname(stats::coef(fit)[2]))
    }) %>%
    dplyr::ungroup()
  
  # Track counts per Species × TimeGroup (for weighted mean)
  track_counts <- data_prepped0 %>%
    dplyr::filter(Time_Adjusted > 0) %>%
    dplyr::group_by(Species, TimeGroup) %>%
    dplyr::summarise(n_tracks = dplyr::n_distinct(TrackID), .groups = "drop")
  
  slopes_joined <- slopes_by_group %>%
    dplyr::left_join(track_counts, by = c("Species", "TimeGroup"))
  
  # Mean slopes per Species
  mean_slopes <- slopes_joined %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      mean_slope_unweighted = mean(Slope, na.rm = TRUE),
      mean_slope_weighted_tracks = stats::weighted.mean(Slope, w = n_tracks, na.rm = TRUE),
      groups_used = sum(!is.na(Slope)),
      .groups = "drop"
    )
  
  list(
    slopes_per_group0 = slopes_joined,  # one slope per TimeGroup
    mean_slopes0      = mean_slopes     # per-species averages
  )
}

compute_mean_msd_slope <- function(data_prepped, tmax = 30) {
  # Ensemble-averaged MSD surface
  eamsd <- eammsd_tbl(data_prepped)
  
  # Slopes α (log10 MSD ~ log10 Time) per Species × TimeGroup in the window (0, tmax]
  slopes_by_group <- eamsd %>%
    dplyr::filter(Time_Adjusted <= tmax, MSD > 0, is.finite(MSD)) %>%
    dplyr::group_by(Species, TimeGroup) %>%
    dplyr::group_modify(~{
      d <- .x
      if (nrow(d) < 3) return(tibble::tibble(Slope = NA_real_))
      fit <- stats::lm(log10(MSD) ~ log10(Time_Adjusted), data = d)
      tibble::tibble(Slope = unname(stats::coef(fit)[2]))
    }) %>%
    dplyr::ungroup()
  
  # Track counts per Species × TimeGroup (for weighted mean)
  track_counts <- data_prepped %>%
    dplyr::filter(Time_Adjusted > 0) %>%
    dplyr::group_by(Species, TimeGroup) %>%
    dplyr::summarise(n_tracks = dplyr::n_distinct(TrackID), .groups = "drop")
  
  slopes_joined <- slopes_by_group %>%
    dplyr::left_join(track_counts, by = c("Species", "TimeGroup"))
  
  # Mean slopes per Species
  mean_slopes <- slopes_joined %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      mean_slope_unweighted = mean(Slope, na.rm = TRUE),
      mean_slope_weighted_tracks = stats::weighted.mean(Slope, w = n_tracks, na.rm = TRUE),
      groups_used = sum(!is.na(Slope)),
      .groups = "drop"
    )
  
  list(
    slopes_per_group = slopes_joined,  # one slope per TimeGroup
    mean_slopes      = mean_slopes     # per-species averages
  )
}

# --- Run it (uses your existing data_prepped) ---
msd_slope_res <- compute_mean_msd_slope(data_prepped, tmax = 30)
msd_slope_res0 <- compute_mean_msd_slope(data_prepped0, tmax = 30)

# View in console
msd_slope_res$mean_slopes
# A tibble with columns: Species, mean_slope_unweighted, mean_slope_weighted_tracks, groups_used

# (Optional) write to Excel
write_xlsx(
  list(
    "MSD_slopes_per_group" = msd_slope_res$slopes_per_group,
    "MSD_mean_slopes"      = msd_slope_res$mean_slopes
  ),
  "MSD_mean_slopes_by_species.xlsx"
)

# Plot 
generate_msd_plot_obj0 <- function(species_name, annotation_label) {
  msd_data <- calculate_msd_curves(data_prepped0, species_name)
  meta_info <- calculate_meta_info(data_prepped0, species_name)
  
  
  slope_info <- msd_data %>%
    filter(Time_Adjusted <= 30) %>%
    group_by(TimeGroup) %>%
    do({
      model <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
      data.frame(Slope = round(coef(model)[2], 2))
    })
  
  meta_merged <- left_join(meta_info, slope_info, by = "TimeGroup")
  
  
  legend_labels <- sapply(levels(meta_merged$TimeGroup), function(g) {
    row <- meta_merged %>% filter(TimeGroup == g)
    paste0(g, ": µ = ", round(row$mean_msd, 2), " µm²",
           ", n = ", row$n_tracks, ", m = ", row$Slope)
  })
  names(legend_labels) <- levels(meta_merged$TimeGroup)
  
  regression_lines <- msd_data %>%
    filter(Time_Adjusted <= 30) %>%
    group_by(TimeGroup) %>%
    do({
      model <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
      Time_vals <- seq(1, 30, length.out = 100)
      preds <- predict(model, newdata = data.frame(Time_Adjusted = Time_vals))
      data.frame(
        Time_Adjusted = Time_vals,
        MSD = 10^preds,
        TimeGroup = unique(.$TimeGroup)
      )
    })
  
  custom_colors <- c("5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
                     "30 min"="darkgreen","60 min"="blue","120 min"="#E75480")
  custom_shapes <- c("5 min"=21,"10 min"=3,"20 min"=23,
                     "30 min"=24,"60 min"=22,"120 min"=4)
  
  p_main <- ggplot(msd_data, aes(x = Time_Adjusted, y = MSD,
                                 color = TimeGroup, shape = TimeGroup)) +
    geom_line(size = 0.75, alpha = ifelse(msd_data$Time_Adjusted > 30, 0.2, 1)) +
    geom_point(size = 2, fill = NA, stroke = 2, alpha = ifelse(msd_data$Time_Adjusted > 30, 0.2, 1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = 4.5, y = 7000, label = annotation_label, hjust = 0.5,
             vjust = 1, size = 12, fontface = "italic") +
    scale_color_manual(
      name = "Time Post-Inoculation",
      values = custom_colors,
      labels = legend_labels
    ) +
    scale_shape_manual(
      name = "Time Post-Inoculation",
      values = custom_shapes,
      labels = legend_labels
    ) +
    scale_y_log10(limits = c(1, 7000), breaks = c(1, 10, 100, 1000)) +
    scale_x_log10(limits = c(1, 7000), breaks = c(1, 10, 100, 1000)) +
    coord_fixed(ratio = 1) +
    labs(x = "Time (s)", y = expression(MSD ~ (µm^2))) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18),
      axis.ticks = element_line(color = "black", size = 1),
      axis.ticks.length = unit(7, "pt"),
      legend.position = c(0.27, 0.78),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18, face = "bold")
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = "blank")),
      shape = guide_legend(override.aes = list(linetype = "blank"))
    )
  
  # Inset plot
  p_inset <- ggplot(regression_lines, aes(x = Time_Adjusted, y = MSD, color = TimeGroup)) +
    geom_line(size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = custom_colors) +
    scale_y_log10(limits = c(1, 200), breaks = c(1, 10, 100)) +
    scale_x_log10(limits = c(1, 200), breaks = c(1, 10, 100)) +
    coord_fixed(ratio = 1) +
    labs(x = "Time (s)", y = expression(MSD ~ (µm^2))) +
    annotate("text", x = 4, y = 100, label = expression(t <= 30~s),
             size = 6, fontface = "bold") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(size = 14),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 1),
      axis.ticks.length = unit(7, "pt"),
      legend.position = "none"
    )
  
  combined <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset, x = 0.66, y = 0.085, width = 0.3, height = 0.3)
  
  return(combined)
}

generate_msd_plot_obj <- function(species_name, annotation_label) {
  msd_data <- calculate_msd_curves(data_prepped, species_name)
  meta_info <- calculate_meta_info(data_prepped, species_name)
  
  
  slope_info <- msd_data %>%
    filter(Time_Adjusted <= 30) %>%
    group_by(TimeGroup) %>%
    do({
      model <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
      data.frame(Slope = round(coef(model)[2], 2))
    })
  
  meta_merged <- left_join(meta_info, slope_info, by = "TimeGroup")
  
  
  legend_labels <- sapply(levels(meta_merged$TimeGroup), function(g) {
    row <- meta_merged %>% filter(TimeGroup == g)
    paste0(g, ": µ = ", round(row$mean_msd, 2), " µm²",
           ", n = ", row$n_tracks, ", m = ", row$Slope)
  })
  names(legend_labels) <- levels(meta_merged$TimeGroup)
  
  regression_lines <- msd_data %>%
    filter(Time_Adjusted <= 30) %>%
    group_by(TimeGroup) %>%
    do({
      model <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
      Time_vals <- seq(1, 30, length.out = 100)
      preds <- predict(model, newdata = data.frame(Time_Adjusted = Time_vals))
      data.frame(
        Time_Adjusted = Time_vals,
        MSD = 10^preds,
        TimeGroup = unique(.$TimeGroup)
      )
    })
  
  custom_colors <- c("5 min"="#008ECC","10 min"="#B63400","20 min"="darkmagenta",
                     "30 min"="darkgreen","60 min"="blue","120 min"="#E75480")
  custom_shapes <- c("5 min"=21,"10 min"=3,"20 min"=23,
                     "30 min"=24,"60 min"=22,"120 min"=4)
  
  p_main <- ggplot(msd_data, aes(x = Time_Adjusted, y = MSD,
                                 color = TimeGroup, shape = TimeGroup)) +
    geom_line(size = 0.75, alpha = ifelse(msd_data$Time_Adjusted > 30, 0.2, 1)) +
    geom_point(size = 2, fill = NA, stroke = 2, alpha = ifelse(msd_data$Time_Adjusted > 30, 0.2, 1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = 4.5, y = 7000, label = annotation_label, hjust = 0.5,
             vjust = 1, size = 12, fontface = "italic") +
    scale_color_manual(
      name = "Time Post-Inoculation",
      values = custom_colors,
      labels = legend_labels
    ) +
    scale_shape_manual(
      name = "Time Post-Inoculation",
      values = custom_shapes,
      labels = legend_labels
    ) +
    scale_y_log10(limits = c(1, 7000), breaks = c(1, 10, 100, 1000)) +
    scale_x_log10(limits = c(1, 7000), breaks = c(1, 10, 100, 1000)) +
    coord_fixed(ratio = 1) +
    labs(x = "Time (s)", y = expression(MSD ~ (µm^2))) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", size = 1),
      axis.text = element_text(size = 18, color = "black"),
      axis.title = element_text(size = 18),
      axis.ticks = element_line(color = "black", size = 1),
      axis.ticks.length = unit(7, "pt"),
      legend.position = c(0.27, 0.78),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 18, face = "bold")
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = "blank")),
      shape = guide_legend(override.aes = list(linetype = "blank"))
    )
  
  # Inset plot
  p_inset <- ggplot(regression_lines, aes(x = Time_Adjusted, y = MSD, color = TimeGroup)) +
    geom_line(size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = custom_colors) +
    scale_y_log10(limits = c(1, 200), breaks = c(1, 10, 100)) +
    scale_x_log10(limits = c(1, 200), breaks = c(1, 10, 100)) +
    coord_fixed(ratio = 1) +
    labs(x = "Time (s)", y = expression(MSD ~ (µm^2))) +
    annotate("text", x = 4, y = 100, label = expression(t <= 30~s),
             size = 6, fontface = "bold") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(size = 14),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks = element_line(color = "black", size = 1),
      axis.ticks.length = unit(7, "pt"),
      legend.position = "none"
    )
  
  combined <- ggdraw() +
    draw_plot(p_main) +
    draw_plot(p_inset, x = 0.66, y = 0.085, width = 0.3, height = 0.3)
  
  return(combined)
}



pf_plot0 <- generate_msd_plot_obj0("PF", "P. falciparum")
py_plot0 <- generate_msd_plot_obj0("PY", "P. yoelii")
pf_plot <- generate_msd_plot_obj("PF",  expression("                      P. falciparum (f" [threshold] * "= 5µm)"))
py_plot <- generate_msd_plot_obj("PY",  expression("                      P. yoelii (f" [threshold] * "= 5µm)"))

# Combine into one panel
combined_panel <- plot_grid(
  pf_plot0,
  py_plot0,
  pf_plot,
  py_plot,
  ncol = 2,
  nrow=2,
  labels = c("A", "B","C","D"),
  label_size = 48
)

ggsave("MSD_PF_PY_combined_panel.pdf", combined_panel, width = 28, height = 20, dpi = 300)
cat("✅ Saved combined panel to: MSD_PF_PY_combined_panel.pdf\n")


################################################################################
## =========================================================
## ==============  QUICK CHECKS: UTILITIES  ================
## =========================================================

library(purrr)
library(tidyr)
library(broom)
library(writexl)
library(ggplot2)
library(cowplot)

# Ensure factors for TimeGroup
lvl_ords <- c(5,10,20,30,60,120)
lbl_ords <- c("5 min","10 min","20 min","30 min","60 min","120 min")

# Helper: per-group (Species x TimeGroup) EAMSD table (already implicit in your plot code)
eammsd_tbl <- function(df) {
  df %>%
    filter(Time_Adjusted > 0) %>%
    group_by(Species, TimeGroup, Time_Adjusted) %>%
    summarise(MSD = mean(Displacement_Squared, na.rm=TRUE), .groups="drop") %>%
    mutate(TimeGroup = factor(
      TimeGroup, levels = lvl_ords, labels = lbl_ords
    ))
}

# Helper: compute per-track time-averaged MSD (TAMSD) over a set of lags
tamsd_one_track <- function(tt, xx, yy, lags) {
  out <- lapply(lags, function(tau){
    # pairs (i,j) with t[j]-t[i] >= tau; choose nearest within tolerance
    # Efficient approach: for each i, find j via binary search
    n <- length(tt)
    vals <- numeric(0)
    for (i in 1:(n-1)) {
      target <- tt[i] + tau
      j <- which.min(abs(tt[(i+1):n] - target))
      j <- j + i
      if (j <= n && abs(tt[j] - target) <= 0.5) { # 0.5 s tolerance; adjust if needed
        dx <- xx[j] - xx[i]
        dy <- yy[j] - yy[i]
        vals <- c(vals, dx*dx + dy*dy)
      }
    }
    if (length(vals) < 2) return(NA_real_)
    mean(vals, na.rm=TRUE)
  })
  tibble(lag = lags, TAMSD = unlist(out))
}

# Helper: velocity autocorrelation for one track (2D dot product form)
vacf_one_track <- function(tt, xx, yy, max_lag = 10) {
  dt <- diff(tt)
  if (!all(is.finite(dt)) || min(dt) <= 0) return(NULL)
  vx <- diff(xx)/dt
  vy <- diff(yy)/dt
  v2 <- vx*vx + vy*vy
  if (length(vx) < (max_lag+2)) return(NULL)
  
  base <- mean(v2, na.rm=TRUE)
  if (!is.finite(base) || base == 0) return(NULL)
  
  L <- min(max_lag, length(vx)-1)
  ac <- sapply(0:L, function(k){
    if (k == 0) return(1.0)
    i1 <- 1:(length(vx)-k); i2 <- (1+k):length(vx)
    num <- mean(vx[i1]*vx[i2] + vy[i1]*vy[i2], na.rm=TRUE)
    num/base
  })
  tibble(lag = 0:L, vacf = ac)
}

# Helper: displacement samples at a specific lag (Van Hove proxy)
displacements_at_lag <- function(tt, xx, yy, lag_s = 5, tol = 0.5) {
  n <- length(tt)
  out <- list(dx = c(), dy = c())
  for (i in 1:(n-1)) {
    target <- tt[i] + lag_s
    jrel <- which.min(abs(tt[(i+1):n] - target))
    j <- jrel + i
    if (j <= n && abs(tt[j] - target) <= tol) {
      out$dx <- c(out$dx, xx[j]-xx[i])
      out$dy <- c(out$dy, yy[j]-yy[i])
    }
  }
  out
}

## =========================================================
## =========== 1) BOOTSTRAP CIs FOR MSD SLOPE m ===========
## =========================================================

# Fit slope m from EAMSD over a time window (default <= 30 s, like your inset)
fit_m_by_group <- function(eamsd, tmax = 30) {
  eamsd %>%
    filter(Time_Adjusted > 0, Time_Adjusted <= tmax, is.finite(MSD), MSD > 0) %>%
    group_by(Species, TimeGroup) %>%
    do({
      mdl <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
      tibble(m_hat = unname(coef(mdl)[2]))
    }) %>%
    ungroup()
}

# Bootstrap m by resampling TrackIDs with replacement inside each (Species, TimeGroup)
boot_m_by_group <- function(df, nboot = 1000, tmax = 30) {
  set.seed(123)
  out <- df %>%
    filter(Time_Adjusted > 0) %>%
    group_by(Species, TimeGroup) %>%
    group_modify(~{
      ids <- unique(.x$TrackID)
      if (length(ids) < 2) return(tibble(ci_lo=NA, ci_hi=NA, p_two_sided=NA))
      boot_m <- replicate(nboot, {
        samp <- sample(ids, replace = TRUE)
        sub  <- .x %>% filter(TrackID %in% samp)
        eam  <- sub %>% group_by(Time_Adjusted) %>%
          summarise(MSD = mean(Displacement_Squared, na.rm=TRUE), .groups="drop")
        eam <- eam %>% filter(Time_Adjusted <= tmax, MSD > 0, is.finite(MSD))
        if (nrow(eam) < 3) return(NA_real_)
        coef(lm(log10(MSD) ~ log10(Time_Adjusted), data = eam))[2]
      })
      boot_m <- boot_m[is.finite(boot_m)]
      if (length(boot_m) < 20) return(tibble(ci_lo=NA, ci_hi=NA, p_two_sided=NA))
      # two-sided bootstrap p-value for H0: m = 1
      p <- 2 * min(mean(boot_m <= 1, na.rm=TRUE), mean(boot_m >= 1, na.rm=TRUE))
      tibble(ci_lo = quantile(boot_m, 0.025), ci_hi = quantile(boot_m, 0.975), p_two_sided = p)
    }) %>% ungroup()
  
  # attach point estimates
  pt <- fit_m_by_group(eammsd_tbl(df), tmax = tmax)
  left_join(pt, out, by = c("Species","TimeGroup")) %>%
    arrange(Species, TimeGroup)
}

msd_boot_results <- boot_m_by_group(data_prepped, nboot = 1000, tmax = 30)
print(msd_boot_results)

## =========================================================
## =============== 2) VELOCITY AUTOCORRELATION =============
## =========================================================

vacf_summary <- data_prepped %>%
  arrange(Species, TimeGroup, TrackID, Time_Adjusted) %>%
  group_by(Species, TimeGroup, TrackID) %>%
  group_modify(~{
    v <- vacf_one_track(.x$Time_Adjusted, .x$X, .x$Y, max_lag = 10)
    if (is.null(v)) return(tibble())
    v
  }) %>%
  ungroup() %>%
  group_by(Species, TimeGroup, lag) %>%
  summarise(mean_vacf = mean(vacf, na.rm=TRUE),
            sd_vacf   = sd(vacf,  na.rm=TRUE),
            n_tr      = n(), .groups="drop")

# Quick plot (saved)
p_vacf <- ggplot(vacf_summary, aes(lag, mean_vacf, color = TimeGroup)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() + geom_point() +
  facet_wrap(~Species, nrow=1) +
  labs(y = "Mean VACF", x = "Lag (frames)") +
  theme_minimal()
ggsave("VACF_by_Species_TimeGroup.pdf", p_vacf, width = 10, height = 4)

## =========================================================
## ================= 3) VAN HOVE (TAILS) ===================
## =========================================================
## We test normality of dx, dy at a chosen lag using KS against N(μ,σ^2).
## This is a pragmatic heavy-tail proxy (Gaussian vs non-Gaussian).

target_lag <- 5     # seconds; change or loop as needed
tol_sec    <- 0.5   # allowable time mismatch

vanhove_tbl <- data_prepped %>%
  arrange(Species, TimeGroup, TrackID, Time_Adjusted) %>%
  group_by(Species, TimeGroup, TrackID) %>%
  group_modify(~{
    d <- displacements_at_lag(.x$Time_Adjusted, .x$X, .x$Y, lag_s = target_lag, tol = tol_sec)
    if (length(d$dx) < 5) return(tibble())
    tibble(dx = d$dx, dy = d$dy)
  }) %>%
  ungroup()

# KS vs Normal for dx, dy per (Species, TimeGroup)
ks_results <- vanhove_tbl %>%
  pivot_longer(cols = c(dx,dy), names_to = "axis", values_to = "d") %>%
  group_by(Species, TimeGroup, axis) %>%
  summarise(
    n = n(),
    mu = mean(d), sd = sd(d),
    ks_p = {
      if (sd(d) == 0 || length(d) < 10) NA_real_
      else ks.test(d, "pnorm", mean(d), sd(d))$p.value
    },
    excess_kurt = {
      m2 <- mean( (d - mean(d))^2 ); m4 <- mean( (d - mean(d))^4 )
      if (is.na(m2) || m2 == 0) NA_real_ else m4/(m2^2) - 3
    },
    .groups = "drop"
  ) %>%
  mutate(TimeGroup = factor(TimeGroup, levels = lvl_ords, labels = lbl_ords))

print(ks_results)

## Quick viz of displacement PDFs (dx) at target lag
p_vanhove <- vanhove_tbl %>%
  group_by(Species, TimeGroup) %>%
  mutate(TimeGroup = factor(TimeGroup, levels = lvl_ords, labels = lbl_ords)) %>%
  ggplot(aes(dx, color = TimeGroup)) +
  geom_density() +
  facet_wrap(~Species, scales = "free_y", nrow=1) +
  labs(title = paste0("Van Hove proxy at ", target_lag, " s (dx)"),
       x = "dx (µm)", y = "Density") +
  theme_minimal()
ggsave("VanHove_dx_density.pdf", p_vanhove, width = 10, height = 4)

## =========================================================
## ======= 4) TIME- VS ENSEMBLE-MSD (ERGODICITY) ==========
## =========================================================

# Choose lags for TAMSD comparison (seconds), within 1..30 to match slope window
lags <- c(1,2,5,10,20,30)

tamsd_all <- data_prepped %>%
  arrange(Species, TimeGroup, TrackID, Time_Adjusted) %>%
  group_by(Species, TimeGroup, TrackID) %>%
  group_modify(~tamsd_one_track(.x$Time_Adjusted, .x$X, .x$Y, lags)) %>%
  ungroup()

# EAMSD at same lags (nearest within tol)
nearest_eamsd <- eammsd_tbl(data_prepped) %>%
  group_by(Species, TimeGroup) %>%
  group_modify(~{
    map_dfr(lags, function(L){
      jj <- which.min(abs(.x$Time_Adjusted - L))
      tibble(lag = L, EAMSD = .x$MSD[jj])
    })
  }) %>% ungroup()

ergod_tbl <- tamsd_all %>%
  group_by(Species, TimeGroup, lag) %>%
  summarise(TAMSD_mean = mean(TAMSD, na.rm=TRUE),
            TAMSD_sd   = sd(TAMSD,  na.rm=TRUE),
            .groups="drop") %>%
  left_join(nearest_eamsd, by = c("Species","TimeGroup","lag")) %>%
  mutate(ratio_TA_over_EA = TAMSD_mean / EAMSD) %>%
  arrange(Species, TimeGroup, lag)

print(ergod_tbl)

# Plot ratio ~ 1 for ergodic Brownian
p_ergod <- ggplot(ergod_tbl, aes(lag, ratio_TA_over_EA, color = TimeGroup)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() + geom_point() +
  facet_wrap(~Species, nrow=1) +
  scale_x_continuous(breaks = lags) +
  labs(y = "TAMSD / EAMSD", x = "Lag (s)") +
  theme_minimal()
ggsave("Ergodicity_TAMSD_over_EAMSD.pdf", p_ergod, width = 10, height = 4)

## =========================================================
## ======= 5) PIECEWISE SLOPES (AND BY HMM STATE) =========
## =========================================================

# Piecewise slope on EAMSD for windows (seconds):
windows <- list(
  early = c(1, 10),
  mid   = c(10, 30),
  late  = c(30, 100)
)

piecewise_slopes <- function(eamsd, wins = windows) {
  bind_rows(lapply(names(wins), function(nm){
    w <- wins[[nm]]
    eamsd %>%
      filter(Time_Adjusted >= w[1], Time_Adjusted <= w[2],
             is.finite(MSD), MSD > 0) %>%
      group_by(Species, TimeGroup) %>%
      do({
        if (nrow(.) < 3) return(tibble(window = nm, m = NA_real_))
        mdl <- lm(log10(MSD) ~ log10(Time_Adjusted), data = .)
        tibble(window = nm, m = unname(coef(mdl)[2]))
      }) %>% ungroup()
  }))
}

pw_eamsd <- piecewise_slopes(eammsd_tbl(data_prepped), windows)
print(pw_eamsd)

# OPTIONAL: If you have per-point HMM states in data_prepped (e.g., column 'state' or 'State'),
# compute state-wise EAMSD slopes over the same windows.
state_col <- intersect(names(data_prepped), c("state","State"))
if (length(state_col) == 1) {
  nm <- state_col[[1]]
  eamsd_state <- data_prepped %>%
    filter(Time_Adjusted > 0) %>%
    group_by(Species, TimeGroup, !!as.name(nm), Time_Adjusted) %>%
    summarise(MSD = mean(Displacement_Squared), .groups="drop") %>%
    rename(State = !!as.name(nm))
  pw_state <- piecewise_slopes(eamsd_state, windows) %>%
    rename(State = TimeGroup)  # harmless rename if needed
  print(pw_state)
}

## =========================================================
## ================== EXPORT ALL SUMMARIES =================
## =========================================================

out_eamsd <- eammsd_tbl(data_prepped)
to_write <- list(
  "MSD_bootstrap_m"       = msd_boot_results,
  "VACF_summary"          = vacf_summary,
  "VanHove_KS"            = ks_results,
  "Ergodicity_ratio"      = ergod_tbl,
  "Piecewise_slopes"      = pw_eamsd,
  "EAMSD_table"           = out_eamsd
)
write_xlsx(to_write, "MSD_quick_checks_summary.xlsx")
cat("✅ Saved: VACF_by_Species_TimeGroup.pdf, VanHove_dx_density.pdf, Ergodicity_TAMSD_over_EAMSD.pdf, MSD_quick_checks_summary.xlsx\n")


################################################################################
# Compare Slopes & Stats

test_slope_difference <- function(msd_data, time_group, time_limit = 30) {
  
  data_subset <- msd_data %>%
    filter(TimeGroup == time_group, Time_Adjusted > 0, Time_Adjusted <= time_limit)
  
  if(nrow(data_subset) < 5) {
    return(list(Slope = NA, Intercept = NA, R_Squared = NA, P_Value = NA, Cohen_f2 = NA))
  }
  
  
  model <- lm(log10(MSD) ~ log10(Time_Adjusted), data = data_subset)
  
  
  hypothesis_test <- linearHypothesis(model, "log10(Time_Adjusted) = 1")
  
  
  r_sq <- summary(model)$r.squared
  
  f2_val <- if(!is.na(r_sq) && r_sq < 1) r_sq / (1 - r_sq) else NA
  
  list(
    Slope      = coef(model)[2],
    Intercept  = coef(model)[1],
    R_Squared  = r_sq,
    P_Value    = hypothesis_test$`Pr(>F)`[2],
    Cohen_f2   = f2_val
  )
}


# RUn Slope Tests

species_list <- unique(data_prepped$Species)
results_list <- list()

for (sp in species_list) {
  
  # Calculate MSD curves
  msd_data <- calculate_msd_curves(data_prepped, sp)
  time_groups_order <- levels(msd_data$TimeGroup)
  
  for (tg in time_groups_order) {
    test_res <- test_slope_difference(msd_data, tg, time_limit = 30)
    results_list[[length(results_list) + 1]] <- data.frame(
      Species   = sp,
      TimeGroup = tg,
      Slope     = round(test_res$Slope, 2),
      Intercept = round(test_res$Intercept, 2),
      R_Squared = round(test_res$R_Squared, 2),
      P_Value   = test_res$P_Value,
      Cohen_f2  = if(!is.na(test_res$Cohen_f2)) round(test_res$Cohen_f2, 2) else NA
    )
  }
}

final_results <- bind_rows(results_list)

# Export results to Excel
write_xlsx(final_results, "Stat_tests_MSD.xlsx")
cat("✅ Slope test results (with R^2 and Cohen's f^2) exported to 'Stat_tests_MSD.xlsx'\n")

