################################################################################
##
## TOA-BASED ACTIVITY ESTIMATION - GROUPER DATA 2021-2023
##
################################################################################

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load libraries
library(data.table)
library(lubridate)
library(suncalc)
library(glmmTMB)
library(gratia)
library(performance)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(patchwork)
library(scales)
library(purrr)
library(emmeans)
library(tidyr)
library(ggplot2)
library(stringr)

# Load tracking data (5 minute intervals)
detect <- fread("./tracking_toa_5min.csv")
fish_ref <- fread("./tracking_fish_metadata.csv")
head(detect)

# Number of observations by fish
table(detect$fish_id)

# Load environmental data
copern <- fread("./environmentals.csv")

# Function to add daynight period to data
getDayNight <- function(times, lat, long) {
  date_range <- range(as.Date(times))
  dates <- seq(date_range[1] - 1, date_range[2] + 1, "days")
  sun <- getSunlightTimes(dates, lat = lat, lon = long,
                          keep = c("sunrise", "sunset"))
  indx_dn <- match(as.Date(times), sun$date)
  return(ifelse(times >= sun$sunrise[indx_dn] & times <= sun$sunset[indx_dn],
                "day", "night"))
}

# 1. PREPARE DETECTION DATA ====================================================

detect[, fish_id := factor(fish_id, levels = sort(unique(fish_id)))]

# Center distances (daily mean for each deployment = 0 m) ----------------------
detect[, norm_dist := dist - mean(dist, na.rm = TRUE),
       by = list(fish_id, deployment_id, as.Date(date_time))]

# Plot one individual
id <- "EM0003"
det_sub <- detect[fish_id == id]
table(det_sub$station)

# Original distances
plot(dist ~ date_time, data = det_sub[station == "TOR05"], type = "l", 
     ylim = range(det_sub$dist), ylab = "Rel. distance (m)")
lines(dist ~ date_time, data = det_sub[station == "TOR04"], col = 2)
lines(dist ~ date_time, data = det_sub[station == "TOR03"], col = 3)

# Centered distances
plot(norm_dist ~ date_time, data = det_sub[station == "TOR05"], type = "l",
     ylim = range(det_sub$norm_dist), ylab = "Norm. dist. (m)")
lines(norm_dist ~ date_time, data = det_sub[station == "TOR04"], col = 2)
lines(norm_dist ~ date_time, data = det_sub[station == "TOR03"], col = 3)

# Steps
plot(step ~ date_time, data = det_sub[station == "TOR05"], type = "l",
     ylab = "Step length (m)")
lines(step ~ date_time, data = det_sub[station == "TOR04"], col = 2)
lines(step ~ date_time, data = det_sub[station == "TOR03"], col = 3)



# 2. HOURLY AND DAILY METRICS ==================================================

# Estimate hourly and daily metrics merging data from all the deployments
#   - n_det: number of 5-min detection bins
#   - temperature: mean temperature
#   - dist_range: range of distances (max_dist - min_dist)
#   - step_mean: mean step size
#   - step_sd: standard deviation of step length
#   - step_max: maximum step length

# Hourly metrics ---------------------------------------------------------------

hour_metrics <- detect[, .(
  n_det = .N,
  temperature = mean(temperature, na.rm = TRUE),
  shift = shift,
  shift_id = shift_id,
  dist_range = diff(range(norm_dist, na.rm = TRUE)),
  step_mean = mean(step, na.rm = TRUE),
  step_sd = sd(step, na.rm = TRUE),
  step_max = quantile(step, probs = 1, na.rm = TRUE)),
  by = list(fish_id, date_time = floor_date(date_time, "hours"))]


# Add day/night times
hour_metrics[, daynight := getDayNight(date_time, 
                                       lat = 39.46256, long = 2.47190)]
hour_metrics[, daynight := factor(daynight, levels = c("day", "night"))]

hour_metrics$sunpos <- getSunlightPosition(hour_metrics$date_time + 30*60, 
                                           lat = 39.46256, lon = 2.47190,
                                           keep = c("altitude"))$altitude


# Add Copernicus environmental data
indx_c <- match(as.Date(detect$date_time), copern$day)
hour_metrics$wave_height <- copern$wave_height[indx_c]
hour_metrics$cloud_cover <- copern$cloud_cover[indx_c]
hour_metrics$moon <- copern$moon[indx_c]

###Add divers data

data_divers <- fread("./divers.csv")

data_divers <- data_divers %>%
  mutate(shift_time = dplyr::recode(hour,
                                    "8-11h"  = 3,
                                    "11-13h" = 2,
                                    "13-18h" = 5,
                                    "19-23h" = 4))

data_divers$divers_h <- data_divers$divers / data_divers$shift_time

hour_metrics <- hour_metrics %>%
  dplyr::left_join(data_divers %>% dplyr::select(shift_id, divers_h), by = "shift_id")
hour_metrics$divers_h[is.na(hour_metrics$divers_h)] <- 0 #missing data means no divers


hour_metrics <- hour_metrics[!is.na(hour_metrics$step_sd), ]#delete incomplete rows


##################
# ==============================================================================
#3.Monthly distribution of day/night metrics (FIGURE 3)
# ==============================================================================

# 3.1. AGGREGATE DATA
# ==============================================================================
month_metrics <- hour_metrics[, .(
  n_days = .N,
  n_detect = sum(n_det),
  temperature = mean(temperature, na.rm = TRUE),
  
  # Mean step
  step_mean_lower  = quantile(step_mean, 0.025, na.rm = TRUE),
  step_mean_q25    = quantile(step_mean, 0.25, na.rm = TRUE),
  step_mean_median = quantile(step_mean, 0.5, na.rm = TRUE),
  step_mean_q75    = quantile(step_mean, 0.75, na.rm = TRUE),
  step_mean_upper  = quantile(step_mean, 0.975, na.rm = TRUE),
  
  # Max step
  step_max_lower   = quantile(step_max, 0.025, na.rm = TRUE),
  step_max_q25     = quantile(step_max, 0.25, na.rm = TRUE),
  step_max_median  = quantile(step_max, 0.5, na.rm = TRUE),
  step_max_q75     = quantile(step_max, 0.75, na.rm = TRUE),
  step_max_upper   = quantile(step_max, 0.975, na.rm = TRUE)),
  
  by = list(fish_id, daynight, date = floor_date(date_time, "month"))]

month_metrics <- month_metrics[!(fish_id == "EM0005" & date == as.Date("2023-08-01"))] ###removed fish with detections only in 01/08/2023, for plotting


# 3.2. PRE-CALCULATE OVERALL MONTHLY TEMPERATURE (WITH BREAKS FOR NAs)
# ==============================================================================
range_dates <- range(month_metrics$date, na.rm = TRUE)
month_seq <- seq(range_dates[1], range_dates[2], by = "months")

# Calculate OVERALL monthly temperature across all fish/records
overall_temp <- hour_metrics[, .(temp = mean(temperature, na.rm = TRUE)), 
                             by = .(m_floor = floor_date(date_time, "month"))]

# Merge with full sequence of months so missing months become NA
full_months <- data.table(m_floor = month_seq, 
                          x_pos = seq_along(month_seq))
overall_temp <- merge(full_months, overall_temp, by = "m_floor", all.x = TRUE)

t_min <- floor(min(overall_temp$temp, na.rm = TRUE)) - 2
t_max <- ceiling(max(overall_temp$temp, na.rm = TRUE)) + 2



# 3.3. PLOTTING SETUP
# ==============================================================================
dir <- "./Figure 3/"
dir.create(dir, showWarnings = FALSE)

col <- rep(c("#577590", "#90be6d", "#f9c74f", "#f3722c"), each = 3)
col_l <- sapply(col, function(c) colorRampPalette(c(c, "black"))(100)[30])
col_m <- sapply(col, function(c) colorRampPalette(c(c, "black"))(100)[70])

month_lab <- format(month_seq, "%b")
month_metrics$month_f <- as.numeric(factor(as.character(month_metrics$date), 
                                           levels = as.character(month_seq)))
month_metrics$month_num <- month(month_metrics$date)

vars <- c("step_mean", "step_max")
ylabs <- c(
  expression("Daily " ~ step[mean] ~ " (m)"),
  expression("Daily " ~ step[max] ~ " (m)")
)


for (i in seq_along(vars)) {
  
  file_name <- paste0(dir, vars[i], ".png")
  xlim <- c(1, length(month_seq))
  
  png(file_name, width = 2600, height = 2600, res = 250)
  
  par(mfrow = c(4, 2), mar = c(3.8, 4.8, 0.5, 3.5))
  
  for (id in sort(unique(month_metrics$fish_id))) {
    
    indx <- month_metrics$fish_id == id & month_metrics$n_days > 5
    
    if(sum(indx) == 0) {
      plot.new()
      next
    }
    
    mid   <- month_metrics[[paste0(vars[i], "_median")]]
    w_low <- month_metrics[[paste0(vars[i], "_lower")]]
    low   <- month_metrics[[paste0(vars[i], "_q25")]]
    upp   <- month_metrics[[paste0(vars[i], "_q75")]]
    w_upp <- month_metrics[[paste0(vars[i], "_upper")]]
    
    ylim <- range(c(0, w_low[indx], w_upp[indx]), na.rm = TRUE)
    ylim <- ylim + c(0, 1) * diff(ylim)/8
    
    x_center <- month_metrics$month_f[indx] + 
      ifelse(month_metrics$daynight[indx] == "day", -0.18, 0.18)
    
    m_num <- month_metrics$month_num[indx]
    fill_col <- ifelse(month_metrics$daynight[indx] == "day", col[m_num], "white")
    
    # --- A) Base Empty Plot ---
    plot(mid[indx] ~ x_center, type = "n", ylim = ylim, main = id, 
         xlim = xlim, axes = FALSE, ann = FALSE)
    
    # --- B) OVERALL MONTHLY TEMPERATURE BACKGROUND ---
    par(new = TRUE)
    plot(overall_temp$x_pos, overall_temp$temp, type = "l", 
         col = "tomato", lwd = 2, 
         xlim = xlim, ylim = c(t_min, t_max),
         axes = FALSE, xlab = "", ylab = "")
    
    axis(4, las = 1, col = "#d90429", col.axis = "#d90429", lwd = 0, lwd.ticks = 1)
    mtext("Temperature (°C)", side = 4, line = 2.2, col = "#d90429", cex = 0.7)
    
    # Reset par to map the boxplots/rectangles back to the primary Y-limits
    par(new = TRUE)
    plot(mid[indx] ~ x_center, type = "n", ylim = ylim, 
         xlim = xlim, axes = FALSE, ann = FALSE)
    
    # --- C) Gray rectangle for gap ---
    # Put back here so it appears for all fish, crossing exactly over missing temperature
    rect(xleft = 17.8 - 0.4, xright = 19.2 + 0.4, ybottom = par("usr")[3], 
         ytop = par("usr")[4], col = "gray90", border = FALSE)
    
    # --- D) DRAW BOXPLOTS ON TOP ---
    arrows(x0 = x_center, x1 = x_center,
           y0 = w_low[indx], y1 = w_upp[indx],
           code = 3, angle = 90, length = 0, col = col_l[m_num], lwd = 1.2)
    
    rect(xleft = x_center - 0.15, xright = x_center + 0.15,
         ybottom = low[indx], ytop = upp[indx], lwd = 1.8,
         col = fill_col, border = col_l[m_num])
    
    rect(xleft = x_center - 0.15, xright = x_center + 0.15,
         ybottom = mid[indx], ytop = mid[indx], lwd = 1.8,
         col = col_m[m_num], border = col_m[m_num])
    
    # --- E) MAIN AXES AND LABELS ---
    axis(2, las = 1, lwd = 0, lwd.ticks = 1)
    lab_at <- pretty(xlim, 15)
    lab_at <- lab_at[lab_at > 0 & lab_at <= xlim[2]]
    year_at <- which(month(month_seq) == 1)
    
    axis(1, at = lab_at, lwd = 0, lwd.ticks = 1, 
         labels = month_lab[lab_at], padj = -0.5)
    axis(1, at = year_at, lwd = 0, lwd.ticks = 1, line = 2, tcl = -1,
         labels = year(month_seq)[year_at], padj = -1.5, hadj = -0.15)
    axis(1, at = seq_along(month_lab)[!seq_along(month_lab) %in% lab_at], 
         lwd = 0, lwd.ticks = 1, labels = FALSE, tcl = -0.25)
    
    title(ylab = ylabs[[i]], line = 2.8)
    box()
    
    text(par("usr")[1], par("usr")[4], adj = c(-0.2, 1.8), labels = id, font = 2)
    
    # --- F) LEGENDS ---
    legend("topright", legend = c("Winter", "Spring", "Summer", "Autumn"),
           fill = unique(col), border = unique(col_l), box.col = "gray60",
           inset = c(0.015, 0.03), cex = 0.8)
    
    legend("topleft", legend = c("Day", "Night"),
           fill = c("gray50", "white"), border = "gray50", box.col = "gray60",
           inset = c(0.15, 0.03), cex = 0.8)
  }
  
  dev.off()
}


##########################################################
# ==============================================================================
#4. Modeled day/night activity (Figure 4)
# ==============================================================================

#4.1. Extract slopes and create the models 
# ==============================================================================

hour_metrics <- hour_metrics %>%
  mutate(
    month_num = lubridate::month(date_time),
    month_f = factor(month.abb[month_num], levels = month.abb)
  )

# Separate the dataset into Day and Night
day_data   <- hour_metrics %>% filter(daynight == "day")
night_data <- hour_metrics %>% filter(daynight == "night")

# Fit intercept-only models on the full dataset to establish baseline (0%)
baseline_model_max <- glmmTMB(step_max ~ 1 + (1|fish_id), 
                              data = hour_metrics, family = Gamma(link = "log"))
baseline_model_mean <- glmmTMB(step_mean ~ 1 + (1|fish_id), 
                               data = hour_metrics, family = Gamma(link = "log"))

global_baseline_max  <- unname(exp(fixef(baseline_model_max)$cond[1]))
global_baseline_mean <- unname(exp(fixef(baseline_model_mean)$cond[1]))

# Daytime Models
day_model_max <- glmmTMB(step_max ~ month_f + temperature + divers_h + 
                           wave_height + moon + cloud_cover + (1|fish_id),
                         data = day_data, family = Gamma(link = "log"))

day_model_mean <- glmmTMB(step_mean ~ month_f + temperature + divers_h + 
                            wave_height + moon + cloud_cover + (1|fish_id),
                          data = day_data, family = Gamma(link = "log"))

# Nighttime Models
night_model_max <- glmmTMB(step_max ~ month_f + temperature + divers_h + 
                             wave_height + moon + cloud_cover + (1|fish_id),
                           data = night_data, family = Gamma(link = "log"))

night_model_mean <- glmmTMB(step_mean ~ month_f + temperature + divers_h + 
                              wave_height + moon + cloud_cover + (1|fish_id),
                            data = night_data, family = Gamma(link = "log"))


# 4.2. Extract predictions and percentage change
# ==============================================================================
# Extract predictions from all 4 models
preds_day_max   <- suppressWarnings(ggpredict(day_model_max, terms = "month_f"))
preds_night_max <- suppressWarnings(ggpredict(night_model_max, terms = "month_f"))

preds_day_mean   <- suppressWarnings(ggpredict(day_model_mean, terms = "month_f"))
preds_night_mean <- suppressWarnings(ggpredict(night_model_mean, terms = "month_f"))

# step_max dataframes
df_day_max <- as.data.frame(preds_day_max) %>%
  mutate(group = "day", Metric = "max")

df_night_max <- as.data.frame(preds_night_max) %>%
  mutate(group = "night", Metric = "max")

# step_mean dataframes
df_day_mean <- as.data.frame(preds_day_mean) %>%
  mutate(group = "day", Metric = "mean")

df_night_mean <- as.data.frame(preds_night_mean) %>%
  mutate(group = "night", Metric = "mean")

# Combine Day/Night and calculate % change
df_combined_max <- bind_rows(df_day_max, df_night_max) %>%
  mutate(
    perc_change = (predicted - global_baseline_max) / global_baseline_max,
    conf.low_perc = (conf.low - global_baseline_max) / global_baseline_max,
    conf.high_perc = (conf.high - global_baseline_max) / global_baseline_max
  ) %>%
  mutate(group = factor(group, levels = c("day", "night"))) 

df_combined_mean <- bind_rows(df_day_mean, df_night_mean) %>%
  mutate(
    perc_change = (predicted - global_baseline_mean) / global_baseline_mean,
    conf.low_perc = (conf.low - global_baseline_mean) / global_baseline_mean,
    conf.high_perc = (conf.high - global_baseline_mean) / global_baseline_mean
  ) %>%
  mutate(group = factor(group, levels = c("day", "night"))) 


# 4.3. Create the plots
# ==============================================================================
season_cols <- rep(c("#577590", "#90be6d", "#f9c74f", "#f3722c"), each = 3)
names(season_cols) <- month.abb

plot_discrete_month_perc <- function(df, y_label) {
  
  ggplot(df, aes(x = x, y = perc_change, group = group)) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    
    geom_line(aes(linetype = group), 
              color = "gray20", position = position_dodge(width = 0.5), 
              alpha = 0.5, linewidth = 0.8) +
    geom_pointrange(aes(ymin = conf.low_perc, ymax = conf.high_perc, 
                        color = x, shape = group),
                    position = position_dodge(width = 0.5), size = 0.8, linewidth = 1) +
    
    scale_y_continuous(labels = scales::percent_format()) +
    scale_color_manual(values = season_cols, guide = "none") + 
    scale_shape_manual(values = c("day" = 16, "night" = 1), name = "Period") +
    scale_linetype_manual(values = c("day" = "solid", "night" = "dotted"), name = "Period") +
    
    labs(x = NULL, y = y_label) +
    
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

p_month_mean <- plot_discrete_month_perc(
  df = df_combined_mean, 
  y_label = expression("Effect Size (%)")
)

p_month_max <- plot_discrete_month_perc(
  df = df_combined_max, 
  y_label = expression("Effect Size (%)")
)

combined_monthly_grid <- (p_month_mean / p_month_max) +
  plot_annotation(title = "Monthly day and night effect size") +
  plot_layout(guides = "collect") & 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "bottom"
  )

print(combined_monthly_grid)



# ==============================================================================
# 5. Modeled day/night effects of environmental variables (Figure 5)
# ==============================================================================
# 5.1. Split datasets by seasons and day and night
# ==============================================================================
dir <- "./Figure 5/"
dir.create(dir, showWarnings = FALSE)

hour_metrics <- hour_metrics %>%
  mutate(
    season = case_when(
      (format(date_time, "%m-%d") >= "12-21" | format(date_time, "%m-%d") < "03-20") ~ "Winter",
      (format(date_time, "%m-%d") >= "03-20" & format(date_time, "%m-%d") < "06-21") ~ "Spring",
      (format(date_time, "%m-%d") >= "06-21" & format(date_time, "%m-%d") < "09-22") ~ "Summer",
      (format(date_time, "%m-%d") >= "09-22" & format(date_time, "%m-%d") < "12-21") ~ "Autumn"
    )
  )


day_data   <- subset(hour_metrics, daynight == "day")
night_data <- subset(hour_metrics, daynight == "night")

list_day <- list(
  Winter  = subset(day_data, season == "Winter"),
  Spring  = subset(day_data, season == "Spring"),
  Summer  = subset(day_data, season == "Summer"),
  Autumn  = subset(day_data, season == "Autumn"),
  Overall = day_data
)

list_night <- list(
  Winter  = subset(night_data, season == "Winter"),
  Spring  = subset(night_data, season == "Spring"),
  Summer  = subset(night_data, season == "Summer"),
  Autumn  = subset(night_data, season == "Autumn"),
  Overall = night_data
)


# 2. Fit 20 models
# ==============================================================================
f_mean_season  <- step_mean ~ temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover + (1|fish_id)
f_mean_overall <- step_mean ~ season + temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover + (1|fish_id)

f_max_season  <- step_max ~ temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover + (1|fish_id)
f_max_overall <- step_max ~ season + temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover + (1|fish_id)

mean_models_day   <- list(Winter=glmmTMB(f_mean_season, data=list_day$Winter, family=Gamma(link="log")), Spring=glmmTMB(f_mean_season, data=list_day$Spring, family=Gamma(link="log")), Summer=glmmTMB(f_mean_season, data=list_day$Summer, family=Gamma(link="log")), Autumn=glmmTMB(f_mean_season, data=list_day$Autumn, family=Gamma(link="log")), Overall=glmmTMB(f_mean_overall, data=list_day$Overall, family=Gamma(link="log")))
max_models_day    <- list(Winter=glmmTMB(f_max_season, data=list_day$Winter, family=Gamma(link="log")), Spring=glmmTMB(f_max_season, data=list_day$Spring, family=Gamma(link="log")), Summer=glmmTMB(f_max_season, data=list_day$Summer, family=Gamma(link="log")), Autumn=glmmTMB(f_max_season, data=list_day$Autumn, family=Gamma(link="log")), Overall=glmmTMB(f_max_overall, data=list_day$Overall, family=Gamma(link="log")))

mean_models_night <- list(Winter=glmmTMB(f_mean_season, data=list_night$Winter, family=Gamma(link="log")), Spring=glmmTMB(f_mean_season, data=list_night$Spring, family=Gamma(link="log")), Summer=glmmTMB(f_mean_season, data=list_night$Summer, family=Gamma(link="log")), Autumn=glmmTMB(f_mean_season, data=list_night$Autumn, family=Gamma(link="log")), Overall=glmmTMB(f_mean_overall, data=list_night$Overall, family=Gamma(link="log")))
max_models_night  <- list(Winter=glmmTMB(f_max_season, data=list_night$Winter, family=Gamma(link="log")),  Spring=glmmTMB(f_max_season, data=list_night$Spring, family=Gamma(link="log")),  Summer=glmmTMB(f_max_season, data=list_night$Summer, family=Gamma(link="log")),  Autumn=glmmTMB(f_max_season, data=list_night$Autumn, family=Gamma(link="log")),  Overall=glmmTMB(f_max_overall, data=list_night$Overall, family=Gamma(link="log")))


# 3. Calculate slopes and curves
# ==============================================================================
seas_vec <- c("Winter", "Spring", "Summer", "Autumn", "Overall")

# --- Extract Linear Slopes for Wave, Cloud, Moon ---
get_subset_slopes <- function(mod_list, metric_label, dn_label, var_name) {
  imap_dfr(mod_list, ~ {
    tr_df <- emtrends(.x, specs = ~ 1, var = var_name, type = "link") %>% as.data.frame()
    trend_col <- paste0(var_name, ".trend")
    tibble(season=.y, metric=metric_label, daynight=dn_label, var=var_name, 
           slope=as.numeric(tr_df[[trend_col]]), lwr=as.numeric(tr_df$asymp.LCL), upr=as.numeric(tr_df$asymp.UCL))
  })
}

vars_to_plot <- c("cloud_cover", "wave_height", "moon")
effects_linear <- map_dfr(vars_to_plot, function(v) {
  bind_rows(
    get_subset_slopes(mean_models_day, "mean", "day", v), 
    get_subset_slopes(max_models_day, "max", "day", v), 
    get_subset_slopes(mean_models_night, "mean", "night", v), 
    get_subset_slopes(max_models_night, "max", "night", v)
  )
}) %>%
  mutate(
    season = factor(season, levels = seas_vec),
    metric = factor(metric, levels = c("mean", "max")),
    perc_change = 100*(exp(slope)-1), 
    lwr_perc = 100*(exp(lwr)-1), 
    upr_perc = 100*(exp(upr)-1)
  )

# --- Extract Prediction Curves for Temperature 
get_temperature_curves <- function(mod_list, data_list, metric_label, dn_label) {
  imap_dfr(mod_list, function(model, season_name) {
    if(season_name == "Overall") return(NULL) 
    
    preds <- suppressWarnings(ggpredict(model, terms = "temperature [all]") %>% as.data.frame())
    
    model_baseline <- mean(preds$predicted)
    
    # Remove predictions above 24°C for spring as they are few (last day of spring only), with a misguiding impact on plot
    if (season_name == "Spring") {
      preds <- preds %>% filter(x <= 24)
    }
    
    preds %>% mutate(
      season = season_name, 
      metric = metric_label, 
      daynight = dn_label, 
      perc_change = (predicted - model_baseline) / model_baseline, 
      lwr_perc = (conf.low - model_baseline) / model_baseline, 
      upr_perc = (conf.high - model_baseline) / model_baseline
    )
  })
}

temp_curves <- bind_rows(
  get_temperature_curves(mean_models_day, list_day, "mean", "day"), 
  get_temperature_curves(max_models_day, list_day, "max", "day"),
  get_temperature_curves(mean_models_night, list_night, "mean", "night"), 
  get_temperature_curves(max_models_night, list_night, "max", "night")
) %>% mutate(season=factor(season, levels=seas_vec), metric=factor(metric, levels=c("mean", "max")))


# 4. Plotting 
# ==============================================================================
base_custom_theme <- function(is_night) {
  t <- theme_minimal() + 
    theme(
      legend.position = "none", 
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      panel.grid.major = element_line(linewidth = 0.5, color = "gray85"), 
      panel.grid.minor = element_line(linewidth = 0.4, color = "gray90")  
    )
  if (is_night) {
    t <- t + theme(
      panel.background = element_rect(fill = "grey90", color = NA),
      panel.grid.major = element_line(linewidth = 0.8, color = "white"),
      panel.grid.minor = element_line(linewidth = 0.4, color = "white")
    )
  }
  return(t)
}

# Plot Temperature Curves
plot_temp <- function(df, target_dn, target_metric) {
  d <- df %>% filter(daynight == target_dn, metric == target_metric)
  m_lab <- ifelse(target_metric == "mean", "Mean", "Max")
  plot_title <- paste0("Temperature (°C)\n(", str_to_title(target_dn), " - ", m_lab, ")")
  
  ggplot(d, aes(x=x, y=perc_change * 100, color=season, fill=season)) +
    geom_hline(yintercept=0, color="gray50", linetype="dashed") +
    geom_ribbon(aes(ymin=lwr_perc * 100, ymax=upr_perc * 100), alpha=0.15, linetype=0) +
    geom_line(linewidth=1) +
    scale_color_manual(values=c("Winter"="#577590", "Spring"="#90be6d", "Summer"="#f9c74f", "Autumn"="#f3722c")) +
    scale_fill_manual(values=c("Winter"="#577590", "Spring"="#90be6d", "Summer"="#f9c74f", "Autumn"="#f3722c")) +
    labs(x="Temperature (°C)", y="Effect size (%)", title=NULL, color="Season", fill="Season") +
    base_custom_theme(target_dn == "night")
}

# Plot Dot-and-Whisker Linear Slopes
plot_linear <- function(df, target_var, target_dn) {
  var_titles <- c(wave_height="Wave height (m)", moon="Lunar illuminated fraction", cloud_cover="Cloud cover fraction")
  d <- df %>% filter(var == target_var, daynight == target_dn)
  plot_title <- paste0(var_titles[[target_var]], "\n(", str_to_title(target_dn), ")")
  
  ggplot(d, aes(x=season, y=perc_change, color=metric)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray50") +
    geom_pointrange(aes(ymin=lwr_perc, ymax=upr_perc), position=position_dodge(width=0.5), fatten=3, linewidth=1) +
    scale_colour_manual(values=c("mean"="blue", "max"="red"), labels=c(mean=expression(step[mean]), max=expression(step[max]))) +
    scale_x_discrete(limits=seas_vec) +
    labs(x=NULL, y="Effect size (%)", title=NULL, color="") +
    base_custom_theme(target_dn == "night") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
}

# 5. Generate Panels & Combine
# ==============================================================================
p_temp_mean_day   <- plot_temp(temp_curves, "day", "mean")
p_temp_mean_night <- plot_temp(temp_curves, "night", "mean")
p_temp_max_day    <- plot_temp(temp_curves, "day", "max")
p_temp_max_night  <- plot_temp(temp_curves, "night", "max")

p_wave_day        <- plot_linear(effects_linear, "wave_height", "day")
p_wave_night      <- plot_linear(effects_linear, "wave_height", "night")
p_cloud_day       <- plot_linear(effects_linear, "cloud_cover", "day")
p_cloud_night     <- plot_linear(effects_linear, "cloud_cover", "night")
p_moon_day        <- plot_linear(effects_linear, "moon", "day")
p_moon_night      <- plot_linear(effects_linear, "moon", "night")

# Combine layout
final_plot <- (p_temp_mean_day | p_temp_mean_night) /
  (p_temp_max_day  | p_temp_max_night) /
  (p_wave_day      | p_wave_night) /
  (p_cloud_day     | p_cloud_night) /
  (p_moon_day      | p_moon_night) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.text = element_text(size=12))

print(final_plot)




# ==============================================================================
# 6. Diver effects per season (FIGURE 7)
# ==============================================================================
# 6.1. Prepare data and create models
# ==============================================================================
dive_hour_metrics <- subset(hour_metrics, shift %in% c("8-11h", "11-13h", "13-18h", "19-23h")) #only diving times are kept for diving analyses

winter <- subset(dive_hour_metrics, season == "Winter")
spring <- subset(dive_hour_metrics, season == "Spring")
summer <- subset(dive_hour_metrics, season == "Summer")
autumn <- subset(dive_hour_metrics, season == "Autumn")

winter <- winter %>%
  filter(shift != "[19,23)") %>%    # keep all but that shift because divers is always = 0 and it interferes in model if not specified
  mutate(shift = droplevels(as.factor(shift)))

winter_mean <- glmmTMB(step_mean ~ divers_h + moon + temperature + I(temperature^2) +  wave_height + cloud_cover + shift +
                         (1 | fish_id),
                       data = winter, family = Gamma(link = "log"))

spring_mean <- glmmTMB(step_mean ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift +
                         (1 | fish_id),
                       data = spring, family = Gamma(link = "log"))

summer_mean <- glmmTMB(step_mean ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift +
                         (1 | fish_id),
                       data = summer, family = Gamma(link = "log"))
autumn_mean <- glmmTMB(step_mean ~ divers_h + moon + temperature + I(temperature^2) +  wave_height + cloud_cover + shift +
                         (1 | fish_id),
                       data = autumn, family = Gamma(link = "log"))
overall_mean <- glmmTMB(step_mean ~ divers_h + moon + temperature + I(temperature^2) +  wave_height + cloud_cover + shift + season +
                          (1 | fish_id),
                        data = hour_metrics, family = Gamma(link = "log"))

winter_max <- glmmTMB(step_max ~ divers_h + moon + temperature + I(temperature^2) +  wave_height + cloud_cover + shift +
                        (1 | fish_id),
                      data = winter, family = Gamma(link = "log"))
spring_max <- glmmTMB(step_max ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id),
                      data = spring, family = Gamma(link = "log"))
summer_max <- glmmTMB(step_max ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id),
                      data = summer, family = Gamma(link = "log"))
autumn_max <- glmmTMB(step_max ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id),
                      data = autumn, family = Gamma(link = "log"))
overall_max <- glmmTMB(step_max ~ divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift + season +
                         (1 | fish_id),
                       data = hour_metrics, family = Gamma(link = "log"))

seas_vec <- c("Winter", "Spring", "Summer", "Autumn", "Overall")

mean_models <- list(
  Winter  = winter_mean,
  Spring  = spring_mean,
  Summer  = summer_mean,
  Autumn  = autumn_mean,
  Overall = overall_mean
)

max_models <- list(
  Winter  = winter_max,
  Spring  = spring_max,
  Summer  = summer_max,
  Autumn  = autumn_max,
  Overall = overall_max
)



# 6.2. Calculate number of days with >= 1 diver 
# ==============================================================================
n_df <- hour_metrics |>
  filter(divers_h > 0) |>
  mutate(date = as.Date(date_time)) |> 
  distinct(season, date) |>
  count(season, name = "n") |>
  complete(season = seas_vec, fill = list(n = 0)) |>
  mutate(
    n = ifelse(
      season == "Overall",
      nrow(distinct(hour_metrics |> filter(divers_h > 0) |> mutate(date = as.Date(date_time)), date)),
      n
    ),
    season = factor(season, levels = seas_vec)
  )


# 6.3. Get slopes 
# ==============================================================================
plot_var <- function(var_name, mean_mods, max_mods, n_data, seas_vec, y_limits = NULL) {
  
  var_titles <- c(
    divers_h    = "Divers",
    temperature = "Temperature",
    moon        = "Moon",
    cloud_cover = "Cloud cover",
    wave_height = "Wave height"
  )
  
  # 1. Helper to extract slopes for this specific variable
  get_slopes <- function(mod_list, metric_label) {
    imap_dfr(mod_list, ~ {
      tr_df <- emtrends(.x, specs = ~ 1, var = var_name, type = "link") |> as.data.frame()
      trend_col <- paste0(var_name, ".trend")
      tibble(
        season = .y, 
        metric = metric_label, 
        slope  = as.numeric(tr_df[[trend_col]]), 
        lwr    = as.numeric(tr_df$asymp.LCL), 
        upr    = as.numeric(tr_df$asymp.UCL)
      )
    })
  }
  
  # 2. Build the data specifically for this plot
  d <- bind_rows(
    get_slopes(mean_mods, "mean"),
    get_slopes(max_mods, "max")
  ) |>
    left_join(n_data, by = "season") |>
    mutate(
      season     = factor(season, levels = seas_vec),
      slope_plot = 100 * (exp(slope) - 1),
      lwr_plot   = 100 * (exp(lwr)   - 1),
      upr_plot   = 100 * (exp(upr)   - 1)
    )
  
  # 3. Create the plot
  offset <- 0.10 * diff(range(c(d$lwr_plot, d$upr_plot), na.rm = TRUE))
  y_lab  <- min(d$lwr_plot, na.rm = TRUE) - offset
  
  p <- ggplot(d, aes(season, slope_plot, colour = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin = lwr_plot, ymax = upr_plot),
                    position = position_dodge(0.6),
                    fatten   = 1.2, linewidth = 1) +
    geom_text(
      data = d |> filter(metric == "mean"),
      aes(x = season, label = paste0("n = ", n)),
      y = y_lab,
      size = 3.5,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = c(mean = "blue", max = "red"),
      name   = "",
      labels = c(mean = expression(step[mean]),
                 max  = expression(step[max]))
    ) +
    scale_x_discrete(limits = seas_vec, labels = seas_vec) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = NULL,
      y = "Effect size (%)",
      title = unname(var_titles[[var_name]])
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}


p_divers      <- plot_var("divers_h",    mean_models, max_models, n_df, seas_vec, y_limits = c(-1.3, 0.5))
p_temperature <- plot_var("temperature", mean_models, max_models, n_df, seas_vec)
p_moon        <- plot_var("moon",        mean_models, max_models, n_df, seas_vec)
p_cloud       <- plot_var("cloud_cover", mean_models, max_models, n_df, seas_vec)
p_wave        <- plot_var("wave_height", mean_models, max_models, n_df, seas_vec)

print(p_divers)
print(p_temperature)
print(p_moon)
print(p_cloud)
print(p_wave)



# ==============================================================================
# 7. Diver effects per season and shift (FIGURE 8)
# ==============================================================================
# 7.1. Prepare data and create models
# ==============================================================================

shift_levs <- c("8-11h", "11-13h", "13-18h", "19-23h")
shift_labs <- c("8-11","11-13","13-18","19-23")


# Models
winter_mean <- glmmTMB(step_mean ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                       data = winter, family = Gamma(link = "log"))
spring_mean <- glmmTMB(step_mean ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                       data = spring, family = Gamma(link = "log"))
summer_mean <- glmmTMB(step_mean ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                       data = summer, family = Gamma(link = "log"))
autumn_mean <- glmmTMB(step_mean ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                       data = autumn, family = Gamma(link = "log"))

winter_max <- glmmTMB(step_max ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                      data = winter, family = Gamma(link = "log"))
spring_max <- glmmTMB(step_max ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                      data = spring, family = Gamma(link = "log"))
summer_max <- glmmTMB(step_max ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                      data = summer, family = Gamma(link = "log"))
autumn_max <- glmmTMB(step_max ~ divers_h*shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover + (1 | fish_id),
                      data = autumn, family = Gamma(link = "log"))

# Pack for looping
mean_models <- list(Winter = winter_mean, Spring = spring_mean, Summer = summer_mean, Autumn = autumn_mean)
max_models  <- list(Winter = winter_max,  Spring = spring_max,  Summer = summer_max,  Autumn = autumn_max)
season_data <- list(Winter = winter,      Spring = spring,      Summer = summer,      Autumn = autumn)


# 7.2. Extract shift-specific slopes (divers_h)
# ==============================================================================
get_shift_slopes <- function(model, season_name, metric, shifts) {
  
  emtrends(model, specs = ~ shift, var = "divers_h", type = "link") |>
    as.data.frame() |>
    rename(slope = divers_h.trend, lwr = asymp.LCL, upr = asymp.UCL) |>
    mutate(
      season = season_name,
      metric = metric,
      shift  = factor(shift, levels = shifts),
      # % change per +1 diver
      slope_plot = 100 * (exp(slope) - 1),
      lwr_plot   = 100 * (exp(lwr)   - 1),
      upr_plot   = 100 * (exp(upr)   - 1)
    )
}

effects_shift <- bind_rows(
  imap_dfr(mean_models, ~ get_shift_slopes(.x, .y, "mean", shift_levs)),
  imap_dfr(max_models,  ~ get_shift_slopes(.x, .y, "max",  shift_levs))
) |>
  mutate(metric = factor(metric, levels = c("mean","max")))


# 7.3. Calculate days with divers > 0
# ==============================================================================
get_shift_n <- function(dat, season_name, shifts) {
  dat |>
    filter(divers_h > 0) |>
    mutate(date = as.Date(date_time)) |> 
    distinct(shift, date) |>
    count(shift, name = "n") |>
    complete(shift = shifts, fill = list(n = 0)) |>
    mutate(
      season = season_name,
      shift  = factor(shift, levels = shifts)
    )
}

n_shift <- imap_dfr(season_data, ~ get_shift_n(.x, .y, shift_levs))

effects_shift <- effects_shift |>
  left_join(n_shift, by = c("season","shift"))


# 7.4. Plotting
# ==============================================================================
plot_season_shift <- function(df, season_name, shift_levs, shift_labs, y_limits = NULL) {
  
  d <- df |>
    filter(season == season_name) |>
    mutate(shift = factor(shift, levels = shift_levs))
  offset <- 0.10 * diff(range(c(d$lwr_plot, d$upr_plot), na.rm = TRUE))
  y_lab  <- min(d$lwr_plot, na.rm = TRUE) - offset
  target_shift <- which(shift_levs == "19-23h")
  
  p <- ggplot(d, aes(x = shift, y = slope_plot, colour = metric)) +
    annotate("rect", 
             xmin = target_shift - 0.5, xmax = target_shift + 0.5, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey80", alpha = 0.4, color = NA) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = lwr_plot, ymax = upr_plot),
                    position = position_dodge(0.6),
                    fatten   = 1.2, linewidth = 1) +
    geom_text(
      data = d |> filter(metric == "mean"),
      aes(x = shift, label = paste0("n = ", n), colour = NULL),
      y = y_lab,
      size = 3.5
    ) +
    scale_colour_manual(
      values = c(mean = "blue", max = "red"),
      name   = "",
      labels = c(mean = expression(step[mean]),
                 max  = expression(step[max]))
    ) +
    scale_x_discrete(limits = shift_levs, labels = shift_labs) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "Shift",
      y = "Effect size (%)",
      title = season_name
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  if (!is.null(y_limits)) {
    p <- p + coord_cartesian(ylim = y_limits)
  }
  
  return(p)
}

p_winter <- plot_season_shift(effects_shift, "Winter", shift_levs, shift_labs, y_limits = c(-12.5, 18))
p_spring <- plot_season_shift(effects_shift, "Spring", shift_levs, shift_labs, y_limits = c(-4, 5))
p_summer <- plot_season_shift(effects_shift, "Summer", shift_levs, shift_labs, y_limits = c(-2, 1.3))
p_autumn <- plot_season_shift(effects_shift, "Autumn", shift_levs, shift_labs, y_limits = c(-6.35, 0.5))

combined_shifts <- (p_winter | p_spring | p_summer | p_autumn) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(combined_shifts)
