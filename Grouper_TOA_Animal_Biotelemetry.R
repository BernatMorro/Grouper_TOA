################################################################################
# TOA-BASED ACTIVITY ESTIMATION — GROUPER DATA
#
# Place this script in the same folder as tracking_detection_data.csv.gz, tracking_sensor_data.csv.gz,
# divers.csv and environmentals.csv from github.com/BernatMorro/Grouper_TOA.
# Run: Rscript Grouper_TOA_Animal_Biotelemetry.R
# Figure 3: observed monthly activity distributions.
# Figure 4: monthly daytime/night-time activity.
# Figure 5: environmental associations.
# Figure 7: diver associations by season.
# Figure 8: diver associations by season and diving shift.
#
################################################################################

# ==============================================================================
# 0. LIBRARIES AND INPUT DATA
# ==============================================================================
# Clean your workspace:remove(list=ls())

# Use the script's folder. 
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (length(script_path) == 1L && nzchar(script_path)) setwd(dirname(script_path))
}

# Install activityTOA once if needed: remotes::install_github("aspillaga/activityTOA")
library(data.table)
library(lubridate)

run_args <- commandArgs(trailingOnly = TRUE)
site_tz <- "Europe/Madrid"

# Keep absolute instants in UTC.
as_utc <- function(x) {
  if (inherits(x, "POSIXt")) {
    stamp <- with_tz(x, "UTC")
  } else if (is.character(x)) {
    stamp <- ymd_hms(x, tz = "UTC", quiet = TRUE)
  } else {
    stop("Expected ISO date-time strings or parsed POSIX timestamps.")
  }
  if (anyNA(stamp)) stop("Missing or unrecognised timestamp in tracking data.")
  stamp
}

# Raw acoustic detections to relative radial distances
# -------------------------------------------------------------------------------
read_tracking <- function() {
  raw <- fread("tracking_detection_data.csv.gz", tz = "UTC",
               select = c("date_time", "fish_id", "station", "deployment_id"))
  sensor <- fread("tracking_sensor_data.csv.gz", tz = "UTC",
                  select = c("date_time", "deployment_id", "temperature"))
  raw[, date_time := as_utc(date_time)]
  sensor[, date_time := as_utc(date_time)]
  start <- ymd("2021-02-13", tz = site_tz)
  end <- ymd("2023-12-31", tz = site_tz)
  raw <- raw[date_time >= start & date_time <= end]
  if (!nrow(raw)) stop("No detections in the study window.")
  if (anyNA(raw[, .(fish_id, station, deployment_id)]) || anyNA(sensor$deployment_id)) {
    stop("Missing fish, station or receiver deployment identifier.")
  }
  if (!is.numeric(sensor$temperature)) stop("Receiver temperature must be numeric.")
  if (!any(as.numeric(raw$date_time) %% 1 > 0)) {
    stop("All detection timestamps are whole seconds. TOA requires fractional seconds.")
  }
  sensor <- unique(sensor[, .(deployment_id, sensor_time = date_time, temperature)])
  if (anyDuplicated(sensor, by = c("deployment_id", "sensor_time"))) {
    stop("Conflicting temperature records for one receiver deployment and time.")
  }
  raw[, sensor_time := floor_date(date_time, "10 minutes")]
  raw[sensor, on = .(deployment_id, sensor_time), temperature := i.temperature]
  dropped <- sum(!is.finite(raw$temperature))
  message("Detections in the study window: ", nrow(raw),
          ". Excluded for missing or non-finite receiver temperature: ", dropped)
  raw <- raw[is.finite(temperature)]
  if (!nrow(raw)) stop("No detections have matched receiver temperatures.")
  
  # Lotek granularity is 0.5 s.
  # Method: Aspillaga et al. (2025), doi:10.1111/2041-210X.70174.
  relative_distance <- function(time, temp, group) {
    seconds <- as.numeric(time)
    if (length(unique(seconds)) < 4L) {
      stop("Too few distinct timestamps for a drift spline in ", group, ".")
    }
    frac <- activityTOA::fractSec(time)
    drift <- activityTOA::unwrapDrift(frac, time, 0.5)
    n_days <- length(unique(as.Date(time, tz = "UTC")))
    spline_fit <- stats::smooth.spline(x = seconds, y = drift, df = ceiling(n_days * 0.5))
    toa <- drift - stats::predict(spline_fit, x = seconds)$y
    distance <- toa * activityTOA::soundSpeed(temp, sal = 38, depth = 20)
    if (any(!is.finite(distance))) stop("Non-finite relative distance in ", group, ".")
    distance - min(distance)
  }
  setorder(raw, fish_id, deployment_id, date_time)
  raw[, dist := relative_distance(date_time, temperature,
                                  paste(.BY$fish_id, .BY$deployment_id, sep = "/")),
      by = .(fish_id, deployment_id)]
  # One row per occupied five-minute bin at one receiver. Steps are differences
  # of successive BIN MEANS. Each step belongs to its later bin.
  bins <- raw[, .(n_det = .N, dist = mean(dist), temperature = mean(temperature)),
              by = .(date_time = floor_date(date_time, "5 minutes"),
                     fish_id, station, deployment_id)]
  setorder(bins, fish_id, station, deployment_id, date_time)
  bin_steps <- function(distance, time, max_h = 12) {
    gap_h <- c(NA_real_, diff(as.numeric(time)) / 3600)
    step <- c(NA_real_, abs(diff(distance)))
    step[is.na(gap_h) | gap_h > max_h] <- NA_real_
    step
  }
  bins[, step := bin_steps(dist, date_time), by = .(fish_id, station, deployment_id)]
  local_hour <- hour(with_tz(bins$date_time, site_tz))
  labels <- c("23-8h", "8-11h", "11-13h", "13-18h", "18-19h", "19-23h", "23-8h")
  bins[, shift := labels[findInterval(local_hour, c(0L, 8L, 11L, 13L, 18L, 19L, 23L))]]
  shift_date <- as.Date(bins$date_time, tz = site_tz) + as.integer(local_hour >= 23L)
  bins[, shift_id := paste(shift_date, shift)]
  setcolorder(bins, c("date_time", "fish_id", "shift", "shift_id", "station",
                      "deployment_id", "n_det", "dist", "temperature", "step"))
  setorder(bins, date_time, fish_id, station, deployment_id)
  bins
}

check_archived_bins <- function(rebuilt) {
  paths <- c("tracking_toa_5min_part1.csv", "tracking_toa_5min_part2.csv")
  if (!all(file.exists(paths))) stop("--check-tracking requires both archived five-minute CSVs.")
  reference <- rbindlist(lapply(paths, fread, tz = "UTC"), use.names = TRUE)
  keys <- c("fish_id", "station", "deployment_id", "date_time")
  ordered <- function(x) {
    x <- copy(x)
    if (!all(names(rebuilt) %in% names(x))) stop("Unexpected archived five-minute columns.")
    for (v in c("fish_id", "station", "deployment_id")) set(x, j = v, value = as.character(x[[v]]))
    x[, date_time := as.numeric(as_utc(date_time))]
    if (anyDuplicated(x, by = keys)) stop("Duplicate five-minute keys.")
    setorderv(x, keys)
    x
  }
  a <- ordered(rebuilt)
  b <- ordered(reference)
  if (nrow(a) != nrow(b)) stop("Five-minute row counts differ: rebuilt ", nrow(a),
                               ", archived ", nrow(b), ".")
  for (v in c(keys, "shift", "shift_id", "n_det")) {
    if (!isTRUE(all.equal(a[[v]], b[[v]], tolerance = 0, check.attributes = FALSE))) {
      stop("Five-minute comparison differs in ", v, ".")
    }
  }
  for (v in c("dist", "temperature", "step")) {
    if (!identical(is.na(a[[v]]), is.na(b[[v]]))) stop("Missing-value pattern differs in ", v, ".")
    keep <- !is.na(a[[v]])
    error <- if (any(keep)) max(abs(a[[v]][keep] - b[[v]][keep])) else 0
    if (!is.finite(error) || error > 1e-6) stop("Five-minute comparison differs in ", v,
                                                ". Maximum absolute difference: ", error)
  }
  message("Archive comparison passed for all ", nrow(a), " five-minute bins.")
}

detect <- read_tracking()
if ("--check-tracking" %in% run_args) check_archived_bins(detect)
dir.create("derived_data", recursive = TRUE, showWarnings = FALSE)
fwrite(detect, "derived_data/tracking_toa_5min.csv")
message("Wrote ", nrow(detect), " occupied bins from ", sum(detect$n_det),
        " detections, with ", sum(is.finite(detect$step)), " valid steps.")
if ("--tracking-only" %in% run_args) quit(save = "no", status = 0L)

# Hourly biological and environmental analysis
# -------------------------------------------------------------------------------
library(suncalc)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(emmeans)

# Assign daytime/night-time using sunrise and sunset.
getDayNight <- function(times, lat, long) {
  local_dates <- as.Date(times, tz = site_tz)
  date_range <- range(local_dates)
  dates <- seq(date_range[1] - 1, date_range[2] + 1, by = "days")
  sun <- getSunlightTimes(dates, lat = lat, lon = long, keep = c("sunrise", "sunset"),
                          tz = site_tz)
  index <- match(local_dates, sun$date)
  
  ifelse(times >= sun$sunrise[index] & times <= sun$sunset[index], "day", "night")
}

# Require a single shift label within each fish-hour.
one_value <- function(x) {
  values <- unique(x[!is.na(x)])
  if (length(values) > 1) stop("Conflicting shift labels within a fish-hour.")
  if (!length(values)) return(NA_character_)
  as.character(values[1])
}


# ==============================================================================
# 1. PREPARE HOURLY ACTIVITY AND ENVIRONMENTAL DATA
# ==============================================================================

detect[, fish_id := factor(fish_id, levels = sort(unique(fish_id)))]

# Centre distances within each fish, deployment and date.
detect[, norm_dist := dist - mean(dist, na.rm = TRUE),
       by = .(fish_id, deployment_id, local_date = as.Date(date_time, tz = site_tz))]

# Combine detections into fish-hour activity metrics.
hour_metrics <- detect[, .(
  n_det = .N,
  n_steps = sum(is.finite(step)),
  temperature = mean(temperature, na.rm = TRUE),
  shift = one_value(shift),
  shift_id = one_value(shift_id),
  dist_range = diff(range(norm_dist, na.rm = TRUE)),
  step_mean = mean(step, na.rm = TRUE),
  step_sd = sd(step, na.rm = TRUE),
  step_max = max(step, na.rm = TRUE)
), by = .(fish_id, date_time = floor_date(date_time, "hour"))]

hour_metrics[, daynight := getDayNight(date_time, lat = 39.46256, long = 2.47190)]
hour_metrics[, daynight := factor(daynight, levels = c("day", "night"))]

# Match daily environmental records to local hourly dates.
copern <- fread("environmentals.csv")
copern$day <- as.Date(copern$day)
if (anyNA(copern$day)) stop("Missing or invalid dates in copern$day.")
if (anyDuplicated(copern$day)) stop("More than one Copernicus record per day.")

hour_metrics[, local_date := as.Date(date_time, tz = site_tz)]
hour_metrics[, date_f := factor(local_date)]
hour_metrics[, year_f := factor(format(local_date, "%Y"))]
index <- match(hour_metrics$local_date, copern$day)

hour_metrics$wave_height <- copern$wave_height[index]
hour_metrics$cloud_cover <- copern$cloud_cover[index]
hour_metrics$moon <- copern$moon[index]

# Join the diver data once, after aggregation to fish-hours.
data_divers <- fread("divers.csv")
data_divers <- data_divers %>%
  mutate(shift_time = recode(as.character(hour), "8-11h" = 3, "11-13h" = 2,
                             "13-18h" = 5, "19-23h" = 4),
         divers_h = divers / shift_time)

if (anyDuplicated(data_divers$shift_id)) {
  stop("Repeated shift_id values in divers.csv; resolve before joining.")
}

hour_metrics <- left_join(
  as.data.frame(hour_metrics),
  data_divers %>% select(shift_id, divers_h),
  by = "shift_id"
)

# Unmatched records represent no divers.
hour_metrics$divers_h[is.na(hour_metrics$divers_h)] <- 0

# Restore data.table format for Figure 3.
setDT(hour_metrics)
hour_metrics <- unique(hour_metrics[!is.na(step_sd)])

if (anyDuplicated(hour_metrics[, .(fish_id, date_time)])) {
  stop("More than one hourly record per fish and timestamp.")
}

message("Hourly rows after preparation: ", nrow(hour_metrics))


# ==============================================================================
# 2. FIGURE 3: OBSERVED MONTHLY DAYTIME/NIGHT-TIME ACTIVITY
# ==============================================================================

dir <- "./Figure 3/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

month_metrics <- hour_metrics[, .(
  n_hours = .N,
  n_detect = sum(n_det),
  temperature = mean(temperature, na.rm = TRUE),
  
  step_mean_lower = quantile(step_mean, 0.025, na.rm = TRUE),
  step_mean_q25 = quantile(step_mean, 0.25, na.rm = TRUE),
  step_mean_median = quantile(step_mean, 0.5, na.rm = TRUE),
  step_mean_q75 = quantile(step_mean, 0.75, na.rm = TRUE),
  step_mean_upper = quantile(step_mean, 0.975, na.rm = TRUE),
  
  step_max_lower = quantile(step_max, 0.025, na.rm = TRUE),
  step_max_q25 = quantile(step_max, 0.25, na.rm = TRUE),
  step_max_median = quantile(step_max, 0.5, na.rm = TRUE),
  step_max_q75 = quantile(step_max, 0.75, na.rm = TRUE),
  step_max_upper = quantile(step_max, 0.975, na.rm = TRUE)
), by = .(fish_id, daynight, date = floor_date(local_date, "month"))]

month_metrics <- month_metrics[
  !(fish_id == "EM0005" & date == as.Date("2023-08-01"))
]

# Complete monthly temperature sequence so missing months appear as gaps.
range_dates <- range(month_metrics$date, na.rm = TRUE)
month_seq <- seq(range_dates[1], range_dates[2], by = "months")

overall_temp <- hour_metrics[, .(temp = mean(temperature, na.rm = TRUE)),
                             by = .(m_floor = floor_date(local_date, "month"))]
full_months <- data.table(m_floor = month_seq, x_pos = seq_along(month_seq))
overall_temp <- merge(full_months, overall_temp, by = "m_floor", all.x = TRUE)

t_min <- floor(min(overall_temp$temp, na.rm = TRUE)) - 2
t_max <- ceiling(max(overall_temp$temp, na.rm = TRUE)) + 2

col <- rep(c("#577590", "#90be6d", "#f9c74f", "#f3722c"), each = 3)
col_l <- sapply(col, function(x) colorRampPalette(c(x, "black"))(100)[30])
col_m <- sapply(col, function(x) colorRampPalette(c(x, "black"))(100)[70])

month_lab <- format(month_seq, "%b")
month_metrics$month_f <- as.numeric(factor(as.character(month_metrics$date),
                                           levels = as.character(month_seq)))
month_metrics$month_num <- month(month_metrics$date)

vars <- c("step_mean", "step_max")
ylabs <- list(expression("Hourly " ~ step[mean] ~ " (m)"),
              expression("Hourly " ~ step[max] ~ " (m)"))

for (i in seq_along(vars)) {
  file_name <- file.path(dir, paste0(vars[i], ".png"))
  xlim <- c(1, length(month_seq))
  
  png(file_name, width = 2600, height = 2600, res = 250)
  par(mfrow = c(4, 2), mar = c(3.8, 4.8, 0.5, 3.5))
  
  for (id in sort(unique(month_metrics$fish_id))) {
    indx <- month_metrics$fish_id == id & month_metrics$n_hours > 5
    
    if (sum(indx) == 0) {
      plot.new()
      next
    }
    
    mid <- month_metrics[[paste0(vars[i], "_median")]]
    w_low <- month_metrics[[paste0(vars[i], "_lower")]]
    low <- month_metrics[[paste0(vars[i], "_q25")]]
    upp <- month_metrics[[paste0(vars[i], "_q75")]]
    w_upp <- month_metrics[[paste0(vars[i], "_upper")]]
    
    ylim <- range(c(0, w_low[indx], w_upp[indx]), na.rm = TRUE)
    ylim <- ylim + c(0, 1) * diff(ylim) / 8
    
    x_center <- month_metrics$month_f[indx] +
      ifelse(month_metrics$daynight[indx] == "day", -0.18, 0.18)
    
    m_num <- month_metrics$month_num[indx]
    fill_col <- ifelse(month_metrics$daynight[indx] == "day", col[m_num], "white")
    
    # Empty activity plot establishes the main coordinates.
    plot(mid[indx] ~ x_center, type = "n", ylim = ylim, main = id,
         xlim = xlim, axes = FALSE, ann = FALSE)
    
    # Temperature background and right-hand axis.
    par(new = TRUE)
    plot(overall_temp$x_pos, overall_temp$temp, type = "l", col = "tomato", lwd = 2,
         xlim = xlim, ylim = c(t_min, t_max), axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1, col = "#d90429", col.axis = "#d90429", lwd = 0, lwd.ticks = 1)
    mtext("Temperature (°C)", side = 4, line = 2.2, col = "#d90429", cex = 0.7)
    
    # Return to activity coordinates.
    par(new = TRUE)
    plot(mid[indx] ~ x_center, type = "n", ylim = ylim, xlim = xlim,
         axes = FALSE, ann = FALSE)
    
    rect(17.8 - 0.4, par("usr")[3], 19.2 + 0.4, par("usr")[4],
         col = "gray90", border = FALSE)
    
    arrows(x0 = x_center, x1 = x_center, y0 = w_low[indx], y1 = w_upp[indx],
           code = 3, angle = 90, length = 0, col = col_l[m_num], lwd = 1.2)
    rect(x_center - 0.15, low[indx], x_center + 0.15, upp[indx],
         lwd = 1.8, col = fill_col, border = col_l[m_num])
    rect(x_center - 0.15, mid[indx], x_center + 0.15, mid[indx],
         lwd = 1.8, col = col_m[m_num], border = col_m[m_num])
    
    axis(2, las = 1, lwd = 0, lwd.ticks = 1)
    lab_at <- pretty(xlim, 15)
    lab_at <- lab_at[lab_at > 0 & lab_at <= xlim[2]]
    year_at <- which(month(month_seq) == 1)
    
    axis(1, at = lab_at, labels = month_lab[lab_at], lwd = 0, lwd.ticks = 1,
         padj = -0.5)
    axis(1, at = year_at, labels = year(month_seq)[year_at], lwd = 0,
         lwd.ticks = 1, line = 2, tcl = -1, padj = -1.5, hadj = -0.15)
    axis(1, at = seq_along(month_lab)[!seq_along(month_lab) %in% lab_at],
         labels = FALSE, lwd = 0, lwd.ticks = 1, tcl = -0.25)
    
    title(ylab = ylabs[[i]], line = 2.8)
    box()
    text(par("usr")[1], par("usr")[4], adj = c(-0.2, 1.8), labels = id, font = 2)
    
    legend("topright", legend = c("Winter", "Spring", "Summer", "Autumn"),
           fill = unique(col), border = unique(col_l), box.col = "gray60",
           inset = c(0.015, 0.03), cex = 0.8)
    legend("topleft", legend = c("Day", "Night"), fill = c("gray50", "white"),
           border = "gray50", box.col = "gray60", inset = c(0.15, 0.03), cex = 0.8)
  }
  
  dev.off()
}


# ==============================================================================
# 3. SHARED SETTINGS FOR FIGURES 4, 5, 7 AND 8
# ==============================================================================

seasons <- c("Winter", "Spring", "Summer", "Autumn")
season_levels <- c(seasons, "Overall")
shift_levels <- c("8-11h", "11-13h", "13-18h", "19-23h")
shift_labels <- c("8-11", "11-13", "13-18", "19-23")
metrics <- c(mean = "step_mean", max = "step_max")

season_colours <- c(Winter = "#577590", Spring = "#90be6d",
                    Summer = "#f9c74f", Autumn = "#f3722c")
month_colours <- setNames(rep(unname(season_colours), each = 3), month.abb)
metric_colours <- c(mean = "blue", max = "red")

nsim <- 500L
lag_hours <- c(1, 2, 3, 6, 24, 48)
model_control <- glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))

# Complete analysis dataset.
required <- c("fish_id", "date_time", "local_date", "date_f", "year_f", "daynight", "shift",
              "shift_id", "step_mean", "step_max", "n_steps",
              "temperature", "divers_h", "wave_height", "moon", "cloud_cover")
missing_columns <- setdiff(required, names(hour_metrics))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}

analysis_data <- unique(as.data.frame(hour_metrics)[required])

if (!inherits(analysis_data$date_time, "POSIXt") || anyNA(analysis_data$date_time)) {
  stop("date_time must contain valid POSIXct/POSIXlt timestamps.")
}
if (anyNA(analysis_data$fish_id)) stop("Missing fish IDs.")
if (anyDuplicated(analysis_data[c("fish_id", "date_time")])) {
  stop("More than one observation per fish and timestamp.")
}
if (anyNA(analysis_data$daynight) ||
    any(!as.character(analysis_data$daynight) %in% c("day", "night"))) {
  stop("daynight must contain only day and night.")
}

numeric_vars <- c("step_mean", "step_max", "n_steps", "temperature", "divers_h",
                  "wave_height", "moon", "cloud_cover")
if (!all(vapply(analysis_data[numeric_vars], is.numeric, logical(1)))) {
  stop("Responses and environmental predictors must be numeric.")
}

complete_rows <- Reduce(`&`, lapply(analysis_data[numeric_vars], is.finite))
message("Removing ", sum(!complete_rows), " incomplete/non-finite modelling rows.")
analysis_data <- analysis_data[complete_rows, , drop = FALSE]

if (!nrow(analysis_data)) stop("No complete modelling observations remain.")
if (any(analysis_data$step_mean <= 0 | analysis_data$step_max <= 0)) {
  stop("Gamma responses must be strictly positive.")
}
if (any(analysis_data$divers_h < 0)) stop("Negative diver values found.")

step_reference <- quantile(analysis_data$n_steps, 0.5, type = 1, names = FALSE)
analysis_data$n_steps_f <- relevel(factor(analysis_data$n_steps),
                                   ref = as.character(step_reference))

analysis_data$fish_id <- factor(analysis_data$fish_id)
analysis_data$month_f <- factor(month.abb[month(analysis_data$local_date)], levels = month.abb)
analysis_data$shift <- as.character(analysis_data$shift)

month_day <- format(analysis_data$local_date, "%m-%d")
analysis_data$season <- case_when(
  month_day >= "12-21" | month_day < "03-20" ~ "Winter",
  month_day >= "03-20" & month_day < "06-21" ~ "Spring",
  month_day >= "06-21" & month_day < "09-22" ~ "Summer",
  TRUE ~ "Autumn"
)
analysis_data$season <- factor(analysis_data$season, levels = seasons)

elapsed <- (as.numeric(analysis_data$date_time) -
              min(as.numeric(analysis_data$date_time))) / 3600
if (any(abs(elapsed - round(elapsed)) > 1e-6)) {
  stop("Timestamps do not lie on a common hourly grid.")
}

analysis_data$hour_index <- as.integer(round(elapsed))
analysis_data <- analysis_data[order(analysis_data$fish_id, analysis_data$hour_index), ]
rownames(analysis_data) <- NULL


# ==============================================================================
# 4. SHARED MODELLING AND PLOTTING FUNCTIONS
# ==============================================================================

# Set up hourly positions within fish/calendar-month blocks.
# Zero-weight padding preserves the hourly factor grid.
prepare_subset <- function(dat) {
  dat <- droplevels(as.data.frame(dat))
  if (!nrow(dat)) stop("Empty model subset.")
  if (nlevels(dat$fish_id) < 2) stop("Fewer than two fish in a model subset.")
  
  dat$calendar_month <- format(dat$local_date, "%Y-%m")
  dat$fish_month <- interaction(dat$fish_id, dat$calendar_month, drop = TRUE, sep = "_")
  
  month_start <- as.POSIXct(paste0(dat$calendar_month, "-01 00:00:00"),
                            tz = site_tz)
  local_hour <- (as.numeric(dat$date_time) - as.numeric(month_start)) / 3600
  if (any(abs(local_hour - round(local_hour)) > 1e-6)) {
    stop("Invalid hourly positions within calendar months.")
  }
  
  dat$local_hour <- as.integer(round(local_hour))
  dat$.observed <- TRUE
  dat$.weight <- 1
  
  full_hours <- seq.int(0L, max(dat$local_hour))
  missing_hours <- setdiff(full_hours, dat$local_hour)
  
  if (length(missing_hours)) {
    padding <- dat[rep(which.max(dat$local_hour), length(missing_hours)), , drop = FALSE]
    padding$local_hour <- missing_hours
    padding$.observed <- FALSE
    padding$.weight <- 0
    dat <- rbind(dat, padding)
  }
  
  dat$time_ar1 <- factor(dat$local_hour, levels = full_hours)
  dat <- dat[order(dat$fish_month, dat$local_hour), ]
  rownames(dat) <- NULL
  dat
}

# Fit a single model and retain the fixed formula for predictions.
fit_model <- function(dat, response, fixed_terms, label, include_year = FALSE,
                      include_shift_id = FALSE) {
  fixed_terms <- paste(fixed_terms, "+ n_steps_f")
  
  random_terms <- c("(1 | fish_id)", "(1 | date_f)")
  if (include_year) {
    random_terms <- c(random_terms, "(1 | year_f)")
  }
  if (include_shift_id) {
    random_terms <- c(random_terms, "(1 | shift_id)")
  }
  random_terms <- c(random_terms, "ar1(time_ar1 + 0 | fish_month)")
  
  model_formula <- as.formula(paste(
    response, "~", fixed_terms, "+", paste(random_terms, collapse = " + ")
  ))
  
  message("Fitting ", label)
  
  mod <- glmmTMB(
    model_formula,
    data = dat,
    weights = .weight,
    family = Gamma(link = "log"),
    dispformula = ~ 1,
    na.action = na.fail,
    control = model_control
  )
  
  if (mod$fit$convergence != 0 || !isTRUE(mod$sdr$pdHess)) {
    message("Refitting ", label, " with orthogonal temperature terms")
    orth_terms <- sub("temperature + I(temperature^2)", "poly(temperature, 2)",
                      fixed_terms, fixed = TRUE)
    orth_formula <- as.formula(paste(
      response, "~", orth_terms, "+", paste(random_terms, collapse = " + ")
    ))
    mod <- update(mod, formula. = orth_formula)
  }
  
  if (mod$fit$convergence != 0 || !isTRUE(mod$sdr$pdHess)) {
    message("Model: ", label)
    message("Optimizer convergence code: ", mod$fit$convergence)
    message("Optimizer message: ", mod$fit$message)
    message("Positive-definite Hessian: ", mod$sdr$pdHess)
    print(VarCorr(mod))
    glmmTMB::diagnose(mod)
    stop("Convergence/Hessian check failed for ", label)
  }
  
  V <- as.matrix(vcov(mod)$cond)
  if (any(!is.finite(V)) || any(diag(V) <= 0)) {
    stop("Unusable fixed-effect covariance for ", label)
  }
  
  mf <- model.frame(mod)
  if (nrow(mf) != nrow(dat) ||
      !identical(levels(mf$time_ar1), levels(dat$time_ar1)) ||
      !isTRUE(all.equal(as.numeric(model.weights(mf)), dat$.weight,
                        check.attributes = FALSE))) {
    stop("Model data/hourly-grid verification failed for ", label)
  }
  
  fixed <- delete.response(terms(mod, component = "cond"))
  list(model = mod, data = dat, fixed = fixed, label = label)
}

fit_diver_model <- function(dat, response, fixed_terms, label) {
  reference <- fit_model(dat, response, fixed_terms, label,
                         include_year = TRUE, include_shift_id = TRUE)
  base <- reference$model
  dat$divers_10 <- dat$divers_h / 10
  model_formula <- update(formula(base), . ~ . + (0 + divers_10 | fish_id))
  environment(model_formula) <- environment(fit_model)
  beta0 <- fixef(base)$cond
  eta0 <- as.numeric(model.matrix(base) %*% beta0)
  theta0 <- getME(base, "theta")
  disp_name <- intersect(c("betadisp", "betad"), names(base$fit$par))
  stopifnot(length(disp_name) == 1L)
  bfgs <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"),
                         optCtrl = list(maxit = 3000L, reltol = 1e-12))
  
  for (attempt in 1:3) {
    if (attempt == 1L) {
      control <- model_control
      beta_start <- beta0
    } else if (attempt == 2L) {
      basis <- poly(dat$temperature, 2) * sqrt(nrow(dat) - 1)
      dat$temp_orth1 <- basis[, 1L]
      dat$temp_orth2 <- basis[, 2L]
      terms_orth <- sub("temperature + I(temperature^2)", "temp_orth1 + temp_orth2",
                        fixed_terms, fixed = TRUE)
      fixed_formula <- as.formula(paste("~", terms_orth, "+ n_steps_f"))
      X <- model.matrix(fixed_formula, dat)
      beta_start <- drop(qr.coef(qr(X), eta0))
      stopifnot(max(abs(drop(X %*% beta_start) - eta0)) < 1e-7)
      eq <- deparse1(model_formula, width.cutoff = 500L)
      eq <- sub("temperature + I(temperature^2)", "temp_orth1 + temp_orth2", eq,
                fixed = TRUE)
      eq <- sub("poly(temperature, 2)", "temp_orth1 + temp_orth2", eq, fixed = TRUE)
      model_formula <- as.formula(eq, env = environment(fit_model))
      stopifnot("temp_orth1" %in% all.vars(model_formula))
      control <- bfgs
    }
    start_sd <- if (attempt == 3L) 0.08 else 0.035
    start <- list(beta = unname(beta_start), theta = c(unname(theta0), log(start_sd)))
    start[[disp_name]] <- unname(fixef(base)$disp)
    candidate <- tryCatch(
      glmmTMB(model_formula, data = dat, weights = .weight, family = Gamma(link = "log"),
              dispformula = ~ 1, na.action = na.fail, control = control, start = start),
      error = function(e) {
        message(label, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(candidate) || candidate$fit$convergence != 0L) next
    stopifnot(identical(names(fixef(candidate)$cond), names(beta_start)),
              identical(rownames(model.frame(candidate)), rownames(model.frame(base))),
              identical(as.numeric(weights(candidate)), as.numeric(weights(base))))
    blocks <- VarCorr(candidate)$cond
    idx <- which(vapply(blocks, function(x) "divers_10" %in% rownames(x), logical(1)))
    stopifnot(length(idx) == 1L)
    sd10 <- unname(attr(blocks[[idx]], "stddev")["divers_10"])
    lr <- 2 * (base$fit$objective - candidate$fit$objective)
    if (!is.finite(sd10) || !is.finite(lr) || lr < -0.002) next
    
    # Compare the limiting model at zero slope variance on the same observations.
    boundary <- sd10 < 1e-4 && abs(lr) < 0.002
    if (abs(lr) < 0.002 && attempt == 1L) next
    if (boundary) {
      entry <- reference
      lr <- sd10 <- 0
    } else {
      if (!isTRUE(candidate$sdr$pdHess)) next
      V <- as.matrix(vcov(candidate)$cond)
      if (any(!is.finite(V)) || any(diag(V) <= 0)) next
      entry <- list(model = candidate, data = dat,
                    fixed = delete.response(terms(candidate, component = "cond")),
                    label = label)
    }
    entry$diver_slope <- data.frame(
      n_fish = nlevels(droplevels(dat$fish_id[dat$.observed])),
      slope_fit = if (boundary) "zero-variance limit" else "estimated",
      slope_sd_per_diver_h = sd10 / 10,
      slope_LR = max(0, lr),
      slope_p = if (lr > 0) pchisq(lr, df = 1, lower.tail = FALSE) / 2 else 1
    )
    return(entry)
  }
  stop("Unable to verify the fish-specific diver slope for ", label)
}


predict_fixed <- function(entry, newdata) {
  beta <- fixef(entry$model)$cond
  V <- as.matrix(vcov(entry$model)$cond)
  X <- model.matrix(entry$fixed, newdata)
  X <- X[, names(beta), drop = FALSE]
  V <- V[names(beta), names(beta), drop = FALSE]
  
  eta <- as.numeric(X %*% beta)
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  
  data.frame(predicted = exp(eta), lower = exp(eta - qnorm(0.975) * se),
             upper = exp(eta + qnorm(0.975) * se))
}

reference_grid <- function(entry, n = 1L) {
  obs <- entry$data[entry$data$.observed, , drop = FALSE]
  vars <- all.vars(entry$fixed)
  nd <- obs[rep(1L, n), vars, drop = FALSE]
  
  for (v in vars) {
    if (is.numeric(obs[[v]])) nd[[v]] <- rep(mean(obs[[v]]), n)
    if (is.factor(obs[[v]])) {
      nd[[v]] <- factor(rep(levels(obs[[v]])[1], n), levels = levels(obs[[v]]))
    }
  }
  
  rownames(nd) <- NULL
  nd
}

linear_effect <- function(entry, variable) {
  beta <- fixef(entry$model)$cond
  V <- as.matrix(vcov(entry$model)$cond)
  b <- unname(beta[variable])
  se <- sqrt(V[variable, variable])
  
  data.frame(estimate = 100 * expm1(b),
             lower = 100 * expm1(b - qnorm(0.975) * se),
             upper = 100 * expm1(b + qnorm(0.975) * se))
}

diving_days <- function(dat) {
  length(unique(dat$local_date[dat$divers_h > 0]))
}

# Plot saving function.
save_plot <- function(p, folder, name, width = 1000, height = 1000) {
  if (!dir.exists(folder)) stop("Output directory does not exist: ", folder)
  print(p)
  ggsave(filename = file.path(folder, paste0(name, ".png")), plot = p,
         width = width, height = height, units = "px", dpi = 300, bg = "white")
}


# ==============================================================================
# 5. SHARED DIAGNOSTIC FUNCTIONS
# ==============================================================================

make_diagnostic_residuals <- function(entry) {
  mod <- entry$model
  dat <- entry$data
  obs <- dat[dat$.observed, , drop = FALSE]
  
  y <- model.response(model.frame(mod))[dat$.observed]
  mu <- as.numeric(predict(mod, type = "response", re.form = NULL))[dat$.observed]
  
  Xd <- model.matrix(formula(mod, component = "disp"), data = dat)
  bd <- fixef(mod)$disp
  shape <- exp(as.numeric(Xd[, names(bd), drop = FALSE] %*% bd))[dat$.observed]
  
  set.seed(123)
  simulated <- matrix(rgamma(length(y) * nsim, shape = rep(shape, nsim),
                             scale = rep(mu / shape, nsim)), nrow = length(y))
  
  dh <- createDHARMa(simulatedResponse = simulated, observedResponse = y,
                     fittedPredictedResponse = mu, integerResponse = FALSE, seed = 123)
  list(dh = dh, obs = obs)
}

diagnose_model <- function(entry) {
  mod <- entry$model
  result <- make_diagnostic_residuals(entry)
  dh <- result$dh
  obs <- result$obs
  
  uniformity <- testUniformity(dh, plot = FALSE)
  dispersion <- testDispersion(dh, plot = FALSE, alternative = "two.sided")
  
  adequacy <- data.frame(
    model = entry$label, n_observed = nrow(obs),
    n_years = nlevels(droplevels(obs$year_f)),
    year_sd = as.numeric(attr(VarCorr(mod)$cond$year_f, "stddev"))[1L],
    convergence = mod$fit$convergence, positive_Hessian = isTRUE(mod$sdr$pdHess),
    uniformity_p = uniformity$p.value,
    dispersion_ratio = as.numeric(dispersion$statistic), dispersion_p = dispersion$p.value
  )
  
  if (!is.null(entry$diver_slope)) {
    adequacy <- cbind(adequacy, entry$diver_slope)
    adequacy$ar1_rho <- attr(VarCorr(mod)$cond$fish_month, "correlation")[1L, 2L]
    V <- as.matrix(vcov(mod)$cond)
    fixed <- rownames(V) != "(Intercept)"
    R <- cov2cor(V[fixed, fixed, drop = FALSE])
    adequacy$divers_h_VIF <- solve(R)["divers_h", "divers_h"]
  }
  r <- dh$scaledResiduals
  
  fish_results <- bind_rows(lapply(unique(as.character(obs$fish_id)), function(id) {
    selected <- as.character(obs$fish_id) == id
    times <- obs$hour_index[selected]
    residual <- r[selected]
    block <- as.character(obs$fish_month[selected])
    
    bind_rows(lapply(lag_hours, function(h) {
      later <- match(times + h, times)
      keep <- !is.na(later)
      positions <- which(keep)
      keep[positions] <- block[positions] == block[later[positions]]
      
      x <- residual[keep]
      z <- residual[later[keep]]
      correlation <- NA_real_
      
      if (length(x) >= 30 && is.finite(sd(x)) && is.finite(sd(z)) &&
          sd(x) > 0 && sd(z) > 0) {
        correlation <- cor(x, z)
      }
      
      data.frame(fish_id = id, lag_h = h, n_pairs = length(x), residual_cor = correlation)
    }))
  }))
  
  temporal <- fish_results %>%
    group_by(lag_h) %>%
    summarise(
      n_fish = sum(is.finite(residual_cor)),
      n_pairs = sum(n_pairs),
      median_cor = if (any(is.finite(residual_cor))) {
        median(residual_cor, na.rm = TRUE)
      } else NA_real_,
      min_cor = if (any(is.finite(residual_cor))) {
        min(residual_cor, na.rm = TRUE)
      } else NA_real_,
      max_cor = if (any(is.finite(residual_cor))) {
        max(residual_cor, na.rm = TRUE)
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(model = entry$label, .before = 1)
  
  list(adequacy = adequacy, temporal = temporal)
}

# Write diagnostic summaries.
run_diagnostics <- function(models, folder) {
  if (!dir.exists(folder)) stop("Output directory does not exist: ", folder)
  
  results <- lapply(models, function(entry) {
    message("Diagnosing ", entry$label)
    diagnose_model(entry)
  })
  
  adequacy <- bind_rows(lapply(results, function(x) x$adequacy))
  temporal <- bind_rows(lapply(results, function(x) x$temporal))
  pooled <- startsWith(adequacy$model, "Overall_")
  if ("slope_p" %in% names(adequacy) && any(pooled)) {
    adequacy$slope_p_holm_pooled <- NA_real_
    adequacy$slope_p_holm_pooled[pooled] <- p.adjust(adequacy$slope_p[pooled], "holm")
  }
  
  write.csv(adequacy, file.path(folder, "adequacy_summary.csv"), row.names = FALSE)
  write.csv(temporal, file.path(folder, "temporal_correlation_summary.csv"), row.names = FALSE)
  
  print(adequacy)
  print(temporal)
  invisible(list(adequacy = adequacy, temporal = temporal))
}


# ==============================================================================
# 6. FIGURE 4: MONTHLY DAYTIME AND NIGHT-TIME ACTIVITY
# ==============================================================================
# Four models.
# Monthly fixed-effect predictions in meters at common covariate values.
# ==============================================================================

dir <- "./Figure 4/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
fig4_dir <- dir

data4 <- setNames(lapply(c("day", "night"), function(period) {
  prepare_subset(analysis_data[analysis_data$daynight == period, ])
}), c("day", "night"))

# Compare day and night at the same supported count, nearest the pooled median.
shared_steps4 <- intersect(data4$day$n_steps, data4$night$n_steps)
if (!length(shared_steps4)) stop("Day and night have no step count in common.")
steps4 <- shared_steps4[which.min(abs(shared_steps4 - step_reference))]

models4 <- list()
predictions4 <- list()

environment_vars <- c("temperature", "divers_h", "wave_height", "moon", "cloud_cover")
pooled_means <- vapply(analysis_data[environment_vars], mean, numeric(1))
terms4 <- "month_f + temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover"

for (period in c("day", "night")) {
  if (!identical(levels(data4[[period]]$month_f), month.abb)) {
    stop("Figure 4 requires observations in all 12 months for each diel period.")
  }
  
  for (metric in names(metrics)) {
    label <- paste(period, metric, sep = "_")
    entry <- fit_model(data4[[period]], metrics[[metric]], terms4, label)
    models4[[label]] <- entry
    
    nd <- reference_grid(entry, 12)
    nd$n_steps_f[] <- as.character(steps4)
    nd$month_f <- factor(month.abb, levels = month.abb)
    for (v in environment_vars) nd[[v]] <- pooled_means[[v]]
    
    predictions4[[label]] <- cbind(
      data.frame(month = factor(month.abb, levels = month.abb),
                 period = period, metric = metric),
      predict_fixed(entry, nd)
    )
  }
}

predictions4 <- bind_rows(predictions4) %>%
  mutate(period = factor(period, levels = c("day", "night")))

plot_monthly <- function(metric_name, y_label) {
  dat <- filter(predictions4, metric == metric_name)
  
  ggplot(dat, aes(month, predicted, group = period)) +
    geom_line(aes(linetype = period), colour = "gray30", alpha = 0.6,
              position = position_dodge(0.5)) +
    geom_pointrange(aes(ymin = lower, ymax = upper, colour = month, shape = period),
                    position = position_dodge(0.5), linewidth = 0.8) +
    scale_colour_manual(values = month_colours, guide = "none") +
    scale_shape_manual(values = c(day = 16, night = 1), name = "Period") +
    scale_linetype_manual(values = c(day = "solid", night = "dotted"), name = "Period") +
    labs(x = NULL, y = y_label) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

p_month_mean <- plot_monthly("mean", expression("Predicted" ~ step[mean] ~ "(m)"))
p_month_max <- plot_monthly("max", expression("Predicted" ~ step[max] ~ "(m)"))
figure4 <- (p_month_mean / p_month_max) + plot_layout(guides = "collect")
figure4 <- figure4 & theme(legend.position = "bottom")

save_plot(figure4, fig4_dir, "Figure4", 2000, 1600)

figure4_manuscript <- figure4 + plot_annotation(tag_levels = "A")
ggsave(file.path(fig4_dir, "Figure4_manuscript.png"), plot = figure4_manuscript,
       width = 1800, height = 1849, units = "px", dpi = 300, bg = "white")


for (panel in c("a", "b")) {
  current_plot <- if (panel == "a") p_month_mean else p_month_max
  ggsave(filename = file.path(fig4_dir, paste0("Fig. 4", panel, "2.png")),
         plot = current_plot, device = "png",
         width = 2000, height = 1200, units = "px", dpi = 250, bg = "white")
}

diagnostics4 <- run_diagnostics(models4, fig4_dir)


# ==============================================================================
# 7. FIGURE 5: ENVIRONMENTAL ASSOCIATIONS
# ==============================================================================
# Twenty models: two metrics x two diel periods x five seasonal subsets.
# ==============================================================================

dir <- "./Figure 5/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
fig5_dir <- dir

models5 <- list()
linear5 <- list()
temperature5 <- list()
terms5 <- "temperature + I(temperature^2) + divers_h + wave_height + moon + cloud_cover"

for (period in c("day", "night")) {
  for (s in season_levels) {
    obs <- analysis_data[analysis_data$daynight == period, ]
    if (s != "Overall") obs <- obs[as.character(obs$season) == s, ]
    dat <- prepare_subset(obs)
    
    for (metric in names(metrics)) {
      label <- paste(period, s, metric, sep = "_")
      fixed_terms <- if (s == "Overall") paste(terms5, "+ season") else terms5
      entry <- fit_model(dat, metrics[[metric]], fixed_terms, label)
      models5[[label]] <- entry
      
      for (v in c("wave_height", "cloud_cover", "moon")) {
        key <- paste(label, v, sep = "_")
        linear5[[key]] <- cbind(
          data.frame(period = period, season = s, metric = metric, variable = v),
          linear_effect(entry, v)
        )
      }
      
      if (s != "Overall") {
        temperatures <- sort(unique(obs$temperature))
        nd <- reference_grid(entry, length(temperatures))
        nd$temperature <- temperatures
        # Each percentage uses the arithmetic mean prediction across this grid.
        beta <- fixef(entry$model)$cond
        X <- model.matrix(entry$fixed, data = nd)[, names(beta), drop = FALSE]
        V <- as.matrix(vcov(entry$model)$cond)[names(beta), names(beta), drop = FALSE]
        predicted <- as.numeric(exp(X %*% beta))
        relative <- predicted / mean(predicted)
        
        # The numerator and the fitted reference share the same coefficients.
        # The rows of gradient are derivatives of log(relative) with respect to beta.
        weights <- predicted / sum(predicted)
        mean_X <- drop(crossprod(weights, X))
        gradient <- sweep(X, 2L, mean_X, "-")
        se_log_relative <- sqrt(pmax(0, rowSums((gradient %*% V) * gradient)))
        half_width <- qnorm(0.975) * se_log_relative
        
        curve <- data.frame(
          temperature = temperatures, predicted = predicted,
          estimate = 100 * (relative - 1),
          lower = 100 * (relative * exp(-half_width) - 1),
          upper = 100 * (relative * exp(half_width) - 1)
        )
        curve$period <- period
        curve$season <- s
        curve$metric <- metric
        
        if (s == "Spring") curve <- curve[curve$temperature <= 24, ]
        temperature5[[label]] <- curve
      }
    }
  }
}

linear5 <- bind_rows(linear5) %>%
  mutate(season = factor(season, levels = season_levels),
         metric = factor(metric, levels = c("mean", "max")))

temperature5 <- bind_rows(temperature5) %>%
  mutate(season = factor(season, levels = seasons))

theme_environment <- function(period) {
  grid_colour <- if (period == "night") "white" else "grey84"
  background <- if (period == "night") "grey89" else "white"
  
  theme_minimal() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = background, colour = NA),
      panel.grid.major = element_line(colour = grid_colour, linewidth = 0.5),
      panel.grid.minor = element_line(colour = grid_colour, linewidth = 0.3)
    )
}

plot_temperature <- function(period_name, metric_name) {
  dat <- filter(temperature5, period == period_name, metric == metric_name)
  
  ggplot(dat, aes(temperature, estimate, colour = season, fill = season)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(name = "Season", values = season_colours, breaks = seasons) +
    scale_fill_manual(name = "Season", values = season_colours, breaks = seasons) +
    labs(x = "Temperature (°C)", y = "Effect size (%)") +
    theme_environment(period_name)
}

plot_environment <- function(variable_name, period_name) {
  dat <- filter(linear5, variable == variable_name, period == period_name)
  
  ggplot(dat, aes(season, estimate, colour = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.5),
                    linewidth = 0.8) +
    scale_colour_manual(name = "Step length", values = metric_colours,
                        breaks = c("mean", "max"),
                        labels = expression(step[mean], step[max])) +
    scale_x_discrete(limits = season_levels) +
    labs(x = NULL, y = "Effect size (%)") +
    theme_environment(period_name) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

panels5 <- list(
  temp_mean_day = plot_temperature("day", "mean"),
  temp_mean_night = plot_temperature("night", "mean"),
  temp_max_day = plot_temperature("day", "max"),
  temp_max_night = plot_temperature("night", "max"),
  wave_day = plot_environment("wave_height", "day"),
  wave_night = plot_environment("wave_height", "night"),
  cloud_day = plot_environment("cloud_cover", "day"),
  cloud_night = plot_environment("cloud_cover", "night"),
  moon_day = plot_environment("moon", "day"),
  moon_night = plot_environment("moon", "night")
)

combined_panels <- panels5
row_titles <- list(expression("Temperature" ~ "(" * step[mean] * ")"),
                   expression("Temperature" ~ "(" * step[max] * ")"),
                   "Wave height", "Cloud cover", "Moon")
for (i in seq_along(combined_panels)) {
  panel <- combined_panels[[i]]
  row <- (i + 1L) %/% 2L
  if (i %% 2L == 1L) panel <- panel + labs(title = row_titles[[row]])
  if (i %% 2L == 0L) panel <- panel + labs(y = NULL)
  if (i <= 2L) panel <- panel + labs(subtitle = c("Day", "Night")[i], x = NULL)
  if (row %in% c(3L, 4L)) {
    panel <- panel + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (!i %in% c(1L, 5L)) panel <- panel + guides(colour = "none", fill = "none")
  combined_panels[[i]] <- panel
}

figure5 <- wrap_plots(combined_panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
figure5 <- figure5 & theme(
  legend.position = "right", legend.text = element_text(size = 12),
  plot.title = element_text(size = 12),
  plot.subtitle = element_text(size = 10),
  plot.tag = element_text(size = 13, face = "bold")
)
save_plot(figure5, fig5_dir, "Figure5", 2600, 2700)

for (name in names(panels5)) {
  ggsave(filename = file.path(fig5_dir, paste0(name, ".png")),
         plot = panels5[[name]], device = "png",
         width = 1000, height = 500, units = "px", dpi = 300, bg = "white")
}

diagnostics5 <- run_diagnostics(models5, fig5_dir)

# ==============================================================================
# 8. SHARED DIVING-HOUR DATA FOR FIGURES 7 AND 8
# ==============================================================================

diver_data <- analysis_data[analysis_data$shift %in% shift_levels, ]
diver_data <- subset(diver_data, !(season == "Winter" & shift == "19-23h"))
diver_data$shift <- factor(diver_data$shift, levels = shift_levels)

season_diver_data <- setNames(lapply(seasons, function(s) {
  dat <- diver_data[as.character(diver_data$season) == s, ]
  dat <- droplevels(dat)
  
  if (!nrow(dat) || nlevels(dat$shift) < 2) stop("Insufficient shifts for ", s)
  if (length(unique(dat$divers_h)) < 2) stop("No diver-count variation for ", s)
  dat
}), seasons)

observed7 <- c(season_diver_data, list(Overall = diver_data))
data7 <- lapply(observed7, prepare_subset)


# ==============================================================================
# 9. FIGURE 7: DIVER ASSOCIATIONS BY SEASON
# ==============================================================================

dir <- "./Figure 7/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
fig7_dir <- dir

models7 <- list()
effects7 <- list()
terms7 <- "divers_h + moon + temperature + I(temperature^2) + wave_height + cloud_cover + shift"

for (s in season_levels) {
  for (metric in names(metrics)) {
    label <- paste(s, metric, sep = "_")
    fixed_terms <- if (s == "Overall") paste(terms7, "+ season") else terms7
    entry <- fit_diver_model(data7[[s]], metrics[[metric]], fixed_terms, label)
    models7[[label]] <- entry
    
    effects7[[label]] <- cbind(
      data.frame(season = s, metric = metric, dd = diving_days(observed7[[s]])),
      linear_effect(entry, "divers_h")
    )
  }
}

effects7 <- bind_rows(effects7) %>%
  mutate(season = factor(season, levels = season_levels),
         metric = factor(metric, levels = c("mean", "max")))

range7 <- range(c(0, effects7$lower, effects7$upper))
span7 <- max(diff(range7), 0.01)
labels7 <- filter(effects7, metric == "mean") %>%
  mutate(label_y = range7[1] - 0.12 * span7)

figure7 <- ggplot(effects7, aes(season, estimate, colour = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6),
                  linewidth = 0.8) +
  geom_text(data = labels7, aes(season, label_y, label = paste0("dd = ", dd)),
            inherit.aes = FALSE, size = 2.5) +
  scale_colour_manual(name = "Step length", values = metric_colours,
                      breaks = c("mean", "max"),
                      labels = expression(step[mean], step[max])) +
  scale_x_discrete(limits = season_levels) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  labs(x = NULL, y = "Effect per additional diver/h (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", axis.title.y = element_text(size = 9))

ggsave(filename = file.path(fig7_dir, "divers.png"), plot = figure7, device = "png",
       width = 1000, height = 500, units = "px", dpi = 300, bg = "white")

figure7_manuscript <- figure7 +
  labs(title = "Divers") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())
ggsave(file.path(fig7_dir, "Figure7_manuscript.png"), plot = figure7_manuscript,
       width = 1350, height = 1509, units = "px", dpi = 350, bg = "white")

diagnostics7 <- run_diagnostics(models7, fig7_dir)


# ==============================================================================
# 10. FIGURE 8: SHIFT-SPECIFIC DIVER ASSOCIATIONS BY SEASON
# ==============================================================================

dir <- "./Figure 8/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
fig8_dir <- dir

models8 <- list()
effects8 <- list()
terms8 <- "divers_h * shift + moon + temperature + I(temperature^2) + wave_height + cloud_cover"

for (s in seasons) {
  obs <- season_diver_data[[s]]
  dat <- data7[[s]]
  
  for (sh in levels(obs$shift)) {
    values <- obs$divers_h[obs$shift == sh]
    if (length(unique(values)) < 2) {
      stop("No diver-count variation for ", s, ", ", sh)
    }
  }
  
  for (metric in names(metrics)) {
    label <- paste(s, metric, sep = "_")
    entry <- fit_diver_model(dat, metrics[[metric]], terms8, label)
    models8[[label]] <- entry
    
    observed_model_data <- entry$data[entry$data$.observed, , drop = FALSE]
    trends <- emtrends(entry$model, specs = ~ shift, var = "divers_h",
                       data = observed_model_data, type = "link")
    tr <- as.data.frame(summary(trends, infer = c(TRUE, FALSE),
                                level = 0.95, adjust = "none"))
    
    lower_name <- intersect(c("asymp.LCL", "lower.CL"), names(tr))
    upper_name <- intersect(c("asymp.UCL", "upper.CL"), names(tr))
    
    if (length(lower_name) != 1 || length(upper_name) != 1) {
      stop("Unexpected confidence interval columns for ", label)
    }
    
    effects8[[label]] <- data.frame(
      season = s, metric = metric, shift = as.character(tr$shift),
      estimate = 100 * expm1(tr$divers_h.trend),
      lower = 100 * expm1(tr[[lower_name]]),
      upper = 100 * expm1(tr[[upper_name]])
    )
  }
}

effects8 <- bind_rows(effects8) %>%
  mutate(shift = factor(shift, levels = shift_levels),
         metric = factor(metric, levels = c("mean", "max")))

dd8 <- bind_rows(lapply(seasons, function(s) {
  obs <- season_diver_data[[s]]
  bind_rows(lapply(shift_levels, function(sh) {
    data.frame(season = s, shift = sh,
               dd = diving_days(obs[obs$shift == sh, , drop = FALSE]))
  }))
})) %>%
  mutate(shift = factor(shift, levels = shift_levels))

plot_shift <- function(season_name, dd_size = 3, base_size = 11) {
  dat <- filter(effects8, season == season_name)
  effect_range <- range(c(0, dat$lower, dat$upper))
  span <- max(diff(effect_range), 0.01)
  labels8 <- filter(dd8, season == season_name) %>%
    mutate(label_y = effect_range[1] - 0.12 * span)
  night_position <- match("19-23h", shift_levels)
  
  ggplot(dat, aes(shift, estimate, colour = metric)) +
    annotate("rect", xmin = night_position - 0.5, xmax = night_position + 0.5,
             ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6),
                    linewidth = 0.8) +
    geom_text(data = labels8, aes(shift, label_y, label = paste0("dd= ", dd)),
              inherit.aes = FALSE, size = dd_size) +
    scale_colour_manual(name = "Step length", values = metric_colours,
                        breaks = c("mean", "max"),
                        labels = expression(step[mean], step[max])) +
    scale_x_discrete(limits = shift_levels, labels = shift_labels, drop = FALSE) +
    labs(title = season_name, x = "Shift (hours)", y = "Effect per additional diver/h (%)") +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom")
}

panels8 <- setNames(lapply(seasons, plot_shift), tolower(seasons))
figure8 <- wrap_plots(panels8, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")
figure8 <- figure8 & theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))

for (season_name in seasons) {
  current_plot <- plot_shift(season_name, dd_size = 1.75, base_size = 9) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = grid::unit(10, "pt"),
          axis.text.x = element_text(size = 7))
  ggsave(filename = file.path(fig8_dir, paste0(tolower(season_name), ".png")),
         plot = current_plot, device = "png",
         width = 500, height = 1000, units = "px", dpi = 300, bg = "white")
}

save_plot(figure8, fig8_dir, "Figure8_combined_AR1", 2600, 1900)
manuscript_panels8 <- lapply(seasons, function(s) {
  p <- plot_shift(s, dd_size = 2.0, base_size = 9.5) +
    labs(x = "Shift (h)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 8))
  if (s != seasons[1L]) p <- p + labs(y = NULL)
  p
})
figure8_manuscript <- wrap_plots(manuscript_panels8, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A")
figure8_manuscript <- figure8_manuscript &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11))
ggsave(file.path(fig8_dir, "Figure8_manuscript.png"), plot = figure8_manuscript,
       width = 2900, height = 1350, units = "px", dpi = 350, bg = "white")

diagnostics8 <- run_diagnostics(models8, fig8_dir)

message("Figures and diagnostic summaries completed.")

# ==============================================================================
# 11. SUPPLEMENTARY MODEL DIAGNOSTIC FIGURES
#
# Generates one six-panel diagnostic figure per model using the fitted models
# from Sections 6, 7, 9 and 10. 
# ==============================================================================

draw_diagnostic_composition <- function(result, model_name) {
  dh <- result$dh
  obs <- result$obs
  
  old_par <- par(mfrow = c(2, 3), mar = c(4.5, 4.2, 4, 1), oma = c(0, 0, 2, 0))
  on.exit(par(old_par))
  
  plot_groups <- function(group, heading, xlab) {
    old_labels <- par(col.lab = "transparent")
    on.exit(par(old_labels))
    DHARMa::testCategorical(dh, catPred = group)
    mtext(heading, side = 3, line = 2.4, font = 2)
    mtext(xlab, side = 1, line = 3)
    mtext("DHARMa residual", side = 2, line = 3)
  }
  
  plot_predictor <- function(predictor, heading, xlab) {
    DHARMa::plotResiduals(dh, form = predictor, quantreg = FALSE,
                          rank = FALSE, asFactor = FALSE, smoothScatter = FALSE,
                          main = "", xlab = xlab)
    title(main = heading)
  }
  
  DHARMa::plotQQunif(dh, testUniformity = TRUE,
                     testDispersion = TRUE, testOutliers = FALSE)
  
  plot_groups(droplevels(obs$fish_id), "Residuals by fish", "Fish identity")
  
  if (length(unique(obs$divers_h)) > 1L) {
    plot_predictor(obs$divers_h, "Residuals vs diver abundance", "Divers per hour")
  } else {
    plot.new()
    title(main = "Residuals vs diver abundance")
    text(0.5, 0.5, "Constant diver abundance")
  }
  
  plot_predictor(obs$temperature, "Residuals vs temperature", "Temperature (°C)")
  
  month_group <- droplevels(factor(format(obs$local_date, "%m"),
                                   levels = sprintf("%02d", 1:12), labels = month.abb))
  if (nlevels(month_group) > 1L) {
    plot_groups(month_group, "Residuals by month", "Month")
  } else {
    hist(dh$scaledResiduals, breaks = seq(0, 1, length.out = 21),
         probability = TRUE, col = "grey85", border = "white",
         main = "Residuals by month", xlab = "Scaled residual")
    abline(h = 1, col = "red", lty = 2)
  }
  
  plot_predictor(obs$hour_index, "Residuals vs elapsed time", "Elapsed hours")
  
  mtext(model_name, outer = TRUE, side = 3, line = 0.4, font = 2)
}

save_diagnostic_compositions <- function(models, folder, figure_name) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  
  for (model_name in names(models)) {
    message("Saving diagnostics: ", figure_name, " — ", model_name)
    result <- make_diagnostic_residuals(models[[model_name]])
    output_file <- file.path(folder, paste0(model_name, "_DHARMa_diagnostics.png"))
    
    png(output_file, width = 2700, height = 1600, res = 180)
    tryCatch(
      draw_diagnostic_composition(result, paste(figure_name, model_name, sep = " — ")),
      finally = dev.off()
    )
  }
}

diagnostic_models <- list("Figure 4" = models4, "Figure 5" = models5,
                          "Figure 7" = models7, "Figure 8" = models8)
for (figure_name in names(diagnostic_models)) {
  save_diagnostic_compositions(diagnostic_models[[figure_name]],
                               file.path("Diagnostics", figure_name), figure_name)
}

message("All diagnostic compositions saved.")


# ==============================================================================
# 12. VARIANCE INFLATION FACTORS FOR THE CONDITIONAL FIXED EFFECTS
# ==============================================================================

model_vifs <- function(entry, figure, model_name) {
  mod <- entry$model
  X <- model.matrix(mod, component = "cond")
  V <- as.matrix(vcov(mod)$cond)
  keep <- rownames(V) != "(Intercept)"
  R <- cov2cor(V[keep, keep, drop = FALSE])
  if (any(!is.finite(R))) stop("Non-finite coefficient correlation in ", model_name)
  
  column <- match(rownames(R), colnames(X))
  assignment <- attr(X, "assign")[column]
  term_labels <- attr(terms(mod, component = "cond"), "term.labels")
  if (anyNA(assignment) || any(assignment < 1L)) {
    stop("Cannot match conditional coefficients to terms in ", model_name)
  }
  term_for_coef <- term_labels[assignment]
  
  # A term and its complement partition the coefficient correlation matrix.
  log_det <- function(M) {
    if (nrow(M) == 0L) return(0)
    2 * sum(log(diag(chol(M))))
  }
  log_det_all <- log_det(R)
  calculate <- function(indices, term, scope) {
    rest <- setdiff(seq_len(nrow(R)), indices)
    log_gvif <- log_det(R[indices, indices, drop = FALSE]) +
      log_det(R[rest, rest, drop = FALSE]) - log_det_all
    df <- length(indices)
    data.frame(figure = figure, model = model_name, scope = scope, term = term,
               df = df, GVIF = exp(log_gvif),
               GVIF_1_over_2df = exp(log_gvif / (2 * df)),
               VIF_scale = exp(log_gvif / df))
  }
  rows <- lapply(unique(term_for_coef), function(term) {
    calculate(which(term_for_coef == term), term, "fitted term")
  })
  
  # The two temperature coefficients span the same quadratic in either basis.
  temperature_terms <- c("temperature", "I(temperature^2)", "poly(temperature, 2)",
                         "temp_orth1", "temp_orth2")
  temp <- which(term_for_coef %in% temperature_terms)
  if (length(temp) == 2L) {
    rows[[length(rows) + 1L]] <- calculate(temp, "temperature (joint quadratic)", "joint term")
  }
  do.call(rbind, rows)
}

vif_models <- list("4" = models4, "5" = models5, "7" = models7, "8" = models8)
model_VIFs <- do.call(rbind, lapply(names(vif_models), function(fig) {
  entries <- vif_models[[fig]]
  do.call(rbind, lapply(names(entries), function(id) model_vifs(entries[[id]], fig, id)))
}))
rownames(model_VIFs) <- NULL
dir.create("./Diagnostics", recursive = TRUE, showWarnings = FALSE)
write.csv(model_VIFs, "./Diagnostics/Model_VIFs.csv", row.names = FALSE)



