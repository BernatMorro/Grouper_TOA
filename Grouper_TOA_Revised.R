################################################################################
# TOA-BASED ACTIVITY ESTIMATION — GROUPER DATA
#
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

# Run with the project folder as your working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(lubridate)
library(suncalc)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(emmeans)

#Load objects
copern <- fread("./environmentals.csv")
data_divers <- fread("./divers.csv")
detect_part1 <- fread("./tracking_toa_5min_part1.csv")
detect_part2 <- fread("./tracking_toa_5min_part2.csv")

# Verify that both files contain the same columns in the same order.
if (!identical(names(detect_part1), names(detect_part2))) {
  stop("The two tracking files must have identical column names and order.")
}

# Reconstruct the complete tracking dataset.
detect <- rbindlist(
  list(detect_part1, detect_part2),
  use.names = TRUE,
  fill = FALSE)

# Remove temporary objects.
rm(detect_part1, detect_part2)


# Assign daytime/night-time using sunrise and sunset.
getDayNight <- function(times, lat, long) {
  date_range <- range(as.Date(times))
  dates <- seq(date_range[1] - 1, date_range[2] + 1, by = "days")
  sun <- getSunlightTimes(dates, lat = lat, lon = long, keep = c("sunrise", "sunset"))
  index <- match(as.Date(times), sun$date)
  
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
       by = .(fish_id, deployment_id, as.Date(date_time))]

# Combine detections into fish-hour activity metrics.
hour_metrics <- detect[, .(
  n_det = .N,
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
copern$day <- as.Date(copern$day)
if (anyNA(copern$day)) stop("Missing or invalid dates in copern$day.")
if (anyDuplicated(copern$day)) stop("More than one Copernicus record per day.")

hour_date <- as.Date(hour_metrics$date_time, tz = "Europe/Madrid")
index <- match(hour_date, copern$day)

hour_metrics$wave_height <- copern$wave_height[index]
hour_metrics$cloud_cover <- copern$cloud_cover[index]
hour_metrics$moon <- copern$moon[index]

# Divers divided by shift duration as variable.
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
  n_days = .N,
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
), by = .(fish_id, daynight, date = floor_date(date_time, "month"))]

month_metrics <- month_metrics[
  !(fish_id == "EM0005" & date == as.Date("2023-08-01"))
]

# Complete monthly temperature sequence so missing months appear as gaps.
range_dates <- range(month_metrics$date, na.rm = TRUE)
month_seq <- seq(range_dates[1], range_dates[2], by = "months")

overall_temp <- hour_metrics[, .(temp = mean(temperature, na.rm = TRUE)),
                             by = .(m_floor = floor_date(date_time, "month"))]
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
ylabs <- list(expression("Daily " ~ step[mean] ~ " (m)"),
              expression("Daily " ~ step[max] ~ " (m)"))

for (i in seq_along(vars)) {
  file_name <- file.path(dir, paste0(vars[i], ".png"))
  xlim <- c(1, length(month_seq))
  
  png(file_name, width = 2600, height = 2600, res = 250)
  par(mfrow = c(4, 2), mar = c(3.8, 4.8, 0.5, 3.5))
  
  for (id in sort(unique(month_metrics$fish_id))) {
    indx <- month_metrics$fish_id == id & month_metrics$n_days > 5
    
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

# Use a separate complete analysis dataset.
required <- c("fish_id", "date_time", "daynight", "shift", "step_mean", "step_max",
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

numeric_vars <- c("step_mean", "step_max", "temperature", "divers_h",
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

analysis_data$fish_id <- factor(analysis_data$fish_id)
analysis_data$month_f <- factor(month.abb[month(analysis_data$date_time)], levels = month.abb)
analysis_data$shift <- as.character(analysis_data$shift)

month_day <- format(analysis_data$date_time, "%m-%d")
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
  
  dat$calendar_month <- format(dat$date_time, "%Y-%m")
  dat$fish_month <- interaction(dat$fish_id, dat$calendar_month, drop = TRUE, sep = "_")
  
  month_start <- as.POSIXct(paste0(dat$calendar_month, "-01 00:00:00"),
                            tz = lubridate::tz(dat$date_time))
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
fit_model <- function(dat, response, fixed_terms, label) {
  fixed_formula <- as.formula(paste("~", fixed_terms))
  model_formula <- as.formula(paste(response, "~", fixed_terms,
                                    "+ (1 | fish_id) + ar1(time_ar1 + 0 | fish_month)"))
  
  message("Fitting ", label)
  mod <- glmmTMB(model_formula, data = dat, weights = .weight,
                 family = Gamma(link = "log"), dispformula = ~ 1,
                 na.action = na.fail, control = model_control)
  
  if (mod$fit$convergence != 0 || !isTRUE(mod$sdr$pdHess)) {
    stop("Convergence/Hessian check failed for ", label)
  }
  
  mf <- model.frame(mod)
  if (nrow(mf) != nrow(dat) ||
      !identical(levels(mf$time_ar1), levels(dat$time_ar1)) ||
      !isTRUE(all.equal(as.numeric(model.weights(mf)), dat$.weight,
                        check.attributes = FALSE))) {
    stop("Model data/hourly-grid verification failed for ", label)
  }
  
  list(model = mod, data = dat, fixed = fixed_formula, label = label)
}

# Predict with fish and temporal random effects set to zero.
# Calculate pointwise 95% CIs on the log scale, then exponentiate.
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

# Construct prediction data using observed subset means.
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

# Percentage change per +1 unit of a linear predictor, with a 95% CI.
linear_effect <- function(entry, variable) {
  beta <- fixef(entry$model)$cond
  V <- as.matrix(vcov(entry$model)$cond)
  b <- unname(beta[variable])
  se <- sqrt(V[variable, variable])
  
  data.frame(estimate = 100 * expm1(b),
             lower = 100 * expm1(b - qnorm(0.975) * se),
             upper = 100 * expm1(b + qnorm(0.975) * se))
}

# Count distinct dates with recorded diving activity, once per date.
diving_days <- function(dat) {
  dates <- as.Date(dat$date_time[dat$divers_h > 0], tz = lubridate::tz(dat$date_time))
  length(unique(dates))
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

diagnose_model <- function(entry) {
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
  
  uniformity <- testUniformity(dh, plot = FALSE)
  dispersion <- testDispersion(dh, plot = FALSE, alternative = "two.sided")
  
  adequacy <- data.frame(
    model = entry$label, n_observed = length(y),
    convergence = mod$fit$convergence, positive_Hessian = isTRUE(mod$sdr$pdHess),
    uniformity_p = uniformity$p.value,
    dispersion_ratio = as.numeric(dispersion$statistic), dispersion_p = dispersion$p.value
  )
  
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
# Reference per metric = mean of the 24 monthly daytime/night-time predictions.
# ==============================================================================

dir <- "./Figure 4/"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
fig4_dir <- dir

data4 <- setNames(lapply(c("day", "night"), function(period) {
  prepare_subset(analysis_data[analysis_data$daynight == period, ])
}), c("day", "night"))

models4 <- list()
predictions4 <- list()

environment_vars <- c("temperature", "divers_h", "wave_height", "moon", "cloud_cover")
pooled_means <- vapply(analysis_data[environment_vars], mean, numeric(1))
terms4 <- "month_f + temperature + divers_h + wave_height + moon + cloud_cover"

for (period in c("day", "night")) {
  if (!identical(levels(data4[[period]]$month_f), month.abb)) {
    stop("Figure 4 requires observations in all 12 months for each diel period.")
  }
  
  for (metric in names(metrics)) {
    label <- paste(period, metric, sep = "_")
    entry <- fit_model(data4[[period]], metrics[[metric]], terms4, label)
    models4[[label]] <- entry
    
    nd <- reference_grid(entry, 12)
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
  group_by(metric) %>%
  mutate(reference = mean(predicted), estimate = 100 * (predicted / reference - 1),
         lower = 100 * (lower / reference - 1), upper = 100 * (upper / reference - 1)) %>%
  ungroup() %>%
  mutate(period = factor(period, levels = c("day", "night")))

plot_monthly <- function(metric_name) {
  dat <- filter(predictions4, metric == metric_name)
  
  ggplot(dat, aes(month, estimate, group = period)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_line(aes(linetype = period), colour = "gray30", alpha = 0.6,
              position = position_dodge(0.5)) +
    geom_pointrange(aes(ymin = lower, ymax = upper, colour = month, shape = period),
                    position = position_dodge(0.5), linewidth = 0.8) +
    scale_colour_manual(values = month_colours, guide = "none") +
    scale_shape_manual(values = c(day = 16, night = 1), name = "Period") +
    scale_linetype_manual(values = c(day = "solid", night = "dotted"), name = "Period") +
    labs(x = NULL, y = "Deviation from reference (%)") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

p4_mean <- plot_monthly("mean")
p4_max <- plot_monthly("max")
figure4 <- (p4_mean / p4_max) + plot_layout(guides = "collect")
figure4 <- figure4 & theme(legend.position = "bottom")

save_plot(figure4, fig4_dir, "Figure4", 2000, 1600)
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
        curve <- cbind(data.frame(temperature = temperatures), predict_fixed(entry, nd))
        reference <- mean(curve$predicted)
        
        curve$estimate <- 100 * (curve$predicted / reference - 1)
        curve$lower <- 100 * (curve$lower / reference - 1)
        curve$upper <- 100 * (curve$upper / reference - 1)
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
  t <- theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank())
  
  if (period == "night") {
    t <- t + theme(panel.background = element_rect(fill = "grey90", colour = NA),
                   panel.grid.major = element_line(colour = "white"))
  }
  
  t
}

plot_temperature <- function(period_name, metric_name) {
  dat <- filter(temperature5, period == period_name, metric == metric_name)
  
  ggplot(dat, aes(temperature, estimate, colour = season, fill = season)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = season_colours) +
    scale_fill_manual(values = season_colours) +
    labs(x = "Temperature (°C)", y = "Effect size (%)") +
    theme_environment(period_name)
}

plot_environment <- function(variable_name, period_name) {
  dat <- filter(linear5, variable == variable_name, period == period_name)
  
  ggplot(dat, aes(season, estimate, colour = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.5),
                    linewidth = 0.8) +
    scale_colour_manual(values = metric_colours) +
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

figure5 <- wrap_plots(panels5, ncol = 2)
save_plot(figure5, fig5_dir, "Figure5", 2000, 3000)

for (name in names(panels5)) {
  save_plot(panels5[[name]], fig5_dir, name, 1000, 600)
}

diagnostics5 <- run_diagnostics(models5, fig5_dir)


# ==============================================================================
# 8. SHARED DIVING-HOUR DATA FOR FIGURES 7 AND 8
# ==============================================================================

diver_data <- analysis_data[analysis_data$shift %in% shift_levels, ]
diver_data$shift <- factor(diver_data$shift, levels = shift_levels)

season_diver_data <- setNames(lapply(seasons, function(s) {
  dat <- diver_data[as.character(diver_data$season) == s, ]
  if (s == "Winter") dat <- dat[as.character(dat$shift) != "19-23h", ]
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
    entry <- fit_model(data7[[s]], metrics[[metric]], fixed_terms, label)
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
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6),
                  linewidth = 0.8) +
  geom_text(data = labels7, aes(season, label_y, label = paste0("dd = ", dd)),
            inherit.aes = FALSE, size = 3) +
  scale_colour_manual(values = metric_colours) +
  scale_x_discrete(limits = season_levels) +
  labs(x = NULL, y = "Effect size (%)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(figure7, fig7_dir, "divers", 1000, 1000)
diagnostics7 <- run_diagnostics(models7, fig7_dir)


# ==============================================================================
# 10. FIGURE 8: SHIFT-SPECIFIC DIVER ASSOCIATIONS BY SEASON
# ==============================================================================
# The interaction estimates separate slopes for each shift.
# Extract point estimates and pointwise 95% CIs only.
# No tests against zero, interaction tests or pairwise slope comparisons.
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
    entry <- fit_model(dat, metrics[[metric]], terms8, label)
    models8[[label]] <- entry
    
    trends <- emtrends(entry$model, specs = ~ shift, var = "divers_h",
                       data = obs, type = "link")
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
    data.frame(
      season = s, shift = sh,
      dd = diving_days(obs[as.character(obs$shift) == sh, , drop = FALSE])
    )
  }))
})) %>%
  mutate(shift = factor(shift, levels = shift_levels))

plot_shift <- function(season_name) {
  dat <- filter(effects8, season == season_name)
  effect_range <- range(c(0, dat$lower, dat$upper))
  span <- max(diff(effect_range), 0.01)
  
  labels <- filter(dd8, season == season_name) %>%
    mutate(label_y = effect_range[1] - 0.12 * span)
  
  ggplot(dat, aes(shift, estimate, colour = metric)) +
    annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray50") +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(0.6),
                    linewidth = 0.8) +
    geom_text(data = labels, aes(shift, label_y, label = paste0("dd= ", dd)),
              inherit.aes = FALSE, size = 3) +
    scale_colour_manual(values = metric_colours) +
    scale_x_discrete(limits = shift_levels, labels = shift_labels, drop = FALSE) +
    labs(x = "Shift", y = "Effect size (%)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")
}

panels8 <- setNames(lapply(seasons, plot_shift), tolower(seasons))
figure8 <- wrap_plots(panels8, nrow = 1)

for (name in names(panels8)) {
  save_plot(panels8[[name]], fig8_dir, name, 800, 1200)
}

save_plot(figure8, fig8_dir, "Figure8_combined_AR1", 2000, 1000)
diagnostics8 <- run_diagnostics(models8, fig8_dir)

message("Figures and diagnostic summaries completed.")