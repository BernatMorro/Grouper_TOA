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
# devtools::install_github("aspillaga/activityTOA", build_vignettes = TRUE)
library(activityTOA)
# devtools::install_github("aspillaga/BTNtools")
library(BTNtools)
library(glmmTMB)
library(effects)
library(gratia)
library(performance)
library(DHARMa)


# Load tracking data (diving shift grouping)
detect <- fread("./detections.csv")
toro <- fread("./divers.csv")
environmentals <- fread("./environmentals.csv")

# Join datasets
dataset <- detect %>%
  dplyr::left_join(toro %>% dplyr::select(shift_id, divers), by = "shift_id")
dataset$divers[is.na(dataset$divers)] <- 0

dataset <- dataset %>%
  mutate(day = as.IDate(day)) %>%
  left_join(
    environmentals %>% select(day, wave_height, cloud_cover, rain, moon),
    by = "day")


# Add month and year
dataset$month <- format(as.Date(dataset$dia), "%m")
dataset$year <- format(as.Date(dataset$dia), "%Y")


# Keep only diving shifts of interest and delete rows without telemetry data
dataset <- subset(dataset, shift %in% c("[8,11)", "[11,13)", "[13,18)", "[19,23)"))
dataset <- dataset[!is.na(dataset$step_sd), ] 

#dataset$date <- as.Date(dataset$dia)

# center temperature for polinomial glmm compatibility 
dataset <- dataset %>%
  mutate(temp_c = scale(temperature, center = TRUE, scale = FALSE))


#We define seasons
dataset$season <- case_when(
  (format(dataset$date, "%m-%d") >= "12-21" | format(dataset$date, "%m-%d") < "03-20") ~ "Winter",
  (format(dataset$date, "%m-%d") >= "03-20" & format(dataset$date, "%m-%d") < "06-21") ~ "Spring",
  (format(dataset$date, "%m-%d") >= "06-21" & format(dataset$date, "%m-%d") < "09-22") ~ "Summer",
  (format(dataset$date, "%m-%d") >= "09-22" & format(dataset$date, "%m-%d") < "12-21") ~ "Autumn"
)

winter <- subset(dataset, season == "Winter")
spring <- subset(dataset, season == "Spring")
summer <- subset(dataset, season == "Summer")
autumn <- subset(dataset, season == "Autumn")

#We remove a diving shift with no dives in winter
winter <- winter %>%
  filter(shift != "[19,23)") %>%  
  mutate(shift = droplevels(as.factor(shift)))


################################## Analysis for Figure 5 and 7 #################################

# Fit GLMM
winter_mean <- glmmTMB(step_mean ~ divers + moon + temp_c + I(temp_c^2) +  wave_height + cloud_cover + shift +
                         (1 | fish_id) + (1 | year),
                       data = winter, family = Gamma(link = "log"))
spring_mean <- glmmTMB(step_mean ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift +
                         (1 | fish_id) + (1 | year),
                       data = spring, family = Gamma(link = "log"))
summer_mean <- glmmTMB(step_mean ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift +
                         (1 | fish_id) + (1 | year),
                       data = summer, family = Gamma(link = "log"))
autumn_mean <- glmmTMB(step_mean ~ divers + moon + temp_c + I(temp_c^2) +  wave_height + cloud_cover + shift +
                         (1 | fish_id) + (1 | year),
                       data = autumn, family = Gamma(link = "log"))
overall_mean <- glmmTMB(step_mean ~ divers + moon + temp_c + I(temp_c^2) +  wave_height + cloud_cover + shift + season +
                          (1 | fish_id) + (1 | year),
                        data = dataset, family = Gamma(link = "log"))

winter_max <- glmmTMB(step_max ~ divers + moon + temp_c + I(temp_c^2) +  wave_height + cloud_cover + shift +
                        (1 | fish_id) + (1 | year),
                      data = winter, family = Gamma(link = "log"))
spring_max <- glmmTMB(step_max ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id) + (1 | year),
                      data = spring, family = Gamma(link = "log"))
summer_max <- glmmTMB(step_max ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id) + (1 | year),
                      data = summer, family = Gamma(link = "log"))
autumn_max <- glmmTMB(step_max ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift +
                        (1 | fish_id) + (1 | year),
                      data = autumn, family = Gamma(link = "log"))
overall_max <- glmmTMB(step_max ~ divers + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover + shift + season +
                         (1 | fish_id) + (1 | year),
                       data = dataset, family = Gamma(link = "log"))

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

## number of days with ≥1 diver per season + Overall
n_df <- dataset |>
  filter(divers > 0) |>
  distinct(season, date) |>
  count(season, name = "n") |>
  complete(season = seas_vec, fill = list(n = 0)) |>
  mutate(
    n = ifelse(
      season == "Overall",
      nrow(distinct(dataset |> filter(divers > 0), date)),
      n
    ),
    season = factor(season, levels = seas_vec)
  )

##---------------------------------------------------------------
## 1. Function to get log-link slopes (per unit) for one variable
##---------------------------------------------------------------
get_var_slopes <- function(mod_list, metric_label, var) {
  imap_dfr(mod_list, ~ {
    tr_df <- emtrends(.x, specs = ~ 1, var = var, type = "link") |> as.data.frame()
    tibble(
      season = .y,
      metric = metric_label,
      var    = var,
      slope  = as.numeric(tr_df[[paste0(var, ".trend")]]),
      lwr    = as.numeric(tr_df$asymp.LCL),
      upr    = as.numeric(tr_df$asymp.UCL)
    )
  })
}

##----------------------------------------------------------
## 2. Compute slopes for all variables and both metrics
##----------------------------------------------------------
vars <- c("divers", "temp_c", "moon", "cloud_cover", "wave_height")

effects_all <- map_dfr(vars, \(v) {
  bind_rows(
    get_var_slopes(mean_models, "mean", v),
    get_var_slopes(max_models,  "max",  v)
  )
}) |>
  mutate(
    season = factor(season, levels = seas_vec),
    metric = factor(metric, levels = c("mean", "max")),
    sig    = ifelse(lwr * upr > 0, "yes", "no")
  ) |>
  left_join(n_df, by = "season")

##----------------------------------------------------------
## 3. Plotting: one plot per variable
##----------------------------------------------------------
plot_var <- function(df, var_name, seas_vec) {
  
  var_titles <- c(
    divers      = "Divers",
    temp_c = "temp_c",
    moon        = "Moon",
    cloud_cover = "Cloud cover",
    wave_height = "Wave height"
  )
  
  d <- df |>
    filter(var == var_name) |>
    mutate(
      season     = factor(season, levels = seas_vec),
      slope_plot = 100 * (exp(slope) - 1),
      lwr_plot   = 100 * (exp(lwr)   - 1),
      upr_plot   = 100 * (exp(upr)   - 1)
    )
  
  offset <- 0.10 * diff(range(c(d$lwr_plot, d$upr_plot), na.rm = TRUE))
  y_lab  <- min(d$lwr_plot, na.rm = TRUE) - offset
  
  ggplot(d, aes(season, slope_plot, colour = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin = lwr_plot, ymax = upr_plot),
                    position = position_dodge(0.6),
                    fatten   = 1.2) +
    geom_text(
      data = d |> filter(metric == "mean"),
      aes(x = season, label = paste0("n = ", n)),
      y = y_lab,
      size = 3,
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
    theme(legend.position = "none")
}

p_divers      <- plot_var(effects_all, "divers",      seas_vec)
p_temp_c <- plot_var(effects_all, "temp_c", seas_vec)
p_moon        <- plot_var(effects_all, "moon",        seas_vec)
p_cloud       <- plot_var(effects_all, "cloud_cover", seas_vec)
p_wave        <- plot_var(effects_all, "wave_height", seas_vec)

print(p_divers + scale_y_continuous(limits = c(-1.03, 0.3), expand = expansion(mult = 0.02)))
print(p_temp_c + ggtitle("Temperature") + scale_y_continuous(limits = c(-9, 10),expand = expansion(mult = 0.02)))
print(p_moon + scale_y_continuous(limits = c(-22, 12), expand = expansion(mult = 0.02)))
print(p_cloud + scale_y_continuous(limits = c(-23.5, 22), expand = expansion(mult = 0.02)))
print(p_wave + scale_y_continuous(limits = c(-36, 10), expand = expansion(mult = 0.02)))


######################################### Analysis for Figure 6 ########################################
##----------------------------------------------------------
## 0) Shift levels
##----------------------------------------------------------
shift_levs <- c("[8,11)", "[11,13)", "[13,18)", "[19,23)")
shift_labs <- c("8-11","11-13","13-18","19-23")

##----------------------------------------------------------
## 1) Fit models (divers * shift) for each season
##----------------------------------------------------------
winter_mean <- glmmTMB(
  step_mean ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = winter, family = Gamma(link = "log")
)
spring_mean <- glmmTMB(
  step_mean ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = spring, family = Gamma(link = "log")
)
summer_mean <- glmmTMB(
  step_mean ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = summer, family = Gamma(link = "log")
)
autumn_mean <- glmmTMB(
  step_mean ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = autumn, family = Gamma(link = "log")
)

winter_max <- glmmTMB(
  step_max ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = winter, family = Gamma(link = "log")
)
spring_max <- glmmTMB(
  step_max ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = spring, family = Gamma(link = "log")
)
summer_max <- glmmTMB(
  step_max ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = summer, family = Gamma(link = "log")
)
autumn_max <- glmmTMB(
  step_max ~ divers*shift + moon + temp_c + I(temp_c^2) + wave_height + cloud_cover +
    (1 | fish_id) + (1 | year),
  data = autumn, family = Gamma(link = "log")
)

## Pack for looping
mean_models <- list(Winter = winter_mean, Spring = spring_mean, Summer = summer_mean, Autumn = autumn_mean)
max_models  <- list(Winter = winter_max,  Spring = spring_max,  Summer = summer_max,  Autumn = autumn_max)
season_data <- list(Winter = winter,      Spring = spring,      Summer = summer,      Autumn = autumn)

##----------------------------------------------------------
## 2) Extract shift-specific slopes (divers)
##----------------------------------------------------------

get_shift_slopes <- function(model, season_name, metric, shifts) {
  
  # shift-specific beta_divers|shift on link scale
  emtrends(model, specs = ~ shift, var = "divers", type = "link") |>
    as.data.frame() |>
    rename(slope = divers.trend, lwr = asymp.LCL, upr = asymp.UCL) |>
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

##----------------------------------------------------------
## 3) n per shift (days with divers > 0) by season
##----------------------------------------------------------

get_shift_n <- function(dat, season_name, shifts) {
  dat |>
    filter(divers > 0) |>
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

##----------------------------------------------------------
## 4) Plot function: one plot per season
##----------------------------------------------------------

plot_season_shift <- function(df, season_name, shift_levs, shift_labs) {
  
  d <- df |>
    filter(season == season_name) |>
    mutate(shift = factor(shift, levels = shift_levs))
  
  # label y position
  offset <- 0.10 * diff(range(c(d$lwr_plot, d$upr_plot), na.rm = TRUE))
  y_lab  <- min(d$lwr_plot, na.rm = TRUE) - offset
  
  p <- ggplot(d, aes(shift, slope_plot, colour = metric)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(ymin = lwr_plot, ymax = upr_plot),
                    position = position_dodge(0.6),
                    fatten   = 1.2) +
    geom_text(
      data = d |> filter(metric == "mean"),
      aes(x = shift, label = paste0("n = ", n), colour = NULL),
      y = y_lab,
      size = 3
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
    theme(legend.position = "none")
  
  return(p)
}


p_winter <- plot_season_shift(effects_shift, "Winter", shift_levs, shift_labs)
p_spring <- plot_season_shift(effects_shift, "Spring", shift_levs, shift_labs)
p_summer <- plot_season_shift(effects_shift, "Summer", shift_levs, shift_labs)
p_autumn <- plot_season_shift(effects_shift, "Autumn", shift_levs, shift_labs)

print(p_winter +  scale_y_continuous(limits = c(-4.2, 4.5), expand = expansion(mult = 0.02)))
print(p_spring +  scale_y_continuous(limits = c(-1.3, 3), expand = expansion(mult = 0.02)))
print(p_summer +  scale_y_continuous(limits = c(-1.22, 1), expand = expansion(mult = 0.02)))
print(p_autumn + scale_y_continuous(limits = c(-2.3, 1), expand = expansion(mult = 0.02)))

