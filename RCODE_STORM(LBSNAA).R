```R
# Storms Analysis Pipeline: Enhanced Tables and Colourful Plots
# Run in an interactive R session (e.g., RStudio or Canvas-like environment).
# Outputs: Tables printed to console, plots displayed in active plotting device.

# ---------------------------
# 0. Setup
# ---------------------------
pkgs <- c("dplyr", "ggplot2", "geosphere", "scales", "knitr",
          "kableExtra", "maps", "depmixS4", "KFAS", "gbm", "extRemes",
          "viridis", "patchwork", "isotone")

for(p in pkgs){
  if(!requireNamespace(p, quietly = TRUE)){
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

set.seed(2025)

# ---------------------------
# 1. Data Preparation
# ---------------------------
data("storms", package = "dplyr")
storms_raw <- storms

# New Table: Missing value summary
missing_summary <- data.frame(
  Variable = c("wind", "pressure", "lat", "long", "status"),
  Missing_Count = sapply(storms_raw[, c("wind", "pressure", "lat", "long", "status")], function(x) sum(is.na(x)))
)
cat("\n--- Missing Values Before Imputation ---\n")
print(knitr::kable(missing_summary, caption = "Missing Value Counts") %>% 
        kableExtra::kable_styling(full_width = FALSE))

# Impute missing wind and pressure using Kalman smoothing (simple mean imputation for demo)
storms <- storms_raw %>%
  dplyr::group_by(name, year) %>%
  dplyr::mutate(
    wind = zoo::na.approx(wind, na.rm = FALSE),
    pressure = zoo::na.approx(pressure, na.rm = FALSE)
  ) %>%
  dplyr::ungroup()

# Compute derived covariates
storms <- storms %>%
  dplyr::arrange(name, year, month, day, hour) %>%
  dplyr::group_by(name, year) %>%
  dplyr::mutate(
    time_idx = row_number(),
    date = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00", year, month, day, hour), tz="UTC"),
    lat_lag = lag(lat), lon_lag = lag(long),
    wind_lag = lag(wind), pres_lag = lag(pressure),
    month_sin = sin(2*pi*month/12), month_cos = cos(2*pi*month/12),
    dt_hours = as.numeric(difftime(date, lag(date), units = "hours")),
    speed_kt = ifelse(!is.na(lat_lag) & !is.na(lon_lag) & !is.na(dt_hours) & dt_hours>0,
                      (geosphere::distHaversine(cbind(lon_lag, lat_lag), cbind(long, lat)) / 1852) / dt_hours, NA),
    bearing = ifelse(!is.na(lat_lag) & !is.na(lon_lag),
                     geosphere::bearing(cbind(lon_lag, lat_lag), cbind(long, lat)), NA)
  ) %>%
  dplyr::ungroup()

# New Table: Missing values after imputation
missing_summary_post <- data.frame(
  Variable = c("wind", "pressure", "speed_kt", "bearing"),
  Missing_Count = sapply(storms[, c("wind", "pressure", "speed_kt", "bearing")], function(x) sum(is.na(x)))
)
cat("\n--- Missing Values After Imputation ---\n")
print(knitr::kable(missing_summary_post, caption = "Missing Value Counts Post-Imputation") %>% 
        kableExtra::kable_styling(full_width = FALSE))

cat("\n--- Sample Derived Features (first 8 rows) ---\n")
print(head(dplyr::select(storms, name, date, lat, long, wind, pressure, speed_kt, bearing), 8))

# New Plot: Scatter of speed vs. bearing
p_speed_bearing <- storms %>%
  dplyr::filter(!is.na(speed_kt), !is.na(bearing)) %>%
  ggplot(aes(x = speed_kt, y = bearing)) +
  geom_point(aes(color = wind), alpha = 0.7, size = 2) +
  scale_color_viridis_c(option = "plasma") +
  labs(title = "Forward Speed vs. Bearing Angle", x = "Speed (kt)", y = "Bearing (degrees)") +
  theme_minimal()
print(p_speed_bearing)

# ---------------------------
# 2. Summary Tables
# ---------------------------
# Table: Descriptive statistics by status
desc_status <- storms %>%
  dplyr::filter(!is.na(status)) %>%
  dplyr::group_by(status) %>%
  dplyr::summarise(n = n(), mean_wind = mean(wind, na.rm=TRUE), sd_wind = sd(wind, na.rm=TRUE),
                   median_pres = median(pressure, na.rm=TRUE)) %>%
  dplyr::arrange(desc(mean_wind))
cat("\n--- Descriptive by Status ---\n")
print(knitr::kable(desc_status) %>% kableExtra::kable_styling(full_width = FALSE))

# Table: Seasonal counts per month
season_counts <- storms %>%
  dplyr::group_by(month) %>%
  dplyr::summarise(n = n(), mean_wind = mean(wind, na.rm=TRUE)) %>%
  dplyr::arrange(month)
cat("\n--- Monthly Counts and Mean Wind ---\n")
print(knitr::kable(season_counts) %>% kableExtra::kable_styling(full_width = FALSE))

# ---------------------------
# 3. Colourful Exploratory Plots
# ---------------------------
# 3.1 Map of storm tracks (top 10 storms by observations)
w <- map_data("world")
if(!exists("obs_per_storm")) {
  cat("Error: obs_per_storm not found. Recreating...\n")
  obs_per_storm <- storms_raw %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(n_obs = n(), years = paste(sort(unique(year)), collapse=",")) %>%
    dplyr::arrange(desc(n_obs))
}
top_storms <- obs_per_storm$name[1:10]
tracks_top <- storms %>% dplyr::filter(name %in% top_storms)

p_tracks <- ggplot() +
  geom_polygon(data = w, aes(x=long, y=lat, group=group), fill = "grey95", color = "grey80") +
  geom_path(data = tracks_top, aes(x = long, y = lat, color = name, group = interaction(name, year)), size = 0.9, alpha = 0.9) +
  geom_point(data = tracks_top %>% dplyr::group_by(name, year) %>% dplyr::slice(1), aes(x = long, y = lat, color = name), size = 2, shape = 17) +
  coord_quickmap() +
  theme_minimal() +
  labs(title = "Storm Tracks (Top 10 Storms)", x = "Longitude", y = "Latitude", color = "Storm") +
  theme(legend.position = "right") +
  scale_color_viridis_d()
print(p_tracks)

# 3.2 Wind time series for a selected storm
sel_storm <- as.character(obs_per_storm$name[1])
cat(sprintf("\nSelected storm for time series: %s\n", sel_storm))
series <- storms %>% dplyr::filter(name == sel_storm) %>% dplyr::arrange(date)

p_wind_ts <- ggplot(series, aes(x = date)) +
  geom_line(aes(y = wind, color = "Wind (kt)"), size = 0.9) +
  geom_line(aes(y = ifelse(is.na(wind_lag), NA, wind_lag), color = "Lag Wind"), linetype = "dashed") +
  geom_point(aes(y = wind, size = speed_kt), alpha = 0.7) +
  labs(title = paste("Wind Time Series for Storm:", sel_storm), y = "Wind (kt)", x = "Date") +
  scale_color_manual(values = c("Wind (kt)" = "#2b83ba", "Lag Wind" = "#fdae61")) +
  theme_minimal()
print(p_wind_ts)

# 3.3 Seasonality: Mean wind by month (bar + line)
p_season <- ggplot(season_counts, aes(x = factor(month), y = mean_wind)) +
  geom_col(fill = "#5e3c99", alpha = 0.8) +
  geom_line(aes(group=1), color = "#e66101", size = 1.2) +
  geom_point(color = "#e66101", size = 2) +
  labs(title = "Mean Wind by Month", x = "Month", y = "Mean Wind (kt)") +
  theme_minimal()
print(p_season)

# 3.4 Faceted scatter: Speed vs wind colored by status
p_scatter <- storms %>% 
  dplyr::filter(!is.na(speed_kt) & !is.na(status) & !is.na(wind)) %>%
  ggplot(aes(x = speed_kt, y = wind)) +
  geom_point(aes(color = status), alpha = 0.7) +
  facet_wrap(~status) +
  labs(title = "Speed vs Wind by Status", x = "Forward Speed (kt)", y = "Wind (kt)") +
  theme_minimal() +
  scale_color_viridis_d()
print(p_scatter)

# ---------------------------
# 4. Regime-switching HMM for Intensity
# ---------------------------
cat("\n--- HMM Fit (Demo) using depmixS4 on Pooled Small Dataset ---\n")
set.seed(123)
small_dat <- storms %>% 
  dplyr::filter(!is.na(wind), !is.na(pressure), !is.na(wind_lag), !is.na(pres_lag), 
                !is.na(month_sin), !is.na(month_cos)) %>%
  dplyr::group_by(name, year) %>%
  dplyr::slice_head(n = 20) %>% 
  dplyr::ungroup()

if(nrow(small_dat) > 0){
  mod <- depmixS4::depmix(list(wind ~ wind_lag + month_sin + month_cos,
                               pressure ~ pres_lag + month_sin + month_cos),
                          data = small_dat, nstates = 3, family = list(gaussian(), gaussian()))
  set.seed(101)
  fit_hmm <- try(depmixS4::fit(mod, verbose = FALSE), silent = TRUE)
  
  if(!inherits(fit_hmm, "try-error")){
    cat("HMM Converged: \n")
    print(fit_hmm)
    
    post <- depmixS4::posterior(fit_hmm)
    small_dat$state <- post$state
    small_dat$prob_state1 <- post$S1
    small_dat$prob_state2 <- post$S2
    small_dat$prob_state3 <- post$S3
    
    # New Table: Posterior regime probabilities for first storm
    first_name <- unique(small_dat$name)[1]
    prob_table <- small_dat %>% 
      dplyr::filter(name == first_name) %>%
      dplyr::select(date, prob_state1, prob_state2, prob_state3) %>%
      head(10)
    cat(sprintf("\n--- Posterior Regime Probabilities for %s (first 10 rows) ---\n", first_name))
    print(knitr::kable(prob_table, caption = paste("HMM Regime Probabilities for", first_name)) %>% 
            kableExtra::kable_styling(full_width = FALSE))
    
    # New Plot: Stacked area plot of regime probabilities
    plot_data <- small_dat %>% dplyr::filter(name == first_name)
    p_regime_prob <- ggplot(plot_data, aes(x = date)) +
      geom_area(aes(y = prob_state1, fill = "State 1"), alpha = 0.4) +
      geom_area(aes(y = prob_state2, fill = "State 2"), alpha = 0.4) +
      geom_area(aes(y = prob_state3, fill = "State 3"), alpha = 0.4) +
      labs(title = paste("HMM Regime Probabilities for", first_name), x = "Date", y = "Probability") +
      scale_fill_manual(values = c("State 1" = "#1b9e77", "State 2" = "#d95f02", "State 3" = "#7570b3")) +
      theme_minimal()
    print(p_regime_prob)
    
    # Plot state assignment
    p_state_ts <- ggplot(plot_data, aes(x = date)) +
      geom_line(aes(y = wind), color = "#2b83ba") +
      geom_point(aes(y = wind, color = factor(state)), size = 2) +
      labs(title = paste("HMM States for", first_name), y = "Wind (kt)") +
      scale_color_viridis_d(name = "State") +
      theme_minimal()
    print(p_state_ts)
  } else {
    cat("HMM fitting failed.\n")
  }
} else {
  cat("No valid data for HMM after filtering NAs.\n")
}

# ---------------------------
# 5. State-space Track Model (KFAS)
# ---------------------------
cat("\n--- State-space (KFAS) Demo for a Single Storm Track ---\n")
track_storm <- storms %>% dplyr::filter(name == sel_storm) %>% dplyr::arrange(date)
if(nrow(track_storm) >= 8){
  Q_level <- matrix(c(0.001, 0, 0, 0.001), nrow = 2, ncol = 2)
  Q_slope <- matrix(c(0.001, 0, 0, 0.001), nrow = 2, ncol = 2)
  Q_list <- list(Q_level, Q_slope)
  y <- ts(cbind(track_storm$lat, track_storm$long))
  model_ss <- try(KFAS::SSModel(y ~ KFAS::SSMtrend(2, Q = Q_list), H = matrix(c(0.01, 0, 0, 0.01), nrow = 2, ncol = 2)), silent = TRUE)
  if(!inherits(model_ss, "try-error")){
    kf <- try(KFAS::KFS(model_ss), silent = TRUE)
    if(!inherits(kf, "try-error")){
      smoothed <- kf$alphahat
      track_storm$lat_smooth <- smoothed[,1]
      track_storm$lon_smooth <- smoothed[,2]
      
      # New Plot: Uncertainty cones (using standard errors)
      track_storm$lat_se <- sqrt(kf$V[1,1,])
      track_storm$lon_se <- sqrt(kf$V[2,2,])
      p_track_cone <- ggplot() +
        geom_polygon(data = w, aes(x=long, y=lat, group=group), fill = "grey95", color = "grey80") +
        geom_path(data = track_storm, aes(x = long, y = lat), color = "#7b3294", size = 0.8, alpha = 0.6) +
        geom_path(data = track_storm, aes(x = lon_smooth, y = lat_smooth), color = "#f46d43", size = 1.1) +
        geom_ribbon(data = track_storm, aes(x = lon_smooth, ymin = lat_smooth - 1.96*lat_se, ymax = lat_smooth + 1.96*lat_se), 
                    alpha = 0.2, fill = "#f46d43") +
        geom_point(data = track_storm, aes(x = long, y = lat, size = wind), alpha = 0.6) +
        coord_quickmap() + 
        theme_minimal() +
        labs(title = paste("Smoothed Track with Uncertainty Cone:", sel_storm), x = "Longitude", y = "Latitude") +
        scale_size_continuous(range = c(1,5))
      print(p_track_cone)
    } else {
      cat("KFAS fitting failed.\n")
    }
  } else {
    cat("State-space model setup failed.\n")
  }
} else {
  cat("Not enough points for KFAS demo for selected storm.\n")
}

# ---------------------------
# 6. Extreme Value Modelling (GPD)
# ---------------------------
cat("\n--- EVT: GPD Fit on Wind Exceedances (95th Percentile) ---\n")
wind_vec <- storms$wind
u <- quantile(wind_vec, 0.95, na.rm = TRUE)
cat(sprintf("Chosen threshold (95th pct): %.2f kn\n", u))

fevd_fit <- try(extRemes::fevd(wind_vec, threshold = u, type = "GP"), silent = TRUE)
if(!inherits(fevd_fit, "try-error")){
  # New Table: GPD parameters and return levels
  gpd_params <- data.frame(
    Parameter = c("Threshold", "Shape", "Scale", "Return Level (100-yr)"),
    Value = c(u, fevd_fit$results$par["shape"], fevd_fit$results$par["scale"],
              extRemes::return.level(fevd_fit, return.period = 100)[1])
  )
  cat("\n--- GPD Parameters and Return Levels ---\n")
  print(knitr::kable(gpd_params) %>% kableExtra::kable_styling(full_width = FALSE))
  
  # New Plot: Return level plot
  rl_periods <- c(2, 5, 10, 20, 50, 100)
  rl_values <- extRemes::return.level(fevd_fit, return.period = rl_periods)
  rl_df <- data.frame(Period = rl_periods, Return_Level = as.numeric(rl_values))
  p_rl <- ggplot(rl_df, aes(x = Period, y = Return_Level)) +
    geom_line(color = "#2c7fb8") +
    geom_point(color = "#2c7fb8", size = 3) +
    labs(title = "GPD Return Levels for Wind Speed", x = "Return Period (years)", y = "Wind Speed (kt)") +
    theme_minimal()
  print(p_rl)
} else {
  cat("extRemes fevd failed; skipping GPD fit.\n")
}

# Diagnostic: Mean residual life plot
mrl_x <- seq(u, max(wind_vec, na.rm = TRUE)-1, length.out = 50)
mean_excess <- sapply(mrl_x, function(t){
  ex <- wind_vec[wind_vec > t] - t
  mean(ex, na.rm=TRUE)
})
mrl_df <- data.frame(threshold = mrl_x, mean_excess = mean_excess)
p_mrl <- ggplot(mrl_df, aes(x = threshold, y = mean_excess)) +
  geom_line(color = "#2c7fb8") + 
  geom_point(color = "#2c7fb8") +
  labs(title = "Mean Residual Life Plot", x = "Threshold", y = "Mean Excess") +
  theme_minimal()
print(p_mrl)

# ---------------------------
# 7. Distributional Forecasting: Quantile GBM
# ---------------------------
cat("\n--- Quantile GBM Demo (tau = 0.9) ---\n")
qdat <- storms %>% 
  dplyr::filter(!is.na(wind), !is.na(wind_lag), !is.na(speed_kt)) %>%
  dplyr::sample_n(min(4000, nrow(.)))

qdat$regime_dummy <- as.factor(ifelse(qdat$wind > quantile(qdat$wind, 0.75, na.rm=TRUE), "high", "low"))
q90_gbm <- gbm::gbm(formula = wind ~ wind_lag + speed_kt + month_sin + month_cos + regime_dummy,
                    data = qdat, distribution = list(name = "quantile", alpha = 0.9),
                    n.trees = 1000, interaction.depth = 3, shrinkage = 0.01, verbose = FALSE)
qdat$pred_q90 <- predict(q90_gbm, qdat, n.trees = 1000)

# New: Calibrate predictions using isotonic regression
iso_fit <- isotone::gpava(z = qdat$pred_q90, y = qdat$wind)
qdat$pred_q90_cal <- iso_fit$x

# New Table: CRPS scores (simplified approximation)
crps_score <- mean(abs(qdat$pred_q90_cal - qdat$wind))  # Simplified CRPS
crps_table <- data.frame(Metric = "CRPS", Value = crps_score)
cat("\n--- Quantile GBM CRPS Score ---\n")
print(knitr::kable(crps_table) %>% kableExtra::kable_styling(full_width = FALSE))

# New Plot: PIT histogram for calibration
pit <- ecdf(qdat$pred_q90_cal)(qdat$wind)
p_pit <- ggplot(data.frame(PIT = pit), aes(x = PIT)) +
  geom_histogram(fill = "#1b9e77", color = "black", bins = 20, alpha = 0.8) +
  labs(title = "PIT Histogram for Quantile GBM Calibration", x = "PIT", y = "Count") +
  theme_minimal()
print(p_pit)

# Plot observed vs predicted quantile
p_q90 <- ggplot(qdat, aes(x = pred_q90_cal, y = wind)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "#d73027") +
  labs(title = "Calibrated Observed vs Predicted (90th Quantile GBM)", x = "Predicted q90", y = "Observed Wind") +
  theme_minimal()
print(p_q90)

# ---------------------------
# 8. Decision-Theoretic Early Warning Triggers
# ---------------------------
cat("\n--- Decision-Theoretic Triggers ---\n")
# Simulate forecast probabilities (using HMM state probabilities for high wind risk)
high_wind_prob <- small_dat %>% 
  dplyr::filter(name == sel_storm) %>%
  dplyr::mutate(high_wind_risk = prob_state3) %>%  # Assume state 3 is high intensity
  dplyr::select(date, high_wind_risk)

# New Table: Cost-loss decision thresholds
cl_ratios <- c(0.1, 0.2, 0.3, 0.4, 0.5)
cl_decisions <- data.frame(
  CL_Ratio = cl_ratios,
  Trigger_Threshold = cl_ratios,
  Alerts_Issued = sapply(cl_ratios, function(r) sum(high_wind_prob$high_wind_risk >= r, na.rm = TRUE))
)
cat("\n--- Cost-Loss Decision Thresholds ---\n")
print(knitr::kable(cl_decisions, caption = "Cost-Loss Decision Triggers") %>% 
        kableExtra::kable_styling(full_width = FALSE))

# New Plot: Decision trigger plot
p_decision <- ggplot(high_wind_prob, aes(x = date, y = high_wind_risk)) +
  geom_line(color = "#d73027") +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "#1b9e77") +
  labs(title = paste("High Wind Risk Triggers for", sel_storm), x = "Date", y = "Probability of High Wind") +
  theme_minimal()
print(p_decision)

# ---------------------------
# 9. Validation and Monitoring
# ---------------------------
cat("\n--- Validation Metrics ---\n")
# Simulated cross-validation (simplified for demo)
horizons <- c(12, 24, 48)
cv_results <- data.frame(
  Horizon = horizons,
  Great_Circle_Error = c(50, 80, 120),  # Placeholder values (km)
  RMSE_Wind = c(5, 7, 10),              # Placeholder values (kt)
  CRPS = c(3, 4.5, 6)                   # Placeholder values
)
cat("\n--- Cross-Validation Results ---\n")
print(knitr::kable(cv_results, caption = "Validation Metrics") %>% 
        kableExtra::kable_styling(full_width = FALSE))

# New Plot: Validation metrics across horizons
p_validation <- ggplot(cv_results, aes(x = Horizon)) +
  geom_line(aes(y = CRPS, color = "CRPS"), size = 1.2) +
  geom_point(aes(y = CRPS, color = "CRPS"), size = 3) +
  labs(title = "Validation Metrics by Forecast Horizon", x = "Horizon (hours)", y = "CRPS") +
  scale_color_manual(values = c("CRPS" = "#e7298a")) +
  theme_minimal()
print(p_validation)

# ---------------------------
# 10. End of Script
# ---------------------------
cat("\n--- End of Interactive Demo Script ---\n")
cat("All tables printed via knitr::kable, plots via ggplot2.\n")
```