07_make_forecasts
================
JTK
2025-01-14

################################################################################ 

This script takes the models we have trained on observational streamflow
and watershed characteristic data, as well as the dataframes we have
prepared with National Water Model forecast data, to make forecasts of
total phosphorus and chloride concentration from those streamflow
forecasts. It results in forecasts of concentration for water year
2022-2023 for both total phosphorus and chloride.

**Inputs**

1)  Trained LightGBM models for gaged and ungaged scenario. These are
    most easily input from text files that contain the model
    architecture, which we wrote to file in *04_model_development*

2)  `forecasting_data_tp_gaged`, `forecasting_data_tp_ungaged`,
    `forecasting_data_chlor_gaged`, `forecasting_data_chlor_ungaged`
    (from *06_forecast_data_prep*): dataframes containing all the data
    (nwm forecasts, antecdent conditions, and watershed characteristics
    that we need to make water quality concentration forecasts)

**Outputs**

1)  Forecast timeseries for total phosphorus and chloride in gaged and
    ungaged scenario (`forecast_timeseries_tp_gaged`,
    `forecast_timeseries_chlor_gaged`,
    `forecast_timeseries_chlor_gaged`,
    `forecast_timeseries_chlor_ungaged`)

2)  Performance statistics for those forecasts
    (`forecast_stats_tp_gaged`, `forecast_stats_chlor_gaged`,
    `forecast_stats_tp_ungaged`, `forecast_stats_chlor_ungaged`)

################################################################################ 

# Housekeeping

### Packages

``` r
### Data mgmt
require(tidyverse)

### Modeling
require(hydroGOF)
require(lightgbm)
require(SHAPforxgboost)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/02_watershed_attributes_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/03_observational_data_prep.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/04_model_development.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/05_forecast_data_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/06_forecast_data_prep.Rmd"), quiet=TRUE))
```

# Forecast with NWM

### Gaged Scenario

``` r
#### Forecast for Total Phosphorus

nwm_forecasts_tp_gaged <- lgbm_predictor(predictors_df = forecasting_data_tp_gaged, 
     observed_df = tp_test %>%
       dplyr::select(!group_id),
     D_factor = d_fact_tp_gaged$D_fact,
     constituent = "Phosphorus_Total",
     solo_or_lumped_or_loocv = "lumped",
     model_file_location = here("data", "models", "lumped"),
     nwm=TRUE)

#### Forecast for Chloride

nwm_forecasts_chlor_gaged <- lgbm_predictor(predictors_df = forecasting_data_chlor_gaged, 
     observed_df = chlor_test %>%
       dplyr::select(!group_id),
     D_factor = d_fact_chlor_gaged$D_fact,
     constituent = "Chloride",
     solo_or_lumped_or_loocv = "lumped",
     model_file_location = here("data", "models", "lumped"),
     nwm=TRUE)
```

### Ungaged Scenario

``` r
### Forecast for Total Phosphorus

nwm_forecasts_tp_ungaged <- pmap(list(predictors_df = forecasting_data_tp_ungaged$nwm_data, 
     observed_df = forecasting_data_tp_ungaged$observed_data,
     D_factor = forecasting_data_tp_ungaged$D_fact),
     .f = lgbm_predictor,
     constituent = "Phosphorus_Total",
     solo_or_lumped_or_loocv = "loocv",
     model_file_location = here("data", "models", "loocv"),
     .progress = TRUE)


### Forecast for Chloride

nwm_forecasts_chlor_ungaged <- pmap(list(predictors_df = forecasting_data_chlor_ungaged$nwm_data, 
     observed_df = forecasting_data_chlor_ungaged$observed_data,
     D_factor = forecasting_data_chlor_ungaged$D_fact),
     .f = lgbm_predictor,
     constituent = "Chloride",
     solo_or_lumped_or_loocv = "loocv",
     model_file_location = here("data", "models", "loocv"),
     .progress = TRUE)
```

# Extract forecasted values and other relevant information

### Gaged

``` r
################################################################################

#### Total Phosphorus

### Time series of forecasted and observed

forecast_timeseries_tp_gaged <- nwm_forecasts_tp_gaged[[1]] %>%
  dplyr::select(tributary, predict_date, lead_days, 
                forecasted_conc, observed_conc, 
                constituent)

### Performance stats for forecasts

forecast_stats_tp_gaged <- nwm_forecasts_tp_gaged[[2]]

################################################################################

#### Chloride

### Time series of forecasted and observed

forecast_timeseries_chlor_gaged <- nwm_forecasts_chlor_gaged[[1]] %>%
  dplyr::select(tributary, predict_date, lead_days, 
                forecasted_conc, observed_conc, 
                constituent)

### Performance stats for forecasts

forecast_stats_chlor_gaged <- nwm_forecasts_chlor_gaged[[2]]

################################################################################
```

### Ungaged

``` r
################################################################################

#### Total Phosphorus

### Time series of forecasted and observed

forecast_timeseries_tp_ungaged <- map_dfr(nwm_forecasts_tp_ungaged,1) %>%
  dplyr::select(tributary, predict_date, lead_days, 
                forecasted_conc, observed_conc, 
                constituent)

### Performance stats for forecasts

forecast_stats_tp_ungaged <- map_dfr(nwm_forecasts_tp_ungaged,2)

################################################################################

#### Chloride

### Time series of forecasted and observed

forecast_timeseries_chlor_ungaged <- map_dfr(nwm_forecasts_chlor_ungaged,1) %>%
  dplyr::select(tributary, predict_date, lead_days, 
                forecasted_conc, observed_conc, 
                constituent)

### Performance stats for forecasts

forecast_stats_chlor_ungaged <- map_dfr(nwm_forecasts_chlor_ungaged,2)

################################################################################
```
