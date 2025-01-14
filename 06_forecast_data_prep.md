06_forecast_data_prep
================
JTK
2025-01-13

################################################################################ 

This script prepares various dataframes to make total phos and chloride
predictions from NWM forecasts in each basin. The most difficult task
with regards to this preparation is the creation of antecedent condition
dataframes, which account for antecedent daily, weekly, monthly, etc.
streamflow. Once these are prepped, we can make the actual forecasts

**Inputs**

1)  NWM streamflow forecast files, which we downloaded in
    *05_forecast_data_download*

2)  `watershed_chars` (from *02_watershed_attributes_download*): a
    dataframe containing static watershed attributes for all the Lake
    Champlain tributary watershed

3)  `picked_vars_tp` and `picked_vars_chlor` (from
    *04_model_development*): dataframes containing the list of features
    that we selected to include in the models to predict total phos and
    chloride concentration in a **gaged** scenario

4)  `picked_vars_tp_ungaged` and `picked_vars_chlor_ungaged` (from
    *04_model_development*): dataframes containing the list of features
    that we selected to include in the models to predict total phos and
    chloride concentration in an **ungaged** scenario

**Outputs**

1)  Dataframes containing all the data we need to forecast concentration
    from NWM streamflow using each of our LightGBM models
    (`forecasting_data_tp_gaged`, `forecasting_data_chlor_gaged`,
    `forecasting_data_tp_ungaged`, and `forecasting_data_tp_ungaged`)

################################################################################ 

# Housekeeping

### Packages

``` r
# Data mgmt
require(tidyverse)

## Parallel computing
require(future)
require(furr)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/02_watershed_attributes_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/03_observational_data_prep.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/04_model_development.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/05_forecast_data_download.Rmd"), quiet=TRUE))
```

# Transform NWM

### Bring in forecast files from where we saved them during download

``` r
##### This is for everything else (Oct 2021-Sept 2022; Nov 2022-Sept 2023)

#### Set forecast file location and list all the files located there

loc <- here("nwm_operational/medium_term")

mt_files <- list.files(loc)

mt_downloads <- list()

#### Import all the downloaded files

for(i in 1:length(mt_files)){
  
  file_name <- paste0(loc, mt_files[i])
  
  mt <- read_csv(file_name)
  
  mt_downloads[[i]] <- mt 
  
  
}

#### Trim down to Lake Champlain gages

mt_downloads_lc <- map(mt_downloads, ~(.x %>%
                                         filter(as.character(comid) %in% 
                                                  lc_sites_metadata_all$comid)),
                       .progress = TRUE)
```

### Transform the NWM to daily values

``` r
#### Transform all NWM forecasts to daily 

mt_downloads_lc_trimmed_daily <- map(mt_downloads_lc, nwm_daily_maker,
                                     .progress = TRUE)

nwm_daily_forecast_lc_tribs <- bind_rows(mt_downloads_lc_trimmed_daily) %>%
  mutate(comid = as.character(comid))


##### Combine all members and take the mean
##### And transform discharge by normalizing by drainage area
##### And taking the log

nwm_daily_mean_forecast_lc_tribs <- nwm_daily_forecast_lc_tribs %>%
  dplyr::group_by(comid, init_date, predict_date) %>%
  summarise(mean_forecasted_q_cms = mean(forecasted_q_cms)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(comid) %>%
  dplyr::arrange(init_date, .by_group = TRUE) %>%
  dplyr::ungroup() %>%
  mutate(lead_days = as.numeric(predict_date - init_date)) %>%
  inner_join(., lc_gages_metadata_clean %>%
               dplyr::select(tributary, site_no, comid) %>%
               mutate(comid = as.character(comid)),
             by = "comid") %>%
  inner_join(., final_streamstat_features %>%
               dplyr::select(tributary, drnarea_km2),
             by = "tributary") %>%
  mutate(mean_forecasted_q_cms_km2 = mean_forecasted_q_cms/drnarea_km2) %>%
  mutate(mean_forecasted_log_q_cms_km2 = log10(mean_forecasted_q_cms_km2))


################################################################################
```

# Calculate various flow variables

### Antecedent conditions

This is mostly challenging because at longer leadtimes, “antecedent”
conditions become a mixture of forecasted and observed values. So for,
say, a three day lead time weekly antecedent discharge day 1 would have
four days of observed discharge (day 4,5,6,7) and three days of
forecasted discharge (1,2,3); on day 2, three days of observed and four
days of forecasted, and so on and so forth. So we need to build a
dataframe that reflects this

``` r
### First, make sure the observed data is arranged by date and station

observed_flow_for_nwm_predictions <- flow_data_transform %>%
  dplyr::group_by(tributary) %>%
  arrange(date, .by_group = TRUE) %>%
  dplyr::ungroup()

#### Create a lookup table to determine which timesteps 
#### For which we need antecedent values 

timesteps_needed <- nwm_daily_mean_forecast_lc_tribs %>%
      ungroup() %>%
      mutate(days_needed_from_gage_weekly_ant = 7 - lead_days) %>%
      mutate(days_needed_from_gage_weekly_ant = 
               ifelse(days_needed_from_gage_weekly_ant < 0, 0, 
                    days_needed_from_gage_weekly_ant)) %>%
      mutate(days_needed_from_nwm_weekly_ant = 7 - days_needed_from_gage_weekly_ant) %>%
      mutate(days_needed_from_gage_monthly_ant = 30 - lead_days) %>%
      mutate(days_needed_from_gage_monthly_ant = 
               ifelse(days_needed_from_gage_monthly_ant < 0, 0, 
                    days_needed_from_gage_monthly_ant)) %>%
      mutate(days_needed_from_nwm_monthly_ant = (30) - 
               days_needed_from_gage_monthly_ant) %>%
          dplyr::select( comid, lead_days, init_date, predict_date, 
                      days_needed_from_gage_weekly_ant, 
                      days_needed_from_nwm_weekly_ant,
                      days_needed_from_gage_monthly_ant, 
                      days_needed_from_nwm_monthly_ant,
                     ) %>%
  inner_join(., lc_gages_metadata_clean %>%
               mutate(comid = as.character(comid)) %>%
               dplyr::select(site_no, comid),
             by = "comid")
  

#### Find antecedent conditions for all LC gages for each forecasted data
#### in water year 2022-2023

flow_ants_lc <- purrr::pmap_dfr(list(predict_date = timesteps_needed$predict_date,
                          init_date = timesteps_needed$init_date,
                          site_no = timesteps_needed$site_no,
                          days_from_gage =
                            timesteps_needed$days_needed_from_gage_weekly_ant,
                          days_from_model =
                            timesteps_needed$days_needed_from_nwm_weekly_ant,
                monthly_days_from_gage =
                  timesteps_needed$days_needed_from_gage_monthly_ant,
                monthly_days_from_model =
                  timesteps_needed$days_needed_from_nwm_monthly_ant
                          ),
                .f = possibly(antecedent_calculator,
                      otherwise = tibble(mean_prior_weekly_log_q_cms_km2 = NA, 
                                         mean_prior_monthly_log_q_cms_km2 = NA)),
                
                obs_df = observed_flow_for_nwm_predictions,
                model_df = nwm_daily_mean_forecast_lc_tribs,
                .progress = TRUE)
```

### Delta Q

``` r
#### Calculate the delta q

delta_q <- limb_getter(nwm_daily_mean_forecast_lc_tribs,
                         observed_flow_for_nwm_predictions)
```

### Combine antecedent and delta Q

``` r
antecedent_flow <- bind_cols(nwm_daily_mean_forecast_lc_tribs, 
                             flow_ants_lc) %>%
  inner_join(., 
             delta_q,
             by = c("predict_date", "init_date", "tributary")) 
```

### Put all the flow conditions together

``` r
nwm_flow_and_antecedents <- antecedent_flow %>%
  dplyr::select(!delta_daily_q_cat) %>%
  mutate(day_of_year = yday(predict_date)) %>%
  dplyr::mutate(season = case_when(day_of_year %in% seq(60,151,1) ~ "Spring",
                                   day_of_year %in% seq(152,243,1) ~ "Summer",
                                   day_of_year %in% seq(244,334,1) ~ "Fall",
                                   day_of_year > 334 | day_of_year < 60 ~ "Winter")) %>%
  mutate(water_year = add_waterYear(predict_date)) %>%
  dplyr::select(!day_of_year) %>%
  rename(log_daily_q = mean_forecasted_log_q_cms_km2,
         mean_prior_weekly_q = mean_prior_weekly_log_q_cms_km2,
         mean_prior_monthly_q = mean_prior_monthly_log_q_cms_km2,
         ) %>%
  dplyr::select(tributary, 
                init_date, predict_date, lead_days,
                log_daily_q, 
                mean_prior_weekly_q, mean_prior_monthly_q, lag_daily_q,
                delta_daily_q,
                season,
                water_year
                )
```

# Join to watershed characteristics

Here we are going to join the watershed characteristics contained in the
models we selected in 04_model_development. These will be the final
dataframes for performing forecasts

### For gaged

``` r
### Total Phosphorus

forecasting_data_tp_gaged <- nwm_flow_and_antecedents %>%
  inner_join(., watershed_chars,
             by = "tributary") %>%
  dplyr::select(tributary, init_date, predict_date, lead_days, 
                picked_vars_tp$Feature) %>%
  dplyr::select(order(colnames(.)))


#### Subset the chloride_test data to just the variables identified earlier
#### And join to 

forecasting_data_chlor_gaged <- nwm_flow_and_antecedents %>%
  inner_join(., watershed_chars,
             by = "tributary") %>%
  dplyr::select(tributary, init_date, predict_date, lead_days, 
                picked_vars_chlor$Feature) %>%
  dplyr::select(order(colnames(.)))
```

### For ungaged

``` r
### For total phosphorus

forecasting_data_tp_ungaged <- nwm_flow_and_antecedents %>%
  inner_join(., watershed_chars,
             by = "tributary") %>%
  dplyr::select(tributary, init_date, predict_date, lead_days, 
                picked_vars_tp_ungaged$Feature) %>%
  dplyr::select(order(colnames(.))) %>%
  dplyr::group_by(tributary, lead_days) %>%
  nest() %>%
  dplyr::ungroup() %>%
  rename(nwm_data = data) %>%
  inner_join(., tp_test %>%
               dplyr::mutate(trib = tributary) %>%
               dplyr::group_by(trib) %>%
               nest(),
                join_by(tributary == trib)) %>%
  rename(observed_data = data) %>%
  inner_join(., d_fact_tp_ungaged, join_by(tributary == trib))

### For chloride

forecasting_data_chlor_ungaged <- nwm_flow_and_antecedents %>%
  inner_join(., watershed_chars,
             by = "tributary") %>%
  dplyr::select(tributary, init_date, predict_date, lead_days, 
                picked_vars_chlor_ungaged$Feature) %>%
  dplyr::select(order(colnames(.))) %>%
  dplyr::group_by(tributary, lead_days) %>%
  nest() %>%
  dplyr::ungroup() %>%
  rename(nwm_data = data) %>%
  inner_join(., chlor_test %>%
               dplyr::mutate(trib = tributary) %>%
               dplyr::group_by(trib) %>%
               nest(),
                join_by(tributary == trib)) %>%
  rename(observed_data = data) %>%
  inner_join(., d_fact_chlor_ungaged, join_by(tributary == trib))
```
