03_observational_data_prep
================
JTK
2025-01-07

################################################################################ 

This script manipulates and cleans the observational data (discharge,
water quality) in order to get it reading for model construction. It
also splits data into training, validation, and testing sets.

################################################################################ 

# Housekeeping

### Packages

``` r
### Data mgmt
require(tidyverse)
require(zoo)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/02_watershed_attributes_download.Rmd"), quiet=TRUE))
```

# Data Cleaning & Prep

### Clean the USGS observed discharge data

``` r
#### Join to tables with site names
#### Transform flow data to m3/s 
#### Take 0 flow days and make them instead a very small number (10^-4 cfs)
#### And finally normalize by watershed area

flow_data_clean <- flow_data %>%
  dplyr::select(site_no, dateTime, Flow, waterYear) %>%
  rename(date = dateTime,
         discharge_cfs = Flow) %>%
  mutate(discharge_cfs = ifelse(discharge_cfs == 0, 1E-4, discharge_cfs)) %>%
  mutate(discharge_cms = discharge_cfs*0.0283168) %>%
  inner_join(lc_gages_metadata_clean %>%
               dplyr::select(tributary, site_no),
             .,
             by = "site_no") %>%
  dplyr::ungroup() %>%
  inner_join(., final_streamstat_features %>%
              dplyr::select(tributary, drnarea_km2),
             by = "tributary") %>%
  mutate(discharge_cms_km2 = discharge_cms/(drnarea_km2)) %>%
  dplyr::select(tributary, site_no,
                date, waterYear, 
                discharge_cms_km2)
```

### Transform to log and add in antecedent features

``` r
#### First, transform normalized discharge (m3/s/km2) to log-scale
#### Then, calculate daily, weekly, and monthly discharge
#### Calculate the the change in discharge from time t to time t-1
#### Calculate that change normalized by discharge at time t-1
#### And then also transform that into a categorical "rise/fall"
#### based on whether it is positive (rise) or negative (fall)
#### And finally determine a "season" categorical variable


flow_data_transform <- flow_data_clean %>% 
  group_by(tributary) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(log_daily_q = log10(discharge_cms_km2)) %>%
  mutate(mean_prior_weekly_q = rollapply(log_daily_q,
                                         width = list(-(7:1)),
                                         FUN = mean, align = "right",
                                         fill= NA)) %>%
  mutate(mean_prior_monthly_q = rollapply(log_daily_q,
                                         width = list(-(30:1)),
                                         FUN = mean, align = "right",
                                         fill= NA)) %>%
  mutate(delta_daily_q = discharge_cms_km2 - lag(discharge_cms_km2)) %>%
  mutate(lag_daily_q = lag(log_daily_q)) %>%
  mutate(dq_q = abs(delta_daily_q/lag(discharge_cms_km2))*sign(delta_daily_q)) %>%
  mutate(delta_daily_q_cat = case_when(dq_q > 0.10 ~ 1,
                                       dq_q < -0.10 ~ -1,
                                       (dq_q <= 0.10 & dq_q >= -0.10) ~ 0)) %>%
  mutate(delta_daily_q = log_daily_q - lag(log_daily_q)) %>% ### Update
  mutate(delta_daily_q_cat = as.factor(delta_daily_q_cat)) %>%
  mutate(day_of_year = yday(date)) %>%
  dplyr::slice(-(1:30)) %>% ### Remove first thirty days of record 
  dplyr::select(tributary, site_no, 
                date, waterYear, 
                day_of_year,
                log_daily_q,
                mean_prior_weekly_q, 
                mean_prior_monthly_q,
                delta_daily_q, 
                delta_daily_q_cat,
                lag_daily_q) %>%
  rename(water_year = waterYear) %>% ### Better follows our naming conventions
  dplyr::ungroup() %>%
  tidyr::drop_na(mean_prior_monthly_q) %>%
  tidyr::drop_na(tributary) %>%
  mutate(date = as_date(date)) %>%
  dplyr::mutate(season = case_when(day_of_year %in% seq(60,151,1) ~ "Spring",
                                   day_of_year %in% seq(152,243,1) ~ "Summer",
                                   day_of_year %in% seq(244,334,1) ~ "Fall",
                                   day_of_year > 334 | day_of_year < 60 ~ "Winter")) %>%
  dplyr::select(!day_of_year) ### Remove day of year
```

# Water quality data

### Clean VTDEC water quality monitoring data

``` r
#### Mostly, select the constituents we want, rename, and join columns to get 
#### a more workable format
#### The most important thing we are doing here is joining the constituent with 
#### its measured fraction to make one variable 
#### so we can distinguish between, say, total and dissolved phosphorus
#### We are also dropping missing data, making sure concentrations are numeric
#### and transforming dates to datetime format
#### Also, if there are multiple observations on the same date at the same site,
#### we take the mean to get one observation per site per date
#### We then want to add in the missing data from Little Otter Creek that we 
#### had to download manually from VTDEC's website
#### We then want to make sure Chloride, which for some reason is sometimes listed 
#### as total OR dissolved, is all under one moniker
#### (Because there is no difference between the two)
#### And finally select only Chloride and Total Phosphorus

lc_tribs_wq_all_clean <- lc_tribs_wq_all %>%
  dplyr::select(where(~!all(is.na(.x)))) %>% ### Drop the columns that have all NAs
  dplyr::select( 
                MonitoringLocationName, MonitoringLocationIdentifier,
                ActivityLocation.LatitudeMeasure, ActivityLocation.LongitudeMeasure,
                CharacteristicName, ResultSampleFractionText,
                ActivityStartDate, ActivityStartTime.Time,
                ResultMeasureValue, ResultMeasure.MeasureUnitCode) %>%
  rename(tributary = MonitoringLocationName,
         storet_site_id = MonitoringLocationIdentifier, 
         latitude = ActivityLocation.LatitudeMeasure,
         longitude = ActivityLocation.LongitudeMeasure,
         constituent = CharacteristicName,
         fraction = ResultSampleFractionText,
         date = ActivityStartDate,
         time = ActivityStartTime.Time,
         conc = ResultMeasureValue,
         units = ResultMeasure.MeasureUnitCode) %>%
  mutate(constituent = ifelse(!is.na(fraction), 
                              paste0(constituent, "_", fraction),
                              constituent)) %>%
  dplyr::select(!fraction) %>%
  as_tibble() %>%
  mutate(conc = as.numeric(conc)) %>%
  drop_na(conc) %>%
  mutate(date = as_date(date)) %>%
  dplyr::group_by(tributary, date, constituent) %>%
  mutate(conc = mean(conc)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::select(tributary, date, constituent, conc) %>%
  bind_rows(., little_otter_missing %>%
              dplyr::select(!time)) %>%
  mutate(log_conc = log10(conc)) %>% ### Transform to log
  mutate(constituent = ifelse(str_detect(constituent, "Chloride"), "Chloride", constituent)) %>%
  dplyr::filter(constituent %in% c("Phosphorus_Total", "Chloride")) %>%
  dplyr::select(!conc) ### Remove raw concentration
```

# Combine datasets

### Join discharge and water quality data

``` r
#### Join them together

tribs_wq_and_q_all <- inner_join(flow_data_transform, 
                                 lc_tribs_wq_all_clean,
                                 by = c("tributary", "date")) %>%
  dplyr::ungroup()
```

### Join combined water quality & discharge data to watershed attributes

``` r
#### Do it 

all_drivers <- tribs_wq_and_q_all %>%
  inner_join(., watershed_chars,
             by = "tributary") %>%
  filter(water_year < 2024) ### Remove data from water year 2024

#### Split into Total Phosphorus and Chloride dataframes
#### Also add a numerical ID field that uniquely represents each tributary
#### This will make it easier to make certain model train/valid/splits

##### For total phosphorus

tp_drivers <- all_drivers %>%
  filter(constituent == "Phosphorus_Total") %>%
  group_by(tributary) %>%
  mutate(group_id = cur_group_id()) %>%
  dplyr::ungroup() 

##### And for chloride

chlor_drivers <- all_drivers %>%
  filter(constituent == "Chloride") %>%
  group_by(tributary) %>%
  mutate(group_id = cur_group_id()) %>%
  dplyr::ungroup() 
```

# Make final modeling dataframes

### Remove unsplittable features

``` r
#### We want check to make sure each feature has at least two unique values across
#### all cross validation data sets 
#### This removes features that would have limited explainability in a leave-one-out
#### cross validation scenario and might instead "identify" a specific basin 
#### Rather than reflect, in some way, physical/chemical process(es)
#### We've written a small function to do so

#### Find the unsplittable features

##### This can really be any of the constituent driver dataframes
##### because these are likely to be static attributes
##### but just to be safe lets do it for both

unsplittable_feats_tp <- remove_unsplittable(tp_drivers)

unsplittable_feats_chlor <- remove_unsplittable(chlor_drivers)

#### Remove those from the drivers data frame

tp_drivers <- tp_drivers %>%
  dplyr::select(!unsplittable_feats_tp)

chlor_drivers <- chlor_drivers %>%
  dplyr::select(!(unsplittable_feats_chlor))

#### Check to make sure it worked 

unsplittable_feats_tp %in% names(tp_drivers)

unsplittable_feats_chlor %in% names(chlor_drivers)
```

### Split into training/valid, testing

``` r
#### For Total Phosphorus

tp_train_valid <- tp_drivers %>%
  filter(water_year < 2019)

tp_test <- tp_drivers %>%
  filter(water_year >= 2019)

#### For Chloride

chlor_train_valid <- chlor_drivers %>%
  filter(water_year < 2019)

chlor_test <- chlor_drivers %>%
  filter(water_year >= 2019)
```
