05_forecast_data_download
================
JTK
2025-01-13

################################################################################ 

This script downloads National Water Model (NWM) hydrological streamflow
forecasts for user-specified sites

It downloads NWM forecasts from an open-source Google Bucket, where
operational forecasts are archived. It uses cloud-based resources to
trim the large files output by the NWM to just our stations of interest.
Users can change the COMIDs to better reflect their locations of
interest. The NWM download relies heavily on the nwmTools package
written by Mike Johnson, but with some rewrites of the core functions,
which I believe were written for past iterations. Please note that this
can take a **VERY LONG** time to run

**Inputs**

1)  `lc_sites_metadata_all` (**from 01_data_discovery_and_download**): a
    list of all the COMIDs for which we are interested in downloading
    NWM forecast data

**Ouputs**

1)  CSV files with forecasts for each user-specified COMID. A single
    file is written for each month of forecasts (for example, forecasts
    for all forecasts issued in Oct 2022 for all COMIDs of interest
    would be ouput in file Oct2022.csv)

################################################################################ 

# Housekeeping

### Packages

``` r
### Data mgmt

require(tidyverse)
require(glue)

### NWM download and NetCDF mgmt

require(terra)
require(tidync)
require(nwmTools)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/02_watershed_attributes_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/03_observational_data_prep.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/04_model_development.Rmd"), quiet=TRUE))
```

# Download NWM Forecasts

### Generate the urls where forecasts might be found

``` r
#### First, get the dates that we are interested in
#### We only want to download forecasts from v2.1 of the NWM
#### Which was operational from from 21 April 2021 to 20 September 2023

download_dates_mt <- tibble(date = seq.Date(from = as.Date("2021-04-21"), 
                                            to = as.Date("2023-09-20"), 
                                            by = "day")) %>%
  .$date

#### Now generate the Google Bucket URLs where each of those forecasts is archived
#### Note that the medium-term forecast is initialized every 6 hours, and forecasts
#### out 204 hrs (8.5 days) in the future
#### We only want to download one daily initialization, which is that made at 
#### midnight UTC
#### There are also seven members of the medium-term NWM
#### We want to download them all

#### URLS for V2.1 ######################################################################
#########################################################################################

#### Make some empty lists to save things

urls_per_member_v21 <- list()

urls_all_days_mt_v21 <- list()

#### Loop over all dates

for(j in 1:length(download_dates_mt)) {
  
  forecast_date_mt <- download_dates_mt[j]
  
  #### Loop over each member
  
  for(i in seq(1, 7, 1)){
    
        #### Print progress
  
        print(paste(forecast_date_mt, "_ ens", i ))
        
        #### Generate URLs using out rewrite of the nwmTools function
      
        urls_v21 <- get_gcp_urls2(domain = "conus",
          output = "channel_rt",
          config = "medium_range",
          ensemble = i,
          date = forecast_date_mt,
          hour = "00",
        minute = "00",
        num = 204)
        
        #### Clean-up the urls
        
        urls_v21 <- urls_v21 %>%
          #mutate(urls = sub("(.*)f ", "\\1f0\\2", urls)) %>%
          mutate(init_time = 0) %>%
          mutate(init_date = forecast_date_mt) %>%
          mutate(member = paste0("mem", i))
        
        #### Store each member for a given date
        
        urls_per_member_v21[[i]] <- urls_v21
    
  }
  
  #### Store all members for each date
  
  urls_all_days_mt_v21[[j]] <- urls_per_member_v21
}





### Join all the urls together

all_urls_mt <- bind_rows(urls_all_days_mt_v21) %>%
  mutate(predict_date = date(dateTime)) %>%
  mutate(lead_time = str_extract(str_extract(urls, "f\\d+"), "\\d+")) %>%
  mutate(init_date_time = as_datetime(init_date)) 



#############################################################################################
#############################################################################################
```

### Format URLs and COMID dataframes to allow for download

``` r
#### First, nest the URL dataframe

all_urls_and_files_mt_nest <- all_urls_mt %>%
  rename(predict_dateTime = dateTime) %>%
  filter(lead_time <= 192) %>% ### Trim off the half-day at the end of each forecast
  mutate(membr = member ) %>%
  nest_by(init_date, membr) %>%
  ungroup() %>%
  mutate(month = format(as.Date(init_date), "%b%Y"))

#### Now, subset the files to download
#### We do this so that we can iteratively save things because downloading 
#### the full years-long record would take forever
#### This allows us to track our progress better 

#### Subset of files to download

#*****************Change this to download month of interest**********************#

trim_urls_and_files_mt_nest <- all_urls_and_files_mt_nest %>%
  filter(month == "Oct2022") 

#********************************************************************************#


#### Unnest the dataframe

trim_urls_and_files_mt_unnest <- trim_urls_and_files_mt_nest %>%
  unnest(data) %>%
  dplyr::select(!membr) 

#### Declare where to write things 

write_file <- paste0(here("nwm_operational/medium_term"), "/",
                     trim_urls_and_files_mt_unnest$month[1],
                     ".csv")
```

### Download NWM

``` r
#### Extract the COMIDs that we are interested in
#### Note that this inherits the lc_sites_metadata_all variable from
#### 01_data_discovery_and_download.Rmd (so it requires that script to be run)

stations <- lc_sites_metadata_all$comid 


#### Do the download
#### Do it
#### Note that this function tries a given URL five times maximum with a 
#### delay of 60 seconds between tries
#### This allows minor blips in internet connection to not break the download
#### If none of those tries returns actual data, it then returns a dataframe
#### with NAs that is formatted (i.e., contains all the fields) of the actual
#### data

walk(trim_urls_and_files_mt_unnest$urls,
     .f = possibly(insistently(get_timeseries3,
                               rate = rate_delay(pause  = 60,
                                                 max_times = 5)),
                   otherwise = tibble(comid = NA, 
                                      init_date = NA, 
                                      lead_time = NA,
                                      member = NA,
                                      predict_dateTime = NA,
                                      modeled_q_cms = NA)),
     stations, ### Constant fed to .f 
     write_file,
     .progress = TRUE)
```

### Check to see if download completed

``` r
#### Since we are downloading a month at a time, the download may break before
#### The entire month is downloaded. To check, we want to see what the last date
#### downloaded was
#### If the entire month was downloaded it should be, for example, 2022-10-31

### Read in the file
nwm_forecasts <- read_csv(write_file)

### See if we downloaded the whole month
tail(nwm_forecasts)

### Check the last date
last_date <- last(nwm_forecasts)

last_date

### Find index of last timestep downloaded in the original df
last_downloaded_index <- which(trim_urls_and_files_mt_unnest$init_date == 
                                 last_date$init_date &  
                               trim_urls_and_files_mt_unnest$predict_dateTime == 
                                 last_date$predict_dateTime &
                               trim_urls_and_files_mt_unnest$member ==
                                 last_date$member)

last_downloaded_index 


### Check for missing
missing_dates_indexes <- which(is.na(nwm_forecasts$predict_dateTime))

missing_dates_indexes

### If they exist, where are they

missing_data <- nwm_forecasts %>%
  mutate(date_mem = paste(init_date, member, sep = "_")) %>%

unique(missing_data$date_mem)

### Check all days are present
unique(nwm_forecasts$init_date)

unique(nwm_forecasts$member)

### Check the time zone
nwm_forecasts$predict_dateTime[1]
```
