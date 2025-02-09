01_data_discovery_and_download
================
JTK
2025-01-03

################################################################################ 

This script downloads flow data from the pre-selected gages of interest.
It also discovers stations with water quality data gathered by Vermont
DEC corresponding to each of the streamflow gages and downloads that
data. We have specifically tailored this approach to the Lake Champlain
basin, but it is flexible enough so that interested users can modify
various lines to download from the gages and sampling sites of interest
to them.

**Required inputs**

1)  NWIS site numbers for sites of interest
    (data/lake_champlain_usgs_gages.csv)

2)  Missing data from Little Otter Creek(data/little_otter_tp.csv &
    data/little_otter_chlor.csv)

\*\*\*\* Note that, for adapting to other sites, only \#1 would be
required

**Outputs/Returns**

1)  Raw dataframe of daily streamflow data for 18 watersheds in Lake
    Champlain

2)  Raw water quality dataframe for total phosphorus and chloride
    concentration as measured in each of those 18 watersheds

################################################################################ 

# Housekeeping

### Packages

``` r
## Data mgmt
require(tidyverse)

## Data download
require(dataRetrieval)
require(nhdplusTools)
require(EPATADA)

## Misc.
require(here)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))
```

# Get metadata for USGS gage of interest

``` r
#### First, retrieve the NWIS IDs of the gages we have selected
#### in the Lake Champlain basin
#### We've pulled these gage IDs from various pubs in the region 
#### (Underwood et al., 2017, WRR; Vaughn, 2019 LCBP Report)

lc_site_ids <- read_csv(here("data/lake_champlain_usgs_gages.csv")) %>%
  dplyr::select(river_basin, site_no)

#### Make sure the NWIS code is eight digits
#### .csv files often remove the leading zeros from ids

lc_site_ids <- lc_site_ids %>%
  mutate(site_no = as.character(sprintf("%08d", site_no)))

#### Get official site names and various other pieces of metadata

##### Specific to gage location and flow monitoring record

lc_gages_metadata <- whatNWISdata(siteNumber = lc_flow_gages$site_no,
             parameterCd = "00060",
             service = "dv") 

##### Drainage area

lc_gages_metadata_da <- readNWISsite(lc_site_ids$site_no) %>%
  dplyr::select(site_no, drain_area_va)

##### And NHDPlus (Medium-res) COMID
##### Which is absolutely essential for working with the NWM
##### To do this we must first get the station IDs into the format needed
##### To query the NHD for the COMID related to each site location

lc_gages_nldi <- lc_site_ids %>%
  mutate(fsrc = "nwissite",
         fid = paste0("USGS-", site_no)) %>%
  mutate(comid = map2_dbl(fsrc, fid,
                      ~discover_nhdplus_id(nldi_feature = list(featureSource = .x, 
                                                featureID = .y)))) %>%
  dplyr::select(!c("river_basin", "fsrc"))
  

#### And clean up the metadata

lc_gages_metadata_clean <- inner_join(lc_site_ids, 
                                      lc_gages_nldi,
                                      by = "site_no") %>%
  inner_join(.,   lc_gages_metadata %>%
               dplyr::select(station_nm, site_no, 
                             dec_lat_va, dec_long_va, dec_coord_datum_cd,
                             begin_date, end_date),
             by = "site_no") %>%
  inner_join(., lc_gages_metadata_da,
             by = "site_no") %>%
  rename(drain_area_mi2 = drain_area_va) %>%
  mutate(drain_area_km2 = drain_area_mi2*2.58999) %>%
  dplyr::select(!drain_area_mi2) %>%
  relocate(drain_area_km2, .after = "station_nm") %>%
  rename(tributary = river_basin) %>%
  relocate(station_nm, .after = "tributary")
```

# Get metadata for water quality monitoring sites

``` r
#### Now, read-in a list of the VTDEC tributary monitoring sites. 
#### To do so, we have to download data monitored by VTDEC 
#### for a subset of the Champlain Tribs water quality monitoring period (1991-present). 
#### Let's take a slice of 2012
#### Otherwise it would take waaaaaay too long to download
#### 2012, according to VTDEC's webpage, should be a period when monitoring in all tribs
#### is active

#### Note that this returns every sample sampled by VTDEC in this slice of 2012

wq_data_profile <- EPATADA::TADA_DataRetrieval(
                           organization = "1VTDECWQ",
                           startDate = "2012-03-01",
                           endDate = "2012-09-30",
                           applyautoclean = FALSE)

#### Map locations of samples and inspect to make sure they fall in both
#### VT & NY and include all the tribs we want

EPATADA::TADA_OverviewMap(wq_data_profile %>%
                         rename(TADA.LatitudeMeasure = ActivityLocation.LatitudeMeasure,
                                TADA.LongitudeMeasure = ActivityLocation.LongitudeMeasure,
                                TADA.CharacteristicName = CharacteristicName) %>%
                         mutate(TADA.LatitudeMeasure = as.numeric(TADA.LatitudeMeasure),
                                TADA.LongitudeMeasure =as.numeric(TADA.LongitudeMeasure)))

#### Select river monitoring sites that fall along the main stem of each tributary
#### Then, slice by the number of samples in our sample monitoring period
#### We will select the site with the most samples, which we know will correspond
#### to the main monitoring site for water quality for each of the tribs
#### (which is the site we want)
#### We will then use these site ids to download all the water quality data for the monitoring period

site_ids_wq <- wq_data_profile %>%
  filter(MonitoringLocationName %in% lc_gages_metadata_clean$tributary) %>%
  dplyr::select(MonitoringLocationIdentifier, MonitoringLocationName, ResultIdentifier) %>%
  dplyr::group_by(MonitoringLocationIdentifier,MonitoringLocationName) %>%
  summarise(sample_count = length(unique(ResultIdentifier))) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(MonitoringLocationName) %>%
  slice_max(sample_count)

#### Now plot again to make sure that we selected the right ones

EPATADA::TADA_OverviewMap(wq_data_profile %>%
                         rename(TADA.LatitudeMeasure = ActivityLocation.LatitudeMeasure,
                                TADA.LongitudeMeasure = ActivityLocation.LongitudeMeasure,
                                TADA.CharacteristicName = CharacteristicName) %>%
                         mutate(TADA.LatitudeMeasure = as.numeric(TADA.LatitudeMeasure),
                                TADA.LongitudeMeasure =as.numeric(TADA.LongitudeMeasure)) %>%
                         filter(MonitoringLocationIdentifier %in% site_ids_wq$MonitoringLocationIdentifier))

#### And finally join the stream gage metadata file so that we have all
#### station identifiers in one location

lc_sites_metadata_all <- inner_join(lc_gages_metadata_clean,
           site_ids_wq %>%
             dplyr::select(MonitoringLocationIdentifier, MonitoringLocationName) %>%
             rename(wq_site_id = MonitoringLocationIdentifier,
                    tributary = MonitoringLocationName),
           by = "tributary") %>%
  relocate(wq_site_id, .after = fid)


################################

#### Remove extraneous variables

rm(wq_data_profile)

################################
```

# Download Data

### Streamflow Data

``` r
#### Downloads USGS flow data 
#### We want to daily values
#### so let's download them

flow_data <- dataRetrieval::readNWISdv(siteNumbers = lc_sites_metadata_all$site_no, 
                          parameterCd = "00060",
                          startDate = "1990-01-01",
                          endDate = "2023-12-31") %>%
  renameNWISColumns() %>%
  addWaterYear() %>%
  dplyr::filter(!str_detect(Flow_cd , "P"))
```

### Water quality data

``` r
#### Download water quality data from each of the eighteen tribs
#### This require looping over each station to avoid breaking the downloader

tribs_wq <- list()

for(i in 1:length(lc_sites_metadata_all$wq_site_id)) {
  
  cat(crayon::cyan("Reading", lc_sites_metadata_all$tributary[i], "\n"))
  
  wq_by_site <- TADA_DataRetrieval(siteid = lc_sites_metadata_all$wq_site_id[i],
                           startDate = "1990-01-01",
                           endDate = "2023-12-31",
                           applyautoclean = FALSE)
  
  
  tribs_wq[[i]] <- wq_by_site
  
}



lc_tribs_wq_all <- bind_rows(tribs_wq)

#### We also need to manually bring in some data for Little Otter Creek
#### For some reason, data from 2020-2023 are not available on STORET for Little
#### Otter Creek, so we have had to download it manually from 
#### https://dec.vermont.gov/watershed/lakes-ponds/monitor/lake-champlain-long-term-monitoring-project#Data
#### We have saved these locally and now import them here

little_otter_total_phos <- read_csv(here("data/little_otter_tp.csv"))

little_otter_chlor <- read_csv(here("data/little_otter_chlor.csv"))

#### Bind together
little_otter_missing <- bind_rows(little_otter_total_phos, little_otter_chlor)
```
