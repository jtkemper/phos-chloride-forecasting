02_watershed_attributes_download
================
JTK
2025-01-04

################################################################################ 

This script downloads various static watershed attributes for each
individual watershed in the Lake Champlain Basin. It draws from a
variety of sources, including the National Hydrography Dataset (high-
and medium-res), US (SSURGO) and Canadian soils datasets, USGS
StreamStats, a USGS-built set of expanded for the NHD (Wieczorek et al.,
2018, <https://doi.org/10.5066/F7765D7V>.), and several other publically
available datasets.

We then compile these into one large dataframe that contains watershed
attributes for each basin. The goal in doing this is to develop “global”
machine learning models that may potentially learn relationships between
dynamic hydrology and static watershed attributes, allowing them to
learn from a diversity of data (i.e., observations in all watersheds) to
make predictions in individual watersheds

**Inputs**

1)  Shapefile that contains outlet points for each watershed

2)  .csv of land cover data derived from Troy et al., 2007

**Outputs**

1)  Wide dataframe of static watershed attributes for each of the 18
    Lake Champlain tributaries

################################################################################ 

# Housekeeping

### Packages

``` r
### Data mgmt 
require(tidyverse)
```

    ## Loading required package: tidyverse

    ## Warning: package 'ggplot2' was built under R version 4.3.2

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.2     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.0
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.2     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
require(tsibble)
```

    ## Loading required package: tsibble
    ## 
    ## Attaching package: 'tsibble'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     interval
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, union

``` r
### Hydrography & watershed characteristics
require(streamstats)
```

    ## Loading required package: streamstats
    ## The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
    ## which was just loaded, will retire in October 2023.
    ## Please refer to R-spatial evolution reports for details, especially
    ## https://r-spatial.org/r/2023/05/15/evolution4.html.
    ## It may be desirable to make the sf package available;
    ## package maintainers should consider adding sf to Suggests:.
    ## The sp package is now running under evolution status 2
    ##      (status 2 uses the sf package in place of rgdal)
    ## Please note that rgdal will be retired during October 2023,
    ## plan transition to sf/stars/terra functions using GDAL and PROJ
    ## at your earliest convenience.
    ## See https://r-spatial.org/r/2023/05/15/evolution4.html and https://github.com/r-spatial/evolution
    ## rgdal: version: 1.6-7, (SVN revision 1203)
    ## Geospatial Data Abstraction Library extensions to R successfully loaded
    ## Loaded GDAL runtime: GDAL 3.6.2, released 2023/01/02
    ## Path to GDAL shared files: C:/Users/jkemper/AppData/Local/R/win-library/4.3/rgdal/gdal
    ##  GDAL does not use iconv for recoding strings.
    ## GDAL binary built with GEOS: TRUE 
    ## Loaded PROJ runtime: Rel. 9.2.0, March 1st, 2023, [PJ_VERSION: 920]
    ## Path to PROJ shared files: C:/Users/jkemper/AppData/Local/R/win-library/4.3/rgdal/proj
    ## PROJ CDN enabled: FALSE
    ## Linking to sp version:1.6-1
    ## To mute warnings of possible GDAL/OSR exportToProj4() degradation,
    ## use options("rgdal_show_exportToProj4_warnings"="none") before loading sp or rgdal.

``` r
require(nhdplusTools)
```

    ## Loading required package: nhdplusTools

``` r
### Spatial 
require(sf)
```

    ## Loading required package: sf
    ## Linking to GEOS 3.11.2, GDAL 3.6.2, PROJ 9.2.0; sf_use_s2() is TRUE

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))
```

# Get location data

### Import watershed outlet points

These points were delienated manually in ArcGIS by selecting the
flowline closest to the watershed outlet and transforming it into a
point feature

``` r
#### Import 'em

start_points <- st_read(here("data/watershed_outlets/watershed_outlets_updated.shp")) %>%
  st_zm(drop = TRUE) %>%
  filter(tributary != "Putnam Creek")
```

### Find COMID (for NHD Medium-Res) for outlet points

``` r
#### Convert the simple feature collection to a set of point features 

start_points2 <- start_points["geometry"] %>% sf::st_as_sfc()

#### Create a list to save the medium-res comids

comid_start <- list()

#### Find COMIDs corresponding to outlet points

for(i in 1:length(start_points2)) {
  
  print(i)
  
   start <- tibble(start = discover_nhdplus_id(start_points2[i]))
    
  
  comid_start[i] <- start
  
}

#### Bind them all together into one dataframe

start_comid_df <- tibble(start_comid = unlist(comid_start)) %>%
  bind_cols(start_points, .) %>%
  as_tibble() %>%
  dplyr::select(!geometry)
```

### Get all NHD Medium-Res comids in a particular basin

This will be necessary if we want to use available USGS datasets to
determine various attributes at each flowline in the entire basin

``` r
#### Get all COMIDs in an individual basin

comids_by_basin <- map2(start_comid_df$tributary,
                        start_comid_df$start_comid,
                        all_comid_getter)

#### Bind them all together

comids_by_basin <- comids_by_basin %>%
  bind_rows() %>%
  unnest(data)
```

### Find outlet COMIDs for NHD high-resolution

``` r
#### First, download the NHD high-resolution

lc_nhd_hr <- download_nhdplushr(here("downloads/nhd_hr"),
                                hu_list = "0430")

#### Get the name of the geodatabase

lc_nhd_hr_gdb <- list.files(lc_nhd_hr, pattern = ".gdb")

#### And create a geopackage to save the data 

lc_hr_gpkg <- here("downloads/nhd_hr/0430_hr.gpkg")

#### Extract the NHD data from the downloaded geodatabase

lc_hr <- get_nhdplushr(lc_nhd_hr,
                       out_gpkg = lc_hr_gpkg)


#### Find high-resolution COMIDs

comid_start_hr <- list()

for(i in 1:length(start_points2)) {
  
  print(i)
  
   start <- get_flowline_index(lc_hr$NHDFlowline,
                   start_points2[i])
    
  
  comid_start_hr[i] <- tibble(COMID = start$COMID)
  
}

#### Bind them all together into one dataframe

start_comid_hr <- tibble(start_comid = unlist(comid_start_hr)) %>%
  bind_cols(start_points, .) %>%
  as_tibble() %>%
  dplyr::select(!geometry)
```

# Get StreamStats

We are downloading the StreamStats parameters available for VT and NY.
In particular, what we are interested in is percent elevation over 1200
ft, which for the LCB is a rough indicator of snow coverage/importance.

``` r
#### Take the shapefile of catchment outlet points 
#### And transform it to the format needed to access the StreamStats API
#### (which requires an individual field for lat and for long)
start_points_df <- start_points %>%
  as_tibble() %>%
  mutate(lon = st_coordinates(geometry)[,1],
         lat = st_coordinates(geometry)[,2])

#### Create an empty list to save returned data

streamstat_features <- list()

#### A loop to extract the workspace IDs for each watershed based on 
#### the outlet points we have imported
#### The workspace IDs, I believe, are the internal StreamStats reference
#### to the watershed that is delineated from any given point
#### Once we have these, we can acquire the watershed features that we want

for(i in 1:length(start_points_df$tributary)) {
  
  #### Track our progress across watersheds
  
  cat(crayon::cyan("\nRetrieving", start_points_df$tributary[i], "\n"))
  
  #### Delineate the watershed in StreamStats and return watershed characteristics
  #### but don't return the actual watershed polygon
  
  #### Do this in a while loop to repeat download if we haven't actually retrieved 
  #### watershed characteristic values
  #### For some reason, sometimes StreamStats returns watershed 
  #### chars without associated values
  #### When this happens, there are less than six columns in the returned
  #### "parameters" dataframe
  #### This while loop simply checks to see how many columns are in the 
  #### parameters dataframe and, if it's less than six (so no values),
  #### repeats the download of watershed characteristics
  
  n_params <- 5
  
  while(n_params < 6){
    
    #### Do the download
    
    watershed <- delineateWatershed(start_points_df$lon[i],
                                 start_points_df$lat[i],
                                 crs = 4269,
                                 includeparameters = "true",
                                 includefeatures = "false")
  
    n_params <- ncol(watershed$parameters)
    
    #### Print to output if we need to repeat 
    
    if(n_params<6){cat(crayon::yellow("\nRepeating\n"))}
   
  } ## End while loop
  

  #### Nest the parameters and append which watershed they are for
  
  streamstat_features[[i]] <- watershed$parameters %>%
    nest() %>%
    mutate(tributary = start_points_df$tributary[i])
  


} ## End for loop

#### Bind all watersheds together 

all_streamstat_features <- bind_rows(streamstat_features)

#### Now trim the features down to just the two that we want

final_streamstat_features <- all_streamstat_features %>%
  unnest(cols = c(data)) %>%
  filter(str_detect(name, "1200 ft|Drainage")) %>%
  dplyr::select(tributary, code, value) %>%
  mutate(code = tolower(code)) %>%
  pivot_wider(names_from = code, values_from = value) %>%
  mutate(drnarea_km2 = drnarea*2.58999) %>% ### transform from mi^2 to km^2
  dplyr::select(!drnarea)
```

# Get land use attributes

We are importing percentage values for land use that we calculated in
ArcGIS. The land use layer we used is one specific to the LCB that was
created to model phosphorus contributions to the lake. You can read more
about its creation in Troy et al., 2007
(<http://www.lcbp.org/techreportPDF/54_LULC-Phosphorus_2007.pdf>) and
you can download the layer here
(<https://www.arcgis.com/home/item.html?id=16043a36e8a64aa79cb1728cf7d98409>)

We determined the percentage of each land cover class by simply totaling
up the area of the pixels in each watershed and dividing by the
watershed area

Similarly, we created a 100 m buffer around all the flowlines in the
basin and totalled up the area for each landcover class within this
buffer, and then divided by total watershed area. This gives us an idea
of riparian land cover.

We import both these datasets here

``` r
#### Get land cover for each basin entirely

lulc_2001 <- lulc_processer(here("data/lulc_2001_by_watershed.csv"),
                            final_streamstat_features) %>%
  drop_na()

#### Get land cover for riparian areas

lulc_all_100m <- lulc_processer(here("data/lulc_2001_by_watershed.csv"),
                            final_streamstat_features) %>%
  drop_na() %>%
  rename_if(is.numeric, ~paste0("all_100m_", .))

#### Join together

lulc_all <- inner_join(lulc_2001, lulc_all_100m, by = "tributary")
```

# Get some hydrography attributes

Here, we want to get a few things from the NHD HR that may impact total
phosphorus and chloride contributions. Namely, these are things like
drainage density (which we derive from flowline lengths), stream slopes,
lengths of different geomorphic stream types, and frequency of
particular stream orders.

``` r
#### Get some stats

flowline_stats_hr <- map2(start_comid_hr$tributary,
                                start_comid_hr$start_comid,
                                .f = hr_flow_length_calculator,
                                hr_data = lc_hr$NHDFlowline,
                                hr_geopack = lc_hr_gpkg,
                                "Trib",
                                .progress = TRUE)

#### Extract the total flowline length in each tributary

flow_length_hr <- bind_cols(start_comid_hr,
                            flowline_length_km = map_dbl(flowline_stats_hr, 
                                                         1))

#### Extract the percentage (by length) of each stream order 

pct_orders <- bind_cols(start_comid_hr, map_dfr(flowline_stats_hr, ~(.[[3]] %>%
                               as_tibble() %>%
                               nest()))) %>%
  unnest(data) %>%
  pivot_wider(names_from = StreamOrde, names_glue = "pct_{StreamOrde}_order",
              values_from = pct_order) %>%
  dplyr::select(!start_comid) %>%
  replace(is.na(.), 0) 

#### Extract the percentage (by length) of each geomorphic stream type
#### Channel type is based on slope thresholds and calculated based on those listed in 
#### Geomorphic Classification of Rivers: An Updated Review
#### (https://www.fs.usda.gov/rm/pubs_journals/2022/rmrs_2022_buffington_j001.pdf))

pct_by_type <- bind_cols(start_comid_hr, map_dfr(flowline_stats_hr, ~(.[[9]] %>%
                               as_tibble() %>%
                               nest()))) %>%
  unnest(data) %>%
  pivot_wider(names_from = channel_type, names_glue = "{channel_type}_pct",
              values_from = pct_by_type) %>%
  dplyr::select(!start_comid) %>%
  replace(is.na(.), 0) 


#### And calculate drainage density based on the total flowline length
#### And the basin area

basin_drainage_density <- full_join(flow_length_hr %>%
                                      dplyr::select(!start_comid),
                                    final_streamstat_features %>%
                                      dplyr::select(tributary, drnarea_km2),
                                    by = "tributary") %>%
  mutate(drain_density_km_km2 = flowline_length_km/drnarea_km2) %>%
  dplyr::select(tributary, drain_density_km_km2)


#### And now get the length of the mainstem for each tributary

##### Calculate it

mainstem_flowline_lengths <- map2(start_comid_hr$tributary,
                                start_comid_hr$start_comid,
                                hr_flow_length_calculator,
                                hr_data = lc_hr$NHDFlowline,
                                hr_geopack = lc_hr_gpkg,
                                "Main",
                                .progress = TRUE)

##### And bind together

mainstem_flowline_lengths <- bind_cols(start_comid_hr %>%
                                          dplyr::select(tributary),
                                       mainstem_length_km = map_dbl(
                                         mainstem_flowline_lengths, 1)) 
```

# Manually calculate some other watershed characteristics

### Hack’s exponent

``` r
hacks <- inner_join(final_streamstat_features %>%
             dplyr::select(tributary, drnarea_km2),
           mainstem_flowline_lengths %>%
             dplyr::select(tributary, mainstem_length_km), 
           by = "tributary") %>%
  mutate(across(where(is.numeric), ~log10(.))) %>%
  rename_with(~paste0("log_", .), where(is.numeric)) %>%
  mutate(h = log_mainstem_length_km/log_drnarea_km2) %>%
  dplyr::select(tributary, h)
```

### Richard-Baker Flashiness Index

Calculates the Richards-Baker Flashiness Index, an indicator of
hydrologic “flashiness”. More details can be found in the original
publication (<https://doi.org/10.1111/j.1752-1688.2004.tb01046.x>)

``` r
#### Make the flow data into a tsibble to check to see if there are 
#### gaps in the flow data

flow_ts <- flow_data %>%
  as_tibble() %>%
  mutate(dateTime = as_date(dateTime)) %>%
    as_tsibble(key = site_no, index = dateTime)



#### Check for gaps

time_gaps <- has_gaps(flow_ts) %>%
  filter(.gaps == TRUE)


#### Extract the sites with gaps

sites_with_gaps <- flow_ts %>%
  filter(site_no %in% time_gaps$site_no)

#### Find where the gaps are

gaps <- count_gaps(flow_ts)

#### If gaps are small
#### We can fill the gaps with linear interpolation

flow_filled_gaps <- flow_ts %>%
  fill_gaps()

#### However, if they are large, we can just trim the flow data
#### to only a period without gaps
#### And assume the flashiness has not changed too drastically outside that period
flow_ts <- flow_ts %>%
  filter(waterYear < 2015)




#### Check to see that the number of gaps we filled with NAs
#### is equal to the number of gaps there are

identical(nrow(flow_filled_gaps %>%
                   filter(is.na(Flow))),
            sum(gaps$.n))


#### Calculate the Richards-Baker Flashiness Index
#### From https://doi.org/10.1111/j.1752-1688.2004.tb01046.x

flashiness <- flow_ts %>%
  as_tibble() %>%
  mutate(Flow = Flow*0.0283168) %>%
  dplyr::group_by(site_no) %>%
  arrange(dateTime, .by_group = TRUE) %>%
  mutate(abs_delta_daily_q = abs(Flow - lag(Flow))) %>%
  drop_na() %>%
  #dplyr::group_by(waterYear, site_no) %>%
  summarise(path_length = sum(abs_delta_daily_q, na.rm = TRUE),
            total_q = sum(Flow, na.rm = TRUE)) %>%
  mutate(rb_flashiness = path_length/total_q) %>%
  #dplyr::group_by(site_no) %>%
  #summarise(mean_rb_flashiness = mean(rb_flashiness)) %>%
  inner_join(., lc_sites_metadata_all %>%
               dplyr::select(tributary, site_no),
             by = "site_no") %>%
  rename(name = tributary) %>%
  dplyr::select(name, 
                rb_flashiness
                )
```
