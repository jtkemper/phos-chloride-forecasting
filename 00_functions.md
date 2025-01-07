00_functions
================
JTK
2025-01-03

## Land Cover Processor

This function takes areas for each land cover class in the Troy et al.,
2007 land cover layer for the Lake Champlain Basin and derives
illustrative classes. In particular, we combine the Urban and Urban Open
class to create a Developed class, combine the Agriculture and Brush
classes to create an Agriculture class, and combine water and wetlands
to create a Water+Wetlands class. We also normalize areas of each class
by watershed area

``` r
#### Here is the function

lulc_processer <- function(file_string,
                           drn_area_df = NULL) {
  
      lulc <- read_csv(file_string) %>%
      rename_all(~tolower(.)) %>%
      dplyr::select(tributary, class, area) %>%
      rename(lulc_class = class,
             area_m2 = area)
    
      lulc <- lulc %>%
      full_join(., drn_area_df %>%
                  dplyr::select(tributary, drnarea_km2),
                by = "tributary") %>%
      mutate(pct_watershed_area = area_m2*1E-6/drnarea_km2*100) %>%
      dplyr::select(tributary, lulc_class, pct_watershed_area) %>%
      pivot_wider(id_cols = tributary, 
                  values_from = pct_watershed_area, names_from = lulc_class) %>%
      mutate(Developed = Urban + `Urban-Open`,
             Agriculture2 = Agriculture + Brush,
             Water_Wetland = Water + Wetland) %>%
      dplyr::select(!c(Urban, `Urban-Open`, 
                       Agriculture, Brush,
                       Water, Wetland)) %>%
      rename(Agriculture = Agriculture2) %>%
      mutate(across(where(is.numeric), ~round(., 1))) %>%
      rename_all(~tolower(.)) %>%
      rename_if(is.numeric, ~paste0("lulc_", .))
      
  
  return(lulc)
  
  
  
}
```

## HR Flowline length calculator

This function calculates the lengths of flowlines in each watershed from
the NHD HR. It also calculates various other flowline stats, like
elevation, basin slope, and the geomorphic type (based on slope)

``` r
hr_flow_length_calculator <- function(tributary = NULL,
                                      begin_comid = NULL,
                                      hr_data = NULL,
                                      hr_geopack = NULL,
                                      main_or_trib = "Trib") {
  
  flowline_stats <- list()
  
   print(paste("Reading ", tributary))
   
              if(main_or_trib == "Trib") {
                
                    ids <- get_UT(hr_data, begin_comid)
                
              } else if(main_or_trib == "Main") {
                
                ids <- get_UM(hr_data, begin_comid)
                
              }
 
              
            
              basin <- subset_nhdplus(comids = ids,
                                      flowline_only = TRUE,
                                    #output_file = outfile,
                                    nhdplus_data = hr_geopack,
                                    return_data = TRUE)

              basin_flowlines <- basin$NHDFlowline
              
              
              total_flowline_length <- sum(basin_flowlines$LENGTHKM)
              
              maximum_stream_order <- max(basin_flowlines$StreamOrde)
              
              pct_by_order <- basin_flowlines %>%
                as_tibble() %>%
                dplyr::group_by(StreamOrde) %>%
                summarise(pct_order = sum(LENGTHKM)/total_flowline_length)
              
              flowlines_and_types <- basin_flowlines %>%
                as_tibble() %>%
                mutate(channel_type = case_when(Slope > 0.26 ~
                                                    "colluvial",
                                                  Slope > 0.08 & Slope <= 0.26 ~
                                                    "cascade",
                                                  Slope > 0.03 & Slope <= 0.08 ~
                                                    "step_pool",
                                                  Slope > 0.01 & Slope <= 0.03 ~
                                                    "plane_bed",
                                                  Slope > 0.001 & Slope <= 0.01 ~
                                                    "pool_riffle",
                                                  Slope <= 0.001 ~
                                                    "dune_ripple")) 
              
              
              pct_by_type <- flowlines_and_types %>%
                dplyr::group_by(channel_type) %>%
                summarise(pct_by_type = sum(LENGTHKM)/total_flowline_length)
              

              median_reach_slope <- median(basin_flowlines$Slope)
              
              max_elev <- max(basin_flowlines$MaxElevSmo)
              
              min_elev <- min(basin_flowlines$MinElevSmo)
              
              total_basin_slope <- ((max(basin_flowlines$MaxElevSmo) - 
                min(basin_flowlines$MinElevSmo))/100)/(total_flowline_length*1000)
              
              
              ### Store them as list elements
              
              flowline_stats[[1]] <- total_flowline_length
              
              flowline_stats[[2]] <- maximum_stream_order
              
              flowline_stats[[3]] <- pct_by_order
              
              flowline_stats[[4]] <- median_reach_slope
              
              flowline_stats[[5]] <- max_elev
              
              flowline_stats[[6]] <- min_elev
              
              flowline_stats[[7]] <- total_basin_slope
              
              flowline_stats[[8]] <- basin_flowlines
              
              flowline_stats[[9]] <- pct_by_type
              
              #### If trib, return everything
              if(main_or_trib == "Trib") {
                
                return(flowline_stats)
                
              } else if(main_or_trib == "Main") {
                  
                  #### But if Main, just return the flowline length
                  #### Which is just the legnth of the mainstem
                  return(flowline_stats[[1]])
                  
                
                }
              
            
            
            
            

  
  
}
```

## Get all COMIDs

``` r
all_comid_getter <- function(tributary = NULL,
                                   begin_comid = NULL,
                                   trib_or_main = "Trib") {
  
  #### Track across tributaries
  
  print(paste("Downloading ", tributary))
      
  #### Get all the comids upstream from the outlet 
  
  comids <- navigate_nldi(list(featureSource = "comid",
                       featureID = begin_comid),
                  mode = "upstreamTributaries",
                  distance_km = 9999)
        
  #### Get just the comids and nest them
  #### And then append with the tributary they are for
  
  basin_comids <- comids$UT_flowlines %>%
    as_tibble() %>%
    dplyr::select(nhdplus_comid) %>%
    nest() %>%
    mutate(tributary = tributary)

      
  return(basin_comids)


  } 
```

## Pivot NHD Chars

This is a short function that basically just pivots the dataframe
retuned by nhdTools::get_catchment_characteristics into wide format

``` r
pivot_nhd_chars_wide <- function(nhd_char_df = NULL,
                                 comids_df = NULL){
  
  pivoted_df <- nhd_char_df %>%
    as_tibble() %>%
    dplyr::select(!contains("nodata")) %>%
    mutate(characteristic_id = tolower(characteristic_id)) %>%
    pivot_wider(names_from = characteristic_id,
                values_from = characteristic_value) %>%
    inner_join(comids_df, .,
             by = "comid") %>%
      dplyr::select(!comid)
    

  
  return(pivoted_df)
  
  
}
```

## Read chars directly from ScienceBase

This is a function that reads NHD flowline characteristic data
(Wieczorek et al., 2018)  
directly from the ScienceBase API

``` r
read_direct_from_sb <- function(sci_base_id = NULL) {
  
  ##### Print error if we are missing ScienceBase ID
  
  if(is.null(sci_base_id)){

    stop("Please enter a ScienceBase ID\n")
    
    }
  
    ##### To do this, first retrieve ScienceBase item
    ##### For this we need the ScienceBase ID, which is just the alphanumeric string
    ##### in the URL for the item of interest
    
    sb_item <- sbtools::item_get(sb_id = sci_base_id)
    
    ##### Download the file(s) from ScienceBase 
    
    sb_file_path <- sbtools::item_file_download(sb_item,
                                                dest_dir = tempdir(),
                                                overwrite_file = TRUE)
    
    ##### Find only the file with the data in it (this will be a zip file)
    
    sb_data_file <- str_subset(sb_file_path, ".zip")
    
    ##### Now we want to unzip this file
    
    sb_data_file_unz <- unzip(sb_data_file,
          exdir = dirname(sb_data_file),
          overwrite = TRUE)
    
    ##### Then, check if there are any zipped files in the newly unzipped file
    ##### And then unzip those
    ##### We want to do this until no zipped files remain
    
    sb_subfiles <- str_subset(sb_data_file_unz, ".zip")
    
    while(length(sb_subfiles) > 0) {
      
      sb_data_file_unz <- unzip(sb_subfiles,
                                     exdir = dirname(sb_data_file),
                                     overwrite = TRUE)
      
      sb_subfiles <- str_subset(sb_data_file_unz, ".zip")
      
      if(length(sb_subfiles) > 0){cat(crayon::yellow("Repeating"))}
      
    }
    
    ##### Finally, read in the unzipped data file(s)
    ##### And retrieve the distance to stream measure
    
    sb_actual_data <- read_csv(sb_data_file_unz)
    ##### And return final data
    
    return(sb_actual_data)

  
  
  
}
```

## Remove unsplittable features

This function essentially removes features that, in at least one
cross-validation split, have only one value for all tributaries that are
“left-in” in the training data. Briefly, this is an issue mostly because
it biases model training to variables that may serve to “identify” a
specific basin rather than reflect process

``` r
remove_unsplittable <- function(drivers_df) {
      
    ### Make an empty list to save things
  
    small_feats <- list()
    
    ### Loop over each basin 
    
      for(i in 1:18) {
        
              #### Track across basins
        
              removed_trib <- drivers_df %>%
                filter(group_id == i) %>%
                .$tributary %>%
                .[1]
              
              
              print(removed_trib)
              
              #### Determine how many unique values there are for each feature
              
              small_feats[[i]] <- drivers_df %>%
                filter(group_id != i) %>%
                dplyr::select(!c(constituent, 
                            site_no,
                            drnarea_km2)) %>%
                dplyr::select(!c(log_conc,
                                 date, 
                                 water_year,
                                 tributary,
                                group_id)) %>%
                summarise(across(everything(), ~length(unique(.x)))) %>% ### How many unique
                pivot_longer(everything(), values_to = "unique_feature_values",
                             names_to = "feature") %>%
                mutate(feature_value_rank = dense_rank(unique_feature_values)) %>%
                filter(feature_value_rank == 1) %>%
                dplyr::select(!feature_value_rank) %>%
                mutate(removed_trib = removed_trib)
                
      
      } ## End for loop
      
      all_small_feats <- bind_rows(small_feats)
      
      #### Get features that only have one unique value across at least one
      #### Cross-validation split 
      
      feats_just_one <- all_small_feats %>%
        filter(unique_feature_values < 2) %>%
        .$feature
        
        
      return(feats_just_one)
  
  
  
}
```
