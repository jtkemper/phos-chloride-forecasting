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
