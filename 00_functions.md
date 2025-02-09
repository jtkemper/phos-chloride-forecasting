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

## LightGBM Variable Selector

This function essentially performs a backwards variable selection on all
of the possible features for predicting total phosphorus or total
chloride concentration. It uses the mean SHAP value and standard
deviation of the SHAP values as criteria for ranking variables, and then
trims off the least important one. In this way, we can see which
variables result in the best performance.

``` r
## Run a loop of different years and different models to pick our best model

lgbm_selector <- function(constituent_df,
                          selection = TRUE,
                          model = NULL,
                          is_tuned = FALSE,
                           tuned_params = NULL,
                          selector = "shap",
                          model_type = "gaged"
                          
                          ) {
  
    
  

      model_stats_lgbm <- list()
      
      var_imp <- list()
      
      vars <- list()
      
      all_model_stats <- list()
      
      shap_values <- list()
      
      shap_value_summary <- list()
      
      model_returns <- list()
      
      predicted_observed_ts <- list()
      
      Ds <- list()
      
      
                ### These are the variables in our training data
                ### Need to refresh this to start
      
      if(selection == TRUE) {
        
                  predictors <- constituent_df %>%
              dplyr::select(!c(
                constituent, 
                site_no,
                drnarea_km2)) %>%
            mutate(water_year2 = water_year) ### Allows water year to be a predictor
                  
          
                  
  
      } else if (selection == FALSE) {
        
        
          predictors <- constituent_df %>%
            dplyr::ungroup() %>%
            dplyr::select(c(model$Feature, water_year,
                            log_conc, date, tributary, group_id)) %>%
            mutate(water_year2 = water_year)
        
        
        
      }

          
      
      
      for(i in 1:ncol(predictors)) {
        
        
        
          run <- paste0("run", i)
        
                    
          cat(crayon::yellow("\nModeling Run", i, "\n"))

      
          if(model_type == "gaged") {
            
            

                
                  for(j in c(2003, 2008, 2013, 2018)) {
                
              
                        test_years <- seq(j-4, j, 1)
                        
                        train_years <- seq(1990, j-5, 1)
                        
                        k <- (2018 - j)/5
                  
                        k <- ifelse(k == 0, 1, k+1)
                        
                        
                        
                        cat(crayon::cyan("\nTraining ", first(train_years), ":",
                                           last(train_years),
                                           "\n")) 
                        cat(crayon::magenta("\nTesting ", first(test_years), ":",
                                           last(test_years),
                                           "\n")) 
                        
                            
                        
                        #### Training data
                            predictor_dataset_train <- predictors %>%
                            filter(water_year %in% train_years)
                            
                            if(nrow(predictor_dataset_train) < 1) next
                     
                          ### Now subset the testing data to that same subset
                          predictor_dataset_test <- predictors %>%
                            filter(water_year %in% test_years) %>%
                            dplyr::select(colnames(predictor_dataset_train))
                          
                          #######################
                        
                          ### Declare the predictor and response variables 
                        preds <- data.matrix(predictor_dataset_train %>%
                                                    dplyr::select(!c(log_conc,
                                                                     water_year,
                                                                     date, 
                                                                     tributary,
                                                                     group_id
                                                                     )))
                        
                        
                        response <- predictor_dataset_train$log_conc
                        
                        ### Set up the environment - this is just preparing the dataset API to be used by lightgbm. This is our training data
                        train_lgbm <- lgb.Dataset(preds, 
                                                         label = response,
                                                         #feature_pre_filter = FALSE,
                                                         ) 
                        
                        ### Declare the test data
                        test_lgbm <- data.matrix(predictor_dataset_test %>%
                                                              dplyr::select(!c(log_conc,
                                                                     water_year,
                                                                     date, 
                                                                     tributary,
                                                                     group_id
                                                                     )))
                    
                        ### Declare the hyperparameters
                    
                        if(is_tuned == FALSE ) { 
                        hyperparams <- list(objective = "regression",
                                                num_leaves = 31L,
                                                learning_rate = 0.1,
                                                min_data_in_leaf = 20L,
                                                num_threads = 10L
                                            )
                        
                        
                        } else if(is_tuned == TRUE ) {
                          
                          
                              hyperparams <- list(objective = "regression",
                                    num_leaves = tuned_params$num_leaves,
                                    min_data_in_leaf = tuned_params$min_n,
                                    bagging_fraction = tuned_params$sample_size,
                                    bagging_freq = 1,
                                    num_iterations = tuned_params$trees,
                                    max_depth = tuned_params$tree_depth
                                    )
                          
                          
                          
                          
                        }
                        
                        
                        #########
                        
                        ### Train the model
                        
                        set.seed(913)
                        
                        nutrient_model_lgbm <- lgb.train(hyperparams,
                                                          data = train_lgbm,
                                                          verbose = 1L,
                                                          nrounds = 100L
                                                         )
                        
                        ### Get model fits on training data
                        nutrient_fits <- predict(nutrient_model_lgbm, 
                                                       data = preds) %>%
                          as_tibble() %>% rename(log_predicted_conc = 1)
                        
                        
                        ### Predict with the model on test data
                        nutrient_predicted <- predict(nutrient_model_lgbm, 
                                                       data = test_lgbm) %>%
                          as_tibble() %>% rename(log_predicted_conc = 1)
                        
                        ### Calculate the SHAP values
                        shap_values[[k]] <- SHAPforxgboost::shap.prep(xgb_model = 
                                                                        nutrient_model_lgbm, 
                                                                   X_train = test_lgbm)
                        
                        shap_value_summary[[j]] <- shap_values[[k]] %>%
                              as_tibble() %>%
                              dplyr::group_by(variable) %>%
                              summarise(sd_shap = sd(value),
                                        feature_importance = mean_value[1]) %>%
                              mutate(sd_plus_imp = sd_shap + feature_importance)
                        
                        #### Bind to training observations 
                        #### and estimate smearing coefficient
                        fitted_observed_smear <- bind_cols(predictor_dataset_train %>%
                                                                dplyr::rename(log_observed_conc =
                                                                                log_conc),
                                                                       nutrient_fits) %>%
                          mutate(exp_log_model_residuals = 10^(log_observed_conc -
                                                                 log_predicted_conc))
                        
                        ### Estimate smearing coefficient  
                        D <- mean(fitted_observed_smear$exp_log_model_residuals)
                        
                        ### Bind predictions on test data
                        ### to observatios of test data
                        predicted_observed <- bind_cols(predictor_dataset_test %>%
                                                                dplyr::rename(log_observed_conc =
                                                                                log_conc),
                                                                       nutrient_predicted) 
                        
                        ### Now use the smearing coefficient to convert back to non-log
                        predicted_observed_plus_err <- predicted_observed %>%
                          mutate(predicted_conc = (10^log_predicted_conc)*D) %>%
                          mutate(observed_conc = 10^log_observed_conc) %>%
                          mutate(raw_err = predicted_conc - observed_conc) %>%
                           mutate(sqrerr = raw_err^2,
                                  abs_error = abs(raw_err),
                                  abs_pct_error = abs_error/predicted_conc
                                 )
                        
                        predicted_observed_ts[[j]] <- predicted_observed_plus_err %>%
                          dplyr::select(tributary, date, predicted_conc, observed_conc)
                        
                        Ds[[j]] <- tibble(d_fac = D,
                                          tributary = predictor_dataset_test$tributary[1])
                        
                    
                        #### Evaluate
                        model_stats_lgbm[[j]] <- predicted_observed_plus_err %>%
                          ungroup() %>%
                          dplyr::group_by(tributary) %>%
                          summarise(rmse_inst_conc = sqrt(mean(sqrerr)),
                                    mae_inst_conc = mean(abs_error),
                                    nse = hydroGOF::NSE(predicted_conc, 
                                                        observed_conc),
                                    kge = hydroGOF::KGE(predicted_conc, 
                                                        observed_conc),
                                    mape = mean(abs_pct_error),
                                    pbias = hydroGOF::pbias(predicted_conc,
                                                            observed_conc),
                                    r = cor(predicted_conc, observed_conc,
                                            method = "pearson", use = "pairwise.complete.obs"),
                                    br2 = hydroGOF::br2(predicted_conc, observed_conc),
                                    bias = mean(predicted_conc)/mean(observed_conc),
                                    variability = sd(predicted_conc)/mean(observed_conc)) %>%
                          dplyr::ungroup() %>%
                          summarise(median_rmse = median(rmse_inst_conc),
                                    median_mae = median(mae_inst_conc),
                                    median_nse = median(nse),
                                    median_kge = median(kge),
                                    median_r = median(r),
                                    median_pbias = median(pbias),
                                    median_br2 = median(br2),
                                    sd_rmse = sd(rmse_inst_conc),
                                    sd_mae = sd(mae_inst_conc),
                                    sd_nse = sd(nse),
                                    sd_kge = sd(kge),
                                    sd_r = sd(r),
                                    sd_pbias = sd(pbias),
                                    sd_br2 = sd(br2)) %>%
                          mutate(test_years = paste(j-1, j, sep = "_"))

                        
                        var_imp[[j]] <- lgb.importance(nutrient_model_lgbm , 
                                                                         percentage = TRUE)
                        
                        
                  }
            #####################################################################
          
          } else if (model_type == "ungaged") {
            
            for(j in 1:18) {
              
                                tr <- constituent_df %>% 
                                  filter(group_id == j) %>%
                                  dplyr::slice(1) %>%
                                  .$tributary
                                
                              cat(crayon::cyan("\nTesting Trib", tr, "\n")) 
                              
                          #### Training data
                          predictor_dataset_train <- predictors %>%
                                                    filter(group_id != j) 
                    
                          if(nrow(predictor_dataset_train) < 1) next
             
                       ### Now subset the testing data to that same subset
                        predictor_dataset_test <- predictors %>%
                            filter(group_id == j) %>%
                            dplyr::select(colnames(predictor_dataset_train))
                          
                                   #######################
                
                  ### Declare the predictor and response variables 
                preds <- data.matrix(predictor_dataset_train %>%
                                            dplyr::select(!c(log_conc,
                                                             water_year,
                                                             date, 
                                                             tributary,
                                                             group_id
                                                             )))
                
                
                response <- predictor_dataset_train$log_conc
                
                ### Set up the environment - 
                ### this is just preparing the dataset API to be used by lightgbm. This is our training data
                train_lgbm <- lgb.Dataset(preds, 
                                                 label = response
                                                 ) 
                
                ### Declare the test data
                test_lgbm <- data.matrix(predictor_dataset_test %>%
                                                      dplyr::select(!c(log_conc,
                                                             water_year,
                                                             date, 
                                                             tributary,
                                                             group_id
                                                             )))
            
                ### Declare the hyperparameters
            
                if(is_tuned == FALSE ) { 
                hyperparams <- list(objective = "regression",
                                        num_leaves = 31L,
                                        learning_rate = 0.1,
                                        min_data_in_leaf = 20L,
                                        num_threads = 10L
                                    )
                
                
                } else if(is_tuned == TRUE ) {
                  
                  
                      hyperparams <- list(objective = "regression",
                            num_leaves = tuned_params$num_leaves,
                            min_data_in_leaf = tuned_params$min_n,
                            bagging_fraction = tuned_params$sample_size,
                            bagging_freq = 1,
                            num_iterations = tuned_params$trees,
                            max_depth = tuned_params$tree_depth
                            )
                  
                  
                  
                  
                }
                
                
                #########
                #########
                
                ### Train the model
                
                set.seed(913)
                
                nutrient_model_lgbm <- lgb.train(hyperparams,
                                                  data = train_lgbm,
                                                  verbose = 1L,
                                                  nrounds = 100L
                                                 )
                
                ### Predict 
                
                ### Get model fits on training data
                nutrient_fits <- predict(nutrient_model_lgbm, 
                                               data = preds) %>%
                  as_tibble() %>% rename(log_predicted_conc = 1)
                
                
                ### Predict with the model on test data
                nutrient_predicted <- predict(nutrient_model_lgbm, 
                                               data = test_lgbm) %>%
                  as_tibble() %>% rename(log_predicted_conc = 1)
                
                #######
                
                ### Calculate the SHAP values
                shap_values[[j]] <- SHAPforxgboost::shap.prep(xgb_model = 
                                                                nutrient_model_lgbm, 
                                                           X_train = test_lgbm)
                
                shap_value_summary[[j]] <- shap_values[[j]] %>%
                      as_tibble() %>%
                      dplyr::group_by(variable) %>%
                      summarise(sd_shap = sd(value),
                                feature_importance = mean_value[1]) %>%
                      mutate(sd_plus_imp = sd_shap + feature_importance)
                
                #########
                
                #### Calcuate smearing factor
                
                #### Bind to training observations 
                #### and estimate smearing coefficient
                fitted_observed_smear <- bind_cols(predictor_dataset_train %>%
                                                        dplyr::rename(log_observed_conc =
                                                                        log_conc),
                                                               nutrient_fits) %>%
                  mutate(exp_log_model_residuals = 10^(log_observed_conc -
                                                         log_predicted_conc))
                
                ### Estimate smearing coefficient  
                D <- mean(fitted_observed_smear$exp_log_model_residuals)
                
                #####
                
                ### Transform predictions on test data
                
                ### Bind predictions on test data
                ### to observations of test data
                predicted_observed <- bind_cols(predictor_dataset_test %>%
                                                        dplyr::rename(log_observed_conc =
                                                                        log_conc),
                                                               nutrient_predicted) 
                
                ### Now use the smearing coefficient to convert back to non-log
                predicted_observed_plus_err <- predicted_observed %>%
                  mutate(predicted_conc = (10^log_predicted_conc)*D) %>%
                  mutate(observed_conc = 10^log_observed_conc) %>%
                  mutate(raw_err = predicted_conc - observed_conc) %>%
                   mutate(sqrerr = raw_err^2,
                          abs_error = abs(raw_err),
                          abs_pct_error = abs_error/predicted_conc,
                         )
                
                ### Save stuff
                
                predicted_observed_ts[[j]] <- predicted_observed_plus_err %>%
                  dplyr::select(tributary, date, predicted_conc, observed_conc)
                
                Ds[[j]] <- tibble(d_fac = D)
                
            
                #### Evaluate
                model_stats_lgbm[[j]] <- predicted_observed_plus_err %>%
                  ungroup() %>%
                  dplyr::group_by(tributary) %>%
                  summarise(rmse_inst_conc = sqrt(mean(sqrerr)),
                            mae_inst_conc = mean(abs_error),
                            nse = hydroGOF::NSE(predicted_conc, 
                                                observed_conc),
                            kge = hydroGOF::KGE(predicted_conc, 
                                                observed_conc),
                            mape = mean(abs_pct_error),
                            pbias = hydroGOF::pbias(predicted_conc,
                                                    observed_conc),
                            r = cor(predicted_conc, observed_conc),
                            br2 = hydroGOF::br2(predicted_conc, observed_conc)
                            ) %>%
                          dplyr::ungroup() %>%
                          summarise(median_rmse = median(rmse_inst_conc),
                                    median_mae = median(mae_inst_conc),
                                    median_nse = median(nse),
                                    median_kge = median(kge),
                                    median_r = median(r),
                                    median_pbias = median(pbias),
                                    median_br2 = median(br2),
                                    sd_rmse = sd(rmse_inst_conc),
                                    sd_mae = sd(mae_inst_conc),
                                    sd_nse = sd(nse),
                                    sd_kge = sd(kge),
                                    sd_r = sd(r),
                                    sd_pbias = sd(pbias),
                                    sd_br2 = sd(br2)) %>%
                          mutate(test_years = paste(j-1, j, sep = "_"))
                
                
                var_imp[[j]] <- lgb.importance(nutrient_model_lgbm , 
                                                                 percentage = TRUE)
                
             
              
            }
            

            
            
            
            
          }
          
          if(selection == TRUE) {
            
            
              all_var_imp <- bind_rows(var_imp) 
              
              all_shap_value_summary <- bind_rows(shap_value_summary) %>%
                dplyr::group_by(variable) %>%
                summarise(mean_sd_plus_imp = mean(sd_plus_imp))
          
          all_model_stats[[i]] <- bind_rows(model_stats_lgbm) %>%
            mutate(model = i)
          
          summary_var_imp <- all_var_imp %>%
            dplyr::group_by(Feature) %>%
            summarise(mean_Gain = mean(Gain)) %>%
            dplyr::ungroup()
          
            if(selector == "shap") {
              
                one_removed_predictors <- all_shap_value_summary %>%
                  dplyr::ungroup() %>%
                  arrange(desc(mean_sd_plus_imp)) %>%
                  dplyr::slice(-nrow(.))
                
          vars[[i]] <- all_shap_value_summary %>%
            ungroup() %>%
            mutate(model = i) %>%
            rename(Feature = variable)
          
          ### See how many we have left
          var_count <- length(one_removed_predictors$variable)
          
          if(var_count == 0) break 
          
          
            ### Update variable list
          predictors <- predictors %>%
            dplyr::select(one_removed_predictors$variable,
                          log_conc,
                          water_year,
                          date, 
                          tributary,
                          group_id)
          
              
            } else if(selector == "gain"){
              
                one_removed_predictors <- summary_var_imp %>%
                  dplyr::ungroup() %>%
                  arrange(desc(mean_Gain)) %>%
                  dplyr::slice(-nrow(.))
                
                          vars[[i]] <- summary_var_imp %>%
            ungroup() %>%
            mutate(model = i)
          
          ### See how many we have left
          var_count <- length(one_removed_predictors$Feature)
          
          if(var_count == 0) break 
          
          
            ### Update variable list
          predictors <- predictors %>%
            dplyr::select(one_removed_predictors$Feature,
                          log_conc,
                          water_year,
                          date, 
                          tributary,
                          group_id)

            }

          
          

          

            
            
            
          } else if(selection == FALSE) {
            
              all_var_imp <- bind_rows(var_imp) 
              
              all_predicted_observed_ts <- bind_rows(predicted_observed_ts)
              
              all_Ds <- bind_rows(Ds)
              
              model_returns[[6]] <- all_Ds
              
              model_returns[[5]] <- all_predicted_observed_ts
              
              
              all_model_stats[[i]] <- bind_rows(model_stats_lgbm) %>%
                mutate(model = model$model[1])
              
              summary_var_imp <- all_var_imp %>%
                    group_by(Feature) %>%
                    summarise(mean_Gain = mean(Gain)) %>%
                    ungroup()
                
              #all_shap_values <- bind_rows(shap_values)
              
              model_returns[[4]] <- shap_values
              
              vars[[i]] <- summary_var_imp %>%
                dplyr::ungroup() %>%
                mutate(model = model$model[1])
            
              break
            
          } 
          
        
      
      }
      
      all_all_model_stats <- bind_rows(all_model_stats)
      
      all_var_imp <- bind_rows(vars)
      
      
      #### Save outputs
      
 
      
      model_returns[[1]] <- all_var_imp
      
      #model_returns[[4]] <- shap_values
      
      
      ## examine the model summary statistics
      summary_model_stats <- all_all_model_stats %>%
        ungroup()
      
      
      collapsed_models <- all_var_imp %>%
        group_by(model) %>%
        arrange(Feature, .by_group = TRUE) %>%
        summarise(all_vars = paste(Feature, collapse = ",")) %>%
        full_join(., summary_model_stats, 
                  by = "model")
      
        
        model_stats <- collapsed_models %>%
              dplyr::group_by(model, all_vars) %>%
                  summarise(#n = n(),
                            mean_kge = mean(median_kge, na.rm = TRUE),
                                  mean_nse = mean(median_nse, na.rm = TRUE),
                                  mean_mae = mean(median_mae, na.rm = TRUE),
                                  mean_rmse = mean(median_rmse, na.rm = TRUE),
                                  mean_pbias = mean(median_pbias, na.rm = TRUE),
                                  mean_r = mean(median_r, na.rm = TRUE),
                                  mean_br2 = mean(median_br2, na.rm = TRUE),
                                  mean_sd_kge = sd(median_kge),
                                  mean_sd_nse = sd(median_nse),
                                  mean_sd_mae = sd(median_mae),
                                  mean_sd_rmse = sd(median_rmse),
                                  mean_sd_pbias = sd(median_pbias),
                                  mean_sd_r = sd(median_r),
                                  mean_sd_br2 = sd(median_br2),
                            sd_kge = sd(median_kge),
                            sd_nse = sd(median_nse),
                            sd_mae = sd(median_mae),
                            sd_rmse = sd(median_rmse),
                            sd_pbias = sd(median_pbias),
                            sd_r = sd(median_r),
                            sd_br2 = sd(median_pbias),
                            ) 
            
              model_returns[[2]] <- model_stats
              
        
        
      
      
      return(model_returns)

}
```

## Plot variable selection stats

This function takes a dataframe of model performance stats generated by
the variable selection process and plots those as a function of model
“number”, where model number 1 has all n variables and model n has 1
variable

``` r
plot_stats <- function(model_stats_df) {
  
  model_stats_df %>%
    dplyr::select(model, 
                  mean_kge, mean_nse, mean_mae, mean_br2, mean_pbias,
                  sd_kge,  sd_nse,  sd_mae, sd_br2, sd_pbias) %>%
    pivot_longer(cols = -model, names_to = c(".value", "metric"), names_sep = "_") %>%
    ggplot() +
      geom_line(aes(x = model, y = mean, color = metric)
                #color = "black"
                ) +
      # geom_errorbar(aes(x = model, ymin = mean-sd, ymax = mean + sd),
      #               color = "gray65") + 
      geom_point(aes(x = model, y = mean, color = metric),
                #color = "black", 
                shape = 19, size = 1) +
      scale_color_brewer(palette = "Set1",
                         guide = "none") + 
      #scale_y_log10() + 
      scale_x_continuous(breaks = seq(0,95,5),
                         minor_breaks = seq(0,95,1)) + 
      labs(y = element_blank(),
           x = "Model Iteration") + 
      theme_few() +
      theme(panel.grid = element_line(color = "gray90"),
            legend.position = "bottom",
            strip.placement = "outside") +
      facet_wrap(~metric, scales = "free", ncol = 1,
                 strip.position = "left")

  
}
```

## Make training-validation splits for tuning

This function makes training-validation splits for both gaged and
ungaged scenarios in the formats necessary to perform hyperparameter
tuning.

``` r
train_valid_splitter <- function(trimmed_df = NULL,
                                 model_features_df = NULL,
                                 loocv = FALSE) {
  
    
    ### Empty list where we are going to save things
  
    split_list <- list()
    
    #### Make splits for a gaged scenario 
    #### We will have four training-validation splits
    #### (1990-99 training, 2000-03 valid; etc.)
    
    if(loocv == FALSE){
      
             for(h in c(2003, 2008, 2013, 2018)) {
              
                  test_years <- seq(h-4, h, 1)
                  
                  train_years <- seq(1990, h-5, 1)
              
                  #### Training data
                  
                  trimmed_train <- trimmed_df %>%
                      filter(water_year %in% train_years) %>%
                      dplyr::select(model_features_df$Feature, log_conc)
                      
                      if(nrow(trimmed_train) < 1) next
              
                    ### Now subset the testing data to that same subset
                      
                    trimmed_test <- trimmed_df %>%
                      filter(water_year %in% test_years) %>%
                      dplyr::select(model_features_df$Feature, log_conc)
                    
                     k <- (2018 - h)/5
                    
                    k <- ifelse(k == 0, 1, k+1)
                          
                    
                    split_list[[k]] <- (rsample::make_splits(trimmed_train, 
                                                   assessment = trimmed_test)) 
                    
                    } ### End training/valid split for loop
    
            return(split_list)
      
    #### Make splits for an ungaged scenario 
    #### We will have 18 training-validation splits
    #### 17 watershed training, 1 validation 

    } else if (loocv == TRUE) {
      
        for(g in 1:18){
          
          
               #### Training data
                      trimmed_train <- trimmed_df %>%
                      filter(group_id != g) %>%
                        dplyr::select(!group_id) %>%
                      dplyr::select(model_features_df$Feature, log_conc)
                      
                      if(nrow(trimmed_train) < 1) next
              
                    ### Now subset the testing data to that same subset
                    trimmed_test <- trimmed_df %>%
                      filter(group_id == g) %>%
                      dplyr::select(!group_id) %>%
                      dplyr::select(model_features_df$Feature, log_conc)
                    
                    split_list[[g]] <- (rsample::make_splits(trimmed_train, 
                                                   assessment = trimmed_test)) 
          
          
          
          
          
          } ### End of loocv split list
          
             return(split_list)
    
      }
    
 
}
```

## Tune LightGBM models

This function tunes the hyperparameters for each LightGBM model and
modeling scenario.

``` r
#### FUNCTION THAT TUNES A CHOSEN LGBM MODEL ##################################

lgbm_tuner <- function(observational_df, model_df, 
                       small = FALSE,
                       loocv = FALSE) {

  
      ### First, trim down the dataset to just our identified predictors
  
      trimmed_data <- observational_df %>%
        filter(water_year < 2019) %>%
        ungroup() 
      
      
      ### Now, pull out the training data
      
      training_data <- trimmed_data %>%
        filter(water_year < 2019) %>%
        dplyr::select(model_df$Feature, group_id, log_conc) %>%
        dplyr::select(!group_id)
      
      
      ### Now, build our training & validation split
      ### This is an expanding window design
      ### Here, we are telling tidymodels which bits of the training data
      ### We want to use as resamples (validation data)
      ### To tune the hyperparameters
      ### This is akin to 4-fold CV but for time series data
      ### We are basically training the model on 1990-1998 and then validating on 1999-2003
      ### Training on 1990-2003 and validating on 2004-2008
      ### And so on and so forth
      
      ### Make splits
      
      if(loocv == FALSE){
        
        splits_list <- train_valid_splitter(trimmed_data,
                                            model_df,
                                            loocv = FALSE)
        
        split_names <- paste0("Split", seq(1:length(splits_list)))
        
        
        
      } else if(loocv == TRUE) {
        
        splits_list <- train_valid_splitter(trimmed_data,
                                            model_df,
                                            loocv = TRUE)
        
        split_names <- paste0("Split", seq(1:length(splits_list)))
        
        
        
      }

      ### Declare the training/validation splits
      
      valid_data_expand <- rsample::manual_rset(splits_list, split_names)
      
      ### Create a model recipe
      
      #### First, find the variables we have selected
      
      rec <- recipe(log_conc ~ .,
                         data = training_data) %>%
        prep()
      
      
      #### Declare the specific parameters we want tuned
      
      lgbm_model <- parsnip::boost_tree(mode = "regression",
                                              min_n = tune(),
                                             tree_depth = tune(),
                                             trees = tune(),
                                             sample_size = tune(),
                                             ) %>%
        set_engine("lightgbm", metric = "mae", verbose = -1, 
                   num_leaves = tune(),
                   counts = FALSE,
                   num_threads = 10)

      #### Supply ranges to examine those parameters
      
      if(small == TRUE) {
        
                  lgbm_params <- dials::parameters(
                                  min_n(c(1L,20L), trans = NULL),
                                           tree_depth(c(3,12L)),
                                           num_leaves(c(2,32)),
                                           sample_prop(c(1/10, 1)),
                                           trees(c(10,1000))
                                           )
        
        
        
        
      } else if(small == FALSE) {
        
        
          lgbm_params <- dials::parameters(min_n(c(5L,100L), trans = NULL),
                            tree_depth(c(3,12L)),
                                           num_leaves(c(8,256)),
                                           sample_prop(c(1/10, 1)),
                                           trees(c(10,1000))
                               )
        
        
        
      }
      
          #### Just a small factor to control the number of models we run
          #### For loocv scenarios (this means we run 50 diff models in that case)
      
          divider <- ifelse(loocv == TRUE, 2, 1)
          
          
          #### Construct the grid over which we are going to do the search 
          
          lgbm_grid <- dials::grid_max_entropy(lgbm_params,
                                                size = ifelse(small == TRUE,
                                                              20, 100/divider))
          
      
          #### Set up the workflow
          
          lgbm_wf <- workflows::workflow() %>%
            add_model(lgbm_model) %>%
            add_formula(log_conc ~ .)
          
          ## DO it
          
          lgbm_tuned <- tune::tune_grid(object = lgbm_wf,
                                        resamples = valid_data_expand,
                                        grid = lgbm_grid,
                                        metrics = yardstick::metric_set(rmse, mae, rsq),
                                        control = tune::control_grid(verbose = TRUE))
          
          
          #### Extract scores  
          
         rmse_scores  <- lgbm_tuned %>%
            tune::show_best(metric = "rmse", n = 100) %>%
            mutate(mean_rmse_rank = dense_rank(mean)) %>%
            mutate(mean_rmse_se_rank = dense_rank(std_err)) %>%
            mutate(rmse_mean_plus_se_rank  = mean_rmse_rank + mean_rmse_se_rank) %>%
            dplyr::rename(mean_rmse = mean,
                   se_rmse = std_err)
         
          mae_scores  <- lgbm_tuned %>%
            tune::show_best(metric = "mae", n = 100) %>%
            mutate(mean_mae_rank = dense_rank(mean)) %>%
            mutate(mean_mae_se_rank = dense_rank(std_err)) %>%
            mutate(mae_mean_plus_se_rank  = 0.7*mean_mae_rank + 0.3*mean_mae_se_rank) %>%
            dplyr::rename(mean_mae = mean,
                   se_mae = std_err)
          
          rsq_scores  <- lgbm_tuned %>%
            tune::show_best(metric = "rsq", n = 100) %>%
            mutate(mean_rsq_rank = dense_rank(mean)) %>%
            mutate(mean_rsq_se_rank = dense_rank(std_err)) %>%
            mutate(rsq_mean_plus_se_rank  = 0.7*mean_rsq_rank + 0.3*mean_rsq_se_rank) %>%
            dplyr::rename(mean_rsq = mean,
                   se_rsq = std_err)
          
          scores <- inner_join(rmse_scores,
                               mae_scores %>%
                                 dplyr::select(.config, mean_mae,
                                               se_mae,
                                               mean_mae_rank,
                                               mean_mae_se_rank,
                                               mae_mean_plus_se_rank),
                               by = ".config") %>%
            inner_join(.,
                       rsq_scores %>%
                                 dplyr::select(.config, 
                                               mean_rsq,
                                               se_rsq,
                                               mean_rsq_rank,
                                               mean_rsq_se_rank,
                                               rsq_mean_plus_se_rank),
                               by = ".config") 
          
      return(scores)
             
} ### End function
```

## Choose hyper parameters

This function essentially runs the top 10 performing hyperparameters
returned by the tuning processes to make a final selection of the hypers

``` r
####  PLOTS TO EVALUATE WHICH OF THE HYPERPARAMETERS IN THE TOP TEN FIT BEST ###

#### Re-run tuned to pick best ####

hyper_chooser <- function(constit_df = NULL,
                          hypers_df = NULL,
                          selected_params = NULL, 
                          loocv = FALSE,
                          constit = NULL) {
  
  
          #### Empty lists to save things
  
          stats_by_trib_tunez <- list()
            

          #hyper_choosey <- list()
          
          #### Loop over the top ten
          
          for(i in 1:10) {
            
          #### Track progress 
            
          print(paste("Running hyperparams #", i))
          
          #### Pick a individual set of parameters
          
          best_params <- hypers_df %>%
            arrange(mae_mean_plus_se_rank) %>%
            dplyr::slice(i)
          
          #### If statements for which scenario we are tuning for
          
          if(loocv == FALSE){
            
              #### Run the LightGBM model  
            
              final_models <- lgbm_selector(constituent_df = constit_df,
                                         selection = FALSE,
                                         model = selected_params,
                                         is_tuned = TRUE,
                                        tuned_params = best_params,
                                        model_type = "gaged"
                                       
                                         )
              
              #### Evaluate the model performance
              
              stats_by_trib_tunez[[i]] <- final_models[[2]] %>%
                            mutate(hyper_param = i) %>%
                            rename(kge = mean_kge,
                                   nse = mean_nse,
                                   mae = mean_mae,
                                   rmse = mean_rmse,
                                   pbias = mean_pbias,
                                   br2 = mean_br2)
          
              

            
          } else if(loocv == TRUE){
            
            #### Run the model again in a loocv scenario 
            
              final_models <- lgbm_selector(constituent_df = constit_df,
                                         selection = FALSE,
                                         model = selected_params,
                                         is_tuned = TRUE,
                                        tuned_params = best_params,
                                        model_type = "ungaged"
                                         )
              
              #### Evaluate
              
              stats_by_trib_tunez[[i]] <- final_models[[2]] %>%
                            mutate(hyper_param = i) %>%
                            rename(kge = mean_kge,
                                   nse = mean_nse,
                                   mae = mean_mae,
                                   rmse = mean_rmse,
                                   pbias = mean_pbias,
                                   br2 = mean_br2)
            
          } ### End if statement
 
          } ### End of for loop over each of the top ten
            
          #### Bind performance stats together
          
          stats_by_trib_tunez_all <- bind_rows(stats_by_trib_tunez)
          
          #### Output them
          
          hyper_choosey <- stats_by_trib_tunez_all %>%
            dplyr::select(hyper_param, 
                          kge, nse,
                          mae, rmse,
                          pbias,
                          br2)
          
           
           return(hyper_choosey)
            
}
```

## Test LightGBM model

This is a function that essentially makes it easier to run LightGBM
models in a testing scenario (so with features we’ve selected and
hyperparameters tuned). It trains models that contain user-idenified
features and evaluates how models performs on test data, which
heretofore has not been involved in any of the model construction steps.
This is, essentially, a robust test on how the model might do on data it
has never seen before. It also saves trained model “architecture” to
file, so that they may be used again in a forecasting capacity in a
later step.

``` r
########################### FUNCTION TO RUN LGBM ##############################
### Importantly, this is really most useful when you only want to run it once
### The runs we have set up with the leave one out are more complex and require 
### Their own custom thing to be written
### Okay here's the function 

#### It takes a bunch of arguments
#### 1) training_data: the training data, with concentration data and all 
#### features
#### 2) testing_data: the test dataset, with concentration data and all 
#### features
#### 3) model: a dataframe that has contains the names of features we have selected in 
#### the backwads variable selection process
#### 4) save: boolean about whether to save final model architecture to file
#### Importantly, this is essential to do if we want to use trained models for future applications 
#### (like forecasting)
#### 5) save_file: a string to help format save files in a way that will make them easy to find
#### 6) tuned: boolean whether to use tuned or default hyperparameters
#### 7) tuned_params: dataframe with values for the tuned hyperparameters
#### 8) small: boolean that guides which default hyperparameters to use
#### 9) solo: boolean that tells the function whether we are running a "global" model
#### (which contains both dynamic hydrology features and our selected watershed attributes
#### and has all basins in the training data)
#### or a basin-specific model which contains on dynamic hydrology features and has only 
#### one basin in training (and testing) data
#### loo: boolean that tells the function whether it as a gaged or ungaged scenario 

lgbm_runner <- function(training_data = NULL, 
                        testing_data = NULL, 
                        model = NULL,
                        save = FALSE,
                        save_file = NULL,
                        tuned = FALSE,
                        tuned_params = NULL,
                        small = FALSE,
                        solo = FALSE,
                        loo = FALSE) {
  
        stats_and_importance <- list()
        
        #save_file <- as.character(save_file)
      
      
          ################## SET UP MODEL ############################################
        

          ### These are the variables in our training data
        chosen_model_train <- training_data %>%
            ungroup() %>%
            dplyr::select(c(model$Feature, log_conc, date, tributary)) %>%
             dplyr::select(order(colnames(.)))
        
    
        
            ### Now subset the testing data to that same subset
           chosen_model_test <- testing_data %>%
              ungroup() %>%
              dplyr::select(colnames(chosen_model_train)) %>%
             dplyr::select(order(colnames(.)))
      
        
              ### Declare the predictor and response variables 
            preds <- data.matrix(chosen_model_train %>%
                                        dplyr::select(!c(log_conc,
                                                 date, 
                                                 tributary
                                                 )))
            
            
            response <- chosen_model_train$log_conc
            
            ### Set up the environment - this is just preparing the dataset API to be used by lightgbm. 
            ### This is our training data
            train_lgbm <- lgb.Dataset(preds, 
                                             label = response,
                                             #feature_pre_filter = FALSE,
                                             ) 
            
            ### Declare the test data
            test_lgbm <- data.matrix(chosen_model_test %>%
                                              dplyr::select(!c(log_conc,
                                                 #water_year,
                                                 date, 
                                                 tributary
                                                 )))
            

            ### Declare the hyperparameters
            ### These are defaults, but we have spelled them out 
            ### To be most clear
            
            ### User can specify whether they want tuned params or the defaults
            
            if(tuned == FALSE) {
              
                params <- list(objective = "regression",
                               metric = "mae",
                                num_leaves = ifelse(small == FALSE, 31L, 10),
                                learning_rate = 0.1,
                                min_data_in_leaf = ifelse(small == FALSE, 20L, 5),
                               num_iterations = 100
                               ) ### These are the defaults
            
            } else if(tuned == TRUE) {
              
                params <- list(objective = "regression",
                               metric = "mae",
                            num_leaves = tuned_params$num_leaves,
                            min_data_in_leaf = tuned_params$min_n,
                            bagging_fraction = tuned_params$sample_size,
                            bagging_freq = 1,
                            #feature_fraction = tuned_params$mtry,
                            num_iterations = tuned_params$trees,
                            max_depth = tuned_params$tree_depth
                            )
              
            }
            
            
    ############################################################################
    
    
    ############################# MAKE MODEL #####################################
            

            ### Train the model
            
            set.seed(913)
            

            
            model_lgbm <- lgb.train(params,
                                    data = train_lgbm,
                                    verbose = 1L, 
                                    nrounds = 100L) ###nrounds is the default 100
            
            ### Get values of splits and stuff
            tree_table <- lgb.model.dt.tree(model_lgbm) %>%
              as_tibble()
            
            
            #### Save the model 
            ### According to user specifications
            ### Default is not to save the model 
            
            if(save == TRUE) {
              
              if(solo == FALSE) {
                
                  if(loo == FALSE){
                    
                     lgb.save(model_lgbm, 
                              filename =
                                here("data", "models", "lumped", 
                                     paste0("final_tuned_lgbm_for_",
                                            save_file, 
                                            ".txt")))

                  } else if(loo == TRUE){
                    
                    
                      lgb.save(model_lgbm, 
                              filename =
                                here("data", "models", "loocv", 
                                     paste0("final_tuned_lgbm_for_",
                                            save_file, 
                                            ".txt")))
                    

                  }       
                
              } else if (solo == TRUE){
                
              
                                    
                      lgb.save(model_lgbm, 
                               filename =
                                here("data", "models", "solo", 
                                     paste0("final_tuned_lgbm_for_",
                                            training_data$tributary[1],
                                            save_file, 
                                            ".txt")))

                } 
              
            }
            
            ### Get shapley values
            
              
            shap_values <- SHAPforxgboost::shap.prep(xgb_model =  model_lgbm, 
                                                           X_train = preds)
            
          stats_and_importance[[4]] <- shap_values
              
              
              
              
            shap_values2 <- SHAPforxgboost::shap.prep(xgb_model =  model_lgbm, 
                                                           X_train = test_lgbm)
            
            stats_and_importance[[6]] <- shap_values2
            
            just_shaps <- shap.values(xgb_model =  model_lgbm, X_train = test_lgbm)
              
              
            stats_and_importance[[7]] <- just_shaps$shap_score

            
            stats_and_importance[[8]] <- tree_table
            
            ### Predict with the model
            fits <- predict(model_lgbm,
                            data = preds) %>%
              as_tibble() %>% rename(log_predicted_conc = 1)
            
            predicted <- predict(model_lgbm, 
                                           data = test_lgbm) %>%
              as_tibble() %>% rename(log_predicted_conc = 1)
            
    
           #### Bind to observations (training data) to estimate residuals
            ### and calculate the smearing factor
          predicted_observed_smear <- bind_cols(chosen_model_train %>%
                                            dplyr::rename(log_observed_conc = log_conc),
                                                fits) %>%
          mutate(exp_log_model_residuals = 10^(log_observed_conc - 
                                                 log_predicted_conc))
          
          stats_and_importance[[9]] <- predicted_observed_smear
          
             ### Estimate smearing coefficient  
        D <- mean(predicted_observed_smear$exp_log_model_residuals)
          
        #### Bind predictions to test values
          predicted_observed <- bind_cols(chosen_model_test %>%
                                            dplyr::rename(log_observed_conc = log_conc),
                                                predicted) 
        
     stats_and_importance[[5]] <- tibble(D_fact = D)
        
        ### Now use the smearing coefficient to convert back to non-log
        predicted_observed_err <- predicted_observed %>%
              mutate(predicted_conc = (10^log_predicted_conc)*D) %>%
              mutate(observed_conc = 10^log_observed_conc) %>%
              mutate(raw_err = predicted_conc - observed_conc) %>%
              mutate(sqrerr = raw_err^2,
                      abs_error = abs(raw_err),
                      abs_pct_error = abs_error/predicted_conc)
        
        
        
    ###########################################################################
    ############################################################################
     
    ################# EVALUATE MODEL ##########################################
        
    model_stats_by_trib_lgbm <- predicted_observed_err %>%
      dplyr::group_by(tributary) %>%
      summarise(rmse = sqrt(mean(sqrerr)),
                mae = mean(abs_error),
                nse = hydroGOF::NSE(predicted_conc, 
                                    observed_conc),
                kge = hydroGOF::KGE(predicted_conc, 
                                    observed_conc),
                pbias = hydroGOF::pbias(predicted_conc,
                                        observed_conc),
                r = cor(predicted_conc, observed_conc),
                br2 = hydroGOF::br2(predicted_conc, observed_conc)) 
       
       
       stats_and_importance[[1]] <- model_stats_by_trib_lgbm
        
        #### Calculate variable    
        chosen_lgbm_var_imp <- lgb.importance(model_lgbm , percentage = TRUE)
        
        stats_and_importance[[2]] <- chosen_lgbm_var_imp
        
        #### Save the output timeseries
        stats_and_importance[[3]] <- predicted_observed_err
        
        #### Return what we want
        return(stats_and_importance)
        
          
} ### End function
```

## Test ungaged models

This function essentially makes it easier to train and evaluate LightGBM
models with our selected features and hyperparamters in an **ungaged**
scenario, using the lgbm_runner function above. It is basically a
wrapper around lgbm_runner for an ungaged scenario.

``` r
#### This function takes a bunch of arguments
#### 1) constituent_df: this is the training data with concentration and 
#### all our selected features
#### 2) test_df: This is the concentration data we have not shown to the model at all yet 
#### 3) picked_vars: these are the features we have selected to be included in our
#### final model
#### 4) do_save: boolean for whether the model architecture should be written to disk
#### this is essential to do if we want to use trained models for future applications 
#### (like forecasting)
#### 5) is_tuned: boolean for whether or not to use tuned or default hyperparameters
#### 6) constit: which water quality constituent it is?

watershed_loo_runner <- function(constituent_df = NULL ,
                                 test_df = NULL,
                                 picked_vars = NULL,
                                 tuned_params = NULL,
                                 do_save = FALSE,
                                 is_tuned = FALSE,
                                 constit = NULL) {
  
       ### Make a bunch of empty lists to save various outputs 
        
        summary_stats_overall <- list()
        
        summary_stats_each <- list()
        
        var_imp_each <- list()
        
        each_pred_obs_ts <- list()
        
        each_trib_test_shap_values <- list()
        
        each_trib_train_shap_values <- list()
        
        each_just_shaps <- list()
        
        each_each_shap_plots <- list()
        
        each_d_fact <- list()
        
        loo_output <- list()
        
        tree_tables <- list()
        
        in_sample_fits <- list()
        
        
        ###### Dunno 
        
        cnst <- str_extract(deparse(substitute(constituent_df)), "^[^_]*")
        
        #### Loop over each watershed to train the model on 17 of 18
        #### And then test it on the remaining 1
        
        for(i in 1:18) {
          
          #### Track our progress
          
            tr <- constituent_df %>% 
                filter(group_id == i) %>%
                dplyr::slice(1) %>%
                .$tributary
              
            cat(crayon::cyan("\nTesting Trib", tr, "\n")) 
              

            ### Now run the final model with the tuned parameters 
            
            tuned_model <- lgbm_runner(training_data = constituent_df %>%
                                           filter(group_id != i) %>%
                                           dplyr::select(!group_id), 
                            testing_data =  test_df %>%
                                           filter(group_id == i) %>%
                              dplyr::select(!group_id), 
                            picked_vars,
                            save = do_save,
                            save_file = paste0(constit, "_", 
                                               str_replace_all(tr, " ", "_")),
                            tuned = is_tuned,
                            tuned_params = tuned_params,
                            loo = TRUE)
            
            cat(crayon::green("Done Training/Testing √ \n")) 
            
            #### Calculate summary stats 

            summary_stats_each[[i]] <-  tuned_model[[1]]
            
            cat(crayon::green("Done Summary Stats √ \n")) 
            
            current_trib <- summary_stats_each[[i]] %>%
              dplyr::ungroup() %>%
              dplyr::slice(1) %>%
              .$tributary
            
            ### And variable importance
            
            var_imp_each[[i]] <- tuned_model[[2]] %>% 
              as_tibble() %>%
              mutate(tributary = current_trib)
            
            cat(crayon::green("Done Var Imp √ \n")) 
              
            #### And in-sample fits (predictions on training data)
            
            in_sample_fits[[i]] <- tuned_model[[9]] %>%
                  mutate(testing_trib = tr)
        
            #### And final time series
            
            each_pred_obs_ts[[i]] <- tuned_model[[3]]  %>%
              mutate(tributary = current_trib)
            
            cat(crayon::green("Done Trib Importance √ \n")) 
            
            ### And shap values
            
            test_shap_values <- tuned_model[[6]] 
            
            train_shap_values <- tuned_model[[4]] 
            
            each_just_shaps[[i]] <- tuned_model[[7]]
            
            each_trib_test_shap_values[[i]] <- test_shap_values %>%
              as_tibble() %>%
              mutate(tributary = current_trib)
            
            each_trib_train_shap_values[[i]] <- train_shap_values %>%
              as_tibble() %>%
              mutate(tributary = current_trib)
            
            cat(crayon::green("Shap Values Calculated √ \n")) 
            
            #### And parameters related to constuction of the trees
            
            tree_tables[[i]] <- tuned_model[[8]] %>%
              mutate(testing_trib = tr)
            
            
            #### And SHAP plots

            each_shap_plots_test <- shap.plot.summary(test_shap_values)
            
            each_shap_plots_test <- each_shap_plots_test +
              labs(title = current_trib)
            
            
            each_shap_plots_train <- shap.plot.summary(train_shap_values)
            
              each_shap_plots_train <- each_shap_plots_train +
              labs(title = paste0("Left Out: ", current_trib))
            
            
            cat(crayon::green("Shap Plots Plotted √ \n")) 
            
            
            ggsave(here("plots",paste0(
            cnst, "_testing_on_", current_trib, "_v2.jpg")),
                   each_shap_plots_test,
                   height = 5, width = 5, dpi = 300)
            
            ggsave(here("plots",paste0(
                          cnst,
                          "_leaving_out_", current_trib, "_v2.jpg")),
                   each_shap_plots_train,
                   height = 5, width = 5, dpi = 300)
            
                  cat(crayon::green("Shap Plots Saved √ \n")) 
        
            #### And smearing factor from final training
                  
            each_d_fact[[i]] <- tuned_model[[5]] %>%
              mutate(trib = tr)
            
                  cat(crayon::green("Smearing factor Saved √ \n")) 
                  
                  #cat(crayon::bgBlack("On to the next one")) 
        
        }
        
        
        ### Bind everything together
        
        summary_stats_each_combined <- bind_rows(summary_stats_each)
        
        loo_output[[1]] <- summary_stats_each_combined
        
        var_imp_each_combined <- bind_rows(var_imp_each)
        
        loo_output[[2]] <- var_imp_each_combined
        
        each_pred_obs_ts_combined <- bind_rows(each_pred_obs_ts)
        
        loo_output[[3]] <- each_pred_obs_ts_combined
        
        each_trib_test_shap_values_combined <- bind_rows(each_trib_test_shap_values)
        
        loo_output[[4]] <- each_trib_test_shap_values_combined
        
        each_trib_train_shap_values_combined <- bind_rows(each_trib_train_shap_values)

        loo_output[[6]] <- each_trib_train_shap_values_combined
        
        each_d_fact_combined <- bind_rows(each_d_fact)
        
        loo_output[[5]] <- each_d_fact_combined
        
        loo_output[[7]] <- (each_just_shaps)
        
        loo_output[[8]] <- bind_rows(tree_tables)
        
        loo_output[[9]] <- bind_rows(in_sample_fits)
        
        return(loo_output)

}
```

## New Get GCP URLS

``` r
### We need to re-write some functions from the nwmTools package
### Here we re-write the geturls command to better reflect 
### the URLS of the medium-term forecasts (the built-in function had syntax
### more reflective of the short-term urls)

get_gcp_urls2 <- function (config = "short_range", 
                                domain = "conus", date, hour = "00",
                                minute = "00", num, ensemble = NULL, 
                                output = "channel_rt") 
{
    meta = nwm_filter(source = "gcp", version = NULL, config = config, 
        date = date, ensemble = ensemble, output = output, domain = domain)
    dates = seq.POSIXt(as.POSIXlt(paste(date, hour, minute), 
        tz = "UTC"), length.out = num+1, by = "1 hour")
    YYYYDDMM = rep(format(as.POSIXct(date, tz = "UTC"), "%Y%m%d"), num)
    forward = sprintf("%03d", seq(1, num))
    urls = glue(meta$http_pattern, YYYYDDMM = YYYYDDMM, config = meta$config[1], 
        HH = hour, foward = forward, output = meta$output[1], 
        domain = meta$domain, ensemble = meta$ensemble[1], prefix = meta$prefix[1])
    dates = ymd_hm(paste0(date[1], hour, "00")) + hours(1:(num))
    data.frame(dateTime = dates, urls = urls, output = output)
} 

#### Must do this so that some functions internal to nwmTools can be used

environment(get_gcp_urls2) <- environment(get_gcp_urls)



###################################################################################
```

## File Reader

``` r
#### Function to read downloaded NetCDF files and retrieve various bits of 
#### information from those files 

file_reader <- function(file = NULL, ids = NULL) {
  
channel_params <- tidync(file) %>%
  activate("D0") %>%
  hyper_tibble() %>%
  filter(feature_id %in% ids)

channel_plus_time <- channel_params %>%
  mutate(feature_id  = as.character(feature_id)) %>%
  as_tibble() 

return(channel_plus_time)
  
  
}
```

## Temp File Maker

``` r
### Function to generate the temporary files to store downloaded NWM forecasts

temp_file_maker <- function(urls = NULL) {
  
 urls %>%
    mutate(date_time_mem = paste0(init_date_time, "_", init_time, "_", 
                                  member)) %>%
      mutate(filename = tempfile(pattern = paste0(date_time_mem,"plot"), 
                             fileext = ".nc")) %>%
   dplyr::select(!c(date_time_mem)) 
  
  
}
```

## NWM Downloader w/ cloud-based operations

``` r
### This function dowloads archived NWM forecasts from the Google Bucket
### where they are stored. It utilizes cloud-based operations to extract
### the reaches we are interested in BEFORE downloading the entire NWM file
### to our local machine. In this way we avoid downloading giant files that 
### contain largely extraneous data. This is perhaps the more "correct" way to 
### do this. 
### Note that this is a rewrite of the nwmTools get_timeseries function within 
### Mike Johnson's nwmTools package. We are heavily indebted to the great work
### Mike has done. Thanks Mike!!!!!! (https://github.com/mikejohnson51)
### Also, *******It is NOT FAST*******
### Perhaps it could be rewritten for efficiency, but hey, it works



get_timeseries3 <- function(urls = NULL, 
                            ids = NULL,
                            outfile = NULL,
                            index_id = "feature_id",
                            varname = "streamflow"

                            
) {
  
            #### Get values function #######################################################
            get_values <- function(url, var, ind = NULL) {
                        v = suppressWarnings(values(terra::rast(glue("{url}://{var}"))))
                        
                        if (!is.null(ind)) {
                            v = v[ind]
                        }
                        else {
                            v = v
                        }
                        
                        return(v)
            }
            ################################################################################


### Get lead time from URL 
### And get init_date from URL
            
lead_time <- str_extract(str_extract(urls, "f\\d+"), "\\d+")

init_date <- as_date(str_extract(str_extract(urls, ".\\d+"), "\\d+"))

member <- str_extract(urls, "mem\\d+")

#### Little updater

print(paste("Downloading", init_date, lead_time, member, " "))

### First, set up a URL that turns the netCDF on the Google Bucket
### Into a HD5 file
### And utilizes external pointer (vsicurl) to access file "on the fly"

nwm_url2 <- glue("HDF5:\"/vsicurl/{urls}\"")


### Now get the feature_ids (comids) from the file

all_ids <- get_values(url = nwm_url2, index_id)


### Now find the indexes of our comids (reaches) of interest
### In the file that contains all the comids
### The thinking here is that the index of a given comid in the feature_id "layer"
### Should be the same index of where the streamflow value corresponding to that comid
### Is located
### We need to put in a bunch of safety if-else statements to keep from breaking 
### if the file is not found for whatever reason

 
   if (!is.null(index_id)) {
     
     
     
      all_ids = get_values(url = nwm_url2, index_id)
                
     ### If no COMIDs are entered, return every streamflow ID
     ### So all 2.7 million reaches in the NWM
     ### But if particular COMIDs are desired
     ### Find the index of those COMIDs
     ### in the feature_id "layer"
     
        if (is.null(ids)) {
            ind = NULL
            ids = all_ids
        }
     
        else {
         
                ind <- which(all_ids %in% ids)
            
        }
   }

    else {
      
        ind = NULL
        ids = NULL
        
    } ### Returns null if error in reading files



#### Now let's get streamflow
#### It comes in without decimals (not sure why)

q_values <- get_values(nwm_url2, varname, ind)

q_values <- q_values/100

#### And time 
#### Which we need to multiply by 60
#### to turn into "true" POSICxt time (seconds from 1970-01-01)

nwm_time <- get_values(nwm_url2, "time")[1]*60


#### Now, we have to actually extract the feature_id in the same way as we did
#### for discharge
#### This gets us a vector ordered in the same order as our discharge vector
#### Without this, we will bind things that our in different orders
#### and our final output dataframe will be meaningless 
#### THIS IS EXTREMELY IMPORTANT

comids <- get_values(nwm_url2, index_id, ind)

### Bind the time, COMIDs, and modeled streamflow values into one tibble
### Make into a tibble

streamflow_by_reach <- tibble(comid = comids, 
                              init_date = init_date, 
                              lead_time = lead_time,
                              member = member,
                              predict_dateTime = as_datetime(nwm_time),
                              modeled_q_cms = q_values
                              )



### And write that to file 

write_csv(streamflow_by_reach, outfile,
          append = TRUE,
          col_names = !file.exists(outfile),
          )


#return(streamflow_by_reach)

}
```

## Transform NWM forecasts to daily

This function takes NWM forecasts, which are hourly, and transforms them
to daily by simply taking a mean of the forecasted discharge

``` r
#### Make daily

nwm_daily_maker <- function(nwm_df) {
  
  nwm_df_daily <- nwm_df %>%
    mutate(predict_date = as_date(predict_dateTime)) %>%
    dplyr::group_by(comid, init_date, member, predict_date) %>%
    summarise(forecasted_q_cms = mean(modeled_q_cms)) %>%
    dplyr::ungroup()
  
  return(nwm_df_daily)
  
  
  
}
```

## Calculate forecast antecedent conditions

This function essentially calculates various antecedent conditions to be
used in a forecasting scenario. Importantly, this requires combining
observations AND forecasted values for a given forecast timestep,
meaning that, for example, forecasted flow for Oct 31st 2022 has a
weekly antecedent discharge of seven observed values for a forecast with
lead time of zero days (so the one issued Oct 31st at midnight); of six
observed and one forecasted for a forecast with a lead time of one day
(so issued at Oct 30th); and so on. This is kinda complicated to create
a dataframe that follows this format

``` r
antecedent_calculator <- function(predict_date,
                                  init_date,
                                  site_no,
                              days_from_gage = NULL, 
                              days_from_model = NULL,
                              monthly_days_from_gage = NULL,
                              monthly_days_from_model = NULL,
                              #member = NULL,
                              obs_df = NULL,
                              model_df = NULL,
                              forecast_model = "NWM"
                              ) {

  
########################## SET UP ##############################################
  
  ### Find the index of the dateTime of interest in the observed data (minus 1)
  
  date_index_obs <- which(obs_df$date == predict_date &
                                obs_df$site_no == site_no) 
  
  ### Find the index of thepredict  date
  
  date_index_model <- which(model_df$predict_date == predict_date & 
                              model_df$site_no == site_no & 
                              model_df$init_date == init_date
                              )
  
  ########## WEEKLY ANTECEDENTS ###################################################
  
  ### Calculate weekly antecedents based on values calculated in the 
  ### timesteps_needed step
  
    if (days_from_gage != 0){
          log_obs_flow <- obs_df %>%
            dplyr::ungroup() %>%
            dplyr::slice(((date_index_obs - days_from_model) - days_from_gage):
                           (date_index_obs - days_from_model - 1)) %>%
            .$log_daily_q
          
        
    } else{
        
      log_obs_flow <- NA
      
      } ### End log observed flow ifelse 

    if(days_from_model != 0 ){
      
          log_modeled_flow <- model_df %>%
            dplyr::ungroup() %>%
            dplyr::slice((date_index_model - days_from_model):(date_index_model -
                                                                  1)) %>%
            .$mean_forecasted_log_q_cms_km2
        
      
        }  else {
        
        log_modeled_flow <- NA
      }

    obs_and_modeled_flow <-  c(log_obs_flow, log_modeled_flow) %>%
      na.omit() 
      
      mean_prior_weekly_flow <- mean(obs_and_modeled_flow) %>%
            as_tibble() %>%
            rename(mean_prior_weekly_log_q_cms_km2 = 1)
  

  ################################################################################
    
  ############# MONTHLY ANTECEDENTS ##############################################

    if (monthly_days_from_gage > 0 & monthly_days_from_model <= 7){

          log_obs_flow_minus3 <- obs_df %>%
            dplyr::slice((date_index_obs - 30):(date_index_obs -
                                                       7)) %>%
            .$log_daily_q

          obs_and_modeled_flow_monthly <- c(obs_and_modeled_flow,
                                            log_obs_flow_minus3) %>%
            na.omit()

    } else if(monthly_days_from_model > 7) {

      log_modeled_flow_monthly <- model_df %>%
            dplyr::slice((date_index_model -
                            monthly_days_from_model):(date_index_model - 1)) %>%
            .$mean_forecasted_log_q_cms_km2

      log_observed_flow_monthly <- obs_df %>%
            dplyr::slice(((date_index_obs - monthly_days_from_model)
                          - monthly_days_from_gage):(date_index_obs -
                                                        monthly_days_from_model - 1)) %>%
            .$log_daily_q

      obs_and_modeled_flow_monthly <- c(log_modeled_flow_monthly,
                                        log_observed_flow_monthly) %>%
        na.omit()

    } else {

      obs_and_modeled_flow_monthly <- NA

    }
      
            mean_prior_monthly_flow <- mean(obs_and_modeled_flow_monthly) %>%
            as_tibble() %>%
            rename(mean_prior_monthly_log_q_cms_km2 = 1)
            
            weekly_and_monthly_ant <- bind_cols(mean_prior_weekly_flow, 
                                                mean_prior_monthly_flow)


}
```

## Calculate Delta Q

This function calculates the change in discharge from time t-1 to time
t. It calculates this as both a raw numeric and then a categorical
variable with three possible values of -1, 0, and 1.

``` r
limb_getter <- function(model_df, observed_df){
  
  
  model_df %>%
    mutate(init_date_minus_one = init_date - days(1)) %>%
    group_by(init_date, tributary, site_no, comid) %>% #### This is to calculate limb
    mutate(delta_q_modeled = mean_forecasted_q_cms_km2  - 
             lag(mean_forecasted_q_cms_km2 )) %>%
    mutate(modeled_lag_daily_q = lag(mean_forecasted_log_q_cms_km2
                                     )) %>%
    mutate(delta_q_q_modeled = 
             abs(delta_q_modeled/lag(mean_forecasted_q_cms_km2))*sign(delta_q_modeled)) %>%
    inner_join(., observed_flow_for_nwm_predictions %>%
                 dplyr::select(tributary, site_no,
                               date, log_daily_q),
               join_by(init_date_minus_one == date,
                       tributary == tributary,
                       site_no == site_no)) %>%
    mutate(delta_q_day_zero = mean_forecasted_q_cms_km2 - 10^log_daily_q) %>%
    mutate(delta_q_q_day_zero = abs(delta_q_day_zero/10^log_daily_q)*sign(delta_q_day_zero)) %>%
    mutate(delta_q_q = coalesce(delta_q_q_modeled, delta_q_q_day_zero)) %>%
    group_by(init_date, tributary, site_no, comid) %>% #### This is to calculate limb
    mutate(delta_q_modeled = mean_forecasted_log_q_cms_km2 - 
             lag(mean_forecasted_log_q_cms_km2)) %>% ## Update
    mutate(delta_q_day_zero = mean_forecasted_log_q_cms_km2 - 
             log_daily_q) %>% ### Replace w/ log
    mutate(delta_q = coalesce(delta_q_modeled, delta_q_day_zero)) %>%
    rename(delta_daily_q = delta_q) %>%
    mutate(lag_daily_q = coalesce(modeled_lag_daily_q, log_daily_q)) %>%
    dplyr::ungroup() %>%
    mutate(delta_daily_q_cat = case_when(delta_q_q > 0.10 ~ 1,
                                       delta_q_q < -0.10 ~ -1,
                                       (delta_q_q <= 0.10 & delta_q_q >= -0.10) ~ 0)) %>%
    mutate(delta_daily_q_cat = as.factor(delta_daily_q_cat)) %>%
    dplyr::select(c(predict_date, init_date, tributary,
                    delta_daily_q, 
                    delta_daily_q_cat,
                    lag_daily_q))
    
  
}
```

## LightGBM Forecaster

This function makes forecasts of water quality concentration using
previously trained LightGBM models and NWM streamflow forecast data. The
function takes a variety of arguments:

1)  `predictors_df` A dataframe containing all the data to be used for
    forecasting. This is, essentially, NWM flow and derived antecedent
    conditions, as well as whatever watershed characteristics we have
    included in the model during the development process

2)  `observed_df`: A dataframe containing observed data, which we will
    use for evaluation

3)  `D_factor`: dataframe containing the smearing factor calculated
    during the development process

4)  `constituent`: which constituent we are forecasting

5)  `solo_or_lumped_or_loocv`: a string indicating which type a model
    scenario the forecasting is to be done for

6)  `nwm`: Boolean indicating in the streamflow forecast source is the
    National Water Model or something else

``` r
lgbm_predictor <- function(predictors_df = NULL ,
                           observed_df = NULL,
                           D_factor = NULL,
                           constituent = NULL,
                           solo_or_lumped_or_loocv = "solo",
                           nwm = TRUE,
                           model_file_location = NULL){
  
  model_files <- list.files(model_file_location)
  
  #### Declare the list we're going to retun things in
  nwm_error_and_full_ts <- list()
  
      #### Get the name of the tributary we are on 
  #### For a given iteration
  
  if(nwm == TRUE) {
    
      trib <- observed_df$tributary[1]
    
  } else if(nwm == FALSE){
    
    trib = observed_df
    
  }

  
  
  #### Addd to the forecasted df
  
  if(solo_or_lumped_or_loocv == "solo" |
     solo_or_lumped_or_loocv == "loocv") {
    
    if(nwm == TRUE){
      
      forecasting_data <- predictors_df %>%
        mutate(tributary = trib) %>%
        mutate(lead_days = as.numeric(predict_date - init_date)) %>%
        mutate(lead_days = as.factor(lead_days))
      
    } else if (nwm == FALSE) {
      
      forecasting_data <- predictors_df %>%
        mutate(tributary = trib)
      
    }
    
  } else if (solo_or_lumped_or_loocv == "lumped") {
    
    if(nwm == TRUE){
      
          forecasting_data <- predictors_df %>%
       mutate(lead_days = as.factor(lead_days))
      
      
    } else if(nwm == FALSE) {
      
      forecasting_data <- predictors_df
      
    }
    

    
  }

  ### Find the relevant file name
  ### For the trained model

  if(solo_or_lumped_or_loocv == "lumped") {
    
      lgbm_model <- tibble(file = model_files) %>%
    filter(str_detect(file, paste0(constituent))) %>%
    .$file
    
    
  } else if(solo_or_lumped_or_loocv == "loocv") {
    
      
    lgbm_model <- tibble(file = model_files) %>%
        filter(str_detect(file, paste0(constituent, "_", 
                                       str_replace_all(trib, " ", "_")))) %>%
        .$file
    
    
  }



  if(solo_or_lumped_or_loocv == "solo") {
    
    
      #### Make the test data matrix
    
  forecast_data <- data.matrix(forecasting_data %>%  
                               dplyr::select(water_year, season, log_daily_q,
                                                      delta_daily_q,
                                             #delta_daily_q_cat,
                                                      mean_prior_weekly_q, 
                                             mean_prior_monthly_q,
                                                      lag_daily_q,
                                                 ))
  } else if(solo_or_lumped_or_loocv == "lumped" |
            solo_or_lumped_or_loocv == "loocv") {
    
     forecast_data <- data.matrix(forecasting_data %>%
                                  dplyr::select(!c(init_date, 
                                                   lead_days, 
                                                   predict_date,
                                                   tributary)) %>%  
                               dplyr::select(order(colnames(.))))
    
    
  } 


  #### Load in the tuned LGBM model
  
        cat(crayon::yellow(paste("Loading model ", lgbm_model, "\n")))

    
        final_tuned_model <- lgb.load(paste0(model_file_location,
                                             "/",
                                        lgbm_model))
      

  
  #### For multiple observations on the same day
  #### Take the mean of those 
  
  if(nwm == TRUE) {
    
        observed_df <- observed_df %>%
            group_by(tributary, date) %>%
            mutate(log_conc = mean(log_conc)) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()

  }
  
  
  ### Forecast with the model, bind to other relevant data
  
        forecast_from_nwm <- predict(final_tuned_model,
                                       data = forecast_data) %>%
          as_tibble() %>%
          rename(log_forecasted_conc = 1) %>%
          bind_cols(., forecasting_data %>%
                      rename_with(~paste0(., "_forc"), where(is.numeric))) %>%
          rename(log_modeled_flow =
                      log_daily_q_forc)

    ### Now bind to observations
        
        if(nwm == TRUE){
          
          
        
          forecast_and_obs <- forecast_from_nwm %>%
            inner_join(., observed_df %>%
                   dplyr::select(tributary,
                                 date,
                                 log_conc),
                   join_by(predict_date == date,
                           tributary == tributary)) %>%
          rename(log_observed_conc = log_conc) 

        ### Smearing

        ### Transform back to linear
            forecast_and_obs <- forecast_and_obs %>%
                    mutate(forecasted_conc = D_factor*10^log_forecasted_conc) %>%
                    mutate(observed_conc = 10^log_observed_conc)

            nwm_error_and_full_ts[[1]] <- forecast_and_obs %>%
              mutate(constituent = constituent)

            ### Calculate some errors
            
           nwm_forecast_error <-  forecast_and_obs %>%
                    mutate(raw_error = forecasted_conc - observed_conc,
                           sqrerr = (raw_error)^2,
                           abs_error = abs(raw_error)) %>%
                    dplyr::group_by(tributary, lead_days) %>%
                    summarise(rmse = sqrt(mean(sqrerr)),
                      mae = mean(abs_error),
                      kge = hydroGOF::KGE(forecasted_conc,
                                          observed_conc),
                     nse = hydroGOF::NSE(forecasted_conc,
                                          observed_conc),
                      pbias = hydroGOF::pbias(forecasted_conc,
                                              observed_conc),
                      br2 = hydroGOF::br2(forecasted_conc, observed_conc),
                     decompose_kge(forecasted_conc, observed_conc)) %>%
             dplyr::select(!c(r_term, variability_term, bias_term))

            nwm_error_and_full_ts[[2]] <- nwm_forecast_error %>%
              mutate(constituent = constituent)
            
            forecasted_contin_ts <- forecast_from_nwm %>%
              mutate(forecasted_conc = 10^log_forecasted_conc*D_factor) %>%
              mutate(lead_days = as.numeric(predict_date - init_date)) %>%
              dplyr::select(tributary,
                            init_date, predict_date, lead_days,
                            forecasted_conc)
            
            
        } else if(nwm == FALSE){
          
            forecasted_contin_ts <- forecast_from_nwm %>%
              mutate(forecasted_conc = 10^log_forecasted_conc*D_factor) %>%
              {if(solo_or_lumped_or_loocv == "solo") dplyr::select(., 
                                                          tributary,
                                                          date,
                                                          forecasted_conc) 
                else dplyr::select(., 
                                   tributary,
                                   init_date,
                                   forecasted_conc) }
          
          
          
        }

            nwm_error_and_full_ts[[3]] <- forecasted_contin_ts %>%
              mutate(constituent = constituent)

            return(nwm_error_and_full_ts)
            
            
}
```

## Decompose KGE

``` r
#### This function decomposes KGE into its constituent parts
#### of bias, variability, and correlation (r)
#### See Gupta et al. (2009) for more information
#### (https://doi.org/10.1016/j.jhydrol.2009.08.003)


decompose_kge <- function(sim, obs) {
  
  #### Calculate the correlation between modeled and observed
  
  r <- cor(sim, obs, method = "pearson", use = "pairwise.complete.obs")
  
  #### And the bias
  
  bias <- mean(sim)/mean(obs)
  
  #### And the variability 

  variability <- sd(sim)/sd(obs)
  
  #### Now calculate the terms as they appear in eqn. 9 in Gupta et al. (2009)
  #### Which are used to calculate KGE 
  
  r_term <- (r-1)^2
  
  variability_term <- (variability-1)^2
  
  bias_term <- (bias-1)^2
  
  #### And the final KGE calculation
  
  kge_manual <- 1 - sqrt(r_term + variability_term + bias_term)
  
  #### Finally, make a summary table of each 
  
  decomposed_kge <- tibble(r, variability, bias,
                           r_term, variability_term, bias_term,
                           kge_manual)
  
  
  return(decomposed_kge)
  
}
```
