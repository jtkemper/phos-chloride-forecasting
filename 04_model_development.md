04_model_development
================
JTK
2025-01-07

################################################################################ 

This script builds boosted regression trees (specifically, the LightGBM
implementation) to forecast total phosphorus and chloride concentration
in 18 tributary watersheds to Lake Champlain (see Fig. 1 in the
[ReadMe](https://github.com/jtkemper/phos-chloride-forecasting/blob/main/README.md))

This script builds two model “types” for each constituent:

1)  a model designed to predict/forecast water quality concentration in
    a **“gaged”** scenario, where data from a given watershed is used to
    train the model and predict in that same watershed

AND

2)  a model designed to predict/forecast water quality concentration in
    an **“ungaged”** scenario, where data from all watershed save one is
    used to train a model to predict concentration in that left-out
    watershed

Models are trained on observational data and then predictions are tested
on “test” datasets which have been wholly witheld from training. This,
in essence, provides an upper benchmark for model performance in terms
of forecasting ability: how well can these models do in predicting
(rather than forecasting) concentration if discharge is “perfect”

The basic workflow in this script is as follows: we use a backwrds
variable selection method to identify the most relevant variables for
each LightGBM model, we then train the model on those variables and tune
the hyperparameters, and then finally train the model with tuned
hyperparameters. The ultimate goal is to develop models that can be fed
streamflow forecast data and forecast future in-stream concentrations.

Model evaluation is done is subsequent scripts, as is forecasting (and
forecast evaluation).

**Inputs**

1)  `tp_drivers` (*from 03_observational_data_prep*): dataframe with
    total phosphorus concentration data and all possible predictive
    features

2)  `chlor_drivers` (*from 03_observational_data_prep*): dataframe with
    chloride concentration data and all possible predictive features

**Outputs**

1)  Model architectures for chloride and total phosphorus prediction in
    gaged & ungaged scenarios (`tuned_tp_lgbm`, `tuned_chlor_lgbm`,
    `tp_loo`, `chlor_loo`). Also output text files with these saved
    architecture to xxx/data/model/lumped and xxx/data/models/loocv.

2)  Performance statistics for total phosphorus and chloride
    concentration predictions on test data
    (`summary_stats_final_model_tp_gaged`,
    `summary_stats_final_model_chlor_gaged`, `sum_stats_tp_ungaged`,
    `sum_stats_chlor_ungaged`).

3)  Predicted time series for total phosphorus and chloride for test
    years (2019-2023) (`tp_pred_obs_ts`, `chlor_pred_obs_ts`,
    `predicted_vs_observed_tp_ungaged`,
    `predicted_vs_observed_chlor_ungaged`)

################################################################################ 

# Housekeeping

### Packages

``` r
### Data mgmt
require(tidyverse)

### Model development
require(hydroGOF)
require(lightgbm)
require(SHAPforxgboost)

### Model tuning
require(recipes)
require(parsnip)
require(tune)
require(dials)
require(workflows)
require(yardstick)
require(bonsai)
require(tidymodels)

### Data viz
require(ggthemes)
```

### Load prior scripts

``` r
source(knitr::purl(here("Rmd-files/00_functions.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/01_data_discovery_and_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/02_watershed_attributes_download.Rmd"), quiet=TRUE))

source(knitr::purl(here("Rmd-files/03_observational_data_prep.Rmd"), quiet=TRUE))
```

# Variable selection

Here we perform a backwards variable selection that essentially
calculates performance for a model with all features, then all features
minus one, then all features minus two, etc. Feature removal is done by
using mean SHAP importance and the distribution of SHAP values to rank
feature importance, removing the least important is removed, re-running
the models, and repeating until only one variable remains. The ultimate
goal is to balance performance and parsimony, so model selection is done
by examining model performance as a function of number of variables
contained within the model. The selected model is the model with the
least amount of features before performance begins to decline
precipitously (if that does indeed occur).

We will use a custom function to do the backwards variable selection for
each constituent. This function requires inputs of: 1) the data the
model should be trained on

2)  a boolean variable of whether the user wishes to perform variable
    selection

3)  a boolean for how performance statistics should be returned (grouped
    by trib or not)

4)  what model is to be used (not relevant here, really)

5)  whether or not the model hyperparameters have been tuned

6)  the “scenario” for which models should be trained - ungaged or gaged

7)  how variable importance should be ranked - using SHAP values or
    traditional feature importance

We’ll do this seperately for models to be deployed in a gaged or ungaged
scenario

### Backwards Variable Selection - Gaged Scenario

#### Total Phosphorus

``` r
###############################################################################

################################################################################

#### Perform variable selection for TP

set.seed(913)


tp_all_models <- lgbm_selector(tp_drivers,
                               selection = TRUE,
                               is_tuned = FALSE,
                               model_type = "gaged",
                               selector = "shap")


#### See the performance stats for each model run

tp_all_models_stats <- tp_all_models[[2]]

#### Plot performance stats as a function of model "number"
#### where model number 1 has all n variables and model n has 1 variable

plot_stats(tp_all_models_stats)

#### Select the model that most balances performance and parsimony 
#### So the one with the least amount of features before performance starts
#### to observably "drop off"

tp_models_we_like <- tp_all_models[[1]] %>%
  filter(model == 70)

#### Rename variables, which will be useful later

picked_vars_tp <- tp_models_we_like %>%
  mutate(Feature = dplyr::if_else(Feature == "water_year2", "water_year", Feature))


################################################################################
################################################################################
```

#### Chloride

``` r
#### Perform variable selection for Chlor

set.seed(913)


chlor_all_models <- lgbm_selector(chlor_drivers,
                               selection = TRUE,
                               is_tuned = FALSE,
                               model_type = "gaged",
                               selector = "shap")

#### See the performance stats for each model run

chlor_all_models_stats <- chlor_all_models[[2]]

#### Plot performance stats as a function of model "number"
#### where model number 1 has all n variables and model n has 1 variable

plot_stats(chlor_all_models_stats)

#### Select the model that most balances performance and parsimony 
#### So the one with the least amount of features before performance starts
#### to observably "drop off"

chlor_models_we_like <- chlor_all_models[[1]] %>%
  filter(model == 68)

#### Rename variables, which will be useful later

picked_vars_chlor <- chlor_models_we_like %>%
      mutate(Feature = if_else(Feature == "water_year2", "water_year", Feature))
```

### Backwards Variable Selection - Ungaged Scenario

#### Total Phosphorus

``` r
#### Perform variable selection for TP
#### Important: here we need to use the tp_train_valid dataframe
#### So that we can keep out our test data

set.seed(913)

tp_ungaged_models_all <- lgbm_selector(tp_train_valid, 
                               selection = TRUE,
                               is_tuned = FALSE,
                               model_type = "ungaged",
                               selector = "shap")

#### See the performance stats for each model run

tp_ungaged_models_stats <- tp_ungaged_models_all[[2]]

#### Plot performance stats as a function of model "number"
#### where model number 1 has all n variables and model n has 1 variable

plot_stats(tp_ungaged_models_stats)

#### Select the model that most balances performance and parsimony 
#### So the one with the least amount of features before performance starts
#### to observably "drop off"

tp_ungaged_models_we_like <- tp_ungaged_models_all[[1]] %>%
  filter(model == 68)

#### Rename variables, which will be useful later

picked_vars_tp_ungaged <- tp_ungaged_models_we_like %>%
      mutate(Feature = if_else(Feature == "water_year2", "water_year", Feature))
```

#### Chloride

``` r
#### Perform variable selection for Chlor

set.seed(913)



chlor_ungaged_models_all <- lgbm_selector(chlor_train_valid,
                               selection = TRUE,
                               is_tuned = FALSE,
                               model_type = "ungaged",
                               selector = "shap")

#### See the performance stats for each model run

chlor_ungaged_models_stats <- chlor_ungaged_models_all[[2]]

#### Plot performance stats as a function of model "number"
#### where model number 1 has all n variables and model n has 1 variable

plot_stats(chlor_ungaged_models_stats)

#### Select the model that most balances performance and parsimony 
#### So the one with the least amount of features before performance starts
#### to observably "drop off"

chlor_ungaged_models_we_like <- chlor_ungaged_models_all[[1]] %>%
  filter(model == 67)

#### Rename variables, which will be useful later

picked_vars_chlor_ungaged <- chlor_ungaged_models_we_like %>%
      mutate(Feature = if_else(Feature == "water_year2", "water_year", Feature))
```

# Tune selected models and choose best hyperparameters

Here, we tune the hyperparameters for the final models chosen in the
variable selection process. Tuned hyperparameters are as follows:

1)  num_leaves: the number of leaves in a given tree (which is the
    primary control on model complexity)

2)  min_data_in_leaf: the minimum number of data points that must fall
    in a single leaf. This is essential to mitigate over-fitting

3)  num_interations: this is the number of trees grown. Also important
    to prevent over-fitting

4)  sample_size: this is the fraction of data chosen to train any
    particular tree. Speeds up training and can prevent overfitting

5)  tree_depth: this explicitly tunes the depth to which trees can grow,
    limiting their maximum complexity and potentially mitigating
    overfitting.

For more about approaches to tuning LightGBM models, see here:
<https://lightgbm.readthedocs.io/en/latest/Parameters-Tuning.html>

For more about LightGBM hyperparameters, see here:
<https://lightgbm.readthedocs.io/en/v3.3.5/Parameters.html#feature_pre_filter>

To tune, we are going to use a few custom functions that create
training/validation splits and then tune the various parameters above.
Then we are going to select which of these is the best for our purposes
by running the top ten hyperparameters again.

### Gaged models

#### Tune the hypers

``` r
################################################################################
#################### TUNE THE MODEL ############################################

############### Total Phosphorus ###############################################


scores_tp <- lgbm_tuner(tp_drivers, picked_vars_tp, loocv = FALSE)


scores_tp_loocv <- lgbm_tuner(tp_drivers, picked_vars_tp, loocv = TRUE)


########## Chloride #############################################################

scores_chlor <- lgbm_tuner(chlor_drivers, picked_vars_chlor,
                           loocv = FALSE)

scores_chlor_loocv <- lgbm_tuner(chlor_drivers, picked_vars_chlor, loocv = TRUE)



################################################################################
```

#### Choose best hyperparameters (second pass)

``` r
############### Total Phosphorus ###############################################

#### Once tuning is done, we want to run each of the top ten returned hyperparameters
#### To see which actually yields the *best* performance in terms of our metrics of interest
#### This is bascially a second-pass filter for determining the best hyper parameters
#### We use a custom function to choose from our top ten

hype_tp <- hyper_chooser(constit_df = tp_drivers,
                         hypers_df = scores_tp,
                         selected_params = picked_vars_tp,
                         loocv = FALSE)



#### Now select the best performing one based on the hyperparameter ID and the 
#### performance metrics

#### View which iteration did best

hype_tp

##### Choose

best_params_tp <- scores_tp %>%
  arrange(mae_mean_plus_se_rank) %>%
  dplyr::slice(3)

best_params_tp




################################################################################

########## Chloride #############################################################

hype_chlor <- hyper_chooser(constit_df = chlor_drivers,
                         hypers_df = scores_chlor,
                         selected_params = picked_vars_chlor,
                         loocv = FALSE)

#### Now select the best performing one based on the hyperparameter ID and the 
#### performance metrics

##### View which iteration did best

hype_chlor

##### Choose

best_params_chlor <- scores_chlor %>%
  arrange(mae_mean_plus_se_rank) %>%
  dplyr::slice(7)
```

### Ungaged models

#### Tune the hypers

``` r
################################################################################
#################### TUNE THE MODEL ############################################

############### Total Phosphorus ###############################################


scores_tp_loocv <- lgbm_tuner(tp_drivers, picked_vars_tp, loocv = TRUE)


########## Chloride #############################################################


scores_chlor_loocv <- lgbm_tuner(chlor_drivers, picked_vars_chlor, loocv = TRUE)



################################################################################
```

#### Choose best hyperparameters (second pass)

``` r
############### Total Phosphorus ###############################################

#### Once tuning is done, we want to run each of the top ten returned hyperparameters
#### To see which actually yields the *best* performance in terms of our metrics of interest
#### This is bascially a second-pass filter for determining the best hyper parameters
#### We use a custom function to choose from our top ten

hype_tp_ungaged <- hyper_chooser(constit_df = tp_train_valid,
                         hypers_df = scores_tp_loocv,
                         selected_params = picked_vars_tp_ungaged,
                         loocv = TRUE)

#### Now select the best performing one based on the hyperparameter ID and the 
#### performance metrics

best_params_tp_ungaged <- scores_tp_loocv %>%
  arrange(mae_mean_plus_se_rank) %>%
  dplyr::slice(4)




################################################################################

########## Chloride #############################################################

hype_chlor_ungaged <- hyper_chooser(constit_df = chlor_train_valid,
                         hypers_df = scores_chlor_loocv,
                         selected_params = picked_vars_chlor_ungaged,
                         loocv = TRUE)

#### Now select the best performing one based on the hyperparameter ID and the 
#### performance metrics

best_params_chlor_ungaged <- scores_chlor_loocv %>%
  arrange(mae_mean_plus_se_rank) %>%
  dplyr::slice(5)
```

# Test models

Here, we are going to test models on data that has been totally left out
of the variable selection, tuning, and training process. These models
will represent the “upper” benchmark for forecast model performance -
how well water quality forecasts might do if the flow forecasts were
perfect.

### Gaged

#### Run the models

``` r
##### For Total Phosphorus

tuned_lgbm_tp <- lgbm_runner(tp_train_valid , 
                tp_test , 
                picked_vars_tp,
                save = TRUE,
                save_file = "Phosphorus_Total",
                tuned = TRUE,
                tuned_params = best_params_tp)

##### For Chloride

tuned_lgbm_chlor <- lgbm_runner(chlor_train_valid , 
                chlor_test , 
                picked_vars_chlor,
                save = FALSE,
                save_file = "Chloride",
                tuned = TRUE,
                tuned_params = best_params_chlor)
```

#### Extract relevant outputs

``` r
#############################################################################

#### Total Phosphorus

#### Calculate summary stats 

summary_stats_final_model_tp_gaged <- tuned_lgbm_tp[[1]]

#### Final time series

tp_pred_obs_ts <- tuned_lgbm_tp[[3]] 

### And shap values and plots

shap_values_tp_gaged <- tuned_lgbm_tp[[6]]

shap.plot.summary(shap_values_tp_gaged) 


#############################################################################

#### Chloride

#### Calculate summary stats 

summary_stats_final_model_chlor_gaged <- tuned_lgbm_chlor[[1]]

#### Final time series

pred_obs_ts_chlor_gaged <- tuned_lgbm_chlor[[3]] 

### And SHAP values and plots

shap_values_chlor_gaged <- tuned_lgbm_chlor[[6]]

shap.plot.summary(shap_values_chlor_gaged) 
```

### Ungaged

#### Run the models

``` r
#### For Total Phosphorus

tp_loo <- watershed_loo_runner(tp_train_valid,
                               tp_test, 
                               picked_vars_tp_ungaged,
                               best_params_tp_ungaged,
                               do_save = TRUE,
                               is_tuned = TRUE,
                               constit = "Phosphorus_Total")

#### For Chloride 

chlor_loo <- watershed_loo_runner(chlor_train_valid,
                               chlor_test,
                               picked_vars_chlor_ungaged,
                               best_params_chlor_ungaged,
                               do_save = TRUE,
                               is_tuned = TRUE,
                               constit = "Chloride")
```

#### Extract relevant outputs

``` r
################################################################################

#### For Total Phosphorus

##### First, summary statistics for model predictive performance

sum_stats_tp_ungaged <- tp_loo[[1]]

##### A predicted and observed time series

predicted_vs_observed_tp_ungaged <- tp_loo[[3]]

##### SHAP values for predictions
###### Note that these are for the *test* data

shap_values_tp_ungaged <- tp_loo[[4]]

raw_shap_values_tp_ungaged <- tp_loo[[7]]

##### Split values and other information about the trees
##### Which can help to understand what the model is doing

tree_tables_tp_ungaged <- tp_loo[[8]]

##### Smearing factor for each tributary

d_fact_tp_ungaged <- tp_loo[[5]]

################################################################################

#### For Chloride

##### First, summary statistics for model predictive performance

sum_stats_chlor_ungaged <- chlor_loo[[1]]


##### A predicted and observed time series

predicted_vs_observed_chlor_ungaged <- chlor_loo[[3]]

##### SHAP values for predictions
###### Note that these are for the *test* data

shap_values_chlor_ungaged <- chlor_loo[[4]]

raw_shap_values_chlor_ungaged <- chlor_loo[[7]]

##### Split values and other information about the trees
##### Which can help to understand what the model is doing

tree_tables_chlor_ungaged <- chlor_loo[[8]]

##### Smearing factor for each tributary

d_fact_chlor_ungaged <- chlor_loo[[5]]
```
