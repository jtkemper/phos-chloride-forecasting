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
    an **ungaged** scenario, where data from all watershed save one is
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

**Outputs**

################################################################################ 
