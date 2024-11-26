---
title: "README"
author: "Bridger Huhn"
date: "2024-11-26"
output: html_document
---


 This folder has three example R scripts from a recent publication I am finalizing. No data is provided.

# My CV labeled: Bridger Huhn CV

 My github respositories are not publically available, but I did create one with WEST when I interned there. If I remember correctly it is under a name similar to Jackson_Lake_Levels. Simon Weller (at WEST) would know where to find it. 



# Sample_bayesian_Parameterization.R
 - Uses photosynthetic light response data from a range of common and endemic plants throughout Wyoming (Data Not Included)
 - Parameterizes these curves using nlsML(). Each plant was fit using this least squared regresssion because we are interested in the parameters of these curves which have a biophysical inportance. (meaning we didn't want to fit the whole population level parameters because we are treating these parameters as a trait of the individual)
 - The parameter outputs of these curves are fit using a bayesian regression with a variety of distributions for each parameter. (we competed vairous distributional fits using leave one out cross validation and found these to be the most parsimonious models)
- I chose iterate over each site in our study and not combine them into one large model for several reasons.

1.  Finding the best (most parsimonious) distribution to model each parameter was no the same for every site.

2. The model had a very difficult time converging when extremely different sites were being used even when each site was run with a fixed mean and variance parameter for each distribution.
3.  Separating out each site made it easier to organize results, plots, and posterior draws.

# Draws_Summaries.R
 - Contains the code to extract posterior draws and summarize model output. It also merges data and predicted model output into one dataframe for easy plotting.

# PLOTS.R
 - This is the code for plotting all of the data and model outputs for a recent publication we are finishing. An example of one of the graphs is incuded as SAMPLEGRAPH.png. It shows the posterior distributions of bayesian parameter estimates (alpha and beta) along with data of light use efficiency curves. In simpler terms, it tells us how different plants convert light to energy at a range of light levels.
