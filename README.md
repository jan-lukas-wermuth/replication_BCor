
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication Material for the Paper “Measuring Dependence between Events”

<!-- badges: start -->
<!-- badges: end -->

The code in this repository generates the results of the simulations and
the empirical applications of the working paper Pohle, M.-O.,
Dimitriadis, T. and Wermuth, J.-L. (2024), Measuring Dependence between
Events available on [arXiv](https://arxiv.org/abs/2403.17580).

## Data Availability

We unfortunately do not have the licence to make the data publicly
available. Therefore, we alienated the data in terms of p, q and r and
perform our analysis with this dataset. The simulations of course work
without any additional files.

## Simulations

The files for the simulations are available in the folder ‘simulations’.
Figure 3 is generated in the file ‘sim_AsyDist.R’.
‘Asymptotic_Skewness_Q.R’ generates Figure 6. Figures 7-9 are a product
of multiple files. The data is simulated in ‘sim_par_BinEvents.R’,
Figure 7 and Figure 8 are produced in ‘sim_eval_par_BinEvents.R’ and
Figure 9 is generated in ‘sim_eval_par_BinEvents_Cdetails.R’.

## Applications

The files for the applications are available in the folder
‘applications’. All the estimators from the introductory example can be
found in the file ‘YuleSmallpox_CIs.R’. The code for Figure 1, Figure 2
and Figure 11 as well as Figure 12 in the appendix is in the file
‘correlations_interdependence.Rmd’. Finally, all the graphics for the
case study on drug use are generated in ‘drugs_interdependence_pqr.Rmd’.
