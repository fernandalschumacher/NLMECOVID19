# NLMECOVID19

This repository contains the R codes used to fit a nonlinear mixed-effects model to COVID-19 deaths data, which is presented in the following work:
Schumacher, F. L., Ferreira, C. S., Prates, M. O., Lachos, A., Lachos, V. H. (2020) A robust nonlinear mixed-effects model for COVID-19 deaths data. Submitted. Preprint available at [arXiv:2007.00848](https://arxiv.org/abs/2007.00848).

If you would like to run the application, please follow instructions in file **"covid-application.R"**, and feel free to change the selected countries, but please be aware that nonlinear mixed models are quite sensible to initial values, and they might need adjustments.

The R package skewlmm (available at CRAN and Github) is the linear version of the model, and its documentation might be helpful.
