# Statistical analysis for "Time Thresholds for using Pressure Reactivity Index in Neuroprognostication for Patients with Severe Traumatic Brain Injury"

This repository contains all R code necessary to reproduce the statistical analysis for the manuscript "Time Thresholds for using Pressure Reactivity Index in Neuroprognostication for Patients with Severe Traumatic Brain Injury".

* *analysis.qmd:* main document creating all figures and result tables.
* *feature_engineering.R:* utility functions to engineer features from the wavelet-filtered trajectories.
* *modeling_utils.R:* utility functions to fit and assess the penalized elastic-net logistic regression models.
* *plotting.R:* utility functions for plotting the results.
* *read_data.R:* R code the read the data, which is not publicly available.
* *wavelet_filtering.R:* utility functions for performing the wavelet filtering.
