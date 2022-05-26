# GRACE-downscaling
Codes needed to implement the model-based downscaling (with depth adjustment) of GRACE ocean bottom pressure near steep bathymetry

This repository consists of a suite of scripts and functions that can be run in MATLAB (or Octave, with possible minor modifications). These scripts use output from the GRACE JPL mascon product (available from https://grace.jpl.nasa.gov/data/get-data/jpl_global_mascons/) and ECCO2 1/4-degree interpolated model output (available from https://ecco.jpl.nasa.gov/drive/files/ECCO2/cube92_latlon_quart_90S90N/). The in-situ data used in validation are available as described in the reference below or you can contact me for them. If you have any questions or need assistance running the code, please do not hesitate to contact me at andrewdelman@ucla.edu.

The reference for these scripts/functions is as follows (please cite if using or adapting this code):

Delman, A. S., and F. Landerer (2022), Process-specific contributions to anomalous Java mixed layer cooling during positive IOD events, Remote Sensing, 14, 1764. https://doi.org/10.3390/rs14071764.

The following codes are included (with a listing of figures in the above paper that each one produced).  In addition, most scripts require the use of MATLAB functions contained in the base_functions.tar file.  Some of these functions also call a file, topo2.grd, which is the 2-minute bathymetry ETOPO2v2g available from https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/.

ECCO_OBP_corr_maps.m: Produces spatial correlation maps (zero-lag and optimum lag) between ocean bottom pressure (OBP) at a given point and OBP in the surrounding region, based on the ECCO2 model output.  Used to produce Fig. 1a-b.

ECCO_OBP_mascon_corr_maps.m: Produces spatial correlation maps (zero-lag and optimum lag) between OBP at a given point and OBP in synthetic mascons in the ECCO2 model output.  Used to produce Fig. 1c-d.

mascon_ECCO_downscale_demo.m: Plots the GRACE mascon OBP and the downscaled OBP (based on ECCO2 spatial correlations/covariances) in a given region and month, along with maps of the number of mascons used in downscaling for each point in the region, and the correlation of a downscaling based on synthetic mascons with synthetic OBP in the model (a sort of theoretical "upper bound" for the accuracy of the downscaling).  Used to produce Fig. 2 and 3.

GRACE_ECCO_downscale_insitu_compare.m: This code generates downscaled OBP at points where in-situ OBP data is also present, and generates statistics comparing/validating the time series.  A global map of downscaled/in-situ correlations is generated, and a .mat file is that is then used in other analysis codes.  Used to produce Fig. 4a.

GRACE_downscale_method_compare.m: This code is used to compare downscaled and non-downscaled GRACE correlations (CDF distributions and global map), and can also be used to compare different downscaling methods.  Uses .mat files produced by GRACE_ECCO_downscale_insitu_compare.m as input.  Used to produce Fig. 4b.

downscale_improv_scatter_plots.m: This code also compares the correlations of downscaled vs. non-downscaled time series (or 2 different downscaling methods), but produces scatter plots comparing the correlations independent of location, and color codes the circles based on some parameter at that location.  Used to produce Fig. 5.

obs_downscale_tseries_compare.m: This code compares the time series of GRACE mascon (non-downscaled), GRACE downscaled, and in-situ OBP at a given point where in-situ data is available.  Used to produce Fig. 6.


