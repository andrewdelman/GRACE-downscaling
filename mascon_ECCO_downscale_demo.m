% demonstrate downscaled GRACE along with mascon data for a given month

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
cd('/indopac/adelman/GRACE/')


lon_bounds = [-83.0 -67.0];
lat_bounds = [20.0 32.0];
time_to_plot = [2010 03];

% temporal filtering parameters
steepness_factor = 5;
low_freq_bound = 1/(1*365.24);
high_freq_bound = 1/((3/12)*365.24);
season_cyc_opt = 0;    % 0 = remove seasonal/annual cycle, 1 = retain seasonal/annual cycle

radius_mascons_deg = 20;    % radius of potential mascons to include, in degrees latitude
hybrid_factor = 0;    % ranging from 0 = just use ECCO cross-correlations, to 1 = apply full adjustment to "GRACE range" cross-correlations
hybrid_factor_stddev = 0;   % ranging from 0 = just use ECCO mascon standard deviations to 1 = just use GRACE mascon standard deviations
n_mascons_max = 10;     % maximum number of mascons to use in objective reconstruction
min_corr_to_include = 0.3;    % minimum correlation coefficient (with given point) to include in objective reconstruction
adjust_corr_max = 0.05;    % maximum size of correlation adjustment

% *Important: when changing downscale_id, also change commenting in lines 1010-1150 as needed*

% % downscale_id = ['cs510_cylindweight1.5_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
% downscale_id = 'coloc';    % identifier of method used to define covariances for downscaling
% downscale_id = ['cs510_cylindweight1.5_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
% % downscale_id = ['cs510_cylindweight1.5_slopeperp',num2str(adjust_corr_max),'_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
downscale_id = ['cs510_cylindweight1.5_depthadj',num2str(adjust_corr_max),'_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
% % downscale_id = ['cs510_cylindweight1.5_foverHadj',num2str(adjust_corr_max),'_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
downscale_opt = 1;     % 0 = no downscaling (just correlate with co-located mascon); 1 = downscaling
downscale_adj_opt = 1;   % 0 = no adjustment to spatial correlation/covariance; 1 = adjust spatial correlation

% normalized error tolerance
norm_err_tolerance = 0.3;



% mascon_lat_separation = 3;    % in degrees


curr_file = 'LAND_MASK.CRIv01.nc';
land_mask = ncread(curr_file,'land_mask');
curr_file = 'GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';
lon = ncread(curr_file,'lon');
lat = ncread(curr_file,'lat');
time = ncread(curr_file,'time') + datenum([2002 1 1 0 0 0]);
in_time_range_ind = (1:1:length(time))';
time_range_ind_span = max(in_time_range_ind) - min(in_time_range_ind) + 1;
time = time(in_time_range_ind);
lwe_thickness = ncread(curr_file,'lwe_thickness',[1 1 min(in_time_range_ind)],[length(lon) length(lat) time_range_ind_span]);
lwe_uncert = ncread(curr_file,'uncertainty',[1 1 min(in_time_range_ind)],[length(lon) length(lat) time_range_ind_span]);


time_plot_ind = find(abs(time - datenum([time_to_plot 16 0 0 0])) < 15);



% remove global ocean mean

nan_mask = ones(size(lwe_thickness));
nan_mask(isnan(lwe_thickness) == 1) = 0;
sum_nan_mask = sum(nan_mask,3);
nan_mask(repmat(sum_nan_mask,[1 1 size(nan_mask,3)]) < 0.8*max(max((1 - land_mask).*sum_nan_mask))) = 0;

lwe_thickness(isnan(lwe_thickness) == 1) = 0;
lwe_uncert(isnan(lwe_uncert) == 1) = 0;

lwe_thickness_ocean_mean = sum(sum(nan_mask.*repmat(1 - land_mask,[1 1 size(lwe_thickness,3)]).*lwe_thickness,2),1)./(sum(sum(nan_mask.*repmat(1 - land_mask,[1 1 size(lwe_thickness,3)]),2),1));
lwe_uncert_ocean_mean = (sum(sum(nan_mask.*repmat(1 - land_mask,[1 1 size(lwe_thickness,3)]).*(lwe_uncert.^2),2),1).^(1/2))./(sum(sum(nan_mask.*repmat(1 - land_mask,[1 1 size(lwe_thickness,3)]),2),1));

lwe_thickness_nomean = lwe_thickness - repmat(lwe_thickness_ocean_mean,[size(lwe_thickness,1) size(lwe_thickness,2) 1]);
lwe_thickness_nomean(abs(lwe_thickness) < 1e-10) = NaN;
lwe_uncert(abs(lwe_uncert) < 1e-10) = NaN;
lwe_thickness = lwe_thickness_nomean;
clear lwe_thickness_nomean


% load bathymetry (2-minute, Smith and Sandwell)

ncid_topo = netcdf.open('topo2.grd','NC_NOWRITE');
z = netcdf.getVar(ncid_topo,5);

spacing = 1/30;     %spacing of bathymetry data, in degrees

lat_bathy = 90 - ((0.5*spacing):spacing:180)';
lon_bathy = ((0.5*spacing):spacing:360)';

z_reshaped = reshape(z,360/spacing,180/spacing);
clear z

z_reshaped = flip(z_reshaped,2);
lat_bathy = flip(lat_bathy,1);


% identify ocean mascons

[lat_grid,lon_grid] = meshgrid(lat,lon);
curr_lwe_thickness = lwe_thickness(:,:,1);
ocean_mascon_curr_values_all = unique(curr_lwe_thickness(abs(land_mask) < 1e-5));

mascon_lon_center_all = NaN([length(ocean_mascon_curr_values_all) 1]);
mascon_lat_center_all = NaN([length(ocean_mascon_curr_values_all) 1]);
mascon_lon_bounds_all = NaN([length(ocean_mascon_curr_values_all) 2]);
mascon_lat_bounds_all = NaN([length(ocean_mascon_curr_values_all) 2]);
lwe_thickness_ocean_mascons_all = NaN([length(ocean_mascon_curr_values_all) length(time)]);
lwe_uncert_ocean_mascons_all = NaN([length(ocean_mascon_curr_values_all) length(time)]);
mascon_depth_avg_all = NaN([length(ocean_mascon_curr_values_all) 1]);
for mascon_ind = 1:length(mascon_lon_center_all)
    curr_mascon_value = ocean_mascon_curr_values_all(mascon_ind);
    
    in_mascon_ind = find(abs(curr_lwe_thickness - curr_mascon_value) < 1e-10);
    lon_grid_spacing = mean(mean(diff(lon_grid,1,1)));
    lat_grid_spacing = mean(mean(diff(lat_grid,1,2)));
    if max(diff(sort(lon_grid(in_mascon_ind),'ascend'))) > 2*lon_grid_spacing
        lon_in_mascon_adj = mod(lon_grid(in_mascon_ind) + 180,360) - 180;
        mascon_lon_center_all(mascon_ind) = mod(mean(lon_in_mascon_adj),360);
        mascon_lon_bounds_all(mascon_ind,:) = mod(min(lon_in_mascon_adj),360) + [(-0.5*lon_grid_spacing) (mod(max(lon_in_mascon_adj) - min(lon_in_mascon_adj),360) + (0.5*lon_grid_spacing))];
        in_mascon_lon_bathy_ind = find((mod(lon_bathy - (min(lon_in_mascon_adj) - (0.5*lon_grid_spacing)) + 180,360) - 180 > 0) & (mod(lon_bathy - (max(lon_in_mascon_adj) + (0.5*lon_grid_spacing)) + 180,360) - 180 < 0));
    else
        mascon_lon_center_all(mascon_ind) = mean(lon_grid(in_mascon_ind));
        mascon_lon_bounds_all(mascon_ind,:) = [(min(lon_grid(in_mascon_ind)) - (0.5*lon_grid_spacing)) (max(lon_grid(in_mascon_ind)) + (0.5*lon_grid_spacing))];
        in_mascon_lon_bathy_ind = find((mod(lon_bathy - (min(lon_grid(in_mascon_ind)) - (0.5*lon_grid_spacing)) + 180,360) - 180 > 0) & (mod(lon_bathy - (max(lon_grid(in_mascon_ind)) + (0.5*lon_grid_spacing)) + 180,360) - 180 < 0));
    end
    mascon_lat_center_all(mascon_ind) = mean(lat_grid(in_mascon_ind));
    mascon_lat_bounds_all(mascon_ind,:) = [(min(lat_grid(in_mascon_ind)) - (0.5*lat_grid_spacing)) (max(lat_grid(in_mascon_ind)) + (0.5*lat_grid_spacing))];
    in_mascon_lat_bathy_ind = find((lat_bathy - (min(lat_grid(in_mascon_ind)) - (0.5*lat_grid_spacing)) > 0) & (lat_bathy - (max(lat_grid(in_mascon_ind)) + (0.5*lat_grid_spacing)) < 0));
    
    curr_mascon_i = mod(in_mascon_ind(1) - 1,size(curr_lwe_thickness,1)) + 1;
    curr_mascon_j = ceil(in_mascon_ind(1)/size(curr_lwe_thickness,1));
    lwe_thickness_ocean_mascons_all(mascon_ind,:) = reshape(lwe_thickness(curr_mascon_i,curr_mascon_j,:),[1 size(lwe_thickness,3)]);
    lwe_uncert_ocean_mascons_all(mascon_ind,:) = reshape(lwe_uncert(curr_mascon_i,curr_mascon_j,:),[1 size(lwe_uncert,3)]);
    
    in_mascon_depth = -z_reshaped(in_mascon_lon_bathy_ind,in_mascon_lat_bathy_ind);
    mascon_depth_avg_all(mascon_ind) = mean(in_mascon_depth(in_mascon_depth > 0));
    
end


% create grid of points for downscaling

grid_spacing_downsc = 0.1;

lon_downsc = ((lon_bounds(1) - grid_spacing_downsc):grid_spacing_downsc:(lon_bounds(2) + grid_spacing_downsc))';
lat_downsc = ((lat_bounds(1) - grid_spacing_downsc):grid_spacing_downsc:(lat_bounds(2) + grid_spacing_downsc));
[lat_downsc_grid,lon_downsc_grid] = meshgrid(lat_downsc,lon_downsc);


% interpolate bathymetry to downscaling grid

in_local_bathy_lat_range_ind = find((lat_bathy - lat_bounds(1) >= -5) & (lat_bathy - lat_bounds(2) <= 5));
local_lat_bathy = lat_bathy(in_local_bathy_lat_range_ind);
if mod(lon_bounds(2) + 5 - max(lon_bathy),360) < mod(lon_bounds(1) - 5 - min(lon_bathy),360)
    in_local_bathy_lon_range_ind_1 = find(lon_bathy - (mod(lon_bounds(1) - 5 - (min(lon_bathy) - 1e-5),360) + (min(lon_bathy) - 1e-5)) >= 0);
    in_local_bathy_lon_range_ind_2 = find(lon_bathy - (mod(lon_bounds(2) + 5 - (min(lon_bathy) - 1e-5),360) + (min(lon_bathy) - 1e-5)) <= 0);
    local_lon_bathy = mod([lon_bathy(in_local_bathy_lon_range_ind_1); (lon_bathy(in_local_bathy_lon_range_ind_2) + 360)] - (mean(lon_bounds) - 180),360) + (mean(lon_bounds) - 180);
    z_local = z_reshaped([in_local_bathy_lon_range_ind_1; in_local_bathy_lon_range_ind_2],in_local_bathy_lat_range_ind);
else
    in_local_bathy_lon_range_ind = find((mod(lon_bathy - lon_bounds(1) + 180,360) - 180 >= -5) & (mod(lon_bathy - lon_bounds(2) + 180,360) - 180 <= 5));
    local_lon_bathy = mod(lon_bathy(in_local_bathy_lon_range_ind) - (mean(lon_bounds) - 180),360) + (mean(lon_bounds) - 180);
    z_local = z_reshaped(in_local_bathy_lon_range_ind,in_local_bathy_lat_range_ind);
end


depth_downsc_grid = -reshape(interp2_fft(mod(local_lon_bathy - mean(lon_bounds) + 180,360) + mean(lon_bounds) - 180,local_lat_bathy,z_local,reshape(lon_downsc_grid,[numel(lon_downsc_grid) 1]),reshape(lat_downsc_grid,[numel(lat_downsc_grid) 1])),[length(lon_downsc) length(lat_downsc)]);



% find mascons near map region

in_region_ind = find(((cosd(mascon_lat_center_all)).*(mod(mascon_lon_center_all - lon_bounds(1) + 180,360) - 180) >= -radius_mascons_deg - 1) & ((cosd(mascon_lat_center_all)).*(mod(mascon_lon_center_all - lon_bounds(2) + 180,360) - 180) <= radius_mascons_deg + 1) & (mascon_lat_center_all - lat_bounds(1) >= -radius_mascons_deg - 1) & (mascon_lat_center_all - lat_bounds(2) <= radius_mascons_deg + 1));

ocean_mascon_curr_values = ocean_mascon_curr_values_all(in_region_ind);
mascon_lon_center = mascon_lon_center_all(in_region_ind);
mascon_lat_center = mascon_lat_center_all(in_region_ind);
mascon_lon_bounds = mascon_lon_bounds_all(in_region_ind,:);
mascon_lat_bounds = mascon_lat_bounds_all(in_region_ind,:);
lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons_all(in_region_ind,:);
lwe_uncert_ocean_mascons = lwe_uncert_ocean_mascons_all(in_region_ind,:);
mascon_depth_avg = mascon_depth_avg_all(in_region_ind);


% if high_freq_bound < 1/365
if abs(season_cyc_opt) < 1e-5

    % remove annual cycle

    nan_mask = ones(size(lwe_thickness_ocean_mascons));
    nan_mask((isnan(lwe_thickness_ocean_mascons) == 1) | (abs(lwe_thickness_ocean_mascons) < 1e-10)) = 0;
    sum_nan_mask = sum(nan_mask,1);
    nan_mask(repmat(sum_nan_mask,[size(nan_mask,1) 1]) < 0.8*max(sum_nan_mask)) = 0;

    lwe_thickness_ocean_mascons(isnan(lwe_thickness_ocean_mascons) == 1) = 0;

    month_centers = datenum([(2002*ones([12 1])) (1:1:12)' repmat([16 0 0 0],[12 1])]);

    lwe_thickness_ocean_mascons_anom = NaN(size(lwe_thickness_ocean_mascons));
    for month_ind = 1:12
        curr_month_center = month_centers(month_ind);
        curr_t_bin_ind = find(abs(mod(time - curr_month_center + 180,365.24) - 180) < 15);

        lwe_thickness_ocean_mascons_anom(:,curr_t_bin_ind) = lwe_thickness_ocean_mascons(:,curr_t_bin_ind) - repmat(sum(nan_mask(:,curr_t_bin_ind).*lwe_thickness_ocean_mascons(:,curr_t_bin_ind),2)./(sum(nan_mask(:,curr_t_bin_ind),2)),[1 length(curr_t_bin_ind)]);
    end

    lwe_thickness_ocean_mascons_anom(abs(lwe_thickness_ocean_mascons) < 1e-10) = NaN;
    lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons_anom;
    clear lwe_thickness_ocean_mascons_anom
end


% remove linear trend (and bandpass if desired)

half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)
edge_handling_opt = 0;    % mask out edges where filtering is less accurate? (0 = no, 1 = yes)

% interpolate to regular grid so that filtering accounts for gaps in data
time_GRACE_interp = ((((365.24/12)*round((min(time) - datenum([2002 4 16 0 0 0]))/(365.24/12))) + datenum([2002 4 16 0 0 0])):(365.24/12):(((365.24/12)*round((max(time) - datenum([2018 4 16 0 0 0]))/(365.24/12))) + datenum([2018 4 16 0 0 0])))';
nan_mask_tinterp = zeros([size(nan_mask,1) length(time_GRACE_interp)]);
curr_tinterp_ind = 1;
while curr_tinterp_ind <= length(time_GRACE_interp)
    within_range_ind = find((time - time_GRACE_interp(curr_tinterp_ind) >= -((365.24/12)/2) + 1e-5) & (time - time_GRACE_interp(curr_tinterp_ind) < ((365.24/12)/2) + 1e-5));
    if isempty(within_range_ind) == 0
        time_GRACE_interp = [time_GRACE_interp(1:(curr_tinterp_ind - 1)); time(within_range_ind); time_GRACE_interp((curr_tinterp_ind + 1):length(time_GRACE_interp))];
        nan_mask_tinterp = [nan_mask_tinterp(:,1:(curr_tinterp_ind - 1)) nan_mask(:,within_range_ind) nan_mask_tinterp(:,(curr_tinterp_ind + 1):size(nan_mask_tinterp,2))]; 
        curr_tinterp_ind = curr_tinterp_ind + length(within_range_ind);
    else
        curr_tinterp_ind = curr_tinterp_ind + 1;
    end
end
[time_GRACE_interp,unique_ind] = unique(time_GRACE_interp);
nan_mask_tinterp = nan_mask_tinterp(:,unique_ind);
good_ind = find(sum(nan_mask,1) >= 0.8*max(sum(nan_mask,1)));
lwe_thickness_tinterp = (interp1(time(good_ind),lwe_thickness_ocean_mascons(:,good_ind)',time_GRACE_interp))';

[lwe_thickness_ocean_mascons_filtered,lwe_thickness_ocean_mascons_trend,~] = bandpass_err_fcn(lwe_thickness_tinterp,2,mean(diff(time_GRACE_interp)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,edge_handling_opt,1);
filter_gain_coeffs_array = bandpass_err_fcn_gain_coeffs(lwe_thickness_tinterp,2,mean(diff(time_GRACE_interp)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1);
% % lwe_uncert_ocean_mascons_filtered = squeeze((sum((abs(filter_gain_coeffs_array).^2).*(repmat(reshape(lwe_uncert_ocean_mascons,[size(lwe_uncert_ocean_mascons,1) 1 size(lwe_uncert_ocean_mascons,2)]).^2,[1 size(lwe_uncert_ocean_mascons,2) 1])),3)).^(1/2));
% % clear filter_gain_coeffs_array
% [lwe_thickness_nan_mask_filtered,~,~] = bandpass_err_fcn(nan_mask,2,mean(diff(time)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,edge_handling_opt,1);


% estimate uncertainty due to missing obs., and propagate through time filter

sum_nan_mask_tinterp = sum(nan_mask_tinterp,2);
% [~,mascons_dof_zero_lag,~,~,~,~,mascons_std_dev,~,~] = regression_linear_scalar_scalar(lwe_thickness_tinterp,lwe_thickness_tinterp,2,mean(diff(time_GRACE_interp)),mean(diff(time_GRACE_interp)),(1/5)*length(time_GRACE_interp)*mean(diff(time_GRACE_interp))*[0 1],0.95,0);
[~,mascons_dof_zero_lag,~,~,~,~,mascons_std_dev,~,lags] = regression_linear_scalar_scalar(lwe_thickness_ocean_mascons_filtered,lwe_thickness_ocean_mascons_filtered,2,mean(diff(time_GRACE_interp)),mean(diff(time_GRACE_interp)),(1/5)*length(time_GRACE_interp)*mean(diff(time_GRACE_interp))*[0 1],0.95,0);
mascons_std_dev = mascons_std_dev(:,abs(lags) < 1e-5);
mascons_decorr_timescale = (mean(diff(time_GRACE_interp))*sum_nan_mask_tinterp)./mascons_dof_zero_lag;
mascons_obs_separation_array = repmat(abs(repmat(time_GRACE_interp',[1 1 length(time_GRACE_interp)]) - repmat(reshape(time_GRACE_interp,[1 1 length(time_GRACE_interp)]),[1 length(time_GRACE_interp) 1])),[size(nan_mask_tinterp,1) 1 1]);
mascons_missing_obs_err_cov = repmat(mascons_std_dev.^2,[1 size(nan_mask_tinterp,2) size(nan_mask_tinterp,2)]).*(1 - (mascons_obs_separation_array./repmat(mascons_decorr_timescale,[1 size(nan_mask_tinterp,2) size(nan_mask_tinterp,2)])));
mascons_missing_obs_err_cov(mascons_missing_obs_err_cov < 0) = 0;
mascons_missing_obs_err = NaN(size(lwe_thickness_tinterp));
for curr_mascon_ind = 1:size(nan_mask_tinterp,1)
    for t = 1:size(nan_mask_tinterp,2)
        curr_in_range_ind = find(abs(filter_gain_coeffs_array(curr_mascon_ind,t,:)) >= 0.01*max(abs(filter_gain_coeffs_array(curr_mascon_ind,t,:))));
        curr_in_range_ind = (min(curr_in_range_ind):1:max(curr_in_range_ind))';
        curr_mascons_err_cross_mask = (1 - repmat(nan_mask_tinterp(curr_mascon_ind,curr_in_range_ind),[1 1 length(curr_in_range_ind)])).*(1 - repmat(reshape(nan_mask_tinterp(curr_mascon_ind,curr_in_range_ind),[1 1 length(curr_in_range_ind)]),[1 length(curr_in_range_ind) 1]));
        curr_mascons_missing_obs_err_cov = curr_mascons_err_cross_mask.*mascons_missing_obs_err_cov(curr_mascon_ind,curr_in_range_ind,curr_in_range_ind);
        curr_cov_with_filter_coeffs = repmat(reshape(filter_gain_coeffs_array(curr_mascon_ind,t,curr_in_range_ind),[1 length(curr_in_range_ind)]),[1 1 length(curr_in_range_ind)]).*repmat(filter_gain_coeffs_array(curr_mascon_ind,t,curr_in_range_ind),[1 length(curr_in_range_ind) 1]).*curr_mascons_missing_obs_err_cov;
        mascons_missing_obs_err(curr_mascon_ind,t) = (sum(sum(curr_cov_with_filter_coeffs,3),2))^(1/2);
    end
end


lwe_thickness_ocean_mascons = (interp1(time_GRACE_interp,lwe_thickness_ocean_mascons_filtered',time))';
% lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons_filtered;
% clear lwe_thickness_ocean_mascons_filtered
% lwe_uncert_ocean_mascons = lwe_uncert_ocean_mascons_filtered;
% clear lwe_uncert_ocean_mascons_filtered

mascons_missing_obs_norm_err = mascons_missing_obs_err./repmat(mascons_std_dev,[1 size(mascons_missing_obs_err,2)]);
lwe_missing_obs_norm_err = (interp1(time_GRACE_interp,mascons_missing_obs_norm_err',time))';

% lwe_thickness_ocean_mascons(abs(lwe_thickness_nan_mask_filtered) > 0.2) = NaN;
lwe_thickness_ocean_mascons(abs(lwe_missing_obs_norm_err) > norm_err_tolerance) = NaN;




if abs(downscale_opt - 1) < 1e-5

    % load ECCO OBP for covariance calculations

    ECCO_nc_file = '/indopac/adelman/ECCO2/PHIBOT.ECCO2.lonlatinterp.1992-2018.nc';

    longitude = ncread(ECCO_nc_file,'LONGITUDE_T');
    latitude = ncread(ECCO_nc_file,'LATITUDE_T');
    time_ECCO = ncread(ECCO_nc_file,'TIME');
    time_ECCO = double(time_ECCO) + datenum([1992 1 1 0 0 0]);

    in_lon_range_ind = find((cosd(mean(lat_bounds)).*(mod(longitude - lon_bounds(1) + 180,360) - 180) >= -radius_mascons_deg - 5) & (cosd(mean(lat_bounds)).*(mod(longitude - lon_bounds(2) + 180,360) - 180) <= radius_mascons_deg + 5));
    in_lat_range_ind = find((latitude - lat_bounds(1) >= -radius_mascons_deg - 5) & (latitude - lat_bounds(2) <= radius_mascons_deg + 5));
    % in_time_range_ind = find((time >= datenum(time_range_start)) & (time < datenum(time_range_end)));
    in_time_range_ind = (1:1:length(time_ECCO))';

    if length(find(ismember([1 length(longitude)],in_lon_range_ind) == 1)) > 1
        gap_ind = find(diff(in_lon_range_ind) > 1.5);

        lon_in_range_ECCO = longitude(in_lon_range_ind([((gap_ind + 1):1:length(in_lon_range_ind)) (1:1:gap_ind)]));
        lat_in_range_ECCO = latitude(in_lat_range_ind);
        time_in_range_ECCO = time_ECCO(in_time_range_ind);

        start_vec = [in_lon_range_ind(gap_ind + 1) min(in_lat_range_ind) min(in_time_range_ind)];
        count_vec = [(max(in_lon_range_ind) - in_lon_range_ind(gap_ind + 1) + 1) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];

        OBP_ECCO = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
        OBP_ECCO = 1027.5*OBP_ECCO(in_lon_range_ind((gap_ind + 1):length(in_lon_range_ind)) - in_lon_range_ind(gap_ind + 1) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1);

        start_vec = [1 min(in_lat_range_ind) min(in_time_range_ind)];
        count_vec = [in_lon_range_ind(gap_ind) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];

        OBP_ECCO_2 = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
        OBP_ECCO = [OBP_ECCO; (1027.5*OBP_ECCO_2(in_lon_range_ind(1:gap_ind) - min(in_lon_range_ind) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1))];
    else

        lon_in_range_ECCO = longitude(in_lon_range_ind);
        lat_in_range_ECCO = latitude(in_lat_range_ind);
        time_in_range_ECCO = time_ECCO(in_time_range_ind);

        start_vec = [min(in_lon_range_ind) min(in_lat_range_ind) min(in_time_range_ind)];
        count_vec = [(max(in_lon_range_ind) - min(in_lon_range_ind) + 1) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];

        OBP_ECCO = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
        OBP_ECCO = 1027.5*OBP_ECCO(in_lon_range_ind - min(in_lon_range_ind) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1);
    end

    OBP_ECCO(OBP_ECCO < -1e26) = NaN;

    diff_lon_in_range_ECCO = diff(lon_in_range_ECCO);
    diff_lon_in_range_ECCO = mod(diff_lon_in_range_ECCO + 180,360) - 180;
    lon_in_range_ECCO = lon_in_range_ECCO(1) + (360*floor((mean(lon_bounds) - lon_in_range_ECCO(1))/360)) + [0; cumsum(diff_lon_in_range_ECCO)];
    diff_lat_in_range_ECCO = diff(lat_in_range_ECCO);

    lon_in_range_ECCO_bounds = [(lon_in_range_ECCO(1) - (0.5*diff_lon_in_range_ECCO(1))); (lon_in_range_ECCO(2:length(lon_in_range_ECCO)) - (diff_lon_in_range_ECCO/2)); (lon_in_range_ECCO(length(lon_in_range_ECCO)) + (0.5*diff_lon_in_range_ECCO(length(diff_lon_in_range_ECCO))))];
    lat_in_range_ECCO_bounds = [(lat_in_range_ECCO(1) - (0.5*diff_lat_in_range_ECCO(1))); (lat_in_range_ECCO(2:length(lat_in_range_ECCO)) - (diff_lat_in_range_ECCO/2)); (lat_in_range_ECCO(length(lat_in_range_ECCO)) + (0.5*diff_lat_in_range_ECCO(length(diff_lat_in_range_ECCO))))];

    
    
    size_array = size(OBP_ECCO);

    time_datevec_ECCO = datevec(time_in_range_ECCO);
    n_months_ECCO = ((12*time_datevec_ECCO(size(time_datevec_ECCO,1),1)) + time_datevec_ECCO(size(time_datevec_ECCO,1),2)) - ((12*time_datevec_ECCO(1,1)) + time_datevec_ECCO(1,2)) + 1;
    nan_mask_ECCO_reshaped = ones([prod(size_array(1:2)) size_array(3)]);
    nan_mask_ECCO_reshaped((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-10)) = 0;
    OBP_ECCO_reshaped_nans_zeroed = reshape(OBP_ECCO,[prod(size_array(1:2)) size_array(3)]);
    OBP_ECCO_reshaped_nans_zeroed(abs(nan_mask_ECCO_reshaped) < 1e-5) = 0;

    time_ECCO_monthavg = NaN([n_months_ECCO 1]);
    OBP_ECCO_monthavg = NaN([prod(size_array(1:2)) n_months_ECCO]);
    for curr_month_ind = 1:n_months_ECCO
        curr_month_num = (12*time_datevec_ECCO(1,1)) + time_datevec_ECCO(1,2) + curr_month_ind - 1;
        curr_yearmonth = [floor((curr_month_num - 1)/12) (mod(curr_month_num - 1,12) + 1)];
        in_month_ind = find((abs(time_datevec_ECCO(:,1) - curr_yearmonth(1)) < 1e-3) & (abs(time_datevec_ECCO(:,2) - curr_yearmonth(2)) < 1e-3));

        enough_ind = find(sum(nan_mask_ECCO_reshaped(:,in_month_ind),2) > (0.8*28)/mean(diff(time_in_range_ECCO)));
        time_ECCO_monthavg(curr_month_ind) = sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*repmat(reshape(time_in_range_ECCO(in_month_ind),[1 length(in_month_ind)]),[length(enough_ind) 1]),2),1)./(sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2),1));
        OBP_ECCO_monthavg(enough_ind,curr_month_ind) = sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*OBP_ECCO_reshaped_nans_zeroed(enough_ind,in_month_ind),2)./(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2));    
    end

    time_in_range_ECCO = time_ECCO_monthavg;
    OBP_ECCO = reshape(OBP_ECCO_monthavg,[size_array(1:2) n_months_ECCO]);
    
    
    size_array = size(OBP_ECCO);



    % define locations of GRACE mascons

    mascon_lat_separation = 3;    % in degrees

    curr_file = 'LAND_MASK.CRIv01.nc';
    land_mask = ncread(curr_file,'land_mask');
    curr_file = 'GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';

    [lat_ECCO_grid,lon_ECCO_grid] = meshgrid(lat_in_range_ECCO,lon_in_range_ECCO);

    % mascon_lon_center = NaN([length(ocean_mascon_curr_values) 1]);
    % mascon_lat_center = NaN([length(ocean_mascon_curr_values) 1]);
    % mascon_lon_bounds = NaN([length(ocean_mascon_curr_values) 2]);
    % mascon_lat_bounds = NaN([length(ocean_mascon_curr_values) 2]);
    OBP_ECCO_ocean_mascons = NaN([length(ocean_mascon_curr_values) size(OBP_ECCO,3)]);
    % mascon_depth_avg = NaN([length(ocean_mascon_curr_values) 1]);
    for mascon_ind = 1:length(mascon_lon_center)
        curr_mascon_value = ocean_mascon_curr_values(mascon_ind);

        in_mascon_ind = find(abs(curr_lwe_thickness - curr_mascon_value) < 1e-10);

    %     lon_GRACE_grid_spacing = mean(mean(diff(lon_grid,1,1)));
    %     lat_GRACE_grid_spacing = mean(mean(diff(lat_grid,1,2)));
    %     if max(diff(sort(lon_grid(in_mascon_ind),'ascend'))) > 2*lon_GRACE_grid_spacing
    %         lon_in_mascon_adj = mod(lon_grid(in_mascon_ind) + 180,360) - 180;
    %         mascon_lon_center(mascon_ind) = mod(mean(lon_in_mascon_adj),360);
    %         mascon_lon_bounds(mascon_ind,:) = mod(min(lon_in_mascon_adj),360) + [(-0.5*lon_GRACE_grid_spacing) (mod(max(lon_in_mascon_adj) - min(lon_in_mascon_adj),360) + (0.5*lon_GRACE_grid_spacing))];
    %     else
    %         mascon_lon_center(mascon_ind) = mean(lon_grid(in_mascon_ind));
    %         mascon_lon_bounds(mascon_ind,:) = [(min(lon_grid(in_mascon_ind)) - (0.5*lon_GRACE_grid_spacing)) (max(lon_grid(in_mascon_ind)) + (0.5*lon_GRACE_grid_spacing))];
    %         
    %     end
    % %     in_mascon_lon_bathy_ind = find((mod(lon_bathy - mascon_lon_bounds(mascon_ind,1) + 180,360) - 180 > 0) & (mod(lon_bathy - mascon_lon_bounds(mascon_ind,2) + 180,360) - 180 < 0));
    %     mascon_lat_center(mascon_ind) = mean(lat_grid(in_mascon_ind));
    %     mascon_lat_bounds(mascon_ind,:) = [(min(lat_grid(in_mascon_ind)) - (0.5*lat_GRACE_grid_spacing)) (max(lat_grid(in_mascon_ind)) + (0.5*lat_GRACE_grid_spacing))];
    % %     in_mascon_lat_bathy_ind = find((lat_bathy - mascon_lat_bounds(mascon_ind,1) > 0) & (lat_bathy - mascon_lat_bounds(mascon_ind,2) < 0));

        if (mod(mascon_lon_bounds(mascon_ind,1) - min(lon_in_range_ECCO_bounds) + 180,360) - 180 > 0) && (mod(mascon_lon_bounds(mascon_ind,2) - max(lon_in_range_ECCO_bounds) + 180,360) - 180 < 0) && (mascon_lat_bounds(mascon_ind,1) - min(lat_in_range_ECCO_bounds) > 0) && (mascon_lat_bounds(mascon_ind,2) - max(lat_in_range_ECCO_bounds) < 0)
            % round region
            radius_cone = 1.5;   % in degrees latitude
            dist_from_center = abs((111100*cosd(mascon_lat_center(mascon_ind)).*(mod(lon_ECCO_grid - mascon_lon_center(mascon_ind) + 180,360) - 180)) + (1i*111100*(lat_ECCO_grid - mascon_lat_center(mascon_ind))));
            weight_matrix = 1 - (1*dist_from_center/(111100*radius_cone));
            weight_matrix(weight_matrix < 0) = 0;

            % for cylindrical (not conical) weight
            weight_matrix(weight_matrix > 1e-5) = 1;

            in_weight_range_ind = find(weight_matrix > 1e-5);
            in_mascon_lon_ECCO_ind = unique(mod(in_weight_range_ind - 1,size(lon_in_range_ECCO,1)) + 1);
            in_mascon_lat_ECCO_ind = unique(ceil(in_weight_range_ind/size(lon_in_range_ECCO,1)));
            weight_matrix = weight_matrix(in_mascon_lon_ECCO_ind,in_mascon_lat_ECCO_ind);
            % % 

    %         % for box region
    %         in_mascon_lon_ECCO_ind = find((mod(lon_in_range_ECCO - mascon_lon_bounds(mascon_ind,1) + 180,360) - 180 >= 0) & (mod(lon_in_range_ECCO - mascon_lon_bounds(mascon_ind,2) + 180,360) - 180 < 0));
    %         in_mascon_lat_ECCO_ind = find((lat_in_range_ECCO - mascon_lat_bounds(mascon_ind,1) >= 0) & (lat_in_range_ECCO - mascon_lat_bounds(mascon_ind,2) < 0));
    %         weight_matrix = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind)]);
    %         % % 

            if max(diff(lon_in_range_ECCO(in_mascon_lon_ECCO_ind))) > 3*median(diff(lon_in_range_ECCO))
                gap_ind = find(diff(lon_in_range_ECCO(in_mascon_lon_ECCO_ind)) > 3*median(diff(lon_in_range_ECCO)));
                in_mascon_lon_width_ECCO = [diff(lon_in_range_ECCO_bounds(in_mascon_lon_ECCO_ind((gap_ind + 1):length(in_mascon_lon_ECCO_ind)))); diff(lon_in_range_ECCO_bounds(in_mascon_lon_ECCO_ind([(1:1:gap_ind)'; (mod(gap_ind + [-1; 0] - 1,length(in_mascon_lon_ECCO_ind)) + 1)])))];
                in_mascon_lon_ECCO_ind = in_mascon_lon_ECCO_ind([((gap_ind + 1):1:length(in_mascon_lon_ECCO_ind))'; (1:1:gap_ind)']);
            else
                in_mascon_lon_width_ECCO = diff(lon_in_range_ECCO_bounds([in_mascon_lon_ECCO_ind; (max(in_mascon_lon_ECCO_ind) + 1)]));
            end
            in_mascon_lat_width_ECCO = diff(lat_in_range_ECCO_bounds([in_mascon_lat_ECCO_ind; (max(in_mascon_lat_ECCO_ind) + 1)]));

            nan_mask_ECCO_in_mascon = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind) size(OBP_ECCO,3)]);
            OBP_ECCO_in_mascon_zeronans = OBP_ECCO(in_mascon_lon_ECCO_ind,in_mascon_lat_ECCO_ind,:);
            nan_mask_ECCO_in_mascon((isnan(OBP_ECCO_in_mascon_zeronans) == 1) | (abs(OBP_ECCO_in_mascon_zeronans) < 1e-15)) = 0;
            nan_mask_ECCO_not_enough = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind)]);
            nan_mask_ECCO_not_enough(sum(nan_mask_ECCO_in_mascon,3) < 0.9*max(max(sum(nan_mask_ECCO_in_mascon,3)))) = 0;
            nan_mask_ECCO_in_mascon(abs(repmat(nan_mask_ECCO_not_enough,[1 1 size(nan_mask_ECCO_in_mascon,3)])) < 1e-5) = 0;
            OBP_ECCO_in_mascon_zeronans(abs(nan_mask_ECCO_in_mascon) < 1e-5) = 0;

            weight_matrix = weight_matrix.*((111100*repmat(cosd(lat_in_range_ECCO(in_mascon_lat_ECCO_ind)'),[length(in_mascon_lon_ECCO_ind) 1]).*repmat(in_mascon_lon_width_ECCO,[1 length(in_mascon_lat_ECCO_ind)])).*(111100*repmat(in_mascon_lat_width_ECCO',[length(in_mascon_lon_ECCO_ind) 1])));
            OBP_ECCO_ocean_mascons(mascon_ind,:) = reshape(sum(sum(repmat(weight_matrix,[1 1 size(OBP_ECCO,3)]).*nan_mask_ECCO_in_mascon.*OBP_ECCO_in_mascon_zeronans,2),1)./(sum(sum(repmat(weight_matrix,[1 1 size(OBP_ECCO,3)]).*nan_mask_ECCO_in_mascon,2),1)),[1 size(OBP_ECCO,3)]);
        end

    %     in_mascon_depth = -z_reshaped(in_mascon_lon_bathy_ind,in_mascon_lat_bathy_ind);
    %     mascon_depth_avg(mascon_ind) = mean(in_mascon_depth(in_mascon_depth > 0));

    end

    nan_mask_ocean_mascons = ones(size(OBP_ECCO_ocean_mascons));
    nan_mask_ocean_mascons((isnan(OBP_ECCO_ocean_mascons) == 1) | (abs(OBP_ECCO_ocean_mascons) < 1e-15)) = 0;
    good_mascon_ind = find(sum(nan_mask_ocean_mascons,2) > 0.9*max(sum(nan_mask_ocean_mascons,2)));
    mascon_lon_center = mascon_lon_center(good_mascon_ind);
    mascon_lat_center = mascon_lat_center(good_mascon_ind);
    mascon_lon_bounds = mascon_lon_bounds(good_mascon_ind,:);
    mascon_lat_bounds = mascon_lat_bounds(good_mascon_ind,:);
    lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons(good_mascon_ind,:);
    lwe_uncert_ocean_mascons = lwe_uncert_ocean_mascons(good_mascon_ind,:);
    OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons(good_mascon_ind,:);
    mascon_depth_avg = mascon_depth_avg(good_mascon_ind);

    % mascon_lon_center = mod(mascon_lon_center - (lon_pt_obs - 180),360) + (lon_pt_obs - 180);
    % mascon_lon_bounds = mod(mascon_lon_bounds - (lon_pt_obs - 180),360) + (lon_pt_obs - 180);


    % temporally filter time series

    steepness_factor = 5;
    half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)


    % ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO));
    % ECCO_nan_mask((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-10)) = -1;
    % [OBP_ECCO_filtered,OBP_ECCO_trend,~] = bandpass_err_fcn(OBP_ECCO,3,mean(diff(time_in_range_ECCO)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    % [OBP_ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,3,mean(diff(time_in_range_ECCO)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    % 
    % OBP_ECCO = OBP_ECCO_filtered;
    % clear OBP_ECCO_filtered


    ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO_ocean_mascons));
    ECCO_nan_mask((isnan(OBP_ECCO_ocean_mascons) == 1) | (abs(OBP_ECCO_ocean_mascons) < 1e-10)) = -1;
    [OBP_ECCO_ocean_mascons_filtered,OBP_ECCO_ocean_mascons_trend,~] = bandpass_err_fcn(OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range_ECCO)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    [ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,2,mean(diff(time_in_range_ECCO)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

    OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons_filtered;
    clear OBP_ECCO_ocean_mascons_filtered
    
    % nan_OBP_ECCO_ind = find(isnan(OBP_ECCO) == 1);

    nan_ECCO_ind = find(isnan(OBP_ECCO_ocean_mascons) == 1);
    if abs(season_cyc_opt) < 1e-5
        % remove annual cycle


    %     OBP_ECCO(nan_OBP_ECCO_ind) = 0;
    %     G = [cos(((2*pi)/365.2425)*time_in_range_ECCO) sin(((2*pi)/365.2425)*time_in_range_ECCO) cos(((2*(2*pi))/365.2425)*time_in_range_ECCO) sin(((2*(2*pi))/365.2425)*time_in_range_ECCO) cos(((3*(2*pi))/365.2425)*time_in_range_ECCO) sin(((3*(2*pi))/365.2425)*time_in_range_ECCO) cos(((4*(2*pi))/365.2425)*time_in_range_ECCO) sin(((4*(2*pi))/365.2425)*time_in_range_ECCO)];
    %     coeffs = (((G')*G)^(-1))*((G')*(permute(reshape(OBP_ECCO,[(size(OBP_ECCO,1)*size(OBP_ECCO,2)) size(OBP_ECCO,3)]),[2 1])));
    %     OBP_ECCO = OBP_ECCO - reshape((G*coeffs)',size(OBP_ECCO));


        OBP_ECCO_ocean_mascons(nan_ECCO_ind) = 0;

        G = [cos(((2*pi)/365.2425)*time_in_range_ECCO) sin(((2*pi)/365.2425)*time_in_range_ECCO) cos(((2*(2*pi))/365.2425)*time_in_range_ECCO) sin(((2*(2*pi))/365.2425)*time_in_range_ECCO) cos(((3*(2*pi))/365.2425)*time_in_range_ECCO) sin(((3*(2*pi))/365.2425)*time_in_range_ECCO) cos(((4*(2*pi))/365.2425)*time_in_range_ECCO) sin(((4*(2*pi))/365.2425)*time_in_range_ECCO)];
        coeffs = (((G')*G)^(-1))*((G')*(OBP_ECCO_ocean_mascons'));
        OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons - reshape((G*coeffs)',size(OBP_ECCO_ocean_mascons));
    end

    % curr_tseries(curr_nan_ind) = NaN;


    % OBP_ECCO(unique([nan_OBP_ECCO_ind; find(abs(OBP_ECCO_nan_mask_filtered) > 0.2)])) = NaN;

    OBP_ECCO_ocean_mascons(unique([nan_ECCO_ind; find(abs(ECCO_nan_mask_filtered) > 0.2)])) = NaN;
        
    OBP_ECCO_ocean_mascons_array_1 = repmat(reshape(OBP_ECCO_ocean_mascons,[size(OBP_ECCO_ocean_mascons,1) 1 size(OBP_ECCO_ocean_mascons,2)]),[1 size(OBP_ECCO_ocean_mascons,1) 1]);
    OBP_ECCO_ocean_mascons_array_2 = repmat(reshape(OBP_ECCO_ocean_mascons,[1 size(OBP_ECCO_ocean_mascons,1) size(OBP_ECCO_ocean_mascons,2)]),[size(OBP_ECCO_ocean_mascons,1) 1 1]);
    lwe_thickness_ocean_mascons_array_1 = repmat(reshape(lwe_thickness_ocean_mascons,[size(lwe_thickness_ocean_mascons,1) 1 size(lwe_thickness_ocean_mascons,2)]),[1 size(lwe_thickness_ocean_mascons,1) 1]);
    lwe_thickness_ocean_mascons_array_2 = repmat(reshape(lwe_thickness_ocean_mascons,[1 size(lwe_thickness_ocean_mascons,1) size(lwe_thickness_ocean_mascons,2)]),[size(lwe_thickness_ocean_mascons,1) 1 1]);

    delta_lag = 365.24/12;
    lag_range_to_test = (365.24/12)*[0 1];
    
    [OBP_ECCO_ocean_mascons_corr_array,~,~,~,~,~,lags_cov] = correlation_scalar_scalar_uncert_bounds(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
    [GRACE_ocean_mascons_corr_array,~,~,~,~,~,~] = correlation_scalar_scalar_uncert_bounds((9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_1,(9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
    [~,~,~,~,~,~,ECCO_std_dev_array_1,ECCO_std_dev_array_2,~] = regression_linear_scalar_scalar(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);
    [~,~,~,~,~,~,GRACE_std_dev_array_1,GRACE_std_dev_array_2,~] = regression_linear_scalar_scalar((9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_1,(9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);
    
    curr_lag = 0;
    OBP_ECCO_ocean_mascons_corr_array = OBP_ECCO_ocean_mascons_corr_array(:,:,abs(lags_cov - curr_lag) < 1e-5);
    GRACE_ocean_mascons_corr_array = GRACE_ocean_mascons_corr_array(:,:,abs(lags_cov - curr_lag) < 1e-5);
    ECCO_std_dev_array_1 = ECCO_std_dev_array_1(:,:,abs(lags_cov - curr_lag) < 1e-5);
    ECCO_std_dev_array_2 = ECCO_std_dev_array_2(:,:,abs(lags_cov - curr_lag) < 1e-5);
    GRACE_std_dev_array_1 = GRACE_std_dev_array_1(:,:,abs(lags_cov - curr_lag) < 1e-5);
    GRACE_std_dev_array_2 = GRACE_std_dev_array_2(:,:,abs(lags_cov - curr_lag) < 1e-5);
    
    
    time_pt = time_in_range_ECCO;    
    
%     time_pt_datevec = datevec(time_pt);
%     n_months_pt = ((12*time_pt_datevec(size(time_pt_datevec,1),1)) + time_pt_datevec(size(time_pt_datevec,1),2)) - ((12*time_pt_datevec(1,1)) + time_pt_datevec(1,2)) + 1;
    

    OBP_ECCO_zeronans = OBP_ECCO;
    OBP_ECCO_zeronans(isnan(OBP_ECCO) == 1) = 0;
    
    
    OBP_pt_downsc = NaN([length(lon_downsc) length(lat_downsc) size(OBP_ECCO_zeronans,3)]);
    for i_sec = 1:ceil(length(lon_downsc)/10)
        in_sec_range_ind = (max([1 ((10*(i_sec - 1)) + 1)]):1:min([length(lon_downsc) (10*i_sec)]))';
        near_lon_range_ind = find(abs(mod(lon_in_range_ECCO - lon_downsc(floor(mean(in_sec_range_ind))) + 180,360) - 180) <= 3);
        OBP_pt_downsc(in_sec_range_ind,:,:) = reshape(interp2_fft(lon_in_range_ECCO(near_lon_range_ind),lat_in_range_ECCO,OBP_ECCO_zeronans(near_lon_range_ind,:,:),repmat(lon_downsc(in_sec_range_ind),[numel(lat_downsc) 1]),reshape(repmat(lat_downsc,[length(in_sec_range_ind) 1]),[(length(in_sec_range_ind)*length(lat_downsc)) 1])),[length(in_sec_range_ind) length(lat_downsc) size(OBP_ECCO_zeronans,3)]);
        
        disp(['i_sec = ',num2str(i_sec)])
    end
    
%     [~,~,~,~,~,~,OBP_pt_downsc_stddev,~,~] = regression_linear_scalar_scalar(OBP_pt_downsc,OBP_pt_downsc,3,1,1,[0 1],0.95,1);
    
end

lwe_thickness_downsc = NaN([length(lon_downsc) length(lat_downsc)]);
n_mascons_in_downsc = NaN([length(lon_downsc) length(lat_downsc)]);
corr_model_downsc = NaN([length(lon_downsc) length(lat_downsc)]);
for i = 1:length(lon_downsc)
    for j = 1:length(lat_downsc)
        lon_pt = lon_downsc(i);
        lat_pt = lat_downsc(j);
        depth_pt = depth_downsc_grid(i,j);

        if abs(downscale_opt - 1) < 1e-5
            OBP_pt = squeeze(OBP_pt_downsc(i,j,:));
            
            if length(find(isnan(OBP_pt) == 0)) < 0.5*length(OBP_pt)
                continue
            end
                

            curr_tseries_time = time_pt;
            curr_ECCO_tseries = OBP_pt;
    
            curr_tseries_nan_mask = (1e-5)*ones(size(curr_ECCO_tseries));
            curr_tseries_nan_mask((isnan(curr_ECCO_tseries) == 1) | (abs(curr_ECCO_tseries) < 1e-10)) = -1;
            [curr_tseries_filtered,curr_tseries_trend,~] = bandpass_err_fcn(curr_ECCO_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
            [curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask,1,mean(diff(curr_tseries_time)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    
            curr_ECCO_tseries = curr_tseries_filtered;
            clear curr_tseries_filtered
    
    
            nan_curr_tseries_ind = find(isnan(curr_ECCO_tseries) == 1);
            if abs(season_cyc_opt) < 1e-5
                % remove annual cycle
    
                curr_ECCO_tseries(nan_curr_tseries_ind) = 0;
    
                G = [cos(((2*pi)/365.2425)*curr_tseries_time) sin(((2*pi)/365.2425)*curr_tseries_time) cos(((2*(2*pi))/365.2425)*curr_tseries_time) sin(((2*(2*pi))/365.2425)*curr_tseries_time) cos(((3*(2*pi))/365.2425)*curr_tseries_time) sin(((3*(2*pi))/365.2425)*curr_tseries_time) cos(((4*(2*pi))/365.2425)*curr_tseries_time) sin(((4*(2*pi))/365.2425)*curr_tseries_time)];
                coeffs = (((G')*G)^(-1))*((G')*curr_ECCO_tseries);
                curr_ECCO_tseries = curr_ECCO_tseries - (G*coeffs);
            end
    
            % curr_tseries(curr_nan_ind) = NaN;
    
    
            curr_ECCO_tseries(unique([nan_curr_tseries_ind; find(abs(curr_tseries_nan_mask_filtered) > 0.2)])) = NaN;
            
            
            % compute covariances
            
            [OBP_tseries_ECCO_corr_array,~,~,~,~,~,~] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_ECCO_tseries,[1 length(curr_ECCO_tseries)]),[size(OBP_ECCO_ocean_mascons,1) 1]),OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
            [~,~,~,~,~,~,~,std_dev_curr_ECCO_tseries,~] = regression_linear_scalar_scalar(squeeze(OBP_ECCO_ocean_mascons_array_1(:,1,:)),repmat(reshape(curr_ECCO_tseries,[1 length(curr_ECCO_tseries)]),[size(OBP_ECCO_ocean_mascons_array_1,1) 1]),2,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);
            
            
            % use zero lag for covariances
    
            curr_lag = 0;
            OBP_tseries_ECCO_corr_array = OBP_tseries_ECCO_corr_array(:,abs(lags_cov - curr_lag) < 1e-5);
            std_dev_curr_ECCO_tseries = std_dev_curr_ECCO_tseries(:,abs(lags_cov - curr_lag) < 1e-5);
    
    
            % only include mascons within radius of given point
            
            mascons_within_radius_ind = find(abs(((cosd(mean([mascon_lat_center (lat_pt*ones([length(mascon_lat_center) 1]))],2))).*(mod(mascon_lon_center - lon_pt + 180,360) - 180)) + (1i*(mascon_lat_center - lat_pt))) <= radius_mascons_deg);
            curr_OBP_ECCO_ocean_mascons_corr = OBP_ECCO_ocean_mascons_corr_array(mascons_within_radius_ind,mascons_within_radius_ind);
            curr_GRACE_ocean_mascons_corr = GRACE_ocean_mascons_corr_array(mascons_within_radius_ind,mascons_within_radius_ind);
            curr_ECCO_std_dev_array_1 = ECCO_std_dev_array_1(mascons_within_radius_ind,mascons_within_radius_ind);
            curr_ECCO_std_dev_array_2 = ECCO_std_dev_array_2(mascons_within_radius_ind,mascons_within_radius_ind);
            curr_GRACE_std_dev_array_1 = GRACE_std_dev_array_1(mascons_within_radius_ind,mascons_within_radius_ind);
            curr_GRACE_std_dev_array_2 = GRACE_std_dev_array_2(mascons_within_radius_ind,mascons_within_radius_ind);
            OBP_tseries_ECCO_corr_array = OBP_tseries_ECCO_corr_array(mascons_within_radius_ind);
            std_dev_curr_ECCO_tseries = std_dev_curr_ECCO_tseries(mascons_within_radius_ind);
            
            % only include mascons whose correlation with point is above a threshold
            [sorted_tseries_corr,~] = sort(OBP_tseries_ECCO_corr_array,'descend');
            sorted_tseries_corr = sorted_tseries_corr(isnan(sorted_tseries_corr) == 0);
            corr_threshold_n_mascons = sorted_tseries_corr(n_mascons_max) - 1e-10;
            high_corr_ind = find(OBP_tseries_ECCO_corr_array >= max([corr_threshold_n_mascons min_corr_to_include]));
%             disp(['corr_threshold_n_mascons = ',num2str(corr_threshold_n_mascons)])
            n_mascons_in_downsc(i,j) = length(high_corr_ind);
            
            if length(high_corr_ind) < 3
        %         keyboard
                continue
            end            
            
            if abs(downscale_adj_opt - 1) < 1e-5
                
                % adjust mascon-point correlations based on mascon depth relative to obs. point depth
    
                depth_radius_adjust = 2500;     % radius of depth to apply (positive) adjustment
                
                adjust_corr_vec = zeros([length(high_corr_ind) 1]);
                mascon_tseries_depth_diff = abs(mascon_depth_avg(mascons_within_radius_ind(high_corr_ind)) - depth_pt);
                in_range_pos_adjust_ind = find(mascon_tseries_depth_diff < depth_radius_adjust);
                adjust_corr_vec(in_range_pos_adjust_ind) = adjust_corr_vec(in_range_pos_adjust_ind) + (adjust_corr_max*(1 - (mascon_tseries_depth_diff(in_range_pos_adjust_ind)/depth_radius_adjust)));
                OBP_tseries_ECCO_highcorr_array = OBP_tseries_ECCO_corr_array(high_corr_ind) + adjust_corr_vec;
                OBP_tseries_ECCO_highcorr_array(OBP_tseries_ECCO_highcorr_array > 1) = 1;
                OBP_tseries_ECCO_highcorr_array(OBP_tseries_ECCO_highcorr_array < -1) = -1;
                
                ocean_mascons_tseries_cov = OBP_tseries_ECCO_highcorr_array.*(curr_ECCO_std_dev_array_1(high_corr_ind,1)).*std_dev_curr_ECCO_tseries(high_corr_ind);
                
            else
                ocean_mascons_tseries_cov = (OBP_tseries_ECCO_corr_array(high_corr_ind)).*(curr_ECCO_std_dev_array_1(high_corr_ind,1)).*std_dev_curr_ECCO_tseries(high_corr_ind);
            end
    
            ocean_mascons_cov = (curr_OBP_ECCO_ocean_mascons_corr(high_corr_ind,high_corr_ind)).*(curr_ECCO_std_dev_array_1(high_corr_ind,high_corr_ind)).*(curr_ECCO_std_dev_array_2(high_corr_ind,high_corr_ind));
            
            gain_vec = (ocean_mascons_cov^(-1))*ocean_mascons_tseries_cov;
            lwe_thickness_downsc(i,j) = (lwe_thickness_ocean_mascons(mascons_within_radius_ind(high_corr_ind),time_plot_ind)')*gain_vec;
            
            % entirely model-based reconstruction
            lwe_thickness_model_downsc = (OBP_ECCO_ocean_mascons(mascons_within_radius_ind(high_corr_ind),:)')*gain_vec;
            [model_reconstr_corr_array,~,~,~,~,~,~] = correlation_scalar_scalar_uncert_bounds(curr_ECCO_tseries,lwe_thickness_model_downsc,1,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
            corr_model_downsc(i,j) = model_reconstr_corr_array(abs(lags_cov - curr_lag) < 1e-5);
            
        else
            colloc_mascon_ind = find((mod(lon_pt - mascon_lon_bounds(:,1) + 180,360) - 180 >= 0) & (mod(lon_pt - mascon_lon_bounds(:,2) + 180,360) - 180 < 0) & (lat_pt - mascon_lat_bounds(:,1) >= 0) & (lat_pt - mascon_lat_bounds(:,2) < 0));
            lwe_thickness_downsc(i,j) = lwe_thickness_ocean_mascons(colloc_mascon_ind,time_plot_ind);
        end
        
    end
    
    disp(['Completed i = ',num2str(i),' (out of ',num2str(length(lon_downsc)),')'])
end


% plot maps

in_plot_orig_grid_lon_ind = find((mod(lon - lon_bounds(1) + 180,360) - 180 >= -0.5) & (mod(lon - lon_bounds(2) + 180,360) - 180 <= 0.5));
in_plot_orig_grid_lat_ind = find((lat - lat_bounds(1) >= -0.5) & (lat - lat_bounds(2) <= 0.5));
lon_orig_grid_SW_corner = lon(in_plot_orig_grid_lon_ind) - (diff(lon((min(in_plot_orig_grid_lon_ind) - 1):max(in_plot_orig_grid_lon_ind)))/2);
lat_orig_grid_SW_corner = lat(in_plot_orig_grid_lat_ind) - (diff(lat((min(in_plot_orig_grid_lat_ind) - 1):max(in_plot_orig_grid_lat_ind)))/2);
lon_orig_grid_SW_corner = lon_orig_grid_SW_corner + (360*round((lon_bounds(1) - lon_orig_grid_SW_corner)/360));

in_plot_mascon_center_ind = find((mod(mascon_lon_center - lon_bounds(1) + 180,360) - 180 >= -0.5) & (mod(mascon_lon_center - lon_bounds(2) + 180,360) - 180 <= 0.5) & (mascon_lat_center - lat_bounds(1) >= -0.5) & (mascon_lat_center - lat_bounds(2) <= 0.5));


lwe_thickness_ocean_time_plot = NaN(size(curr_lwe_thickness));
for curr_mascon_ind = 1:length(good_mascon_ind)
    in_curr_mascon_ind = find(abs(curr_lwe_thickness - ocean_mascon_curr_values(good_mascon_ind(curr_mascon_ind))) < 1e-5);
    lwe_thickness_ocean_time_plot(in_curr_mascon_ind) = lwe_thickness_ocean_mascons(curr_mascon_ind,time_plot_ind);
end


c_levels = [(-50) ((-20):1:20) 50];     % specified levels for contours
% cmap = colormap(0.85*bcyr(40));    % colormap for contours
cmap = colormap(0.85*bcyr(42));    % colormap for contours
% % cmap = flip(cmap([(1:1:10) (12:2:18) (20:1:31) (33:2:39) (41:1:50)],:),1);
% cmap = flip(cmap,1);

% plot pcolor map with land mask

[output_array,output_clabels] = irregular_clevels_plot(10*lwe_thickness_ocean_time_plot(in_plot_orig_grid_lon_ind,in_plot_orig_grid_lat_ind),c_levels);

fig1000 = figure(1000);
close(figure(1))
colormap(cmap)
fig_paper_pos = get(fig1000,'PaperPosition');
fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
p = pcolor(lon_orig_grid_SW_corner,lat_orig_grid_SW_corner',output_array');
caxis([0 size(cmap,1)])
set(p,'edgecolor','none')
set(gca,'xlim',lon_bounds,'ylim',lat_bounds)
daspect([1 cosd(mean(lat_bounds)) 1])
set(gca,'position',[0 0 1 1],'units','normalized','Visible','off')
pause(5)
plotboxaspectratio = get(gca,'PlotBoxAspectRatio');
fig_pos = get(gcf,'Position');
ratio_plotbox = plotboxaspectratio(2)/plotboxaspectratio(1);
ratio_fig = fig_pos(4)/fig_pos(3);
if ratio_fig > ratio_plotbox
    fig_pos(4) = ratio_plotbox*fig_pos(3);
else
    fig_pos(3) = fig_pos(4)/ratio_plotbox;
end
set(gcf,'Position',fig_pos)

print(gcf,'contour_plot_temporary.png','-dpng','-r300')
close(fig1000)

% get contour plot and overlay land mask

rgb_array = imread('contour_plot_temporary.png');
rgb_array = flip(permute(rgb_array,[2 1 3]),2);

delete('contour_plot_temporary.png')

size_rgb_array = size(rgb_array);
landmask_ind = landfind_indices(lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2),size_rgb_array(1:2));
rgb_array_reshaped = reshape(rgb_array,[prod(size_rgb_array(1:2)) 3]);

black_shaded_ind = find((rgb_array_reshaped(:,1) == 0) & (rgb_array_reshaped(:,2) == 0) & (rgb_array_reshaped(:,3) == 0));
rgb_array_reshaped(black_shaded_ind,:) = 255*ones(length(black_shaded_ind),3);

% put black mask on land areas
rgb_array_reshaped(landmask_ind,:) = zeros(length(landmask_ind),3);

rgb_array_masked = reshape(rgb_array_reshaped,size_rgb_array);

fig1 = figure(1);
h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
set(gca,'ydir','normal')
daspect([1 cosd(mean(lat_bounds)) 1])
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
marker_radius = 0.2;
for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
    curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
    rectangle('Position',[(mascon_lon_center(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon_center(curr_mascon_ind))/360)) - marker_radius) (mascon_lat_center(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
end
hold off
xlabel('Longitude','FontSize',10)
ylabel('Latitude','FontSize',10)
title({['Equivalent water thickness anomaly from JPL ocean mascons, ',datestr(datenum([time_to_plot 16 0 0 0]),'mmm yyyy'),',']; 'bathymetry CI = 1000 m'; ' '},'FontSize',10)
colormap(cmap)
cbar = colorbar('location','southoutside');
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'xtick',(1:10:size(cmap,1)),'xticklabel',output_clabels(1 + (1:10:size(cmap,1))),'FontSize',10)
set(get(cbar,'xlabel'),'String','Equivalent water thickness anomaly (mm)','FontSize',10)

print(fig1,['GRACE_mascon_lwe_thickness_',datestr(datenum([time_to_plot 16 0 0 0]),'yyyymm'),'_time_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat.pdf'],'-dpdf','-r800')
close(fig1)


% contour plot of interpolated field

[output_array,output_clabels] = irregular_clevels_plot(10*lwe_thickness_downsc,c_levels);

fig2 = figure(2);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig2,lon_downsc_grid',lat_downsc_grid',output_array',0:1:(length(c_levels) - 1),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)])
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% marker_radius = 0.2;
% for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
%     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
%     rectangle('Position',[(mascon_lon_center(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon_center(curr_mascon_ind))/360)) - marker_radius) (mascon_lat_center(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% end
hold off
xlabel('Longitude','FontSize',10)
ylabel('Latitude','FontSize',10)
title({['Downscaled GRACE water thickness anomaly, ',datestr(datenum([time_to_plot 16 0 0 0]),'mmm yyyy'),',']; ' '},'FontSize',10)
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'xtick',(1:10:size(cmap,1)),'xticklabel',output_clabels(1 + (1:10:size(cmap,1))),'FontSize',10)
set(get(cbar,'xlabel'),'String','Equivalent water thickness anomaly (mm)','FontSize',10)

print(fig2,['GRACE_downsc_lwe_thickness_',datestr(datenum([time_to_plot 16 0 0 0]),'yyyymm'),'_time_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',downscale_id,'.pdf'],'-dpdf','-r800')
close(fig2)


% plot number of mascons involved in downscaling at each point

c_levels = ((-0.5):1:10.5)';
[output_array,~] = irregular_clevels_plot(n_mascons_in_downsc,c_levels);
cmap = 0.85*bcyr(11);
cmap = [cmap(1:10,:); [0.85 0 0]];

fig3 = figure(3);
close(figure(1))
colormap(cmap)
[h,cbar] = contour_plot_with_landmask(fig3,lon_downsc_grid',lat_downsc_grid',output_array',0:1:(length(c_levels) - 1),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)])
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% marker_radius = 0.2;
% for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
%     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
%     rectangle('Position',[(mascon_lon_center(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon_center(curr_mascon_ind))/360)) - marker_radius) (mascon_lat_center(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% end
hold off
xlabel('Longitude','FontSize',10)
ylabel('Latitude','FontSize',10)
title({'Number of mascons involved in downscaling'; ' '},'FontSize',10)
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'xtick',0.5:1:10.5,'xticklabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'},'FontSize',10)
set(get(cbar,'xlabel'),'String','Number of mascons','FontSize',10)

print(fig3,['GRACE_downsc_n_mascons_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',downscale_id,'.bmp'],'-dbmp','-r300')
close(fig3)



% plot correlation of reconstruction based on model point OBP (no actual GRACE data)

c_levels = [-1 (0:0.1:1)]';
[output_array,~] = irregular_clevels_plot(corr_model_downsc,c_levels);

cmap = 0.85*bcyr(28);
cmap = [0.4 0.4 0.85; cmap([16; 17; 18; 19; 21; 24; 26; 28],:)];

fig4 = figure(4);
close(figure(1))
colormap(cmap)
[~,cbar] = contour_plot_with_landmask(fig4,lon_downsc_grid',lat_downsc_grid',output_array',0:1:(length(c_levels) - 1),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)])
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% marker_radius = 0.2;
% for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
%     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
%     rectangle('Position',[(mascon_lon_center(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon_center(curr_mascon_ind))/360)) - marker_radius) (mascon_lat_center(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% end
hold off
xlabel('Longitude','FontSize',10)
ylabel('Latitude','FontSize',10)
title({'Correlation of model-based downscaling reconstruction'; ' '},'FontSize',10)
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'xtick',1:1:11,'xticklabel',{'0' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'},'FontSize',10)
set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',10)

print(fig4,['GRACE_downsc_model_reconstr_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',downscale_id,'.bmp'],'-dbmp','-r300')
close(fig4)