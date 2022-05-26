% compute volume flux across a given latitude and within a certain layer, with downscaled GRACE

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
cd('/indopac/adelman/GRACE/')


lat_transect = 26.5;
lon_bounds = [-81.0 -8.0];
lon_transect_id = 'Atlantic';
depth_bounds = [1100 3000];
time_range_start = [2002 1 1];
time_range_end = [2018 1 1];

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

% downscale_id = 'coloc';    % identifier of method used to define covariances for downscaling
% downscale_id = ['cs510_cylindweight1.5_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
downscale_id = ['cs510_cylindweight1.5_depthadj',num2str(adjust_corr_max),'_hybrid',num2str(hybrid_factor),'_maxmascons',num2str(n_mascons_max),'_mincorr',num2str(min_corr_to_include)];    % identifier of method used to define covariances for downscaling
downscale_opt = 1;     % 0 = no downscaling (just correlate with co-located mascon); 1 = downscaling
downscale_adj_opt = 1;   % 0 = no adjustment to spatial correlation/covariance; 1 = adjust spatial correlation

% normalized error tolerance
norm_err_tolerance = 0.3;


mascon_lat_separation = 3;    % in degrees

% curr_file = 'CLM4.SCALE_FACTOR.JPL.MSCNv01CRIv01.nc';
% lon = ncread(curr_file,'lon');
% lat = ncread(curr_file,'lat');
% scale_factor = ncread(curr_file,'scale_factor');
% curr_file = 'JPL_MSCNv01_PLACEMENT.nc';
% mascon_lon = ncread(curr_file,'mascon_lon');
% mascon_lat = ncread(curr_file,'mascon_lat');
% mascon_rad = ncread(curr_file,'mascon_rad');
curr_file = 'LAND_MASK.CRIv01.nc';
land_mask = ncread(curr_file,'land_mask');
curr_file = 'GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';
lon = ncread(curr_file,'lon');
lat = ncread(curr_file,'lat');
time = ncread(curr_file,'time') + datenum([2002 1 1 0 0 0]);
in_time_range_ind = find((time >= datenum(time_range_start)) & (time < datenum(time_range_end)));
time = time(in_time_range_ind);
lwe_thickness = ncread(curr_file,'lwe_thickness',[1 1 min(in_time_range_ind)],[length(lon) length(lat) length(in_time_range_ind)]);
lwe_uncert = ncread(curr_file,'uncertainty',[1 1 min(in_time_range_ind)],[length(lon) length(lat) length(in_time_range_ind)]);
% lon_bounds = ncread(curr_file,'lon_bounds');
% lat_bounds = ncread(curr_file,'lat_bounds');
% time_bounds = ncread(curr_file,'time_bounds');


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



% find mascons near map region

in_region_ind = find(((cosd(mascon_lat_center_all)).*(mod(mascon_lon_center_all - lon_bounds(1) + 180,360) - 180) >= -radius_mascons_deg - 1) & ((cosd(mascon_lat_center_all)).*(mod(mascon_lon_center_all - lon_bounds(2) + 180,360) - 180) <= radius_mascons_deg + 1) & (mascon_lat_center_all - lat_transect >= -radius_mascons_deg - 1) & (mascon_lat_center_all - lat_transect <= radius_mascons_deg + 1));

ocean_mascon_curr_values = ocean_mascon_curr_values_all(in_region_ind);
mascon_lon_center = mascon_lon_center_all(in_region_ind);
mascon_lat_center = mascon_lat_center_all(in_region_ind);
mascon_lon_bounds = mascon_lon_bounds_all(in_region_ind,:);
mascon_lat_bounds = mascon_lat_bounds_all(in_region_ind,:);
lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons_all(in_region_ind,:);
lwe_uncert_ocean_mascons = lwe_uncert_ocean_mascons_all(in_region_ind,:);
mascon_depth_avg = mascon_depth_avg_all(in_region_ind);



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
% lwe_uncert_ocean_mascons_filtered = squeeze((sum((abs(filter_gain_coeffs_array).^2).*(repmat(reshape(lwe_uncert_ocean_mascons,[size(lwe_uncert_ocean_mascons,1) 1 size(lwe_uncert_ocean_mascons,2)]).^2,[1 size(lwe_uncert_ocean_mascons,2) 1])),3)).^(1/2));
% clear filter_gain_coeffs_array


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



% interpolate bathymetry to transect latitude, and focus on specified longitudes

transect_lat_bathy_closest_before_ind = find(lat_bathy <= lat_transect,1,'last');
transect_lat_bathy_closest_after_ind = find(lat_bathy > lat_transect,1,'first');
weight_after = (lat_transect - lat_bathy(transect_lat_bathy_closest_before_ind))/(diff(lat_bathy([transect_lat_bathy_closest_before_ind transect_lat_bathy_closest_after_ind])));
weight_before = 1 - weight_after;

if ((mod(lon_bounds(2) + (2*spacing) - lon_bounds(1) - 1e-5,360) + 1e-5 > mod(max(lon_bathy) - lon_bounds(1) - 1e-5,360) + 1e-5) || (mod(lon_bounds(2) - (lon_bounds(1) - (2*spacing)) - 1e-5,360) + 1e-5 > mod(lon_bounds(2) - min(lon_bathy) - 1e-5,360) + 1e-5))
    in_bathy_lon_range_ind_1 = find(lon_bathy - (mod(lon_bounds(1) - min(lon_bathy),360) + min(lon_bathy)) >= (-2*spacing));
    in_bathy_lon_range_ind_2 = find(lon_bathy - (mod(lon_bounds(2) - min(lon_bathy),360) + min(lon_bathy)) <= (2*spacing));
    lon_bathy_along_transect = [lon_bathy(in_bathy_lon_range_ind_1); (lon_bathy(in_bathy_lon_range_ind_2) + 360)];
    z_along_transect = (weight_before*z_reshaped([in_bathy_lon_range_ind_1; in_bathy_lon_range_ind_2],transect_lat_bathy_closest_before_ind)) + (weight_after*z_reshaped([in_bathy_lon_range_ind_1; in_bathy_lon_range_ind_2],transect_lat_bathy_closest_after_ind));
else
    in_bathy_lon_range_ind = find(((lon_bathy - (mod(lon_bounds(1) - min(lon_bathy),360) + min(lon_bathy))) >= (-2*spacing)) & (lon_bathy - (mod(lon_bounds(2) - min(lon_bathy),360) + min(lon_bathy)) <= (2*spacing)));
    lon_bathy_along_transect = lon_bathy(in_bathy_lon_range_ind);
    z_along_transect = (weight_before*z_reshaped(in_bathy_lon_range_ind,transect_lat_bathy_closest_before_ind)) + (weight_after*z_reshaped(in_bathy_lon_range_ind,transect_lat_bathy_closest_after_ind));
end


% define depth levels to use in integration

depth_level_bounds_all = [(0:10:100)'; (120:20:200)'; (240:40:400)'; (450:50:1000)'; (1100:100:6000)'];

depth_ind_inrange = find((depth_level_bounds_all(2:length(depth_level_bounds_all)) > depth_bounds(1)) & (depth_level_bounds_all(1:(length(depth_level_bounds_all) - 1)) < depth_bounds(2)));

depth_level_bounds = [depth_bounds(1); depth_level_bounds_all((min(depth_ind_inrange) + 1):max(depth_ind_inrange)); depth_bounds(2)];
depth_level_centers = depth_level_bounds(2:length(depth_level_bounds)) - (diff(depth_level_bounds)/2);

lon_endpoints_at_depth = cell([length(depth_level_centers) 1]);
for depth_ind = 1:length(depth_level_centers)
    curr_depth_center = depth_level_centers(depth_ind);
    
    depth_mask = zeros(size(z_along_transect));
    depth_mask(z_along_transect < -curr_depth_center) = 1;
    
    lon_start_ind = find(diff(depth_mask) > 1e-5);
    lon_end_ind = find(diff(depth_mask) < -1e-5);
    
    weight_after = (curr_depth_center - (-z_along_transect(lon_start_ind)))./((-z_along_transect(lon_start_ind + 1)) - (-z_along_transect(lon_start_ind)));
    weight_before = 1 - weight_after;
    lon_start = (weight_before.*lon_bathy_along_transect(lon_start_ind)) + (weight_after.*lon_bathy_along_transect(lon_start_ind + 1));
    
    weight_after = (curr_depth_center - (-z_along_transect(lon_end_ind)))./((-z_along_transect(lon_end_ind + 1)) - (-z_along_transect(lon_end_ind)));
    weight_before = 1 - weight_after;
    lon_end = (weight_before.*lon_bathy_along_transect(lon_end_ind)) + (weight_after.*lon_bathy_along_transect(lon_end_ind + 1));
    
    if isempty(lon_start_ind) == 1
        lon_start = mod(lon_bounds(1),360);
    end
    if isempty(lon_end_ind) == 1
        lon_end = mod(lon_bounds(2) - (max(lon_bathy_along_transect) + 1e-5 - 360),360) + (max(lon_bathy_along_transect) + 1e-5 - 360);
    end
    try
        if lon_end_ind(1) < lon_start_ind(1)
            disp(['Transect starts in open water at depth = ',num2str(curr_depth_center),' m'])
            lon_start = [mod(lon_bounds(1),360); lon_start];
        end
        if lon_start_ind(length(lon_start_ind)) > lon_end_ind(length(lon_end_ind))
            disp(['Transect ends in open water at depth = ',num2str(curr_depth_center),' m'])
            lon_end = [lon_end; (mod(lon_bounds(2) - (max(lon_bathy_along_transect) + 1e-5 - 360),360) + (max(lon_bathy_along_transect) + 1e-5 - 360))];
        end
    catch
    end
    
    gap_low_width_threshold = 1;    % minimum width of gap, in degrees longitude
    long_enough_ind = find(diff([lon_start lon_end],1,2) > 1);
    
    lon_endpoints_at_depth{depth_ind} = [lon_start(long_enough_ind) lon_end(long_enough_ind)];        

end


f = 2*((2*pi)/86164)*sind(lat_transect);


if abs(downscale_opt - 1) < 1e-5

    % load ECCO OBP for covariance calculations

    ECCO_nc_file = '/indopac/adelman/ECCO2/PHIBOT.ECCO2.lonlatinterp.1992-2018.nc';

    longitude = ncread(ECCO_nc_file,'LONGITUDE_T');
    latitude = ncread(ECCO_nc_file,'LATITUDE_T');
    time_ECCO = ncread(ECCO_nc_file,'TIME');
    time_ECCO = double(time_ECCO) + datenum([1992 1 1 0 0 0]);

    in_lon_range_ind = find((cosd(lat_transect).*(mod(longitude - lon_bounds(1) + 180,360) - 180) >= -radius_mascons_deg - 5) & (cosd(mean(lat_transect)).*(mod(longitude - lon_bounds(2) + 180,360) - 180) <= radius_mascons_deg + 5));
    in_lat_range_ind = find((latitude - lat_transect >= -radius_mascons_deg - 5) & (latitude - lat_transect <= radius_mascons_deg + 5));
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


    % temporally filter time series

    steepness_factor = 5;
    half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)


    ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO_ocean_mascons));
    ECCO_nan_mask((isnan(OBP_ECCO_ocean_mascons) == 1) | (abs(OBP_ECCO_ocean_mascons) < 1e-10)) = -1;
    [OBP_ECCO_ocean_mascons_filtered,OBP_ECCO_ocean_mascons_trend,~] = bandpass_err_fcn(OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range_ECCO)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    [ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,2,mean(diff(time_in_range_ECCO)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

    OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons_filtered;
    clear OBP_ECCO_ocean_mascons_filtered

    nan_ECCO_ind = find(isnan(OBP_ECCO_ocean_mascons) == 1);
    if abs(season_cyc_opt) < 1e-5
        % remove annual cycle


        OBP_ECCO_ocean_mascons(nan_ECCO_ind) = 0;

        G = [cos(((2*pi)/365.2425)*time_in_range_ECCO) sin(((2*pi)/365.2425)*time_in_range_ECCO) cos(((2*(2*pi))/365.2425)*time_in_range_ECCO) sin(((2*(2*pi))/365.2425)*time_in_range_ECCO) cos(((3*(2*pi))/365.2425)*time_in_range_ECCO) sin(((3*(2*pi))/365.2425)*time_in_range_ECCO) cos(((4*(2*pi))/365.2425)*time_in_range_ECCO) sin(((4*(2*pi))/365.2425)*time_in_range_ECCO)];
        coeffs = (((G')*G)^(-1))*((G')*(OBP_ECCO_ocean_mascons'));
        OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons - reshape((G*coeffs)',size(OBP_ECCO_ocean_mascons));
    end


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
    

    OBP_ECCO_zeronans = OBP_ECCO;
    OBP_ECCO_zeronans(isnan(OBP_ECCO) == 1) = 0;
    
    
end
    
    
vol_flux_segments = cell([length(depth_level_centers) 1]);
vol_flux_in_layers = NaN([length(depth_level_centers) length(time)]);
vol_flux_stddev_segments = cell([length(depth_level_centers) 1]);
% vol_flux_uncert_segments = cell([length(depth_level_centers) 1]);
% vol_flux_uncert_in_layers = NaN([length(depth_level_centers) length(time)]);
cum_n_segments = 0;
% cum_mascon_ind = [];
max_stddev_segment = [];
for depth_ind = 1:length(depth_level_centers)
    curr_depth_center = depth_level_centers(depth_ind);
    curr_depth_thickness = diff(depth_level_bounds(depth_ind + [0 1]));
    curr_lon_endpoints = lon_endpoints_at_depth{depth_ind};
    
%     cum_n_segments = cum_n_segments + size(curr_lon_endpoints,1);
    cum_n_segments = cum_n_segments + sum(min([ones([size(curr_lon_endpoints,1) 1]) (diff(curr_lon_endpoints,1,2)/(mascon_lat_separation/cosd(lat_transect)))],[],2),1);
    
    curr_endpoint_lwe_thickness = NaN([size(curr_lon_endpoints) size(lwe_thickness_ocean_mascons,2)]);
    curr_endpoint_lwe_uncert = NaN([size(curr_lon_endpoints) size(lwe_uncert_ocean_mascons,2)]);
    for curr_endpt = 1:numel(curr_lon_endpoints)
        curr_lon = curr_lon_endpoints(curr_endpt);
        
        if abs(downscale_opt - 1) < 1e-5
            near_lat_range_ind = find(abs(lat_in_range_ECCO - lat_transect) <= 3);
            near_lon_range_ind = find(abs(mod(lon_in_range_ECCO - curr_lon + 180,360) - 180) <= 3);
            OBP_pt = squeeze(interp2_fft(lon_in_range_ECCO(near_lon_range_ind),lat_in_range_ECCO(near_lat_range_ind),OBP_ECCO_zeronans(near_lon_range_ind,near_lat_range_ind,:),mod(curr_lon - (mean(lon_in_range_ECCO(near_lon_range_ind)) - 180),360) + (mean(lon_in_range_ECCO(near_lon_range_ind)) - 180),lat_transect));
            
            
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
    
    
            % only include mascons within radius of given transect
            
            mascons_within_radius_ind = find(abs(((cosd(mean([mascon_lat_center (lat_transect*ones([length(mascon_lat_center) 1]))],2))).*(mod(mascon_lon_center - curr_lon + 180,360) - 180)) + (1i*(mascon_lat_center - lat_transect))) <= radius_mascons_deg);
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
            if length(high_corr_ind) < 3
        %         keyboard
                continue
            end            
            
            if abs(downscale_adj_opt - 1) < 1e-5
                
                % adjust mascon-point correlations based on mascon depth relative to obs. point depth
    
                depth_radius_adjust = 2500;     % radius of depth to apply (positive) adjustment
                
                adjust_corr_vec = zeros([length(high_corr_ind) 1]);
                mascon_tseries_depth_diff = abs(mascon_depth_avg(mascons_within_radius_ind(high_corr_ind)) - curr_depth_center);
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
            curr_i = mod(curr_endpt - 1,size(curr_lon_endpoints,1)) + 1;
            curr_j = ceil(curr_endpt/size(curr_lon_endpoints,1));
            curr_endpoint_lwe_thickness(curr_i,curr_j,:) = (gain_vec')*(lwe_thickness_ocean_mascons(mascons_within_radius_ind(high_corr_ind),:));
            
        else
            colloc_mascon_ind = find((mod(curr_lon - mascon_lon_bounds(:,1) + 180,360) - 180 >= 0) & (mod(curr_lon - mascon_lon_bounds(:,2) + 180,360) - 180 < 0) & (lat_transect - mascon_lat_bounds(:,1) >= 0) & (lat_transect - mascon_lat_bounds(:,2) < 0));
            curr_i = mod(curr_endpt - 1,size(curr_lon_endpoints,1)) + 1;
            curr_j = ceil(curr_endpt/size(curr_lon_endpoints,1));
            curr_endpoint_lwe_thickness(curr_i,curr_j,:) = reshape(lwe_thickness_ocean_mascons(colloc_mascon_ind,:),[1 1 size(lwe_thickness_ocean_mascons,2)]);
        end
        
    end
    
    % estimated density at current depth level (not based on actual water characteristics)
    curr_depth_est_density = 1027 + ((0.0049 - ((6e-8)*curr_depth_center))*curr_depth_center);
    
    vol_flux_segments{depth_ind} = reshape(diff(1000*9.81*(0.01*curr_endpoint_lwe_thickness),1,2)*curr_depth_thickness,[size(curr_lon_endpoints,1) length(time)])/(curr_depth_est_density*f);
    vol_flux_in_layers(depth_ind,:) = sum(vol_flux_segments{depth_ind},1);
    
    
    nan_mask_vol_flux_segments = ones(size(vol_flux_segments{depth_ind}));
    nan_mask_vol_flux_segments(isnan(vol_flux_segments{depth_ind}) == 1) = 0;
    sum_nan_mask = sum(nan_mask_vol_flux_segments,2);
    nan_mask_vol_flux_segments(repmat(sum_nan_mask,[1 size(vol_flux_segments{depth_ind},2)]) < 0.8*max(sum_nan_mask)) = 0;
    
    vol_flux_segments{depth_ind}(isnan(vol_flux_segments{depth_ind}) == 1) = 0;
    
    vol_flux_segments_tmean = sum(nan_mask_vol_flux_segments.*vol_flux_segments{depth_ind},2)./sum(nan_mask_vol_flux_segments,2);
    vol_flux_stddev_segments{depth_ind} = (sum(nan_mask_vol_flux_segments.*((vol_flux_segments{depth_ind} - repmat(vol_flux_segments_tmean,[1 size(vol_flux_segments{depth_ind},2)])).^2),2)./(sum(nan_mask_vol_flux_segments,2) - 1)).^(1/2);
    vol_flux_segments{depth_ind}(abs(vol_flux_segments{depth_ind}) < 1e-10) = NaN;
    vol_flux_stddev_segments{depth_ind}(isnan(vol_flux_stddev_segments{depth_ind}) == 1) = 0;
    max_stddev_segment = max([max_stddev_segment; vol_flux_stddev_segments{depth_ind}]);
    
%     vol_flux_uncert_segments{depth_ind} = squeeze(1000*9.81*(0.01*((sum(curr_endpoint_lwe_uncert.^2,2)).^(1/2)))*curr_depth_thickness)/(curr_depth_est_density*f);
%     vol_flux_uncert_in_layers(depth_ind,:) = (sum(vol_flux_uncert_segments{depth_ind}.^2,1)).^(1/2);
    
end

% fill in NaN layers
nan_mask_layer_flux = ones(size(vol_flux_in_layers));
nan_mask_layer_flux(isnan(vol_flux_in_layers) == 1) = 0;
good_depth_ind = find(sum(nan_mask_layer_flux,2) > 0.8*max(sum(nan_mask_layer_flux,2)));
if isempty(setdiff((1:1:size(nan_mask_layer_flux,1))',good_depth_ind)) == 0
    vol_flux_in_layers_interp = interp1(depth_level_centers(good_depth_ind),vol_flux_in_layers(good_depth_ind,:),depth_level_centers);
    [dist_to_closest_good_depth,closest_ind] = min(repmat((1:1:size(vol_flux_in_layers,1))',[1 length(good_depth_ind)]) - repmat(good_depth_ind',[size(vol_flux_in_layers,1) 1]),[],2);
    vol_flux_in_layers(abs(abs(dist_to_closest_good_depth) - 1.5) < 1,:) = vol_flux_in_layers_interp(abs(abs(dist_to_closest_good_depth) - 1.5) < 1,:);
    
    nan_mask_layer_flux = ones(size(vol_flux_in_layers_interp));
    nan_mask_layer_flux(isnan(vol_flux_in_layers_interp) == 1) = 0;
    close_outside_ind = intersect(find(sum(nan_mask_layer_flux,2) < (0.8*max(sum(nan_mask_layer_flux,2)))),find(abs(dist_to_closest_good_depth) < 2.5));
    vol_flux_in_layers(close_outside_ind) = vol_flux_in_layers(good_depth_ind(closest_ind(close_outside_ind)));
    vol_flux_in_layers(isnan(vol_flux_in_layers) == 1) = 0;
end

vol_flux_total = sum(vol_flux_in_layers,1);


nan_mask_vol_flux_in_layers = ones(size(vol_flux_in_layers));
nan_mask_vol_flux_in_layers(isnan(vol_flux_in_layers) == 1) = 0;
sum_nan_mask = sum(nan_mask_vol_flux_in_layers,2);
nan_mask_vol_flux_in_layers(repmat(sum_nan_mask,[1 size(vol_flux_in_layers,2)]) < 0.8*max(sum_nan_mask)) = 0;
    
vol_flux_in_layers(isnan(vol_flux_in_layers) == 1) = 0;
% vol_flux_uncert_in_layers(isnan(vol_flux_uncert_in_layers) == 1) = 0;

vol_flux_in_layers_tmean = sum(nan_mask_vol_flux_in_layers.*vol_flux_in_layers,2)./sum(nan_mask_vol_flux_in_layers,2);
% vol_flux_in_layers_tmean_uncert = (sum(nan_mask_vol_flux_in_layers.*(vol_flux_uncert_in_layers.^2),2)./sum(nan_mask_vol_flux_in_layers,2)).^(1/2);
vol_flux_in_layers_tmean_perdepth = vol_flux_in_layers_tmean./diff(depth_level_bounds);
% vol_flux_in_layers_tmean_uncert_perdepth = vol_flux_in_layers_tmean_uncert./diff(depth_level_bounds);

vol_flux_in_layers(abs(vol_flux_in_layers) < 1e-10) = NaN;
% vol_flux_uncert_in_layers(abs(vol_flux_uncert_in_layers) < 1e-10) = NaN;


% % scaling factor for uncertainty in depth integration: the square root of the ratio between
% % (2*cum_n_segments) and an estimate of the actual dof based on the number
% % of mascons used in the calculation
% 
% % mascon_in_closest_lat_band_ind = find((mascon_lat - lat_transect >= -(mascon_lat_separation/2)) & (mascon_lat - lat_transect < mascon_lat_separation/2));
% % est_actual_dof_mascon = (1/(2*0.8))*length(mascon_in_closest_lat_band_ind);
% est_actual_dof_mascon = (1/(2*0.8))*length(cum_mascon_ind);
% uncert_scaling_factor = ((2*cum_n_segments)/est_actual_dof_mascon)^(1/2);
% 
% vol_flux_uncert_total = uncert_scaling_factor*((sum(vol_flux_uncert_in_layers.^2,1)).^(1/2));



% % identifiers for obs. time series
% obs_tseries_text_id = 'SAMOC upper branch';
% obs_tseries_id = 'SAMOC_upper';


% load RAPID data

curr_file = '/indopac/adelman/GRACE/RAPID/moc_transports.nc';
RAPID_time = ncread(curr_file,'time') + datenum([2004 04 01 0 0 0]);
RAPID_transport_0_800 = ncread(curr_file,'t_therm10');
RAPID_transport_800_1100 = ncread(curr_file,'t_aiw10');
RAPID_transport_1100_3000 = ncread(curr_file,'t_ud10');
RAPID_transport_3000_5000 = ncread(curr_file,'t_ld10');
RAPID_transport_gt5000 = ncread(curr_file,'t_bw10');
RAPID_transport_FloridaStr = ncread(curr_file,'t_gs10');
RAPID_transport_Ekman = ncread(curr_file,'t_ek10');
RAPID_transport_upper_midocean = ncread(curr_file,'t_umo10');
RAPID_transport_overturning = ncread(curr_file,'moc_mar_hc10');


% % remove linear trend (and bandpass if desired)
% 
% % steepness_factor = 5;
% % low_freq_bound = 1/(1.5*mean(diff(time))*length(time));
% % high_freq_bound = 1/426;
% half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)
% edge_handling_opt = 0;    % mask out edges where filtering is less accurate? (0 = no, 1 = yes)
% 
% [RAPID_transport_3000_5000_filtered,RAPID_transport_3000_5000_trend,~] = bandpass_err_fcn(RAPID_transport_3000_5000,1,mean(diff(RAPID_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,edge_handling_opt,1);
% 
% RAPID_transport_3000_5000 = RAPID_transport_3000_5000_filtered;
% clear RAPID_transport_3000_5000_filtered
% 
% 
% % remove annual cycle
% 
% curr_nan_ind = find(isnan(RAPID_transport_3000_5000) == 1);
% RAPID_transport_3000_5000(curr_nan_ind) = 0;
% 
% G = [cos(((2*pi)/365.2425)*RAPID_time) sin(((2*pi)/365.2425)*RAPID_time) cos(((2*(2*pi))/365.2425)*RAPID_time) sin(((2*(2*pi))/365.2425)*RAPID_time) cos(((3*(2*pi))/365.2425)*RAPID_time) sin(((3*(2*pi))/365.2425)*RAPID_time) cos(((4*(2*pi))/365.2425)*RAPID_time) sin(((4*(2*pi))/365.2425)*RAPID_time)];
% coeffs = (((G')*G)^(-1))*((G')*RAPID_transport_3000_5000);
% RAPID_transport_3000_5000 = RAPID_transport_3000_5000 - (G*coeffs);
% 
% RAPID_transport_3000_5000(curr_nan_ind) = NaN;



% plot bathymetry profile, depth vs. flux, and time series of flux

bcyr_40 = colormap(bcyr(40));
hsv_40 = hsv(40);
cmap = colormap([bcyr_40(20:(-1):14,:); hsv_40(23:1:40,:); bcyr_40(36:1:40,:)]);

fig1 = figure(1);
% h = plot(lon_bathy_along_transect,-z_along_transect,lon_endpoints_at_depth{1},depth_bounds(1)*[1 1],lon_endpoints_at_depth{length(lon_endpoints_at_depth)},depth_bounds(2)*[1 1]);
h = plot(lon_bathy_along_transect,-z_along_transect);
set(gca,'xlim',lon_bounds + (360*round((lon_bathy_along_transect(1) - lon_bounds(1))/360)),'ydir','reverse','clim',[0 (0.01*ceil((1e-6)*max_stddev_segment/0.01))])
curr_ytick = get(gca,'ytick');
curr_yticklabel = get(gca,'yticklabel');
curr_yticklabel{curr_ytick < 0} = '';
set(gca,'yticklabel',curr_yticklabel)
set(h(1),'Color',[0.6 0.3 0.2],'LineWidth',1)
% set(h(1 + (1:1:size(lon_endpoints_at_depth{1},1))),'Color',[0 0 0],'LineWidth',1.5)
% set(h((1 + size(lon_endpoints_at_depth{1},1)) + (1:1:size(lon_endpoints_at_depth{length(lon_endpoints_at_depth)},1))),'Color',[0 0 0],'LineWidth',1.5)
hold on
line(lon_bounds + (360*round((lon_bathy_along_transect(1) - lon_bounds(1))/360)),[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
for depth_ind = 1:length(depth_level_centers)
%     patch('XData',[(lon_endpoints_at_depth{depth_ind}'); flip(lon_endpoints_at_depth{depth_ind}',1)],'YData',repmat([(depth_level_bounds(depth_ind)*[1; 1]); (depth_level_bounds(depth_ind + 1)*[1; 1])],[1 size(lon_endpoints_at_depth{depth_ind},1)]),'EdgeColor',[0 0 0.8],'FaceColor','none')
    p = patch('XData',[(lon_endpoints_at_depth{depth_ind}'); flip(lon_endpoints_at_depth{depth_ind}',1)],'YData',repmat([(depth_level_bounds(depth_ind)*[1; 1]); (depth_level_bounds(depth_ind + 1)*[1; 1])],[1 size(lon_endpoints_at_depth{depth_ind},1)]),'FaceColor','flat','FaceVertexCData',(1e-6)*vol_flux_stddev_segments{depth_ind});
end
hold off
xlabel('Longitude')
ylabel('Depth (meters)')
% title(['Bathymetry profile in the ',lon_transect_id,' at ',num2str(lat_transect),' ^{\circ} lat, and integrated area for flux calculation'],'FontSize',10)
title({['Bathymetry profile in the ',lon_transect_id,' at ',num2str(lat_transect),' ^{o} lat, and standard deviation in segments']; ' '},'FontSize',10)
cbar = colorbar('southoutside');
set(cbar,'FontSize',10)
set(get(cbar,'xlabel'),'String','Volume flux standard deviation (Sv)')

% saveas(fig1,['Bathymetry_profile_',num2str(lat_transect),'_lat_transect_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(depth_bounds(1)),'_',num2str(depth_bounds(2)),'_depth_',datestr(min(time),'yyyymm'),'_',datestr(max(time),'yyyymm'),'_time_bounds.pdf'])
saveas(fig1,['Bathy_profile_std_dev_',num2str(lat_transect),'_lat_transect_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(depth_bounds(1)),'_',num2str(depth_bounds(2)),'_depth_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',downscale_id,'.pdf'])
close(fig1)


% fig2 = figure(2);
% h = plot((1e-6)*vol_flux_in_layers_tmean_perdepth,depth_level_centers,(1e-6)*(vol_flux_in_layers_tmean_perdepth - vol_flux_in_layers_tmean_uncert_perdepth),depth_level_centers,(1e-6)*(vol_flux_in_layers_tmean_perdepth + vol_flux_in_layers_tmean_uncert_perdepth),depth_level_centers);
% set(gca,'ydir','reverse')
% set(h(1),'Color',[0 0 0.8],'LineWidth',2)
% set(h(2),'Color',[0.5 0.5 0.5],'LineWidth',1)
% set(h(3),'Color',[0.5 0.5 0.5],'LineWidth',1)
% hold on
% line([0 0],depth_bounds,[0 0],'Color',[0 0 0],'LineWidth',1)
% hold off
% xlabel('Volume flux per meter depth (Sv m^{-1})')
% ylabel('Depth (meters)')
% title({['Time mean volume flux per meter depth in the ',lon_transect_id,' at ',num2str(lat_transect),' ^{o} lat, with standard error bounds']; ' '},'FontSize',10)
% 
% saveas(fig2,['Vol_flux_depth_profile_tmean_',num2str(lat_transect),'_lat_transect_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(depth_bounds(1)),'_',num2str(depth_bounds(2)),'_depth_',datestr(min(time),'yyyymm'),'_',datestr(max(time),'yyyymm'),'_time_bounds.pdf'])
% close(fig2)




% remove linear trend of obs. time series (and bandpass if desired)

curr_obs_time = RAPID_time;
curr_obs_tseries = RAPID_transport_1100_3000;


good_ind = find(isnan(curr_obs_tseries) == 0);
[~,unique_good_ind,~] = unique(curr_obs_tseries(good_ind));
good_ind = good_ind(unique_good_ind);

if abs(season_cyc_opt) < 1e-5
    % remove annual cycle

    G = [ones([length(good_ind) 1]) cos(((2*pi)/365.2425)*curr_obs_time(good_ind)) sin(((2*pi)/365.2425)*curr_obs_time(good_ind)) cos(((2*(2*pi))/365.2425)*curr_obs_time(good_ind)) sin(((2*(2*pi))/365.2425)*curr_obs_time(good_ind)) cos(((3*(2*pi))/365.2425)*curr_obs_time(good_ind)) sin(((3*(2*pi))/365.2425)*curr_obs_time(good_ind)) cos(((4*(2*pi))/365.2425)*curr_obs_time(good_ind)) sin(((4*(2*pi))/365.2425)*curr_obs_time(good_ind))];
    coeffs = (((G')*G)^(-1))*((G')*curr_obs_tseries(good_ind));
    curr_obs_tseries(good_ind) = curr_obs_tseries(good_ind) - (G(:,2:size(G,2))*coeffs(2:size(G,2)));
end

curr_tseries_nan_mask = ones(size(curr_obs_tseries));
curr_tseries_nan_mask((isnan(curr_obs_tseries) == 1) | (abs(curr_obs_tseries) < 1e-10)) = 0;
% interpolate to NaN data points and clip NaN ends
in_obs_range_ind = (min(good_ind):1:max(good_ind))';
[~,unique_obs_range_ind] = unique(curr_obs_time(in_obs_range_ind));
in_obs_range_ind = in_obs_range_ind(unique_obs_range_ind);
curr_obs_tseries = interp1(curr_obs_time(good_ind),curr_obs_tseries(good_ind),curr_obs_time(in_obs_range_ind));
curr_obs_time = curr_obs_time(in_obs_range_ind);
curr_tseries_nan_mask = curr_tseries_nan_mask(in_obs_range_ind);

% fill in gaps (if any exist)
gap_inds = find(diff(curr_obs_time) > 1.5*mean(diff(curr_obs_time)));
for gap_ind = 1:length(gap_inds)
    curr_gap_ind = gap_inds(gap_ind);
    gap_length = round(diff(curr_obs_time(curr_gap_ind + [0 1]))/(mean(diff(curr_obs_time)))) - 1;
    curr_obs_time = [curr_obs_time(1:curr_gap_ind); (curr_obs_time(curr_gap_ind) + ((mean(diff(curr_obs_time)))*(1:1:gap_length)')); curr_obs_time((curr_gap_ind + 1):length(curr_obs_time))];
    curr_obs_tseries = [curr_obs_tseries(1:curr_gap_ind); NaN([gap_length 1]); curr_obs_tseries((curr_gap_ind + 1):length(curr_obs_tseries))];
    curr_tseries_nan_mask = [curr_tseries_nan_mask(1:curr_gap_ind); zeros([gap_length 1]); curr_tseries_nan_mask((curr_gap_ind + 1):length(curr_tseries_nan_mask))];
end
[~,unique_ind] = unique(curr_obs_time);
curr_obs_time = curr_obs_time(unique_ind);
curr_obs_tseries = curr_obs_tseries(unique_ind);
curr_tseries_nan_mask = curr_tseries_nan_mask(unique_ind);
curr_obs_tseries = interp1(curr_obs_time(isnan(curr_obs_tseries) == 0),curr_obs_tseries(isnan(curr_obs_tseries) == 0),curr_obs_time);

[curr_obs_tseries_filtered,~,~] = bandpass_err_fcn(curr_obs_tseries,1,mean(diff(curr_obs_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,edge_handling_opt,1);
filter_gain_coeffs_array = squeeze(bandpass_err_fcn_gain_coeffs(curr_obs_tseries,1,mean(diff(curr_obs_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1));


% estimate uncertainty due to missing obs., and propagate through time filter

sum_nan_mask = sum(curr_tseries_nan_mask,1);
% [~,curr_tseries_dof_zero_lag,~,~,~,~,curr_tseries_std_dev,~,~] = regression_linear_scalar_scalar(curr_tseries_filtered,curr_tseries_filtered,1,mean(diff(curr_tseries_time)),mean(diff(curr_tseries_time)),(1/5)*length(curr_tseries_time)*mean(diff(curr_tseries_time))*[0 1],0.95,0);
[~,curr_obs_tseries_dof_zero_lag,~,~,~,~,curr_obs_tseries_std_dev,~,lags] = regression_linear_scalar_scalar(curr_obs_tseries_filtered,curr_obs_tseries_filtered,1,mean(diff(curr_obs_time)),mean(diff(time)),(1/5)*length(curr_obs_time)*mean(diff(curr_obs_time))*[0 1],0.95,0);
curr_obs_tseries_std_dev = curr_obs_tseries_std_dev(abs(lags) < 1e-5);
curr_obs_tseries_decorr_timescale = (mean(diff(curr_obs_time))*sum_nan_mask)./curr_obs_tseries_dof_zero_lag;
curr_tseries_missing_obs_err = NaN(size(curr_obs_tseries_filtered));
for t = 1:length(curr_tseries_nan_mask)
    curr_in_range_ind = find(abs(filter_gain_coeffs_array(t,:)) >= 0.01*max(abs(filter_gain_coeffs_array(t,:))));
    curr_in_range_ind = (min(curr_in_range_ind):1:max(curr_in_range_ind))';

    curr_tseries_obs_separation_array = abs(repmat(curr_obs_time(curr_in_range_ind),[1 length(curr_in_range_ind)]) - repmat(curr_obs_time(curr_in_range_ind)',[length(curr_in_range_ind) 1]));
    curr_tseries_missing_obs_err_cov = (curr_obs_tseries_std_dev.^2)*(1 - (curr_tseries_obs_separation_array/curr_obs_tseries_decorr_timescale));
    curr_tseries_missing_obs_err_cov(curr_tseries_missing_obs_err_cov < 0) = 0;

    curr_obs_tseries_err_cross_mask = (1 - repmat(curr_tseries_nan_mask(curr_in_range_ind),[1 length(curr_in_range_ind)])).*(1 - repmat(curr_tseries_nan_mask(curr_in_range_ind)',[length(curr_in_range_ind) 1]));
    curr_tseries_missing_obs_err_cov = curr_obs_tseries_err_cross_mask.*curr_tseries_missing_obs_err_cov;
    curr_cov_with_filter_coeffs = repmat(filter_gain_coeffs_array(t,curr_in_range_ind)',[1 length(curr_in_range_ind)]).*repmat(filter_gain_coeffs_array(t,curr_in_range_ind),[length(curr_in_range_ind) 1]).*curr_tseries_missing_obs_err_cov;
    curr_tseries_missing_obs_err(t) = (sum(sum(curr_cov_with_filter_coeffs,2),1))^(1/2);
end

curr_obs_tseries = curr_obs_tseries_filtered;
clear curr_obs_tseries_filtered

curr_tseries_missing_obs_norm_err = curr_tseries_missing_obs_err./repmat(curr_obs_tseries_std_dev,[length(curr_tseries_missing_obs_err) 1]);

curr_obs_tseries(abs(curr_tseries_missing_obs_norm_err) > norm_err_tolerance) = NaN;


% % monthly-average obs. time series prior to interpolation
% 
% obs_time_datevec = datevec(curr_obs_time);
% obs_month_num = (12*(obs_time_datevec(:,1) - min(obs_time_datevec(:,1)))) + obs_time_datevec(:,2);
% n_months_obs = max(obs_month_num) - min(obs_month_num) + 1;
% 
% obs_time_monthavg = NaN([n_months_obs 1]);
% obs_tseries_monthavg = NaN([n_months_obs 1]);
% for curr_month_num = min(obs_month_num):1:max(obs_month_num)
%     in_month_ind = find(abs(obs_month_num - curr_month_num) < 1e-3);
% 
%     if length(find(isnan(curr_obs_tseries(in_month_ind)) == 0)) > (0.8*28)/mean(diff(curr_obs_time))
%         obs_time_monthavg(curr_month_num - min(obs_month_num) + 1) = mean(curr_obs_time(in_month_ind(isnan(curr_obs_tseries(in_month_ind)) == 0)));
%         obs_tseries_monthavg(curr_month_num - min(obs_month_num) + 1) = mean(curr_obs_tseries(in_month_ind(isnan(curr_obs_tseries(in_month_ind)) == 0)));
%     end
% end
% 
% curr_obs_time = obs_time_monthavg;
% curr_obs_tseries = obs_tseries_monthavg;


curr_obs_tseries_tinterp = interp1(curr_obs_time,curr_obs_tseries,time);


[GRACE_RAPID_tseries_corr,~,~,~,GRACE_RAPID_tseries_lowmag_bound,GRACE_RAPID_tseries_highmag_bound,lags_tseries] = correlation_scalar_scalar_uncert_bounds(vol_flux_total',curr_obs_tseries_tinterp,1,mean(diff(time)),365.24*(1/12),365.24*[0 1],0.95);
GRACE_RAPID_tseries_corr = GRACE_RAPID_tseries_corr(abs(lags_tseries) < 1e-5);
GRACE_RAPID_tseries_lowmag_bound = GRACE_RAPID_tseries_lowmag_bound(abs(lags_tseries) < 1e-5);
GRACE_RAPID_tseries_highmag_bound = GRACE_RAPID_tseries_highmag_bound(abs(lags_tseries) < 1e-5);


fig3 = figure(3);
% h = plot(time,(1e-6)*vol_flux_total);
h = plot(time,(1e-6)*vol_flux_total,curr_obs_time,curr_obs_tseries);
tick_year_spacing = 2;
years_to_plot_ticks = ((ceil(time_range_start(1)/tick_year_spacing)*tick_year_spacing):tick_year_spacing:(ceil((time_range_end(1) - 1)/tick_year_spacing)*tick_year_spacing))';
xtick_datenums_plot = datenum([years_to_plot_ticks ones(length(years_to_plot_ticks),2)]);
xtick_labels_plot = cell(length(years_to_plot_ticks),1);
for xtick_ind = 1:length(years_to_plot_ticks)
    xtick_labels_plot{xtick_ind} = num2str(years_to_plot_ticks(xtick_ind));
end
set(gca,'xlim',[datenum(time_range_start) datenum(time_range_end)],'xtick',xtick_datenums_plot,'xticklabel',xtick_labels_plot)
set(h(1),'Color',[0 0 0.8],'LineWidth',2)
set(h(2),'Color',[0.8 0 0],'LineWidth',2)
hold on
line([datenum(time_range_start) datenum(time_range_end)],[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
% line(time,(1e-6)*(vol_flux_total - vol_flux_uncert_total),zeros([length(time) 1]),'Color',[0.5 0.5 0.5],'LineWidth',1)
% line(time,(1e-6)*(vol_flux_total + vol_flux_uncert_total),zeros([length(time) 1]),'Color',[0.5 0.5 0.5],'LineWidth',1)
hold off
ylabel('Volume flux anomaly (Sv)')
title({['Volume flux anomaly in the ',lon_transect_id,' at ',num2str(lat_transect),' ^{o} lat between ',num2str(depth_bounds(1)),' and ',num2str(depth_bounds(2)),' m depth']; ['Correlation coefficient: ',num2str(GRACE_RAPID_tseries_corr),', 95% CI: ',num2str(GRACE_RAPID_tseries_lowmag_bound),' to ',num2str(GRACE_RAPID_tseries_highmag_bound),'.']; ' '},'FontSize',10)
leg = legend('GRACE','RAPID','location','northeast');
keyboard
% saveas(fig3,['Vol_flux_anom_time_series_GRACE_',num2str(lat_transect),'_lat_transect_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(depth_bounds(1)),'_',num2str(depth_bounds(2)),'_depth_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',num2str(local_radius),'_horizrad_',num2str(depth_radius),'_depthrad_',datestr(min(time),'yyyymm'),'_',datestr(max(time),'yyyymm'),'_time_bounds.pdf'])
saveas(fig3,['Vol_flux_anom_time_series_GRACE_RAPID_',num2str(lat_transect),'_lat_transect_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(depth_bounds(1)),'_',num2str(depth_bounds(2)),'_depth_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',downscale_id,'.pdf'])
close(fig3)


% save .mat file (change file name and variable names as needed)

% save('GRACE_RAPID_tseries_method_comparisons.mat','time','*time','vol_flux*')
save('GRACE_RAPID_tseries_method_comparisons.mat','time','*time','vol_flux*','-append')