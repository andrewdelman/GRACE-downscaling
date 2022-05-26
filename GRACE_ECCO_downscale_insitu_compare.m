% downscale OBP from GRACE data, using covariances from ECCO with adjustments, and compare with available DART time series

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
cd('/indopac/adelman/GRACE/')


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

% % normalized noise parameter
% noise_param = 0;

% normalized error tolerance
norm_err_tolerance = 0.3;


% map bounds
longwest = -340;
longeast = 20;
latsouth = -90;
latnorth = 90;
projection = 'robinson';

DART_filenames = dir('/indopac/adelman/GRACE/DART/DART_obp_*_all.nc');
filenames_cellarray = struct2cell(DART_filenames);
filenames_cellarray = filenames_cellarray(1,:);


time_range_start = [2002 4 1];
time_range_end = [2017 7 1];

% mascon_lat_separation = 3;    % in degrees

curr_file = 'LAND_MASK.CRIv01.nc';
land_mask = ncread(curr_file,'land_mask');
curr_file = 'GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';
lon = ncread(curr_file,'lon');
lat = ncread(curr_file,'lat');
time = ncread(curr_file,'time') + datenum([2002 1 1 0 0 0]);
in_time_range_ind = find((time >= datenum(time_range_start)) & (time < datenum(time_range_end)));
time_range_ind_span = max(in_time_range_ind) - min(in_time_range_ind) + 1;
time = time(in_time_range_ind);
lwe_thickness = ncread(curr_file,'lwe_thickness',[1 1 min(in_time_range_ind)],[length(lon) length(lat) time_range_ind_span]);
lwe_uncert = ncread(curr_file,'uncertainty',[1 1 min(in_time_range_ind)],[length(lon) length(lat) time_range_ind_span]);


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


% add path to DART filenames
for curr_file = 1:length(filenames_cellarray)
    filenames_cellarray{curr_file} = ['/indopac/adelman/GRACE/DART/',filenames_cellarray{curr_file}];
end


% add non-DART data to list

n_pts = length(filenames_cellarray);
end_dataset_ind = n_pts;

% SAMOC PIES
for curr_ind = 1:1:4
    filenames_cellarray{n_pts + curr_ind} = '/indopac/adelman/GRACE/SAMOC/SAM_PIES_data.txt';
end
n_pts = n_pts + 4;
end_dataset_ind = [end_dataset_ind; n_pts];

% ANT PIES
filenames_cellarray{n_pts + 1} = '/indopac/adelman/GRACE/ANT/PIES-3.3-prs-anomaly-daily.dat';
filenames_cellarray{n_pts + 2} = '/indopac/adelman/GRACE/ANT/PIES-17-1-prs-anomaly-daily.dat';
n_pts = n_pts + 2;
end_dataset_ind = [end_dataset_ind; n_pts];

% NOAC moorings and PIES
filenames_cellarray{n_pts + 1} = '/indopac/adelman/GRACE/NOAC/WAtl/datasets/PIES_BP27-1_pressure_TWT.tab';
filenames_cellarray{n_pts + 2} = '/indopac/adelman/GRACE/NOAC/WAtl/datasets/PIES_BP29-1_pressure_TWT.tab';
n_pts = n_pts + 2;
end_dataset_ind = [end_dataset_ind; n_pts];

% RAPID data
RAPID_filenames = dir('/indopac/adelman/GRACE/RAPID/RAPID_obp_first50removed*.nc');
RAPID_filenames_cellarray = struct2cell(RAPID_filenames);
RAPID_filenames_cellarray = RAPID_filenames_cellarray(1,:);
for curr_file = 1:length(RAPID_filenames)
    RAPID_filenames_cellarray{curr_file} = ['/indopac/adelman/GRACE/RAPID/',RAPID_filenames_cellarray{curr_file}];
end

filenames_cellarray = [filenames_cellarray RAPID_filenames_cellarray];
n_pts = n_pts + length(RAPID_filenames_cellarray);
end_dataset_ind = [end_dataset_ind; n_pts];


obs_lon = NaN([n_pts 1]);
obs_lat = NaN([n_pts 1]);
obs_depth = NaN([n_pts 1]);
obs_n_days_equiv = NaN([n_pts 1]);
obs_stddev_unfiltered = NaN([n_pts 1]);
obs_stddev_filtered = NaN([n_pts 1]);
ECCO_stddev_unfiltered = NaN([n_pts 1]);
ECCO_stddev_filtered = NaN([n_pts 1]);
obs_tseries_reconstr_corrcoeff = NaN([n_pts 1]);
obs_tseries_reconstr_lowmag_bound = NaN([n_pts 1]);
obs_tseries_reconstr_highmag_bound = NaN([n_pts 1]);
for curr_pt = 1:n_pts
    
    % load bottom pressure observations
    if curr_pt <= end_dataset_ind(1)
        curr_filename = filenames_cellarray{curr_pt};
        
        lon_pt_obs = ncread(curr_filename,'longitude');
        lat_pt_obs = ncread(curr_filename,'latitude');
        obs_lon(curr_pt) = ncread(curr_filename,'longitude');
        obs_lat(curr_pt) = ncread(curr_filename,'latitude');
        time_obs = ncread(curr_filename,'time') + 0.5;
        n_samples = ncread(curr_filename,'n_samples');
        obs_n_days_equiv(curr_pt) = (1/(4*24))*sum(n_samples);
        
        OBP_obs = (1e2)*ncread(curr_filename,'obp_nodrift');
        OBP_obs(OBP_obs < -99000) = NaN;
        frac_good_obs = double((1/96)*n_samples);
        frac_good_obs(isnan(OBP_obs) == 1) = 0;
    elseif curr_pt <= end_dataset_ind(2)
        curr_file = filenames_cellarray{curr_pt};
        curr_fid = fopen(curr_file);
        
        fseek(curr_fid,931,'bof');
        data_array = (fscanf(curr_fid,'%d %d %d %f %f %f %f %f %f %f %f',[11 Inf]))';
        time_obs = datenum(data_array(:,1:3)) + 0.5;
        lat_pt_obs = -34.5;
        if curr_pt == end_dataset_ind(1) + 1
            lon_pt_obs = -51.5;
            OBP_obs = (1e4)*data_array(:,5);
        elseif curr_pt == end_dataset_ind(1) + 2
            lon_pt_obs = -49.5;
            OBP_obs = (1e4)*data_array(:,7);
        elseif curr_pt == end_dataset_ind(1) + 3
            lon_pt_obs = -47.5;
            OBP_obs = (1e4)*data_array(:,9);
        elseif curr_pt == end_dataset_ind(1) + 4
            lon_pt_obs = -44.5;
            OBP_obs = (1e4)*data_array(:,11);
        end
        obs_lon(curr_pt) = lon_pt_obs;
        obs_lat(curr_pt) = lat_pt_obs;
        frac_good_obs = ones(size(OBP_obs));
        frac_good_obs(isnan(OBP_obs) == 1) = 0;
        obs_n_days_equiv(curr_pt) = sum(1*frac_good_obs);
        fclose(curr_fid);
    elseif curr_pt <= end_dataset_ind(3)
        curr_file = filenames_cellarray{curr_pt};
        curr_fid = fopen(curr_file);
        
        fseek(curr_fid,0,'bof');
        data_array = (fscanf(curr_fid,'%f %f',[2 Inf]))';
        if curr_pt == end_dataset_ind(2) + 1
            lon_pt_obs = 12 + (45.23/60);
            lat_pt_obs = -(37 + (05.84/60));
            time_obs = data_array(:,1) + datenum([2009 12 31 09 46 16]);
        elseif curr_pt == end_dataset_ind(2) + 2
            lon_pt_obs = -(0 + (02.72/60));
            lat_pt_obs = -(64 + (00.70/60));
            time_obs = data_array(:,1) + datenum([2009 12 31 14 56 01]);
        end
        OBP_obs = (9.81*1000)*data_array(:,2);
        obs_lon(curr_pt) = lon_pt_obs;
        obs_lat(curr_pt) = lat_pt_obs;
        frac_good_obs = ones(size(OBP_obs));
        frac_good_obs(isnan(OBP_obs) == 1) = 0;
        obs_n_days_equiv(curr_pt) = sum(1*frac_good_obs);
        fclose(curr_fid);
    elseif curr_pt <= end_dataset_ind(4)
        curr_file = filenames_cellarray{curr_pt};
        curr_fid = fopen(curr_file);
        
        if curr_pt == end_dataset_ind(3) + 1
            fseek(curr_fid,1725,'bof');
            lon_pt_obs = -40.8755;
            lat_pt_obs = 47.0973;
            time_obs = [];
            OBP_obs = [];
            curr_line = fgetl(curr_fid);
            while curr_line ~= -1
                time_obs = [time_obs; datenum([str2double(curr_line(1:4)) str2double(curr_line(6:7)) str2double(curr_line(9:10)) str2double(curr_line(12:13)) str2double(curr_line(15:16)) 00])];
                OBP_obs = [OBP_obs; ((1e4)*str2double(curr_line(18:25)))];
                curr_line = fgetl(curr_fid);
            end
            OBP_obs(OBP_obs < 4.582e7) = NaN;
        elseif curr_pt == end_dataset_ind(3) + 2
            fseek(curr_fid,1725,'bof');
            lon_pt_obs = -38.5182;
            lat_pt_obs = 47.2087;
            time_obs = [];
            OBP_obs = [];
            curr_line = fgetl(curr_fid);
            while curr_line ~= -1
                time_obs = [time_obs; datenum([str2double(curr_line(1:4)) str2double(curr_line(6:7)) str2double(curr_line(9:10)) str2double(curr_line(12:13)) str2double(curr_line(15:16)) 00])];
                OBP_obs = [OBP_obs; ((1e4)*str2double(curr_line(18:25)))];
                curr_line = fgetl(curr_fid);
            end
            OBP_obs(OBP_obs < 4.714e7) = NaN;
        end
        obs_lon(curr_pt) = lon_pt_obs;
        obs_lat(curr_pt) = lat_pt_obs;
        frac_good_obs = ones(size(OBP_obs));
        frac_good_obs(isnan(OBP_obs) == 1) = 0;
        obs_n_days_equiv(curr_pt) = sum(1*frac_good_obs);
        fclose(curr_fid);
    elseif curr_pt <= end_dataset_ind(5)
        curr_file = filenames_cellarray{curr_pt};
        
        lon_str_ind = strfind(curr_file,'lon');
        lat_str_ind = strfind(curr_file,'lat');
%         lon_pt_obs = str2double(curr_file(40:(lon_str_ind - 2)));
        lon_pt_obs = str2double(curr_file(55:(lon_str_ind - 2)));
        lat_pt_obs = str2double(curr_file((lon_str_ind + 4):(lat_str_ind - 2)));
        time_obs = ncread(curr_file,'time');
        OBP_obs = (1e4)*ncread(curr_file,'obp_nodrift');
    %     OBP_obs(OBP_obs < -99000) = NaN;
        frac_good_obs = ones(size(OBP_obs));
        frac_good_obs(isnan(OBP_obs) == 1) = 0;

        good_obs_ind = find(frac_good_obs > 0.5);
        [sorted_diff,~] = sort(diff(sort(time_obs(good_obs_ind),'ascend')),'ascend');
        median_diff = sorted_diff(ceil(length(good_obs_ind)/2));
        if median_diff < 0.4
            time_obs_daily = unique(floor(time_obs(good_obs_ind))) + 0.5;
            time_obs_daily = (min(time_obs_daily):1:max(time_obs_daily))';
            curr_dt = mode((1e3)*((round((1e3)*((diff(time_obs(isnan(time_obs) == 0))).^(-1)))).^(-1)));
            [OBP_obs_filtered,~,~] = bandpass_err_fcn(OBP_obs,1,curr_dt,1/(4*length(time_obs)*curr_dt),1/2,steepness_factor,1,1,0,0);
            OBP_obs_daily = NaN(size(time_obs_daily));
            frac_good_obs_daily = NaN(size(time_obs_daily));
            for curr_t_ind = 1:length(time_obs_daily)
                at_curr_t_ind = find(abs(floor(time_obs) + 0.5 - time_obs_daily(curr_t_ind)) < 1e-5);
                at_curr_t_good_ind = at_curr_t_ind(isnan(OBP_obs(at_curr_t_ind)) == 0);

                OBP_obs_daily(curr_t_ind) = mean(OBP_obs_filtered(at_curr_t_good_ind));
                if isempty(at_curr_t_good_ind) == 1
                    frac_good_obs_daily(curr_t_ind) = 0;
                else
%                     frac_good_obs_daily(curr_t_ind) = 0.01*round((length(at_curr_t_good_ind)/(1/(mean(diff(time_obs(at_curr_t_ind))))))/0.01);
                    frac_good_obs_daily(curr_t_ind) = max([min([1 (length(at_curr_t_good_ind)/24)]) (0.01*round((length(at_curr_t_good_ind)/(1/(mean(diff(time_obs(at_curr_t_ind))))))/0.01))]);
                end
            end
            time_obs = time_obs_daily;
            OBP_obs = OBP_obs_daily;
            frac_good_obs = frac_good_obs_daily;
            obs_n_days_equiv(curr_pt) = sum(1*frac_good_obs_daily);
        else
            obs_n_days_equiv(curr_pt) = sum(median_diff*frac_good_obs);
        end
        obs_lon(curr_pt) = lon_pt_obs;
        obs_lat(curr_pt) = lat_pt_obs;
    end
    
    
    % calculations with bathymetry in region surrounding obs. point
    
    in_local_bathy_lat_range_ind = find(abs(lat_bathy - lat_pt_obs) <= radius_mascons_deg + 5);
    local_lat_bathy = lat_bathy(in_local_bathy_lat_range_ind);
    if mod(lon_pt_obs + radius_mascons_deg + 5 - max(lon_bathy),360) < mod(lon_pt_obs - radius_mascons_deg - 5 - min(lon_bathy),360)
        in_local_bathy_lon_range_ind_1 = find(lon_bathy - (mod(lon_pt_obs - radius_mascons_deg - 5 - (min(lon_bathy) - 1e-5),360) + (min(lon_bathy) - 1e-5)) >= 0);
        in_local_bathy_lon_range_ind_2 = find(lon_bathy - (mod(lon_pt_obs + radius_mascons_deg + 5 - (min(lon_bathy) - 1e-5),360) + (min(lon_bathy) - 1e-5)) <= 0);
        local_lon_bathy = mod([lon_bathy(in_local_bathy_lon_range_ind_1); (lon_bathy(in_local_bathy_lon_range_ind_2) + 360)] - (lon_pt_obs - 180),360) + (lon_pt_obs - 180);
        z_local = z_reshaped([in_local_bathy_lon_range_ind_1; in_local_bathy_lon_range_ind_2],in_local_bathy_lat_range_ind);
    else
        in_local_bathy_lon_range_ind = find(abs(mod(lon_bathy - lon_pt_obs + 180,360) - 180) <= radius_mascons_deg + 5);
        local_lon_bathy = mod(lon_bathy(in_local_bathy_lon_range_ind) - (lon_pt_obs - 180),360) + (lon_pt_obs - 180);
        z_local = z_reshaped(in_local_bathy_lon_range_ind,in_local_bathy_lat_range_ind);
    end
    
    
    bathy_near_obs_lon_ind = find(abs(mod(lon_bathy - lon_pt_obs + 180,360) - 180) <= 5);
    gap_ind = find(diff(bathy_near_obs_lon_ind) > 1.5);
    if isempty(gap_ind) == 0
        bathy_near_obs_lon_ind = bathy_near_obs_lon_ind([((gap_ind + 1):1:length(bathy_near_obs_lon_ind))'; (1:1:gap_ind)']);
    end
    bathy_near_obs_lat_ind = find(abs(lat_bathy - lat_pt_obs) <= 5);
    depth_pt_obs = -interp2_fft(mod(lon_bathy(bathy_near_obs_lon_ind) - lon_pt_obs + 180,360) + lon_pt_obs - 180,lat_bathy(bathy_near_obs_lat_ind),z_reshaped(bathy_near_obs_lon_ind,bathy_near_obs_lat_ind),lon_pt_obs,lat_pt_obs);
    
    obs_depth(curr_pt) = depth_pt_obs;
    
    
    % planar fit of local bathymetry within certain radius
    bathy_radius = 6;    % in degrees latitude
    range_depth_threshold = 2000;    % threshold for inclusion in analysis (to include only points near large bathymetry slopes)

    bathy_near_pt_ind = find(abs(((cosd(lat_pt_obs))*(mod(repmat(local_lon_bathy,[1 length(local_lat_bathy)]) - lon_pt_obs + 180,360) - 180)) + (1i*(repmat(local_lat_bathy',[length(local_lon_bathy) 1]) - lat_pt_obs))) <= bathy_radius);
    lon_bathy_near_pt = local_lon_bathy(mod(bathy_near_pt_ind - 1,length(local_lon_bathy)) + 1);
    lat_bathy_near_pt = local_lat_bathy(ceil(bathy_near_pt_ind/length(local_lon_bathy)));

    G = [ones([length(bathy_near_pt_ind) 1]) lon_bathy_near_pt lat_bathy_near_pt];
    m_depth_planar_fit = (((G')*G)^(-1))*((G')*z_local(bathy_near_pt_ind));
    local_depth_planar_fit = G*m_depth_planar_fit;

    mean_depth_within_radius = -mean(z_local(bathy_near_pt_ind));
    range_depth_planar_fit = abs(max(local_depth_planar_fit) - min(local_depth_planar_fit));
    if ((range_depth_planar_fit < range_depth_threshold) || (depth_pt_obs < mean_depth_within_radius) || (obs_n_days_equiv(curr_pt) < 730))
        continue
    else
        disp(['lon_pt_obs = ',num2str(lon_pt_obs),', lat_pt_obs = ',num2str(lat_pt_obs)])
    end
    angle_steepest_descent = imag(log(-m_depth_planar_fit(2) - (1i*m_depth_planar_fit(3))));



    % find mascons near obs. point

    in_radius_ind = find(abs((cosd(mean([mascon_lat_center_all (lat_pt_obs*ones([length(mascon_lat_center_all) 1]))],2)).*(mod(mascon_lon_center_all - lon_pt_obs + 180,360) - 180)) + (1i*(mascon_lat_center_all - lat_pt_obs))) < radius_mascons_deg);

    ocean_mascon_curr_values = ocean_mascon_curr_values_all(in_radius_ind);
    mascon_lon_center = mascon_lon_center_all(in_radius_ind);
    mascon_lat_center = mascon_lat_center_all(in_radius_ind);
    mascon_lon_bounds = mascon_lon_bounds_all(in_radius_ind,:);
    mascon_lat_bounds = mascon_lat_bounds_all(in_radius_ind,:);
    lwe_thickness_ocean_mascons = lwe_thickness_ocean_mascons_all(in_radius_ind,:);
    lwe_uncert_ocean_mascons = lwe_uncert_ocean_mascons_all(in_radius_ind,:);
    mascon_depth_avg = mascon_depth_avg_all(in_radius_ind);
    
    
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



    % remove linear trend of obs. (and bandpass if desired)

    curr_tseries_time = time_obs;
    curr_tseries = OBP_obs;


    good_ind = find(isnan(curr_tseries) == 0);
    [~,unique_good_ind,~] = unique(curr_tseries_time(good_ind));
    good_ind = good_ind(unique_good_ind);

    if abs(season_cyc_opt) < 1e-5
        % remove annual cycle

    %     curr_nan_ind = find(isnan(curr_tseries) == 1);
    %     curr_tseries(curr_nan_ind) = 0;

        G = [ones([length(good_ind) 1]) cos(((2*pi)/365.2425)*curr_tseries_time(good_ind)) sin(((2*pi)/365.2425)*curr_tseries_time(good_ind)) cos(((2*(2*pi))/365.2425)*curr_tseries_time(good_ind)) sin(((2*(2*pi))/365.2425)*curr_tseries_time(good_ind)) cos(((3*(2*pi))/365.2425)*curr_tseries_time(good_ind)) sin(((3*(2*pi))/365.2425)*curr_tseries_time(good_ind)) cos(((4*(2*pi))/365.2425)*curr_tseries_time(good_ind)) sin(((4*(2*pi))/365.2425)*curr_tseries_time(good_ind))];
        coeffs = (((G')*G)^(-1))*((G')*curr_tseries(good_ind));
        curr_tseries(good_ind) = curr_tseries(good_ind) - (G(:,2:size(G,2))*coeffs(2:size(G,2)));
    end

    % curr_tseries_nan_mask = ones(size(curr_tseries));
    % curr_tseries_nan_mask((isnan(curr_tseries) == 1) | (abs(curr_tseries) < 1e-10)) = 0;
    curr_tseries_nan_mask = frac_good_obs;
    curr_tseries_nan_mask(curr_tseries_nan_mask > 1) = 1;
    % interpolate to NaN data points and clip NaN ends
    in_obs_range_ind = (min(good_ind):1:max(good_ind))';
    [~,unique_obs_range_ind] = unique(curr_tseries_time(in_obs_range_ind));
    in_obs_range_ind = in_obs_range_ind(unique_obs_range_ind);
    curr_tseries = interp1(curr_tseries_time(good_ind),curr_tseries(good_ind),curr_tseries_time(in_obs_range_ind));
    curr_tseries_time = curr_tseries_time(in_obs_range_ind);
    curr_tseries_nan_mask = curr_tseries_nan_mask(in_obs_range_ind);

    % fill in gaps (if any exist)
    gap_inds = find(diff(curr_tseries_time) > 1.5*mean(diff(curr_tseries_time)));
    for gap_ind = 1:length(gap_inds)
        curr_gap_ind = gap_inds(gap_ind);
        gap_length = round(diff(curr_tseries_time(curr_gap_ind + [0 1]))/(mean(diff(curr_tseries_time)))) - 1;
        curr_tseries_time = [curr_tseries_time(1:curr_gap_ind); (curr_tseries_time(curr_gap_ind) + ((mean(diff(curr_tseries_time)))*(1:1:gap_length)')); curr_tseries_time((curr_gap_ind + 1):length(curr_tseries_time))];
        curr_tseries = [curr_tseries(1:curr_gap_ind); NaN([gap_length 1]); curr_tseries((curr_gap_ind + 1):length(curr_tseries))];
        curr_tseries_nan_mask = [curr_tseries_nan_mask(1:curr_gap_ind); zeros([gap_length 1]); curr_tseries_nan_mask((curr_gap_ind + 1):length(curr_tseries_nan_mask))];
    end
    [~,unique_ind] = unique(curr_tseries_time);
    curr_tseries_time = curr_tseries_time(unique_ind);
    curr_tseries = curr_tseries(unique_ind);
    curr_tseries_nan_mask = curr_tseries_nan_mask(unique_ind);
    curr_tseries = interp1(curr_tseries_time(isnan(curr_tseries) == 0),curr_tseries(isnan(curr_tseries) == 0),curr_tseries_time);

    [curr_tseries_filtered,~,~] = bandpass_err_fcn(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,edge_handling_opt,1);
    % [curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask + 1e-5,1,mean(diff(curr_tseries_time)),(1/(4*mean(diff(curr_tseries_time))*length(curr_tseries_time)))/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,1,1,edge_handling_opt,1);
    filter_gain_coeffs_array = squeeze(bandpass_err_fcn_gain_coeffs(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1));
    
    
    % estimate uncertainty due to missing obs., and propagate through time filter

    sum_nan_mask = sum(curr_tseries_nan_mask,1);
    % [~,curr_tseries_dof_zero_lag,~,~,~,~,curr_tseries_std_dev,~,~] = regression_linear_scalar_scalar(curr_tseries_filtered,curr_tseries_filtered,1,mean(diff(curr_tseries_time)),mean(diff(curr_tseries_time)),(1/5)*length(curr_tseries_time)*mean(diff(curr_tseries_time))*[0 1],0.95,0);
    [~,curr_tseries_dof_zero_lag,~,~,~,~,curr_tseries_std_dev,~,lags] = regression_linear_scalar_scalar(curr_tseries_filtered,curr_tseries_filtered,1,mean(diff(curr_tseries_time)),mean(diff(time_GRACE_interp)),(1/5)*length(curr_tseries_time)*mean(diff(curr_tseries_time))*[0 1],0.95,0);
    curr_tseries_std_dev = curr_tseries_std_dev(abs(lags) < 1e-5);
    curr_tseries_decorr_timescale = (mean(diff(curr_tseries_time))*sum_nan_mask)./curr_tseries_dof_zero_lag;
    curr_tseries_missing_obs_err = NaN(size(curr_tseries_filtered));
    midpt = ceil(size(filter_gain_coeffs_array,1)/2);
    if length(find(abs(filter_gain_coeffs_array(midpt,:)) >= 0.01*max(abs(filter_gain_coeffs_array(midpt,:))))) <= 0.2*size(filter_gain_coeffs_array,1)
        for t = 1:length(curr_tseries_nan_mask)
            curr_in_range_ind = find(abs(filter_gain_coeffs_array(t,:)) >= 0.01*max(abs(filter_gain_coeffs_array(t,:))));
            curr_in_range_ind = (min(curr_in_range_ind):1:max(curr_in_range_ind))';

            curr_tseries_obs_separation_array = abs(repmat(curr_tseries_time(curr_in_range_ind),[1 length(curr_in_range_ind)]) - repmat(curr_tseries_time(curr_in_range_ind)',[length(curr_in_range_ind) 1]));
            curr_tseries_missing_obs_err_cov = (curr_tseries_std_dev.^2)*(1 - (curr_tseries_obs_separation_array/curr_tseries_decorr_timescale));
            curr_tseries_missing_obs_err_cov(curr_tseries_missing_obs_err_cov < 0) = 0;

            curr_tseries_err_cross_mask = (1 - repmat(curr_tseries_nan_mask(curr_in_range_ind),[1 length(curr_in_range_ind)])).*(1 - repmat(curr_tseries_nan_mask(curr_in_range_ind)',[length(curr_in_range_ind) 1]));
            curr_tseries_missing_obs_err_cov = curr_tseries_err_cross_mask.*curr_tseries_missing_obs_err_cov;
            curr_cov_with_filter_coeffs = repmat(filter_gain_coeffs_array(t,curr_in_range_ind)',[1 length(curr_in_range_ind)]).*repmat(filter_gain_coeffs_array(t,curr_in_range_ind),[length(curr_in_range_ind) 1]).*curr_tseries_missing_obs_err_cov;
            curr_tseries_missing_obs_err(t) = (sum(sum(curr_cov_with_filter_coeffs,2),1))^(1/2);
        end
    else
        curr_tseries_obs_separation_array = abs(repmat(curr_tseries_time,[1 length(curr_tseries_time)]) - repmat(curr_tseries_time',[length(curr_tseries_time) 1]));
        curr_tseries_missing_obs_err_cov = (curr_tseries_std_dev.^2)*(1 - (curr_tseries_obs_separation_array/curr_tseries_decorr_timescale));
        curr_tseries_missing_obs_err_cov(curr_tseries_missing_obs_err_cov < 0) = 0;

        curr_tseries_err_cross_mask = (1 - repmat(curr_tseries_nan_mask,[1 length(curr_tseries_nan_mask)])).*(1 - repmat(curr_tseries_nan_mask',[length(curr_tseries_nan_mask) 1]));
        curr_tseries_missing_obs_err_cov = curr_tseries_err_cross_mask.*curr_tseries_missing_obs_err_cov;
        nonzero_err_ind = find(abs(curr_tseries_err_cross_mask) > 0.001);

        sum_curr_cov_with_filter_coeffs = zeros(size(curr_tseries_filtered));
        for curr_err_ind = 1:length(nonzero_err_ind)
            curr_e_ind = nonzero_err_ind(curr_err_ind);
            curr_i_ind = mod(curr_e_ind - 1,size(filter_gain_coeffs_array,1)) + 1;
            curr_j_ind = ceil(curr_e_ind/size(filter_gain_coeffs_array,1));
            sum_curr_cov_with_filter_coeffs = sum_curr_cov_with_filter_coeffs + ((filter_gain_coeffs_array(:,curr_i_ind).*filter_gain_coeffs_array(:,curr_j_ind))*curr_tseries_missing_obs_err_cov(curr_e_ind));
        end
        curr_tseries_missing_obs_err = (sum_curr_cov_with_filter_coeffs).^(1/2);
        
    end
    
    curr_tseries = curr_tseries_filtered;
    % curr_tseries_nan_mask = curr_tseries_nan_mask_filtered - 1e-5;
    clear curr_tseries_filtered curr_tseries_nan_mask_filtered

    % curr_tseries(curr_nan_ind) = NaN;

    curr_tseries_missing_obs_norm_err = curr_tseries_missing_obs_err./repmat(curr_tseries_std_dev,[length(curr_tseries_missing_obs_err) 1]);

    % curr_tseries(abs(curr_tseries_nan_mask) > 0.2) = NaN;
    curr_tseries(abs(curr_tseries_missing_obs_norm_err) > norm_err_tolerance) = NaN;
    
    
    % standard deviations of original and filtered time series
    [~,~,~,~,~,~,OBP_obs_stddev,~,~] = regression_linear_scalar_scalar(OBP_obs,OBP_obs,1,1,1,[0 1],0.95,0);
    [~,~,~,~,~,~,curr_tseries_stddev,~,~] = regression_linear_scalar_scalar(curr_tseries,curr_tseries,1,1,1,[0 1],0.95,0);
    obs_stddev_unfiltered(curr_pt) = OBP_obs_stddev(1);
    obs_stddev_filtered(curr_pt) = curr_tseries_stddev(1);
    
    obs_tseries_tinterp = interp1(curr_tseries_time,curr_tseries,time);



    if abs(downscale_opt - 1) < 1e-5

        % load ECCO OBP for covariance calculations

        ECCO_nc_file = '/indopac/adelman/ECCO2/PHIBOT.ECCO2.lonlatinterp.1992-2018.nc';

        longitude = ncread(ECCO_nc_file,'LONGITUDE_T');
        latitude = ncread(ECCO_nc_file,'LATITUDE_T');
        time_ECCO = ncread(ECCO_nc_file,'TIME');
        time_ECCO = double(time_ECCO) + datenum([1992 1 1 0 0 0]);

        in_lon_range_ind = find(abs(mod(longitude - lon_pt_obs + 180,360) - 180) <= radius_mascons_deg + 5);
        in_lat_range_ind = find(abs(latitude - lat_pt_obs) <= radius_mascons_deg + 5);
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
        lon_in_range_ECCO = lon_in_range_ECCO(1) + (360*floor((lon_pt_obs - lon_in_range_ECCO(1))/360)) + [0; cumsum(diff_lon_in_range_ECCO)];
        diff_lat_in_range_ECCO = diff(lat_in_range_ECCO);

        lon_in_range_ECCO_bounds = [(lon_in_range_ECCO(1) - (0.5*diff_lon_in_range_ECCO(1))); (lon_in_range_ECCO(2:length(lon_in_range_ECCO)) - (diff_lon_in_range_ECCO/2)); (lon_in_range_ECCO(length(lon_in_range_ECCO)) + (0.5*diff_lon_in_range_ECCO(length(diff_lon_in_range_ECCO))))];
        lat_in_range_ECCO_bounds = [(lat_in_range_ECCO(1) - (0.5*diff_lat_in_range_ECCO(1))); (lat_in_range_ECCO(2:length(lat_in_range_ECCO)) - (diff_lat_in_range_ECCO/2)); (lat_in_range_ECCO(length(lat_in_range_ECCO)) + (0.5*diff_lat_in_range_ECCO(length(diff_lat_in_range_ECCO))))];


        OBP_ECCO_zeronans = OBP_ECCO;
        OBP_ECCO_zeronans(isnan(OBP_ECCO) == 1) = 0;
        time_pt = time_in_range_ECCO;
        OBP_pt = squeeze(interp2_fft(lon_in_range_ECCO,lat_in_range_ECCO,OBP_ECCO_zeronans,lon_pt_obs,lat_pt_obs));
        
        [~,~,~,~,~,~,OBP_pt_stddev,~,~] = regression_linear_scalar_scalar(OBP_pt,OBP_pt,1,1,1,[0 1],0.95,1);
        
        
        clear OBP_ECCO_zeronans

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

            enough_ind = find(sum(nan_mask_ECCO_reshaped(:,in_month_ind),2) > (0.8*28)/mean(diff(time_pt)));
            time_ECCO_monthavg(curr_month_ind) = sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*repmat(reshape(time_in_range_ECCO(in_month_ind),[1 length(in_month_ind)]),[length(enough_ind) 1]),2),1)./(sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2),1));
            OBP_ECCO_monthavg(enough_ind,curr_month_ind) = sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*OBP_ECCO_reshaped_nans_zeroed(enough_ind,in_month_ind),2)./(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2));    
        end

        time_in_range_ECCO = time_ECCO_monthavg;
        OBP_ECCO = reshape(OBP_ECCO_monthavg,[size_array(1:2) n_months_ECCO]);


        time_pt_datevec = datevec(time_pt);
        n_months_pt = ((12*time_pt_datevec(size(time_pt_datevec,1),1)) + time_pt_datevec(size(time_pt_datevec,1),2)) - ((12*time_pt_datevec(1,1)) + time_pt_datevec(1,2)) + 1;
        nan_mask_pt = ones(size(OBP_pt));
        nan_mask_pt((isnan(OBP_pt) == 1) | (abs(OBP_pt) < 1e-10)) = 0;
        OBP_pt_nans_zeroed = OBP_pt;
        OBP_pt_nans_zeroed(abs(nan_mask_pt) < 1e-5) = 0;

        time_pt_monthavg = NaN([n_months_pt 1]);
        OBP_pt_monthavg = NaN([n_months_pt 1]);
        for curr_month_ind = 1:n_months_pt
            curr_month_num = (12*time_pt_datevec(1,1)) + time_pt_datevec(1,2) + curr_month_ind - 1;
            curr_yearmonth = [floor((curr_month_num - 1)/12) (mod(curr_month_num - 1,12) + 1)];
            in_month_ind = find((abs(time_pt_datevec(:,1) - curr_yearmonth(1)) < 1e-3) & (abs(time_pt_datevec(:,2) - curr_yearmonth(2)) < 1e-3));

            if sum(nan_mask_pt(in_month_ind)) > (0.8*28)/mean(diff(time_pt))
                time_pt_monthavg(curr_month_ind) = sum(nan_mask_pt(in_month_ind).*time_pt(in_month_ind),1)./(sum(nan_mask_pt(in_month_ind),1));
                OBP_pt_monthavg(curr_month_ind) = sum(nan_mask_pt(in_month_ind).*OBP_pt_nans_zeroed(in_month_ind),1)./(sum(nan_mask_pt(in_month_ind),1));
            end
        end

        time_pt = time_pt_monthavg(isnan(OBP_pt_monthavg) == 0);
        OBP_pt = OBP_pt_monthavg(isnan(OBP_pt_monthavg) == 0);

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


        curr_tseries_time = time_pt;
        curr_ECCO_tseries = OBP_pt;

        curr_tseries_nan_mask = (1e-5)*ones(size(curr_ECCO_tseries));
        curr_tseries_nan_mask((isnan(curr_ECCO_tseries) == 1) | (abs(curr_ECCO_tseries) < 1e-10)) = -1;
        [curr_tseries_filtered,curr_tseries_trend,~] = bandpass_err_fcn(curr_ECCO_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
        [curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask,1,mean(diff(curr_tseries_time)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

        curr_ECCO_tseries = curr_tseries_filtered;
        clear curr_tseries_filtered


        % nan_OBP_ECCO_ind = find(isnan(OBP_ECCO) == 1);

        nan_ECCO_ind = find(isnan(OBP_ECCO_ocean_mascons) == 1);
        nan_curr_tseries_ind = find(isnan(curr_ECCO_tseries) == 1);
        if abs(season_cyc_opt) < 1e-5
            % remove annual cycle


        %     OBP_ECCO(nan_OBP_ECCO_ind) = 0;
        %     G = [cos(((2*pi)/365.2425)*time_in_range_ECCO) sin(((2*pi)/365.2425)*time_in_range_ECCO) cos(((2*(2*pi))/365.2425)*time_in_range_ECCO) sin(((2*(2*pi))/365.2425)*time_in_range_ECCO) cos(((3*(2*pi))/365.2425)*time_in_range_ECCO) sin(((3*(2*pi))/365.2425)*time_in_range_ECCO) cos(((4*(2*pi))/365.2425)*time_in_range_ECCO) sin(((4*(2*pi))/365.2425)*time_in_range_ECCO)];
        %     coeffs = (((G')*G)^(-1))*((G')*(permute(reshape(OBP_ECCO,[(size(OBP_ECCO,1)*size(OBP_ECCO,2)) size(OBP_ECCO,3)]),[2 1])));
        %     OBP_ECCO = OBP_ECCO - reshape((G*coeffs)',size(OBP_ECCO));


            OBP_ECCO_ocean_mascons(nan_ECCO_ind) = 0;
            curr_ECCO_tseries(nan_curr_tseries_ind) = 0;

            G = [cos(((2*pi)/365.2425)*time_in_range_ECCO) sin(((2*pi)/365.2425)*time_in_range_ECCO) cos(((2*(2*pi))/365.2425)*time_in_range_ECCO) sin(((2*(2*pi))/365.2425)*time_in_range_ECCO) cos(((3*(2*pi))/365.2425)*time_in_range_ECCO) sin(((3*(2*pi))/365.2425)*time_in_range_ECCO) cos(((4*(2*pi))/365.2425)*time_in_range_ECCO) sin(((4*(2*pi))/365.2425)*time_in_range_ECCO)];
            coeffs = (((G')*G)^(-1))*((G')*(OBP_ECCO_ocean_mascons'));
            OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons - reshape((G*coeffs)',size(OBP_ECCO_ocean_mascons));
            G = [cos(((2*pi)/365.2425)*curr_tseries_time) sin(((2*pi)/365.2425)*curr_tseries_time) cos(((2*(2*pi))/365.2425)*curr_tseries_time) sin(((2*(2*pi))/365.2425)*curr_tseries_time) cos(((3*(2*pi))/365.2425)*curr_tseries_time) sin(((3*(2*pi))/365.2425)*curr_tseries_time) cos(((4*(2*pi))/365.2425)*curr_tseries_time) sin(((4*(2*pi))/365.2425)*curr_tseries_time)];
            coeffs = (((G')*G)^(-1))*((G')*curr_ECCO_tseries);
            curr_ECCO_tseries = curr_ECCO_tseries - (G*coeffs);
        end

        % curr_tseries(curr_nan_ind) = NaN;


        % OBP_ECCO(unique([nan_OBP_ECCO_ind; find(abs(OBP_ECCO_nan_mask_filtered) > 0.2)])) = NaN;

        OBP_ECCO_ocean_mascons(unique([nan_ECCO_ind; find(abs(ECCO_nan_mask_filtered) > 0.2)])) = NaN;
        curr_ECCO_tseries(unique([nan_curr_tseries_ind; find(abs(curr_tseries_nan_mask_filtered) > 0.2)])) = NaN;



        % compute covariances


        OBP_ECCO_ocean_mascons_array_1 = repmat(reshape(OBP_ECCO_ocean_mascons,[size(OBP_ECCO_ocean_mascons,1) 1 size(OBP_ECCO_ocean_mascons,2)]),[1 size(OBP_ECCO_ocean_mascons,1) 1]);
        OBP_ECCO_ocean_mascons_array_2 = repmat(reshape(OBP_ECCO_ocean_mascons,[1 size(OBP_ECCO_ocean_mascons,1) size(OBP_ECCO_ocean_mascons,2)]),[size(OBP_ECCO_ocean_mascons,1) 1 1]);
        lwe_thickness_ocean_mascons_array_1 = repmat(reshape(lwe_thickness_ocean_mascons,[size(lwe_thickness_ocean_mascons,1) 1 size(lwe_thickness_ocean_mascons,2)]),[1 size(lwe_thickness_ocean_mascons,1) 1]);
        lwe_thickness_ocean_mascons_array_2 = repmat(reshape(lwe_thickness_ocean_mascons,[1 size(lwe_thickness_ocean_mascons,1) size(lwe_thickness_ocean_mascons,2)]),[size(lwe_thickness_ocean_mascons,1) 1 1]);

        delta_lag = 365.24/12;
        lag_range_to_test = (365.24/12)*[0 1];

        [OBP_tseries_ECCO_corr_array,~,~,~,~,~,~] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_ECCO_tseries,[1 length(curr_ECCO_tseries)]),[size(OBP_ECCO_ocean_mascons,1) 1]),OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
        [OBP_ECCO_ocean_mascons_corr_array,~,~,~,~,~,lags_cov] = correlation_scalar_scalar_uncert_bounds(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
        [GRACE_ocean_mascons_corr_array,~,~,~,~,~,~] = correlation_scalar_scalar_uncert_bounds((9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_1,(9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95);
        [~,~,~,~,~,~,ECCO_std_dev_array_1,ECCO_std_dev_array_2,~] = regression_linear_scalar_scalar(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);
        [~,~,~,~,~,~,GRACE_std_dev_array_1,GRACE_std_dev_array_2,~] = regression_linear_scalar_scalar((9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_1,(9.81*1000)*0.01*lwe_thickness_ocean_mascons_array_2,3,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);
        [~,~,~,~,~,~,~,std_dev_curr_ECCO_tseries,~] = regression_linear_scalar_scalar(squeeze(OBP_ECCO_ocean_mascons_array_1(:,1,:)),repmat(reshape(curr_ECCO_tseries,[1 length(curr_ECCO_tseries)]),[size(OBP_ECCO_ocean_mascons_array_1,1) 1]),2,mean(diff(time_in_range_ECCO)),max([delta_lag (mean(diff(time_in_range_ECCO)) - 1e-2)]),lag_range_to_test,0.95,1);


        % use zero lag for covariances

        curr_lag = 0;
        OBP_tseries_ECCO_corr_array = OBP_tseries_ECCO_corr_array(:,abs(lags_cov - curr_lag) < 1e-5);
        OBP_ECCO_ocean_mascons_corr_array = OBP_ECCO_ocean_mascons_corr_array(:,:,abs(lags_cov - curr_lag) < 1e-5);
        GRACE_ocean_mascons_corr_array = GRACE_ocean_mascons_corr_array(:,:,abs(lags_cov - curr_lag) < 1e-5);
        ECCO_std_dev_array_1 = ECCO_std_dev_array_1(:,:,abs(lags_cov - curr_lag) < 1e-5);
        ECCO_std_dev_array_2 = ECCO_std_dev_array_2(:,:,abs(lags_cov - curr_lag) < 1e-5);
        GRACE_std_dev_array_1 = GRACE_std_dev_array_1(:,:,abs(lags_cov - curr_lag) < 1e-5);
        GRACE_std_dev_array_2 = GRACE_std_dev_array_2(:,:,abs(lags_cov - curr_lag) < 1e-5);
        std_dev_curr_ECCO_tseries = std_dev_curr_ECCO_tseries(:,abs(lags_cov - curr_lag) < 1e-5);
        
        % standard deviations of original and filtered time series
        ECCO_stddev_unfiltered(curr_pt) = OBP_pt_stddev(1);
        ECCO_stddev_filtered(curr_pt) = std_dev_curr_ECCO_tseries(1);


    %     % adjustments to cross-correlations
    %     
    %     mascon_depth_array = repmat(mascon_depth_avg,[1 length(mascon_depth_avg)]);
    %     mascon_depth_crossarray_avg = (0.5*mascon_depth_array) + (0.5*(mascon_depth_array'));
    %     mascon_depth_crossarray_diff = abs(mascon_depth_array - mascon_depth_array');
    %     mascon_tseries_depth_diff = abs(mascon_depth_avg - depth_pt_obs);
    %     
    %     mascon_row_id_array = repmat((1:1:length(mascon_depth_avg))',[1 length(mascon_depth_avg)]);
    %     mascon_col_id_array = repmat(1:1:length(mascon_depth_avg),[length(mascon_depth_avg) 1]);
    %     row_gt_col_ind = find(mascon_row_id_array - mascon_col_id_array > 1e-5);
    %     
    %     
    %     % average GRACE - ECCO correlation values in depth difference range
    %     depth_diff_min = 0;
    %     depth_diff_max = 1500;
    %     in_curr_range_ind = find((mascon_depth_crossarray_diff > depth_diff_min) & (mascon_depth_crossarray_diff < depth_diff_max));
    %     in_curr_range_unique_ind = intersect(in_curr_range_ind,row_gt_col_ind);
    %     GRACE_minus_ECCO_mean_in_range = mean(GRACE_ocean_mascons_corr_array(in_curr_range_unique_ind) - OBP_ECCO_ocean_mascons_corr_array(in_curr_range_unique_ind));
    %     curr_adjust = repmat(GRACE_minus_ECCO_mean_in_range,[length(in_curr_range_ind) 1]);
    %     hybrid_ocean_mascons_corr_array = OBP_ECCO_ocean_mascons_corr_array;  % use ECCO cross-correlation array as starting point
    %     % taper adjustment as correlations approach 1
    %     close_to_1_ind = find(OBP_ECCO_ocean_mascons_corr_array(in_curr_range_ind) - (1 - (2*curr_adjust)) > 0);
    %     curr_adjust(close_to_1_ind) = 0.5*(1 - OBP_ECCO_ocean_mascons_corr_array(in_curr_range_ind(close_to_1_ind)));
    %     
    %     hybrid_ocean_mascons_corr_array(in_curr_range_ind) = hybrid_ocean_mascons_corr_array(in_curr_range_ind) + (hybrid_factor*curr_adjust);
    %     
    %     in_curr_range_tseries_ind = find((mascon_tseries_depth_diff > depth_diff_min) & (mascon_tseries_depth_diff < depth_diff_max));
    %     hybrid_tseries_corr_array = OBP_tseries_ECCO_corr_array;
    %     curr_adjust = repmat(GRACE_minus_ECCO_mean_in_range,[length(in_curr_range_tseries_ind) 1]);
    %     % taper adjustment as correlations approach 1
    %     close_to_1_ind = find(OBP_tseries_ECCO_corr_array(in_curr_range_tseries_ind) - (1 - (2*curr_adjust)) > 0);
    %     curr_adjust(close_to_1_ind) = 0.5*(1 - OBP_tseries_ECCO_corr_array(in_curr_range_tseries_ind(close_to_1_ind)));
    %     
    %     hybrid_tseries_corr_array(in_curr_range_tseries_ind) = hybrid_tseries_corr_array(in_curr_range_tseries_ind) + (hybrid_factor*curr_adjust);
    %     
    %     
    %     % apply tapered adjustment to avoid discontinuity with previous adjustment
    %     depth_diff_min = 1500;
    %     adjust_at_min = GRACE_minus_ECCO_mean_in_range;
    %     depth_diff_max = 2000;
    %     adjust_at_max = 0;
    %     in_curr_range_ind = find((mascon_depth_crossarray_diff > depth_diff_min) & (mascon_depth_crossarray_diff < depth_diff_max));
    %     % in_curr_range_unique_ind = intersect(in_curr_range_ind,row_gt_col_ind);
    %     min_weight = (depth_diff_max - mascon_depth_crossarray_diff(in_curr_range_ind))/(depth_diff_max - depth_diff_min);
    %     max_weight = 1 - min_weight;
    %     curr_taper_adjust = (adjust_at_min*min_weight) + (adjust_at_max*max_weight);
    %     % taper adjustment as correlations approach 1
    %     close_to_1_ind = find(OBP_ECCO_ocean_mascons_corr_array(in_curr_range_ind) - (1 - (2*curr_taper_adjust)) > 0);
    %     curr_taper_adjust(close_to_1_ind) = 0.5*(1 - OBP_ECCO_ocean_mascons_corr_array(in_curr_range_ind(close_to_1_ind)));
    %     
    %     % hybrid_ocean_mascons_corr_array = OBP_ECCO_ocean_mascons_corr_array;  % use ECCO cross-correlation array as starting point
    %     hybrid_ocean_mascons_corr_array(in_curr_range_ind) = hybrid_ocean_mascons_corr_array(in_curr_range_ind) + (hybrid_factor*curr_taper_adjust);
    %     
    %     in_curr_range_tseries_ind = find((mascon_tseries_depth_diff > depth_diff_min) & (mascon_tseries_depth_diff < depth_diff_max));
    %     
    %     min_weight = (depth_diff_max - mascon_tseries_depth_diff(in_curr_range_tseries_ind))/(depth_diff_max - depth_diff_min);
    %     max_weight = 1 - min_weight;
    %     curr_taper_adjust = (adjust_at_min*min_weight) + (adjust_at_max*max_weight);
    %     % taper adjustment as correlations approach 1
    %     close_to_1_ind = find(OBP_tseries_ECCO_corr_array(in_curr_range_tseries_ind) - (1 - (2*curr_taper_adjust)) > 0);
    %     curr_taper_adjust(close_to_1_ind) = 0.5*(1 - OBP_tseries_ECCO_corr_array(in_curr_range_tseries_ind(close_to_1_ind)));
    %     
    %     % hybrid_tseries_corr_array = OBP_tseries_ECCO_corr_array;
    %     hybrid_tseries_corr_array(in_curr_range_tseries_ind) = hybrid_tseries_corr_array(in_curr_range_tseries_ind) + (hybrid_factor*curr_taper_adjust);
    %     
    %     
    %     % only include mascons whose correlation with obs. point is above a threshold
    %     [sorted_hybrid_tseries_corr,~] = sort(hybrid_tseries_corr_array,'descend');
    %     sorted_hybrid_tseries_corr = sorted_hybrid_tseries_corr(isnan(sorted_hybrid_tseries_corr) == 0);
    %     corr_threshold_n_mascons = sorted_hybrid_tseries_corr(n_mascons_max) - 1e-10;
    %     high_corr_ind = find(hybrid_tseries_corr_array >= max([corr_threshold_n_mascons min_corr_to_include]));
    %     disp(['corr_threshold_n_mascons = ',num2str(corr_threshold_n_mascons)])
    %     if length(high_corr_ind) < 3
    % %         keyboard
    %         continue
    %     end
    %     
    % %     ocean_mascons_cov = hybrid_ocean_mascons_corr_array(high_corr_ind,high_corr_ind).*GRACE_std_dev_array_1(high_corr_ind,high_corr_ind).*GRACE_std_dev_array_2(high_corr_ind,high_corr_ind);
    %     ocean_mascons_cov = hybrid_ocean_mascons_corr_array(high_corr_ind,high_corr_ind).*(((1 - hybrid_factor_stddev)*ECCO_std_dev_array_1(high_corr_ind,high_corr_ind)) + (hybrid_factor_stddev*GRACE_std_dev_array_1(high_corr_ind,high_corr_ind))).*(((1 - hybrid_factor_stddev)*ECCO_std_dev_array_2(high_corr_ind,high_corr_ind)) + (hybrid_factor_stddev*GRACE_std_dev_array_2(high_corr_ind,high_corr_ind)));
    %     % ocean_mascons_cov_diag = diag(ocean_mascons_cov,0);
    %     % ocean_mascons_cov = ocean_mascons_cov + (noise_param*diag(ocean_mascons_cov_diag,0));
    %     
    % %     ocean_mascons_tseries_cov = hybrid_tseries_corr_array(high_corr_ind).*GRACE_std_dev_array_1(high_corr_ind,1).*std_dev_curr_ECCO_tseries(high_corr_ind);
    %     ocean_mascons_tseries_cov = hybrid_tseries_corr_array(high_corr_ind).*(((1 - hybrid_factor_stddev)*ECCO_std_dev_array_1(high_corr_ind,1)) + (hybrid_factor_stddev*GRACE_std_dev_array_1(high_corr_ind,1))).*std_dev_curr_ECCO_tseries(high_corr_ind);


        % only include mascons whose correlation with obs. point is above a threshold
        [sorted_tseries_corr,~] = sort(OBP_tseries_ECCO_corr_array,'descend');
        sorted_tseries_corr = sorted_tseries_corr(isnan(sorted_tseries_corr) == 0);
        corr_threshold_n_mascons = sorted_tseries_corr(n_mascons_max) - 1e-10;
        high_corr_ind = find(OBP_tseries_ECCO_corr_array >= max([corr_threshold_n_mascons min_corr_to_include]));
        disp(['corr_threshold_n_mascons = ',num2str(corr_threshold_n_mascons)])
        if length(high_corr_ind) < 3
    %         keyboard
            continue
        end
        
        
        if abs(downscale_adj_opt - 1) < 1e-5
            
%             % adjust mascon-point correlations based on local bathymetry slope, orientation and proximity
% 
%             dist_max_adjust = 2;     % distance of maximum adjustment, in degrees latitude
%             dist_radius_adjust = 2;    % radius of adjustment region from max adjust point, in degrees latitude
%             
%             % apply adjustment only if slope is steep enough, and depth at pt_obs is greater than the average of the surrounding region
%             adjust_corr_vec = zeros([length(high_corr_ind) 1]);
%             if ((range_depth_planar_fit > 1000) && (depth_pt_obs > mean_depth_within_radius))
%                 lon_center_pos_adjust = lon_pt_obs + (cos(angle_steepest_descent)*((1/cosd(lat_pt_obs))*dist_max_adjust));
%                 lat_center_pos_adjust = lat_pt_obs + (sin(angle_steepest_descent)*dist_max_adjust);
%                 mascon_dist_from_pos_center = abs(((cosd(lat_center_pos_adjust))*(mod(mascon_lon_center(high_corr_ind) - lon_center_pos_adjust + 180,360) - 180)) + (1i*(mascon_lat_center(high_corr_ind) - lat_center_pos_adjust)));
%                 in_range_pos_adjust_ind = find(mascon_dist_from_pos_center < dist_radius_adjust);
%                 adjust_corr_vec(in_range_pos_adjust_ind) = adjust_corr_vec(in_range_pos_adjust_ind) + (adjust_corr_max*(1 - ((mascon_dist_from_pos_center(in_range_pos_adjust_ind)/dist_radius_adjust).^2)));
% 
%                 lon_center_neg_adjust = lon_pt_obs - (cos(angle_steepest_descent)*((1/cosd(lat_pt_obs))*dist_max_adjust));
%                 lat_center_neg_adjust = lat_pt_obs - (sin(angle_steepest_descent)*dist_max_adjust);
%                 mascon_dist_from_neg_center = abs(((cosd(lat_center_neg_adjust))*(mod(mascon_lon_center(high_corr_ind) - lon_center_neg_adjust + 180,360) - 180)) + (1i*(mascon_lat_center(high_corr_ind) - lat_center_neg_adjust)));
%                 in_range_neg_adjust_ind = find(mascon_dist_from_neg_center < dist_radius_adjust);
%                 adjust_corr_vec(in_range_neg_adjust_ind) = adjust_corr_vec(in_range_neg_adjust_ind) - (adjust_corr_max*(1 - ((mascon_dist_from_neg_center(in_range_neg_adjust_ind)/dist_radius_adjust).^2)));
%                 
%             end
            

            
            % adjust mascon-point correlations based on mascon depth relative to obs. point depth

            depth_radius_adjust = 2500;     % radius of depth to apply (positive) adjustment
            
%             adjust_corr_array = zeros([length(high_corr_ind) length(high_corr_ind)]);
%             mascon_depth_array = repmat(mascon_depth_avg(high_corr_ind),[1 length(high_corr_ind)]);
%             mascon_depth_crossarray_diff = abs(mascon_depth_array - mascon_depth_array');
%             in_range_pos_adjust_ind = find(mascon_depth_crossarray_diff < depth_radius_adjust);
%             adjust_corr_array(in_range_pos_adjust_ind) = adjust_corr_array(in_range_pos_adjust_ind) + (adjust_corr_max*(1 - (mascon_depth_crossarray_diff(in_range_pos_adjust_ind)/depth_radius_adjust)));
%             OBP_ECCO_ocean_mascons_highcorr_array = OBP_ECCO_ocean_mascons_corr_array(high_corr_ind,high_corr_ind) + adjust_corr_array;
%             OBP_ECCO_ocean_mascons_highcorr_array(OBP_ECCO_ocean_mascons_highcorr_array > 1) = 1;
%             OBP_ECCO_ocean_mascons_highcorr_array(OBP_ECCO_ocean_mascons_highcorr_array < -1) = -1;
%             
%             ocean_mascons_cov = OBP_ECCO_ocean_mascons_highcorr_array.*(ECCO_std_dev_array_1(high_corr_ind,high_corr_ind)).*(ECCO_std_dev_array_2(high_corr_ind,high_corr_ind));
            
            adjust_corr_vec = zeros([length(high_corr_ind) 1]);
            mascon_tseries_depth_diff = abs(mascon_depth_avg(high_corr_ind) - depth_pt_obs);
            in_range_pos_adjust_ind = find(mascon_tseries_depth_diff < depth_radius_adjust);
            adjust_corr_vec(in_range_pos_adjust_ind) = adjust_corr_vec(in_range_pos_adjust_ind) + (adjust_corr_max*(1 - (mascon_tseries_depth_diff(in_range_pos_adjust_ind)/depth_radius_adjust)));
            OBP_tseries_ECCO_highcorr_array = OBP_tseries_ECCO_corr_array(high_corr_ind) + adjust_corr_vec;
            OBP_tseries_ECCO_highcorr_array(OBP_tseries_ECCO_highcorr_array > 1) = 1;
            OBP_tseries_ECCO_highcorr_array(OBP_tseries_ECCO_highcorr_array < -1) = -1;
            
            ocean_mascons_tseries_cov = OBP_tseries_ECCO_highcorr_array.*(ECCO_std_dev_array_1(high_corr_ind,1)).*std_dev_curr_ECCO_tseries(high_corr_ind);
            
                        
%             % adjust mascon-point correlations based on mascon depth relative to obs. point depth
% 
%             depth_radius_adjust = 1000;     % radius of depth equivalent (with f = value at obs. lat.) to apply (positive) adjustment
%             
%             adjust_corr_vec = zeros([length(high_corr_ind) 1]);
%             f_pt_obs = 2*(2*pi/86164)*(sind(lat_pt_obs));
%             f_mascon = 2*(2*pi/86164)*(sind(mascon_lat_center(high_corr_ind)));
%             f_over_H_mascon_tseries_diff = (f_mascon./mascon_depth_avg(high_corr_ind)) - (f_pt_obs/depth_pt_obs);
%             in_range_pos_adjust_deeper_ind = find((f_over_H_mascon_tseries_diff > (f_pt_obs/(depth_pt_obs + depth_radius_adjust)) - (f_pt_obs/depth_pt_obs)) & (f_over_H_mascon_tseries_diff < 0));
%             adjust_corr_vec(in_range_pos_adjust_deeper_ind) = adjust_corr_vec(in_range_pos_adjust_deeper_ind) + (adjust_corr_max*(1 - (f_over_H_mascon_tseries_diff(in_range_pos_adjust_deeper_ind)/((f_pt_obs/(depth_pt_obs + depth_radius_adjust)) - (f_pt_obs/depth_pt_obs)))));
%             in_range_pos_adjust_shallower_ind = find((f_over_H_mascon_tseries_diff < (f_pt_obs/(depth_pt_obs - depth_radius_adjust)) - (f_pt_obs/depth_pt_obs)) & (f_over_H_mascon_tseries_diff >= 0));
%             adjust_corr_vec(in_range_pos_adjust_shallower_ind) = adjust_corr_vec(in_range_pos_adjust_shallower_ind) + (adjust_corr_max*(1 - (f_over_H_mascon_tseries_diff(in_range_pos_adjust_shallower_ind)/((f_pt_obs/(depth_pt_obs - depth_radius_adjust)) - (f_pt_obs/depth_pt_obs)))));
%             
%             ocean_mascons_tseries_cov = (OBP_tseries_ECCO_corr_array(high_corr_ind) + adjust_corr_vec).*(ECCO_std_dev_array_1(high_corr_ind,1)).*std_dev_curr_ECCO_tseries(high_corr_ind);
            
        else
%             ocean_mascons_cov = (OBP_ECCO_ocean_mascons_corr_array(high_corr_ind,high_corr_ind)).*(ECCO_std_dev_array_1(high_corr_ind,high_corr_ind)).*(ECCO_std_dev_array_2(high_corr_ind,high_corr_ind));
            ocean_mascons_tseries_cov = (OBP_tseries_ECCO_corr_array(high_corr_ind)).*(ECCO_std_dev_array_1(high_corr_ind,1)).*std_dev_curr_ECCO_tseries(high_corr_ind);
        end

        ocean_mascons_cov = (OBP_ECCO_ocean_mascons_corr_array(high_corr_ind,high_corr_ind)).*(ECCO_std_dev_array_1(high_corr_ind,high_corr_ind)).*(ECCO_std_dev_array_2(high_corr_ind,high_corr_ind));
        
        gain_vec = (ocean_mascons_cov^(-1))*ocean_mascons_tseries_cov;
        obs_tseries_reconstr = ((9.81*1000)*0.01*(lwe_thickness_ocean_mascons(high_corr_ind,:)'))*gain_vec;

    else
        colloc_mascon_ind = find((mod(lon_pt_obs - mascon_lon_bounds(:,1) + 180,360) - 180 >= 0) & (mod(lon_pt_obs - mascon_lon_bounds(:,2) + 180,360) - 180 < 0) & (lat_pt_obs - mascon_lat_bounds(:,1) >= 0) & (lat_pt_obs - mascon_lat_bounds(:,2) < 0));
        obs_tseries_reconstr = (9.81*1000)*0.01*(lwe_thickness_ocean_mascons(colloc_mascon_ind,:)');
    end
    
    
    % compute correlation of original and reconstructed OBP
    
    obs_tseries_tinterp(abs(obs_tseries_tinterp) < 1e-15) = NaN;
    obs_tseries_reconstr(abs(obs_tseries_reconstr) < 1e-15) = NaN;
    
    
    delta_lag = 365.24/12;
    lag_range_to_test = (365.24/12)*[0 12];
    conf_level = 0.95;
    
    [corr_obs_tseries_reconstr,~,~,~,corr_lowmag_bound,corr_highmag_bound,lags] = correlation_scalar_scalar_uncert_bounds(obs_tseries_tinterp,obs_tseries_reconstr,1,mode(diff(time)),delta_lag,lag_range_to_test,conf_level);
    zero_lag_ind = find(abs(lags) < 1e-10);
    corr_obs_tseries_reconstr_zerolag = corr_obs_tseries_reconstr(zero_lag_ind);
    corr_lowmag_bound_zerolag = corr_lowmag_bound(zero_lag_ind);
    corr_highmag_bound_zerolag = corr_highmag_bound(zero_lag_ind);
    
    
    obs_tseries_reconstr_corrcoeff(curr_pt) = corr_obs_tseries_reconstr_zerolag;
    obs_tseries_reconstr_lowmag_bound(curr_pt) = corr_lowmag_bound_zerolag;
    obs_tseries_reconstr_highmag_bound(curr_pt) = corr_highmag_bound_zerolag;
    
    disp(['Computed correlation for site number ',num2str(curr_pt),' (out of ',num2str(length(filenames_cellarray)),')'])
    
end


% if ((hybrid_factor > 1e-5) && (hybrid_factor < 1 - 1e-5))
%     hybrid_str = ['p',num2str(round(hybrid_factor/0.1))];
% else
%     hybrid_str = num2str(hybrid_factor);
% end
% if ((hybrid_factor_stddev > 1e-5) && (hybrid_factor_stddev < 1 - 1e-5))
%     hybrid_std_str = ['p',num2str(round(hybrid_factor_stddev/0.1))];
% else
%     hybrid_std_str = num2str(hybrid_factor_stddev);
% end
if abs(downscale_opt) < 1e-5
    corrcoeff_coloc_mascons = obs_tseries_reconstr_corrcoeff;
    lowmag_bound_coloc_mascons = obs_tseries_reconstr_lowmag_bound;
    highmag_bound_coloc_mascons = obs_tseries_reconstr_highmag_bound;
elseif abs(downscale_opt - 1) < 1e-5
    if abs(downscale_adj_opt) < 1e-5
        eval(['corrcoeff_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_corrcoeff;'])
        eval(['lowmag_bound_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_lowmag_bound;'])
        eval(['highmag_bound_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_highmag_bound;'])
    elseif abs(downscale_adj_opt - 1) < 1e-5
        if mod(adjust_corr_max,0.1) == 0
            adjust_corr_str = num2str(round(adjust_corr_max/0.1));
        else
            if adjust_corr_max < 0.1
                adjust_corr_str = ['0',num2str(round(adjust_corr_max/0.01))];
            else
                adjust_corr_str = num2str(round(adjust_corr_max/0.01));
            end
        end
%         eval(['corrcoeff_adjcorrp',adjust_corr_str,'_distadj',num2str(dist_max_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_corrcoeff;'])
%         eval(['lowmag_bound_adjcorrp',adjust_corr_str,'_distadj',num2str(dist_max_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_lowmag_bound;'])
%         eval(['highmag_bound_adjcorrp',adjust_corr_str,'_distadj',num2str(dist_max_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_highmag_bound;'])

        eval(['corrcoeff_adjcorrp',adjust_corr_str,'_depthadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_corrcoeff;'])
        eval(['lowmag_bound_adjcorrp',adjust_corr_str,'_depthadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_lowmag_bound;'])
        eval(['highmag_bound_adjcorrp',adjust_corr_str,'_depthadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_highmag_bound;'])
%         eval(['corrcoeff_adjcorrp',adjust_corr_str,'_foverHadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_corrcoeff;'])
%         eval(['lowmag_bound_adjcorrp',adjust_corr_str,'_foverHadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_lowmag_bound;'])
%         eval(['highmag_bound_adjcorrp',adjust_corr_str,'_foverHadj',num2str(depth_radius_adjust),'_maxmascons',num2str(n_mascons_max),'_mincorrp',num2str(10*min_corr_to_include),' = obs_tseries_reconstr_highmag_bound;'])
    end
end


save(['DARTplus_downscale_corrcoeff_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods.mat'],'corrcoeff_*','lowmag_bound_*','highmag_bound_*','-append')




fig1000 = figure(1000);
fig_paper_pos = get(fig1000,'PaperPosition');
fig_paper_pos(4) = ((latnorth - latsouth)/((longeast - longwest)*(cosd(mean([latsouth latnorth])))))*fig_paper_pos(3);
fig_paper_pos(3:4) = 4*fig_paper_pos(3:4);
fig_pos = get(fig1000,'Position');
fig_pos(4) = ((latnorth - latsouth)/((longeast - longwest)*(cosd(mean([latsouth latnorth])))))*fig_pos(3);
fig_pos(3:4) = 4*fig_pos(3:4);
set(fig1000,'PaperPosition',fig_paper_pos,'Position',fig_pos,'PaperSize',[44 34])
h = plot([-1000 -900],[0 0]);
set(h,'Color',[1 1 1],'LineWidth',0.5)
set(gca,'xlim',[longwest longeast],'ylim',[latsouth latnorth])
daspect([1 cosd(mean([latsouth latnorth])) 1])
set(gca,'position',[0 0 1 1],'units','normalized','Visible','off')
pause(20)
fig_pos = get(fig1000,'Position');
fig_pos(4) = ((latnorth - latsouth)/((longeast - longwest)*(cosd(mean([latsouth latnorth])))))*fig_pos(3);
set(fig1000,'Position',fig_pos)
print(fig1000,'contour_plot_temporary.png','-dpng','-r300')

close(fig1000)


% get contour plot and overlay land mask

rgb_array = imread('contour_plot_temporary.png');
rgb_array = flip(permute(rgb_array,[2 1 3]),2);

delete('contour_plot_temporary.png')

size_rgb_array = size(rgb_array);
pixel_dim = size_rgb_array(1:2);


% load bathymetry data

%begin bathymetry calculations

ncid_topo = netcdf.open('~/plotting_scripts/topo2.grd','NC_NOWRITE');
z = netcdf.getVar(ncid_topo,5);

spacing = 1/30;     %spacing of bathymetry data, in degrees

%masking bathymetry threshold (in meters elevation)

mask_threshold = 0;


lat_bathy = 90 - ((0.5*spacing):spacing:180)';

long_bathy = ((0.5*spacing):spacing:360)';


z_reshaped = reshape(z,360/spacing,180/spacing);

clear z


z_reshaped = flip(z_reshaped,2);

lat_bathy = flip(lat_bathy,1);
% resolutions of imported image

x_res_image = (mod(longeast - longwest - 1e-10,360) + 1e-10)/(pixel_dim(1));
y_res_image = (latnorth - latsouth)/(pixel_dim(2));

x_res_ratio = x_res_image/spacing;
y_res_ratio = y_res_image/spacing;


% % apply smoothing filter
% x_smooth_filter_vec = zeros(size(z_land_mask,1),1);
% y_smooth_filter_vec = zeros(size(z_land_mask,2),1);
x_smooth_filter_vec = zeros(size(z_reshaped,1),1);
y_smooth_filter_vec = zeros(size(z_reshaped,2),1);

% edge_full_ind = floor((x_res_ratio + 1)/2);
% x_smooth_filter_vec(1:edge_full_ind) = 1;
% x_smooth_filter_vec(edge_full_ind + 1) = ((x_res_ratio + 1)/2) - edge_full_ind;
% x_smooth_filter_vec((length(x_smooth_filter_vec) - edge_full_ind + 1):length(x_smooth_filter_vec)) = flipdim(x_smooth_filter_vec(2:(edge_full_ind + 1)),1);
% 
% y_smooth_filter_vec(1:edge_full_ind) = 1;
% y_smooth_filter_vec(edge_full_ind + 1) = ((y_res_ratio + 1)/2) - edge_full_ind;
% y_smooth_filter_vec((length(y_smooth_filter_vec) - edge_full_ind + 1):length(y_smooth_filter_vec)) = flipdim(y_smooth_filter_vec(2:(edge_full_ind + 1)),1);

% lower-resolution smoothing (for quicker contouring)
edge_full_ind = floor(x_res_ratio + 1);
x_smooth_filter_vec(1:edge_full_ind) = 1;
x_smooth_filter_vec(edge_full_ind + 1) = (x_res_ratio + 1) - edge_full_ind;
x_smooth_filter_vec((length(x_smooth_filter_vec) - edge_full_ind + 1):length(x_smooth_filter_vec)) = flip(x_smooth_filter_vec(2:(edge_full_ind + 1)),1);
y_smooth_filter_vec(1:edge_full_ind) = 1;
y_smooth_filter_vec(edge_full_ind + 1) = (y_res_ratio + 1) - edge_full_ind;
y_smooth_filter_vec((length(y_smooth_filter_vec) - edge_full_ind + 1):length(y_smooth_filter_vec)) = flip(y_smooth_filter_vec(2:(edge_full_ind + 1)),1);

% x_smooth_filter_array = repmat(x_smooth_filter_vec,[1 size(z_land_mask,2)]);
% y_smooth_filter_array = repmat(y_smooth_filter_vec',[size(z_land_mask,1) 1]);
x_smooth_filter_array = repmat(x_smooth_filter_vec,[1 size(z_reshaped,2)]);
y_smooth_filter_array = repmat(y_smooth_filter_vec',[size(z_reshaped,1) 1]);

smooth_filter_array = x_smooth_filter_array.*y_smooth_filter_array;
smooth_filter_array = smooth_filter_array/(sum(sum(smooth_filter_array)));

fft_smooth_filter_array = fft(fft(smooth_filter_array,[],1),[],2);

% z_land_mask_smoothed = ifft(ifft(fft_smooth_filter_array.*fft(fft(z_land_mask,[],1),[],2),[],2),[],1);
z_reshaped_smoothed = ifft(ifft(fft_smooth_filter_array.*fft(fft(z_reshaped,[],1),[],2),[],2),[],1);


% "interpolate" to the same grid as the image

if mod(longwest,360) >= mod(longeast,360)
    in_bathy_lon_range_ind_1 = find(long_bathy >= mod(longwest,360));
    in_bathy_lon_range_ind_2 = find(long_bathy <= mod(longeast,360));
    in_bathy_lon_range_ind = [in_bathy_lon_range_ind_1; in_bathy_lon_range_ind_2];
else
    in_bathy_lon_range_ind = find((long_bathy >= mod(longwest,360)) & (long_bathy <= mod(longeast,360)));
end
in_bathy_lat_range_ind = find((lat_bathy >= latsouth) & (lat_bathy <= latnorth));
SW_corner_x_weight_inner = ((mod(longwest - long_bathy(in_bathy_lon_range_ind(1)) + 180,360) - 180)/spacing) + 1;
SW_corner_y_weight_inner = ((latsouth - lat_bathy(min(in_bathy_lat_range_ind)))/spacing) + 1;

x_rel_placing_vec = SW_corner_x_weight_inner + (x_res_ratio*((0:1:(pixel_dim(1) - 1))'));

% land_mask_x_interp = NaN(pixel_dim(1),size(z_land_mask_smoothed,2));
z_reshaped_x_interp = NaN(pixel_dim(1),size(z_reshaped_smoothed,2));
for x_ind = 1:length(x_rel_placing_vec)
    curr_x_rel_placing = x_rel_placing_vec(x_ind);
    
    x_bathy_lower_ind = mod(in_bathy_lon_range_ind(1) - 1 + floor(curr_x_rel_placing) - 1,length(long_bathy)) + 1;
    x_bathy_upper_ind = mod(x_bathy_lower_ind + 1 - 1,length(long_bathy)) + 1;
    
    x_upper_weight = curr_x_rel_placing - floor(curr_x_rel_placing);
    x_lower_weight = 1 - x_upper_weight;
    
%     land_mask_x_interp(x_ind,:) = (x_lower_weight*z_land_mask_smoothed(x_bathy_lower_ind,:)) + (x_upper_weight*z_land_mask_smoothed(x_bathy_upper_ind,:)); 
    z_reshaped_x_interp(x_ind,:) = (x_lower_weight*z_reshaped_smoothed(x_bathy_lower_ind,:)) + (x_upper_weight*z_reshaped_smoothed(x_bathy_upper_ind,:)); 
    
end

y_rel_placing_vec = SW_corner_y_weight_inner + (y_res_ratio*((0:1:(pixel_dim(2) - 1))'));

% land_mask = NaN(pixel_dim);
z_reshaped_interp = NaN(pixel_dim);
for y_ind = 1:length(y_rel_placing_vec)
    curr_y_rel_placing = y_rel_placing_vec(y_ind);
    
    y_bathy_lower_ind = max([1 (min(in_bathy_lat_range_ind) - 1 + floor(curr_y_rel_placing))]);
    y_bathy_upper_ind = min([length(lat_bathy) (y_bathy_lower_ind + 1)]);
    
    y_upper_weight = curr_y_rel_placing - floor(curr_y_rel_placing);
    y_lower_weight = 1 - y_upper_weight;
    
%     land_mask(:,y_ind) = (y_lower_weight*land_mask_x_interp(:,y_bathy_lower_ind)) + (y_upper_weight*land_mask_x_interp(:,y_bathy_upper_ind)); 
    z_reshaped_interp(:,y_ind) = (y_lower_weight*z_reshaped_x_interp(:,y_bathy_lower_ind)) + (y_upper_weight*z_reshaped_x_interp(:,y_bathy_upper_ind)); 
    
end


% find indices of land points
% landmask_ind = find(land_mask > 0.5);
landmask_ind = find(z_reshaped_interp >= mask_threshold);


% frac_first_part = (min([longeast_offset 0]) - max([(-360) longwest_offset]))/(min([longeast_offset 360]) - max([(-360) longwest_offset]));
% if frac_first_part <= 0
%     landmask_ind = landfind_indices(longwest_offset,longeast_offset,latsouth,latnorth,size_rgb_array(1:2));
% elseif frac_first_part >= 1
%     landmask_ind = landfind_indices(longwest_offset + 360,longeast_offset + 360,latsouth,latnorth,size_rgb_array(1:2));
% else
%     landmask_i_ind = [(mod(landfind_indices(longwest_offset + 360,(longwest_offset*(1 - (floor(frac_first_part*(size_rgb_array(1)))/(frac_first_part*(size_rgb_array(1)))))) + 360,latsouth,latnorth,floor([frac_first_part 1].*(size_rgb_array(1:2)))) - 1,floor(frac_first_part*(size_rgb_array(1)))) + 1); (mod(landfind_indices(longwest_offset*(1 - (floor(frac_first_part*(size_rgb_array(1)))/(frac_first_part*(size_rgb_array(1))))),longeast_offset,latsouth,latnorth,size_rgb_array(1:2) - floor([frac_first_part 0].*(size_rgb_array(1:2)))) - 1,size_rgb_array(1) - floor(frac_first_part*(size_rgb_array(1)))) + 1 + (floor(frac_first_part*(size_rgb_array(1)))))];
%     landmask_j_ind = [(ceil(landfind_indices(longwest_offset + 360,(longwest_offset*(1 - (floor(frac_first_part*(size_rgb_array(1)))/(frac_first_part*(size_rgb_array(1)))))) + 360,latsouth,latnorth,floor([frac_first_part 1].*(size_rgb_array(1:2))))/(floor(frac_first_part*(size_rgb_array(1)))))); (ceil(landfind_indices(longwest_offset*(1 - (floor(frac_first_part*(size_rgb_array(1)))/(frac_first_part*(size_rgb_array(1))))),longeast_offset,latsouth,latnorth,size_rgb_array(1:2) - floor([frac_first_part 0].*(size_rgb_array(1:2))))/(size_rgb_array(1) - floor(frac_first_part*(size_rgb_array(1))))))];
%     landmask_ind = double(landmask_i_ind) + ((size_rgb_array(1))*(double(landmask_j_ind) - 1));
% end
rgb_array_reshaped = reshape(rgb_array,[prod(size_rgb_array(1:2)) 3]);

% convert black undefined areas in SSH field to white shading
black_shaded_ind = find((rgb_array_reshaped(:,1) == 0) & (rgb_array_reshaped(:,2) == 0) & (rgb_array_reshaped(:,3) == 0));
rgb_array_reshaped(black_shaded_ind,:) = 255*ones(length(black_shaded_ind),3);

% put black mask on land areas
rgb_array_reshaped(landmask_ind,:) = zeros(length(landmask_ind),3);

rgb_array_masked = reshape(rgb_array_reshaped,size_rgb_array);


% plot contour map with land mask

% bathy_CI = 2000;
% bathy_contours = ((-bathy_CI):(-bathy_CI):(-6000))';
% bathy_contours = [-2000; -2000];

% color levels for correlation between obs. time series and downscale estimate
c_levels = ((-1):0.05:1)';
cmap = colormap(0.85*bcyr(length(c_levels) - 1));

[irregular_clevels_vec,~] = irregular_clevels_plot(obs_tseries_reconstr_corrcoeff,c_levels);

fig3 = figure(3);
fig_pos = get(fig3,'Position');
fig_pos(3:4) = 4*fig_pos(3:4);
set(fig3,'Position',fig_pos,'PaperSize',[11 8.5])
close(figure(1))
axm = axesm(projection,'origin',[0 mean([longwest longeast]) 0]);
h = geoshow(repmat((latsouth:((latnorth - latsouth)/(size(rgb_array,2) - 1)):latnorth)',[1 size(rgb_array,1)]),repmat(longwest:((longeast - longwest)/(size(rgb_array,1) - 1)):longeast,[size(rgb_array,2) 1]),permute(rgb_array_masked,[2 1 3]));
setm(axm,'frame','on','grid','on','flonlimit',[longwest longeast] - mean([longwest longeast]),'flatlimit',[latsouth latnorth],'mlinelocation',30*round(((longeast - longwest)/6)/30),'plinelocation',30,'glinestyle','-','meridianlabel','on','parallellabel','on','mlabellocation',60,'mlabelparallel',latnorth - 0.01,'plabellocation',30,'plabelmeridian',longwest - 10,'FontSize',12)
set(axm,'XColor','none','YColor','none')
children_handles = get(axm,'Children');
set(children_handles(1),'String','90^{\circ} N        ')
set(children_handles(2),'String','60^{\circ} N  ')
set(children_handles(6),'String','60^{\circ} S  ')
set(children_handles(7),'String','90^{\circ} S        ')
set(children_handles(8),'String','  0^{\circ} ')
set(children_handles(9),'String','  60^{\circ} W')
set(children_handles(10),'String','  120^{\circ} W')
set(children_handles(11),'String','  180^{\circ}')
set(children_handles(12),'String','  120^{\circ} E')
set(children_handles(13),'String','  60^{\circ} E')
for curr_handle = 8:13
    try
        curr_pos = get(children_handles(curr_handle),'Position');
        set(children_handles(curr_handle),'Position',curr_pos + [0 0.07 0])
    catch
    end
end
hold on
% contourm_overlay(longwest:((longeast - longwest)/(size(rgb_array,1) - 1)):longeast,(latsouth:((latnorth - latsouth)/(size(rgb_array,2) - 1)):latnorth)',real(z_reshaped_interp)',bathy_contours,[0.5 0.5 0.5],'-',0.5)
for pt_ind = 1:length(obs_lon)
    if isnan(irregular_clevels_vec(pt_ind)) == 1
        continue
    end
    pt_radius = 0.03;
    center_lon_lat = [((mod(obs_lon(pt_ind) - longwest,360) + longwest)) obs_lat(pt_ind)];
    [center_x,center_y] = projfwd(getm(axm),center_lon_lat(2),center_lon_lat(1));
%     h = rectangle('Position',[((mod(obs_lon(pt_ind) - longwest,360) + longwest) - (pt_radius/cosd(obs_lat(pt_ind)))) (obs_lat(pt_ind) - pt_radius) (2*pt_radius/cosd(obs_lat(pt_ind))) (2*pt_radius)],'Curvature',[1 1]);
    h = rectangle('Position',[(center_x - pt_radius) (center_y - pt_radius) (2*pt_radius) (2*pt_radius)],'Curvature',[1 1]);
    set(h,'facecolor',cmap(ceil(irregular_clevels_vec(pt_ind)),:),'edgecolor','none')
end
hold off
title({'Zero-lag correlations of in-situ OBP with co-located downscaling estimate,'; [num2str(round(1/high_freq_bound)),'-',num2str(round(1/low_freq_bound)),' day periods, with hybrid factor = ',num2str(hybrid_factor),', max mascons = ',num2str(n_mascons_max),', min corr. = ',num2str(min_corr_to_include)]; ' '},'FontSize',10)
cbar = colorbar('location','southoutside');
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'ticks',0:10:40,'xticklabel',{'-1' '-0.5' '0' '0.5' '1'},'FontSize',12)
set(get(cbar,'label'),'String','Correlation coefficient','FontSize',14)

print(fig3,['Insitu_OBP_corrcoeff_downscale_',downscale_id,'_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(datenum(time_range_start),'yyyymm'),'_',datestr(datenum(time_range_end) - 1,'yyyymm'),'_time_bounds.pdf'],'-dpdf','-r300')
close(fig3)