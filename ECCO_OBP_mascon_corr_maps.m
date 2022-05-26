% correlate between OBP at a point (in ECCO and/or obs.) and ECCO OBP in surrounding region

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
path(path,'~/plotting_scripts/gsw_matlab/')
cd('/indopac/adelman/GRACE/')


% OBP observation points from SAMOC PIES
% site A: -51.5, -34.5
% site B: -49.5, -34.5
% site C: -47.5, -34.5
% site D: -44.5, -34.5


ECCO_id = 'cs510_cylindweight1.5';   % ECCO simulation configuration used

lon_pt = -70.59;    % longitude of point time series
lat_pt = 39.49;     % latitude of point time series
lon_bounds = [-90 -30];    % longitude bounds for map plots
lat_bounds = [20 60];    % latitude bounds for map plots

pt_ECCO_obs_opt = 0;     % 0 = use point from ECCO, 1 = use point from obs.
pt_text_id = 'ECCO cs510 -70.59 lon, 39.49 lat';
pt_id = 'ECCO_cs510_-70.59_lon_39.49_lat';

time_range_start = [1992 1 1];
time_range_end = [2019 1 1];

low_freq_bound = 1/365.24;
high_freq_bound = 1/((3/12)*365.24);
season_cyc_opt = 0;    % 0 = remove seasonal/annual cycle, 1 = retain seasonal/annual cycle
monthavg_opt = 1;     % 0 = no monthly averaging, 1 = monthly averaging
reconstr_opt = 0;     % 0 = plot time series-mascon correlations; 1 = plot time series-reconstruction correlations

% lag parameters for correlations
% delta_lag = (1/72)*365.24;
% delta_lag = 2;
delta_lag = 365.24/12;
lag_range_to_test = (365.24/12)*[-1 1];


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

interp_spacing = 0.1;
lon_z_interp = ((lon_bounds(1) - (0.5*interp_spacing)):interp_spacing:(lon_bounds(2) + (0.5*interp_spacing)))';
lat_z_interp = ((lat_bounds(1) - (0.5*interp_spacing)):interp_spacing:(lat_bounds(2) + (0.5*interp_spacing)))';

warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')

in_local_bathy_lat_range_ind = find((lat_bathy - lat_bounds(1) >= (-10*spacing)) & (lat_bathy - lat_bounds(2) <= (10*spacing)));
local_lat_bathy = lat_bathy(in_local_bathy_lat_range_ind);
if ((mod(lon_bounds(2) + (10*spacing) - lon_bounds(1) - 1e-5,360) + 1e-5 > mod(max(lon_bathy) - lon_bounds(1) - 1e-5,360) + 1e-5) || (mod(lon_bounds(2) - (lon_bounds(1) - (10*spacing)) - 1e-5,360) + 1e-5 > mod(lon_bounds(2) - min(lon_bathy) - 1e-5,360) + 1e-5))
    in_local_bathy_lon_range_ind_1 = find(lon_bathy - (mod(lon_bounds(1) - min(lon_bathy),360) + min(lon_bathy)) >= (-10*spacing));
    in_local_bathy_lon_range_ind_2 = find(lon_bathy - (mod(lon_bounds(2) - min(lon_bathy),360) + min(lon_bathy)) <= (10*spacing));
    local_lon_bathy = mod([lon_bathy(in_local_bathy_lon_range_ind_1); (lon_bathy(in_local_bathy_lon_range_ind_2) + 360)] - (lon_bounds(1) - (10*spacing)),360) + (lon_bounds(1) - (10*spacing));
    z_local = z_reshaped([in_local_bathy_lon_range_ind_1; in_local_bathy_lon_range_ind_2],in_local_bathy_lat_range_ind);
else
    in_local_bathy_lon_range_ind = find(((lon_bathy - (mod(lon_bounds(1) - min(lon_bathy),360) + min(lon_bathy))) >= (-10*spacing)) & (lon_bathy - (mod(lon_bounds(2) - min(lon_bathy),360) + min(lon_bathy)) <= (10*spacing)));
    local_lon_bathy = mod(lon_bathy(in_local_bathy_lon_range_ind) - (lon_bounds(1) - (10*spacing)),360) + (lon_bounds(1) - (10*spacing));
    z_local = z_reshaped(in_local_bathy_lon_range_ind,in_local_bathy_lat_range_ind);
end

% low-pass filter and interpolate bathymetry

% z_local = Gaussian_2D_filter(z_local,spacing,(2*0.1)/4,(2*0.1)/2,0,0,1);
% z_interp = interp2_fft(local_lon_bathy,local_lat_bathy,z_local,lon_z_interp,lat_z_interp);
% local_lon_bathy = lon_z_interp;
% local_lat_bathy = lat_z_interp;
% z_local = z_interp;



% load ECCO OBP

ECCO_nc_file = '/indopac/adelman/ECCO2/PHIBOT.ECCO2.lonlatinterp.1992-2018.nc';

longitude = ncread(ECCO_nc_file,'LONGITUDE_T');
latitude = ncread(ECCO_nc_file,'LATITUDE_T');
time = ncread(ECCO_nc_file,'TIME');
time = double(time) + datenum([1992 1 1 0 0 0]);

in_lon_range_ind = find((mod(longitude - (lon_bounds(1) - 4) + 180,360) - 180 >= 0) & (mod(longitude - (lon_bounds(2) + 4) + 180,360) - 180 <= 0));
in_lat_range_ind = find((latitude - (lat_bounds(1) - 3) >= 0) & (latitude - (lat_bounds(2) + 3) <= 0));
in_time_range_ind = find((time >= datenum(time_range_start)) & (time < datenum(time_range_end)));

if length(find(ismember([1 length(longitude)],in_lon_range_ind) == 1)) > 1
    gap_ind = find(diff(in_lon_range_ind) > 1.5);
    
    lon_in_range = longitude(in_lon_range_ind([((gap_ind + 1):1:length(in_lon_range_ind)) (1:1:gap_ind)]));
    lat_in_range = latitude(in_lat_range_ind);
    time_in_range = time(in_time_range_ind);
    
    start_vec = [in_lon_range_ind(gap_ind + 1) min(in_lat_range_ind) min(in_time_range_ind)];
    count_vec = [(max(in_lon_range_ind) - in_lon_range_ind(gap_ind + 1) + 1) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];
    
    OBP_ECCO = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
    OBP_ECCO = 1027.5*OBP_ECCO(in_lon_range_ind((gap_ind + 1):length(in_lon_range_ind)) - in_lon_range_ind(gap_ind + 1) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1);
    
    start_vec = [1 min(in_lat_range_ind) min(in_time_range_ind)];
    count_vec = [in_lon_range_ind(gap_ind) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];
    
    OBP_ECCO_2 = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
    OBP_ECCO = [OBP_ECCO; (1027.5*OBP_ECCO_2(in_lon_range_ind(1:gap_ind) - min(in_lon_range_ind) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1))];
else
    
    lon_in_range = longitude(in_lon_range_ind);
    lat_in_range = latitude(in_lat_range_ind);
    time_in_range = time(in_time_range_ind);
    
    start_vec = [min(in_lon_range_ind) min(in_lat_range_ind) min(in_time_range_ind)];
    count_vec = [(max(in_lon_range_ind) - min(in_lon_range_ind) + 1) (max(in_lat_range_ind) - min(in_lat_range_ind) + 1) (max(in_time_range_ind) - min(in_time_range_ind) + 1)];
    
    OBP_ECCO = ncread(ECCO_nc_file,'PHIBOT',start_vec,count_vec);
    OBP_ECCO = 1027.5*OBP_ECCO(in_lon_range_ind - min(in_lon_range_ind) + 1,in_lat_range_ind - min(in_lat_range_ind) + 1,in_time_range_ind - min(in_time_range_ind) + 1);
end

OBP_ECCO(OBP_ECCO < -1e26) = NaN;

diff_lon_in_range = diff(lon_in_range);
diff_lon_in_range = mod(diff_lon_in_range + 180,360) - 180;
lon_in_range = lon_in_range(1) + (360*round((lon_bounds(1) - lon_in_range(1))/360)) + [0; cumsum(diff_lon_in_range)];
diff_lat_in_range = diff(lat_in_range);

lon_in_range_bounds = [(lon_in_range(1) - (0.5*diff_lon_in_range(1))); (lon_in_range(2:length(lon_in_range)) - (diff_lon_in_range/2)); (lon_in_range(length(lon_in_range)) + (0.5*diff_lon_in_range(length(diff_lon_in_range))))];
lat_in_range_bounds = [(lat_in_range(1) - (0.5*diff_lat_in_range(1))); (lat_in_range(2:length(lat_in_range)) - (diff_lat_in_range/2)); (lat_in_range(length(lat_in_range)) + (0.5*diff_lat_in_range(length(diff_lat_in_range))))];


if abs(pt_ECCO_obs_opt) < 1e-5
    OBP_ECCO_zeronans = OBP_ECCO;
    OBP_ECCO_zeronans(isnan(OBP_ECCO) == 1) = 0;
    time_pt = time_in_range;
    OBP_pt = squeeze(interp2_fft(lon_in_range,lat_in_range,OBP_ECCO_zeronans,lon_pt,lat_pt));
    clear OBP_ECCO_zeronans
    frac_good_obs = ones(size(OBP_pt));
    frac_good_obs((isnan(OBP_pt) == 1) | (abs(OBP_pt) < 1e-5)) = 0;
elseif abs(pt_ECCO_obs_opt - 1) < 1e-5
    
    if strcmp(pt_id(1:5),'SAMOC') == 1
        % load bottom pressure observations (from PIES)

        curr_file = '/indopac/adelman/GRACE/SAMOC/SAM_PIES_data.txt';
        curr_fid = fopen(curr_file);

        fseek(curr_fid,931,'bof');
        data_array = (fscanf(curr_fid,'%d %d %d %f %f %f %f %f %f %f %f',[11 Inf]))';
        time_pt = datenum(data_array(:,1:3));
        if strcmp('SAMOC_A',pt_id) == 1
            OBP_pt = (1e4)*data_array(:,5);
        elseif strcmp('SAMOC_B',pt_id) == 1
            OBP_pt = (1e4)*data_array(:,7);
        elseif strcmp('SAMOC_C',pt_id) == 1
            OBP_pt = (1e4)*data_array(:,9);
        elseif strcmp('SAMOC_D',pt_id) == 1
            OBP_pt = (1e4)*data_array(:,11);
        end
        frac_good_obs = ones(size(OBP_pt));
        frac_good_obs(isnan(OBP_pt) == 1) = 0;
    elseif strcmp(pt_id(1:4),'DART') == 1
        curr_id = pt_id(6:10);
        curr_file = ['/indopac/adelman/GRACE/DART/DART_obp_',curr_id,'_all.nc'];

        lon_pt_obs = ncread(curr_file,'longitude');
        lat_pt_obs = ncread(curr_file,'latitude');
        time_obs = ncread(curr_file,'time') + 0.5;
        n_samples = ncread(curr_file,'n_samples');
        OBP_obs = (1e2)*ncread(curr_file,'obp_nodrift');
        OBP_obs(OBP_obs < -99000) = NaN;
        frac_good_obs = double((1/96)*n_samples);
        frac_good_obs(isnan(OBP_obs) == 1) = 0;
    end
end

size_array = size(OBP_ECCO);

if abs(monthavg_opt - 1) < 1e-5
    % compute monthly-averaged versions of time series
    
    time_datevec_ECCO = datevec(time_in_range);
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
        time_ECCO_monthavg(curr_month_ind) = sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*repmat(reshape(time_in_range(in_month_ind),[1 length(in_month_ind)]),[length(enough_ind) 1]),2),1)./(sum(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2),1));
        OBP_ECCO_monthavg(enough_ind,curr_month_ind) = sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind).*OBP_ECCO_reshaped_nans_zeroed(enough_ind,in_month_ind),2)./(sum(nan_mask_ECCO_reshaped(enough_ind,in_month_ind),2));    
    end
    
    time_in_range = time_ECCO_monthavg;
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
    
end


% define locations of GRACE mascons

mascon_lat_separation = 3;    % in degrees

curr_file = 'LAND_MASK.CRIv01.nc';
land_mask = ncread(curr_file,'land_mask');
curr_file = 'GRCTellus.JPL.200204_201706.GLO.RL06M.MSCNv01CRIv01.nc';
lon_GRACE = ncread(curr_file,'lon');
lat_GRACE = ncread(curr_file,'lat');
curr_lwe_thickness = ncread(curr_file,'lwe_thickness',[1 1 1],[length(lon_GRACE) length(lat_GRACE) 1]);

[lat_GRACE_grid,lon_GRACE_grid] = meshgrid(lat_GRACE,lon_GRACE);
ocean_mascon_curr_values = unique(curr_lwe_thickness(abs(land_mask) < 1e-5));

[lat_ECCO_grid,lon_ECCO_grid] = meshgrid(lat_in_range,lon_in_range);

mascon_lon_center = NaN([length(ocean_mascon_curr_values) 1]);
mascon_lat_center = NaN([length(ocean_mascon_curr_values) 1]);
mascon_lon_bounds = NaN([length(ocean_mascon_curr_values) 2]);
mascon_lat_bounds = NaN([length(ocean_mascon_curr_values) 2]);
OBP_ECCO_ocean_mascons = NaN([length(ocean_mascon_curr_values) size(OBP_ECCO,3)]);
mascon_depth_avg = NaN([length(ocean_mascon_curr_values) 1]);
for mascon_ind = 1:length(mascon_lon_center)
    curr_mascon_value = ocean_mascon_curr_values(mascon_ind);
    
    in_mascon_ind = find(abs(curr_lwe_thickness - curr_mascon_value) < 1e-10);
    lon_GRACE_grid_spacing = mean(mean(diff(lon_GRACE_grid,1,1)));
    lat_GRACE_grid_spacing = mean(mean(diff(lat_GRACE_grid,1,2)));
    if max(diff(sort(lon_GRACE_grid(in_mascon_ind),'ascend'))) > 2*lon_GRACE_grid_spacing
        lon_in_mascon_adj = mod(lon_GRACE_grid(in_mascon_ind) + 180,360) - 180;
        mascon_lon_center(mascon_ind) = mod(mean(lon_in_mascon_adj),360);
        mascon_lon_bounds(mascon_ind,:) = mod(min(lon_in_mascon_adj),360) + [(-0.5*lon_GRACE_grid_spacing) (mod(max(lon_in_mascon_adj) - min(lon_in_mascon_adj),360) + (0.5*lon_GRACE_grid_spacing))];
    else
        mascon_lon_center(mascon_ind) = mean(lon_GRACE_grid(in_mascon_ind));
        mascon_lon_bounds(mascon_ind,:) = [(min(lon_GRACE_grid(in_mascon_ind)) - (0.5*lon_GRACE_grid_spacing)) (max(lon_GRACE_grid(in_mascon_ind)) + (0.5*lon_GRACE_grid_spacing))];
        
    end
    in_mascon_lon_bathy_ind = find((mod(lon_bathy - mascon_lon_bounds(mascon_ind,1) + 180,360) - 180 > 0) & (mod(lon_bathy - mascon_lon_bounds(mascon_ind,2) + 180,360) - 180 < 0));
    mascon_lat_center(mascon_ind) = mean(lat_GRACE_grid(in_mascon_ind));
    mascon_lat_bounds(mascon_ind,:) = [(min(lat_GRACE_grid(in_mascon_ind)) - (0.5*lat_GRACE_grid_spacing)) (max(lat_GRACE_grid(in_mascon_ind)) + (0.5*lat_GRACE_grid_spacing))];
    in_mascon_lat_bathy_ind = find((lat_bathy - mascon_lat_bounds(mascon_ind,1) > 0) & (lat_bathy - mascon_lat_bounds(mascon_ind,2) < 0));
    
    if (mod(mascon_lon_bounds(mascon_ind,1) - min(lon_in_range_bounds) + 180,360) - 180 > 0) && (mod(mascon_lon_bounds(mascon_ind,2) - max(lon_in_range_bounds) + 180,360) - 180 < 0) && (mascon_lat_bounds(mascon_ind,1) - min(lat_in_range_bounds) > 0) && (mascon_lat_bounds(mascon_ind,2) - max(lat_in_range_bounds) < 0)
        % round region
        radius_cone = 1.5;   % in degrees latitude
        dist_from_center = abs((111100*cosd(mascon_lat_center(mascon_ind)).*(mod(lon_ECCO_grid - mascon_lon_center(mascon_ind) + 180,360) - 180)) + (1i*111100*(lat_ECCO_grid - mascon_lat_center(mascon_ind))));
        weight_matrix = 1 - (1*dist_from_center/(111100*radius_cone));
        weight_matrix(weight_matrix < 0) = 0;
        
        % for cylindrical (not conical) weight
        weight_matrix(weight_matrix > 1e-5) = 1;
        
        in_weight_range_ind = find(weight_matrix > 1e-5);
        in_mascon_lon_ECCO_ind = unique(mod(in_weight_range_ind - 1,size(lon_in_range,1)) + 1);
        in_mascon_lat_ECCO_ind = unique(ceil(in_weight_range_ind/size(lon_in_range,1)));
        weight_matrix = weight_matrix(in_mascon_lon_ECCO_ind,in_mascon_lat_ECCO_ind);
        % % 
        
%         % for box region
%         in_mascon_lon_ECCO_ind = find((mod(lon_in_range_ECCO - mascon_lon_bounds(mascon_ind,1) + 180,360) - 180 >= 0) & (mod(lon_in_range_ECCO - mascon_lon_bounds(mascon_ind,2) + 180,360) - 180 < 0));
%         in_mascon_lat_ECCO_ind = find((lat_in_range_ECCO - mascon_lat_bounds(mascon_ind,1) >= 0) & (lat_in_range_ECCO - mascon_lat_bounds(mascon_ind,2) < 0));
%         weight_matrix = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind)]);
%         % %         
        
%         in_mascon_lon_ECCO_ind = find((mod(lon_in_range_bounds(2:length(lon_in_range_bounds)) - mascon_lon_bounds(mascon_ind,1) + 180,360) - 180 > 0) & (mod(lon_in_range_bounds(1:length(lon_in_range)) - mascon_lon_bounds(mascon_ind,2) + 180,360) - 180 < 0));
        if max(diff(lon_in_range(in_mascon_lon_ECCO_ind))) > 3*median(diff(lon_in_range))
            gap_ind = find(diff(lon_in_range(in_mascon_lon_ECCO_ind)) > 3*median(diff(lon_in_range)));
            in_mascon_lon_width_ECCO = [diff(lon_in_range_bounds(in_mascon_lon_ECCO_ind((gap_ind + 1):length(in_mascon_lon_ECCO_ind)))); diff(lon_in_range_bounds(in_mascon_lon_ECCO_ind([(1:1:gap_ind)'; (mod(gap_ind + [-1; 0] - 1,length(in_mascon_lon_ECCO_ind)) + 1)])))];
            in_mascon_lon_ECCO_ind = in_mascon_lon_ECCO_ind([((gap_ind + 1):1:length(in_mascon_lon_ECCO_ind))'; (1:1:gap_ind)']);
        else
            in_mascon_lon_width_ECCO = diff(lon_in_range_bounds([in_mascon_lon_ECCO_ind; (max(in_mascon_lon_ECCO_ind) + 1)]));
        end
%         in_mascon_lat_ECCO_ind = find((lat_in_range_bounds(2:length(lat_in_range_bounds)) - mascon_lat_bounds(mascon_ind,1) > 0) & (lat_in_range_bounds(1:length(lat_in_range)) - mascon_lat_bounds(mascon_ind,2) < 0));
        in_mascon_lat_width_ECCO = diff(lat_in_range_bounds([in_mascon_lat_ECCO_ind; (max(in_mascon_lat_ECCO_ind) + 1)]));
        
        nan_mask_ECCO_in_mascon = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind) size(OBP_ECCO,3)]);
        OBP_ECCO_in_mascon_zeronans = OBP_ECCO(in_mascon_lon_ECCO_ind,in_mascon_lat_ECCO_ind,:);
        nan_mask_ECCO_in_mascon((isnan(OBP_ECCO_in_mascon_zeronans) == 1) | (abs(OBP_ECCO_in_mascon_zeronans) < 1e-15)) = 0;
        nan_mask_ECCO_not_enough = ones([length(in_mascon_lon_ECCO_ind) length(in_mascon_lat_ECCO_ind)]);
        nan_mask_ECCO_not_enough(sum(nan_mask_ECCO_in_mascon,3) < 0.9*max(max(sum(nan_mask_ECCO_in_mascon,3)))) = 0;
        nan_mask_ECCO_in_mascon(abs(repmat(nan_mask_ECCO_not_enough,[1 1 size(nan_mask_ECCO_in_mascon,3)])) < 1e-5) = 0;
        OBP_ECCO_in_mascon_zeronans(abs(nan_mask_ECCO_in_mascon) < 1e-5) = 0;
        
        weight_matrix = weight_matrix.*((111100*repmat(cosd(lat_in_range(in_mascon_lat_ECCO_ind)'),[length(in_mascon_lon_ECCO_ind) 1]).*repmat(in_mascon_lon_width_ECCO,[1 length(in_mascon_lat_ECCO_ind)])).*(111100*repmat(in_mascon_lat_width_ECCO',[length(in_mascon_lon_ECCO_ind) 1])));
        OBP_ECCO_ocean_mascons(mascon_ind,:) = reshape(sum(sum(repmat(weight_matrix,[1 1 size(OBP_ECCO,3)]).*nan_mask_ECCO_in_mascon.*OBP_ECCO_in_mascon_zeronans,2),1)./(sum(sum(repmat(weight_matrix,[1 1 size(OBP_ECCO,3)]).*nan_mask_ECCO_in_mascon,2),1)),[1 size(OBP_ECCO,3)]);
    end
    
    in_mascon_depth = -z_reshaped(in_mascon_lon_bathy_ind,in_mascon_lat_bathy_ind);
    mascon_depth_avg(mascon_ind) = mean(in_mascon_depth(in_mascon_depth > 0));
    
end

nan_mask_ocean_mascons = ones(size(OBP_ECCO_ocean_mascons));
nan_mask_ocean_mascons((isnan(OBP_ECCO_ocean_mascons) == 1) | (abs(OBP_ECCO_ocean_mascons) < 1e-15)) = 0;
mascon_ind_in_range = find(sum(nan_mask_ocean_mascons,2) > 0.9*max(sum(nan_mask_ocean_mascons,2)));
mascon_lon_center = mascon_lon_center(mascon_ind_in_range);
mascon_lat_center = mascon_lat_center(mascon_ind_in_range);
mascon_lon_bounds = mascon_lon_bounds(mascon_ind_in_range,:);
mascon_lat_bounds = mascon_lat_bounds(mascon_ind_in_range,:);
OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons(mascon_ind_in_range,:);
mascon_depth_avg = mascon_depth_avg(mascon_ind_in_range);

mascon_lon_center = mod(mascon_lon_center - (mean(lon_bounds) - 180),360) + (mean(lon_bounds) - 180);
mascon_lon_bounds = mod(mascon_lon_bounds - (mean(lon_bounds) - 180),360) + (mean(lon_bounds) - 180);


% temporally filter time series

steepness_factor = 5;
half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)


% if abs(reconstr_opt - 1) < 1e-5
    ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO));
    ECCO_nan_mask((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-10)) = -1;
    [OBP_ECCO_filtered,OBP_ECCO_trend,~] = bandpass_err_fcn(OBP_ECCO,3,mean(diff(time_in_range)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
    [OBP_ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,3,mean(diff(time_in_range)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

    OBP_ECCO = OBP_ECCO_filtered;
    clear OBP_ECCO_filtered
% end

ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO_ocean_mascons));
ECCO_nan_mask((isnan(OBP_ECCO_ocean_mascons) == 1) | (abs(OBP_ECCO_ocean_mascons) < 1e-10)) = -1;
[OBP_ECCO_ocean_mascons_filtered,OBP_ECCO_ocean_mascons_trend,~] = bandpass_err_fcn(OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
[ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,2,mean(diff(time_in_range)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons_filtered;
clear OBP_ECCO_ocean_mascons_filtered


curr_tseries_time = time_pt;
curr_tseries = OBP_pt;

% curr_tseries_nan_mask = (1e-5)*ones(size(curr_tseries));
% curr_tseries_nan_mask((isnan(curr_tseries) == 1) | (abs(curr_tseries) < 1e-10)) = -1;
% [curr_tseries_filtered,curr_tseries_trend,~] = bandpass_err_fcn(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
% [curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask,1,mean(diff(curr_tseries_time)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
% 
% curr_tseries = curr_tseries_filtered;
% clear curr_tseries_filtered

if abs(reconstr_opt - 1) < 1e-5
    nan_OBP_ECCO_ind = find(isnan(OBP_ECCO) == 1);
end
nan_ECCO_ind = find(isnan(OBP_ECCO_ocean_mascons) == 1);
% nan_curr_tseries_ind = find(isnan(curr_tseries) == 1);
if abs(season_cyc_opt) < 1e-5
    % remove annual cycle
    
    if abs(reconstr_opt - 1) < 1e-5
        OBP_ECCO(nan_OBP_ECCO_ind) = 0;
        G = [cos(((2*pi)/365.2425)*time_in_range) sin(((2*pi)/365.2425)*time_in_range) cos(((2*(2*pi))/365.2425)*time_in_range) sin(((2*(2*pi))/365.2425)*time_in_range) cos(((3*(2*pi))/365.2425)*time_in_range) sin(((3*(2*pi))/365.2425)*time_in_range) cos(((4*(2*pi))/365.2425)*time_in_range) sin(((4*(2*pi))/365.2425)*time_in_range)];
        coeffs = (((G')*G)^(-1))*((G')*(permute(reshape(OBP_ECCO,[(size(OBP_ECCO,1)*size(OBP_ECCO,2)) size(OBP_ECCO,3)]),[2 1])));
        OBP_ECCO = OBP_ECCO - reshape((G*coeffs)',size(OBP_ECCO));
    end
    
    OBP_ECCO_ocean_mascons(nan_ECCO_ind) = 0;
%     curr_tseries(nan_curr_tseries_ind) = 0;
    
    G = [cos(((2*pi)/365.2425)*time_in_range) sin(((2*pi)/365.2425)*time_in_range) cos(((2*(2*pi))/365.2425)*time_in_range) sin(((2*(2*pi))/365.2425)*time_in_range) cos(((3*(2*pi))/365.2425)*time_in_range) sin(((3*(2*pi))/365.2425)*time_in_range) cos(((4*(2*pi))/365.2425)*time_in_range) sin(((4*(2*pi))/365.2425)*time_in_range)];
    coeffs = (((G')*G)^(-1))*((G')*(OBP_ECCO_ocean_mascons'));
    OBP_ECCO_ocean_mascons = OBP_ECCO_ocean_mascons - reshape((G*coeffs)',size(OBP_ECCO_ocean_mascons));
%     G = [cos(((2*pi)/365.2425)*curr_tseries_time) sin(((2*pi)/365.2425)*curr_tseries_time) cos(((2*(2*pi))/365.2425)*curr_tseries_time) sin(((2*(2*pi))/365.2425)*curr_tseries_time) cos(((3*(2*pi))/365.2425)*curr_tseries_time) sin(((3*(2*pi))/365.2425)*curr_tseries_time) cos(((4*(2*pi))/365.2425)*curr_tseries_time) sin(((4*(2*pi))/365.2425)*curr_tseries_time)];
%     coeffs = (((G')*G)^(-1))*((G')*curr_tseries);
%     curr_tseries = curr_tseries - (G*coeffs);
end

% curr_tseries(curr_nan_ind) = NaN;

if abs(reconstr_opt - 1) < 1e-5
    OBP_ECCO(unique([nan_OBP_ECCO_ind; find(abs(OBP_ECCO_nan_mask_filtered) > 0.2)])) = NaN;
end
OBP_ECCO_ocean_mascons(unique([nan_ECCO_ind; find(abs(ECCO_nan_mask_filtered) > 0.2)])) = NaN;
% curr_tseries(unique([nan_curr_tseries_ind; find(abs(curr_tseries_nan_mask_filtered) > 0.2)])) = NaN;



% filter point time series

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

[curr_tseries_filtered,~,~] = bandpass_err_fcn(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
% [curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask + 1e-5,1,mean(diff(curr_tseries_time)),(1/(4*mean(diff(curr_tseries_time))*length(curr_tseries_time)))/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,1,1,edge_handling_opt,1);
filter_gain_coeffs_array = squeeze(bandpass_err_fcn_gain_coeffs(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1));


% estimate uncertainty due to missing obs., and propagate through time filter

sum_nan_mask = sum(curr_tseries_nan_mask,1);
% [~,curr_tseries_dof_zero_lag,~,~,~,~,curr_tseries_std_dev,~,~] = regression_linear_scalar_scalar(curr_tseries_filtered,curr_tseries_filtered,1,mean(diff(curr_tseries_time)),mean(diff(curr_tseries_time)),(1/5)*length(curr_tseries_time)*mean(diff(curr_tseries_time))*[0 1],0.95,0);
[~,curr_tseries_dof_zero_lag,~,~,~,~,curr_tseries_std_dev,~,lags] = regression_linear_scalar_scalar(curr_tseries_filtered,curr_tseries_filtered,1,mean(diff(curr_tseries_time)),(1/12)*365.24,(1/5)*length(curr_tseries_time)*mean(diff(curr_tseries_time))*[0 1],0.95,0);
curr_tseries_std_dev = curr_tseries_std_dev(abs(lags) < 1e-5);
curr_tseries_decorr_timescale = (mean(diff(curr_tseries_time))*sum_nan_mask)./curr_tseries_dof_zero_lag;
curr_tseries_missing_obs_err = NaN(size(curr_tseries_filtered));
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

curr_tseries = curr_tseries_filtered;
% curr_tseries_nan_mask = curr_tseries_nan_mask_filtered - 1e-5;
clear curr_tseries_filtered curr_tseries_nan_mask_filtered

% curr_tseries(curr_nan_ind) = NaN;

curr_tseries_missing_obs_norm_err = curr_tseries_missing_obs_err./repmat(curr_tseries_std_dev,[length(curr_tseries_missing_obs_err) 1]);

% curr_tseries(abs(curr_tseries_nan_mask) > 0.2) = NaN;
norm_err_tolerance = 0.3;
curr_tseries(abs(curr_tseries_missing_obs_norm_err) > norm_err_tolerance) = NaN;



% compute correlation of point time series with ECCO OBP

curr_tseries_tinterp = interp1(curr_tseries_time,curr_tseries,time_in_range);

[corr_tseries_ECCO_array,~,~,~,corr_tseries_ECCO_lowmag_bound_array,~,lags] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_tseries_tinterp,[1 length(curr_tseries_tinterp)]),[size(OBP_ECCO_ocean_mascons,1) 1]),OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95);

corr_tseries_ECCO_array(isnan(corr_tseries_ECCO_array) == 1) = 0;
corr_tseries_ECCO_lowmag_bound_array(isnan(corr_tseries_ECCO_lowmag_bound_array) == 1) = 0;
signal_noise_ratio_array = abs(corr_tseries_ECCO_array./(corr_tseries_ECCO_array - corr_tseries_ECCO_lowmag_bound_array));
signal_noise_ratio_array(isnan(signal_noise_ratio_array) == 1) = 0;

size_array_corr = size(corr_tseries_ECCO_array);
[~,max_corr_ind] = max(abs(corr_tseries_ECCO_array),[],2);

from_max_corr = abs(repmat(1:1:size_array_corr(2),[size_array_corr(1) 1]) - repmat(max_corr_ind,[1 size_array_corr(2)]));
max_corr_mask = zeros(size(corr_tseries_ECCO_array));
max_corr_mask(from_max_corr < 1e-5) = 1;
corr_tseries_ECCO_opt = sum(max_corr_mask.*corr_tseries_ECCO_array,2);
corr_tseries_ECCO_opt_lag = lags(max_corr_ind);
corr_tseries_ECCO_opt_SNR = sum(max_corr_mask.*signal_noise_ratio_array,2);


% create local grid to map mascons

lon_mascon_grid = unique(reshape(mascon_lon_bounds,[numel(mascon_lon_bounds) 1]));
lat_mascon_grid = unique(reshape(mascon_lat_bounds,[numel(mascon_lat_bounds) 1]));
corr_tseries_ECCO_mascon_grid_opt = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
corr_tseries_ECCO_mascon_grid_opt_lag = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
corr_tseries_ECCO_mascon_grid_opt_SNR = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
for mascon_ind = 1:size(OBP_ECCO_ocean_mascons,1)
    in_mascon_lon_ind = find((lon_mascon_grid - mascon_lon_bounds(mascon_ind,1) > -1e-5) & (lon_mascon_grid - mascon_lon_bounds(mascon_ind,2) < -1e-5));
    in_mascon_lat_ind = find((lat_mascon_grid - mascon_lat_bounds(mascon_ind,1) > -1e-5) & (lat_mascon_grid - mascon_lat_bounds(mascon_ind,2) < -1e-5));
    
    corr_tseries_ECCO_mascon_grid_opt(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_opt(mascon_ind);
    corr_tseries_ECCO_mascon_grid_opt_lag(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_opt_lag(mascon_ind);
    corr_tseries_ECCO_mascon_grid_opt_SNR(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_opt_SNR(mascon_ind);
end


% optimum correlation and lag plots

if abs(monthavg_opt - 1) < 1e-5
    ECCO_file_id = [ECCO_id,'_monthavg'];
else
    ECCO_file_id = ECCO_id;
end

if abs(reconstr_opt) < 1e-5
    % plot maps

    % in_plot_mascon_center_ind = find((mod(mascon_lon - lon_bounds(1) + 180,360) - 180 >= -0.5) & (mod(mascon_lon - lon_bounds(2) + 180,360) - 180 <= 0.5) & (mascon_lat - lat_bounds(1) >= -0.5) & (mascon_lat - lat_bounds(2) <= 0.5));

    c_levels = (-1):0.05:1;     % specified levels for contours
    cmap = colormap(0.85*bcyr(40));    % colormap for contours

    % plot pcolor map with land mask

    fig1000 = figure(1000);
    close(figure(1))
    colormap(cmap)
    fig_paper_pos = get(fig1000,'PaperPosition');
    fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
    fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
    set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
    p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_tseries_ECCO_mascon_grid_opt');
    caxis([min(c_levels) max(c_levels)])
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
    % marker_radius = 0.2;
    % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
    %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
    %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
    % end
    plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    hold off
    xlabel('Longitude','FontSize',12)
    ylabel('Latitude','FontSize',12)
    title({['Optimum correlation of ',pt_text_id,' OBP time series to ECCO ',ECCO_id]; ['OBP in mascons, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
    colormap(cmap)
    cbar = colorbar('location','southoutside');
    set(gca,'clim',[0 size(cmap,1)] - 0.5)
    set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',10)
    set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',10)

    print(fig1,['ECCO_map_corr_opt_',pt_id,'_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
    close(fig1)


    c_levels = [((2*lags(1)) - lags(2)) (lags(1:(length(lags) - 1)) + (diff(lags)/2))' ((2*lags(length(lags))) - lags(length(lags) - 1))];     % specified levels for contours
    cmap = colormap(0.85*bcyr(length(lags)));    % colormap for contours

    % plot pcolor map with land mask

    fig1000 = figure(1000);
    close(figure(1))
    colormap(cmap)
    fig_paper_pos = get(fig1000,'PaperPosition');
    fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
    fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
    set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
    p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_tseries_ECCO_mascon_grid_opt_lag');
    caxis([min(c_levels) max(c_levels)])
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

    fig2 = figure(2);
    close(figure(1))
    h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
    set(gca,'ydir','normal')
    daspect([1 cosd(mean(lat_bounds)) 1])
    hold on
    contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
    % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
    % marker_radius = 0.2;
    % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
    %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
    %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
    % end
    plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    hold off
    xlabel('Longitude','FontSize',12)
    ylabel('Latitude','FontSize',12)
    title({['Optimum lag of ECCO ',ECCO_id,' OBP in mascons relative to ',pt_text_id]; ['OBP time series, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
    colormap(cmap)
    cbar = colorbar('location','southoutside');
    set(gca,'clim',[0 size(cmap,1)] - 0.5)
    if abs(monthavg_opt - 1) < 1e-5
        curr_xtick_spacing = 1;
    else
        curr_xtick_spacing = 6;
    end
    curr_xtick_labels = cell(1,ceil(length(lags)/curr_xtick_spacing));
    for curr_label_ind = 1:length(curr_xtick_labels)
        curr_lag_ind = (curr_xtick_spacing*(curr_label_ind - 1)) + 1;
        if max(abs(lags)) >= 2*(365.24/12)
            curr_xtick_labels{curr_label_ind} = num2str(0.1*round((lags(curr_lag_ind)/(365.24/12))/0.1));
            set(get(cbar,'xlabel'),'String','Optimum lag (months)','FontSize',12)
        else
            curr_xtick_labels{curr_label_ind} = num2str(0.01*round(lags(curr_lag_ind)/0.01));
            set(get(cbar,'xlabel'),'String','Optimum lag (days)','FontSize',12)
        end
    end
    set(cbar,'xtick',(0:curr_xtick_spacing:size(cmap,1)),'xticklabel',curr_xtick_labels,'FontSize',12)

    print(fig2,['ECCO_map_corr_opt_lag_',pt_id,'_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
    close(fig2)


    % plot map for correlation at specific lag

    curr_lag_plot = 0;

    curr_lag_plot_ind = find(abs(lags - curr_lag_plot) < 1e-5);
    corr_tseries_ECCO_curr_lag = corr_tseries_ECCO_array(:,curr_lag_plot_ind);
    corr_tseries_ECCO_curr_lag_SNR = signal_noise_ratio_array(:,curr_lag_plot_ind);

    corr_tseries_ECCO_mascon_curr_lag = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
    corr_tseries_ECCO_mascon_curr_lag_SNR = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
    for mascon_ind = 1:size(OBP_ECCO_ocean_mascons,1)
        in_mascon_lon_ind = find((lon_mascon_grid - mascon_lon_bounds(mascon_ind,1) > -1e-5) & (lon_mascon_grid - mascon_lon_bounds(mascon_ind,2) < -1e-5));
        in_mascon_lat_ind = find((lat_mascon_grid - mascon_lat_bounds(mascon_ind,1) > -1e-5) & (lat_mascon_grid - mascon_lat_bounds(mascon_ind,2) < -1e-5));

        corr_tseries_ECCO_mascon_curr_lag(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_curr_lag(mascon_ind);
        corr_tseries_ECCO_mascon_curr_lag_SNR(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_curr_lag_SNR(mascon_ind);
    end


    c_levels = (-1):0.05:1;     % specified levels for contours
    cmap = colormap(0.85*bcyr(40));    % colormap for contours

    % plot pcolor map with land mask

    fig1000 = figure(1000);
    close(figure(1))
    colormap(cmap)
    fig_paper_pos = get(fig1000,'PaperPosition');
    fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
    fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
    set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
    p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_tseries_ECCO_mascon_curr_lag');
    caxis([min(c_levels) max(c_levels)])
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

    fig3 = figure(3);
    close(figure(1))
    h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
    set(gca,'ydir','normal')
    daspect([1 cosd(mean(lat_bounds)) 1])
    hold on
    contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
    % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
    % marker_radius = 0.2;
    % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
    %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
    %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
    % end
    plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    hold off
    xlabel('Longitude','FontSize',12)
    ylabel('Latitude','FontSize',12)
    title({['Correlation of ',pt_text_id,' OBP time series to ECCO ',ECCO_id,' OBP in mascons,']; ['at ',num2str(curr_lag_plot),' days lag, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
    colormap(cmap)
    cbar = colorbar('location','southoutside');
    set(gca,'clim',[0 size(cmap,1)] - 0.5)
    set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
    set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)

    print(fig3,['ECCO_map_corr_',num2str(curr_lag_plot),'_lag_',pt_id,'_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
    close(fig3)


%     % plot map of OBP standard deviation
% 
%     nan_mask_ECCO = ones(size(OBP_ECCO));
%     nan_mask_ECCO((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-15)) = 0;
%     OBP_ECCO_nans_zeroed = OBP_ECCO;
%     OBP_ECCO_nans_zeroed(isnan(OBP_ECCO) == 1) = 0;
% 
%     OBP_ECCO_mean = sum(nan_mask_ECCO.*OBP_ECCO_nans_zeroed,3)./(sum(nan_mask_ECCO,3));
%     OBP_ECCO_std_dev = (sum(nan_mask_ECCO.*((OBP_ECCO_nans_zeroed - repmat(OBP_ECCO_mean,[1 1 size_array(3)])).^2),3)./(sum(nan_mask_ECCO,3))).^(1/2);
% 
% 
%     c_levels = [(0:0.005:0.1) 1];     % specified levels for contours
%     cmap = colormap(0.85*flip(hot(length(c_levels) - 1),1));    % colormap for contours
% 
%     [irregular_clevels_array,~] = irregular_clevels_plot((1/(0.1*101325))*OBP_ECCO_std_dev,c_levels);
% 
%     fig4 = figure(4);
%     close(figure(1))
%     [~,cbar] = contour_plot_with_landmask(fig4,lon_in_range',lat_in_range,irregular_clevels_array',0:1:(length(c_levels) - 1),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
%     set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
%     hold on
%     contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
%     % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
%     hold off
%     xlabel('Longitude','FontSize',12)
%     ylabel('Latitude','FontSize',12)
%     title({'Standard deviation of ECCO OBP (dbar),'; [num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',10)
%     colormap(cmap)
%     set(gca,'clim',[0 size(cmap,1)] - 0.5)
%     set(cbar,'xtick',(0:4:size(cmap,1)) - 0.5,'xticklabel',{'0' '0.02' '0.04' '0.06' '0.08' '0.1'},'FontSize',12)
%     set(get(cbar,'xlabel'),'String','OBP standard deviation (dbar)','FontSize',12)
% 
%     print(fig4,['ECCO_map_std_dev_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
%     close(fig4)



%     % plot time series for pt time series and another point in ECCO
% 
%     ECCO_pt_lon = -51.5;
%     ECCO_pt_lat = -34.5;
% 
%     OBP_ECCO(isnan(OBP_ECCO) == 1) = 0;
%     curr_ECCO_tseries = squeeze(interp2_fft(lon_in_range,lat_in_range,OBP_ECCO,ECCO_pt_lon,ECCO_pt_lat));
%     ECCO_nan_mask_interp = squeeze(interp2_fft(lon_in_range,lat_in_range,OBP_ECCO_nan_mask_filtered,ECCO_pt_lon,ECCO_pt_lat));
%     curr_ECCO_tseries((abs(ECCO_nan_mask_interp) > 0.2) | (abs(curr_ECCO_tseries) < 1e-15)) = NaN;
% 
%     ECCO_tseries_text_id = [ECCO_id,' OBP at ',num2str(0.01*round(ECCO_pt_lon/0.01)),' ^{o} lon, ',num2str(0.01*round(ECCO_pt_lat/0.01)),' ^{o} lat'];
%     ECCO_text_id = [ECCO_id,' OBP, ',num2str(0.01*round(ECCO_pt_lon/0.01)),' ^{o} lon, ',num2str(0.01*round(ECCO_pt_lat/0.01)),' ^{o} lat'];
%     ECCO_tseries_id = [ECCO_file_id,'_',num2str(0.01*round(ECCO_pt_lon/0.01)),'_lon_',num2str(0.01*round(ECCO_pt_lat/0.01)),'_lat'];
% 
% 
%     fig5 = figure(5);
%     h = plot(curr_tseries_time,(1/(0.1*101325))*curr_tseries,time_in_range,(1/(0.1*101325))*curr_ECCO_tseries);
%     tick_year_spacing = 2;
%     years_to_plot_ticks = ((ceil(time_range_start(1)/tick_year_spacing)*tick_year_spacing):tick_year_spacing:(ceil((time_range_end(1))/tick_year_spacing)*tick_year_spacing))';
%     xtick_datenums_plot = datenum([years_to_plot_ticks ones(length(years_to_plot_ticks),2)]);
%     xtick_labels_plot = cell(length(years_to_plot_ticks),1);
%     for xtick_ind = 1:length(years_to_plot_ticks)
%         xtick_labels_plot{xtick_ind} = num2str(years_to_plot_ticks(xtick_ind));
%     end
%     set(gca,'xlim',[datenum(time_range_start) datenum(time_range_end)],'xtick',xtick_datenums_plot,'xticklabel',xtick_labels_plot,'xgrid','on')
%     set(h(1),'Color',[0.8 0 0],'LineWidth',2)
%     set(h(2),'Color',[0 0 0.8],'LineWidth',2)
%     hold on
%     line([datenum(time_range_start) datenum(time_range_end)],[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
%     hold off
%     ylabel('OBP anomaly (dbar)')
%     title({['Time series of OBP anomalies from ',pt_text_id,' and ECCO ',ECCO_tseries_text_id]; ['filtered for ',num2str(1/high_freq_bound),' to ',num2str(round(1/low_freq_bound)),' day periods']; ' '},'FontSize',8)
%     leg = legend(pt_text_id,ECCO_text_id,'location','northeast');
% 
%     saveas(fig5,['OBP_anom_time_series_',pt_id,'_',ECCO_tseries_id,'_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(datenum(time_range_start),'yyyymm'),'_',datestr(datenum(time_range_end) - 1,'yyyymm'),'_time_bounds.pdf'])
%     close(fig5)


%     % plot bathymetry map
% 
%     c_levels = 0:100:6000;     % specified levels for contours
%     cmap = flip(colormap(0.85*bcyr(length(c_levels) - 1)),1);    % colormap for contours
%     [irregular_clevels_array,~] = irregular_clevels_plot(-z_local,c_levels);
% 
%     contour_line_interval = 200;
% 
%     fig7 = figure(7);
%     close(figure(1))
%     [~,cbar] = contour_plot_with_landmask(fig7,local_lon_bathy',local_lat_bathy,irregular_clevels_array',0:1:length(c_levels),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
%     set(gca,'ydir','normal','clim',[0 size(cmap,1)],'FontSize',12)
%     hold on
%     % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
%     contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),contour_line_interval:contour_line_interval:10000,[0 0 0],'-',0.5)
%     hold off
%     xlabel('Longitude','FontSize',12)
%     ylabel('Latitude','FontSize',12)
%     title({['Bathymetry map, contour lines every ',num2str(contour_line_interval),' meters']; ' '},'FontSize',10)
%     colormap(cmap)
%     set(gca,'clim',[0 size(cmap,1)] - 0.5)
%     xlabel_spacing = 1000;
%     xlabel_ind = find(mod(c_levels,xlabel_spacing) == 0);
%     curr_xtick_labels = cell([1 length(xlabel_ind)]);
%     for curr_label_ind = 1:length(xlabel_ind)
%         curr_xlabel_ind = xlabel_ind(curr_label_ind);
%         curr_xtick_labels{curr_label_ind} = num2str(c_levels(curr_xlabel_ind));
%     end
%     set(cbar,'xtick',xlabel_ind - 1.5,'xticklabel',curr_xtick_labels,'FontSize',12)
%     set(get(cbar,'xlabel'),'String','Depth (meters)','FontSize',12)
% 
%     print(fig7,['Bathymetry_map_CI_',num2str(contour_line_interval),'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat.bmp'],'-dbmp','-r800')
%     close(fig7)
    
    
    keyboard    
end



% reconstruct point time series using objective analysis

% range of mascons to use (corners specified in this order: SW NW SE NE)
lon_corners = [-65 -50 -50 -35];
lat_corners = [-45 -25 -45 -25];

m_west = diff(lon_corners(1:2))/(diff(lat_corners(1:2)));
b_west = lon_corners(1) - (m_west*lat_corners(1));
m_east = diff(lon_corners(3:4))/(diff(lat_corners(3:4)));
b_east = lon_corners(3) - (m_east*lat_corners(3));

mascon_curr_ind = find((mascon_lon_center - ((m_west*mascon_lat_center) + b_west) >= 0) & (mascon_lon_center - ((m_east*mascon_lat_center) + b_east) <= 0));

m_south = diff(lat_corners([1 3]))/(diff(lon_corners([1 3])));
b_south = lat_corners(1) - (m_south*lon_corners(1));
m_north = diff(lat_corners([2 4]))/(diff(lon_corners([2 4])));
b_north = lat_corners(2) - (m_north*lon_corners(2));

mascon_curr_ind_2 = find((mascon_lon_center(mascon_curr_ind) - ((m_west*mascon_lat_center(mascon_curr_ind)) + b_west) >= 0) & (mascon_lon_center(mascon_curr_ind) - ((m_east*mascon_lat_center(mascon_curr_ind)) + b_east) <= 0));

mascon_in_reconstr_range_ind = mascon_curr_ind(mascon_curr_ind_2);


% normalized noise parameter
noise_param = 0;

OBP_ECCO_ocean_mascons_array_1 = repmat(reshape(OBP_ECCO_ocean_mascons(mascon_in_reconstr_range_ind,:),[length(mascon_in_reconstr_range_ind) 1 size(OBP_ECCO_ocean_mascons,2)]),[1 length(mascon_in_reconstr_range_ind) 1]);
OBP_ECCO_ocean_mascons_array_2 = repmat(reshape(OBP_ECCO_ocean_mascons(mascon_in_reconstr_range_ind,:),[1 length(mascon_in_reconstr_range_ind) size(OBP_ECCO_ocean_mascons,2)]),[length(mascon_in_reconstr_range_ind) 1 1]);


[OBP_ECCO_ocean_mascons_corr_array,~,~,~,~,~,lags_cov] = correlation_scalar_scalar_uncert_bounds(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95);
[~,~,~,~,~,~,std_dev_array_1,std_dev_array_2,~] = regression_linear_scalar_scalar(OBP_ECCO_ocean_mascons_array_1,OBP_ECCO_ocean_mascons_array_2,3,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95,1);
[~,~,~,~,~,~,~,std_dev_curr_tseries,~] = regression_linear_scalar_scalar(squeeze(OBP_ECCO_ocean_mascons_array_1(:,1,:)),repmat(reshape(curr_tseries_tinterp,[1 length(curr_tseries_tinterp)]),[length(mascon_in_reconstr_range_ind) 1]),2,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95,1);

% use zero lag for covariances
curr_lag = 0;
OBP_ECCO_ocean_mascons_corr_array = OBP_ECCO_ocean_mascons_corr_array(:,:,abs(lags_cov - curr_lag) < 1e-5);
std_dev_array_1 = std_dev_array_1(:,:,abs(lags_cov - curr_lag) < 1e-5);
std_dev_array_2 = std_dev_array_2(:,:,abs(lags_cov - curr_lag) < 1e-5);
std_dev_curr_tseries = std_dev_curr_tseries(:,abs(lags_cov - curr_lag) < 1e-5);
ocean_mascons_cov = OBP_ECCO_ocean_mascons_corr_array.*std_dev_array_1.*std_dev_array_2;
ocean_mascons_cov_diag = diag(ocean_mascons_cov,0);
ocean_mascons_cov = ocean_mascons_cov + (noise_param*diag(ocean_mascons_cov_diag,0));

ocean_mascons_tseries_cov = corr_tseries_ECCO_array(mascon_in_reconstr_range_ind,abs(lags - curr_lag) < 1e-5).*std_dev_array_1(:,1).*std_dev_curr_tseries;

gain_vec = (ocean_mascons_cov^(-1))*ocean_mascons_tseries_cov;
curr_tseries_reconstr = (OBP_ECCO_ocean_mascons(mascon_in_reconstr_range_ind,:)')*gain_vec;


% compute correlation of reconstructed time series with mascon OBP

% [corr_reconstr_ECCO_array,~,~,~,corr_reconstr_ECCO_lowmag_bound_array,~,lags] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_tseries_reconstr,[1 length(curr_tseries_reconstr)]),[size(OBP_ECCO_ocean_mascons,1) 1]),OBP_ECCO_ocean_mascons,2,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95);
[corr_reconstr_ECCO_array,~,~,~,corr_reconstr_ECCO_lowmag_bound_array,~,lags] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_tseries_reconstr,[1 1 length(curr_tseries_reconstr)]),[size(OBP_ECCO,1) size(OBP_ECCO,2) 1]),OBP_ECCO,3,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95);

corr_reconstr_ECCO_array(isnan(corr_reconstr_ECCO_array) == 1) = 0;
corr_reconstr_ECCO_lowmag_bound_array(isnan(corr_reconstr_ECCO_lowmag_bound_array) == 1) = 0;
signal_noise_ratio_array = abs(corr_reconstr_ECCO_array./(corr_reconstr_ECCO_array - corr_reconstr_ECCO_lowmag_bound_array));
signal_noise_ratio_array(isnan(signal_noise_ratio_array) == 1) = 0;

size_array_corr = size(corr_reconstr_ECCO_array);
% [max_corr,max_corr_ind] = max(abs(corr_reconstr_ECCO_array),[],2);
% 
% from_max_corr = abs(repmat(1:1:size_array_corr(2),[size_array_corr(1) 1]) - repmat(max_corr_ind,[1 size_array_corr(2)]));
% max_corr_mask = zeros(size(corr_reconstr_ECCO_array));
% max_corr_mask(from_max_corr < 1e-5) = 1;
% corr_reconstr_ECCO_opt = sum(max_corr_mask.*corr_reconstr_ECCO_array,2);
% corr_reconstr_ECCO_opt_lag = lags(max_corr_ind);
% corr_reconstr_ECCO_opt_SNR = sum(max_corr_mask.*signal_noise_ratio_array,2);

[max_corr,max_corr_ind] = max(abs(corr_reconstr_ECCO_array),[],3);

from_max_corr = abs(repmat(reshape(1:1:size_array_corr(3),[1 1 size_array_corr(3)]),[size_array_corr(1:2) 1]) - repmat(max_corr_ind,[1 1 size_array_corr(3)]));
max_corr_mask = zeros(size(corr_reconstr_ECCO_array));
max_corr_mask(from_max_corr < 1e-5) = 1;
corr_reconstr_ECCO_opt = sum(max_corr_mask.*corr_reconstr_ECCO_array,3);
corr_reconstr_ECCO_opt_lag = lags(max_corr_ind);
corr_reconstr_ECCO_opt_SNR = sum(max_corr_mask.*signal_noise_ratio_array,3);


% % create local grid to map mascons
% 
% % lon_mascon_grid = unique(reshape(mascon_lon_bounds,[numel(mascon_lon_bounds) 1]));
% % lat_mascon_grid = unique(reshape(mascon_lat_bounds,[numel(mascon_lat_bounds) 1]));
% corr_reconstr_ECCO_mascon_grid_opt = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
% corr_reconstr_ECCO_mascon_grid_opt_lag = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
% corr_reconstr_ECCO_mascon_grid_opt_SNR = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
% for mascon_ind = 1:size(OBP_ECCO_ocean_mascons,1)
%     in_mascon_lon_ind = find((lon_mascon_grid - mascon_lon_bounds(mascon_ind,1) > -1e-5) & (lon_mascon_grid - mascon_lon_bounds(mascon_ind,2) < -1e-5));
%     in_mascon_lat_ind = find((lat_mascon_grid - mascon_lat_bounds(mascon_ind,1) > -1e-5) & (lat_mascon_grid - mascon_lat_bounds(mascon_ind,2) < -1e-5));
%     
%     corr_reconstr_ECCO_mascon_grid_opt(in_mascon_lon_ind,in_mascon_lat_ind) = corr_reconstr_ECCO_opt(mascon_ind);
%     corr_reconstr_ECCO_mascon_grid_opt_lag(in_mascon_lon_ind,in_mascon_lat_ind) = corr_reconstr_ECCO_opt_lag(mascon_ind);
%     corr_reconstr_ECCO_mascon_grid_opt_SNR(in_mascon_lon_ind,in_mascon_lat_ind) = corr_reconstr_ECCO_opt_SNR(mascon_ind);
% end
% 
% 
% % plot maps
% 
% c_levels = (-1):0.05:1;     % specified levels for contours
% cmap = colormap(0.85*bcyr(40));    % colormap for contours
% 
% % plot pcolor map with land mask
% 
% fig1000 = figure(1000);
% close(figure(1))
% colormap(cmap)
% fig_paper_pos = get(fig1000,'PaperPosition');
% fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
% fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
% set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
% p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_reconstr_ECCO_mascon_grid_opt');
% caxis([min(c_levels) max(c_levels)])
% set(p,'edgecolor','none')
% set(gca,'xlim',lon_bounds,'ylim',lat_bounds)
% daspect([1 cosd(mean(lat_bounds)) 1])
% set(gca,'position',[0 0 1 1],'units','normalized','Visible','off')
% pause(5)
% plotboxaspectratio = get(gca,'PlotBoxAspectRatio');
% fig_pos = get(gcf,'Position');
% ratio_plotbox = plotboxaspectratio(2)/plotboxaspectratio(1);
% ratio_fig = fig_pos(4)/fig_pos(3);
% if ratio_fig > ratio_plotbox
%     fig_pos(4) = ratio_plotbox*fig_pos(3);
% else
%     fig_pos(3) = fig_pos(4)/ratio_plotbox;
% end
% set(gcf,'Position',fig_pos)
% 
% print(gcf,'contour_plot_temporary.png','-dpng','-r300')
% close(fig1000)
% 
% % get contour plot and overlay land mask
% 
% rgb_array = imread('contour_plot_temporary.png');
% rgb_array = flip(permute(rgb_array,[2 1 3]),2);
% 
% delete('contour_plot_temporary.png')
% 
% size_rgb_array = size(rgb_array);
% landmask_ind = landfind_indices(lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2),size_rgb_array(1:2));
% rgb_array_reshaped = reshape(rgb_array,[prod(size_rgb_array(1:2)) 3]);
% 
% black_shaded_ind = find((rgb_array_reshaped(:,1) == 0) & (rgb_array_reshaped(:,2) == 0) & (rgb_array_reshaped(:,3) == 0));
% rgb_array_reshaped(black_shaded_ind,:) = 255*ones(length(black_shaded_ind),3);
% 
% % put black mask on land areas
% rgb_array_reshaped(landmask_ind,:) = zeros(length(landmask_ind),3);
% 
% rgb_array_masked = reshape(rgb_array_reshaped,size_rgb_array);
% 
% fig1 = figure(1);
% h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
% set(gca,'ydir','normal')
% daspect([1 cosd(mean(lat_bounds)) 1])
% hold on
% contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% % marker_radius = 0.2;
% % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
% %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
% %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% % end
% plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
% hold off
% xlabel('Longitude','FontSize',12)
% ylabel('Latitude','FontSize',12)
% title({['Optimum correlation of ',pt_text_id,' reconstructed OBP time series to ECCO ',ECCO_id]; ['OBP in mascons, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',10)
% colormap(cmap)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)] - 0.5)
% set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',10)
% set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',10)
% 
% print(fig1,['ECCO_map_corr_opt_',pt_id,'_reconstr_',num2str(noise_param),'_noiseparam_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
% close(fig1)
% 
% 
% c_levels = [((2*lags(1)) - lags(2)) (lags(1:(length(lags) - 1)) + (diff(lags)/2))' ((2*lags(length(lags))) - lags(length(lags) - 1))];     % specified levels for contours
% cmap = colormap(0.85*bcyr(length(lags)));    % colormap for contours
% 
% % plot pcolor map with land mask
% 
% fig1000 = figure(1000);
% close(figure(1))
% colormap(cmap)
% fig_paper_pos = get(fig1000,'PaperPosition');
% fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
% fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
% set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
% p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_reconstr_ECCO_mascon_grid_opt_lag');
% caxis([min(c_levels) max(c_levels)])
% set(p,'edgecolor','none')
% set(gca,'xlim',lon_bounds,'ylim',lat_bounds)
% daspect([1 cosd(mean(lat_bounds)) 1])
% set(gca,'position',[0 0 1 1],'units','normalized','Visible','off')
% pause(5)
% plotboxaspectratio = get(gca,'PlotBoxAspectRatio');
% fig_pos = get(gcf,'Position');
% ratio_plotbox = plotboxaspectratio(2)/plotboxaspectratio(1);
% ratio_fig = fig_pos(4)/fig_pos(3);
% if ratio_fig > ratio_plotbox
%     fig_pos(4) = ratio_plotbox*fig_pos(3);
% else
%     fig_pos(3) = fig_pos(4)/ratio_plotbox;
% end
% set(gcf,'Position',fig_pos)
% 
% print(gcf,'contour_plot_temporary.png','-dpng','-r300')
% close(fig1000)
% 
% % get contour plot and overlay land mask
% 
% rgb_array = imread('contour_plot_temporary.png');
% rgb_array = flip(permute(rgb_array,[2 1 3]),2);
% 
% delete('contour_plot_temporary.png')
% 
% size_rgb_array = size(rgb_array);
% landmask_ind = landfind_indices(lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2),size_rgb_array(1:2));
% rgb_array_reshaped = reshape(rgb_array,[prod(size_rgb_array(1:2)) 3]);
% 
% black_shaded_ind = find((rgb_array_reshaped(:,1) == 0) & (rgb_array_reshaped(:,2) == 0) & (rgb_array_reshaped(:,3) == 0));
% rgb_array_reshaped(black_shaded_ind,:) = 255*ones(length(black_shaded_ind),3);
% 
% % put black mask on land areas
% rgb_array_reshaped(landmask_ind,:) = zeros(length(landmask_ind),3);
% 
% rgb_array_masked = reshape(rgb_array_reshaped,size_rgb_array);
% 
% fig2 = figure(2);
% close(figure(1))
% h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
% set(gca,'ydir','normal')
% daspect([1 cosd(mean(lat_bounds)) 1])
% hold on
% contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% % marker_radius = 0.2;
% % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
% %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
% %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% % end
% plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
% hold off
% xlabel('Longitude','FontSize',12)
% ylabel('Latitude','FontSize',12)
% title({['Optimum lag of ECCO ',ECCO_id,' OBP in mascons relative to reconstructed ',pt_text_id]; ['OBP time series, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',10)
% colormap(cmap)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)] - 0.5)
% if abs(monthavg_opt - 1) < 1e-5
%     curr_xtick_spacing = 1;
% else
%     curr_xtick_spacing = 6;
% end
% curr_xtick_labels = cell(1,ceil(length(lags)/curr_xtick_spacing));
% for curr_label_ind = 1:length(curr_xtick_labels)
%     curr_lag_ind = (curr_xtick_spacing*(curr_label_ind - 1)) + 1;
%     if max(abs(lags)) >= 2*(365.24/12)
%         curr_xtick_labels{curr_label_ind} = num2str(0.1*round((lags(curr_lag_ind)/(365.24/12))/0.1));
%         set(get(cbar,'xlabel'),'String','Optimum lag (months)','FontSize',12)
%     else
%         curr_xtick_labels{curr_label_ind} = num2str(0.01*round(lags(curr_lag_ind)/0.01));
%         set(get(cbar,'xlabel'),'String','Optimum lag (days)','FontSize',12)
%     end
% end
% set(cbar,'xtick',(0:curr_xtick_spacing:size(cmap,1)),'xticklabel',curr_xtick_labels,'FontSize',12)
% 
% print(fig2,['ECCO_map_corr_opt_lag_',pt_id,'_reconstr_',num2str(noise_param),'_noiseparam_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
% close(fig2)
% 
% 
% % plot map for correlation at specific lag
% 
% curr_lag_plot = 0;
% 
% curr_lag_plot_ind = find(abs(lags - curr_lag_plot) < 1e-5);
% corr_tseries_ECCO_curr_lag = corr_tseries_ECCO_array(:,curr_lag_plot_ind);
% corr_tseries_ECCO_curr_lag_SNR = signal_noise_ratio_array(:,curr_lag_plot_ind);
% 
% corr_tseries_ECCO_mascon_curr_lag = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
% corr_tseries_ECCO_mascon_curr_lag_SNR = NaN([length(lon_mascon_grid) length(lat_mascon_grid)]);
% for mascon_ind = 1:size(OBP_ECCO_ocean_mascons,1)
%     in_mascon_lon_ind = find((lon_mascon_grid - mascon_lon_bounds(mascon_ind,1) > -1e-5) & (lon_mascon_grid - mascon_lon_bounds(mascon_ind,2) < -1e-5));
%     in_mascon_lat_ind = find((lat_mascon_grid - mascon_lat_bounds(mascon_ind,1) > -1e-5) & (lat_mascon_grid - mascon_lat_bounds(mascon_ind,2) < -1e-5));
%     
%     corr_tseries_ECCO_mascon_curr_lag(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_curr_lag(mascon_ind);
%     corr_tseries_ECCO_mascon_curr_lag_SNR(in_mascon_lon_ind,in_mascon_lat_ind) = corr_tseries_ECCO_curr_lag_SNR(mascon_ind);
% end
% 
% 
% c_levels = (-1):0.05:1;     % specified levels for contours
% cmap = colormap(0.85*bcyr(40));    % colormap for contours
% 
% % plot pcolor map with land mask
% 
% fig1000 = figure(1000);
% close(figure(1))
% colormap(cmap)
% fig_paper_pos = get(fig1000,'PaperPosition');
% fig_paper_pos(4) = ((lat_bounds(2) - lat_bounds(1))/(cosd(mean(lat_bounds))*(lon_bounds(2) - lon_bounds(1))))*fig_paper_pos(3);
% fig_paper_pos(3:4) = 2*fig_paper_pos(3:4);
% set(fig1000,'PaperPosition',fig_paper_pos,'PaperSize',[22 17])
% p = pcolor(lon_mascon_grid,lat_mascon_grid',corr_reconstr_ECCO_mascon_curr_lag');
% caxis([min(c_levels) max(c_levels)])
% set(p,'edgecolor','none')
% set(gca,'xlim',lon_bounds,'ylim',lat_bounds)
% daspect([1 cosd(mean(lat_bounds)) 1])
% set(gca,'position',[0 0 1 1],'units','normalized','Visible','off')
% pause(5)
% plotboxaspectratio = get(gca,'PlotBoxAspectRatio');
% fig_pos = get(gcf,'Position');
% ratio_plotbox = plotboxaspectratio(2)/plotboxaspectratio(1);
% ratio_fig = fig_pos(4)/fig_pos(3);
% if ratio_fig > ratio_plotbox
%     fig_pos(4) = ratio_plotbox*fig_pos(3);
% else
%     fig_pos(3) = fig_pos(4)/ratio_plotbox;
% end
% set(gcf,'Position',fig_pos)
% 
% print(gcf,'contour_plot_temporary.png','-dpng','-r300')
% close(fig1000)
% 
% % get contour plot and overlay land mask
% 
% rgb_array = imread('contour_plot_temporary.png');
% rgb_array = flip(permute(rgb_array,[2 1 3]),2);
% 
% delete('contour_plot_temporary.png')
% 
% size_rgb_array = size(rgb_array);
% landmask_ind = landfind_indices(lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2),size_rgb_array(1:2));
% rgb_array_reshaped = reshape(rgb_array,[prod(size_rgb_array(1:2)) 3]);
% 
% black_shaded_ind = find((rgb_array_reshaped(:,1) == 0) & (rgb_array_reshaped(:,2) == 0) & (rgb_array_reshaped(:,3) == 0));
% rgb_array_reshaped(black_shaded_ind,:) = 255*ones(length(black_shaded_ind),3);
% 
% % put black mask on land areas
% rgb_array_reshaped(landmask_ind,:) = zeros(length(landmask_ind),3);
% 
% rgb_array_masked = reshape(rgb_array_reshaped,size_rgb_array);
% 
% fig3 = figure(3);
% close(figure(1))
% h = image(lon_bounds(1):((lon_bounds(2) - lon_bounds(1))/(size(rgb_array,1) - 1)):lon_bounds(2),(lat_bounds(1):((lat_bounds(2) - lat_bounds(1))/(size(rgb_array,2) - 1)):lat_bounds(2))',permute(rgb_array_masked,[2 1 3]));
% set(gca,'ydir','normal')
% daspect([1 cosd(mean(lat_bounds)) 1])
% hold on
% contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% % marker_radius = 0.2;
% % for curr_mascon_center = 1:length(in_plot_mascon_center_ind)
% %     curr_mascon_ind = in_plot_mascon_center_ind(curr_mascon_center);
% %     rectangle('Position',[(mascon_lon(curr_mascon_ind) + (360*round((mean(lon_bounds) - mascon_lon(curr_mascon_ind))/360)) - marker_radius) (mascon_lat(curr_mascon_ind) - marker_radius) (2*marker_radius) (2*marker_radius)],'Curvature',[1 1],'EdgeColor',[0 0 0],'FaceColor',[0.7 0 0.7]);
% % end
% plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
% hold off
% xlabel('Longitude','FontSize',12)
% ylabel('Latitude','FontSize',12)
% title({['Correlation of ',pt_text_id,' reconstructed OBP time series to ECCO ',ECCO_id,' OBP in mascons,']; ['at ',num2str(curr_lag_plot),' days lag, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',10)
% colormap(cmap)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)] - 0.5)
% set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
% set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)
% 
% print(fig3,['ECCO_map_corr_',num2str(curr_lag_plot),'_lag_',pt_id,'_reconstr_',num2str(noise_param),'_noiseparam_mascons_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
% close(fig3)



c_levels = (-1):0.05:1;     % specified levels for contours
cmap = colormap(0.85*bcyr(40));    % colormap for contours

fig1 = figure(1);
% close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig1,lon_in_range',lat_in_range,corr_reconstr_ECCO_opt',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% contour_overlay(double(lon_in_range'),double(lat_in_range),double(corr_reconstr_ECCO_opt_SNR'),[1 1],[0 0 0],'-',1)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
patch('xdata',lon_corners([1 2 4 3]),'ydata',lat_corners([1 2 4 3]),'FaceColor','none','EdgeColor',[0.4 0.4 0.4],'LineWidth',1.5)
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Optimum correlation of ',pt_text_id,' OBP reconstructed time series to ECCO']; ['OBP in region, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)

print(fig1,['ECCO_map_corr_opt_',pt_id,'_reconstr_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig1)


c_levels = [((2*lags(1)) - lags(2)) (lags(1:(length(lags) - 1)) + (diff(lags)/2))' ((2*lags(length(lags))) - lags(length(lags) - 1))];     % specified levels for contours
cmap = colormap(0.85*bcyr(length(lags)));    % colormap for contours

fig2 = figure(2);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig2,lon_in_range',lat_in_range',corr_reconstr_ECCO_opt_lag',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
patch('xdata',lon_corners([1 2 4 3]),'ydata',lat_corners([1 2 4 3]),'FaceColor','none','EdgeColor',[0.4 0.4 0.4],'LineWidth',1.5)
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Optimum lag of ECCO OBP in region relative to ',pt_text_id]; ['OBP reconstructed time series, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
if abs(monthavg_opt - 1) < 1e-5
    curr_xtick_spacing = 1;
else
    curr_xtick_spacing = 6;
end
curr_xtick_labels = cell(1,ceil(length(lags)/curr_xtick_spacing));
for curr_label_ind = 1:length(curr_xtick_labels)
    curr_lag_ind = (curr_xtick_spacing*(curr_label_ind - 1)) + 1;
    if max(abs(lags)) >= 2*(365.24/12)
        curr_xtick_labels{curr_label_ind} = num2str(0.1*round((lags(curr_lag_ind)/(365.24/12))/0.1));
        set(get(cbar,'xlabel'),'String','Optimum lag (months)','FontSize',12)
    else
        curr_xtick_labels{curr_label_ind} = num2str(0.01*round(lags(curr_lag_ind)/0.01));
        set(get(cbar,'xlabel'),'String','Optimum lag (days)','FontSize',12)
    end
end
set(cbar,'xtick',(0:curr_xtick_spacing:size(cmap,1)),'xticklabel',curr_xtick_labels,'FontSize',12)

print(fig2,['ECCO_map_corr_opt_lag_',pt_id,'_reconstr_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig2)


% plot map for correlation at specific lag

curr_lag_plot = 0;

curr_lag_plot_ind = find(abs(lags - curr_lag_plot) < 1e-5);
corr_reconstr_ECCO_curr_lag = corr_reconstr_ECCO_array(:,:,curr_lag_plot_ind);
corr_reconstr_ECCO_curr_lag_SNR = signal_noise_ratio_array(:,:,curr_lag_plot_ind);


c_levels = (-1):0.05:1;     % specified levels for contours
cmap = colormap(0.85*bcyr(40));    % colormap for contours

fig3 = figure(3);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig3,lon_in_range',lat_in_range,corr_reconstr_ECCO_curr_lag',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
contour_overlay(double(lon_in_range'),double(lat_in_range),double(corr_reconstr_ECCO_curr_lag_SNR'),[1 1],[0 0 0],'-',1)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
patch('xdata',lon_corners([1 2 4 3]),'ydata',lat_corners([1 2 4 3]),'FaceColor','none','EdgeColor',[0.4 0.4 0.4],'LineWidth',1.5)
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Correlation of ',pt_text_id,' OBP reconstructed time series to ECCO OBP in region,']; ['at ',num2str(curr_lag_plot),' days lag, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)

print(fig3,['ECCO_map_corr_',num2str(curr_lag_plot),'_lag_',pt_id,'_reconstr_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig3)




% plot time series comparison of original and reconstructed OBP


curr_tseries_tinterp(abs(curr_tseries_tinterp) < 1e-15) = NaN;
curr_tseries_reconstr(abs(curr_tseries_reconstr) < 1e-15) = NaN;


fig5 = figure(5);
h = plot(time_in_range,(1/(0.1*101325))*curr_tseries_tinterp,time_in_range,(1/(0.1*101325))*curr_tseries_reconstr);
tick_year_spacing = 2;
years_to_plot_ticks = ((ceil(time_range_start(1)/tick_year_spacing)*tick_year_spacing):tick_year_spacing:(ceil((time_range_end(1))/tick_year_spacing)*tick_year_spacing))';
xtick_datenums_plot = datenum([years_to_plot_ticks ones(length(years_to_plot_ticks),2)]);
xtick_labels_plot = cell(length(years_to_plot_ticks),1);
for xtick_ind = 1:length(years_to_plot_ticks)
    xtick_labels_plot{xtick_ind} = num2str(years_to_plot_ticks(xtick_ind));
end
set(gca,'xlim',[datenum(time_range_start) datenum(time_range_end)],'xtick',xtick_datenums_plot,'xticklabel',xtick_labels_plot,'xgrid','on')
set(h(1),'Color',[0.8 0 0],'LineWidth',1.5)
set(h(2),'Color',[0 0 0.8],'LineWidth',1.5)
hold on
line([datenum(time_range_start) datenum(time_range_end)],[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
hold off
ylabel('OBP anomaly (dbar)')
title({['Original and reconstructed time series of OBP anomaly from ',pt_text_id]; ['filtered for ',num2str(1/high_freq_bound),' to ',num2str(round(1/low_freq_bound)),' day periods']; ' '},'FontSize',10)
leg = legend(pt_text_id,[ECCO_id,' reconstruction'],'location','northeast');

saveas(fig5,['OBP_anom_time_series_',pt_id,'_reconstr_',ECCO_file_id,'_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(datenum(time_range_start),'yyyymm'),'_',datestr(datenum(time_range_end) - 1,'yyyymm'),'_time_bounds.pdf'])
close(fig5)