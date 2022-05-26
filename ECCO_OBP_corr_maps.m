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

ECCO_id = 'cs510';   % ECCO simulation configuration used

lon_pt = -63.91;    % longitude of point time series
lat_pt = 23.40;     % latitude of point time series
lon_bounds = [-90 -40];    % longitude bounds for map plots
lat_bounds = [5 45];    % latitude bounds for map plots

pt_ECCO_obs_opt = 0;     % 0 = use point from ECCO, 1 = use point from obs.
pt_text_id = 'ECCO cs510 -63.91 lon, 23.40 lat';
pt_id = 'ECCO_cs510_-63.91_lon_23.40_lat_5850_depth';

time_range_start = [1992 1 1];
time_range_end = [2019 1 1];

low_freq_bound = 1/365.24;
high_freq_bound = 1/((3/12)*365.24);
season_cyc_opt = 0;    % 0 = remove seasonal/annual cycle, 1 = retain seasonal/annual cycle
monthavg_opt = 0;     % 0 = no monthly averaging, 1 = monthly averaging

% lag parameters for correlations
% delta_lag = (1/72)*365.24;
delta_lag = 2;
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

in_lon_range_ind = find((mod(longitude - (lon_bounds(1) - 1) + 180,360) - 180 >= 0) & (mod(longitude - (lon_bounds(2) + 1) + 180,360) - 180 <= 0));
in_lat_range_ind = find((latitude - lat_bounds(1) >= 0) & (latitude - lat_bounds(2) <= 0));
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


if abs(pt_ECCO_obs_opt) < 1e-5
    OBP_ECCO_zeronans = OBP_ECCO;
    OBP_ECCO_zeronans(isnan(OBP_ECCO) == 1) = 0;
    time_pt = time_in_range;
    OBP_pt = squeeze(interp2_fft(lon_in_range,lat_in_range,OBP_ECCO_zeronans,lon_pt,lat_pt));
    clear OBP_ECCO_zeronans
elseif abs(pt_ECCO_obs_opt - 1) < 1e-5
    
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


% temporally filter time series

steepness_factor = 5;
half_power_adj = exp(erfinv((2^(1/2)) - 1)/steepness_factor);   % adjustment factor to set bounds at half-power (rather than half-amplitude)


ECCO_nan_mask = (1e-5)*ones(size(OBP_ECCO));
ECCO_nan_mask((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-10)) = -1;
[OBP_ECCO_filtered,OBP_ECCO_trend,~] = bandpass_err_fcn(OBP_ECCO,3,mean(diff(time_in_range)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
[ECCO_nan_mask_filtered,~,~] = bandpass_err_fcn(ECCO_nan_mask,3,mean(diff(time_in_range)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

OBP_ECCO = OBP_ECCO_filtered;
clear OBP_ECCO_filtered


curr_tseries_time = time_pt;
curr_tseries = OBP_pt;

curr_tseries_nan_mask = (1e-5)*ones(size(curr_tseries));
curr_tseries_nan_mask((isnan(curr_tseries) == 1) | (abs(curr_tseries) < 1e-10)) = -1;
[curr_tseries_filtered,curr_tseries_trend,~] = bandpass_err_fcn(curr_tseries,1,mean(diff(curr_tseries_time)),low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);
[curr_tseries_nan_mask_filtered,~,~] = bandpass_err_fcn(curr_tseries_nan_mask,1,mean(diff(curr_tseries_time)),0.5*low_freq_bound/half_power_adj,high_freq_bound*half_power_adj,steepness_factor,0,1,1,1);

curr_tseries = curr_tseries_filtered;
clear curr_tseries_filtered


nan_ECCO_ind = find(isnan(OBP_ECCO) == 1);
nan_curr_tseries_ind = find(isnan(curr_tseries) == 1);
if abs(season_cyc_opt) < 1e-5
    % remove annual cycle

    OBP_ECCO(nan_ECCO_ind) = 0;
    curr_tseries(nan_curr_tseries_ind) = 0;
    
    G = [cos(((2*pi)/365.2425)*time_in_range) sin(((2*pi)/365.2425)*time_in_range) cos(((2*(2*pi))/365.2425)*time_in_range) sin(((2*(2*pi))/365.2425)*time_in_range) cos(((3*(2*pi))/365.2425)*time_in_range) sin(((3*(2*pi))/365.2425)*time_in_range) cos(((4*(2*pi))/365.2425)*time_in_range) sin(((4*(2*pi))/365.2425)*time_in_range)];
    coeffs = (((G')*G)^(-1))*((G')*(reshape(OBP_ECCO,[prod(size_array(1:2)) size_array(3)])'));
    OBP_ECCO = OBP_ECCO - reshape((G*coeffs)',size_array);
    G = [cos(((2*pi)/365.2425)*curr_tseries_time) sin(((2*pi)/365.2425)*curr_tseries_time) cos(((2*(2*pi))/365.2425)*curr_tseries_time) sin(((2*(2*pi))/365.2425)*curr_tseries_time) cos(((3*(2*pi))/365.2425)*curr_tseries_time) sin(((3*(2*pi))/365.2425)*curr_tseries_time) cos(((4*(2*pi))/365.2425)*curr_tseries_time) sin(((4*(2*pi))/365.2425)*curr_tseries_time)];
    coeffs = (((G')*G)^(-1))*((G')*curr_tseries);
    curr_tseries = curr_tseries - (G*coeffs);
end

% curr_tseries(curr_nan_ind) = NaN;

OBP_ECCO(unique([nan_ECCO_ind; find(abs(ECCO_nan_mask_filtered) > 0.2)])) = NaN;
curr_tseries(unique([nan_curr_tseries_ind; find(abs(curr_tseries_nan_mask_filtered) > 0.2)])) = NaN;



% compute correlation of point time series with ECCO OBP

curr_tseries_tinterp = interp1(curr_tseries_time,curr_tseries,time_in_range);

[corr_tseries_ECCO_array,~,~,~,corr_tseries_ECCO_lowmag_bound_array,~,lags] = correlation_scalar_scalar_uncert_bounds(repmat(reshape(curr_tseries_tinterp,[1 1 length(curr_tseries_tinterp)]),[size_array(1:2) 1]),OBP_ECCO,3,mean(diff(time_in_range)),max([delta_lag (mean(diff(time_in_range)) - 1e-2)]),lag_range_to_test,0.95);

corr_tseries_ECCO_array(isnan(corr_tseries_ECCO_array) == 1) = 0;
corr_tseries_ECCO_lowmag_bound_array(isnan(corr_tseries_ECCO_lowmag_bound_array) == 1) = 0;
signal_noise_ratio_array = abs(corr_tseries_ECCO_array./(corr_tseries_ECCO_array - corr_tseries_ECCO_lowmag_bound_array));
signal_noise_ratio_array(isnan(signal_noise_ratio_array) == 1) = 0;

size_array_corr = size(corr_tseries_ECCO_array);
corr_tseries_ECCO_reshaped = reshape(corr_tseries_ECCO_array,[prod(size_array_corr(1:2)) size_array_corr(3)]);
signal_noise_ratio_reshaped = reshape(signal_noise_ratio_array,[prod(size_array_corr(1:2)) size_array_corr(3)]);
[max_corr,max_corr_ind] = max(abs(corr_tseries_ECCO_reshaped),[],2);

from_max_corr = abs(repmat(1:1:size_array_corr(3),[prod(size_array_corr(1:2)) 1]) - repmat(max_corr_ind,[1 size_array_corr(3)]));
max_corr_mask = zeros(size(corr_tseries_ECCO_reshaped));
max_corr_mask(from_max_corr < 1e-5) = 1;
corr_tseries_ECCO_opt = reshape(sum(max_corr_mask.*corr_tseries_ECCO_reshaped,2),size_array_corr(1:2));
corr_tseries_ECCO_opt_lag = reshape(lags(max_corr_ind),size_array_corr(1:2));
corr_tseries_ECCO_opt_SNR = reshape(sum(max_corr_mask.*signal_noise_ratio_reshaped,2),size_array_corr(1:2));



% optimum correlation and lag plots

if abs(monthavg_opt - 1) < 1e-5
    ECCO_file_id = [ECCO_id,'_monthavg'];
else
    ECCO_file_id = ECCO_id;
end


c_levels = (-1):0.05:1;     % specified levels for contours
cmap = colormap(0.85*bcyr(40));    % colormap for contours

fig1 = figure(1);
% close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig1,lon_in_range',lat_in_range,corr_tseries_ECCO_opt',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% contour_overlay(double(lon_in_range'),double(lat_in_range),double(corr_tseries_ECCO_opt_SNR'),[1 1],[0 0 0],'-',1)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Optimum correlation of ',pt_text_id,' OBP time series to ECCO ',ECCO_id]; ['OBP in region, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)

print(fig1,['ECCO_map_corr_opt_',pt_id,'_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig1)


c_levels = [((2*lags(1)) - lags(2)) (lags(1:(length(lags) - 1)) + (diff(lags)/2))' ((2*lags(length(lags))) - lags(length(lags) - 1))];     % specified levels for contours
cmap = colormap(0.85*bcyr(length(lags)));    % colormap for contours

fig2 = figure(2);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig2,lon_in_range',lat_in_range',corr_tseries_ECCO_opt_lag',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Optimum lag of ECCO ',ECCO_id,' OBP in region relative to ',pt_text_id]; ['OBP time series, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
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

print(fig2,['ECCO_map_corr_opt_lag_',pt_id,'_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig2)


% plot map for correlation at specific lag

curr_lag_plot = 0;

curr_lag_plot_ind = find(abs(lags - curr_lag_plot) < 1e-5);
corr_tseries_ECCO_curr_lag = corr_tseries_ECCO_array(:,:,curr_lag_plot_ind);
corr_tseries_ECCO_curr_lag_SNR = signal_noise_ratio_array(:,:,curr_lag_plot_ind);


c_levels = (-1):0.05:1;     % specified levels for contours
cmap = colormap(0.85*bcyr(40));    % colormap for contours

fig3 = figure(3);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig3,lon_in_range',lat_in_range,corr_tseries_ECCO_curr_lag',c_levels,cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
contour_overlay(double(lon_in_range'),double(lat_in_range),double(corr_tseries_ECCO_curr_lag_SNR'),[1 1],[0 0 0],'-',1)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
plot(lon_pt,lat_pt,'p','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({['Correlation of ',pt_text_id,' OBP time series to ECCO ',ECCO_id,' OBP in region,']; ['at ',num2str(curr_lag_plot),' days lag, ',num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
set(cbar,'xtick',(0:5:size(cmap,1)) - 0.5,'xticklabel',{'-1' '' '-0.5' '' '0' '' '0.5' '' '1'},'FontSize',12)
set(get(cbar,'xlabel'),'String','Correlation coefficient','FontSize',12)

print(fig3,['ECCO_map_corr_',num2str(curr_lag_plot),'_lag_',pt_id,'_',ECCO_file_id,'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig3)


% plot map of OBP standard deviation

nan_mask_ECCO = ones(size(OBP_ECCO));
nan_mask_ECCO((isnan(OBP_ECCO) == 1) | (abs(OBP_ECCO) < 1e-15)) = 0;
OBP_ECCO_nans_zeroed = OBP_ECCO;
OBP_ECCO_nans_zeroed(isnan(OBP_ECCO) == 1) = 0;

OBP_ECCO_mean = sum(nan_mask_ECCO.*OBP_ECCO_nans_zeroed,3)./(sum(nan_mask_ECCO,3));
OBP_ECCO_std_dev = (sum(nan_mask_ECCO.*((OBP_ECCO_nans_zeroed - repmat(OBP_ECCO_mean,[1 1 size_array(3)])).^2),3)./(sum(nan_mask_ECCO,3))).^(1/2);


c_levels = [(0:0.005:0.1) 1];     % specified levels for contours
cmap = colormap(0.85*flip(hot(length(c_levels) - 1),1));    % colormap for contours

[irregular_clevels_array,~] = irregular_clevels_plot((1/(0.1*101325))*OBP_ECCO_std_dev,c_levels);

fig4 = figure(4);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig4,lon_in_range',lat_in_range,irregular_clevels_array',0:1:(length(c_levels) - 1),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)] - 0.5,'FontSize',12)
hold on
contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),1000:1000:10000,[0.5 0.3 0.2],'-',0.5)
% line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
title({'Standard deviation of ECCO OBP (dbar),'; [num2str(round(1/high_freq_bound)),' to ',num2str(round(1/low_freq_bound)),' day periods, bathymetry CI = 1000 m']; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
set(cbar,'xtick',(0:4:size(cmap,1)) - 0.5,'xticklabel',{'0' '0.02' '0.04' '0.06' '0.08' '0.1'},'FontSize',12)
set(get(cbar,'xlabel'),'String','OBP standard deviation (dbar)','FontSize',12)

print(fig4,['ECCO_map_std_dev_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(min(time_in_range),'yyyymm'),'_',datestr(max(time_in_range),'yyyymm'),'_time_bounds.bmp'],'-dbmp','-r300')
close(fig4)



% plot time series for pt time series and another point in ECCO

ECCO_pt_lon = -51.5;
ECCO_pt_lat = -34.5;

OBP_ECCO(isnan(OBP_ECCO) == 1) = 0;
curr_ECCO_tseries = squeeze(interp2_fft(lon_in_range,lat_in_range,OBP_ECCO,ECCO_pt_lon,ECCO_pt_lat));
ECCO_nan_mask_interp = squeeze(interp2_fft(lon_in_range,lat_in_range,ECCO_nan_mask_filtered,ECCO_pt_lon,ECCO_pt_lat));
curr_ECCO_tseries((abs(ECCO_nan_mask_interp) > 0.2) | (abs(curr_ECCO_tseries) < 1e-15)) = NaN;

ECCO_tseries_text_id = [ECCO_id,' OBP at ',num2str(0.01*round(ECCO_pt_lon/0.01)),' ^{o} lon, ',num2str(0.01*round(ECCO_pt_lat/0.01)),' ^{o} lat'];
ECCO_text_id = [ECCO_id,' OBP, ',num2str(0.01*round(ECCO_pt_lon/0.01)),' ^{o} lon, ',num2str(0.01*round(ECCO_pt_lat/0.01)),' ^{o} lat'];
ECCO_tseries_id = [ECCO_file_id,'_',num2str(0.01*round(ECCO_pt_lon/0.01)),'_lon_',num2str(0.01*round(ECCO_pt_lat/0.01)),'_lat'];


fig5 = figure(5);
h = plot(curr_tseries_time,(1/(0.1*101325))*curr_tseries,time_in_range,(1/(0.1*101325))*curr_ECCO_tseries);
tick_year_spacing = 2;
years_to_plot_ticks = ((ceil(time_range_start(1)/tick_year_spacing)*tick_year_spacing):tick_year_spacing:(ceil((time_range_end(1))/tick_year_spacing)*tick_year_spacing))';
xtick_datenums_plot = datenum([years_to_plot_ticks ones(length(years_to_plot_ticks),2)]);
xtick_labels_plot = cell(length(years_to_plot_ticks),1);
for xtick_ind = 1:length(years_to_plot_ticks)
    xtick_labels_plot{xtick_ind} = num2str(years_to_plot_ticks(xtick_ind));
end
set(gca,'xlim',[datenum(time_range_start) datenum(time_range_end)],'xtick',xtick_datenums_plot,'xticklabel',xtick_labels_plot,'xgrid','on')
set(h(1),'Color',[0.8 0 0],'LineWidth',2)
set(h(2),'Color',[0 0 0.8],'LineWidth',2)
hold on
line([datenum(time_range_start) datenum(time_range_end)],[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
hold off
ylabel('OBP anomaly (dbar)')
title({['Time series of OBP anomalies from ',pt_text_id,' and ECCO ',ECCO_tseries_text_id]; ['filtered for ',num2str(1/high_freq_bound),' to ',num2str(round(1/low_freq_bound)),' day periods']; ' '},'FontSize',8)
leg = legend(pt_text_id,ECCO_text_id,'location','northeast');

saveas(fig5,['OBP_anom_time_series_',pt_id,'_',ECCO_tseries_id,'_',num2str(round(1/high_freq_bound)),'_',num2str(round(1/low_freq_bound)),'_periods_',datestr(datenum(time_range_start),'yyyymm'),'_',datestr(datenum(time_range_end) - 1,'yyyymm'),'_time_bounds.pdf'])
close(fig5)


% plot bathymetry map

c_levels = 0:100:6000;     % specified levels for contours
cmap = flip(colormap(0.85*bcyr(length(c_levels) - 1)),1);    % colormap for contours
[irregular_clevels_array,~] = irregular_clevels_plot(-z_local,c_levels);

contour_line_interval = 200;

fig7 = figure(7);
close(figure(1))
[~,cbar] = contour_plot_with_landmask(fig7,local_lon_bathy',local_lat_bathy,irregular_clevels_array',0:1:length(c_levels),cmap,'southoutside',lon_bounds(1),lon_bounds(2),lat_bounds(1),lat_bounds(2));
set(gca,'ydir','normal','clim',[0 size(cmap,1)],'FontSize',12)
hold on
% % line(lon_bounds,[0 0],'Color','k','LineStyle','-','LineWidth',1)
% contour_overlay(double(local_lon_bathy'),double(local_lat_bathy),double((-z_local)'),contour_line_interval:contour_line_interval:10000,[0 0 0],'-',0.5)
hold off
xlabel('Longitude','FontSize',12)
ylabel('Latitude','FontSize',12)
% title({['Bathymetry map, contour lines every ',num2str(contour_line_interval),' meters']; ' '},'FontSize',8)
title({'Bathymetry map'; ' '},'FontSize',8)
colormap(cmap)
set(gca,'clim',[0 size(cmap,1)] - 0.5)
xlabel_spacing = 1000;
xlabel_ind = find(mod(c_levels,xlabel_spacing) == 0);
curr_xtick_labels = cell([1 length(xlabel_ind)]);
for curr_label_ind = 1:length(xlabel_ind)
    curr_xlabel_ind = xlabel_ind(curr_label_ind);
    curr_xtick_labels{curr_label_ind} = num2str(c_levels(curr_xlabel_ind));
end
set(cbar,'xtick',xlabel_ind - 1.5,'xticklabel',curr_xtick_labels,'FontSize',12)
set(get(cbar,'xlabel'),'String','Depth (meters)','FontSize',12)

% print(fig7,['Bathymetry_map_CI_',num2str(contour_line_interval),'_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat.bmp'],'-dbmp','-r800')
print(fig7,['Bathymetry_map_',num2str(lon_bounds(1)),'_',num2str(lon_bounds(2)),'_lon_',num2str(lat_bounds(1)),'_',num2str(lat_bounds(2)),'_lat.bmp'],'-dbmp','-r800')
close(fig7)