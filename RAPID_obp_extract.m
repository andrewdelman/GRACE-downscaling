% extract OBP data from specified DART bottom pressure recorder

path(path,'~/GRACE/')
path(path,'/indopac/adelman/GRACE/RAPID/')
cd('/indopac/adelman/GRACE/RAPID/')


curr_filenames = dir('/indopac/adelman/GRACE/RAPID/b*.nc');
filenames_cellarray = struct2cell(curr_filenames);
filenames_cellarray = filenames_cellarray(1,:);


longitude = NaN([length(filenames_cellarray) 1]);
latitude = NaN([length(filenames_cellarray) 1]);
depth = NaN([length(filenames_cellarray) 1]);
time_datenum_begin = NaN([length(filenames_cellarray) 1]);
time_datenum_end = NaN([length(filenames_cellarray) 1]);
time_data_length = NaN([length(filenames_cellarray) 1]);
for curr_filenum = 1:length(filenames_cellarray)
    curr_filename = filenames_cellarray{curr_filenum};
    
    longitude(curr_filenum) = ncread(curr_filename,'LONGITUDE');
    latitude(curr_filenum) = ncread(curr_filename,'LATITUDE');
    depth(curr_filenum) = ncread(curr_filename,'DEPTH');
    
    curr_time_vec = ncread(curr_filename,'TIME');
    curr_datenum_vec = datenum(datetime(curr_time_vec - 0.5,'convertfrom','juliandate'));
    time_datenum_begin(curr_filenum) = min(curr_datenum_vec);
    time_datenum_end(curr_filenum) = max(curr_datenum_vec);
    time_data_length(curr_filenum) = length(curr_datenum_vec);
end


cluster_radius = 0.02;      % radius (in deg lat) of points included in each cluster
cluster_depth_radius = 100;   % radius (in depth) of points included in eah cluster
cluster_first_ind = [];
in_cluster_ind_cellarray = cell([0 1]);
remaining_filenums = (1:1:length(filenames_cellarray))';
curr_filenum = 1;
curr_cluster = 1;
while curr_filenum <= length(remaining_filenums)
    if ismember(curr_filenum,remaining_filenums) == 1
        cluster_first_ind = [cluster_first_ind; curr_filenum];
        in_cluster_ind_cellarray{curr_cluster} = find((abs(((cosd(latitude(curr_filenum))).*(longitude - longitude(curr_filenum))) + (1i*(latitude - latitude(curr_filenum)))) <= cluster_radius) & (depth - depth(curr_filenum) <= cluster_depth_radius));
        remaining_filenums = setdiff(remaining_filenums,in_cluster_ind_cellarray{curr_cluster});
        curr_cluster = curr_cluster + 1;
    end
    curr_filenum = curr_filenum + 1;
end


lon_cluster = NaN([length(cluster_first_ind) 1]);
lat_cluster = NaN([length(cluster_first_ind) 1]);
depth_cluster = NaN([length(cluster_first_ind) 1]);
for curr_cluster_ind = 1:length(cluster_first_ind)
    in_curr_cluster_ind = in_cluster_ind_cellarray{curr_cluster_ind};
    
    lon_in_cluster = NaN([1 length(in_curr_cluster_ind)]);
    lat_in_cluster = NaN([1 length(in_curr_cluster_ind)]);
    depth_in_cluster = NaN([1 length(in_curr_cluster_ind)]);
    time_range_all_in_cluster = NaN([1 length(in_curr_cluster_ind)]);
    time_datenum_all_in_cluster = NaN([max(time_data_length(in_curr_cluster_ind)) length(in_curr_cluster_ind)]);
    pres_detrended_all_in_cluster = NaN([max(time_data_length(in_curr_cluster_ind)) length(in_curr_cluster_ind)]);
    for curr_tseries_ind = 1:length(in_curr_cluster_ind)
        curr_filenum = in_curr_cluster_ind(curr_tseries_ind);
        curr_filename = filenames_cellarray{curr_filenum};
        
        lon_in_cluster(curr_tseries_ind) = ncread(curr_filename,'LONGITUDE');
        lat_in_cluster(curr_tseries_ind) = ncread(curr_filename,'LATITUDE');
        depth_in_cluster(curr_tseries_ind) = ncread(curr_filename,'DEPTH');
        
        curr_time_vec = ncread(curr_filename,'TIME');
        curr_datenum_vec = datenum(datetime(curr_time_vec - 0.5,'convertfrom','juliandate'));
        time_datenum_all_in_cluster(1:length(curr_datenum_vec),curr_tseries_ind) = curr_datenum_vec;
        
        try
            pres_corrected_qcflag = ncread(curr_filename,'PRSTDR01_SEADATANET_QC');
            pres_corrected = ncread(curr_filename,'PRSTDR01');
        catch
            try
                pres_corrected_qcflag = ncread(curr_filename,'PRSTPS01_SEADATANET_QC');
                pres_corrected = ncread(curr_filename,'PRSTPS01');
            catch
                pres_corrected_qcflag = ncread(curr_filename,'PRSTRS01_SEADATANET_QC');
                pres_corrected = ncread(curr_filename,'PRSTRS01');   
            end
        end
        % remove first part of each time series (to try to exclude exponential drift)
        n_days_to_remove = 50;
        pres_corrected(curr_datenum_vec - min(curr_datenum_vec(curr_datenum_vec > 0.5)) < n_days_to_remove) = NaN;
        
        good_pres_ind = find(isnan(pres_corrected) == 0);
        [pres_detrended,~,~] = bandpass_err_fcn_no_pad(pres_corrected(good_pres_ind),1,1,1/(1.5*length(good_pres_ind)),1/(0.5*2),5,1,0,1,0);
        pres_detrended_vec = pres_corrected;
        pres_detrended_vec(good_pres_ind) = pres_detrended;
        pres_detrended_all_in_cluster(1:length(pres_corrected_qcflag),curr_tseries_ind) = pres_detrended_vec;
        time_range_all_in_cluster(curr_tseries_ind) = (mode(diff(sort(curr_datenum_vec(isnan(curr_datenum_vec) == 0),'ascend'))))*(length(good_pres_ind));
    end
    
    % averages of lon, lat, depth, weighted by time range in each time series
    lon_cluster(curr_cluster_ind) = sum(time_range_all_in_cluster.*(mod(lon_in_cluster - (lon_in_cluster(1) - 180),360) + (lon_in_cluster(1) - 180)))./(sum(time_range_all_in_cluster));
    lat_cluster(curr_cluster_ind) = sum(time_range_all_in_cluster.*lat_in_cluster)./(sum(time_range_all_in_cluster));
    depth_cluster(curr_cluster_ind) = sum(time_range_all_in_cluster.*depth_in_cluster)./(sum(time_range_all_in_cluster));
    
    time_datenum_all_in_cluster(isnan(pres_detrended_all_in_cluster) == 1) = NaN;
    sorted_time_datenum_in_cluster = sort(time_datenum_all_in_cluster,1,'ascend');
    good_sorted_ind = find(isnan(sorted_time_datenum_in_cluster) == 0);
    diff_time_datenum_in_cluster = diff(sorted_time_datenum_in_cluster(good_sorted_ind),1,1);
    mode_diff_cluster = mode((1/1440)*round(reshape(diff_time_datenum_in_cluster,[(length(good_sorted_ind) - 1) 1])/(1/1440)));
    
    mod_time_datenum_in_cluster = mod(sorted_time_datenum_in_cluster + (mode_diff_cluster/55),mode_diff_cluster) - (mode_diff_cluster/55);
    mode_mod_time_datenum = mode((1/1440)*round(reshape(mod_time_datenum_in_cluster(good_sorted_ind),[length(good_sorted_ind) 1])/(1/1440)));
    
    time_datenum_cluster_start = min(sorted_time_datenum_in_cluster(1,:)) + (mod(mode_mod_time_datenum - min(sorted_time_datenum_in_cluster(1,:)) + (mode_diff_cluster/2),mode_diff_cluster) - (mode_diff_cluster/2));
    time_datenum_cluster_end = max(max(sorted_time_datenum_in_cluster)) + (mod(mode_mod_time_datenum - max(max(sorted_time_datenum_in_cluster)) + (mode_diff_cluster/2),mode_diff_cluster) - (mode_diff_cluster/2));
    time_datenum_cluster = (time_datenum_cluster_start:mode_diff_cluster:time_datenum_cluster_end)';
    
    n_samples_cluster = NaN([length(time_datenum_cluster) 1]);
    obp_detrended_cluster = NaN([length(time_datenum_cluster) 1]);
    for curr_time_ind = 1:length(time_datenum_cluster)
        curr_time_center = time_datenum_cluster(curr_time_ind);
        in_time_bin_ind = find((time_datenum_all_in_cluster - curr_time_center >= (-mode_diff_cluster/2)) & (time_datenum_all_in_cluster - curr_time_center < (mode_diff_cluster/2)));
        
        good_in_time_bin_ind = find(isnan(pres_detrended_all_in_cluster(in_time_bin_ind)) == 0);
        good_in_time_bin_ind = in_time_bin_ind(good_in_time_bin_ind);
        n_samples_cluster(curr_time_ind) = length(good_in_time_bin_ind);
        time_datenum_cluster(curr_time_ind) = mean(time_datenum_all_in_cluster(good_in_time_bin_ind));
        obp_detrended_cluster(curr_time_ind) = mean(pres_detrended_all_in_cluster(good_in_time_bin_ind));
    end    
    
    
    
    % archive data in netCDF file
    
    output_filename = ['/indopac/adelman/GRACE/RAPID/RAPID_obp_first',num2str(n_days_to_remove),'removed_',num2str(0.01*round(lon_cluster(curr_cluster_ind)/0.01)),'_lon_',num2str(0.01*round(lat_cluster(curr_cluster_ind)/0.01)),'_lat_',num2str(round(depth_cluster(curr_cluster_ind))),'_depth.nc'];
    
    nccreate(output_filename,'longitude','Dimensions',{'longitude',1},'Datatype','double')
    ncwrite(output_filename,'longitude',lon_cluster(curr_cluster_ind))
    ncwriteatt(output_filename,'longitude','long_name','longitude of observation')
    ncwriteatt(output_filename,'longitude','units','degrees East')
    nccreate(output_filename,'latitude','Dimensions',{'latitude',1},'Datatype','double')
    ncwrite(output_filename,'latitude',lat_cluster(curr_cluster_ind))
    ncwriteatt(output_filename,'latitude','long_name','latitude of observation')
    ncwriteatt(output_filename,'latitude','units','degrees North')
    nccreate(output_filename,'depth','Dimensions',{'depth',1},'Datatype','double')
    ncwrite(output_filename,'depth',depth_cluster(curr_cluster_ind))
    ncwriteatt(output_filename,'depth','long_name','depth below surface')
    ncwriteatt(output_filename,'depth','units','meters')
    ncwriteatt(output_filename,'depth','positive','down')
    nccreate(output_filename,'time','Dimensions',{'time',Inf},'Datatype','double')
    ncwrite(output_filename,'time',time_datenum_cluster)
    ncwriteatt(output_filename,'time','long_name','time of observation')
    ncwriteatt(output_filename,'time','units','date number (day 1 = 0000-01-01 00:00:00)')
    nccreate(output_filename,'obp_nodrift','Dimensions',{'time',Inf},'Datatype','double')
    ncwrite(output_filename,'obp_nodrift',obp_detrended_cluster)
    ncwriteatt(output_filename,'obp_nodrift','long_name','Ocean bottom pressure minus tidal prediction, with corrections for drift')
    ncwriteatt(output_filename,'obp_nodrift','units','decibars')
    ncwriteatt(output_filename,'obp_nodrift','missing_value',-999)
    
    disp(['Processed cluster ',num2str(curr_cluster_ind),' (out of ',num2str(length(cluster_first_ind)),')'])
    
end
