% compare distributions of DART-downscale estimate correlations

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
cd('/indopac/adelman/GRACE/')


% load correlation values
lowperiod = 91;
highperiod = 365;
load(['DARTplus_downscale_corrcoeff_',num2str(lowperiod),'_',num2str(highperiod),'_periods.mat'])

case1_compare = 'coloc_mascons';
case1_name = 'Not downscaled';
% case2_compare = 'maxmascons10_mincorrp3';
% case2_name = 'Model cov.';
case2_compare = 'adjcorrp05_depthadj2500_maxmascons10_mincorrp3';
case2_name = 'Model cov. + depth adj.';


corrcoeff_case1 = eval(['corrcoeff_',case1_compare]);
corrcoeff_case2 = eval(['corrcoeff_',case2_compare]);

good_ind = find(isnan(corrcoeff_case1 + corrcoeff_case2) == 0);


% % plot correlation improvement, color-coded by site depth
% 
% c_levels = (2000:200:6000)';
% [irregular_clevels_vec,~] = irregular_clevels_plot(obs_depth,c_levels);
% 
% cmap = 0.85*colormap(flip(bcyr(length(c_levels) - 1),1));
% 
% 
% fig1 = figure(1);
% colormap(cmap)
% circle_radius = 0.01;
% for curr_pt = 1:length(good_ind)
%     h = rectangle('Position',[([corrcoeff_case1(good_ind(curr_pt)) corrcoeff_case2(good_ind(curr_pt))] - (circle_radius*[1 1])) (circle_radius*[2 2])],'Curvature',[1 1]);
%     set(h,'FaceColor',cmap(ceil(irregular_clevels_vec(good_ind(curr_pt))),:),'EdgeColor','none')
% %     set(h,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
% end
% hold on
% line([-1 1],[-1 1],'Color',[0 0 0],'LineStyle','--')
% hold off
% set(gca,'xlim',[-0.5 0.8],'xtick',(-0.4):0.2:0.8,'xgrid','on','ylim',[-0.5 0.8],'ygrid','on','DataAspectRatio',[1 1 1],'FontSize',12)
% xlabel([case1_name,' correlation'])
% ylabel('Adjusted correlation')
% title({'Correlation improvement plot, color-coded by site depth'; ' '},'FontSize',8)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)])
% set(cbar,'ticks',0:5:20,'ticklabels',{'2000' '3000' '4000' '5000' '6000'},'FontSize',12)
% set(get(cbar,'label'),'String','Ocean bottom depth (meters)','FontSize',12)
% 
% print(fig1,['Corrcoeff_obs_downscale_obsdepth_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
% close(fig1)


% plot correlation improvement, color-coded by ECCO2 unfiltered (daily average) std. dev.

c_levels = [(200:20:300)'; (330:30:450)'; (500:50:700)'; (800:100:1000)'];
% c_levels = [(35:5:60)'; (70:10:100)'; (120:20:200)'; 250; 300; 400];
[irregular_clevels_vec,~] = irregular_clevels_plot(ECCO_stddev_unfiltered,c_levels);

cmap = 0.85*colormap(bcyr(length(c_levels) - 1));

fig2 = figure(2);
close(figure(1))
colormap(cmap)
circle_radius = 0.02;
for curr_pt = 1:length(good_ind)
    h = rectangle('Position',[([corrcoeff_case1(good_ind(curr_pt)) corrcoeff_case2(good_ind(curr_pt))] - (circle_radius*[1 1])) (circle_radius*[2 2])],'Curvature',[1 1]);
    set(h,'FaceColor',cmap(ceil(irregular_clevels_vec(good_ind(curr_pt))),:),'EdgeColor','none')
end
hold on
line([-1 1],[-1 1],'Color',[0 0 0],'LineStyle','--')
hold off
set(gca,'xlim',[-0.5 0.8],'xtick',(-0.4):0.2:0.8,'xgrid','on','ylim',[-0.5 0.8],'ygrid','on','DataAspectRatio',[1 1 1],'FontSize',12)
xlabel([case1_name,' correlation'])
ylabel('Downscaled+adjusted correlation')
title({'Correlation improvement plot, color-coded by ECCO2 unfiltered std. dev.'; ' '},'FontSize',8)
cbar = colorbar('location','southoutside');
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'ticks',[0 5 11 16],'ticklabels',{'0.02' '0.03' '0.05' '0.08'},'FontSize',12)
% set(cbar,'ticks',[1 5 9 14],'ticklabels',{'0.004' '0.006' '0.01' '0.02'},'FontSize',12)
set(get(cbar,'label'),'String','OBP daily std. dev. (dbar)','FontSize',12)

print(fig2,['Corrcoeff_obs_downscale_ECCOstdunfilt_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
close(fig2)


% compute time mean surface EKE at each obs. site

nc_file_name_satobs = '/indopac/adelman/Global_mesoscale/SSH/SSH_lp_elong_eddy_SSALTO_DUACS_5day_global.nc';
grid_res_satobs = 1/4;
lon_satobs = ncread(nc_file_name_satobs,'longitude');
lat_satobs = ncread(nc_file_name_satobs,'latitude');
time_datenum_satobs = ncread(nc_file_name_satobs,'time');

% obs_EKE_surface_tmean = NaN([length(obs_lon) 1]);
% for curr_pt = 1:length(obs_lon)
%     if ismember(curr_pt,good_ind) == 0
%         continue
%     end
%     
%     % load satellite SSH and compute EKE
%     lon_lat_radius_from_pt = 2;
%     in_lon_range_sat_ind = find(abs(mod(lon_satobs - obs_lon(curr_pt) + 180,360) - 180) <= lon_lat_radius_from_pt);
%     in_lat_range_sat_ind = find(abs(lat_satobs - obs_lat(curr_pt)) <= lon_lat_radius_from_pt);
%     lat_in_range_satobs = lat_satobs(in_lat_range_sat_ind);
%     if max(diff(in_lon_range_sat_ind)) > 5
%         gap_ind = find(diff(in_lon_range_sat_ind) > 5);
%         in_lon_range_sat_ind = [(in_lon_range_sat_ind(gap_ind + 1):1:max(in_lon_range_sat_ind))'; (min(in_lon_range_sat_ind):1:in_lon_range_sat_ind(gap_ind))'];
%         lon_in_range_satobs = lon_satobs(in_lon_range_sat_ind);
%         start_vec = [in_lon_range_sat_ind(1) in_lat_range_sat_ind(1) 1];
%         count_vec = [(max(in_lon_range_sat_ind) - in_lon_range_sat_ind(1) + 1) length(in_lat_range_sat_ind) length(time_datenum_satobs)];
%         SSH_satobs = ncread(nc_file_name_satobs,'SSH',start_vec,count_vec);
%         start_vec_2 = [min(in_lon_range_sat_ind) in_lat_range_sat_ind(1) 1];
%         count_vec_2 = [(in_lon_range_sat_ind(length(in_lon_range_sat_ind)) - min(in_lon_range_sat_ind) + 1) length(in_lat_range_sat_ind) length(time_datenum_satobs)];
%         SSH_satobs_2 = ncread(nc_file_name_satobs,'SSH',start_vec_2,count_vec_2);
%         SSH_satobs = [SSH_satobs; SSH_satobs_2];
%     else
%         lon_in_range_satobs = lon_satobs(in_lon_range_sat_ind);
%         start_vec = [in_lon_range_sat_ind(1) in_lat_range_sat_ind(1) 1];
%         count_vec = [length(in_lon_range_sat_ind) length(in_lat_range_sat_ind) length(time_datenum_satobs)];
%         SSH_satobs = ncread(nc_file_name_satobs,'SSH',start_vec,count_vec);
%     end
%     
%     lon_in_range_satobs = mod(lon_in_range_satobs - (lon_in_range_satobs(1) - 0.5),360) + (lon_in_range_satobs(1) - 0.5);
%     diff_lon_in_range_satobs = diff(lon_in_range_satobs);
%     
% 
%     size_array_satobs = size(SSH_satobs);
%     
%     nan_mask = ones(size(SSH_satobs));
%     nan_mask(isnan(SSH_satobs) == 1) = 0;
%     sum_nan_mask = sum(nan_mask,3);
%     not_enough_ind = sum_nan_mask < (0.9*max(max(sum_nan_mask)));
%     
%     
%     
%     % compute EKE associated with each SSH component
%     
%     nan_mask = ones(size(SSH_satobs));
%     nan_mask(isnan(SSH_satobs) == 1) = 0;
%     SSH_satobs_nans_zeroed = SSH_satobs;
%     SSH_satobs_nans_zeroed(isnan(SSH_satobs) == 1) = 0;
%     sum_nan_mask = sum(nan_mask,3);
%     nan_mask_not_enough = ones(size(sum_nan_mask));
%     nan_mask_not_enough(sum_nan_mask < (0.9*max(max(sum_nan_mask)))) = 0;
%     nan_mask_not_enough_3D = repmat(nan_mask_not_enough,[1 1 size(nan_mask,3)]);
%     nan_mask_SSH = nan_mask & nan_mask_not_enough_3D;
%     
%     clear nan_mask nan_mask_not_enough*
%     
%     SSH_mean_satobs = sum(nan_mask_SSH.*SSH_satobs_nans_zeroed,3)./sum(nan_mask_SSH,3);
%     
%     SSH_nomean_satobs = SSH_satobs - repmat(SSH_mean_satobs,[1 1 size(SSH_satobs,3)]);
%     
%     clear SSH_*_nans_zeroed
%     
%     
%     % compute EKE from geostrophic gradients
%     
%     lon_satobs_midpts = lon_in_range_satobs(2:length(lon_in_range_satobs)) - (diff(lon_in_range_satobs)/2);
%     lat_satobs_midpts = lat_in_range_satobs(2:length(lat_in_range_satobs)) - (diff(lat_in_range_satobs)/2);
%     dx = 111100*(repmat(cosd(lat_satobs_midpts'),[(size_array_satobs(1) - 1) 1]).*repmat(diff(lon_in_range_satobs),[1 (size_array_satobs(2) - 1)]));
%     dy = 111100*repmat(diff(lat_in_range_satobs'),[(size_array_satobs(1) - 1) 1]);
%     
%     diff_x_SSH_nomean_in_r = diff(SSH_nomean_satobs,1,1);
%     diff_x_SSH_nomean_satobs = diff_x_SSH_nomean_in_r(:,2:size(diff_x_SSH_nomean_in_r,2),:) - (diff(diff_x_SSH_nomean_in_r,1,2)/2);
%     diff_y_SSH_nomean_in_r = diff(SSH_nomean_satobs,1,2);
%     diff_y_SSH_nomean_satobs = diff_y_SSH_nomean_in_r(2:size(diff_y_SSH_nomean_in_r,1),:,:) - (diff(diff_y_SSH_nomean_in_r,1,1)/2);
%     
%     
%     f = repmat(2*((2*pi)/86164)*sind(lat_satobs_midpts'),[length(lon_satobs_midpts) 1]);
%         
%     uvel_geostr_nomean_satobs = -9.81*repmat((f.^(-1)),[1 1 size(SSH_nomean_satobs,3)]).*(diff_y_SSH_nomean_satobs./repmat(dy,[1 1 size(SSH_nomean_satobs,3)]));
%     vvel_geostr_nomean_satobs = 9.81*repmat((f.^(-1)),[1 1 size(SSH_nomean_satobs,3)]).*(diff_x_SSH_nomean_satobs./repmat(dx,[1 1 size(SSH_nomean_satobs,3)]));
%     
%     clear diff_*
%     
%     
%     EKE_satobs = (0.5*((uvel_geostr_nomean_satobs.^2) + (vvel_geostr_nomean_satobs.^2)));
%     
%     clear uvel_geostr* vvel_geostr*
%     
%     % mask out low-latitude points
%     lat_mask = ones(size(lat_satobs_midpts));
%     lat_mask(abs(lat_satobs_midpts) < 5) = 0;
%     lat_mask_3d = repmat(lat_mask',[size(EKE_satobs,1) 1 size(EKE_satobs,3)]);
%     lat_mask_ind = find(abs(lat_mask_3d) < 1e-5);
%     
%     EKE_satobs(lat_mask_ind) = NaN;
%     
%     nan_mask_EKE = ones(size(EKE_satobs));
%     nan_mask_EKE(isnan(EKE_satobs) == 1) = 0;
%     EKE_satobs_nans_zeroed = EKE_satobs;
%     EKE_satobs_nans_zeroed(isnan(EKE_satobs) == 1) = 0;
%     
%     EKE_mean_satobs = sum(nan_mask_EKE.*EKE_satobs_nans_zeroed,3)./sum(nan_mask_EKE,3);
%     EKE_mean_satobs(isnan(EKE_mean_satobs) == 1) = 0;
%     
%     obs_EKE_surface_tmean(curr_pt) = interp2_fft(mod(lon_satobs_midpts - (obs_lon(curr_pt) - 180),360) + (obs_lon(curr_pt) - 180),lat_satobs_midpts,EKE_mean_satobs,obs_lon(curr_pt),obs_lat(curr_pt));
%     
%     if abs(obs_EKE_surface_tmean(curr_pt)) < 1e-5
%         obs_EKE_surface_tmean(curr_pt) = NaN;
%     end
%     
% end



% plot correlation improvement, color-coded by surface EKE time mean
c_levels = (1e-4)*[(0:10:50)'; 65; 80; 100; (130:30:250)'; (300:50:500)'; 560; 630; 700; 800; 2000];
[irregular_clevels_vec,~] = irregular_clevels_plot(obs_EKE_surface_tmean,c_levels);

cmap = 0.85*colormap(bcyr(length(c_levels) - 1));

fig3 = figure(3);
close(figure(1))
colormap(cmap)
circle_radius = 0.02;
for curr_pt = 1:length(good_ind)
    if isnan(irregular_clevels_vec(good_ind(curr_pt))) == 1
        continue
    end
    h = rectangle('Position',[([corrcoeff_case1(good_ind(curr_pt)) corrcoeff_case2(good_ind(curr_pt))] - (circle_radius*[1 1])) (circle_radius*[2 2])],'Curvature',[1 1]);
    set(h,'FaceColor',cmap(ceil(irregular_clevels_vec(good_ind(curr_pt))),:),'EdgeColor','none')
end
hold on
line([-1 1],[-1 1],'Color',[0 0 0],'LineStyle','--')
hold off
set(gca,'xlim',[-0.5 0.8],'xtick',(-0.4):0.2:0.8,'xgrid','on','ylim',[-0.5 0.8],'ygrid','on','DataAspectRatio',[1 1 1],'FontSize',12)
xlabel([case1_name,' correlation'])
ylabel('Downscaled+adjusted correlation')
title({'Correlation improvement plot, color-coded by time mean surface EKE'; ' '},'FontSize',8)
cbar = colorbar('location','southoutside');
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'ticks',[0 5 8 13 18 22],'ticklabels',{'0' '50' '100' '250' '500' '800'},'FontSize',12)
set(get(cbar,'label'),'String','Surface EKE time mean (cm^2 s^{-2})','FontSize',12)

print(fig3,['Corrcoeff_obs_downscale_surfaceEKE_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
close(fig3)


% plot standard deviation (normalized by obs. std. dev.)

obs_stddev = obs_stddev_tinterp;
if strcmp(case1_compare,'coloc_mascons') == 1
    stddev_case1 = obs_stddev_coloc;
else
    stddev_case1 = eval(['obs_stddev_',case1_compare]);
end
stddev_case2 = eval(['obs_stddev_',case2_compare]);

good_ind = find(isnan(obs_stddev + stddev_case1 + stddev_case2) == 0);

stddev_ratio_case1 = stddev_case1./obs_stddev;
stddev_ratio_case2 = stddev_case2./obs_stddev;


% plot normalized standard deviations, color-coded

% c_levels = (2000:200:6000)';
% [irregular_clevels_vec,~] = irregular_clevels_plot(obs_depth,c_levels);
% c_levels = (1e-4)*[(0:10:50)'; 65; 80; 100; (130:30:250)'; (300:50:500)'; 560; 630; 700; 800; 2000];
% [irregular_clevels_vec,~] = irregular_clevels_plot(obs_EKE_surface_tmean,c_levels);
c_levels = [(200:20:300)'; (330:30:450)'; (500:50:700)'; (800:100:1000)'];
[irregular_clevels_vec,~] = irregular_clevels_plot(ECCO_stddev_unfiltered,c_levels);

% cmap = 0.85*colormap(flip(bcyr(length(c_levels) - 1),1));
cmap = 0.85*colormap(bcyr(length(c_levels) - 1));

fig1 = figure(1);
colormap(cmap)
circle_radius = 0.035;
for curr_pt = 1:length(good_ind)
    if isnan(irregular_clevels_vec(good_ind(curr_pt))) == 1
        continue
    end
    h = rectangle('Position',[([stddev_ratio_case1(good_ind(curr_pt)) stddev_ratio_case2(good_ind(curr_pt))] - (circle_radius*[1 1])) (circle_radius*[2 2])],'Curvature',[1 1]);
    set(h,'FaceColor',cmap(ceil(irregular_clevels_vec(good_ind(curr_pt))),:),'EdgeColor','none')
%     set(h,'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
end
hold on
line([0 2],[0 2],'Color',[0 0 0],'LineStyle','--')
line([1 1],[0 2],'Color',[0 0 0],'LineStyle','-.')
line([0 2],[1 1],'Color',[0 0 0],'LineStyle','-.')
hold off
set(gca,'xlim',[0 2],'xtick',0:0.25:2,'xticklabel',{'0' '' '0.5' '' '1' '' '1.5' '' '2'},'xgrid','on','ylim',[0 2],'ytick',0:0.25:2,'yticklabel',{'0' '' '0.5' '' '1' '' '1.5' '' '2'},'ygrid','on','DataAspectRatio',[1 1 1],'FontSize',12)
xlabel([case1_name,' normalized std. dev.'])
ylabel('Downscaled+adj. norm. std. dev.')
% title({'Normalized standard deviation plot, color-coded by site depth'; ' '},'FontSize',8)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)])
% set(cbar,'ticks',0:5:20,'ticklabels',{'2000' '3000' '4000' '5000' '6000'},'FontSize',12)
% set(get(cbar,'label'),'String','Ocean bottom depth (meters)','FontSize',12)
% title({'Normalized standard deviation plot, color-coded by time mean surface EKE'; ' '},'FontSize',8)
% cbar = colorbar('location','southoutside');
% set(gca,'clim',[0 size(cmap,1)])
% set(cbar,'ticks',[0 5 8 13 18 22],'ticklabels',{'0' '50' '100' '250' '500' '800'},'FontSize',12)
% set(get(cbar,'label'),'String','Surface EKE time mean (cm^2 s^{-2})','FontSize',12)
title({'Normalized standard deviation plot, color-coded by ECCO2 unfiltered std. dev.'; ' '},'FontSize',8)
cbar = colorbar('location','southoutside');
set(gca,'clim',[0 size(cmap,1)])
set(cbar,'ticks',[0 5 11 16],'ticklabels',{'0.02' '0.03' '0.05' '0.08'},'FontSize',12)
set(get(cbar,'label'),'String','OBP daily std. dev. (dbar)','FontSize',12)

% print(fig1,['Stddev_norm_obs_downscale_obsdepth_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
% print(fig1,['Stddev_norm_obs_downscale_surfaceEKE_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
print(fig1,['Stddev_norm_obs_downscale_ECCOstdunfilt_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.bmp'],'-dbmp','-r300')
close(fig1)


% plot Taylor diagrams

RMS_diff_case1 = (1 + (stddev_ratio_case1.^2) - (2*corrcoeff_case1.*stddev_ratio_case1)).^(1/2);
RMS_diff_case2 = (1 + (stddev_ratio_case2.^2) - (2*corrcoeff_case2.*stddev_ratio_case2)).^(1/2);

good_poscorr_ind = good_ind(corrcoeff_case1(good_ind) > 0);

fig1 = figure(1);
[hp axl] = taylordiag([1; stddev_ratio_case1(good_poscorr_ind)],[0; RMS_diff_case1(good_poscorr_ind)],[1; corrcoeff_case1(good_poscorr_ind)],'limstd',2,'colRMS',[0.4 0.4 0.4],'tickRMSangle',115,'markersize',8,'markerstyle','o','markercolor',[0.8 0 0]);
set(hp(1),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor','none')
set(axl(1).handle,'String','Normalized std. dev.')
title([case1_name,' Taylor diagram, normalized by in-situ obs. std. dev.'],'FontSize',8)

print(fig1,['Taylor_diag_obs_',case1_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.pdf'],'-dpdf','-r300')
close(fig1)


good_poscorr_ind = good_ind(corrcoeff_case2(good_ind) > 0);

fig1 = figure(1);
[hp axl] = taylordiag([1; stddev_ratio_case2(good_poscorr_ind)],[0; RMS_diff_case2(good_poscorr_ind)],[1; corrcoeff_case2(good_poscorr_ind)],'limstd',2,'colRMS',[0.4 0.4 0.4],'tickRMSangle',115,'markersize',8,'markerstyle','o','markercolor',[0.8 0 0]);
set(hp(1),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerFaceColor','none')
set(axl(1).handle,'String','Normalized std. dev.')
title([case2_name,' Taylor diagram, normalized by in-situ obs. std. dev.'],'FontSize',8)

print(fig1,['Taylor_diag_obs_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.pdf'],'-dpdf','-r300')
close(fig1)


% plot Taylor diagram with arrows showing change from one case to another

good_poscorr_ind = good_ind((corrcoeff_case1(good_ind) > 0) & (corrcoeff_case2(good_ind) > 0));

fig1 = figure(1);
[hp axl] = taylordiag_change([[1 1]; [stddev_ratio_case1(good_poscorr_ind) stddev_ratio_case2(good_poscorr_ind)]],[[0 0]; [RMS_diff_case1(good_poscorr_ind) RMS_diff_case2(good_poscorr_ind)]],[[1 1]; [corrcoeff_case1(good_poscorr_ind) corrcoeff_case2(good_poscorr_ind)]],'limstd',2,'colRMS',[0.4 0.4 0.4],'tickRMSangle',115,'linewidth',0.5,'linecolor',[0.8 0 0]);
set(axl(1).handle,'String','Normalized std. dev.')
hold on
rectangle('Position',[0.96 -0.04 0.08 0.08],'Curvature',[1 1],'EdgeColor',[0.4 0.4 0.4],'FaceColor','none')
hold off
title([case1_name,'-',case2_name,' change Taylor diagram, normalized by in-situ obs. std. dev.'],'FontSize',8)

print(fig1,['Taylor_diag_obs_change_',case1_compare,'_',case2_compare,'_',num2str(lowperiod),'_',num2str(highperiod),'_periods.pdf'],'-dpdf','-r300')
close(fig1)