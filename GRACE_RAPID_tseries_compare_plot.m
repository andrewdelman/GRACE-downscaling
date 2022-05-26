% plot non-downscaled and downscaled GRACE time series of volume transport, along with RAPID

path(path,'~/GRACE/')
path(path,'~/plotting_scripts/')
cd('/indopac/adelman/GRACE/')


load('GRACE_RAPID_tseries_method_comparisons.mat')

time_range_start = [2004 4 1];
time_range_end = [2017 2 1];

freq_range_id = '426_8351';
downscale_id = 'cs510_cylindweight1.5_depthadj2500_hybrid0_maxmascons10_mincorr0.3';
GRACE_flux_total_nondownsc = (1e-6)*vol_flux_total_coloc_interannual_only';
GRACE_flux_total_downsc = (1e-6)*vol_flux_total_depthadj_interannual_only';
RAPID_tseries_interp = RAPID_tseries_interp_interannual_only;


time_mid_interp = NaN([((2*length(time)) + 1) 1]);
time_mid_interp(1) = ((1.5)*time(1)) + ((-0.5)*time(2));
time_mid_interp(2:2:(2*length(time))) = time;
time_mid_interp(3:2:((2*length(time)) - 1)) = time(2:length(time)) - (diff(time)/2);
time_mid_interp((2*length(time)) + 1) = ((-0.5)*time(length(time) - 1)) + ((1.5)*time(length(time)));
GRACE_flux_total_nondownsc_only = NaN([((2*length(time)) + 1) 1]);
GRACE_flux_total_nondownsc_only(1) = ((1.5)*GRACE_flux_total_nondownsc(1)) + ((-0.5)*GRACE_flux_total_nondownsc(2));
GRACE_flux_total_nondownsc_only(2:2:(2*length(time))) = GRACE_flux_total_nondownsc;
GRACE_flux_total_nondownsc_only(3:2:((2*length(time)) - 1)) = GRACE_flux_total_nondownsc(2:length(time)) - (diff(GRACE_flux_total_nondownsc)/2);
GRACE_flux_total_nondownsc_only((2*length(time)) + 1) = ((-0.5)*GRACE_flux_total_nondownsc(length(time) - 1)) + ((1.5)*GRACE_flux_total_nondownsc(length(time)));
improv_ind = find(abs(GRACE_flux_total_downsc - RAPID_tseries_interp) - abs(GRACE_flux_total_nondownsc - RAPID_tseries_interp) < 0);
GRACE_flux_total_nondownsc_only(2*improv_ind) = NaN;
GRACE_flux_total_downsc_improv_only = NaN([((2*length(time)) + 1) 1]);
GRACE_flux_total_downsc_improv_only(1) = ((1.5)*GRACE_flux_total_downsc(1)) + ((-0.5)*GRACE_flux_total_downsc(2));
GRACE_flux_total_downsc_improv_only(2:2:(2*length(time))) = GRACE_flux_total_downsc;
GRACE_flux_total_downsc_improv_only(3:2:((2*length(time)) - 1)) = GRACE_flux_total_downsc(2:length(time)) - (diff(GRACE_flux_total_downsc)/2);
GRACE_flux_total_downsc_improv_only((2*length(time)) + 1) = ((-0.5)*GRACE_flux_total_downsc(length(time) - 1)) + ((1.5)*GRACE_flux_total_downsc(length(time)));
not_improv_ind = find(abs(GRACE_flux_total_downsc - RAPID_tseries_interp) - abs(GRACE_flux_total_nondownsc - RAPID_tseries_interp) > 0);
GRACE_flux_total_downsc_improv_only(2*not_improv_ind) = NaN;


fig5 = figure(5);
h = plot(time,GRACE_flux_total_nondownsc,time,GRACE_flux_total_downsc,time,RAPID_tseries_interp);
tick_year_spacing = 2;
years_to_plot_ticks = ((ceil(time_range_start(1)/tick_year_spacing)*tick_year_spacing):tick_year_spacing:(ceil((time_range_end(1) - 1)/tick_year_spacing)*tick_year_spacing))';
xtick_datenums_plot = datenum([years_to_plot_ticks ones(length(years_to_plot_ticks),2)]);
xtick_labels_plot = cell(length(years_to_plot_ticks),1);
for xtick_ind = 1:length(years_to_plot_ticks)
    xtick_labels_plot{xtick_ind} = num2str(years_to_plot_ticks(xtick_ind));
end
set(gca,'xlim',[datenum([2004 4 1]) datenum([2017 2 1])],'xtick',xtick_datenums_plot,'xticklabel',xtick_labels_plot,'xgrid','on')
set(h(1),'Color',[0 0 0.3],'LineWidth',2)
set(h(2),'Color',[0 0.5 1],'LineWidth',2)
set(h(3),'Color',[0.8 0 0],'LineWidth',2)
hold on
line([datenum(time_range_start) datenum(time_range_end)],[0 0],[0 0],'Color',[0 0 0],'LineWidth',1)
% highlight_line_1 = line(time_mid_interp,GRACE_flux_total_nondownsc_only,ones([length(time_mid_interp) 1]),'Color',[0 0 0.3],'LineWidth',2);
% highlight_line_2 = line(time_mid_interp,GRACE_flux_total_downsc_improv_only,ones([length(time_mid_interp) 1]),'Color',[0 0.5 1],'LineWidth',2);
lower_bars_vec_1 = -3.5*ones(size(GRACE_flux_total_nondownsc_only));
lower_bars_vec_1(isnan(GRACE_flux_total_nondownsc_only) == 1) = NaN;
highlight_lower_1 = line(time_mid_interp,lower_bars_vec_1,ones(length(time_mid_interp)),'Color',[0 0 0.3],'LineWidth',3);
lower_bars_vec_2 = -3.5*ones(size(GRACE_flux_total_downsc_improv_only));
lower_bars_vec_2(isnan(GRACE_flux_total_downsc_improv_only) == 1) = NaN;
highlight_lower_2 = line(time_mid_interp,lower_bars_vec_2,ones(length(time_mid_interp)),'Color',[0 0.5 1],'LineWidth',3);
hold off
ylabel('Volume flux anomaly (Sv)')
title({'Volume flux anomaly in the Atlantic at 26.5 ^{o} lat, 3000-5000 m depth.'; 'Comparison of non-downscaled and downscaled GRACE with RAPID time series.'; ' '},'FontSize',10)
% leg = legend([highlight_line_1; highlight_line_2; h(3)],'Non-downscaled GRACE','Downscaled GRACE','RAPID','location','northeast');
leg = legend(h,'Non-downscaled GRACE','Downscaled GRACE','RAPID','location','northeast');
keyboard
% saveas(fig5,['Vol_flux_anom_time_series_GRACE_RAPID_compare_26.5_lat_transect_Atlantic_3000_5000_depth_',freq_range_id,'_periods_',downscale_id,'.pdf'])
print(fig5,['Vol_flux_anom_time_series_GRACE_RAPID_compare_26.5_lat_transect_Atlantic_3000_5000_depth_',freq_range_id,'_periods_',downscale_id,'.bmp'],'-dbmp','-r300')
close(fig5)