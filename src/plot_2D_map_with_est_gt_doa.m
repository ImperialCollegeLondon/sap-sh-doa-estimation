function[handles, cax] = plot_2D_map_with_est_gt_doa(az_grid,inc_grid,counts,est_doa,gt_doa,contour_levels)

if nargin < 6 || isempty(contour_levels)
    imagesc(180/pi*az_grid(1,:),180/pi*inc_grid(:,1),counts)
    set(gca,'ydir','normal')
else
    contourf(180/pi*az_grid,180/pi*inc_grid,counts,contour_levels,'showtext','on')
end
hold all
legstr = {};
handles = [];
if nargin > 3 && ~isempty(est_doa)
    [est_az,est_inc,~] = mycart2sph(est_doa);
    handles = [handles plot_az_inc(est_az,est_inc,'+k')];
    plot_az_inc(est_az,est_inc,'dw');
    legstr = {legstr{:};'Estimated'};
end
if nargin > 4 && ~isempty(gt_doa)
    [gt_az,gt_inc,~] = mycart2sph(gt_doa);
    handles = [handles plot_az_inc(gt_az,gt_inc,'xk')];
    plot_az_inc(gt_az,gt_inc,'sw');
    legstr = {legstr{:};'Actual'};
end

if ~isempty(legstr)
    legend(handles,legstr,'visible','off');
end
cax = colorbar;
colormap(v_colormap('v_thermliny'));
xlabel('Azimuth [deg]');
ylabel('Inclination [deg]');