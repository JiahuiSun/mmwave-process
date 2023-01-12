% plot range-angle heatmap
function [axh] = plot_rangeDop(Dopdata_sum,rng_grid,vel_grid, cnt, Resl_indx)
% disp(size(Dopdata_sum));
Dopdata_sum(1:10,:)=0;
% maximum = max(max(Dopdata_sum));
% [x,y]=find(Dopdata_sum==maximum);
% disp(maximum);
% disp([x,y]);
% plot 2D(range-Doppler)
figure('visible','off')
% set(gcf,'Position',[10,10,530,420])

surf(vel_grid,rng_grid(20:100,:),Dopdata_sum(20:100,:));
hold on;
scatter3(vel_grid(Resl_indx(1,:)), rng_grid(Resl_indx(2,:)), Resl_indx(3,:).^(0.5), 50, 'red', 'o');
hold off;

ax = gca;
ax.XDir = 'reverse';
view(0,90)
axis([-2 2 1 4]);
shading interp;
xlabel('Doppler Velocity (m/s)');
ylabel('Range(meters)');
colorbar;
% clim([0 max(Dopdata_sum)]);
% caxis([0,1e6]);
grid off;
% axis off;
% box off;
% colorbar off;
% set(gca,'looseInset',[0 0 0 0]);
title(['Range-Doppler heatmap no.', num2str(cnt), ', obj: ', num2str(length(Resl_indx))]);
end