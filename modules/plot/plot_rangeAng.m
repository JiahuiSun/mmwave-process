% plot range-angle heatmap
function [axh] = plot_rangeAng(Xcube,rng_grid,agl_grid, cnt)
Nr = size(Xcube,1);   %%%length of Chirp(num of rangeffts)
Ne = size(Xcube,2);   %%%number of angleffts
Nd = size(Xcube,3);   %%%length of chirp loop

% polar coordinates
for i = 1:size(agl_grid,2)
    yvalue(i,:) = (sin(agl_grid(i)*pi/180 )).*rng_grid;
end
for i=1:size(agl_grid,2)
    xvalue(i,:) = (cos(agl_grid(i)*pi/180)).*rng_grid;
end

%% plot 2D(range-angle)

Xpow = abs(Xcube);
Xpow = squeeze(sum(Xpow,3)/size(Xpow,3));
% disp(size(Xpow));
% noisefloor = db2pow(-15);
Xsnr=Xpow;
% Xsnr = pow2db(Xpow/noisefloor);
Xsnr(1:20,:)=0;
figure('visible','off')
% figure()
% set(gcf,'Position',[10,10,530,420])
[axh] = surf(agl_grid,rng_grid,Xsnr);
ax = gca;
% ax.XDir = 'reverse';
view(0,90)
% axis([-90 90 0 25]);
% axis([-60 60 2 25]);
colorbar
% caxis([0 50])
axis([-90 90 0 7]);
grid off
shading interp

xlabel('Angle of arrive(degrees)')
ylabel('Range(meters)')

% caxis([0,0.2])
title(['Range-Angle heatmap no.', num2str(cnt)])




end