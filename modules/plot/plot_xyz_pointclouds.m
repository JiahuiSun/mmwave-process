% plot 3D point clouds
function [axh] = plot_xyz_pointclouds(detout, cnt, ds)
% detout format: % [range bin, velocity bin, angle bin, power, range(m), ...
% velocity (m/s), angle(degree)]

figure('visible','off')
% x-direction: Doppler, y-direction: angle, z-direction: 
% frame num,power,v,x,y,z

if ds==1  % mHomeGes
    power = detout(:,7);
    x_value=detout(:,3);
    y_value=detout(:,4);
    z_value=detout(:,5);
elseif ds==0  % processed point cloud
    power = detout(:,4);
    x_value=detout(:,8);
    y_value=detout(:,9);
    z_value=detout(:,10);
elseif ds==2  % our saved point cloud
    power = detout(:,2);
    x_value=detout(:,4);
    y_value=detout(:,5);
    z_value=detout(:,6);
elseif ds==3  % pantomime
    power = detout(:,5);
    x_value=detout(:,2);
    y_value=detout(:,3);
    z_value=detout(:,4);
else
    power = abs(detout(:,4));
    x_value=detout(:,1);
    y_value=detout(:,2);
    z_value=detout(:,3);
end

[axh] = scatter3(x_value,y_value,z_value,20,power, 'filled');
set(gca,'position',[0.08,0.11,0.77,0.8]);
% ax = gca;
% ax.XDir = 'reverse';
view(0,0);
% view(0,0); 
cb = colorbar('position',[0.87 0.11 0.04 0.8]);
% colorbar('position',[0.15 0.15 0.04 0.2])
colormap turbo;
cb.Label.String = 'veloctity';
caxis([0,1.5]);

%============================
% hold on
% chest=[0,0.65,0];
% head=[-0.35,0.65,0];
% bottom=[0.43,0.65,0];
% left_shoulder=[0,0.65,0.2];
% right_shoulder=[0,0.65,-0.2];
% 
% plot3([chest(1) head(1)],[chest(2) head(2)],[chest(3) head(3)],'r--','LineWidth',2);
% hold on
% plot3([chest(1) bottom(1)],[chest(2) bottom(2)],[chest(3) bottom(3)],'r--','LineWidth',2);
% hold on
% plot3([chest(1) left_shoulder(1)],[chest(2) left_shoulder(2)],[chest(3) left_shoulder(3)],'r--','LineWidth',2);
% hold on
% plot3([chest(1) right_shoulder(1)],[chest(2) right_shoulder(2)],[chest(3) right_shoulder(3)],'r--','LineWidth',2);
% hold on
%============================

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-1 1 2 4 -1 1]);
% title='xyz point cloud '+num2str(cnt);
title(['xyz point cloud no.', num2str(cnt)]);
grid on

end