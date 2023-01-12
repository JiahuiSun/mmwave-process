function [x_vector, y_vector, z_vector] = advance_music(Xcube,detout,num_tx,num_rx,fft_Ang,Is_Windowed)
Nr=size(Xcube,1);   %%%length of Chirp: 512
Ne=size(Xcube,2);   %%%length of channel: 12=4*3
Nd=size(Xcube,3);   %%%length of chirp loop: 120

num_detected_obj=size(detout,2);

virtual_ant=zeros(Ne,num_detected_obj);
for i=1:num_detected_obj
    virtual_ant(:,i)=Xcube(detout(2,i),:,detout(1,i));
end

azimuth_ant = virtual_ant(1:2 * num_rx, :);
azimuth_music=[];
for i=1:num_detected_obj
    [values,~]=pmusic(azimuth_ant(:,i),1);
    azimuth_music(:,i) = values;
end

[~,k_max] = max(abs(azimuth_music));  % shape = (1,num_detected_obj)
% disp(["Here" length(azimuth_music(1,:)) length(azimuth_music(:,1)) length(k_max)]);
for i=1:length(k_max)
    k_max(i)=k_max(i)/256*2;
    if k_max(i)>1
        k_max(i)=k_max(i)-2;
    end
end
x_vector = k_max;

elevation_music=zeros([256 num_detected_obj]);
elevation_ant=zeros([5 num_detected_obj]);
elevation_ant(1,:)=virtual_ant(12,:);
elevation_ant(2,:)=virtual_ant(6,:);

for j=3:5
    elevation_ant(j,:) = elevation_ant(j-1,:)-virtual_ant(14-j, :)+virtual_ant(8-j, :);
end

for i=1:num_detected_obj
        [values,~]=pmusic(elevation_ant(:,i),1);
        elevation_music(:,i) = values;
end

[~,k2_max] = max(abs(elevation_music));  % shape = (1,num_detected_obj)
   

for i=1:length(k2_max)
    k2_max(i)=k2_max(i)/256*2;
    if k2_max(i)>1
        k2_max(i)=k2_max(i)-2;
    end
end

z_vector=k2_max;

% Calculate elevation phase shift

y_vector = sqrt(1 - x_vector .^ 2 - z_vector .^ 2);
% disp(["music" x_vector]);
end