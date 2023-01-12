function [point_all,sep_list] = generate_point_cloud(bin_filename,Pfa,data_path)
addpath(genpath('.\config'));
addpath(genpath('.\utils'));
addpath(genpath('.\modules\plot'));
addpath(genpath('.\modules\fft'));
addpath(genpath('.\modules\detection'));
addpath(genpath('.\modules\cluster'));

%% parameter setting
params = get_params_value();
% constant parameters
c = params.c; % Speed of light in air (m/s)
fc = params.fc; % Center frequency (Hz)
lambda = params.lambda;
Rx = params.Rx;
Tx = params.Tx;

% configuration parameters
Fs = params.Fs;
sweepSlope = params.sweepSlope;
samples = params.samples;
loop = params.loop;

Tc = params.Tc; % us 
fft_Rang = params.fft_Rang;
fft_Vel = params.fft_Vel;
fft_Ang = params.fft_Ang;
num_crop = params.num_crop;
max_value = params.max_value; % normalization the maximum of data WITH 1843

% Creat grid table
rng_grid = params.rng_grid;
agl_grid = params.agl_grid;
vel_grid = params.vel_grid;

% Algorithm parameters
data_each_frame = samples*loop*Tx;

 % 1==> Windowing before doing range and angle fft
Is_Windowed = 1;

%% 循环读取数据
filename=data_path + bin_filename + '.bin';
% Frame_num = data_length/data_each_frame;
fid = fopen(filename);
sep_list=[];
point_all=[];
frame_cnt = 1;
while 1
    data_frame = get_next_frame(fid);
    data_chirp = [];
    if length(data_frame) < data_each_frame
        break
    end
    for cj = 1:Tx*loop
        temp_data = data_frame(:, (cj-1)*samples+1:cj*samples);
        data_chirp(:,:,cj) = temp_data;
    end
    
    % separate the chirps for TDM-MIMO with 3 TXs
    chirp_tx1 = data_chirp(:,:,1:3:end);  % RX*samples*loops
    chirp_tx3 = data_chirp(:,:,2:3:end);
    chirp_tx2 = data_chirp(:,:,3:3:end);

    % permutation with the format [samples, Rx, chirp]
    chirp_tx1 = permute(chirp_tx1, [2,1,3]);
    chirp_tx2 = permute(chirp_tx2, [2,1,3]);
    chirp_tx3 = permute(chirp_tx3, [2,1,3]);

    %% Range FFT for chirps from 3txs
    [Rangedata_tx1] = fft_range(chirp_tx1, fft_Rang, Is_Windowed);
    [Rangedata_tx2] = fft_range(chirp_tx2, fft_Rang, Is_Windowed);
    [Rangedata_tx3] = fft_range(chirp_tx3, fft_Rang, Is_Windowed);

    %% Doppler FFT
    Dopplerdata_tx1 = fft_doppler(Rangedata_tx1, fft_Vel, 0);
    Dopplerdata_tx2 = fft_doppler(Rangedata_tx2, fft_Vel, 0);
    Dopplerdata_tx3 = fft_doppler(Rangedata_tx3, fft_Vel, 0);
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));
    
    %% CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), cell power]
%     Pfa = 5e-3; % probability of false alarm
    [Resl_indx] = cfar_RV(Dopdata_sum, fft_Rang, num_crop, Pfa);
    detout = peakGrouping(Resl_indx);
    
    % doppler compensation on Rangedata_even using the max-intensity peak
    % on each range bin
    for ri = num_crop+1:fft_Rang-num_crop
        find_idx = find(detout(2, :) == ri);
        if isempty(find_idx)
            continue
        else
            % pick the first larger velocity
            pick_idx = find_idx(1);
            % phase compensation for virtual elements
            pha_comp_term_tx3 = exp(-1i * 2*pi * (detout(1,pick_idx)-fft_Vel/2-1) / (fft_Vel*3));
            Dopplerdata_tx3(ri, :, :) = Dopplerdata_tx3(ri, :, :) * pha_comp_term_tx3;
            pha_comp_term_tx2 = exp(-1i * 4*pi * (detout(1,pick_idx)-fft_Vel/2-1) / (fft_Vel*3));
            Dopplerdata_tx2(ri, :, :) = Dopplerdata_tx2(ri, :, :) * pha_comp_term_tx2;
        end
    end

    Dopdata_merge=[Dopplerdata_tx1, Dopplerdata_tx3, Dopplerdata_tx2];
    
    %% point cloud generation
%     [x_vector, y_vector, z_vector]=naive_xyz(Dopdata_merge,detout,Tx,Rx,fft_Ang,Is_Windowed); 
    [x_vector, y_vector, z_vector] = advance_music(Dopdata_merge, detout, Tx, Rx, fft_Ang, Is_Windowed); 
    
%     % Angle estimation for detected point clouds
%     Dopplerdata_merge = permute([Dopplerdata_tx1, Dopplerdata_tx3], [2, 1, 3]);
%     [Resel_agl, ~, rng_excd_list] = angle_estim_dets(detout, Dopplerdata_merge, fft_Vel, ...
%         fft_Ang, Rx, 2, num_crop);
    
%     % Transform bin index to range/velocity/angle
%     Resel_agl_deg = agl_grid(1, Resel_agl)';
    Resel_vel = vel_grid(detout(1,:), 1);
    Resel_rng = rng_grid(detout(2,:), 1);
    
    x_vec=x_vector'.*Resel_rng;
    y_vec=y_vector'.*Resel_rng;
    z_vec=z_vector'.*Resel_rng;
    
    % save_det data format below
%     % [1 range bin, 2 velocity bin, 3 power, 4 range(m), 5 velocity (m/s), 6 x ,7 y ,8 z]
%     save_det_data = [detout(2,:)', detout(1,:)', detout(3,:)', ...
%         Resel_rng, Resel_vel, x_vec,y_vec,z_vec];
%     disp(save_det_data);

%     % delete complex y value
%     valid_idx=find(imag(y_vec)==0);
%     save_det_data=save_det_data(valid_idx,:);
    save_det_data = real([x_vec,y_vec,z_vec,Resel_vel, detout(3,:)']);

%     if ~isempty(rng_excd_list)
%         save_det_data(rng_excd_list, :) = [];
%     end


%     % 按坐标、速度等过滤点
%     save_det_data_crop=save_det_data(find(save_det_data(:,1)>-1.5 & save_det_data(:,1)<1.5 &...  % x
%         save_det_data(:,2)>1.5 & save_det_data(:,2)<4.5 &...  % y
%         save_det_data(:,3)>-1.5 & save_det_data(:,3)<1.5),:);  % z
%         (save_det_data(:,5)<-0.02 | save_det_data(:,5)>0.02)), :);
%     save_det_data_crop=save_det_data;

%     % 去除噪点
%     idx = myDBScan(save_det_data_crop(:,1:3),0.5,2,1);
%     idx_set=unique(idx);
%     disp(length(idx_set)-1)
%     pc_body=save_det_data_crop(find(idx~=0),:);
    
    pc_body = save_det_data;

    point_cloud_feat=pc_body(:,[5,4,1,2,3]);  % power,velocity,x,y,z
    
    % 在第一列加上frame序号
    frame_vec=ones(size(point_cloud_feat,1),1)*frame_cnt;
    point_cloud_feat=[frame_vec point_cloud_feat];% frame num,power,v,x,y,z
    point_all=[point_all; point_cloud_feat];
    
    frame_cnt=frame_cnt+1;
end
% point_all(:,2)=point_all(:,2)/max(point_all(:,2)); % normalize energy
if ~exist("point_clouds","dir")
    mkdir("point_clouds");
end
writematrix(point_all, 'point_clouds/'+bin_filename+'.csv');
fclose(fid);
% clearvars;
