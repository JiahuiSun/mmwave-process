function [save_det_data_crop_all,sep_list] = generate_gifs(bin_filename,Pfa,data_path)

addpath(genpath('.\config'));
addpath(genpath('.\utils'));
addpath(genpath('.\modules\plot'));
addpath(genpath('.\modules\fft'));
addpath(genpath('.\modules\detection'));
addpath(genpath('.\modules\cluster'));

% parameter setting
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
Is_plot_rangeDop = 0;
Is_plot_rangeAng = 0;
Is_plot_pointCloud = 1;
gif=1;
frame_dur = 0.1;

filename=data_path + bin_filename + '.bin';

data_each_frame = samples*loop*Tx;

% Frame_num = data_length/data_each_frame;
% Frame_num = 364;

fid = fopen(filename);
save_det_data_crop_all=[];
sep_list=[];

pic_num=1;

for cnt=1:400
    data_frame=get_next_frame(fid);
    if length(data_frame) < data_each_frame
        break
    end
    data_chirp = [];

    for cj = 1:Tx*loop
        temp_data = data_frame(:, (cj-1)*samples+1:cj*samples);
        data_chirp(:,:,cj) = temp_data;
    end
    
    % separate the chirps for TDM-MIMO with 3 TXs
%     chirp_tx1 = data_chirp(:,:,1:3:end);%4*1024*128
    chirp_tx1 = data_chirp(:,:,1:3:end);
    chirp_tx3 = data_chirp(:,:,2:3:end);
    chirp_tx2 = data_chirp(:,:,3:3:end);

    % permutation with the format [samples, Rx, chirp]
    chirp_tx1 = permute(chirp_tx1, [2,1,3]);
    chirp_tx2 = permute(chirp_tx2, [2,1,3]);
    chirp_tx3 = permute(chirp_tx3, [2,1,3]);

    % Range FFT for chirps from 3txs
    [Rangedata_tx1] = fft_range(chirp_tx1,fft_Rang,Is_Windowed);
    [Rangedata_tx2] = fft_range(chirp_tx2,fft_Rang,Is_Windowed);
    [Rangedata_tx3] = fft_range(chirp_tx3,fft_Rang,Is_Windowed);

    % Doppler FFT
    Dopplerdata_tx1 = fft_doppler(Rangedata_tx1, fft_Vel, 0);
    Dopplerdata_tx2 = fft_doppler(Rangedata_tx2, fft_Vel, 0);
    Dopplerdata_tx3 = fft_doppler(Rangedata_tx3, fft_Vel, 0);
    Dopdata_sum = squeeze(mean(abs(Dopplerdata_tx1), 2));
    
    % Plot range-Doppler image
    % CFAR detector on Range-Velocity to detect targets 
    % Output format: [doppler index, range index(start from index 1), ...
    % cell power]
%     Pfa = 5e-3; % probability of false alarm
    [Resl_indx] = cfar_RV(Dopdata_sum, fft_Rang, num_crop, Pfa);
    detout = peakGrouping(Resl_indx);
    
    if Is_plot_rangeDop
%         saveas(gcf,'range-doppler\'+bin_filename+"_"+num2str(cnt)+'.jpg');
        if gif
            plot_rangeDop(Dopdata_sum,rng_grid,vel_grid, pic_num, detout(:, detout(2,:)>=20 & detout(2,:)<=100));
%             f=plot_xyz_pointclouds(save_det_data_crop, pic_num);
            F=getframe(gcf);
            I=frame2im(F);
            [I,map]=rgb2ind(I,256);
            if pic_num == 1
                imwrite(I,map,'gif\'+bin_filename+'_3-3_rangeDop.gif','gif','Loopcount',inf,'DelayTime',frame_dur);
            else
                imwrite(I,map,'gif\'+bin_filename+'_3-3_rangeDop.gif','gif','WriteMode','append','DelayTime',frame_dur);
            end
        else
            dir_name = "range-doppler\1009\"+bin_filename+"\";
            if ~exist(dir_name+pic_num+".jpg","file")
                if pic_num==1 && ~exist(dir_name,"dir") 
                    mkdir(dir_name);
                end
                plot_rangeDop(Dopdata_sum,rng_grid,vel_grid, pic_num);
%               saveas(gcf,dir_name+pic_num+".jpg");
                exportgraphics(gcf,dir_name+pic_num+".jpg");
                clear gcf;
                clear Dopdata_sum;
            else
                break;
            end
        end
    end
    if ~Is_plot_rangeAng && ~Is_plot_pointCloud
        pic_num = pic_num + 1;
        continue;
    end
    
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
    
%     Rangedata_merge = [Rangedata_tx1, Rangedata_tx3];
    Dopdata_merge=[Dopplerdata_tx1, Dopplerdata_tx3, Dopplerdata_tx2];
    
%     disp(size(Dopplerdata_tx1));
%     disp(size(Dopdata_merge));
%     % Angle FFT (Azimuth angle estimation)
%     Angdata = fft_angle(Dopdata_merge(:,1:2 * Rx, :),fft_Ang,Is_Windowed);
%     Angdata_crop = Angdata(num_crop + 1:fft_Rang - num_crop, :, :);
%     [Angdata] = Normalize(Angdata, max_value);
%     [Angdata_crop] = Normalize(Angdata_crop, max_value);

    % Plot range-angle (RA) image
    if Is_plot_rangeAng
        Angdata = fft_angle(Dopdata_merge(:,1:2 * Rx, :),fft_Ang,Is_Windowed);
        [Angdata] = Normalize(Angdata, max_value);
        if gif
            f=plot_rangeAng(Angdata,rng_grid,agl_grid, pic_num);
%             f=plot_xyz_pointclouds(save_det_data_crop, pic_num);
            F=getframe(gcf);
            I=frame2im(F);
            [I,map]=rgb2ind(I,256);
            if pic_num == 1
                imwrite(I,map,'gif\'+bin_filename+'_rangeAngle.gif','gif','Loopcount',inf,'DelayTime',frame_dur);
            else
                imwrite(I,map,'gif\'+bin_filename+'_rangeAngle.gif','gif','WriteMode','append','DelayTime',frame_dur);
            end
        end
%         saveas(gcf,'E:\大四下学期\毕设\figs_0223\side45\'+action+'_angle.jpg')
%         plot_rangeAng(Angdata_crop,rng_grid(num_crop+1:fft_Rang-num_crop),agl_grid);
    end
    
    %% point cloud generation
    if Is_plot_pointCloud
%     [x_vector, y_vector, z_vector]=naive_xyz(Dopdata_merge,detout,Tx,Rx,fft_Ang,Is_Windowed); 
    [x_vector, y_vector, z_vector]=advance_music(Dopdata_merge,detout,Tx,Rx,fft_Ang,Is_Windowed); 
    
    % Transform bin index to range/velocity
    Resel_vel = vel_grid(detout(1,:), 1);
    Resel_rng = rng_grid(detout(2,:), 1);

    x_vec=x_vector'.*Resel_rng;
    y_vec=y_vector'.*Resel_rng;
    z_vec=z_vector'.*Resel_rng;

    cur_frame_data=real([x_vec,y_vec,z_vec,Resel_vel]);

    [idx, isnoise] = myDBScan(cur_frame_data(:,1:3),0.7,2,1);
    idx_set=unique(idx);
%         disp(length(idx_set)-1)
    pc_body=cur_frame_data(find(idx~=0),:);
%         pc_body=cur_frame_data;

    pc_body=pc_body(pc_body(:,1)>-2 & pc_body(:,1)<2 &...  % x
            pc_body(:,2)>2 & pc_body(:,2)<4 &...  % y
            pc_body(:,3)>-2 & pc_body(:,3)<2, :);  % z
        
%     % Angle estimation for detected point clouds
%     Dopplerdata_merge = permute([Dopplerdata_tx1, Dopplerdata_tx3], [2, 1, 3]);
%     [Resel_agl, ~, rng_excd_list] = angle_estim_dets(detout, Dopplerdata_merge, fft_Vel, ...
%         fft_Ang, Rx, 2, num_crop);
    
    % Transform bin index to range/velocity/angle
%     Resel_agl_deg = agl_grid(1, Resel_agl)';
%     Resel_vel = vel_grid(detout(1,:), 1);
%     Resel_rng = rng_grid(detout(2,:), 1);
%     
%     x_vec=x_vector'.*Resel_rng;
%     y_vec=y_vector'.*Resel_rng;
%     z_vec=z_vector'.*Resel_rng;
%         
%     % save_det data format below
%     % [1 range bin, 2 velocity bin, 3 angle bin, 4 power, 5 range(m), 6 velocity (m/s), 7 angle(degree),8 x ,9 y ,10 z]
    save_det_data = [detout(2,:)', detout(1,:)', Resel_agl', detout(3,:)', ...
        Resel_rng, Resel_vel, Resel_agl_deg,x_vec,y_vec,z_vec];
% %     disp(save_det_data)
%     
%     % delete complex y value
%     valid_idx=find(imag(y_vec)==0);
%     save_det_data=save_det_data(valid_idx,:);

%     % filter out the points with range_bin within the crop region
%     if ~isempty(rng_excd_list)
%         save_det_data(rng_excd_list, :) = [];
%     end
    
%     save_det_data_crop=save_det_data(find(save_det_data(:,8)>-1.5 & save_det_data(:,8)<1.5 &...  % x
%         save_det_data(:,9)>2 & save_det_data(:,9)<4 &...  % y
%         save_det_data(:,10)>-2 & save_det_data(:,10)<2), :);  % z
%         save_det_data(:,6)<-0.05 & save_det_data(:,6)>0.05),:);
%     save_det_data_crop=save_det_data;

%     save_det_data_crop_all=[save_det_data_crop_all;save_det_data_crop];
%     sep_list=[sep_list,size(save_det_data_crop,1)];
% 
%     point_cloud_feat=save_det_data_crop(:,[4,8,9,10]);% power,x,y,z
%     point_cloud_feat=sortrows(point_cloud_feat,1,'descend');
%     point_num=size(point_cloud_feat,1);
%     if point_num>=30
%         point_cloud_feat_crop=point_cloud_feat(1:30,:); %30*4
%     else
%         point_cloud_feat_crop=zeros(30,4);
%         point_cloud_feat_crop(1:point_num,:)=point_cloud_feat;
%     end
%     
%     point_cloud_feat_crop=reshape( point_cloud_feat_crop,1,[]); %1*120
%     point_cloud_feat_all=[point_cloud_feat_all;point_cloud_feat_crop];
    
    if gif
        f=plot_xyz_pointclouds(pc_body, pic_num, 4);
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if pic_num == 1
            imwrite(I,map,'gif\'+bin_filename+'_3-15_PC.gif','gif','Loopcount',inf,'DelayTime',frame_dur);
        else
            imwrite(I,map,'gif\'+bin_filename+'_3-15_PC.gif','gif','WriteMode','append','DelayTime',frame_dur);
        end
    end
    
    end
    pic_num = pic_num + 1;
    clear data_frame;
%     break
end
fclose(fid);