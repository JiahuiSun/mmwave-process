clc
clear
close all

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

% raw_data=readDCA1000(data_path+"sitting_posture_static_front_0223\"+action+".bin");

% Algorithm parameters
data_each_frame = samples*loop*Tx;
% set_frame_number = size(raw_data,2)/(data_each_frame)
% frame_start = 1;
% frame_end = set_frame_number;
Is_Windowed = 1;% 1==> Windowing before doing range and angle fft
Is_plot_rangeDop = 0;
Is_plot_rangeAng = 1;

data_path="D:\ti\mmwave_RawData\";
action_list=["center","bend","nod","back","left","right","leftleg","rightleg"];
distance_list=[60,80]; % 2 domains
pkg_idx=[1 2];


distance=distance_list(1);
bin_filename=action_list(1)+"_d"+num2str(distance_list(1))+"_"+num2str(pkg_idx(1));
filename = data_path+"sitting_posture_static_front_0223\"+bin_filename+".bin";
disp(filename) %100 frames
Angdata_all = generate_ra_3dfft(filename,0);
Angdata_all=reshape(Angdata_all,40,128); %1*5120

disp(size(Angdata_all)) %(100,5120)
if Is_plot_rangeAng
%     plot_rangeAng(Angdata_all,rng_grid(1:40),agl_grid);
    imagesc(abs(Angdata_all))
end
title('Ô­Ê¼Í¼Ïñ');
