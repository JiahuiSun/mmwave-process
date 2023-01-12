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


% Algorithm parameters
data_each_frame = samples*loop*Tx;

data_path="D:\ti\mmwave_RawData\";
i=8;
action_list=["center","bend","nod","back","left","right","leftleg","rightleg"];


Pfa=1e-2;
i = 1;
bin_filename=action_list(i)+"_d60_1";
filename = data_path+"sitting_posture_static_front_0223\"+bin_filename+".bin";
detout = point_cloud_generation(filename,Pfa);
% Plot xyz points

power = detout(:,4);
x_value=detout(:,8);
y_value=detout(:,9);
z_value=detout(:,10);
X = [x_value,y_value,z_value]; % Noisy 2-D circular data set

color_list=['r','g','b','y','m'];

idx = dbscan(X,0.1,5); % The default distance metric is Euclidean distance
idx_set=unique(idx);
for i=1:length(idx_set)
    X_local=X(find(idx==idx_set(i)),:);
    scatter3(X_local(:,1),X_local(:,2),X_local(:,3),20,color_list(i), 'filled');
    hold on
end
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis([-1.5, 1.5 -1.5 1.5 -1.5 1.5]);
title('DBSCAN Using Euclidean Distance Metric')