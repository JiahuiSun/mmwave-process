clc;
clear;
close all;

radar_type='Infineon';

%% constants
c=3e8;

if strcmp(radar_type,'TI')
    fc=77e9;
    B=4e9;
    num_tx=2;
    num_rx=4;
    slope=36.837e12;
    chirp_duration=118.5e-6;
    fs=1e7;
    loops=128;
    N_elevation=2;
    theta=0;
else
    fc=24e9;
    B=200e6;
    c=3e8;
    num_tx=1;
    num_rx=2;
    up_chirp_duration=3e-4;
    chirp_duration=5e-4;
    fs=426666;
    loops=16;
    N_elevation=1;
    theta=0;
end

%% Derived Params
lambda=c/fc;
antenna_spacing=lambda/2;
% up_chirp_duration=B/slope;
slope=B/up_chirp_duration;
frame_duration=loops*chirp_duration;
N_azimuth=num_tx*num_rx;

%% Formulas
range_max=fs*c/(2*slope);
range_res=c/(2*B);

vel_max=lambda/(4*chirp_duration);
vel_res=lambda/(2*frame_duration);

angle_max=asin(lambda/(2*antenna_spacing))*180/pi;
azimuth_angle_res=lambda/(N_azimuth*antenna_spacing*cos(theta))*180/pi;
elevation_angle_res=lambda/(N_elevation*antenna_spacing*cos(theta))*180/pi;

%% Print Results
fprintf('Maximum range is %.2f m\n',range_max)
fprintf('Range resolution is %.2f cm\n',range_res*100)
fprintf('Maximum velocity is %.2f m/s\n',vel_max)
fprintf('Velocity resolution is %.2f m/s\n',vel_res)
fprintf('Maximum angle is %.2f degrees\n',angle_max)
fprintf('Azimuth angular resolution is %.2f degrees\n',azimuth_angle_res)
fprintf('Elevation angular resolution is %.2f degrees\n',elevation_angle_res)

