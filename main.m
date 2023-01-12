clc
clearvars
close all
%% 生成每个数据集的点云图
data_path="data/gesture-dataset/";
file_list = ["mid_pull_Raw_0"];
pfar=0.001;
for i=1:length(file_list)     
    bin_filename = file_list(i);
    filename = data_path+bin_filename+".bin";
    disp(filename);
    points_all = generate_point_cloud(bin_filename, pfar, data_path);
    disp(size(points_all,1));
    % TODO: 都处理完了，就是画不出图来
    plot_xyz_pointclouds(points_all, i, 2);
end
