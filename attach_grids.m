clear; close all; clc;
logging_func("Start attaching grids");
%% Explanation

%
% attach_grids
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Attach grids on the picture
% 
% Warning
%  This script needs following files in the same folder.
%   * logging_func.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1/Truck_data.mat";

% Image folder path
image_folder_path = "../PIV_xo350W30固定壁L1H1/データ/データ1/data1";
image_title = "data1";

% Output folder
%output_folder = '../PIV_xo350W30固定壁L1H1定常/result';

% Start time step of image
image_start_step = 529;

% Search t
search_t = [1,2];

% Search reynolds number
search_reynolds = [500,700];

% Ratio between real length [mm] and pixcel
standard_real_length = 10; % [mm]
standard_pixcel = 56.99; % [pxs]

% Range
standard_point_x = 1;
standard_point_y = 216;
standard_x_length = 1023;
standard_y_length = 170;

% Grid interval
grid_interval = 10; % [mm]

%% Check Inputs
if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (~exist(image_folder_path,'dir'))
    error("Cannot find %s",image_folder_path);
    return;
end

%{
if (~exist(output_folder,'dir'))
    error("Cannot find %s",output_folder);
    return;
end
%}

if (isempty(search_t) && isempty(search_reynolds))
    error("Please input search t or search reynolds");
    return;
end

if (~isempty(search_t))
    if (~isvector(search_t))
        error("search_t should be vector!!");
        return;
    end
end

if (~isempty(search_reynolds))
    if (~isvector(search_reynolds))
        error("search_t should be vector!!");
        return;
    end
end

if (~isscalar(standard_real_length))
    error("Please input scalar!!");
    return;
end

if (~isscalar(standard_pixcel))
    error("Please input scalar!!");
    return;
end

if (~isscalar(standard_point_x))
    error("Please input scalar!!");
    return;
end

if (~isscalar(standard_point_y))
    error("Please input scalar!!");
    return;
end

if (~isscalar(standard_x_length))
    error("Please input scalar!!");
    return;
end

if (~isscalar(standard_y_length))
    error("Please input scalar!!");
    return;
end

if (length(search_t(1,:)) > 1)
    search_t = search_t';
end

if (length(search_reynolds(1,:)) > 1)
    search_reynolds = search_reynolds';
end

%% Load
load(PIV_data_filename, ...
        "width","object_width","object_length", ...
        "mean_x","mean_y", ...
        "meanMap_u_filtered","meanMap_v_filtered",...
        "meanMap_vorticity");
load(Truck_data_filename,"sfreq");
time_step = length(meanMap_u_filtered);

%% Search t
[~,ind_t_t] = min(abs((0:time_step-1)./sfreq - search_t),[],2);

%% Search reynolds (Re_H)
[Re_D,Re_H,Re_L] = solve_reynolds(PIV_data_filename,Truck_data_filename);
[~,ind_t_reynolds] = min(abs(Re_H' - search_reynolds),[],2);

%% Convert to image file path
ind_t = sort(unique([ind_t_t,ind_t_reynolds]));
numFiles = length(ind_t);
logging_func(sprintf("Output %d files",numFiles));

image_file_path = string(zeros(1,numFiles));
output_file_path = string(zeros(1,numFiles));
for ii =1:numFiles
    image_file_path(ii)  = image_folder_path + "/" + image_title + sprintf("%06d",round(ind_t(ii))+image_start_step) + ".bmp";
    output_file_path(ii) = image_folder_path + "/" + sprintf("t%03d_reynolds%03d",round(ind_t(ii)),round(Re_H(ind_t(ii)))) + ".bmp";
end

%% Attach grids

converter = standard_pixcel / standard_real_length;

for ii = 1:numFiles
    logging_func(sprintf("Input %s",image_file_path(ii)));
    I = imread(image_file_path(ii));
    sz = size(I);
    line_list_x = (standard_point_x:converter*grid_interval:standard_point_x+standard_x_length) - standard_point_x;
    line_list_y = (standard_point_y+standard_y_length:-converter*grid_interval:standard_point_y) - standard_point_y;

    figure
    imshow(I(standard_point_y:standard_point_y+standard_y_length,standard_point_x:standard_point_x+standard_x_length));
    hold on;
    xline(line_list_x,"k-");
    yline(flip(line_list_y),"k-");
    hold off
    axis on
    xticks(line_list_x);
    xticklabels({0:grid_interval:grid_interval*(length(line_list_x)-1)});
    yticks(flip(line_list_y));
    yticklabels({flip(0:grid_interval:grid_interval*(length(line_list_y)-1))});
    saveas(gcf,output_file_path(ii));

    logging_func(sprintf("Output %s",output_file_path(ii)));
end


