clear; close all; clc;
logging_func("Visualize continuity");
%% Explanation

%
% visualize_continutiy
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create figure to check continuity
% 
% Warning
%  This script needs following files in the same folder.
%   * logging_func.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";

% Search time
search_t = [1,3];

%% Check inputs

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (~isvector(search_t))
    error("search_t should be vector!!");
    return;
end

if (length(search_t(1,:)) > 1)
    search_t = search_t';
end
num_search_t = length(search_t);

%% Load data
load(PIV_data_filename, ...
    "object_length", "object_width", ...
    "mean_x","mean_y", ...
    "meanMap_u_filtered","meanMap_v_filtered");

load(Truck_data_filename,"sfreq");
time_step = length(meanMap_u_filtered);

%% Get index of search time
[~,ind_t] = min(abs((0:time_step-1)./sfreq - search_t),[],2);
real_t = (ind_t-1)./sfreq;

%% Compute continuity

for ii = 1:num_search_t
    % Get interval
    x_int = abs(mean_x{ind_t(ii)}(2) - mean_x{ind_t(ii)}(1));
    y_int = abs(mean_y{ind_t(ii)}(2) - mean_y{ind_t(ii)}(1));
    [u_continuity_map,v_continuity_map,continuity_map] = compute_continuity(meanMap_u_filtered{ind_t(ii)},meanMap_v_filtered{ind_t(ii)},x_int,y_int);

    % Visualize
    figure
    pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,continuity_map);
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    caxis([min(continuity_map,[],'all') max(continuity_map,[],'all')]);
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title(sprintf("Velocity u (t=%.3f[s])",real_t(ii)));
end