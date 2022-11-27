clear; close all; clc;
logging_func("Analyze X line data");
%% Explanation

%
% analyze_Xline
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Analyze a x-line data in mean map
% 
% Warning
%  This script needs following files in the same folder.
%   * logging_func.m
%   * compute_Xline.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";

% Search pixel
search_x = [360, 500]; % [mm]

% Search time
search_t = [1.5, 3.0]; % [s]

% Normalized option
NORMALIZE_OPTION = 1; % 0-object_width, 1-width

%% Check input

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (length(search_x) ~= length(search_t))       
    error("Search-x and Search-t should be same size");
    return;
end
numSearch = length(search_x);

%% Get normalize Parameter
load(PIV_data_filename,"object_width");

%% Visualize data
[search_u_filtered,normalized_y,real_t,real_t_ind,real_x,real_x_ind] = compute_Xline(PIV_data_filename,Truck_data_filename,search_x,search_t);

if (NORMALIZE_OPTION == 0)
    load(PIV_data_filename,"object_width");
    normalized_y = normalized_y./object_width;
else
    load(PIV_data_filename,"width");
    normalized_y = normalized_y./width;
end


for ii = 1:numSearch
    figure
    plot(search_u_filtered{ii},normalized_y{ii}./object_width,"-o");
    grid on
    xlabel("u/U");
    ylabel("y/H");
    ylim([min(normalized_y{ii}),max(normalized_y{ii})]);
    title(sprintf("u in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
    if (NORMALIZE_OPTION == 1)
        ylim([0,1]);
    end
end