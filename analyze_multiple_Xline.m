clear; close all; clc;
logging_func("Analyze multiple X line data");
%% Explanation

%
% analyze_multiple_Xline
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
PIV_data_filename(1) = "../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat";
PIV_data_filename(2) = "../PIV_xo350W60固定壁L1H1定常/result/PIV_data.mat";

% Truck data file name
Truck_data_filename(1) = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";
Truck_data_filename(2) = "../加速度データ_xo350W60固定壁L1H1定常/Truck_data.mat";

% Search pixel
search_x = [360, 500]; % [mm]

% Search time
search_t = [1.5, 3.0]; % [s]

%% Check input

if (length(PIV_data_filename) < 2 || length(PIV_data_filename) ~= length(Truck_data_filename))
    error("Please input multiple PIV data file name");
    return;
else
    searchFiles = length(PIV_data_filename);
end

for ii = 1:searchFiles
    if (~exist(PIV_data_filename(ii),'file'))
        error("Cannot find %s",PIV_data_filename(ii));
        return;
    end
end

for ii = 1:searchFiles
    if (~exist(Truck_data_filename(ii),'file'))
        error("Cannot find %s",Truck_data_filename(ii));
        return;
    end
end

if (length(search_x) ~= length(search_t))
    error("Search-x and Search-t should be same size");
    return;
end
numSearch = length(search_x);

%% Visualize data
% Get all u_filtered data
search_u_filtered = cell(searchFiles,1);
normalized_y = cell(searchFiles,1);
for ii = 1:numSearch
    for jj = 1:searchFiles
        [search_u_filtered(jj),normalized_y(jj),real_t,real_t_ind,real_x,real_x_ind] = compute_Xline(PIV_data_filename(jj),Truck_data_filename(jj),search_x(ii),search_t(ii));
    end

    % Visualize
    figure
    hold on
    for jj = 1:searchFiles
        plot(search_u_filtered{jj},normalized_y{jj},"-o");
    end
    hold off
    grid on
    xlabel("u/U");
    ylabel("y/H");
    %ylim([min(normalized_y{ii}),max(normalized_y{ii})]);
    title(sprintf("u in x line (t=%.3f[s], x=%.3f[m])",real_t,real_x));
end
