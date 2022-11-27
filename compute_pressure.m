clear; close all; clc;
logging_func("Compute pressure");
%% Explanation

%
% compute_pressure
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%   * solve_Poisson_rectangle.m
%   * compute_dudx.m
%   * compute_dvdy.m
%   * logging_func.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W30固定壁L1H1定常/result';

% Search time
search_t = 0:0.5:3;

% Vorticity existance
VORTICITY = true;

% Constant
rho = 1.201;

calculation_check = false;

%% Check input

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (~exist(output_folder,'dir'))
    error("Cannot find %s",output_folder);
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

%% Import data

load(PIV_data_filename, ...
"object_length","object_width","width", ...
"mean_x","mean_y", ...
"meanMap_u_filtered","meanMap_v_filtered");

load(Truck_data_filename,"sfreq");
time_step = length(meanMap_u_filtered);

%% Get index of search time
[~,ind_t] = min(abs((0:time_step-1)./sfreq - search_t),[],2);
real_t = (ind_t-1)./sfreq;

%% Compute pressure
for ii = 1:num_search_t
    logging_func(sprintf("search_t = %.3f (%d)",real_t(ii),ind_t(ii)));
    % Get interval
    x_int = abs(mean_x{ind_t(ii)}(2) - mean_x{ind_t(ii)}(1));
    y_int = abs(mean_y{ind_t(ii)}(2) - mean_y{ind_t(ii)}(1));

    % Compute Poisson equation
    pressure = solve_Poisson_rectangle(meanMap_u_filtered{ind_t(ii)},meanMap_v_filtered{ind_t(ii)},x_int,y_int,calculation_check);
    normalize_pressure = pressure./(1/2*rho*mean(meanMap_u_filtered{ind_t(ii)},"all","omitnan"));

    % Visualize
    figure
    pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,normalize_pressure);
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    %caxis([min(normalize_pressure,[],'all') max(normalize_pressure,[],'all')]);
    caxis([-0.5,0.5]);
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title(sprintf("Pressure (t=%.3f[s])",real_t(ii)));
end