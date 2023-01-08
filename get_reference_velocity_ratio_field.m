clear; close all; clc;
logging_func("Get characteristic velocity ratio field");
%% Explanation

%
% get_characteristic_velocity_ratio_field
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%

%% Setting

% Experiment PIV files
experiment_PIV_folder = "../PIV_xo350W60固定壁L1H1/result";


% Reference PIV data file
reference_PIV_filename = "../基準データ_W60固定壁/result/PIV_data.mat";

% Reference Truck velocity data file
reference_truck_filename = "../基準加速度_W60固定壁/truck_data.mat";

%% Check Inputs

if (~exist(experiment_PIV_folder,'dir'))
    error("Error: Cannot find %s",experiment_PIV_folder);
    return;
end


if (~exist(reference_PIV_filename,'file'))
    error("Error: Cannot find %s",reference_PIV_filename);
    return;
end

if (~exist(reference_truck_filename,'file'))
    error("Error: Cannot find %s",reference_truck_filename);
    return;
end

%% Compute characteristic velocity

compute_velocity_ratio_field(reference_truck_filename,reference_PIV_filename,experiment_PIV_folder);