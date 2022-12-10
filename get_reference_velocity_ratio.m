clear; close all; clc;
logging_func("Get reference velocity ratio");
%% Explanation

%台車速度と基準PIV流速の比を算出するコード
% Compare Truck and PIV data
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Compare truck velocity and PIV velocity and compute characteristic velocity
% 
% Warning
%  This script needs following files in the same folder.
%   * compute_velocity_ratio.m
%   * logging_func.m
%

%% Setting
% Experiment PIV data
experiment_PIV_data_folder = "../PIV_xo350W45固定壁L1H1定常/result";

% Reference Truck data
reference_truck_data = "../加速度データ_xo350W100固定壁_新/Truck_data.mat";  %代表速度の基準にする速度データ

% Reference PIV data
reference_PIV_data = '../PIV_xo350W100固定壁/新データ/最適化/result/PIV_data.mat';

%% Compute velocity ratio of reference solution
compute_velocity_ratio(reference_truck_data,reference_PIV_data,experiment_PIV_data_folder);
visualize_velocity_ratio(reference_truck_data,reference_PIV_data,experiment_PIV_data_folder);
