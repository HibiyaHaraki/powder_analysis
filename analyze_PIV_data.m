clear; close all; clc;
logging_func("Analyze PIV data");
%% Explanation

%
% analyze_PIV_data
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%   * compute_meanMap.m
%     * get_PIV_XAxis.m
%     * get_PIV_YAxis.m
%     * get_PIV_Data_meanMap.m
%     * logging_func.m
%   * compute_meanVlue.m
%

%% Setting

% Input files

file_name(1) = "../PIV_xo350W45固定壁L1H1定常/データ/データ1/data1combine.mat";
file_name(2) = "../PIV_xo350W45固定壁L1H1定常/データ/データ2/data2combine.mat";
file_name(3) = "../PIV_xo350W45固定壁L1H1定常/データ/データ3/data3combine.mat";
file_name(4) = "../PIV_xo350W45固定壁L1H1定常/データ/データ4/data4combine.mat";
file_name(5) = "../PIV_xo350W45固定壁L1H1定常/データ/データ5/data5combine.mat";

% Is there typevector_filtered?
typevector_existance = 1;

% Sampling frequency
sfreq = 1000; % [Hz]

% Front distance
frontDistance = 350; % [mm]

% Width
width = 45; % [mm]  % 縦軸の正規化値

% Object data
object_width = 15; % [mm]   % 横軸の正規化値
object_length = 15; % [mm]  % 横軸原点のずらし値

% Output folder name
output_folder_name = '../PIV_xo350W45固定壁L1H1定常/result';


%% Analyze
if (typevector_existance)
    compute_meanMap(file_name,output_folder_name,frontDistance,width,object_width,object_length);
else
    compute_meanMap_without_typevector(file_name,output_folder_name,frontDistance,width,object_width,object_length);
end
compute_meanValue(output_folder_name);
visualize_meanValue(output_folder_name,sfreq);
