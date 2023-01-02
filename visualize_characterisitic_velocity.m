clear; close all; clc;
logging_func("Visualize characteristic velocity");
%% Explanation

%
% Visualize characteristic velocity
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Check and visualize truck velocity data
% 
% Warning
%  This script needs determinePixel.m
%

%% Setting

% PIV data file name
PIV_data_filename = "../PIV_xo350W100固定壁/新データ/最適化/result/PIV_data.mat";

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W100固定壁_新/Truck_data.mat";

%% Check Inputs

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

%% Load data
load(PIV_data_filename,'reference_velocity_ratio');
load(Truck_data_filename,'sfreq','truck_acceleration_mean_vx');

%% Determine Pixel
[neccesary_pixel_size,determined_pixel_size,convert_point] = determinePixel(sfreq,truck_acceleration_mean_vx);

%% Compute features
for ii = 1:length(convert_point)
    start_index = convert_point(ii);
    if (ii < length(convert_point))
        stop_index = convert_point(ii + 1) - 1;
    else
        stop_index = length(reference_velocity_ratio);
    end
    fprintf("Pixel Size: %d (%.3f - %.3f [s])\n",determined_pixel_size(start_index),(start_index-1)./sfreq,(stop_index)./sfreq);
    fprintf("  Max : %.3f\n", max(reference_velocity_ratio(start_index:stop_index)));
    fprintf("  Min : %.3f\n", min(reference_velocity_ratio(start_index:stop_index)));
    fprintf("  Mean: %.3f\n",mean(reference_velocity_ratio(start_index:stop_index)));
end

%% Visualize
figure
hold on
plot((0:(length(reference_velocity_ratio)-1))./sfreq,reference_velocity_ratio);
for ii = 1:length(convert_point)
    xline((convert_point(ii)-1)./sfreq,'k-',{sprintf('Pixel Size: %d',determined_pixel_size(convert_point(ii))),sprintf('Time: %.3f [s]',(convert_point(ii)-1)./sfreq)});
end
hold off
grid on
ylabel("characteristic velocity ratio");
xlabel("Time [s]");

figure
hold on
plot((0:(length(reference_velocity_ratio)-1))./sfreq,reference_velocity_ratio);
for ii = 1:length(convert_point)
    xline((convert_point(ii)-1)./sfreq,'k-',{sprintf('Pixel Size: %d',determined_pixel_size(convert_point(ii))),sprintf('Time: %.3f [s]',(convert_point(ii)-1)./sfreq)});
end
hold off
grid on
ylim([-inf max(reference_velocity_ratio(convert_point(2):end))]);
ylabel("characteristic velocity ratio");
xlabel("Time [s]");


