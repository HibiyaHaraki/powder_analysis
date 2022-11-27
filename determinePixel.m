function [neccesary_pixel_size,determined_pixel_size,convert_point] = determinePixel(sFreq,truck_mean_vx)

%
% determinePixel(sFreq,truck_mean_vx)
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% This script shows the appropriate PIV pixel size
%
% Inputs
%  - sFreq          : sampling frequency (scalar)
%  - truck_mean_vx  : Truck x-velocity (vector)
%

%% Check Inputs
% Sampling frequency shoud be scalar
if (~isscalar(sFreq))
    error("Sampling frequency shoud be scalar!!");
    return;
end

% Truck x-velocity data shoud be vector
if (~isvector(truck_mean_vx))
    error("Truck x-velocity data shoud be vector!!");
    return;
end
num_data = length(truck_mean_vx);

%% Constants
x_size = 0.1767; % [mm]
num_width_pixels = 1024; % [pixel]
size_factor = 2;
min_analysis_size = 2;

%% Determine appropriate pixel size
truck_movement = diff(cumtrapz(1/sFreq,[0;truck_mean_vx]));
length_per_pixel = x_size / num_width_pixels;
neccesary_pixel_size = size_factor.*truck_movement./length_per_pixel;

determined_pixel_size = 2.^(ceil(log(neccesary_pixel_size) ./ log(2) - 0.5));
determined_pixel_size(determined_pixel_size < min_analysis_size) = min_analysis_size;

%% Search converting point
diff_pixel_size = determined_pixel_size(2:end) - determined_pixel_size(1:end-1);
diff_pixel_size(diff_pixel_size ~= 0) = 1;
convert_point = [1;find(logical(diff_pixel_size))+1];

end