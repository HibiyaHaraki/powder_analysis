clear; close all; clc;
logging_func("Animate continuity");
%% Explanation

%
% animate_continutiy
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

% Record (1-on, 0-off)
RECORD_MODE = 0;
output_filename = "colormap_animation_W30固定壁L1H1定常.mp4";
frame_rate = 100; % Frame Rate

%% Check inputs

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

%% Load data
load(PIV_data_filename, ...
    "object_length", "object_width", ...
    "mean_x","mean_y", ...
    "meanMap_u_filtered","meanMap_v_filtered");

load(Truck_data_filename,"sfreq");
time_step = length(meanMap_u_filtered);


%% Create Animation
if (RECORD_MODE)
    video = VideoWriter(output_folder+"/"+output_filename,"MPEG-4");
    video.FrameRate = frame_rate;
    open(video);
end

ii = 1;
% Get interval
x_int = abs(mean_x{ii}(2) - mean_x{ii}(1));
y_int = abs(mean_y{ii}(2) - mean_y{ii}(1));
continuity_map = compute_continuity(meanMap_u_filtered{ii},meanMap_u_filtered{ii},x_int,y_int);

% Visualize
fig = figure;
fig.Position = [10,10,1900,500];
pfig = pcolor((mean_x{ii}*1000-object_length)/object_width,mean_y{ii}*1000/object_width,continuity_map);
pfig.EdgeColor = 'none';
pfig.FaceColor = 'interp';
colormap("jet");
colorbar
xlim([min((mean_x{end}*1000-object_length)/object_width,[],"all"),max((mean_x{end}*1000-object_length)/object_width,[],"all")]);
ylim([min(mean_y{end}*1000/object_width,[],"all"),max(mean_y{end}*1000/object_width,[],"all")]);
caxis([-300,300]);
axis equal
xlabel("x/H");
ylabel("y/H");
title(sprintf("t=%.3f[s]",0));

for ii = 1:time_step
    % Update figure
    pfig.XData = (mean_x{ii}*1000-object_length)/object_width;
    pfig.YData = mean_y{ii}*1000/object_width;
    pfig.ZData = ones(length(mean_y{ii}),length(mean_x{ii}));

    % Get interval
    x_int = abs(mean_x{ii}(2) - mean_x{ii}(1));
    y_int = abs(mean_y{ii}(2) - mean_y{ii}(1));

    % Compute continuity
    pfig.CData = compute_continuity(meanMap_u_filtered{ii},meanMap_v_filtered{ii},x_int,y_int);

    title(sprintf("t=%.3f[s]",(ii-1)/sfreq));
    drawnow;
    if (RECORD_MODE)
        frame = getframe(gcf);
        writeVideo(video,frame);
    end
end

if (RECORD_MODE)
    close(video);
    logging_func(sprintf("Output %s",output_folder+"/"+output_filename));
end