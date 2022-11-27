clear; close all; clc;
logging_func("Create animation");
%% Explanation

%
% crate_animation
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W30固定壁L1H1定常/result';

% Record (1-on, 0-off)
RECORD_MODE = 0;
output_filename = "colormap_animation_W30固定壁L1H1定常.mp4"
frame_rate = 100; % Frame Rate

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

%% Import data
load(PIV_data_filename, ...
    "object_width","object_length", ...
    "mean_x","mean_y", ...
    "meanMap_u_filtered","meanMap_v_filtered");
load(Truck_data_filename,"sfreq");

time_step = (0:length(mean_x)-1)./sfreq;
logging_func(sprintf("Time step : %d",length(time_step)));

%% Search Max Value
maxValue = 0;
for ii = 1:length(time_step)
    tmp_maxValue = max(sqrt(meanMap_u_filtered{ii}.^2 + meanMap_v_filtered{ii}.^2),[],"all");
    if (tmp_maxValue > maxValue)
        maxValue = tmp_maxValue;
    end
end
logging_func(sprintf("Max magnitude flow velocity of PIV : %.3f [m/s]",maxValue));

%% Create animation and record

if (RECORD_MODE)
    video = VideoWriter(output_folder+"/"+output_filename,"MPEG-4");
    video.FrameRate = frame_rate;
    open(video);
end

fig = figure;
fig.Position = [10,10,1900,500];
pfig = pcolor((mean_x{1}*1000-object_length)/object_width,mean_y{1}*1000/object_width,sqrt(meanMap_u_filtered{1}.^2 + meanMap_v_filtered{1}.^2));
pfig.EdgeColor = 'none';
pfig.FaceColor = 'interp';
colormap("jet");
colorbar
axis equal
xlim([min((mean_x{end}*1000-object_length)/object_width,[],"all"),max((mean_x{end}*1000-object_length)/object_width,[],"all")]);
ylim([min(mean_y{end}*1000/object_width,[],"all"),max(mean_y{end}*1000/object_width,[],"all")]);
xlabel("x/H");
ylabel("y/H");
titletext = title(sprintf("t=%.3f [s]",time_step(1)));
caxis([0 maxValue]);
set(gca,"FontSize",15);

for ii = 1:length(time_step)
    pfig.XData = (mean_x{ii}*1000-object_length) / object_width;
    pfig.YData = mean_y{ii}*1000 / object_width;
    pfig.ZData = ones(length(mean_y{ii}),length(mean_x{ii}));
    pfig.CData = sqrt(meanMap_u_filtered{ii}.^2 + meanMap_v_filtered{ii}.^2);
    titletext.String = sprintf("t=%.3f [s]",time_step(ii));
    drawnow
    if (RECORD_MODE)
        frame = getframe(gcf);
        writeVideo(video,frame);
    end
end

if (RECORD_MODE)
    close(video);
    logging_func(sprintf("Output %s",output_folder+"/"+output_filename));
end

