clear; close all; clc;
logging_func("Analyze 1 Pixel data");
%% Explanation

%速度変動の算出とフーリエ解析を実行するコード
% analyze_a_pixel
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%   * compute_velocity_u_flucutuation.m
%   * compute_velocity_v_flucutuation.m
%   * compute_characteristicVelocity.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W30固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W30固定壁L1H1定常/Truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W30固定壁L1H1定常/result';

% Search pixel
search_x = [400, 400, 500]; % [mm]
search_y = [5  ,20, 20]; % [mm]

% Time range for FFT
acceleration_fft_time_range = false;
acceleration_fft_start_time = 0.5; % [s]
acceleration_fft_stop_time = 1; % [s]

constVelocity_fft_time_range = false;
constVelocity_fft_start_time = 0.5; % [s]
constVelocity_fft_stop_time = 1; % [s]

% x-range in FFT result
x_range = [0 100];

%% Convert unit of search pixel
load(PIV_data_filename,"frontDistance");
search_x = search_x - frontDistance;

search_x = search_x ./ 1000;
search_y = search_y ./ 1000;

%% Compute characteristic velocity
logging_func("Compute characteristic velocity");
[characteristicVelocity_acceleration,characteristicVelocity_constVelocity] = compute_characteristicVelocity(PIV_data_filename,Truck_data_filename);

%% Compute Velocity fluctuation
logging_func("Compute Velocity fluctuation");
[velocity_u_fluctuation_acceleration,velocity_u_fluctuation_constVelocity,real_x,real_y] = compute_velocity_u_fluctuation(search_x,search_y,characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename);
[velocity_v_fluctuation_acceleration,velocity_v_fluctuation_constVelocity,     ~,     ~] = compute_velocity_v_fluctuation(search_x,search_y,characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename);

real_x = real_x .* 1000; % convert unit
real_y = real_y .* 1000; % convert unit
%% Show search point
num_search_pixel = length(search_x);
load(PIV_data_filename,"mean_x","mean_y","meanMap_u_filtered");
figure
pfig = pcolor(mean_x{end},mean_y{end},meanMap_u_filtered{end});
pfig.EdgeColor = 'none';
pfig.FaceColor = 'interp';
hold on
for ii = 1:num_search_pixel
    scatter(search_x,search_y,100,"r.");
end
hold off
axis equal
xlabel("x [m]");
ylabel("y [m]");

%% Visualize velocity fluctuation
load(Truck_data_filename,"sfreq");

% Acceleration
if (~isnan(velocity_u_fluctuation_acceleration))
    time_step = length(velocity_u_fluctuation_acceleration(1,:));
    [acceleration_fft_start_ind,acceleration_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        acceleration_fft_start_time, ...
        acceleration_fft_stop_time, ...
        acceleration_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_u_fluctuation_acceleration(ii,:));
        hold on
        xline((acceleration_fft_start_ind-1)/sfreq,"r-");
        xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in acceleration part");
    saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_acceleration.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"velocity_u_fluctuation_acceleration.png"));
end

if (~isnan(velocity_v_fluctuation_acceleration))
    time_step = length(velocity_v_fluctuation_acceleration(1,:));
    [acceleration_fft_start_ind,acceleration_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        acceleration_fft_start_time, ...
        acceleration_fft_stop_time, ...
        acceleration_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_v_fluctuation_acceleration(ii,:));
        hold on
        xline((acceleration_fft_start_ind-1)/sfreq,"r-");
        xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("v filtered in acceleration part");
    saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_acceleration.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"velocity_v_fluctuation_acceleration.png"));
end

% Const velocity
if (~isnan(velocity_u_fluctuation_constVelocity))
    time_step = length(velocity_u_fluctuation_constVelocity(1,:));
    [constVelocity_fft_start_ind,constVelocity_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        constVelocity_fft_start_time, ...
        constVelocity_fft_stop_time, ...
        constVelocity_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_u_fluctuation_constVelocity(ii,:));
        hold on
        xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
        xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"velocity_u_fluctuation_constVelocity.png"));
end

if (~isnan(velocity_v_fluctuation_constVelocity))
    time_step = length(velocity_v_fluctuation_constVelocity(1,:));
    [constVelocity_fft_start_ind,constVelocity_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        constVelocity_fft_start_time, ...
        constVelocity_fft_stop_time, ...
        constVelocity_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_v_fluctuation_constVelocity(ii,:));
        hold on
        xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
        xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"velocity_v_fluctuation_constVelocity.png"));
end

%% Frequency Analysis
% Acceleration
if (~isnan(velocity_u_fluctuation_acceleration))
    time_step = length(velocity_u_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    fft_acceleration = fft(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind),[],2);
    fft_acceleration = abs(fft_acceleration(:,2:end-1) ./ time_step);
    fft_frequency = (1:time_step-2)./time_step.*sfreq;

    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(fft_frequency,fft_acceleration(ii,:));
        %semilogx(fft_frequency,fft_acceleration(ii,:));
        grid on
        xlim(x_range);
        title(sprintf("FFT u result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of u filtered in acceleration part");
    saveas(gcf,output_folder+"/"+"FFT_u_acceleration.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_u_acceleration.png"));
end

if (~isnan(velocity_v_fluctuation_acceleration))
    time_step = length(velocity_v_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    fft_acceleration = fft(velocity_v_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind),[],2);
    tmp_show = abs(fft_acceleration ./ time_step);
    fft_acceleration = abs(fft_acceleration(:,2:end-1) ./ time_step);
    fft_frequency = (1:time_step-2)./time_step.*sfreq;

    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(fft_frequency,fft_acceleration(ii,:));
        %semilogx(fft_frequency,fft_acceleration(ii,:));
        grid on
        xlim(x_range);
        title(sprintf("FFT v result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of v filtered in acceleration part");
    saveas(gcf,output_folder+"/"+"FFT_v_acceleration.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_acceleration.png"));
end

% const velocity
if (~isnan(velocity_u_fluctuation_constVelocity))
    time_step = length(velocity_u_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
    fft_constVelocity = fft(velocity_u_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind),[],2);
    fft_constVelocity = abs(fft_constVelocity(:,2:end-1) ./ time_step);
    fft_frequency = (1:time_step-2)./time_step.*sfreq;

    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(fft_frequency,fft_constVelocity(ii,:));
        %semilogx(fft_frequency,fft_constVelocity(ii,:));
        grid on
        xlim(x_range);
        title(sprintf("FFT u result (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of u filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_u_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_u_constVelocity.png"));
end

if (~isnan(velocity_v_fluctuation_constVelocity))
    time_step = length(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
    fft_constVelocity = fft(velocity_v_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind),[],2);
    fft_constVelocity = abs(fft_constVelocity(:,2:end-1) ./ time_step);
    fft_frequency = (1:time_step-2)./time_step.*sfreq;

    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(fft_frequency,fft_constVelocity(ii,:));
        %semilogx(fft_frequency,fft_constVelocity(ii,:));
        grid on
        xlim(x_range);
        title(sprintf("FFT v result (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of v filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_v_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_constVelocity.png"));
end

%% Function
function [start_ind,stop_ind] = search_index_range(time,start_time,stop_time,ch)
if (ch)
    [~,start_ind] = min(abs(time - start_time));
    [~,stop_ind] = min(abs(time - stop_time));
else
    start_ind = 1;
    stop_ind = length(time);
end
end

