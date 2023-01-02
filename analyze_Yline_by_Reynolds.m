clear; close all; clc;
logging_func("Analyze Y line data");
%% Explanation

%
% analyze_Yline
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Analyze a y-line data in mean map
% 
% Warning
%  This script needs following files in the same folder.
%   * logging_func.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W45固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W45固定壁L1H1定常/Truck_data.mat";

% Search pixel
search_y = [355, 500, 500]; % [mm]

% Search reynolds number
search_reynolds = [250, 500, 700]; % [s]

% Object width
object_width = 15; % [mm]

% Vorticity existance
VORTICITY = true;

%% Check input

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (length(search_y) ~= length(search_reynolds))
    error("Search-x and Search-t should be same size");
    return;
end
numSearch = length(search_y);

%% Import data
if (VORTICITY)
    load(PIV_data_filename, ...
        "mean_x","mean_y", ...
        "reference_velocity_ratio",...
        "meanMap_u_filtered","meanMap_v_filtered","meanMap_vorticity");
else
    load(PIV_data_filename, ...
        "mean_x","mean_y", ...
        "reference_velocity_ratio",...
        "meanMap_u_filtered","meanMap_v_filtered");
end
load(Truck_data_filename,"sfreq", ...
    'truck_acceleration_time_step','truck_constVelocity_time_step', ...
    'truck_acceleration_mean_vx','truck_constVelocity_mean_vx');

time_step = (0:length(mean_x)-1)./sfreq;

%% Convert unit of search pixel
search_y = search_y ./ 1000;

%% Get index of search reynolds
[Re_D,Re_H,Re_L] = solve_reynolds(PIV_data_filename,Truck_data_filename);
if (~isempty(search_reynolds))
    [~,real_t_ind] = min(abs(Re_H' - search_reynolds'),[],2);
    search_t = (real_t_ind' - 1)./sfreq;
    real_t = search_t;
end

%% Search y
real_y_ind = zeros(numSearch,1);
real_y = zeros(numSearch,1);
for ii = 1:numSearch
    [~,real_y_ind(ii)] = min(abs(mean_y{real_t_ind(ii)}' - search_y(ii)),[],2);
    real_y(ii) = mean_y{real_t_ind(ii)}(real_y_ind(ii));
end

%% Visualize meanMap_u_filtered in search y-line and FFT
for ii = 1:numSearch
    if (truck_acceleration_time_step < real_t_ind(ii))
        search_u_filtered = meanMap_u_filtered{real_t_ind(ii)}(real_y_ind(ii),:) ./ max(truck_constVelocity_mean_vx);
    else
        search_u_filtered = meanMap_u_filtered{real_t_ind(ii)}(real_y_ind(ii),:) ./ (reference_velocity_ratio(real_t_ind(ii)) * truck_acceleration_mean_vx(real_t_ind(ii)));
    end

    % Visualize u 
    figure
    plot(mean_x{real_t_ind(ii)}*1000/object_width,search_u_filtered,"-o");
    grid on
    xlabel("x/H");
    ylabel("u/U");
    xlim([0 max(mean_x{real_t_ind(ii)}*1000/object_width)]);
    %ylim([min(mean_x{real_t_ind(ii)}*1000/object_width),max(mean_x{real_t_ind(ii)}*1000/object_width)]);
    title(sprintf("u in y line (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));

    % FFT
    search_u_filtered_FFT = abs(fft(rmmissing(search_u_filtered))./length(search_u_filtered));
    search_u_filtered_FFT = search_u_filtered_FFT(1:round(length(search_u_filtered)/2)+1);
    search_u_filtered_FFT(2:end-1) = 2*search_u_filtered_FFT(2:end-1);
    freq = 1/(mean_x{real_t_ind(ii)}(2)*1000/object_width);
    f = freq*(0:(length(search_u_filtered_FFT)-1)) ./ length(search_u_filtered_FFT) ;

    figure
    plot(f(2:end),search_u_filtered_FFT(2:end));
    grid on
    xlabel("Freq [1/mm]");
    ylabel("Amplitude");
    xlim([min(f),max(f)]);
    title(sprintf("FFT u result (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));
end

%% Visualize meanMap_v_filtered in search y-line and FFT
for ii = 1:numSearch
    if (truck_acceleration_time_step < real_t_ind(ii))
        search_v_filtered = meanMap_v_filtered{real_t_ind(ii)}(real_y_ind(ii),:) ./ max(truck_constVelocity_mean_vx);
    else
        search_v_filtered = meanMap_v_filtered{real_t_ind(ii)}(real_y_ind(ii),:) ./ (reference_velocity_ratio(real_t_ind(ii)) * truck_acceleration_mean_vx(real_t_ind(ii)));
    end

    % Visualize v 
    figure
    plot(mean_x{real_t_ind(ii)}*1000/object_width,search_v_filtered,"-o");
    grid on
    xlabel("x/H");
    ylabel("u/U");
    xlim([0 max(mean_x{real_t_ind(ii)}*1000/object_width)]);
    %ylim([min(mean_x{real_t_ind(ii)}*1000/object_width),max(mean_x{real_t_ind(ii)}*1000/object_width)]);
    title(sprintf("v in y line (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));

    % FFT
    search_v_filtered_FFT = abs(fft(rmmissing(search_v_filtered))./length(search_v_filtered));
    search_v_filtered_FFT = search_v_filtered_FFT(1:round(length(search_v_filtered)/2)+1);
    search_v_filtered_FFT(2:end-1) = 2*search_v_filtered_FFT(2:end-1);
    freq = 1/(mean_x{real_t_ind(ii)}(2)*1000/object_width);
    f = freq*(0:(length(search_v_filtered_FFT)-1)) ./ length(search_v_filtered_FFT) ;

    figure
    plot(f(2:end),search_v_filtered_FFT(2:end));
    grid on
    xlabel("Freq [1/mm]");
    ylabel("Amplitude");
    xlim([min(f),max(f)]);
    title(sprintf("FFT v result (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));
end

%% Visualize meanMap_vorticity in search y-line and FFT
if (VORTICITY)
    for ii = 1:numSearch
        %{
        if (truck_acceleration_time_step < real_t_ind)
            search_vorticity = meanMap_vorticity{real_t_ind(ii)}(real_y_ind(ii),:) ./ max(truck_constVelocity_mean_vx);
        else
            search_vorticity = meanMap_vorticity{real_t_ind(ii)}(real_y_ind(ii),:) ./ (reference_velocity_ratio(real_t_ind(ii)) * truck_acceleration_mean_vx(real_t_ind(ii)));
        end
        %}
        search_vorticity = meanMap_vorticity{real_t_ind(ii)}(real_y_ind(ii),:);
    
        % Visualize vorticity 
        figure
        plot(mean_x{real_t_ind(ii)}*1000/object_width,search_vorticity,"-o");
        grid on
        xlabel("x/H");
        ylabel("u/U");
        xlim([0 max(mean_x{real_t_ind(ii)}*1000/object_width)]);
        %ylim([min(mean_x{real_t_ind(ii)}*1000/object_width),max(mean_x{real_t_ind(ii)}*1000/object_width)]);
        title(sprintf("vorticity in y line (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));
    
        % FFT
        search_vorticity_FFT = abs(fft(rmmissing(search_vorticity))./length(search_vorticity));
        search_vorticity_FFT = search_vorticity_FFT(1:round(length(search_vorticity)/2)+1);
        search_vorticity_FFT(2:end-1) = 2*search_vorticity_FFT(2:end-1);
        freq = 1/(mean_x{real_t_ind(ii)}(2)*1000/object_width);
        f = freq*(0:(length(search_vorticity_FFT)-1)) ./ length(search_vorticity_FFT) ;
    
        figure
        plot(f(2:end),search_vorticity_FFT(2:end));
        grid on
        xlabel("Freq [1/mm]");
        ylabel("Amplitude");
        xlim([min(f),max(f)]);
        title(sprintf("FFT vorticity result (t=%.3f[s], y=%.3f[m], Re=%.3f)",real_t(ii),real_y(ii),Re_H(real_t_ind(ii))));
    end
end