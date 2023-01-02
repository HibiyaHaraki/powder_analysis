clear; close all; clc;
logging_func("Analyze frequency data in y-line");
%% Explanation

% 
% analyze_frequency_yline
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Analyze frequency
% 
% Warning
%  This script needs following files in the same folder.
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W45固定壁L1H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W45固定壁L1H1定常/Truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W45固定壁L1H1定常/result';

% Search pixel
search_x = [20,30]; % [mm]

% Search frequency
search_freq = 12; % [Hz]

% Time range for FFT
acceleration_fft_time_range = false;
acceleration_fft_start_time = 0.5; % [s]
acceleration_fft_stop_time = 1; % [s]

constVelocity_fft_time_range = true;
constVelocity_fft_start_time = 1; % [s]
constVelocity_fft_stop_time = 3; % [s]

% x-range in FFT result
x_range = [0 100];

%% Convert unit of search pixel
search_x = search_x ./ 1000;

load(PIV_data_filename,"mean_y");
search_y = mean_y{1};

num_search_x = length(search_x);

%% Compute characteristic velocity
logging_func("Compute characteristic velocity");
[characteristicVelocity_acceleration,characteristicVelocity_constVelocity] = compute_characteristicVelocity(PIV_data_filename,Truck_data_filename);

%% Compute Velocity fluctuation
load(Truck_data_filename,"sfreq");

% Create array for data
num_search_y = length(search_y);
amplitude_record_u_acceleration  = nan(num_search_x,num_search_y);
amplitude_record_v_acceleration  = nan(num_search_x,num_search_y);
amplitude_record_u_constVelocity = nan(num_search_x,num_search_y);
amplitude_record_v_constVelocity = nan(num_search_x,num_search_y);

power_spectrum_u_acceleration  = nan(num_search_x,num_search_y);
power_spectrum_v_acceleration  = nan(num_search_x,num_search_y);
power_spectrum_u_constVelocity = nan(num_search_x,num_search_y);
power_spectrum_v_constVelocity = nan(num_search_x,num_search_y);

for jj = 1:num_search_x
    logging_func(sprintf("Iteration step %d / %d",jj,num_search_x));

    parfor ii = 1:num_search_y
    
        acceleration_fft_start_ind  = NaN;
        acceleration_fft_stop_ind   = NaN;
        constVelocity_fft_start_ind = NaN;
        constVelocity_fft_stop_ind  = NaN;
    
        [velocity_u_fluctuation_acceleration,velocity_u_fluctuation_constVelocity,real_x,real_y] = compute_velocity_u_fluctuation(search_x(jj),search_y(ii),characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename);
        [velocity_v_fluctuation_acceleration,velocity_v_fluctuation_constVelocity,     ~,     ~] = compute_velocity_v_fluctuation(search_x(jj),search_y(ii),characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename);
    
        real_x = real_x .* 1000; % convert unit
        real_y = real_y .* 1000; % convert unit
        
        % Acceleration
        if (~isnan(velocity_u_fluctuation_acceleration))
            time_step = length(velocity_u_fluctuation_acceleration(1,:));
            [acceleration_fft_start_ind,acceleration_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
                acceleration_fft_start_time, ...
                acceleration_fft_stop_time, ...
                acceleration_fft_time_range);
        end
        
        if (~isnan(velocity_v_fluctuation_acceleration))
            time_step = length(velocity_v_fluctuation_acceleration(1,:));
            [acceleration_fft_start_ind,acceleration_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
                acceleration_fft_start_time, ...
                acceleration_fft_stop_time, ...
                acceleration_fft_time_range);
        end
        
        % Const velocity
        if (~isnan(velocity_u_fluctuation_constVelocity))
            time_step = length(velocity_u_fluctuation_constVelocity(1,:));
            [constVelocity_fft_start_ind,constVelocity_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
                constVelocity_fft_start_time, ...
                constVelocity_fft_stop_time, ...
                constVelocity_fft_time_range);
        end
        
        if (~isnan(velocity_v_fluctuation_constVelocity))
            time_step = length(velocity_v_fluctuation_constVelocity(1,:));
            [constVelocity_fft_start_ind,constVelocity_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
                constVelocity_fft_start_time, ...
                constVelocity_fft_stop_time, ...
                constVelocity_fft_time_range);
        end
        
        % Frequency Analysis
        % Acceleration
        if (~isnan(velocity_u_fluctuation_acceleration))
            % FFT
            time_step = length(velocity_u_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind));
            fft_acceleration = fft(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind),[],2);
            fft_acceleration = abs(fft_acceleration(:,2:end-1) ./ time_step);
            fft_frequency = (1:time_step-2)./time_step.*sfreq;
            [~,tmp_ind] = min(abs(fft_frequency - search_freq));
            amplitude_record_u_acceleration(jj,ii) = fft_acceleration(tmp_ind);
    
            % Power Spectrum
            [pxx,ps_freq] = periodogram(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind), ...
                rectwin(length(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind))), ...
                length(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind)), ...
                sfreq)
            [~,tmp_ind] = min(abs(ps_freq - search_freq));
            power_spectrum_u_acceleration(jj,ii) = pxx(tmp_ind);
        end
        
        if (~isnan(velocity_v_fluctuation_acceleration))
            time_step = length(velocity_v_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind));
            fft_acceleration = fft(velocity_v_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind),[],2);
            tmp_show = abs(fft_acceleration ./ time_step);
            fft_acceleration = abs(fft_acceleration(:,2:end-1) ./ time_step);
            fft_frequency = (1:time_step-2)./time_step.*sfreq;
            [~,tmp_ind] = min(abs(fft_frequency - search_freq));
            amplitude_record_v_acceleration(jj,ii) = fft_acceleration(tmp_ind);
    
            % Power Spectrum
            [pxx,ps_freq] = periodogram(velocity_v_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind), ...
                rectwin(length(velocity_v_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind))), ...
                length(velocity_v_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind)), ...
                sfreq)
            [~,tmp_ind] = min(abs(ps_freq - search_freq));
            power_spectrum_v_acceleration(jj,ii) = pxx(tmp_ind);
        end
        
        % const velocity
        if (~isnan(velocity_u_fluctuation_constVelocity))
            time_step = length(velocity_u_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
            fft_constVelocity = fft(velocity_u_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind),[],2);
            fft_constVelocity = abs(fft_constVelocity(:,2:end-1) ./ time_step);
            fft_frequency = (1:time_step-2)./time_step.*sfreq;
            [~,tmp_ind] = min(abs(fft_frequency - search_freq));
            amplitude_record_u_constVelocity(jj,ii) = fft_constVelocity(tmp_ind);
    
            % Power Spectrum
            [pxx,ps_freq] = periodogram(velocity_u_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind), ...
                rectwin(length(velocity_u_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind))), ...
                length(velocity_u_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind)), ...
                sfreq)
            [~,tmp_ind] = min(abs(ps_freq - search_freq));
            power_spectrum_u_constVelocity(jj,ii) = pxx(tmp_ind);
        end
        
        if (~isnan(velocity_v_fluctuation_constVelocity))
            time_step = length(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
            fft_constVelocity = fft(velocity_v_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind),[],2);
            fft_constVelocity = abs(fft_constVelocity(:,2:end-1) ./ time_step);
            fft_frequency = (1:time_step-2)./time_step.*sfreq;
            [~,tmp_ind] = min(abs(fft_frequency - search_freq));
            amplitude_record_v_constVelocity(jj,ii) = fft_constVelocity(tmp_ind);
    
            % Power Spectrum
            [pxx,ps_freq] = periodogram(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind), ...
                rectwin(length(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind))), ...
                length(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind)), ...
                sfreq)
            [~,tmp_ind] = min(abs(ps_freq - search_freq));
            power_spectrum_v_constVelocity(jj,ii) = pxx(tmp_ind);
        end
    end
end

%% Visualize FFT result
% Create legend
for ii = 1:num_search_x
    legend_str(ii) = sprintf("x=%.3f mm",search_x(ii)*1000); 
end

% u acceleration
figure
hold on
for ii = 1:num_search_x
    plot(search_y.*1000,amplitude_record_u_acceleration(ii,:));
end
hold off
grid on
legend(legend_str);
xlabel("y");
ylabel("Amplitude");
title(sprintf("accelerating u FFT (f=%.3f)",search_freq));

% v acceleration
figure
hold on
for ii = 1:num_search_x
    plot(search_y.*1000,amplitude_record_v_acceleration(ii,:));
end
hold off
grid on
legend(legend_str);
xlabel("y");
ylabel("Amplitude");
title(sprintf("accelerating v FFT (f=%.3f)",search_freq));

% u constVelocity
figure
hold on
for ii = 1:num_search_x
    semilogy(search_y.*1000,amplitude_record_u_constVelocity(ii,:));
end
hold off
grid on
legend(legend_str);
gca.YScale='log';
xlabel("y");
ylabel("Amplitude");
title(sprintf("Constant Velocity u FFT (f=%.3f)",search_freq));

% v constVelocity
figure
hold on
for ii = 1:num_search_x
    semilogy(search_y.*1000,amplitude_record_v_constVelocity(ii,:));
end
hold off
grid on
legend(legend_str);
gca.YScale='log';
xlabel("y");
ylabel("Amplitude");
title(sprintf("Constant Velocity v FFT (f=%.3f)",search_freq));

%% Visualize power spectrum
% u acceleration
figure
hold on
for ii = 1:num_search_x
    plot(search_y.*1000,power_spectrum_u_acceleration(ii,:));
end
hold off
grid on
legend(legend_str);
xlabel("y");
ylabel("Amplitude");
title(sprintf("accelerating u Power spectrum (f=%.3f)",search_freq));

% v acceleration
figure
hold on
for ii = 1:num_search_x
    plot(search_y.*1000,power_spectrum_v_acceleration(ii,:));
end
hold off
grid on
legend(legend_str);
xlabel("y");
ylabel("Amplitude");
title(sprintf("accelerating v Power spectrum (f=%.3f)",search_freq));

% u constVelocity
figure
hold on
for ii = 1:num_search_x
    plot(search_y.*1000,power_spectrum_u_constVelocity(ii,:));
end
hold off
grid on
legend(legend_str);
gca.YScale='log';
xlabel("y");
ylabel("Amplitude");
title(sprintf("Constant Velocity u Power spectrum (f=%.3f)",search_freq));

% v constVelocity
figure
hold on
for ii = 1:num_search_x
    semilogy(search_y.*1000,power_spectrum_v_constVelocity(ii,:));
end
hold off
ax = gca;
grid on
legend(legend_str);
ax.YScale='log';
xlabel("y");
ylabel("Amplitude");
title(sprintf("Constant Velocity v Power spectrum (f=%.3f)",search_freq));

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
