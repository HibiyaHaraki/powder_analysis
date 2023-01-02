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
PIV_data_filename = '../PIV_xo350W60移動壁L3H1定常/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W60移動壁L3H1定常/truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W60移動壁L3H1定常/result';

% Search pixel
search_x = [405, 435, 465, 495]; % [mm]
search_y = [58, 58, 58, 58]; % [mm]

% Time range for FFT
acceleration_fft_time_range = false;
acceleration_fft_start_time = 1.3; % [s]
acceleration_fft_stop_time = 3.2; % [s]

constVelocity_fft_time_range = false;
constVelocity_fft_start_time = 1.0; % [s]
constVelocity_fft_stop_time = 2.0; % [s]

% x-range in FFT result
x_range = [0 100];

% Number of frequencies to compute Strouhal number
strouhal_freqs = 5;

% autocorrelation
autocoorelation_lag = 1.0*10^(-1);

% Seting for saving figures
% Acceleration
SAVE_position_acceleration = [ 
    true, true, true, true,
    true, true, true, true
    ];

SAVE_u_fluctuation_acceleration = true;
SAVE_v_fluctuation_acceleration = true;
SAVE_absVelocity_fluctuation_acceleration = true;

SAVE_u_fft_acceleration = false;
SAVE_v_fft_acceleration = false;
SAVE_absVelocity_fft_acceleration = false;

% Constant Velocity
SAVE_position_constVelocity = [ ...
    true, true, true, true
    ];

SAVE_u_fluctuation_constVelocity = false;
SAVE_v_fluctuation_constVelocity = false;
SAVE_absVelocity_fluctuation_constVelocity = false;

SAVE_u_fft_constVelocity = false;
SAVE_v_fft_constVelocity = false;
SAVE_absVelocity_fft_constVelocity = false;

%% Check Inputs
if (length(search_x)  ~= length(SAVE_position_acceleration))
    error("Please input same number of the data for saving figures of acceleration part");
end

if (length(search_x)  ~= length(SAVE_position_constVelocity))
    error("Please input same number of the data for saving figures of constant velocity part");
end

%% Convert unit of search pixel
load(PIV_data_filename,"frontDistance","");
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
[velocity_absVelocity_fluctuation_acceleration,velocity_absVelocity_fluctuation_constVelocity,     ~,     ~] = compute_velocity_absVelocity_fluctuation(search_x,search_y,characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename);

real_x = real_x .* 1000; % convert unit
real_y = real_y .* 1000; % convert unit

%% Show search point
num_search_pixel = length(search_x);
load(PIV_data_filename,"mean_x","mean_y","meanMap_u_filtered","object_width","reference_velocity_ratio");
figure
pfig = pcolor(1000*mean_x{end},1000*mean_y{end},meanMap_u_filtered{end});
pfig.EdgeColor = 'none';
pfig.FaceColor = 'interp';
hold on
for ii = 1:num_search_pixel
    scatter(1000*search_x,1000*search_y,100,"r.");
end
hold off
axis equal
xlabel("x [mm]");
ylabel("y [mm]");

%% Visualize velocity fluctuation
load(Truck_data_filename,"sfreq",'truck_constVelocity_mean_vx');

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
        yline(0,"b-");
        yline(1,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        ylim([-0.6, 1.5]);
        xlabel("t[s]");
        ylabel("u/U");
        title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in acceleration part");
    if (SAVE_u_fluctuation_acceleration)
        saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_u_fluctuation_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_u_fluctuation_acceleration(ii,:));
                hold on
                xline((acceleration_fft_start_ind-1)/sfreq,"r-");
                xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
                yline(0,"b-");
                yline(1,"r-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([-0.6, 1.5]);
                xlabel("x/H");
                ylabel("u/U");
                title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
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
        yline(0,"b-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        ylim([-0.7, 0.7]);
        xlabel("t[s]");
        ylabel("v/U");
        title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("v filtered in acceleration part");
    if (SAVE_v_fluctuation_acceleration)
        saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_v_fluctuation_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_v_fluctuation_acceleration(ii,:));
                hold on
                xline((acceleration_fft_start_ind-1)/sfreq,"r-");
                xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
                yline(0,"b-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([-0.7, 0.7]);
                xlabel("x/H");
                ylabel("v/U");
                title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
end

if (~isnan(velocity_absVelocity_fluctuation_acceleration))
    time_step = length(velocity_absVelocity_fluctuation_acceleration(1,:));
    [acceleration_fft_start_ind,acceleration_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        acceleration_fft_start_time, ...
        acceleration_fft_stop_time, ...
        acceleration_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_absVelocity_fluctuation_acceleration(ii,:));
        hold on
        xline((acceleration_fft_start_ind-1)/sfreq,"r-");
        xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
        yline(1,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        ylim([0, 1.5]);
        xlabel("t[s]");
        ylabel("Velocity/U");
        title(sprintf("Velocity absolute velocity fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("Absolute velocity filtered in acceleration part");
    if (SAVE_absVelocity_fluctuation_acceleration)
        saveas(gcf,output_folder+"/"+"velocity_absVelocity_fluctuation_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_absVelocity_fluctuation_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_absVelocity_fluctuation_acceleration(ii,:));
                hold on
                xline((acceleration_fft_start_ind-1)/sfreq,"r-");
                xline((acceleration_fft_stop_ind-1)/sfreq,"r-");
                yline(1,"r-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([0, 1.5]);
                xlabel("t[s]");
                ylabel("Velocity/U");
                title(sprintf("Velocity absVelocity fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_absVelocity_fluctuation_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
end
%------------------------------------------------------------------------------------------------------------------------------------------
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
        xlabel("x/H");
        ylabel("u/U");
        title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in constant velocity part");
    if (SAVE_u_fluctuation_constVelocity)
        saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_constVelocity.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_u_fluctuation_constVelocity.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_constVelocity(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_u_fluctuation_constVelocity(ii,:));
                hold on
                xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
                xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
                yline(0,"b-");
                yline(1,"r-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([-0.6, 1.5]);
                xlabel("t[s]");
                ylabel("u/U");
                title(sprintf("Velocity u fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_u_fluctuation_constVelocity_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
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
        xlabel("t[s]");
        ylabel("u/U");
        title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("u filtered in constant velocity part");
    if (SAVE_v_fluctuation_constVelocity)
        saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_constVelocity.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_v_fluctuation_constVelocity.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_constVelocity(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_v_fluctuation_constVelocity(ii,:));
                hold on
                xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
                xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
                yline(0,"b-");
                yline(1,"r-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([-0.6, 1.5]);
                xlabel("t[s]");
                ylabel("u/U");
                title(sprintf("Velocity v fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_v_fluctuation_constVelocity_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
end

if (~isnan(velocity_absVelocity_fluctuation_constVelocity))
    time_step = length(velocity_absVelocity_fluctuation_constVelocity(1,:));
    [constVelocity_fft_start_ind,constVelocity_fft_stop_ind] = search_index_range((0:time_step-1)./sfreq, ...
        constVelocity_fft_start_time, ...
        constVelocity_fft_stop_time, ...
        constVelocity_fft_time_range);
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot((0:time_step-1)./sfreq,velocity_absVelocity_fluctuation_constVelocity(ii,:));
        hold on
        xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
        xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
        hold off
        grid on
        xlim([0,(time_step-1)./sfreq]);
        xlabel("t[s]");
        ylabel("v/U");
        title(sprintf("Velocity absolute velocity fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("Absolute velocity filtered in constant velocity part");
    if (SAVE_absVelocity_fluctuation_constVelocity)
        saveas(gcf,output_folder+"/"+"velocity_absVelocity_fluctuation_constVelocity.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"velocity_absVelocity_fluctuation_constVelocity.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_constVelocity(ii))
                figure("visible","off");
                plot((0:time_step-1)./sfreq,velocity_absVelocity_fluctuation_constVelocity(ii,:));
                hold on
                xline((constVelocity_fft_start_ind-1)/sfreq,"r-");
                xline((constVelocity_fft_stop_ind-1)/sfreq,"r-");
                yline(0,"b-");
                yline(1,"r-");
                hold off
                grid on
                xlim([0,(time_step-1)./sfreq]);
                ylim([-0.6, 1.5]);
                xlabel("x/H");
                ylabel("Velocity/U");
                title(sprintf("Absolute Velocity fluctuation (x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_absVelocity_fluctuation_constVelocity_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT u result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of u filtered in acceleration part");
    if (SAVE_u_fft_acceleration)
        saveas(gcf,output_folder+"/"+"FFT_u_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"FFT_u_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot(fft_frequency,fft_acceleration(ii,:));
                %semilogx(fft_frequency,fft_acceleration(ii,:));
                grid on
                xlim(x_range);
                xlabel("f[Hz]");
                ylabel("amplitude");
                title(sprintf("FFT u result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_u_fft_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (acceleration u)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_acceleration(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        U = mean(characteristicVelocity_acceleration(acceleration_fft_start_ind:acceleration_fft_stop_ind),"all");
        strouhal_number = locs(index) .* (object_width ./ U) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n",jj,strouhal_number(jj),locs(index(jj)),U);
        end
     end
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT v result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of v filtered in acceleration part")
    if (SAVE_v_fft_acceleration)
        saveas(gcf,output_folder+"/"+"FFT_v_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot(fft_frequency,fft_acceleration(ii,:));
                %semilogx(fft_frequency,fft_acceleration(ii,:));
                grid on
                xlim(x_range);
                xlabel("f[Hz]");
                ylabel("amplitude");
                title(sprintf("FFT v result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_v_fft_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end

     % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (acceleration v)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_acceleration(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        U = mean(characteristicVelocity_acceleration(acceleration_fft_start_ind:acceleration_fft_stop_ind),"all");
        strouhal_number = locs(index) .* (object_width ./ U) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n",jj,strouhal_number(jj),locs(index(jj)),U);
        end
     end
end

if (~isnan(velocity_u_fluctuation_acceleration) & ~isnan(velocity_v_fluctuation_acceleration))
    time_step = length(velocity_u_fluctuation_acceleration(1,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    fft_acceleration = fft(sqrt(velocity_u_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind).^2+velocity_v_fluctuation_acceleration(:,acceleration_fft_start_ind:acceleration_fft_stop_ind).^2),[],2);
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT absolute velocity result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of absolute velocity filtered in acceleration part");
    if (SAVE_absVelocity_fft_acceleration)
        saveas(gcf,output_folder+"/"+"FFT_absVelocity_acceleration.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"FFT_absVelocity_acceleration.png"));
        for ii =1:num_search_pixel
            if (SAVE_position_acceleration(ii))
                figure("visible","off");
                plot(fft_frequency,fft_acceleration(ii,:));
                %semilogx(fft_frequency,fft_acceleration(ii,:));
                grid on
                xlim(x_range);
                xlabel("f[Hz]");
                ylabel("amplitude");
                title(sprintf("FFT absolute Velocity result (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
                saveas(gcf,output_folder+"/"+"velocity_absVelocity_fft_acceleration_"+real_x(ii,1)+"_"+real_y(ii,1)+".png");
            end
        end
    end

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (acceleration absolute velocity)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_acceleration(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        U = mean(characteristicVelocity_acceleration(acceleration_fft_start_ind:acceleration_fft_stop_ind),"all");
        strouhal_number = locs(index) .* (object_width ./ U) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n",jj,strouhal_number(jj),locs(index(jj)),U);
        end
     end
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT u result (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of u filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_u_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_u_constVelocity.png"));

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (constant velocity u)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_constVelocity(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        strouhal_number = locs(index) .* ...
            (object_width ./ (reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1))) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n", ...
                jj,strouhal_number(jj),locs(index(jj)),reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1));
        end
     end
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT v result (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of v filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_v_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_constVelocity.png"));

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (constant velocity v)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_constVelocity(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        strouhal_number = locs(index) .* ...
            (object_width ./ (reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1))) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n", ...
                jj,strouhal_number(jj),locs(index(jj)),reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1));
        end
     end
end

if (~isnan(velocity_u_fluctuation_constVelocity) & ~isnan(velocity_v_fluctuation_constVelocity))
    time_step = length(velocity_v_fluctuation_constVelocity(1,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
    fft_constVelocity = fft(sqrt(velocity_v_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind).^2+velocity_u_fluctuation_constVelocity(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind).^2),[],2);
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
        xlabel("f[Hz]");
        ylabel("amplitude");
        title(sprintf("FFT absolute velocity result (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of absolute velocity filtered in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_absVelocity_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_absVelocity_constVelocity.png"));

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (constant velocity absolute velocity)\n"));
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_constVelocity(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        strouhal_number = locs(index) .* ...
            (object_width ./ (reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1))) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n", ...
                jj,strouhal_number(jj),locs(index(jj)),reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1));
        end
     end
end

%% Analyze vorticity
load(PIV_data_filename,"meanMap_vorticity");
acceleration_timeStep = length(velocity_u_fluctuation_acceleration);
vorticity_acceleration = meanMap_vorticity(1:acceleration_timeStep);
if (~isnan(velocity_u_fluctuation_constVelocity))
    vorticity_constVelocity = meanMap_vorticity(acceleration_timeStep+1:end);
    constVelocity_timeStep = length(vorticity_constVelocity);
end

% Acceleration
logging_func("Analyze Vorticity (acceleration)");
time_step = length(meanMap_vorticity(acceleration_fft_start_ind:acceleration_fft_stop_ind));
if (~isnan(acceleration_timeStep))
    figure
    tiledlayout('flow')
    vorticity_data = zeros(num_search_pixel,length(vorticity_acceleration));
    for ii = 1:num_search_pixel
        ind_x = zeros(1,length(vorticity_acceleration));
        ind_y = zeros(1,length(vorticity_acceleration));
        % Search index
        for jj = 1:length(vorticity_acceleration)
            [~,ind_x(jj)] = min(abs(mean_x{jj} - search_x(ii)));
            [~,ind_y(jj)] = min(abs(mean_y{jj} - search_y(ii)));
            vorticity_data(ii,jj) = vorticity_acceleration{jj}(ind_y(jj),ind_x(jj));
        end

        nexttile
        hold on
        plot((0:(length(vorticity_data)-1))./sfreq,vorticity_data(ii,:));
        xline(([acceleration_fft_start_ind,acceleration_fft_stop_ind]-1)./sfreq,"r-");
        hold off
        grid on
    end

    fft_vorticity = fft(vorticity_data(:,acceleration_fft_start_ind:acceleration_fft_stop_ind),[],2);
    fft_vorticity = abs(fft_vorticity(:,2:end-1) ./ time_step);
    fft_frequency = (1:time_step-2)./time_step.*sfreq;

    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(fft_frequency,fft_vorticity(ii,:));
        %semilogx(fft_frequency,fft_constVelocity(ii,:));
        grid on
        xlim(x_range);
        title(sprintf("FFT vorticity (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    end
    sgtitle("FFT result of vorticity in constant velocity part");
    saveas(gcf,output_folder+"/"+"FFT_v_constVelocity.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_constVelocity.png"));

    % Compute Strouhal number
     logging_func(sprintf("Compute Strouhal number (acceleration vorticity)\n"));
     U = mean(characteristicVelocity_acceleration(acceleration_fft_start_ind:acceleration_fft_stop_ind),"all");
     for ii = 1:num_search_pixel
        fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
        [pks,locs] = findpeaks(fft_vorticity(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
            fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
        [~,index] = maxk(pks,strouhal_freqs);
        strouhal_number = locs(index) .* (object_width ./ U) ./ 1000;
        for jj = 1:strouhal_freqs
            fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n",jj,strouhal_number(jj),locs(index(jj)),U);
        end
     end
end

% constant velocity
if (~isnan(velocity_u_fluctuation_constVelocity))
    logging_func("Analyze Vorticity (constant Velocity)");
    time_step = length(meanMap_vorticity(constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
    if (~isnan(constVelocity_timeStep))
        figure
        tiledlayout('flow')
        vorticity_data = zeros(num_search_pixel,length(vorticity_constVelocity));
        for ii = 1:num_search_pixel
            ind_x = zeros(1,length(vorticity_constVelocity));
            ind_y = zeros(1,length(vorticity_constVelocity));
            % Search index
            for jj = 1:length(vorticity_constVelocity)
                [~,ind_x(jj)] = min(abs(mean_x{acceleration_timeStep+jj} - search_x(ii)));
                [~,ind_y(jj)] = min(abs(mean_y{acceleration_timeStep+jj} - search_y(ii)));
                vorticity_data(ii,jj) = vorticity_constVelocity{jj}(ind_y(jj),ind_x(jj));
            end
    
            nexttile
            hold on
            plot((0:(length(vorticity_data)-1))./sfreq,vorticity_data(ii,:));
            xline(([constVelocity_fft_start_ind,constVelocity_fft_stop_ind]-1)./sfreq,"r-");
            hold off
            grid on
        end
    
        fft_vorticity = fft(vorticity_data(:,constVelocity_fft_start_ind:constVelocity_fft_stop_ind),[],2);
        fft_vorticity = abs(fft_vorticity(:,2:end-1) ./ time_step);
        fft_frequency = (1:time_step-2)./time_step.*sfreq;
    
        figure
        tiledlayout('flow')
        for ii = 1:num_search_pixel
            nexttile
            plot(fft_frequency,fft_vorticity(ii,:));
            %semilogx(fft_frequency,fft_constVelocity(ii,:));
            grid on
            xlim(x_range);
            title(sprintf("FFT vorticity (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        end
        sgtitle("FFT result of vorticity in constant velocity part");
        saveas(gcf,output_folder+"/"+"FFT_v_constVelocity.png");
        logging_func(sprintf("Save %s",output_folder+"/"+"FFT_v_constVelocity.png"));
    
        % Compute Strouhal number
         logging_func(sprintf("Compute Strouhal number (constant velocity vorticity)\n"));
         for ii = 1:num_search_pixel
            fprintf(" Position %d (%.3f,%.3f)\n",ii,real_x(ii,1),real_y(ii,1));
            [pks,locs] = findpeaks(fft_vorticity(ii,x_range(1) <= fft_frequency & x_range(2) >= fft_frequency), ...
                fft_frequency(x_range(1) <= fft_frequency & x_range(2) >= fft_frequency));
            [~,index] = maxk(pks,strouhal_freqs);
            strouhal_number = locs(index) .* (object_width ./ reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1)) ./ 1000;
            for jj = 1:strouhal_freqs
                fprintf("  Strouhal %d : %.3f (f = %.3f [Hz], U = %.3f [m/s])\n", ...
                    jj,strouhal_number(jj),locs(index(jj)),reference_velocity_ratio(end) .* truck_constVelocity_mean_vx(1));
            end
         end
    end
end

%% Autocorrelation
logging_func("Compute autocorrelation");
acceleration_timeStep = length(velocity_u_fluctuation_acceleration);
meanMap_u_filtered_acceleration = meanMap_u_filtered(1:acceleration_timeStep);
meanMap_v_filtered_acceleration = meanMap_v_filtered(1:acceleration_timeStep);
%meanMap_absVelocity_filtered_acceleration = sqrt(meanMap_u_filtered_acceleration{:}.^2 + meanMap_v_filtered_acceleration{:}.^2);
if (~isnan(velocity_v_fluctuation_constVelocity))
    meanMap_u_filtered_constVelocity = meanMap_u_filtered(acceleration_timeStep+1:end);
    meanMap_v_filtered_constVelocity = meanMap_v_filtered(acceleration_timeStep+1:end);
    %meanMap_absVelocity_filtered_constVelocity = sqrt(meanMap_u_filtered_constVelocity{:}.^2 + meanMap_v_filtered_constVelocity{:}.^2);
    constVelocity_timeStep = length(meanMap_u_filtered_constVelocity);
end

% Acceleration
% u
logging_func("Compute autocorrelation of u (acceleration)");
u_acceleration = zeros(num_search_pixel,acceleration_timeStep);
for ii = 1:num_search_pixel
    for jj = 1:acceleration_timeStep
        [~,ind_x(jj)] = min(abs(mean_x{jj} - search_x(ii)));
        [~,ind_y(jj)] = min(abs(mean_y{jj} - search_y(ii)));
        u_acceleration(ii,jj) = meanMap_u_filtered_acceleration{jj}(ind_y(jj),ind_x(jj)) - characteristicVelocity_acceleration(jj);
    end
    u_acceleration(ii,:) =  u_acceleration(ii,:) - mean(u_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    [autoCorrelation_u_acceleration(:,ii),autoCorrelation_u_lag] = ...
        xcorr(u_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot((0:(length(u_acceleration(ii,:))-1))./sfreq,u_acceleration(ii,:));
    hold off
    title(sprintf("u for Autocorrelation (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot(autoCorrelation_u_lag./sfreq,autoCorrelation_u_acceleration(:,ii));
    hold off
    title(sprintf("Autocorrelation u (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

% v
logging_func("Compute autocorrelation of v (acceleration)");
v_acceleration = zeros(num_search_pixel,acceleration_timeStep);
for ii = 1:num_search_pixel
    for jj = 1:acceleration_timeStep
        [~,ind_x(jj)] = min(abs(mean_x{jj} - search_x(ii)));
        [~,ind_y(jj)] = min(abs(mean_y{jj} - search_y(ii)));
        v_acceleration(ii,jj) = meanMap_v_filtered_acceleration{jj}(ind_y(jj),ind_x(jj)) - characteristicVelocity_acceleration(jj);
    end
    v_acceleration(ii,:) =  v_acceleration(ii,:) - mean(v_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    [autoCorrelation_v_acceleration(:,ii),autoCorrelation_v_lag] = ...
        xcorr(v_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot((0:(length(v_acceleration(ii,:))-1))./sfreq,v_acceleration(ii,:));
    hold off
    title(sprintf("v for Autocorrelation (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot(autoCorrelation_v_lag./sfreq,autoCorrelation_v_acceleration(:,ii));
    hold off
    title(sprintf("Autocorrelation v (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

% absVelocity
logging_func("Compute autocorrelation of absolute velocity (acceleration)");
absVelocity_acceleration = zeros(num_search_pixel,acceleration_timeStep);
for ii = 1:num_search_pixel
    for jj = 1:acceleration_timeStep
        [~,ind_x(jj)] = min(abs(mean_x{jj} - search_x(ii)));
        [~,ind_y(jj)] = min(abs(mean_y{jj} - search_y(ii)));
        absVelocity_acceleration(ii,jj) = sqrt(meanMap_u_filtered_acceleration{jj}(ind_y(jj),ind_x(jj)).^2 + meanMap_v_filtered_acceleration{jj}(ind_y(jj),ind_x(jj)).^2) - characteristicVelocity_acceleration(jj);
    end
    absVelocity_acceleration(ii,:) =  absVelocity_acceleration(ii,:) - mean(absVelocity_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind));
    [autoCorrelation_absVelocity_acceleration(:,ii),autoCorrelation_absVelocity_lag] = ...
        xcorr(absVelocity_acceleration(ii,acceleration_fft_start_ind:acceleration_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot((0:(length(absVelocity_acceleration(ii,:))-1))./sfreq,absVelocity_acceleration(ii,:));
    hold off
    title(sprintf("v for Autocorrelation (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

figure
tiledlayout('flow')
for ii = 1:num_search_pixel
    nexttile
    hold on
    plot(autoCorrelation_absVelocity_lag./sfreq,autoCorrelation_absVelocity_acceleration(:,ii));
    hold off
    title(sprintf("Autocorrelation v (acceleration, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
    grid on
end

% Constant Velocity
% u
if (~isnan(velocity_u_fluctuation_constVelocity))
    logging_func("Compute autocorrelation of u (constant velocity)");
    u_constVelocity = zeros(num_search_pixel,constVelocity_timeStep);
    for ii = 1:num_search_pixel
        for jj = 1:constVelocity_timeStep
            [~,ind_x(jj)] = min(abs(mean_x{acceleration_timeStep+jj} - search_x(ii)));
            [~,ind_y(jj)] = min(abs(mean_y{acceleration_timeStep+jj} - search_y(ii)));
            u_constVelocity(ii,jj) = meanMap_u_filtered_constVelocity{jj}(ind_y(jj),ind_x(jj));
        end
        u_constVelocity(ii,:) =  u_constVelocity(ii,:) - mean(u_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
        %autoCorr(u_constVelocity(ii,:),sfreq);
        [autoCorrelation_u_constVelocity(:,ii),autoCorrelation_u_lag] = xcorr(u_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        hold on
        plot((0:(constVelocity_timeStep-1))./sfreq,u_constVelocity(ii,:));
        xline(([constVelocity_fft_start_ind,constVelocity_fft_stop_ind]-1)./sfreq,"r-");
        hold off
        title(sprintf("u (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(autoCorrelation_u_lag./sfreq,autoCorrelation_u_constVelocity(:,ii));
        title(sprintf("Autocorrelation u (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
end

% v
if (~isnan(velocity_v_fluctuation_constVelocity))
    logging_func("Compute autocorrelation of u (constant velocity)");
    v_constVelocity = zeros(num_search_pixel,constVelocity_timeStep);
    for ii = 1:num_search_pixel
        for jj = 1:constVelocity_timeStep
            [~,ind_x(jj)] = min(abs(mean_x{acceleration_timeStep+jj} - search_x(ii)));
            [~,ind_y(jj)] = min(abs(mean_y{acceleration_timeStep+jj} - search_y(ii)));
            v_constVelocity(ii,jj) = meanMap_v_filtered_constVelocity{jj}(ind_y(jj),ind_x(jj));
        end
        v_constVelocity(ii,:) =  v_constVelocity(ii,:) - mean(v_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
        %autoCorr(u_constVelocity(ii,:),sfreq);
        [autoCorrelation_v_constVelocity(:,ii),autoCorrelation_v_lag] = xcorr(v_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        hold on
        plot((0:(constVelocity_timeStep-1))./sfreq,v_constVelocity(ii,:));
        xline(([constVelocity_fft_start_ind,constVelocity_fft_stop_ind]-1)./sfreq,"r-");
        hold off
        title(sprintf("v (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(autoCorrelation_v_lag./sfreq,autoCorrelation_v_constVelocity(:,ii));
        title(sprintf("Autocorrelation v (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
end

% absVelocity
if (~isnan(velocity_u_fluctuation_constVelocity) & ~isnan(velocity_v_fluctuation_constVelocity))
    logging_func("Compute autocorrelation of u (constant velocity)");
    absVelocity_constVelocity = zeros(num_search_pixel,constVelocity_timeStep);
    for ii = 1:num_search_pixel
        for jj = 1:constVelocity_timeStep
            [~,ind_x(jj)] = min(abs(mean_x{acceleration_timeStep+jj} - search_x(ii)));
            [~,ind_y(jj)] = min(abs(mean_y{acceleration_timeStep+jj} - search_y(ii)));
            absVelocity_constVelocity(ii,jj) = sqrt(meanMap_v_filtered_constVelocity{jj}(ind_y(jj),ind_x(jj)).^2 + meanMap_v_filtered_constVelocity{jj}(ind_y(jj),ind_x(jj)).^2);
        end
        absVelocity_constVelocity(ii,:) =  absVelocity_constVelocity(ii,:) - mean(absVelocity_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind));
        %autoCorr(u_constVelocity(ii,:),sfreq);
        [autoCorrelation_absVelocity_constVelocity(:,ii),autoCorrelation_absVelocity_lag] = xcorr(absVelocity_constVelocity(ii,constVelocity_fft_start_ind:constVelocity_fft_stop_ind)',autocoorelation_lag*sfreq,'coeff'); 
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        hold on
        plot((0:(constVelocity_timeStep-1))./sfreq,absVelocity_constVelocity(ii,:));
        xline(([constVelocity_fft_start_ind,constVelocity_fft_stop_ind]-1)./sfreq,"r-");
        hold off
        title(sprintf("Absolute velocity (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
    
    figure
    tiledlayout('flow')
    for ii = 1:num_search_pixel
        nexttile
        plot(autoCorrelation_absVelocity_lag./sfreq,autoCorrelation_absVelocity_constVelocity(:,ii));
        title(sprintf("Autocorrelation absolute velosity (constant velocity, x:%.3f [mm], y:%.3f [mm])",real_x(ii,1),real_y(ii,1)));
        grid on
    end
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

%{
function autoCorr(data,Fs)
M = length(data);
N = M / 2;
dt = 1 / Fs;

v = data(:,2);
v_norm = v - mean(v);

tau = 6 * 10^-2;    % 遅れ時間

[autocor, lags] = xcorr(v_norm, tau*Fs, 'coeff');

[pks,locs] = findpeaks(autocor);

x = lags/Fs;

% plot(x,autocor)   % 自己相関関数の実の出力の場合
plot(x, autocor,x(locs),pks,'or')   % ピークの算出を行う場合
xlabel('lag(sec)')
ylabel('Autocorrelation')

cycles = diff(x(locs));
meanCycle = mean(cycles);

f = 1 / meanCycle;
end
%}
