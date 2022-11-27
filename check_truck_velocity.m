clear; close all; clc;
tic
logging_func("Check truck velocity");
%% Explanation

%
% check_truck_velocity
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

% Specify input files (truck acceleration data)
input_file_name = "../加速度データ_xo350W45固定壁L1H1/maindata";
input_file_ID = [1:5];

% Specify max number of data
max_num_data = 2*10^4; % 最大20s計測で1000fpsだから、データの最大数は2*10^4

% Point of starting acceleration in PIV data
PIV_start_accleration_point = [590, 591, 680, 632, 560];

% Specify velocity analysis mode (1-analysis 0-not analsis)
ACCELERATION_MODE  = 1;
CONSTVELOCITY_MODE = 0;

% Specify threshold for accelerating start
accelerating_start_threshold = 3E-3;

% Specify threshold for constant velocity stop
constVelocity_stop_threshold = 0.7;
offset = 2;

% Save mat file (1-save, 0-unsave)
SAVE_MODE = 1;

% Specify output folder and file
output_folder_name = '../加速度データ_xo350W45固定壁L1H1';
output_file_name = 'truck_data';

%% Check specified files and folders

% Create input file name string
num_input_files = length(input_file_ID);
input_filenames = string(zeros(1,num_input_files));

% Check input files
for ii = 1:num_input_files
    input_filenames(ii) = sprintf("%s%04d.CSV",input_file_name,input_file_ID(ii));
    if (~exist(output_folder_name,'file'))
        error("Error: Cannot find %s !!",input_filenames(ii));
        return;
    end
end

% Check output folder
if (~exist(output_folder_name,'dir'))
    mkdir(output_folder_name);
    logging_func(sprintf("Created output folder : %s",output_folder_name));
end

% Check PIV starting point
if (num_input_files ~= length(PIV_start_accleration_point))
    error("Error: number of PIV start data should be %d !!",num_input_files);
    return;
end

% Check Analysis Mode
if (~ACCELERATION_MODE && ~CONSTVELOCITY_MODE)
    error("Error: Choose analysis mode!!");
    return;
end

%% Create Matrixes
timesteps = zeros(1,num_input_files); % Number of ime step
sf = zeros(1,num_input_files); % Sampling Frequency

t  = nan(max_num_data,num_input_files); % Time
ax = nan(max_num_data,num_input_files); % Truck x-acceleration
ay = nan(max_num_data,num_input_files); % Truck y-acceleration
az = nan(max_num_data,num_input_files); % Truck z-acceleration
%vx = nan(max_num_data,num_input_files); % Truck x-velocity
%vy = nan(max_num_data,num_input_files); % Truck y-velocity
%vz = nan(max_num_data,num_input_files); % Truck z-velocity
%x  = nan(max_num_data,num_input_files); % Truck x-position
%y  = nan(max_num_data,num_input_files); % Truck y-position
%z  = nan(max_num_data,num_input_files); % Truck z-position

start_acceleration  = nan(1,num_input_files); % Point of starting acceleration
stop_acceleration   = nan(1,num_input_files); % Point of stopping acceleration
start_constVelocity = nan(1,num_input_files); % Point of starting constant velocity move
stop_constVelocity  = nan(1,num_input_files); % Point of stopping constant velocity move

%% Import truck acceleration data
for ii = 1:num_input_files
    logging_func(sprintf("Load %s",input_filenames(ii)));

    % Import .csv data
    csv_data = readcell(input_filenames(ii));

    % Get Sampling Frequency
    sf(ii) = cell2mat(csv_data(6,2));

    % Get number of time step
    timesteps(ii) = length(csv_data(15:end,1)); % CSVファイルの1列目の15行目から最後までのデータを読み込んで格納、15行目から時間データが格納されている

    % Get time data
    t(1:timesteps(ii),ii) = cell2mat(csv_data(15:end,1));

    % Get Acceleration data
    ax(1:timesteps(ii),ii) = cell2mat(csv_data(15:end,2));
    ay(1:timesteps(ii),ii) = cell2mat(csv_data(15:end,3));
    az(1:timesteps(ii),ii) = cell2mat(csv_data(15:end,4));
end

if (length(unique(sf)) ~= 1)
    error("Error: This script can only treat same sampling frequency data!!");
    return;
end
sfreq = sf(1);

%% Compute velocity and position by integration

% Compute Velocity
vx = cumtrapz(t(:,1),ax,1);
vy = cumtrapz(t(:,1),ay,1);
vz = cumtrapz(t(:,1),az,1);

% Compute position
x = cumtrapz(t(:,1),vx,1);
y = cumtrapz(t(:,1),vy,1);
z = cumtrapz(t(:,1),vz,1);

%% Specify start and stop point of accelration and constant velocity
for ii = 1:num_input_files
    % Specify stop point of acceleration
    [~,stop_acceleration(ii)] = max(vx(:,ii));

    % Specify start point of acceleration
    check_variable = 1;
    for jj = 1:stop_acceleration(ii)-1
        if (vx(stop_acceleration(ii)-jj,ii) < accelerating_start_threshold)
            check_variable = 0;
            start_acceleration(ii) = stop_acceleration(ii)-jj;
            break;  % 見つかった時点で繰り返しを強制終了
        end
    end

    if (check_variable)
        error("Error!! Cannot find start point of acceleration in %s.",input_filenames(ii)); %閾値を全てのデータが超えていたらエラーを吐く
        return;
    end

    % Specify constant velocity part
    if (CONSTVELOCITY_MODE)
        % Specify start point of constant velocity
        start_constVelocity(ii) = stop_acceleration(ii) + 1;
    
        % Specify stop point of constant velocity
        check_variable = 1;
        [pks,locs] = findpeaks(vx(:,ii));
        for jj = 1:length(locs)
            if (pks(end-jj+1) > constVelocity_stop_threshold)
                stop_constVelocity(ii) = locs(end-jj+1-offset);
                check_variable = 0;
                break
            end
        end

        if (check_variable)
            error("Error!! Cannot find stop point of constant velocity in %s.",input_filenames(ii)); %閾値を全てのデータが超えていたらエラーを吐く
            return;
        end

        if (stop_constVelocity(ii) < start_constVelocity(ii))
            error("Error!! Cannot specify the constant velocity part in %s.",input_filenames(ii));
            return;
        end
    end
end

%% Compute mean data

truck_acceleration_time_step  = NaN;
truck_constVelocity_time_step = NaN;

truck_acceleration_mean_vx  = NaN;
truck_constVelocity_mean_vx = NaN;

% Compute mean vx in acceleration
if (ACCELERATION_MODE)
    % Get minimum time step in acceleration
    truck_acceleration_time_step = min(stop_acceleration - start_acceleration + 1);
    
    % Compute truck mean x-velocity
    truck_acceleration_vx = zeros(truck_acceleration_time_step,num_input_files);
    for ii = 1:num_input_files
        truck_acceleration_vx(:,ii) = vx(start_acceleration(ii):start_acceleration(ii)+truck_acceleration_time_step-1,ii);
    end
    truck_acceleration_mean_vx = mean(truck_acceleration_vx,2);
    logging_func("Compute truck mean vx in acceleration");
end

% Compute mean vx in acceleration
if (CONSTVELOCITY_MODE)
    % Get minimum time step in acceleration
    truck_constVelocity_time_step = min(stop_constVelocity - start_constVelocity + 1);
    
    % Compute truck mean x-velocity
    truck_constVelocity_vx = zeros(truck_constVelocity_time_step,num_input_files);
    for ii = 1:num_input_files
        truck_constVelocity_vx(:,ii) = vx(start_constVelocity(ii):start_constVelocity(ii)+truck_constVelocity_time_step-1,ii);
    end
    truck_constVelocity_mean_vx = mean(truck_constVelocity_vx,2);
    logging_func("Compute truck mean vx in constant velocity move");
end

% Create all mean data
truck_all_mean_vx = [];
if (ACCELERATION_MODE)
    truck_all_mean_vx = [truck_all_mean_vx;truck_acceleration_mean_vx];
end

if (CONSTVELOCITY_MODE)
    truck_all_mean_vx = [truck_all_mean_vx;truck_constVelocity_mean_vx];
end


%% Compute PIV data order
PIV_truck_offset = zeros(1,num_input_files);

%% Determine Pixel
% Compute pixel size (determinePixel.m)
[neccesary_pixel_size,determined_pixel_size,convert_point] = determinePixel(sfreq,truck_all_mean_vx);

% Visualize
figure
hold on
plot((0:length(truck_all_mean_vx)-1)./sfreq,neccesary_pixel_size,'LineWidth',2.0);
stairs((0:length(truck_all_mean_vx)-1)./sfreq,determined_pixel_size,'LineWidth',2.5);
for ii = 1:length(convert_point)
    xline(convert_point(ii)/sfreq,'k-', ...
        {sprintf("Size: %d ,Time:%.3f [s] (%d)", ...
        determined_pixel_size(convert_point(ii)), ...
        convert_point(ii)/sfreq, ...
        convert_point(ii))}, ...
        'LineWidth',1.1);
end
hold off
grid on
legend("Necessary Pixel size","Determined PIVlab pixel size","Location","best");
xlabel("t [s]");
ylabel("Pixel size");
xlim([0,(length(determined_pixel_size)-1)./sfreq]);
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"determied_pixel_size.png");
end

%% Show information
% Show pixel size for each experiment
fprintf("\n");
for ii = 1:num_input_files
    % Show PIV cut information
    fprintf("Ex %d) \n",input_file_ID(ii));
    fprintf(" PIV acceleration start   : %4d\n",start_acceleration(ii));
    fprintf(" Truck Acceleration start : %4d\n",PIV_start_accleration_point(ii));
    PIV_truck_offset(ii) = PIV_start_accleration_point(ii) - start_acceleration(ii);
    fprintf(" Offset : %d\n",PIV_truck_offset(ii));
    fprintf(" Acceleration start     : (Truck,PIV) = (%4d,%4d)\n", start_acceleration(ii), start_acceleration(ii) + PIV_truck_offset(ii));
    fprintf(" Acceleration stop      : (Truck,PIV) = (%4d,%4d)\n", start_acceleration(ii) + truck_acceleration_time_step - 1, start_acceleration(ii) + PIV_truck_offset(ii) +  truck_acceleration_time_step - 1);
    if (CONSTVELOCITY_MODE)
        fprintf(" Const velocity start   : (Truck,PIV) = (%4d,%4d)\n",start_constVelocity(ii),start_constVelocity(ii) + PIV_truck_offset(ii));
        fprintf(" Const velocity start   : (Truck,PIV) = (%4d,%4d)\n",start_constVelocity(ii) + truck_constVelocity_time_step - 1,start_constVelocity(ii) + PIV_truck_offset(ii) + truck_constVelocity_time_step - 1);
    end

    % Show pixel size data
    truck_all_timestep = [];
    if (ACCELERATION_MODE)
        truck_all_timestep = [truck_all_timestep,start_acceleration(ii):start_acceleration(ii) + truck_acceleration_time_step - 1];
    end
    
    if (CONSTVELOCITY_MODE)
        truck_all_timestep = [truck_all_timestep,start_constVelocity(ii):start_constVelocity(ii) + truck_constVelocity_time_step - 1];
    end

    PIV_all_timestep = truck_all_timestep + PIV_truck_offset(ii);

    for jj = 1:length(convert_point)-1
        fprintf(" Pixel size %2d : (Truck,PIV) = (%4d~%4d,%4d~%4d) \n", ...
            determined_pixel_size(convert_point(jj)), ...
            convert_point(jj),convert_point(jj+1)-1, ...
            PIV_all_timestep(convert_point(jj)),PIV_all_timestep(convert_point(jj+1)-1));
    end
    fprintf(" Pixel size %2d : (Truck,PIV) = (%4d~%4d,%4d~%4d) \n", ...
            determined_pixel_size(convert_point(end)), ...
            convert_point(end),length(truck_all_mean_vx), ...
            PIV_all_timestep(convert_point(end)),PIV_all_timestep(end));
    fprintf("\n");
end
for ii = 1:num_input_files
    truck_all_t = [];
    if (ACCELERATION_MODE)
        truck_all_t = [truck_all_t,start_acceleration(ii):start_acceleration(ii) + truck_acceleration_time_step - 1];
    end

    if (CONSTVELOCITY_MODE)
        truck_all_t = [truck_all_t,start_constVelocity(ii):start_constVelocity(ii) + truck_constVelocity_time_step - 1];
    end


end

%% Output data
if (SAVE_MODE)
    % Output following data
    %  - sampling frequency
    %  - number of time step in acceleration
    %  - number of time step in constant velocity move
    %  - mean x-velocity in acceleration
    %  - mean x-velocity in constant velocity move

    save(output_folder_name+"/"+output_file_name, ...
        'sfreq', ...
        'determined_pixel_size','convert_point', ...
        'truck_acceleration_time_step','truck_constVelocity_time_step', ...
        'truck_acceleration_mean_vx','truck_constVelocity_mean_vx');
    logging_func(sprintf("Output %s",output_folder_name+"/"+output_file_name));
end

%% Visualize
% Visualize each experiment acceleration data
%{
% a_x
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),ax(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('a_x [m/s^2]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("a_x");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_ax.png");
end

% a_y
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),ay(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('a_y [m/s^2]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("a_y");if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_ay.png");
end

% a_z
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),az(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('a_z [m/s^2]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("a_z");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_az.png");
end
%}

% Visualize each experiment velocity data
% v_x
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),vx(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('v_x [m/s]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("v_x");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_vx.png");
end

%{
% v_y
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),vy(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('v_y [m/s]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("v_y");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_vy.png");
end

% v_z
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),vz(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('v_x [m/s]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("v_z");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_vz.png");
end
%}

% Visualize each experiment position data
%{
% x
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),x(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('x [m]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("x");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_x.png");
end

% y
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),y(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('y [m]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("y");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_y.png");
end

% z
figure
tiledlayout('flow') % tiledlayoutは1枚のfigureにグラフを複数作るときに使う、flowは配置を画面サイズに応じて調整するという意味
for ii = 1:num_input_files
    nexttile % 新しいグラフを作る
    hold on
    plot(t(1:timesteps(ii),ii),z(1:timesteps(ii),ii),"k-");
    xline(t(start_acceleration(ii)),"b-");
    xline(t( stop_acceleration(ii)),"b-");
    if (CONSTVELOCITY_MODE)
        xline(t(start_constVelocity(ii)),"r-");
        xline(t( stop_constVelocity(ii)),"r-");
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('z [m]');
    title(sprintf("Experiment %d",input_file_ID(ii)));
end
sgtitle("z");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_z.png");
end
%}

%% Visualize all data during accelerating in 1 figure
if (ACCELERATION_MODE)
    visualize_start_point = start_acceleration;
else
    visualize_start_point = start_constVelocity;
end

if (CONSTVELOCITY_MODE)
    visualize_stop_point = stop_constVelocity;
else
    visualize_stop_point = stop_acceleration;
end

visualize_timestep = visualize_stop_point - visualize_start_point + 1;

% Create figure
figure
legend_str = []; % 判例を入れる変数
hold on
for ii = 1:num_input_files
    plot(t(1:visualize_timestep(ii),1),vx(visualize_start_point(ii):visualize_stop_point(ii),ii));
    legend_str = [legend_str,sprintf("Experiment %d",input_file_ID(ii))];
end
hold off
grid on
legend(legend_str,'Location','southoutside');
ylabel('Velocity x [m/s]');
xlabel('t [s]');
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_all_vx.png");
end

%% Visualize truck mean data
figure
plot((0:length(truck_all_mean_vx)-1)./sfreq,truck_all_mean_vx);
hold on
if (ACCELERATION_MODE && CONSTVELOCITY_MODE)
    xline(truck_acceleration_time_step./sfreq,"k-")
end
hold off
grid on
xlim([0 (length(truck_all_mean_vx)-1)./sfreq]);
ylabel('Velocity x [m/s]');
xlabel('t [s]');
title("Truck mean x-velocity");
if (SAVE_MODE)
    saveas(gcf,output_folder_name+"/"+"truck_mean_vx.png");
end
toc

%% Function for log and time-stamp
function logging_func(message_str)
str_datetime = string(datetime('now','Format','MM/dd HH:mm:ss.SSS'));
fprintf("[%s] %s\n",str_datetime,message_str);
end
