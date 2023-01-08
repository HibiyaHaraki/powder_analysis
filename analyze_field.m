clear; close all; clc;
logging_func("Analyze field data");
%% Explanation

%
% analyze_field
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%

%% Setting

% Experiment PIV data file name
experiment_PIV_data_filename = '../PIV_xo350W60固定壁L1H1/result/PIV_data.mat';

% Truck data file name
experiment_Truck_data_filename = "../加速度データ_xo350W60固定壁L1H1/truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W60固定壁L1H1/result';

% Search pixel
search_x = [405, 435, 465, 495, 405, 435, 465, 495]; % [mm]
search_y = [10, 10, 10, 10, 50, 50, 50, 50]; % [mm]

% Integration time
int_start = 0; % [s]
int_stop  = 3; % [s]

%% Check Inputs

if (~exist(experiment_PIV_data_filename,'file'))
    error("Cannot find %s!!",experiment_PIV_data_filename);
end

if (~exist(experiment_Truck_data_filename,'file'))
    error("Cannot find %s!!",experiment_Truck_data_filename);
end

if (length(search_x)  ~= length(search_y))
    error("Please input same number of the position in x and y");
end
numSearchPosition = length(search_x);

%% Load data
experiment_PIV_data = load(experiment_PIV_data_filename, ...
    "meanMap_u_filtered", ...
    "mean_x", "mean_y", "object_width");

Truck_data = load(experiment_Truck_data_filename, ...
    'truck_acceleration_mean_vx', 'sfreq');

%% Convert unit of search pixel
load(experiment_PIV_data_filename,"frontDistance");
search_x = search_x - frontDistance;

search_x = search_x ./ 1000;
search_y = search_y ./ 1000;

%% Compute characteristic velocity
logging_func("Compute characteristic velocity");
characteristicVelocity = compute_characteristicVelocity_field(experiment_PIV_data_filename,experiment_Truck_data_filename);

%% Compare experiment and reference PIV data
time_step = length(characteristicVelocity);
compare_result = zeros(numSearchPosition,time_step);
for ii = 1:numSearchPosition
    for jj = 1:time_step
        experiment_x = experiment_PIV_data.mean_x{jj};
        experiment_y = experiment_PIV_data.mean_y{jj};
        [~,search_x_ind] = min(abs(experiment_x - search_x(ii)),[],'all');
        [~,search_y_ind] = min(abs(experiment_y - search_y(ii)),[],'all');
        experiment_data = experiment_PIV_data.meanMap_u_filtered{jj}(search_y_ind,search_x_ind);
        reference_data  = characteristicVelocity{jj}(search_y_ind,search_x_ind);
        compare_result(ii,jj) = compare_datas(reference_data,experiment_data);
    end
end

%% Compute integration
int_start_ind = round(int_start*Truck_data.sfreq)+1;
int_stop_ind  = round(int_stop*Truck_data.sfreq)+1;

if (int_start_ind <= 0)
    int_start_ind = 1;
end

if (int_stop_ind > time_step)
    int_stop_ind = time_step;
    warning("Max time is %.3f [s]",time_step/Truck_data.sfreq);
end

real_int_start = int_start_ind ./ Truck_data.sfreq;
real_int_stop = int_stop_ind ./ Truck_data.sfreq;

for ii = 1:numSearchPosition
    fprintf("Data %d (x=%.3f, y=%.3f)\n",ii,search_x(ii),search_y(ii));
    fprintf(" Time step : %d\n",time_step);
    fprintf(" Integration : %.3f\n",sum(compare_result(ii,int_start_ind:int_stop_ind)));
end

%% Visualize

for ii = 1:numSearchPosition
    figure
    hold on
    plot((0:(time_step-1))./Truck_data.sfreq,compare_result(ii,:));
    xline([real_int_start,real_int_stop],"r-");
    hold off
    grid on
    xlabel("Time [s]");
    ylabel("Difference");
    xlim([0,(time_step-1)./Truck_data.sfreq]);
    title(sprintf("Compare Reference and Experiment data (x=%.3f, y=%.3f)",search_x(ii),search_y(ii)));
end

%% Function for comparing data
function answer = compare_datas(reference,experiment)
    answer = abs(experiment-reference);  % 出力したい数式
end
