clear; close all; clc;
logging_func("Analyze a map data");
%% Explanation

%
% analyze_a_map
%
% created by Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Create mean map and mean value from PIV data
% 
% Warning
%  This script needs following files in the same folder.
%   * 
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W60移動壁L1H1/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W60移動壁L1H1/Truck_data.mat";

% Output folder
output_folder = '../PIV_xo350W60移動壁L1H1/result';

% Search time
search_t = [2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2];

% Search reynolds
search_reynolds = [];

% Quiver size
quiverSize = 3;

% Vorticity existance
VORTICITY = true;

% Constant
nu = 1.5*10^(-5);

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

if (isempty(search_t) && isempty(search_reynolds))
    error("Please input search conditions!!");
    return;
end

if (~isempty(search_t))
    if (~isvector(search_t))
        error("search_t should be vector!!");
        return;
    end

    if (length(search_t(1,:)) > 1)
        search_t = search_t';
    end
end

if (~isempty(search_reynolds))
    if (~isvector(search_reynolds))
        error("search_reynolds should be vector!!");
        return;
    end
    
    if (length(search_reynolds(1,:)) > 1)
        search_reynolds = search_reynolds';
    end
end

num_search_t = length(search_t);

%% Import data
if (VORTICITY)
    load(PIV_data_filename, ...
        "width","object_width","object_length", ...
        "mean_x","mean_y", ...
        "meanMap_u_filtered","meanMap_v_filtered",...
        "meanMap_vorticity");
else
    load(PIV_data_filename, ...
     "width","object_width","object_length",  ...
     "mean_x","mean_y", ...
     "meanMap_u_filtered","meanMap_v_filtered");
end

load(Truck_data_filename,"sfreq");
time_step = length(meanMap_u_filtered);

%% Get index of search time
ind_t_t = [];
if (~isempty(search_t))
    [~,ind_t_t] = min(abs((0:time_step-1)./sfreq - search_t),[],2);
end

%% Get index of search reynolds
ind_t_reynolds = [];
[Re_D,Re_H,Re_L] = solve_reynolds(PIV_data_filename,Truck_data_filename);
if (~isempty(search_reynolds))
    [~,ind_t_reynolds] = min(abs(Re_H' - search_reynolds),[],2);
end

%% Time step data
ind_t = sort(unique([ind_t_t;ind_t_reynolds]));
real_t = (ind_t-1)./sfreq;
num_search_t = length(ind_t);

for ii = 1:num_search_t
    fprintf("t=%.3f\n",real_t(ii));
    fprintf("  Re_D = %.3f\n",Re_D(ind_t(ii)));
    fprintf("  Re_H = %.3f\n",Re_H(ind_t(ii)));
    fprintf("  Re_L = %.3f\n",Re_L(ind_t(ii)));
end

%% Search Max and Min value
max_u = 0; max_v = 0;
min_u = 0; min_v = 0;
max_magnitude = 0;
min_magnitude = 0;
max_vorticity = 0;
min_vorticity = 0;
for ii = 1:length(meanMap_u_filtered)
    if (max(meanMap_u_filtered{ii},[],"all") > max_u)
        max_u = max(meanMap_u_filtered{ii},[],"all");
    end

    if (min(meanMap_u_filtered{ii},[],"all") < min_u)
        min_u = min(meanMap_u_filtered{ii},[],"all");
    end

    if (max(meanMap_v_filtered{ii},[],"all") > max_v)
        max_v = max(meanMap_v_filtered{ii},[],"all");
    end

    if (min(meanMap_v_filtered{ii},[],"all") < min_v)
        min_v = min(meanMap_v_filtered{ii},[],"all");
    end

    if (max(sqrt(meanMap_v_filtered{ii}.^2 + meanMap_u_filtered{ii}.^2),[],"all") > max_magnitude)
        max_magnitude = max(sqrt(meanMap_v_filtered{ii}.^2 + meanMap_u_filtered{ii}.^2),[],"all");
    end

    if (min(sqrt(meanMap_v_filtered{ii}.^2 + meanMap_u_filtered{ii}.^2),[],"all") < min_magnitude)
        min_magnitude = min(sqrt(meanMap_v_filtered{ii}.^2 + meanMap_u_filtered{ii}.^2),[],"all");
    end

    if (VORTICITY)
        if (max(meanMap_vorticity{ii},[],"all") > max_vorticity)
            max_vorticity = max(meanMap_vorticity{ii},[],"all");
        end
    
        if (min(meanMap_vorticity{ii},[],"all") < min_vorticity)
            min_vorticity = min(meanMap_vorticity{ii},[],"all");
        end
    end
end

%% Create color map
logging_func("Create color map");
for ii = 1:num_search_t
    % u
    figure
    pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,meanMap_u_filtered{ind_t(ii)});
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    caxis([min_u max_u]);
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title(sprintf("Velocity u (t=%.3f[s])",real_t(ii)));

    % v
    figure
    pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,meanMap_v_filtered{ind_t(ii)});
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    caxis([min_v max_v]);
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title(sprintf("Velocity v (t=%.3f[s])",real_t(ii)));

    % magnitude
    figure
    pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,sqrt(meanMap_u_filtered{ind_t(ii)}.^2+meanMap_v_filtered{ind_t(ii)}.^2));
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    caxis([min_magnitude max_magnitude]);
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title(sprintf("Velocity magnitude (t=%.3f[s])",real_t(ii)));

    % vorticity
    if (VORTICITY)
        figure
        pfig = pcolor((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,meanMap_vorticity{ind_t(ii)});
        pfig.EdgeColor = 'none';
        pfig.FaceColor = 'interp';
        colormap("jet");
        colorbar
        caxis([min_vorticity max_vorticity]);
        axis equal
        xlabel("x/H");
        ylabel("y/H");
        title(sprintf("Vorticity (t=%.3f[s])",real_t(ii)));
    end
end

%% Create quiver map
logging_func("Create quiver map");
for ii = 1:num_search_t
    figure
    quiver((mean_x{ind_t(ii)}*1000-object_length)/object_width,mean_y{ind_t(ii)}*1000/object_width,meanMap_u_filtered{ind_t(ii)},meanMap_v_filtered{ind_t(ii)},quiverSize);
    axis equal
    xlabel("x [m]");
    ylabel("y [m]");
    title(sprintf("Flow velocity (t=%.3f[s])",real_t(ii)));
end

%% Create constVelocity meanMap
load(Truck_data_filename,'truck_constVelocity_time_step');

if (truck_constVelocity_time_step)
    x_length = length(mean_x{end});
    y_length = length(mean_y{end});
    all_PIV_length = length(meanMap_u_filtered);
    constVelocity_PIV_data = nan(truck_constVelocity_time_step,y_length,x_length);
    for ii = 1:truck_constVelocity_time_step
        constVelocity_PIV_data(ii,:,:) = meanMap_u_filtered{all_PIV_length-truck_constVelocity_time_step+ii};
    end
    constVelocity_meanMap = squeeze(mean(constVelocity_PIV_data,1));

    %Visualize
    figure
    pfig = pcolor((mean_x{end}*1000-object_length)/object_width,mean_y{end}*1000/object_width,constVelocity_meanMap);
    pfig.EdgeColor = 'none';
    pfig.FaceColor = 'interp';
    colormap("jet");
    colorbar
    axis equal
    xlabel("x/H");
    ylabel("y/H");
    title("Mean Velocity map in constant truck velocity");
end

