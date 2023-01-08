clear; close all; clc;
logging_func("Analyze X line data");
%% Explanation

%
% analyze_Xline
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% Analyze a x-line data in mean map
% 
% Warning
%  This script needs following files in the same folder.
%   * logging_func.m
%   * compute_Xline.m
%

%% Setting

% PIV data file name
PIV_data_filename = '../PIV_xo350W60固定壁L3H1/result/PIV_data.mat';

% Truck data file name
Truck_data_filename = "../加速度データ_xo350W60固定壁L3H1/Truck_data.mat";

% Search pixel
search_x = [525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525, 525]; % [mm]

% Search time
search_t = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0]; % [s]

% Normalized option
NORMALIZE_OPTION = 0; % 0-object_width, 1-width


%% Check input

if (~exist(PIV_data_filename,'file'))
    error("Cannot find %s",PIV_data_filename);
    return;
end

if (~exist(Truck_data_filename,'file'))
    error("Cannot find %s",Truck_data_filename);
    return;
end

if (length(search_x) ~= length(search_t))       
    error("Search-x and Search-t should be same size");
    return;
end
numSearch = length(search_x);

%% Get normalize Parameter
load(PIV_data_filename,"object_width");

%% Visualize data
[search_u_filtered,search_v_filtered,normalized_y,real_t,real_t_ind,real_x,real_x_ind] = compute_Xline(PIV_data_filename,Truck_data_filename,search_x,search_t);

% Visualize u
if (NORMALIZE_OPTION == 0)
    load(PIV_data_filename,"object_width");
    for ii = 1:numSearch
        figure
        hold on
        plot(search_u_filtered{ii},normalized_y{ii}./object_width,"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("u/U");
        xlim([-0.6, 1.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        ylim([min(normalized_y{ii}./object_width),max(normalized_y{ii}./object_width)]);
        title(sprintf("u in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        %saveas(gcf,sprintf("XLine_u_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
else
    load(PIV_data_filename,"width");
    for ii = 1:numSearch
        figure
        hold on
        plot(search_u_filtered{ii},normalized_y{ii}./max(normalized_y{ii},[],"all"),"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("u/U");
        xlim([-0.6, 1.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        title(sprintf("u in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        if (NORMALIZE_OPTION == 1)
            ylim([0,1.0]);
        end
        %saveas(gcf,sprintf("XLine_normalize_u_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
end

% Visualize v
if (NORMALIZE_OPTION == 0)
    load(PIV_data_filename,"object_width");
    for ii = 1:numSearch
        figure
        hold on
        plot(search_v_filtered{ii},normalized_y{ii}./object_width,"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("v/U");
        xlim([-0.6, 0.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        ylim([min(normalized_y{ii}./object_width),max(normalized_y{ii}./object_width)]);
        title(sprintf("v in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        %saveas(gcf,sprintf("XLine_v_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
else
    load(PIV_data_filename,"width");
    for ii = 1:numSearch
        figure
        hold on
        plot(search_v_filtered{ii},normalized_y{ii}./max(normalized_y{ii},[],"all"),"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("v/U");
        xlim([-0.6, 0.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        title(sprintf("v in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        if (NORMALIZE_OPTION == 1)
            ylim([-0.6,0.6]);
        end
        %saveas(gcf,sprintf("XLine_normalize_v_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
end

% Visualize absolute velocity
if (NORMALIZE_OPTION == 0)
    load(PIV_data_filename,"object_width");
    for ii = 1:numSearch
        figure
        hold on
        plot(sqrt(search_v_filtered{ii}.^2 + search_u_filtered{ii}.^2),normalized_y{ii}./object_width,"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("Velocity/U");
        xlim([0, 1.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        ylim([min(normalized_y{ii}./object_width),max(normalized_y{ii}./object_width)]);
        title(sprintf("Absolute velocity in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        %saveas(gcf,sprintf("XLine_absVelocity_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
else
    load(PIV_data_filename,"width");
    for ii = 1:numSearch
        figure
        hold on
        plot(sqrt(search_v_filtered{ii}.^2 + search_u_filtered{ii}.^2),normalized_y{ii}./max(normalized_y{ii},[],"all"),"-o");
        xline(0,"b-");
        xline(1,"r-");
        hold off
        grid on
        xlabel("Velocity/U");
        xlim([0, 1.6]);
        if (NORMALIZE_OPTION == 0)
            ylabel("y/H");
        else
            ylabel("y/W");
        end
        title(sprintf("Absolute velocity in x line (t=%.3f[s], x=%.3f[m])",real_t(ii),real_x(ii)));
        if (NORMALIZE_OPTION == 1)
            ylim([0,1.0]);
        end
        %saveas(gcf,sprintf("XLine_normalize_absVelocity_%.3f_%.3f.png",real_t(ii),real_x(ii)));
    end
end
