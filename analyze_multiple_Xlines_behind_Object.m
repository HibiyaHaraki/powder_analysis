clear; close all; clc;
logging_func("Analyze multiple X line data behind the object");
%% Explanation

%
% analyze_multiple_Xlines_behind_Object
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
PIV_data_filename(1) = "../PIV_xo350W60固定壁L1H1/result/PIV_data.mat";
PIV_data_filename(2) = "../PIV_xo350W60移動壁L1H1/result/PIV_data.mat";
PIV_data_filename(3) = "../PIV_xo350W60固定壁L1H1/result/PIV_data.mat";

% Truck data file name
Truck_data_filename(1) = "../加速度データ_xo350W60固定壁L1H1/Truck_data.mat";
Truck_data_filename(2) = "../加速度データ_xo350W60移動壁L1H1/Truck_data.mat";
Truck_data_filename(3) = "../加速度データ_xo350W60固定壁L1H1/Truck_data.mat";

% Search pixel
search_x_behind_Object = [-10, 20, 40]; % [mm] %物体の後ろからの距離

% Search time
search_t = [1, 1, 1]; % [s]

% Figure legend
figure_legend(1) = "W30";
figure_legend(2) = "W45";
figure_legend(3) = "W60";

% Normalized option
NORMALIZE_OPTION = 1; % 0-object_width, 1-width

% x-lim
XLIM = true;
xlim_range = [-0.4 1.5];

%% Check input

if (length(PIV_data_filename) < 2 || length(PIV_data_filename) ~= length(Truck_data_filename))
    error("Please input multiple PIV data file name");
    return;
else
    searchFiles = length(PIV_data_filename);
end

for ii = 1:searchFiles
    if (~exist(PIV_data_filename(ii),'file'))
        error("Cannot find %s",PIV_data_filename(ii));
        return;
    end
end

for ii = 1:searchFiles
    if (~exist(Truck_data_filename(ii),'file'))
        error("Cannot find %s",Truck_data_filename(ii));
        return;
    end
end

if (searchFiles ~= length(figure_legend))
    error("Number of input files and figure_legend should be same size");
    return;
end

if (length(search_x_behind_Object) ~= length(search_t))
    error("Search-x_behind_Object and Search-t should be same size");
    return;
end
numSearch = length(search_x_behind_Object);

%% Compute search_x
search_x = zeros(searchFiles,length(search_x_behind_Object));
for ii = 1:searchFiles
    search_x(ii,:) = compute_x_behind_object(PIV_data_filename(ii),Truck_data_filename(ii),search_x_behind_Object,search_t);
    load(PIV_data_filename(ii),"frontDistance");
    search_x(ii,:) = search_x(ii,:) .* 1000 + frontDistance;
end

%% Visualize data
% Get all u_filtered data
for ii = 1:numSearch
    search_u_filtered = cell(searchFiles,1);
    search_v_filtered = cell(searchFiles,1);
    normalized_y = cell(searchFiles,1);
    for jj = 1:searchFiles
        [search_u_filtered(jj),search_v_filtered(jj),normalized_y(jj),real_t,real_t_ind,real_x,real_x_ind] = compute_Xline(PIV_data_filename(jj),Truck_data_filename(jj),search_x(jj,ii),search_t(ii));
    end

    % Visualize u
    figure
    hold on
    for jj = 1:searchFiles
        if (NORMALIZE_OPTION == 0)
            load(PIV_data_filename(jj),"object_width");
            plot(search_u_filtered{jj},normalized_y{jj}./object_width,"-o");
        else
            load(PIV_data_filename(jj),"width");
            plot(search_u_filtered{jj},normalized_y{jj}./max(normalized_y{jj},[],"all"),"-o");
            ylim([0,1]);
        end
    end
    xline(0,"b-");
    xline(1,"r-");
    hold off
    grid on
    xlabel("u/U");
    if (NORMALIZE_OPTION == 0)
        ylabel("y/H");
    else
        ylabel("y/W");
    end
    if (XLIM)
        xlim(xlim_range);
    end
    legend(figure_legend);
    %ylim([min(normalized_y{ii}),max(normalized_y{ii})]);
    title(sprintf("u in x line (t=%.3f[s], x=%.3f[m])",real_t,real_x(1)));
    %saveas(gcf,sprintf("multiple_XLine_u_%.3f_%.3f.png",real_t(ii),real_x(ii)));

    % Visualize v
    figure
    hold on
    for jj = 1:searchFiles
        if (NORMALIZE_OPTION == 0)
            load(PIV_data_filename(jj),"object_width");
            plot(search_v_filtered{jj},normalized_y{jj}./object_width,"-o");
        else
            load(PIV_data_filename(jj),"width");
            plot(search_v_filtered{jj},normalized_y{jj}./max(normalized_y{jj},[],"all"),"-o");
            ylim([0,1]);
        end
    end
    xline(0,"b-");
    xline(1,"r-");
    hold off
    grid on
    xlabel("v/U");
    if (NORMALIZE_OPTION == 0)
        ylabel("y/H");
    else
        ylabel("y/W");
    end
    if (XLIM)
        xlim(xlim_range);
    end
    legend(figure_legend);
    %ylim([min(normalized_y{ii}),max(normalized_y{ii})]);
    title(sprintf("v in x line (t=%.3f[s], x=%.3f[m])",real_t,real_x(1)));
    %saveas(gcf,sprintf("multiple_XLine_reynolds_v_%.3f_%.3f.png",real_t(ii),real_x(ii)));

    % Visualize absolute velocity
    figure
    hold on
    for jj = 1:searchFiles
        if (NORMALIZE_OPTION == 0)
            load(PIV_data_filename(jj),"object_width");
            plot(sqrt(search_u_filtered{jj}.^2+search_v_filtered{jj}.^2),normalized_y{jj}./object_width,"-o");
        else
            load(PIV_data_filename(jj),"width");
            plot(sqrt(search_u_filtered{jj}.^2+search_v_filtered{jj}.^2),normalized_y{jj}./max(normalized_y{jj},[],"all"),"-o");
            ylim([0,1]);
        end
    end
    xline(0,"b-");
    xline(1,"r-");
    hold off
    grid on
    xlabel("V/U");
    if (NORMALIZE_OPTION == 0)
        ylabel("y/H");
    else
        ylabel("y/W");
    end
    if (XLIM)
        xlim(xlim_range);
    end
    legend(figure_legend);
    %ylim([min(normalized_y{ii}),max(normalized_y{ii})]);
    title(sprintf("Absolute velocity in x line (t=%.3f[s], x=%.3f[m])",real_t,real_x(1)));
    %saveas(gcf,sprintf("multiple_XLine_reynolds_absVelocity_%.3f_%.3f.png",real_t(ii),real_x(ii)));
end
