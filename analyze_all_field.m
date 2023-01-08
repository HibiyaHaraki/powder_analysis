clear; close all; clc;
logging_func("Analyze all field data");
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

% Leader PIV data file name
leader_PIV_data_filename = "../PIV_xo350W60固定壁L1H1/result/PIV_data.mat";

% Leader Truck data file name
leader_Truck_data_filename = "../加速度データ_xo350W60固定壁L1H1/truck_data.mat";

% Member PIV data file name
member_PIV_data_filename(1) = "../PIV_xo350W60固定壁L3H1/result/PIV_data.mat";
member_PIV_data_filename(2) = "../PIV_xo350W60固定壁L1H1/result/PIV_data.mat";

% Member Truck data file name
member_Truck_data_filename(1) = "../加速度データ_xo350W60固定壁L3H1/truck_data.mat";
member_Truck_data_filename(2) = "../加速度データ_xo350W60固定壁L1H1/truck_data.mat";

% Output folder
output_folder = "../PIV_xo350W60固定壁L1H1/result";

% Integration time
int_start = 0; % [s]
int_stop  = 3; % [s]

% Comparing range
comparing_range = 50; % [mm]

% Figure legend
figure_legend = ["xo350W60固定壁L1H1", "xo350W60固定壁L3H1", "xo350W60固定壁L1H1"];

%% Check Inputs

% Check leader file
if (~exist(leader_PIV_data_filename,'file'))
    error("Cannot find %s!!",leader_PIV_data_filename);
    return;
end

if (~exist(leader_Truck_data_filename,'file'))
    error("Cannot find %s!!",leader_Truck_data_filename);
    return;
end

% Check member file
if (length(member_PIV_data_filename) ~= length(member_Truck_data_filename))
    error("Please input corresponding member PIV data and member Truck data!!");
    return;
end
num_member_files = length(member_PIV_data_filename);

for ii = 1:num_member_files
    if (~exist(member_PIV_data_filename(ii),'file'))
        error("Cannot find %s!!",member_PIV_data_filename(ii));
        return;
    end
    
    if (~exist(member_Truck_data_filename(ii),'file'))
        error("Cannot find %s!!",member_Truck_data_filename(ii));
        return;
    end
end

% Output folder
if (~exist(output_folder,'dir'))
    error("Cannot find %s!!",output_folder);
    return;
end

% Integration time
if (~isscalar(int_start) || ~isscalar(int_stop) || ~isscalar(comparing_range))
    error("Please input scalar value in integration time!!");
    return;
end

load(leader_Truck_data_filename,"sfreq");

%% Compute characteristic velocity and difference
logging_func("Compute characteristic velocity");

% Leader characteristic velocity
[leader_characteristicVelocity_u,leader_characteristicVelocity_v,leader_characteristicVelocity_absVelocity] = compute_characteristicVelocity_field(leader_PIV_data_filename,leader_Truck_data_filename);
leader_time_step = length(leader_characteristicVelocity_u);

% Member characteristic velocity
member_characteristicVelocity_u = cell(1,num_member_files);
member_characteristicVelocity_v = cell(1,num_member_files);
member_characteristicVelocity_absVelocity = cell(1,num_member_files);
member_time_step = zeros(1,num_member_files);
for ii = 1:num_member_files
    [tmp_member_characteristicVelocity_u,tmp_member_characteristicVelocity_v,tmp_member_characteristicVelocity_absVelocity] = compute_characteristicVelocity_field(member_PIV_data_filename(ii),member_Truck_data_filename(ii));
    member_characteristicVelocity_u(ii) = {tmp_member_characteristicVelocity_u};
    member_characteristicVelocity_v(ii) = {tmp_member_characteristicVelocity_v};
    member_characteristicVelocity_absVelocity(ii) = {tmp_member_characteristicVelocity_absVelocity};
    member_time_step(ii) = length(tmp_member_characteristicVelocity_u);
end

min_time_step = min([leader_time_step,member_time_step]);
logging_func(sprintf("Minimum time step is %d",min_time_step));

logging_func("Compute difference between reference and experiments");

% Leader diffrence between reference and experiment
leader_difference_u = cell(1,leader_time_step);
leader_difference_v = cell(1,leader_time_step);
leader_difference_absVelocity = cell(1,leader_time_step);
leader_PIV_data = load(leader_PIV_data_filename,"meanMap_u_filtered","meanMap_v_filtered","mean_x","mean_y");
for ii = 1:leader_time_step
    tmp_difference_u = zeros(size(leader_PIV_data.meanMap_u_filtered{ii}));
    tmp_difference_v = zeros(size(leader_PIV_data.meanMap_v_filtered{ii}));
    tmp_difference_absVelocity = zeros(size(leader_PIV_data.meanMap_v_filtered{ii}));
    for jj = 1:length(tmp_difference_u(:,1))
        for kk = 1:length(tmp_difference_u(1,:))
            tmp_difference_u(jj,kk) = compare_datas(leader_characteristicVelocity_u{ii}(jj,kk), leader_PIV_data.meanMap_u_filtered{ii}(jj,kk));
            tmp_difference_v(jj,kk) = compare_datas(leader_characteristicVelocity_v{ii}(jj,kk), leader_PIV_data.meanMap_v_filtered{ii}(jj,kk));
            tmp_difference_absVelocity(jj,kk) = compare_datas(leader_characteristicVelocity_absVelocity{ii}(jj,kk), sqrt(leader_PIV_data.meanMap_u_filtered{ii}(jj,kk).^2+leader_PIV_data.meanMap_v_filtered{ii}(jj,kk).^2));
        end
    end
    leader_difference_u(ii) = {tmp_difference_u}; 
    leader_difference_v(ii) = {tmp_difference_v}; 
    leader_difference_absVelocity(ii) = {tmp_difference_absVelocity}; 
end

% Member difference between reference and experiment
member_difference_u = cell(num_member_files,min_time_step);
member_difference_v = cell(num_member_files,min_time_step);
member_difference_absVelocity = cell(num_member_files,min_time_step);
for ll = 1:num_member_files
    member_PIV_data(ll) = load(member_PIV_data_filename(ll),"meanMap_u_filtered","meanMap_v_filtered","mean_x","mean_y");
    for ii = 1:min_time_step
        tmp_difference_u = zeros(size(member_PIV_data(ll).meanMap_u_filtered{ii}));
        tmp_difference_v = zeros(size(member_PIV_data(ll).meanMap_u_filtered{ii}));
        tmp_difference_absVelocity = zeros(size(member_PIV_data(ll).meanMap_u_filtered{ii}));
        for jj = 1:length(tmp_difference_u(:,1))
            for kk = 1:length(tmp_difference_u(1,:))
                tmp_difference_u(jj,kk) = compare_datas(member_characteristicVelocity_u{ll}{ii}(jj,kk), member_PIV_data(ll).meanMap_u_filtered{ii}(jj,kk));
                tmp_difference_v(jj,kk) = compare_datas(member_characteristicVelocity_v{ll}{ii}(jj,kk), member_PIV_data(ll).meanMap_v_filtered{ii}(jj,kk));
                tmp_difference_absVelocity(jj,kk) = compare_datas(member_characteristicVelocity_absVelocity{ll}{ii}(jj,kk), sqrt(member_PIV_data(ll).meanMap_u_filtered{ii}(jj,kk).^2+member_PIV_data(ll).meanMap_v_filtered{ii}(jj,kk).^2));
            end
        end
        member_difference_u(ll,ii) = {tmp_difference_u}; 
        member_difference_v(ll,ii) = {tmp_difference_v}; 
        member_difference_absVelocity(ll,ii) = {tmp_difference_absVelocity}; 
    end
end

%% Interpolation and comparing some experiment data

% Compute standard coordinates
leader_object_right_edge = zeros(1,min_time_step);
max_x_behind_object = zeros(1,min_time_step);
min_x_behind_object_ind = zeros(1,min_time_step);
max_x_behind_object_ind = zeros(1,min_time_step);
leader_x_range = cell(1,min_time_step);
leader_y_range = cell(1,min_time_step);
matrix_size = zeros(1,min_time_step);
for ii = 1:min_time_step
    leader_object_right_edge(ii) = find_right_edge_of_object(leader_PIV_data.meanMap_u_filtered{ii});
    min_x_behind_object_ind(ii) = leader_object_right_edge(ii) + 1;
    max_x_behind_object(ii) = leader_PIV_data.mean_x{ii}(leader_object_right_edge(ii)) + comparing_range;
    [~,max_x_behind_object_ind(ii)] = min(abs(leader_PIV_data.mean_x{ii} - max_x_behind_object(ii)),[],'all');
    leader_x_range(ii) = {leader_PIV_data.mean_x{ii}(min_x_behind_object_ind(ii):max_x_behind_object_ind(ii)) - leader_PIV_data.mean_x{ii}(leader_object_right_edge(ii))};
    leader_y_range(ii) = {leader_PIV_data.mean_y{ii} ./ max(leader_PIV_data.mean_y{ii})};
    matrix_size(ii) = length(leader_x_range{ii}) .* length(leader_y_range{ii});
end

[~,minimum_matrix_ind] = min(matrix_size);
interped_leader_difference_u = cell(1,leader_time_step);
interped_leader_difference_v = cell(1,leader_time_step);
interped_leader_difference_absVelocity = cell(1,leader_time_step);
for ii = 1:min_time_step
    leader_x_range(ii) = leader_x_range(minimum_matrix_ind);
    leader_y_range(ii) = leader_y_range(minimum_matrix_ind);

    leader_x_original = leader_PIV_data.mean_x{ii} - leader_PIV_data.mean_x{ii}(leader_object_right_edge(ii));
    leader_y_original = leader_PIV_data.mean_y{ii} ./ max(leader_PIV_data.mean_y{ii});
    interped_leader_difference_u(ii) = {interp2(leader_x_original,leader_y_original,leader_difference_u{ii},leader_x_range{ii},leader_y_range{ii})}; 
    interped_leader_difference_v(ii) = {interp2(leader_x_original,leader_y_original,leader_difference_v{ii},leader_x_range{ii},leader_y_range{ii})}; 
    interped_leader_difference_absVelocity(ii) = {interp2(leader_x_original,leader_y_original,leader_difference_absVelocity{ii},leader_x_range{ii},leader_y_range{ii})}; 
end

% Interpolation
interped_member_difference_u = cell(num_member_files,min_time_step);
for ll = 1:num_member_files
    for ii = 1:min_time_step
        member_object_right_edge = find_right_edge_of_object(member_PIV_data(ll).meanMap_u_filtered{ii});
        member_x = member_PIV_data(ll).mean_x{ii} - member_PIV_data(ll).mean_x{ii}(member_object_right_edge);
        member_y = member_PIV_data(ll).mean_y{ii} ./ max(member_PIV_data(ll).mean_y{ii});
        interped_member_difference_u(ll,ii) = {interp2(member_x,member_y,member_difference_u{ll,ii},leader_x_range{ii},leader_y_range{ii})};
        interped_member_difference_v(ll,ii) = {interp2(member_x,member_y,member_difference_v{ll,ii},leader_x_range{ii},leader_y_range{ii})};
        interped_member_difference_absVelocity(ll,ii) = {interp2(member_x,member_y,member_difference_absVelocity{ll,ii},leader_x_range{ii},leader_y_range{ii})};
    end
end

%% Compute Integration

int_start_ind = round(int_start*sfreq)+1;
int_stop_ind  = round(int_stop*sfreq)+1;

if (int_start_ind <= 0)
    int_start_ind = 1;
end

if (int_stop_ind > min_time_step)
    int_stop_ind = min_time_step;
end

real_int_start = int_start_ind ./ sfreq;
real_int_stop = int_stop_ind ./ sfreq;

logging_func("List integration result_from here");
% Leader
fprintf("Leader data (%s)\n",leader_PIV_data_filename);
% u
leader_sum_time_u = zeros(1,min_time_step);
leader_sum_result_u = 0;
for ii = 1:min_time_step
    leader_sum_time_u(ii) = sum(interped_leader_difference_u{ii},'all','omitnan');
    leader_sum_result_u = sum(leader_sum_time_u(int_start_ind:int_stop_ind));
end
fprintf(" Integration u : %.3e \n",leader_sum_result_u);

% v
leader_sum_time_v = zeros(1,min_time_step);
leader_sum_result_v = 0;
for ii = 1:min_time_step
    leader_sum_time_v(ii) = sum(interped_leader_difference_v{ii},'all','omitnan');
    leader_sum_result_v = sum(leader_sum_time_v(int_start_ind:int_stop_ind));
end
fprintf(" Integration v : %.3e \n",leader_sum_result_v);

% absolute velocity
leader_sum_time_absVelocity = zeros(1,min_time_step);
leader_sum_result_absVelocity = 0;
for ii = 1:min_time_step
    leader_sum_time_absVelocity(ii) = sum(interped_leader_difference_absVelocity{ii},'all','omitnan');
    leader_sum_result_absVelocity = sum(leader_sum_time_absVelocity(int_start_ind:int_stop_ind));
end
fprintf(" Integration absolute velocity : %.3e \n",leader_sum_result_u);

% Member
member_sum_time_u = zeros(num_member_files,min_time_step);
member_sum_result_u = zeros(1,num_member_files);
member_sum_time_v = zeros(num_member_files,min_time_step);
member_sum_result_v = zeros(1,num_member_files);
member_sum_time_absVelocity = zeros(num_member_files,min_time_step);
member_sum_result_absVelocity = zeros(1,num_member_files);
for ii = 1:num_member_files
    fprintf("Member data %d (%s)\n",ii,member_PIV_data_filename(ii));
    % u
    for jj = 1:min_time_step
        member_sum_time_u(ii,jj) = sum(interped_member_difference_u{ii,jj},'all','omitnan');
        member_sum_result_u(ii) = sum(member_sum_time_u(ii,int_start_ind:int_stop_ind),'all','omitnan');
    end
    fprintf(" Integration u : %.3e \n",member_sum_result_u(ii));

    % v
    for jj = 1:min_time_step
        member_sum_time_v(ii,jj) = sum(interped_member_difference_v{ii,jj},'all','omitnan');
        member_sum_result_v(ii) = sum(member_sum_time_v(ii,int_start_ind:int_stop_ind),'all','omitnan');
    end
    fprintf(" Integration v : %.3e \n",member_sum_result_v(ii));

    % absolute velocity
    for jj = 1:min_time_step
        member_sum_time_absVelocity(ii,jj) = sum(interped_member_difference_absVelocity{ii,jj},'all','omitnan');
        member_sum_result_absVelocity(ii) = sum(member_sum_time_absVelocity(ii,int_start_ind:int_stop_ind),'all','omitnan');
    end
    fprintf(" Integration absolute velocity : %.3e \n",member_sum_result_u(ii));
end

%% Visualize
time_domain = (0:(min_time_step-1)) ./ sfreq;

figure
hold on
plot(time_domain,leader_sum_time_u);
for ii = 1:num_member_files
    plot(time_domain,member_sum_time_u(ii,:))
end
xline([real_int_start,real_int_stop],"r-");
hold off
grid on
xlabel("Time [s]");
ylabel("Difference u");
legend(figure_legend,'Location','best');

figure
hold on
plot(time_domain,leader_sum_time_v);
for ii = 1:num_member_files
    plot(time_domain,member_sum_time_v(ii,:))
end
xline([real_int_start,real_int_stop],"r-");
hold off
grid on
xlabel("Time [s]");
ylabel("Difference v");
legend(figure_legend,'Location','best');

figure
hold on
plot(time_domain,leader_sum_time_absVelocity);
for ii = 1:num_member_files
    plot(time_domain,member_sum_time_absVelocity(ii,:))
end
xline([real_int_start,real_int_stop],"r-");
hold off
grid on
xlabel("Time [s]");
ylabel("Difference absolute velocity");
legend(figure_legend,'Location','best');

%% Function for comparing data
function answer = compare_datas(reference,experiment)
    answer = (experiment - reference)^2;  % 出力したい数式
end
