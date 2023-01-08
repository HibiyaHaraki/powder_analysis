clear; close all ; clc;
logging_func("Start combinePIVdata.m");
%% Setting
% Input File name
input_filename(1) = "../基準データ_W60固定壁/データ/データ3/data3_4pix_459to1334.mat";
input_filename(2) = "../基準データ_W60固定壁/データ/データ3/data3_8pix_1333to4087.mat";
%input_filename(3) = "../PIV_xo350W45固定壁L1H1/データ/データ1/data1_8pix_2470to5158.mat";

% Unnecessary Data
unnecessary_data(1) = {[]};
unnecessary_data(2) = {[1]};
%unnecessary_data(3) = {[]};

% Output file name
output_filename = "../基準データ_W60固定壁/データ/データ3/data3combine.mat";
%% Input check
if (length(input_filename) < 2)
    error("Error!!: Input file should be more than 2");
    return
end

if (length(input_filename) ~= length(unnecessary_data))
    error("Error: Set same number of unnecessary data");
end

if (output_filename == '')
    error("Please set output file name!!");
    return
end

num_files = length(input_filename);
logging_func(sprintf("Import %d files",num_files));

% Vorticity check
vorticity_check = true;
%{
for ii = 1:num_files
    var_names = who('-file',input_filename(ii),'vorticity');
    if (~any(string(var_names) == "vorticity"))
        warning("%s do not have vorticity data.",input_filename(ii));
        vorticity_check = false;
    end
end
%}

%% Load and combine data
% Ready for variables
u_filtered = {}; u_original = {};
v_filtered = {}; v_original = {};
typevector_filtered = {}; typevector_original = {};
velocity_magnitude = {};
vorticity = {};
x = {}; y = {};

% Load data
for ii = 1:num_files
    logging_func(sprintf("Start loading %s",input_filename(ii)));
    file_struct = load(input_filename(ii));

    num_data = length(file_struct.u_filtered);
    logging_func(sprintf("%s has %d data",input_filename(ii),num_data));
    importing_index = 1:num_data;
    if (~isempty(cell2mat(unnecessary_data(ii))))
        importing_index(ismember(importing_index,cell2mat(unnecessary_data(ii)))) = [];
    end

    logging_func(sprintf("Start to combine data of %s",input_filename(ii)));

    %u_filtered
    u_filtered = [u_filtered;file_struct.u_filtered(importing_index)];

    %u_original
    u_original = [u_original;file_struct.u_original(importing_index)];

    %v_filtered
    v_filtered = [v_filtered;file_struct.v_filtered(importing_index)];

    %v_original
    v_original = [v_original;file_struct.v_original(importing_index)];

    %velocity_magnitude
    velocity_magnitude = [velocity_magnitude;file_struct.velocity_magnitude(importing_index)];

    %typevector_filtered
    for jj = 1:num_data
        typevector_filtered = [typevector_filtered;file_struct.typevector_filtered(1)];
    end

    % vorticity
    if (vorticity_check)
        vorticity = [vorticity;file_struct.vorticity(importing_index)];
    end

    %typevector_original
    for jj = 1:num_data
        typevector_original = [typevector_original;file_struct.typevector_original(1)];
    end

    %x
    x = [x;file_struct.x(importing_index)];

    %y
    y = [y;file_struct.y(importing_index)];

    logging_func(sprintf("Finish loading %s",input_filename(ii)));
end
logging_func(sprintf("Output %d time steps data",length(x)));

%% Output
save(output_filename,"y","x","velocity_magnitude","v_original","v_filtered","u_original","u_filtered","vorticity","typevector_filtered","typevector_original");
logging_func(sprintf("Output %s",output_filename));

%% Functions
function logging_func(message_str)
str_datetime = string(datetime('now','Format','MM/dd HH:mm:ss.SSS'));
fprintf("[%s] %s\n",str_datetime,message_str);
end