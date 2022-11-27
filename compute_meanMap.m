function compute_meanMap(input_file_names,output_folder_name,frontDistance,width,object_width,object_length)
    logging_func("Compute PIV mean map");
    %% Explanation
    
    %
    % compute_meanMap
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %
    % Compute mean map from PIV data
    %
    % This script needs following function files in same folder
    %   * get_PIV_XAxis.m
    %   * get_PIV_YAxis.m
    %   * get_PIV_Data_meanMap.m
    %   * logging_func.m
    %
    
    %% Check Inputs
    % Check input files
    num_input_files = length(input_file_names);
    for ii = 1:num_input_files
        if (~exist(input_file_names(ii),'file'))
            error("Error: Cannot find %s",input_file_names(ii));
            return;
        end
    end
    
    % Check output folder
    if (~exist(output_folder_name,'dir'))
        mkdir(output_folder_name);
        logging_func(sprintf("Created output folder : %s",output_folder_name));
    end
    
    % Check front distance
    if (~isscalar(frontDistance))
        error("Error: front distance should be scalar!!");
        return;
    end
    
    % Check width
    if (~isscalar(width))
        error("Error: width should be scalar!!");
        return;
    end
    
    %% Import Data
    
    % x data
    logging_func("Compute mean of x data");
    [mean_x,std_x] = get_PIV_XAxis(input_file_names);

    % y data
    logging_func("Compute mean of y data");
    [mean_y,std_y] = get_PIV_YAxis(input_file_names);

    % u_filtered data
    logging_func("Compute mean map of u_filtered data");
    [meanMap_u_filtered,stdMap_u_filtered] = get_PIV_Data_meanMap(input_file_names,"u_filtered");

    % v_filtered data
    logging_func("Compute mean map of v_filtered data");
    [meanMap_v_filtered,stdMap_v_filtered] = get_PIV_Data_meanMap(input_file_names,"v_filtered");

    % u_original data
    logging_func("Compute mean map of u_original data");
    [meanMap_u_original,stdMap_u_original] = get_PIV_Data_meanMap(input_file_names,"u_original");

    % v_original data
    logging_func("Compute mean map of v_original data");
    [meanMap_v_original,stdMap_v_original] = get_PIV_Data_meanMap(input_file_names,"v_original");

    % velocity_magnitude
    %logging_func("Load velocity magnitude");
    %[mean_velocity_magnitude,std_velocity_magnitude] = get_PIV_Data(input_file_names,"velocity_magnitude");

    % Vorticity
    logging_func("Compute mean map of vorticity data");
    [meanMap_vorticity,stdMap_vorticity] = get_PIV_Data_meanMap(input_file_names,"vorticity");

    %% Output
    save(output_folder_name+"/"+"PIV_data", ...
        'width', ...
        'object_width', ...
        'object_length', ...
        'frontDistance', ...
        'mean_x', ...
        'mean_y', ...
        'meanMap_u_filtered', ...
        'meanMap_v_filtered', ...
        'meanMap_u_original', ...
        'meanMap_v_original', ...
        'meanMap_vorticity');
    logging_func(sprintf("Output %s",output_folder_name+"/"+"PIV_data"));
end