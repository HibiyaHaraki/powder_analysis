function compute_meanValue(output_folder_name)
    logging_func("Compute PIV mean value");
    %% Explanation
    
    %
    % compute_meanValue
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %
    % Compute mean value from mean map
    %
    % Warning
    %  This script needs following files.
    %   * output_folder/PIV_data.mat
    %
    
    %% Check Inputs

    % Check output folder
    if (~exist(output_folder_name,'dir'))
        error("Cannot find folder : %s",output_folder_name);
        return;
    end

    % Check input/output file
    if (~exist(output_folder_name+"/"+"PIV_data.mat",'file'))
        error("Cannot find file : %s",output_folder_name+"/"+"PIV_data.mat");
        return;
    end

    %% Load mean map
    % Get axis data
    mean_x = load(output_folder_name+"/"+"PIV_data",'mean_x');
    mean_y = load(output_folder_name+"/"+"PIV_data",'mean_y');

    % Get mean map
    meanMap = load(output_folder_name+"/"+"PIV_data", ...
        'meanMap_u_filtered', ...
        'meanMap_v_filtered');

    %% Compute mean value
    % u_filtered data
    logging_func("Compute mean value of u_filtered data");
    [meanValue_u_filtered,stdValue_u_filtered] = get_PIV_Data_meanValue(meanMap.meanMap_u_filtered);

    % v_filtered data
    logging_func("Compute mean value of v_filtered data");
    [meanValue_v_filtered,stdValue_v_filtered] = get_PIV_Data_meanValue(meanMap.meanMap_v_filtered);

    %% Save data
    save(output_folder_name+"/"+"PIV_data", ...
        'meanValue_u_filtered', ...
        'meanValue_v_filtered', ...
        '-append');
    logging_func(sprintf("Append variables in %s",output_folder_name+"/"+"PIV_data"));

end