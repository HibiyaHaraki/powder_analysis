function visualize_meanValue(output_folder_name,sfreq)
    logging_func("Visualize PIV mean value");
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

    % Check sampling frequency
    if (~isscalar(sfreq))
        error("Error: Sampling frequency should be scalar!!");
        return;
    end

    %% Load mean value
    meanValue = load(output_folder_name+"/"+"PIV_data", ...
        'meanValue_u_filtered', ...
        'meanValue_v_filtered');
    time_step = length(meanValue.meanValue_u_filtered);

    %% Visualize
    % mean value of u_filtered
    figure
    plot((0:time_step-1)./sfreq,meanValue.meanValue_u_filtered);
    grid on
    xlabel("t [s]");
    ylabel("Mean flow velocity x [m/s]");
    xlim([0 (time_step-1)/sfreq]);
    title("PIV mean flow velocity x");
    saveas(gcf,output_folder_name+"/"+"meanValue_u_filtered.png");
    logging_func(sprintf("Output %s",output_folder_name+"/"+"meanValue_u_filtered.png"));

    % mean value of v_filtered
    figure
    plot((0:time_step-1)./sfreq,meanValue.meanValue_v_filtered);
    grid on
    xlabel("t [s]");
    ylabel("Mean flow velocity y [m/s]");
    xlim([0 (time_step-1)/sfreq]);
    title("PIV mean flow velocity y");
    saveas(gcf,output_folder_name+"/"+"meanValue_v_filtered.png");
    logging_func(sprintf("Output %s",output_folder_name+"/"+"meanValue_v_filtered.png"));
end