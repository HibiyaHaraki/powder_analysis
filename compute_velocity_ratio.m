function compute_velocity_ratio(truck_data_filename,PIV_data_filename,output_folder)
    %
    % compute_velocity_ratio
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %

    % Input check
    if (~exist(truck_data_filename,'file'))
        error("Error: Cannot find %s",truck_data_filename);
        return;
    end

    if (~exist(PIV_data_filename,'file'))
        error("Error: Cannot find %s",PIV_data_filename);
        return;
    end

    if (~exist(output_folder,'dir'))
        error("Error: Cannot find %s",output_folder);
        return;
    end

    % Load data
    truck_data = load(truck_data_filename,'truck_acceleration_mean_vx');
    PIV_data = load(PIV_data_filename,'meanValue_u_filtered');

    % Define time step
    time_step = min([length(truck_data.truck_acceleration_mean_vx),length(PIV_data.meanValue_u_filtered)]);

    % Compute velocity ratio
    reference_velocity_ratio = PIV_data.meanValue_u_filtered(1:time_step) ./ truck_data.truck_acceleration_mean_vx(1:time_step);

    % Save
    save(output_folder+"/"+"PIV_data", ...
        "reference_velocity_ratio", ...
        '-append');
    logging_func(sprintf("Append variable in %s",output_folder+"/"+"PIV_data"));
end