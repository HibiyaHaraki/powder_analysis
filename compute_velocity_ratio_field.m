function compute_velocity_ratio_field(truck_data_filename,PIV_data_filename,output_folder)
    %
    % compute_velocity_ratio_field
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
    PIV_data = load(PIV_data_filename,'meanMap_u_filtered','meanMap_v_filtered','mean_x','mean_y');

    % Define time step
    time_step = min([length(truck_data.truck_acceleration_mean_vx),length(PIV_data.meanMap_u_filtered)]);

    % Compute velocity ratio
    logging_func("Compute velocity ratio field u");
    reference_velocity_ratio_field = cell(1,time_step);
    for ii = 1:time_step
        reference_velocity_ratio_field(ii) = {PIV_data.meanMap_u_filtered{ii} ./ truck_data.truck_acceleration_mean_vx(ii)};
    end

    logging_func("Compute velocity ratio field v");
    reference_velocity_ratio_field_v = cell(1,time_step);
    for ii = 1:time_step
        reference_velocity_ratio_field_v(ii) = {PIV_data.meanMap_v_filtered{ii} ./ truck_data.truck_acceleration_mean_vx(ii)};
    end

    logging_func("Compute velocity ratio field absolute velocity");
    reference_velocity_ratio_field_absVelocity = cell(1,time_step);
    for ii = 1:time_step
        reference_velocity_ratio_field_absVelocity(ii) = {sqrt(PIV_data.meanMap_u_filtered{ii}.^2 + PIV_data.meanMap_v_filtered{ii}.^2) ./ truck_data.truck_acceleration_mean_vx(ii)};
    end

    % Reference coordinates
    reference_mean_x = PIV_data.mean_x(1:time_step);
    reference_mean_y = PIV_data.mean_y(1:time_step);

    % Save
    save(output_folder+"/"+"PIV_data", ...
        "reference_velocity_ratio_field", "reference_velocity_ratio_field_v", "reference_velocity_ratio_field_absVelocity", ...
        "reference_mean_x", "reference_mean_y", ...
        '-append');
    logging_func(sprintf("Append variable in %s",output_folder+"/"+"PIV_data"));
end