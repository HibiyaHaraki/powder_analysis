function [characteristicVelocity_u,characteristicVelocity_v,characteristicVelocity_absVelocity]= compute_characteristicVelocity_field(PIV_data_filename,Truck_data_filename)
    %% Explanation
    
    %
    % compute_characteristicVelocity_field
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %
    
    %% Check Input

    % Check PIV data file
    if (~exist(PIV_data_filename,'file'))
        error("Cannot find %s",PIV_data_filename);
        return;
    end

    % Check Truck data file
    if (~exist(Truck_data_filename,'file'))
        error("Cannot find %s",Truck_data_filename);
        return;
    end

    %% Import data

    Truck_data = load(Truck_data_filename, ...
        'truck_acceleration_mean_vx');

    %% Create reference velocity ratio for experiment data by interpolation
    [interped_reference_velocity_ratio_field_u,interped_reference_velocity_ratio_field_v,interped_reference_velocity_ratio_field_absVelocity] = interp_reference_velocity_ratio(PIV_data_filename,Truck_data_filename);

    %% Compute characteristic velocity
    time_step = length(interped_reference_velocity_ratio_field_u);
    characteristicVelocity_u = cell(1,time_step);
    characteristicVelocity_v = cell(1,time_step);
    characteristicVelocity_absVelocity = cell(1,time_step);
    for ii = 1:time_step
        characteristicVelocity_u(ii) = {interped_reference_velocity_ratio_field_u{ii} .* Truck_data.truck_acceleration_mean_vx(ii)};
        characteristicVelocity_v(ii) = {interped_reference_velocity_ratio_field_v{ii} .* Truck_data.truck_acceleration_mean_vx(ii)};
        characteristicVelocity_absVelocity(ii) = {interped_reference_velocity_ratio_field_absVelocity{ii} .* Truck_data.truck_acceleration_mean_vx(ii)};
    end

end