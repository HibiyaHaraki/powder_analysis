function [interped_reference_velocity_ratio_field_u,interped_reference_velocity_ratio_field_v,interped_reference_velocity_ratio_field_absVelocity] = interp_reference_velocity_ratio(PIV_data_filename,Truck_data_filename)
    %% Explanation
    
    %
    % interp_reference_velocity_ratio
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
    PIV_data = load(PIV_data_filename, ...
        "meanMap_u_filtered", ...
        "mean_x", "mean_y", ...
        "reference_velocity_ratio_field", ...
        "reference_velocity_ratio_field_v", ...
        "reference_velocity_ratio_field_absVelocity", ...
        "reference_mean_x", "reference_mean_y");

    Truck_data = load(Truck_data_filename, ...
        'truck_acceleration_mean_vx');

    %% Compute chracteristic velocity
    if (~isnan(Truck_data.truck_acceleration_mean_vx))
        time_step = min([length(PIV_data.meanMap_u_filtered),length(PIV_data.reference_velocity_ratio_field),length(Truck_data.truck_acceleration_mean_vx)]);
        interped_reference_velocity_ratio_field_u = cell(1,time_step);
        interped_reference_velocity_ratio_field_v = cell(1,time_step);
        interped_reference_velocity_ratio_field_absVelocity = cell(1,time_step);
        for ii = 1:time_step
            % Get PIV coordinates
            experiment_x = PIV_data.mean_x{ii};
            experiment_y = PIV_data.mean_y{ii};
            reference_x  = PIV_data.reference_mean_x{ii};
            reference_y  = PIV_data.reference_mean_y{ii};

            % Scaling
            scale_factor_x = max(experiment_x) ./ max(reference_x);
            scale_factor_y = max(experiment_y) ./ max(reference_y);
            scaled_reference_x = scale_factor_x .* reference_x;
            scaled_reference_y = scale_factor_y .* reference_y;

            % Interpolation
            scaled_reference_velocity_ratio_field = interp2(scaled_reference_x,scaled_reference_y, ...
                PIV_data.reference_velocity_ratio_field{ii}, ...
                experiment_x,experiment_y);
            interped_reference_velocity_ratio_field_u(ii) = {scaled_reference_velocity_ratio_field};

            scaled_reference_velocity_ratio_field = interp2(scaled_reference_x,scaled_reference_y, ...
                PIV_data.reference_velocity_ratio_field_v{ii}, ...
                experiment_x,experiment_y);
            interped_reference_velocity_ratio_field_v(ii) = {scaled_reference_velocity_ratio_field};

            scaled_reference_velocity_ratio_field = interp2(scaled_reference_x,scaled_reference_y, ...
                PIV_data.reference_velocity_ratio_field_absVelocity{ii}, ...
                experiment_x,experiment_y);
            interped_reference_velocity_ratio_field_absVelocity(ii) = {scaled_reference_velocity_ratio_field};
        end
    end
end