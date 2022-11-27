function [characteristicVelocity_acceleration,characteristicVelocity_constVelocity] = compute_characteristicVelocity(PIV_data_filename,Truck_data_filename)
    %% Explanation
    
    %
    % compute_characteristicVelocity
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
        "reference_velocity_ratio");

    Truck_data = load(Truck_data_filename, ...
        'truck_acceleration_mean_vx','truck_constVelocity_mean_vx');

    %% Compute chracteristic velocity
    % Acceleration
    characteristicVelocity_acceleration = NaN;
    if (~isnan(Truck_data.truck_acceleration_mean_vx))
        time_step = min([length(PIV_data.meanMap_u_filtered),length(PIV_data.reference_velocity_ratio),length(Truck_data.truck_acceleration_mean_vx)]);
        characteristicVelocity_acceleration = PIV_data.reference_velocity_ratio(1:time_step) .* Truck_data.truck_acceleration_mean_vx(1:time_step);
    end

    % constant velocity
    characteristicVelocity_constVelocity = NaN;
    if (~isnan(Truck_data.truck_constVelocity_mean_vx))
        characteristicVelocity_constVelocity = PIV_data.reference_velocity_ratio(end) .* Truck_data.truck_constVelocity_mean_vx;
    end
end