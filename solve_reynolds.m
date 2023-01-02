function [Re_D,Re_H,Re_L] = solve_reynolds(PIV_data_filename,Truck_data_filename)

    %
    % solve_reynolds
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %
    % Compute Reynold's number
    % 
    % Warning
    %  This script needs following files in the same folder.
    %    * logging_func.m
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

    % Constant
    nu = 1.512*10^(-5);

    %% Compute characteristic velocity
    
    % Compute characteristic velocity
    [characteristicVelocity_acceleration,characteristicVelocity_constVelocity] = compute_characteristicVelocity(PIV_data_filename,Truck_data_filename);

    %% Load data
    load(PIV_data_filename,"width","object_width","frontDistance");

    %% Compute Reynold's number (Acceleration)
    Re_D = []; Re_H = []; Re_L = [];
    Re_D_acceleration = NaN; Re_H_acceleration = NaN; Re_L_acceleration = NaN;
    if (~isnan(characteristicVelocity_acceleration))
        % Compute Re_D
        Re_D_acceleration = (characteristicVelocity_acceleration .* width/1000) ./ nu;
        
        % Compute Re_H
        Re_H_acceleration = (characteristicVelocity_acceleration .* object_width/1000) ./ nu;
        
        % Compute Re_L
        Re_L_acceleration = (characteristicVelocity_acceleration .* frontDistance/1000) ./ nu;
        
        Re_D = [Re_D;Re_D_acceleration];
        Re_H = [Re_H;Re_H_acceleration];
        Re_L = [Re_L;Re_L_acceleration];
    end

    %% Compute Reynold's number (constant Velocity)
    Re_D_constVelocity = NaN; Re_H_constVelocity = NaN; Re_L_constVelocity = NaN;
    if (~isnan(characteristicVelocity_constVelocity))
        % Compute Re_D
        Re_D_constVelocity = (characteristicVelocity_constVelocity .* width/1000) ./ nu;
        
        % Compute Re_H
        Re_H_constVelocity = (characteristicVelocity_constVelocity .* object_width/1000) ./ nu;
        
        % Compute Re_L
        Re_L_constVelocity = (characteristicVelocity_constVelocity .* frontDistance/1000) ./ nu;

        Re_D = [Re_D;Re_D_constVelocity];
        Re_H = [Re_H;Re_H_constVelocity];
        Re_L = [Re_L;Re_L_constVelocity];
    end
end