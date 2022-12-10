function [velocity_fluctuation_acceleration,velocity_fluctuation_constVelocity,real_x,real_y] = compute_velocity_v_fluctuation(search_x,search_y,characteristicVelocity_acceleration,characteristicVelocity_constVelocity,PIV_data_filename,Truck_data_filename)
    %% Explanation
    
    %
    % compute_velocity_v_fluctuation
    %
    % created by Hibiya Haraki 2022
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

    % Check characteristic velocity
    if (~isvector(characteristicVelocity_acceleration) || ~isvector(characteristicVelocity_constVelocity))
        error("Characteristic velocity should be vector!!");
        return;
    end

    % Check search pixel
    if (~isvector(search_x) || ~isvector(search_y))
        error("search_x and search_y should be vector.");
        return;
    end

    if (length(search_x) ~= length(search_y))
        error("search_x and search_y should have same number.");
        return;
    end

    if (length(search_x) < 1 || length(search_y) < 1)
        error("Please input searching pixels!!");
        return;
    end

    if (length(search_x(1,:)) > 1)
        search_x = search_x';
    end

    if (length(search_y(:,1)) > 1)
        search_y = search_y';
    end

    %% Import data
    num_search_pixel = length(search_x);

    PIV_data = load(PIV_data_filename, ...
        "meanMap_v_filtered", ...
        "mean_x","mean_y");

    Truck_data = load(Truck_data_filename, ...
        "sfreq",'truck_acceleration_mean_vx','truck_constVelocity_mean_vx');

    %% Compute velocity fluctuation
    time_step = 0;
    if (~isnan(characteristicVelocity_acceleration))
        time_step = time_step + length(characteristicVelocity_acceleration);
    end

    if (~isnan(characteristicVelocity_constVelocity))
        time_step = time_step + length(characteristicVelocity_constVelocity);
    end

    ind_x = zeros(num_search_pixel,time_step);
    ind_y = zeros(num_search_pixel,time_step);

    real_x = zeros(num_search_pixel,time_step);
    real_y = zeros(num_search_pixel,time_step);

    velocity_fluctuation_acceleration = NaN;
    velocity_fluctuation_constVelocity = NaN;

    % acceleration
    if (~any(isnan(characteristicVelocity_acceleration)))
        velocity_fluctuation_acceleration = nan(num_search_pixel,length(characteristicVelocity_acceleration));
        for ii = 1:length(characteristicVelocity_acceleration)
            % Get coordinates and indexes of search pixel
            % x
            [~,ind_x(:,ii)] = min(abs(PIV_data.mean_x{ii} - search_x),[],2);
            PIV_x = PIV_data.mean_x{ii};
            real_x(:,ii) = PIV_x(1,ind_x(:,ii));
        
            % y
            [~,ind_y(:,ii)] = min(abs(PIV_data.mean_y{ii} - search_y),[],1);
            PIV_y = PIV_data.mean_y{ii};
            real_y(:,ii) = PIV_y(ind_y(:,ii),1);
        
            % Compute velocity fluctuation
            PIV_u_matrix = PIV_data.meanMap_v_filtered{ii};
            velocity_fluctuation_acceleration(:,ii) = diag(PIV_u_matrix(ind_y(:,ii),ind_x(:,ii))) ./ characteristicVelocity_acceleration(ii);
        end
    end

    % constant velocity
    if (~any(isnan(characteristicVelocity_constVelocity)))
        constVelocity_start = 1;
        if (~isnan(characteristicVelocity_acceleration))
            constVelocity_start = constVelocity_start + length(characteristicVelocity_acceleration);
        end
        
        constVelocity_stop = length(PIV_data.meanMap_v_filtered);
        %{
        if (constVelocity_stop > constVelocity_start + length(characteristicVelocity_acceleration))
            constVelocity_stop = constVelocity_start + length(characteristicVelocity_acceleration);
        end
        %}

        velocity_fluctuation_constVelocity = nan(num_search_pixel,length(constVelocity_start:constVelocity_stop));

        for ii = constVelocity_start:constVelocity_stop
            % Get coordinates and indexes of search pixel
            % x
            [~,ind_x(:,ii)] = min(abs(PIV_data.mean_x{ii} - search_x),[],2);
            PIV_x = PIV_data.mean_x{ii};
            real_x(:,ii) = PIV_x(1,ind_x(:,ii));
        
            % y
            [~,ind_y(:,ii)] = min(abs(PIV_data.mean_y{ii} - search_y),[],1);
            PIV_y = PIV_data.mean_y{ii};
            real_y(:,ii) = PIV_y(ind_y(:,ii),1);
        
            % Compute velocity fluctuation
            PIV_u_matrix = PIV_data.meanMap_v_filtered{ii};
            velocity_fluctuation_constVelocity(:,ii-constVelocity_start+1) = diag(PIV_u_matrix(ind_y(:,ii),ind_x(:,ii))) ./ max(Truck_data.truck_constVelocity_mean_vx,[],"all");
        end
    end
end