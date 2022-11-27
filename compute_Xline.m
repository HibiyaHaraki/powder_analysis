function [search_u_filtered,normalized_y,real_t,real_t_ind,real_x,real_x_ind] = compute_Xline(PIV_data_filename,Truck_data_filename,search_x,search_t)
    % Check Inputs
    if (~exist(PIV_data_filename,'file'))
        error("Cannot find %s",PIV_data_filename);
        return;
    end
    
    if (~exist(Truck_data_filename,'file'))
        error("Cannot find %s",Truck_data_filename);
        return;
    end
    
    if (length(search_x) ~= length(search_t))
        error("Search-x and Search-t should be same size");
        return;
    end
    numSearch = length(search_x);

    % Import data

    load(PIV_data_filename, ...
        "frontDistance", ...
        "mean_x","mean_y", ...
        "reference_velocity_ratio",...
        "meanMap_u_filtered");
    load(Truck_data_filename,"sfreq", ...
        'truck_acceleration_time_step','truck_constVelocity_time_step', ...
        'truck_acceleration_mean_vx','truck_constVelocity_mean_vx');
    
    time_step = (0:length(mean_x)-1)./sfreq;

    % Convert unit of search pixel
    search_x = search_x - frontDistance;
    search_x = search_x ./ 1000;

    % Search t
    [~,real_t_ind] = min(abs(time_step-search_t'),[],2);
    real_t = time_step(real_t_ind)';

    % Search x
    real_x_ind = zeros(numSearch,1);
    real_x = zeros(numSearch,1);
    for ii = 1:numSearch
        [~,real_x_ind(ii)] = min(abs(mean_x{real_t_ind(ii)} - search_x(ii)),[],2);
        real_x(ii) = mean_x{real_t_ind(ii)}(real_x_ind(ii));
        %mean_x{real_t_ind(ii)};
    end

    search_u_filtered = cell(numSearch,1);
    normalized_y = cell(numSearch,1);
    for ii = 1:numSearch
        if (truck_acceleration_time_step < real_t_ind(ii))
            search_u_filtered(ii) = {meanMap_u_filtered{real_t_ind(ii)}(:,real_x_ind(ii)) ./ max(truck_constVelocity_mean_vx)};
        else
            search_u_filtered(ii) = {meanMap_u_filtered{real_t_ind(ii)}(:,real_x_ind(ii)) ./ (reference_velocity_ratio(real_t_ind(ii)) * truck_acceleration_mean_vx(real_t_ind(ii)))};
        end
        normalized_y(ii) = {mean_y{real_t_ind(ii)}*1000};
    end
end