function [mean_value,std_value] = get_PIV_Data_meanValue(cell_data)
    
    %
    % get_PIV_Data_meanValue
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %

    % Check input
    if (~isvector(cell_data))
        error("You should input cell vector to get_PIV_Data_meanValue!!");
        return;
    end
    
    time_step = length(cell_data);
    mean_value = nan(time_step,1);
    std_value  = nan(time_step,1);
    parfor ii = 1:time_step
        % Get mean map
        map_data = cell_data{ii};

        % Remove outliers
        map_data = filloutliers(map_data,NaN,'quartiles');

        % Compute mean and standard deviation
        mean_value(ii) = mean(map_data,"all","omitnan");
        std_value(ii) = std(map_data,0,'all',"omitnan");
    end
    
end