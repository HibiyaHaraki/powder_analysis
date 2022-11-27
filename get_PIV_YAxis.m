function [mean_y,std_y] = get_PIV_YAxis(input_file_names)

    %
    % get_PIV_YAxis
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %

    % Get File information via datastore
    fds = fileDatastore(input_file_names,"ReadFcn",@load_variables,"FileExtensions",".mat");
    num_data = length(fds.Files);
    
    % Import Data
    Data = readall(fds,'UseParallel',true);

    % Check time step
    time_step = length(Data{1}.y);
    for kk = 1:num_data
        if (time_step > length(Data{kk}.y))
            time_step = length(Data{kk}.y);
        end
    end
    
    mean_y = cell(time_step,1);
    std_y  = cell(time_step,1);
    for kk = 1:time_step
        all_y = [];
        for jj = 1:num_data
            original_y_matrix = Data{jj}.y{kk};
            % Check y data order
            if (original_y_matrix(1,1) < original_y_matrix(end,1))
                all_y = [all_y,flipud(original_y_matrix(:,1))];
            else
                all_y = [all_y,original_y_matrix(:,1)];
            end
        end
        all_y = all_y - all_y(end,:);
        mean_y(kk) = {mean(all_y,2)};
        std_y(kk)  = {std(all_y,0,2)};
    end


    function S = load_variables(filename)
        S = load(filename,'y');
    end
end