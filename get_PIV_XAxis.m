function [mean_x,std_x] = get_PIV_XAxis(input_file_names)
    
    %
    % get_PIV_XAxis
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
    time_step = length(Data{1}.x);
    for ii = 1:num_data
        if (time_step > length(Data{ii}.x))
            time_step = length(Data{ii}.x);
        end
    end
    
    % Compute mean and standard deviation of x data
    mean_x = cell(time_step,1);
    std_x  = cell(time_step,1);
    for ii = 1:time_step
        all_x = [];
        for jj = 1:num_data
            original_x_matrix = Data{jj}.x{ii};
            % Check x data order
            if (original_x_matrix(1,1) > original_x_matrix(1,end))
                all_x = [all_x;fliplr(original_x_matrix(1,:))];
            else
                all_x = [all_x;original_x_matrix(1,:)];
            end
        end
        all_x = all_x - all_x(:,1);
        mean_x(ii) = {mean(all_x,1)};
        std_x(ii)  = {std(all_x,0,1)};
    end


    function S = load_variables(filename)
        S = load(filename,'x');
    end
end