function [mean_data,std_data] = get_PIV_Data_meanMap(input_file_names,data_label)

    %
    % get_PIV_Data_meanMap
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
    
    mean_data = cell(time_step,1);
    std_data  = cell(time_step,1);
    for ii = 1:time_step
        % Get PIV data
        data_size = size(cell2mat(getfield(Data{1},data_label,{ii})));
        all_data = zeros(num_data,data_size(1),data_size(2));

        for jj = 1:num_data

            % Get data
            original_data_matrix = cell2mat(getfield(Data{jj},data_label,{ii}));
            original_x_matrix = Data{jj}.x{ii};
            original_y_matrix = Data{jj}.y{ii};
            typevector_matrix = Data{jj}.typevector_filtered{ii};

            % Check object shape
            object_typevector_y_max = find(typevector_matrix(:,1) == 0,1,'last');
            object_typevector_x_max = find(typevector_matrix(1,:) == 0,1,'last');

            % Check x data order
            if (original_x_matrix(1,1) > original_x_matrix(1,end))
                original_data_matrix = fliplr(original_data_matrix);
                %typevector_matrix = fliplr(type_vector_matrix);
            end
            % Check y data order
            if (original_y_matrix(1,1) < original_y_matrix(end,1))
                original_data_matrix = flipud(original_data_matrix);
                %typevector_matrix = flipud(type_vector_matrix);
            end
            all_data(jj,:,:) = original_data_matrix;
        end

        % Remove Outlier data
        all_data = filloutliers(all_data,nan,'mean',1);

        % Calculate mean and standard deviation
        tmp_mean_data = squeeze(mean(all_data,1,"omitnan"));
        tmp_std_data = squeeze(std(all_data,1,"omitnan"));

        % Make NaN by typevector
        %tmp_mean_data(1:object_typevector_y_max,1:object_typevector_x_max) = NaN;
        %tmp_std_data(1:object_typevector_y_max,1:object_typevector_x_max)  = NaN;
        tmp_mean_data(typevector_matrix == 0) = NaN;
        tmp_std_data(typevector_matrix == 0)  = NaN;

        % Convert mat to cell
        mean_data(ii) = {tmp_mean_data};
        std_data(ii)  = {tmp_std_data};  
    end


    function S = load_variables(filename)
        S = load(filename,'x','y',data_label,'typevector_filtered');
    end
end