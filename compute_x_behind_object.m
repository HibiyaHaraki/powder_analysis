function search_x = compute_x_behind_object(PIV_data_filename,Truck_data_filename,search_x_behind_object,search_t)
    % Check Inputs
    if (~exist(PIV_data_filename,'file'))
        error("Cannot find %s",PIV_data_filename);
        return;
    end

    if (~exist(Truck_data_filename,'file'))
        error("Cannot find %s",Truck_data_filename(ii));
        return;
    end
    
    if (~isvector(search_x_behind_object))
        error("Search-x_behind_Object should be vector");
        return;
    end

    if (length(search_x_behind_object) ~= length(search_t))
        error("Search-x_behind_Object and Search-t should be same size");
        return;
    end
    numSearch = length(search_x_behind_object);

    % Import data
    load(PIV_data_filename, ...
        "mean_x", ...
        "meanMap_u_filtered");
    
    load(Truck_data_filename,"sfreq");
    
    time_step = (0:length(mean_x)-1)./sfreq;

    % Convert unit of search pixel
    search_x_behind_object = search_x_behind_object ./ 1000;

    % Search t
    [~,real_t_ind] = min(abs(time_step-search_t'),[],2);

    % Search right edge of the Object
    search_x = zeros(size(search_x_behind_object));
    for ii = 1:numSearch
        u_map = meanMap_u_filtered{real_t_ind(ii)};
        tmp_x_ind = find_right_edge_of_object(u_map);
        right_edge_x = mean_x{real_t_ind(ii)}(tmp_x_ind);
        search_x(ii) = right_edge_x + search_x_behind_object(ii);
    end
end