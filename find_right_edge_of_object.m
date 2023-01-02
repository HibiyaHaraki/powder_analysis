function edge_ind = find_right_edge_of_object(u_map)
    nan_low = 0;
    nan_col = 0;
    for ii = 1:length(u_map(:,1))
        if (~isnan(u_map(ii,1)))
            nan_low = ii-1;
            break;
        end
    end

    for ii = 1:nan_low
        for jj = 1:length(u_map(1,:))
            if (~isnan(u_map(ii,jj)))
                nan_col = nan_col + (jj - 1);
                break;
            end
        end
    end
    edge_ind = round(nan_col / nan_low);
end