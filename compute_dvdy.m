function dvdy = compute_dvdy(v_map,y_int)

%
% compute_dvdy
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%

logging_func("Compute dv/dy");

% Check inputs
map_size = size(v_map);
if (length(map_size) ~= 2)
    error("Error!! Input shoukd be 2-dimentional matrix!!");
    return;
end

if (y_int < 0)
    error("Interval should be positive number!!");
    return;
end

% Compute dvdy
index_counter = [2,1,-1,-2];
dvdy = nan(map_size);
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        % Check (ii,jj) value
        if (isnan(v_map(ii,jj)))
            continue;
        end

        % Check neighbor existance
        value_existance_check = false(length(index_counter),1);
        for kk = 1:length(index_counter)
            if (ii+index_counter(kk) < 1 || ii+index_counter(kk) > map_size(1))
                continue;
            else
                if (~isnan(v_map(ii+index_counter(kk),jj)))
                    value_existance_check(kk) = true;
                end
            end
        end

        % Calculate v continuity
        if (value_existance_check(1) && value_existance_check(2) && value_existance_check(3) && value_existance_check(4))
            dvdy(ii,jj) = (2*v_map(ii+index_counter(1),jj) + v_map(ii+index_counter(2),jj) - v_map(ii+index_counter(3),jj) - 2*v_map(ii+index_counter(4),jj)) / (10*y_int);
        elseif (value_existance_check(2) && value_existance_check(3))
            dvdy(ii,jj) = (v_map(ii+index_counter(2),jj) - v_map(ii+index_counter(3),jj)) / (2*y_int);
        elseif (value_existance_check(2))
            dvdy(ii,jj) = (v_map(ii+index_counter(2),jj) - v_map(ii,jj)) / (y_int);
        elseif (value_existance_check(3))
            dvdy(ii,jj) = (v_map(ii,jj) - v_map(ii+index_counter(3),jj)) / (y_int);
        else
            dvdy(ii,jj) = NaN;
        end
    end
end

end