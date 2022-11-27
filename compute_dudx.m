function dudx = compute_dudx(u_map,x_int)

%
% compute_dudx
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%

logging_func("Compute du/dx");

% Check inputs
map_size = size(u_map);
if (length(map_size) ~= 2)
    error("Error!! Input shoukd be 2-dimentional matrix!!");
    return;
end

if (x_int < 0)
    error("Interval should be positive number!!");
    return;
end

% Compute du/dx
index_counter = [2,1,-1,-2];
dudx = nan(map_size);
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        % Check (ii,jj) value
        if (isnan(u_map(ii,jj)))
            continue;
        end

        % Check neighbor existance
        value_existance_check = false(1,length(index_counter));
        for kk = 1:length(index_counter)
            if (jj+index_counter(kk) < 1 || jj+index_counter(kk) > map_size(2))
                continue;
            else
                if (~isnan(u_map(ii,jj+index_counter(kk))))
                    value_existance_check(kk) = true;
                end
            end
        end

        % Calculate u continuity
        if (value_existance_check(1) && value_existance_check(2) && value_existance_check(3) && value_existance_check(4))
            dudx(ii,jj) = (2*u_map(ii,jj+index_counter(1)) + u_map(ii,jj+index_counter(2)) - u_map(ii,jj+index_counter(3)) - 2*u_map(ii,jj+index_counter(4))) / (10*x_int);
        elseif (value_existance_check(2) && value_existance_check(3))
            dudx(ii,jj) = (u_map(ii,jj+index_counter(2)) - u_map(ii,jj+index_counter(3))) / (2*x_int);
        elseif (value_existance_check(2))
            dudx(ii,jj) = (u_map(ii,jj+index_counter(2)) - u_map(ii,jj)) / (x_int);
        elseif (value_existance_check(3))
            dudx(ii,jj) = (u_map(ii,jj) - u_map(ii,jj+index_counter(3))) / (x_int);
        else
            dudx(ii,jj) = NaN;
        end
    end
end
end