function p = solve_Poisson_rectangle(u_map,v_map,x_int,y_int,CHECK)

%
% solve_Poisson_rectangle
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% This script needs following scripts
%  * compute_dudx.m
%  * compute_dvdy.m
%  * logging_func.m
%
% Warning
%  This script assume that object shape is rectangle and the position is upper-left!!
%

logging_func("Solve Poisson equation");

% Neuman Condition
p_top = 0;
p_bottom = 0;
p_surface = 0;

% Constants
rho = 1.201;
map_size = size(u_map);

% Output matrix
p = nan(map_size);

% Compute du/dx
dudx = compute_dudx(u_map,x_int);

% Compute du/dy
dudy = compute_dvdy(u_map,y_int);

% Compute dv/dx
dvdx = compute_dudx(v_map,x_int);

% Compute dv/dy
dvdy = compute_dvdy(v_map,y_int);

% Compute Q
Q = -2*rho*(dudx.*dvdy - dudy.*dvdx);

% Check object boundary (Object should be upper-left)
logging_func("Check object");
object_x_start = 1;
object_x_stop  = find(isnan(u_map(1,:)),1,"last");
object_y_start = 1;
object_y_stop  = find(isnan(u_map(:,1)),1,"last");
logging_func(sprintf("Object x start : %d",object_x_start));
logging_func(sprintf("Object x stop  : %d",object_x_stop ));
logging_func(sprintf("Object y start : %d",object_y_start));
logging_func(sprintf("Object y stop  : %d",object_y_stop ));

% Compute A*p=b
logging_func("Generate matrixes");
A = sparse(zeros(map_size(1)*map_size(2)));
b = zeros(map_size(1)*map_size(2),1);
h = x_int;

% Specify relationship between A and index 
matrix_index = zeros(map_size(1)*map_size(2),2);
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        matrix_index((ii-1)*map_size(2)+jj,1) = jj;
        matrix_index((ii-1)*map_size(2)+jj,2) = ii;
    end
end

% Set coefficients
logging_func("Set coefficients");
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        % Set -1/4*h^2*Q to b
        b((ii-1)*map_size(2)+jj,1) = -h^2*Q(ii,jj);

        % Set Dirichlet boundary condition
        % Left
        if (ii > object_y_stop && jj == 1)
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj  )) =  1;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj+1)) = -1;
            b((ii-1)*map_size(2)+jj,1) = 0;
            continue;
        end

        % Right
        if (jj == map_size(2))
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj  )) =  1;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj-1)) = -1;
            b((ii-1)*map_size(2)+jj,1) = 0;
            continue;
        end

        % Set Equation
        if (ii > 1 && ii < map_size(1) && jj > object_x_start + 1 && jj < map_size(2))
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj  )) = 4;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj-1)) = -1;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj+1)) = -1;
            A((ii-1)*map_size(2)+jj,(ii  )*map_size(2)+(jj  )) = -1;
            A((ii-1)*map_size(2)+jj,(ii-2)*map_size(2)+(jj  )) = -1;
            continue;
        end

        if (ii > object_y_stop && ii < map_size(1) && jj > 1 && jj <= object_x_start + 1)
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj  )) = 4;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj-1)) = -1;
            A((ii-1)*map_size(2)+jj,(ii-1)*map_size(2)+(jj+1)) = -1;
            A((ii-1)*map_size(2)+jj,(ii  )*map_size(2)+(jj  )) = -1;
            A((ii-1)*map_size(2)+jj,(ii-2)*map_size(2)+(jj  )) = -1;
            continue;
        end
    end
end

% Set Neuman boundary condition
logging_func("Set Neuman boundary condition");
delete_count = 0;
delete_list = [];
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        % Set Neuman boundary condition
        % Top
        if (ii == 1 && jj > object_x_start)
            delete_count = delete_count + 1;
            delete_list(delete_count) = (ii-1)*map_size(2)+jj;
            b(:,1) = b(:,1) - p_top*A(:,(ii-1)*map_size(2)+jj);
            continue;
        end
        % Bottom
        if (ii == map_size(1))
            delete_count = delete_count + 1;
            delete_list(delete_count) = (ii-1)*map_size(2)+jj;
            b(:,1) = b(:,1) - p_bottom*A(:,(ii-1)*map_size(2)+jj);
            continue;
        end
        % Object Surface
        if (ii == object_y_stop + 1 && jj >= object_x_start && jj <= object_x_stop + 1)
            delete_count = delete_count + 1;
            delete_list(delete_count) = (ii-1)*map_size(2)+jj;
            b(:,1) = b(:,1) - p_surface*A(:,(ii-1)*map_size(2)+jj);
            continue;
        end
        if (ii >= object_y_start && ii <= object_y_stop + 1 && jj == object_x_stop + 1)
            delete_count = delete_count + 1;
            delete_list(delete_count) = (ii-1)*map_size(2)+jj;
            b(:,1) = b(:,1) - p_surface*A(:,(ii-1)*map_size(2)+jj);
            continue;
        end
    end
end


% Find Object
logging_func("Find object");
for ii = 1:map_size(1)
    for jj = 1:map_size(2)
        % Delete object
        if (ii >= object_y_start && ii <= object_y_stop && jj >= object_x_start && jj <= object_x_stop)
            delete_count = delete_count + 1;
            delete_list(delete_count) = (ii-1)*map_size(2)+jj;
        end
    end
end

% Check deleting points
if (CHECK)
    figure
    plot(matrix_index(:,1),matrix_index(:,2),'k.');
    hold on
    plot(matrix_index(delete_list,1),matrix_index(delete_list,2),'r.');
    plot(matrix_index(isnan(b),1),matrix_index(isnan(b),2),'ro');
    hold off
    grid on
    axis equal
    legend("Nodes","Deleting nodes","b is NaN");
    set(gca,'Ydir','reverse');
    title("Deleting elements");

    figure
    spy(A);
end



% Delete unnecessary element
logging_func("Delete unnecessary element");
delete_list = unique(delete_list);
logging_func(sprintf("Delete %d elements",length(delete_list)));
A(delete_list,:)  = [];
A(:,delete_list) = [];
b(delete_list) = [];
matrix_index(delete_list,:) = [];

% Check result
logging_func("Compute linear system");
logging_func(sprintf("Matrix A > Size: %d, Rank: %d",length(A(:,1)),sprank(A)));
logging_func(sprintf("Vector b includes %d NaN",nnz(isnan(b))));

% Compute linear system
result_p = A\b;

logging_func("Output pressure");
for ii = 1:length(result_p)
    p(matrix_index(ii,2),matrix_index(ii,1)) = result_p(ii);
end
end