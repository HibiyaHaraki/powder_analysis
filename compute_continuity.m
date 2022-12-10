function [u_continuity_map,v_continuity_map,continuity_map] = compute_continuity(u_map,v_map,x_int,y_int)

%
% compute_continuity
%
% created Hibiya Haraki 2022
% All risks of running this script is always with you.
%
% This script needs following scripts
%  * compute_dudx.m
%  * compute_dvdy.m
%  * logging_func.m
%

% Check inputs
map_size = size(u_map);
if (length(map_size) ~= 2)
    error("Error!! Input shoukd be 2-dimentional matrix!!");
    return;
end

if (x_int < 0 || y_int < 0)
    error("Interval should be positive number!!");
    return;
end

% Calulate du/dx
u_continuity_map = compute_dudx(u_map,x_int);

% Calculte dv/dy
v_continuity_map = compute_dvdy(v_map,y_int);

% Calculate dw/dz = -(du/dx + dv/dy)
continuity_map = -(u_continuity_map + v_continuity_map);
end

