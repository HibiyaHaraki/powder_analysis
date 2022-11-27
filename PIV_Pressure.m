clear all
close all
clc

% (1) Prism 15mm 1.0m/s
load Prism_15mm_10ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.0;
num = 1;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_15mm_10ms_dir_p.mat

clear all
close all
clc

%% (2) Prism 30mm 1.0m/s
load Prism_30mm_10ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.0;
num = 2;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_30mm_10ms_dir_p.mat

clear all
close all
clc

%% (3) Prism 60mm 1.0m/s
load Prism_60mm_10ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.0;
num = 3;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_60mm_10ms_dir_p.mat

clear all
close all
clc

%% (4) Prism 15mm 1.5m/s
load Prism_15mm_15ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.5;
num = 4;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_15mm_15ms_dir_p.mat

clear all
close all
clc

%% (5) Prism 30mm 1.5m/s
load Prism_30mm_15ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.5;
num = 5;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_30mm_15ms_dir_p.mat

clear all
close all
clc

%% (6) Prism 60mm 1.5m/s
load Prism_60mm_15ms.mat u_component v_component typevector_original x y
L = 15 * 10^-3;
U = 1.5;
num = 6;
[Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
clearvars u_component v_component typevector_original
save Prism_60mm_15ms_dir_p.mat

clear all
close all
clc

% %% (7) cylinder
% load cylinder_10ms.mat u_component v_component typevector_original x y
% L = 20 * 10^-3;
% U = 1.0;
% num = 7;
% [Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
%     cPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
% clearvars cylinder_10ms.mat u_component v_component typevector_original
% save cylinder_10msleft_edge_0_p.mat
% 
% clear all
% close all
% clc
% 
% % %% (8) new cylinder
% % load cylinder_10ms.mat u_component v_component typevector_original x y
% % L = 20 * 10^-3;
% % U = 1.0;
% % num = 8;
% % [Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
% %     sPoisson(x, y, u_component, v_component, typevector_original, L, U, num);
% % clearvars cylinder_10ms.mat u_component v_component typevector_original
% % save cylinder_10ms_slipwall_p.mat
% % 
% % clear all
% % close all
% % clc

disp("Program is end")
%% FUNCTION pPoisson
function [Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    pPoisson(x, y, u_component, v_component, typevector_original, L, U, num)

%     Nt = 10;
    Nt = length(u_component);

    % flipする
    for k = 1:Nt
        u{k} = flip(u_component{k},2);
        v{k} = flip(v_component{k},2);
    end
    x = flip(x{1},2);
    y = flip(y{1},2);

    % エッジ検知のためまだフリップしない
    vectorType = typevector_original{1};

    Nx = size(x,2);
    Ny = size(x,1);

    % 物体表面検出
    [row,col] = find(~vectorType,1);
    Ex = Nx-(col-1)+1;
    Ey = row-1;
    % vectorTypeのフリップ
    vectorType = flip(vectorType,2);

    clearvars u_component v_component typevector_original

    for k = 1:Nt
        for j = 1:Ny
            for i = 1:Nx
                if vectorType(j,i) == 0
                    u{k}(j,i) = 0;
                    v{k}(j,i) = 0;
                end
            end
        end
    end

    dx = abs(x(1,2) - x(1,1));
    dy = abs(y(2,1) - y(1,1));

    steps = 1e+6;
    terms = 1e-6;
    omega = 1.98;
    rho = 1.201;
    mu = 18.2 * 10^-6;
    Re = rho * U * L / mu;

    u_x = cell(Nt,1);
    u_y = cell(Nt,1);
    v_x = cell(Nt,1);
    v_y = cell(Nt,1);
    RHS = cell(Nt,1);
    p = cell(Nt,1);
    pd = cell(Nt,1);
    Cp = cell(Nt,1);
    ul = cell(Nt,1);
    vl = cell(Nt,1);

    % 無次元化
    xl = x / L;
    yl = y / L;
    for k = 1:Nt
        ul{k} = u{k} / U;
        vl{k} = v{k} / U;
    end

    err = 0;

    for k = 1:Nt
        % 最小自乗法（十字型）
        % 速度勾配の算出
        u_x{k}(:,1) = (u{k}(:,2) - u{k}(:,1)) / dx;
        u_x{k}(:,2) = (u{k}(:,3) - u{k}(:,1)) / (2 * dx);
        u_x{k}(:,Nx) = - (u{k}(:,Nx-1) - u{k}(:,Nx)) / dx;
        u_x{k}(:,Nx-1) = - (u{k}(:,Nx-2) - u{k}(:,Nx)) / (2 * dx);

        u_y{k}(1,:) = (u{k}(2,:) - u{k}(1,:)) / dy;
        u_y{k}(2,:) = (u{k}(3,:) - u{k}(1,:)) / (2 * dy);
        u_y{k}(Ny,:) = - (u{k}(Ny-1,:) - u{k}(Ny,:)) / dy;
        u_y{k}(Ny-1,:) = - (u{k}(Ny-2,:) - u{k}(Ny,:)) / (2 * dy);

        v_x{k}(:,1) = (v{k}(:,2) - v{k}(:,1)) / dx;
        v_x{k}(:,2) = (v{k}(:,3) - v{k}(:,1)) / (2 * dx);
        v_x{k}(:,Nx) = - (v{k}(:,Nx-1) - v{k}(:,Nx)) / dx;
        v_x{k}(:,Nx-1) = - (v{k}(:,Nx-2) - v{k}(:,Nx)) / (2 * dx);

        v_y{k}(1,:) = (v{k}(2,:) - v{k}(1,:)) / dy;
        v_y{k}(2,:) = (v{k}(3,:) - v{k}(1,:)) / (2 * dy);
        v_y{k}(Ny,:) = - (v{k}(Ny-1,:) - v{k}(Ny,:)) / dy;
        v_y{k}(Ny-1,:) = - (v{k}(Ny-2,:) - v{k}(Ny,:)) / (2 * dy);

        for j = 1:Ny
            for i = 3:Nx-2
                u_x{k}(j,i) = (2 * u{k}(j,i+2) + u{k}(j,i+1) - u{k}(j,i-1) - 2 * u{k}(j,i-2)) / (10 * dx);
                v_x{k}(j,i) = (2 * v{k}(j,i+2) + v{k}(j,i+1) - v{k}(j,i-1) - 2 * v{k}(j,i-2)) / (10 * dx);
            end
        end


        for j = 3:Ny-2
            for i = 1:Nx
                u_y{k}(j,i) = (2 * u{k}(j+2,i) + u{k}(j+1,i) - u{k}(j-1,i) - 2 * u{k}(j-2,i)) / (10 * dy);
                v_y{k}(j,i) = (2 * v{k}(j+2,i) + v{k}(j+1,i) - v{k}(j-1,i) - 2 * v{k}(j-2,i)) / (10 * dy);
            end
        end

        % 物体内部の速度勾配
        u_x{k}(Ey+1:end, 1:Ex-1) = 0;
        u_y{k}(Ey+1:end, 1:Ex-1) = 0;
        v_x{k}(Ey+1:end, 1:Ex-1) = 0;
        v_y{k}(Ey+1:end, 1:Ex-1) = 0;

        % 物体表面の速度勾配
        u_x{k}(Ey+1:end, Ex) = (u{k}(Ey+1:end, Ex+1) - u{k}(Ey+1:end, Ex)) / dx;
        u_x{k}(Ey+1:end, Ex+1) = (u{k}(Ey+1:end, Ex+2) - u{k}(Ey+1:end, Ex)) / (2 * dx);

        u_y{k}(Ey, 1:Ex-1) = - (u{k}(Ey-1, 1:Ex-1) - u{k}(Ey, 1:Ex-1)) / dy;
        u_y{k}(Ey-1, 1:Ex-1) = - (u{k}(Ey-2, 1:Ex-1) - u{k}(Ey, 1:Ex-1)) / (2 * dy);

        v_x{k}(Ey+1:end, Ex) = (v{k}(Ey+1:end, Ex+1) - v{k}(Ey+1:end, Ex)) / dx;
        v_x{k}(Ey+1:end, Ex+1) = (v{k}(Ey+1:end, Ex+2) - v{k}(Ey+1:end, Ex)) / (2 * dx);

        v_y{k}(Ey, 1:Ex-1) = - (v{k}(Ey-1, 1:Ex-1) - v{k}(Ey, 1:Ex-1)) / dy;
        v_y{k}(Ey-1, 1:Ex-1) = - (v{k}(Ey-2, 1:Ex-1) - v{k}(Ey, 1:Ex-1)) / (2 * dy);

        if mod(k, 100) == 0
            result = sprintf('ポアソン方程式 %d / %d', k, Nt);
            disp(result)
        end
    end

    for k = 1:Nt
        RHS{k} = - 2 .* rho .* (u_x{k} .* v_y{k} - u_y{k} .* v_x{k});
        p{k} = zeros(Ny,Nx);
        pd{k} = zeros(Ny,Nx);
    end

    for k = 1:Nt
        for n = 1:steps
            pd{k} = p{k};
            err2 = err;
            norm2 = 0;
            for j = 2:Ny-1
                for i = 2:Nx-1
                    if (isempty(Ex)~=1) && (isempty(Ey)~=1)
                        if (j == Ey) && (i == Ex)
                            p{k}(j,i) = p{k}(j,i)+omega*((p{k}(j,i+1)+p{k}(j,i-1)+p{k}(j+1,i)+p{k}(j-1,i))/4 - p{k}(j,i) - dx*dx*RHS{k}(j,i)/4);
                        end
                        if (j>Ey-1) && (i<Ex+1)
                            continue
                        end
                    end
                    p{k}(j,i) = p{k}(j,i)+omega*((p{k}(j,i+1)+p{k}(j,i-1)+p{k}(j+1,i)+p{k}(j-1,i))/4 - p{k}(j,i) - dx*dx*RHS{k}(j,i)/4);
                end
            end


    %% 境界条件
    % ノイマン条件
            % 上
    %         p{k}(1,:) = p{k}(2,:);
            % 下
%             p{k}(Ny,:) = p{k}(Ny-1,:);
            % 左
%             p{k}(:,1) = p{k}(:,2);
            % 右
%             p{k}(:,Nx) = p{k}(:,Nx-1);

    %  物体表面（ノイマン条件）
            if (isempty(Ex)~=1) && (isempty(Ey)~=1)
                p{k}(Ey,1:Ex-1) = p{k}(Ey-1,1:Ex-1);
                p{k}(Ey+1:end,Ex) = p{k}(Ey+1:end,Ex+1);
            end

    % ディレクレ条件
            % 上
            p{k}(1,:) = 0;
            % 下
            p{k}(Ny,:) = 0;
            % 左
            p{k}(:,1) = 0;
            % 右
            p{k}(:,Nx) = 0;



            for j = 1:Ny
                for i = 1:Nx
                    norm2 = norm2 + (p{k}(j,i) - pd{k}(j,i))^2;
                    err = sqrt(norm2);
                end
            end

    %         if mod(n, 1000) == 0
    %             result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
    %             disp(result)
    %         end

            if (isnan(err))
                result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
                disp(result)
                disp('err is NaN')
                return
            end

            if err == err2
                result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
                disp(result)
                disp('err is constant')
                return
            end

            if err < terms
    %             if mod(k,100) == 0
                    result = sprintf('(%d) %d / %d, 反復回数 %d, 誤差 %e', num, k, Nt, n, err);
                    disp(result)
    %             end
                break;
            end
        end

        for j = 1:Ny
            for i = 1:Nx
                Cp{k} = p{k} / (1/2 * rho * U^2);
            end
        end

        for j = 1:Ny
            for i = 1:Nx
                if vectorType(j,i) == 0
                    u{k}(j,i) = NaN;
                    v{k}(j,i) = NaN;
                    p{k}(j,i) = NaN;
                    Cp{k}(j,i) = NaN;
                end
            end
        end        
    end

    clearvars err err2 i j k n norm2 pd result RHS steps terms u_x u_y v_x v_y
end




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------




%% Cylinder Poisson
function [Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
    cPoisson(x, y, u_component, v_component, typevector_original, L, U, num)

%     Nt = 10;
    Nt = length(u_component);

    % flipする
    for k = 1:Nt
        u{k} = flip(u_component{k},2);
        v{k} = flip(v_component{k},2);
    end
    x = flip(x{1},2);
    y = flip(y{1},2);

    % エッジ検知のためまだフリップしない
    vectorType = typevector_original{1};

    Nx = size(x,2);
    Ny = size(x,1);

    % 物体表面検出
    [row,col] = find(~vectorType,1);
    Ex = Nx-(col-1)+1;
    Ey = row-1;
    % vectorTypeのフリップ
    vectorType = flip(vectorType,2);

    clearvars u_component v_component typevector_original

    for k = 1:Nt
        for j = 1:Ny
            for i = 1:Nx
                if vectorType(j,i) == 0
                    u{k}(j,i) = 0;
                    v{k}(j,i) = 0;
                end
            end
        end
    end

    dx = abs(x(1,2) - x(1,1));
    dy = abs(y(2,1) - y(1,1));

    steps = 1e+6;
    terms = 1e-6;
    omega = 1.98;
    rho = 1.201;
    mu = 18.2 * 10^-6;
    Re = rho * U * L / mu;

    u_x = cell(Nt,1);
    u_y = cell(Nt,1);
    v_x = cell(Nt,1);
    v_y = cell(Nt,1);
    RHS = cell(Nt,1);
    p = cell(Nt,1);
    pd = cell(Nt,1);
    Cp = cell(Nt,1);
    ul = cell(Nt,1);
    vl = cell(Nt,1);

    % 無次元化
    xl = x / L;
    yl = y / L;
    for k = 1:Nt
        ul{k} = u{k} / U;
        vl{k} = v{k} / U;
    end

    err = 0;

    for k = 1:Nt
        % 最小自乗法（十字型）
        % 速度勾配の算出
        u_x{k}(:,1) = (u{k}(:,2) - u{k}(:,1)) / dx;
        u_x{k}(:,2) = (u{k}(:,3) - u{k}(:,1)) / (2 * dx);
        u_x{k}(:,Nx) = - (u{k}(:,Nx-1) - u{k}(:,Nx)) / dx;
        u_x{k}(:,Nx-1) = - (u{k}(:,Nx-2) - u{k}(:,Nx)) / (2 * dx);

        u_y{k}(1,:) = (u{k}(2,:) - u{k}(1,:)) / dy;
        u_y{k}(2,:) = (u{k}(3,:) - u{k}(1,:)) / (2 * dy);
        u_y{k}(Ny,:) = - (u{k}(Ny-1,:) - u{k}(Ny,:)) / dy;
        u_y{k}(Ny-1,:) = - (u{k}(Ny-2,:) - u{k}(Ny,:)) / (2 * dy);

        v_x{k}(:,1) = (v{k}(:,2) - v{k}(:,1)) / dx;
        v_x{k}(:,2) = (v{k}(:,3) - v{k}(:,1)) / (2 * dx);
        v_x{k}(:,Nx) = - (v{k}(:,Nx-1) - v{k}(:,Nx)) / dx;
        v_x{k}(:,Nx-1) = - (v{k}(:,Nx-2) - v{k}(:,Nx)) / (2 * dx);

        v_y{k}(1,:) = (v{k}(2,:) - v{k}(1,:)) / dy;
        v_y{k}(2,:) = (v{k}(3,:) - v{k}(1,:)) / (2 * dy);
        v_y{k}(Ny,:) = - (v{k}(Ny-1,:) - v{k}(Ny,:)) / dy;
        v_y{k}(Ny-1,:) = - (v{k}(Ny-2,:) - v{k}(Ny,:)) / (2 * dy);

        for j = 1:Ny
            for i = 3:Nx-2
                u_x{k}(j,i) = (2 * u{k}(j,i+2) + u{k}(j,i+1) - u{k}(j,i-1) - 2 * u{k}(j,i-2)) / (10 * dx);
                v_x{k}(j,i) = (2 * v{k}(j,i+2) + v{k}(j,i+1) - v{k}(j,i-1) - 2 * v{k}(j,i-2)) / (10 * dx);
            end
        end


        for j = 3:Ny-2
            for i = 1:Nx
                u_y{k}(j,i) = (2 * u{k}(j+2,i) + u{k}(j+1,i) - u{k}(j-1,i) - 2 * u{k}(j-2,i)) / (10 * dy);
                v_y{k}(j,i) = (2 * v{k}(j+2,i) + v{k}(j+1,i) - v{k}(j-1,i) - 2 * v{k}(j-2,i)) / (10 * dy);
            end
        end

        % 物体内部の速度勾配
        u_x{k}(Ey+1:end, 1:Ex-1) = 0;
        u_y{k}(Ey+1:end, 1:Ex-1) = 0;
        v_x{k}(Ey+1:end, 1:Ex-1) = 0;
        v_y{k}(Ey+1:end, 1:Ex-1) = 0;

        % 物体表面の速度勾配
        u_x{k}(Ey+1:end, Ex) = (u{k}(Ey+1:end, Ex+1) - u{k}(Ey+1:end, Ex)) / dx;
        u_x{k}(Ey+1:end, Ex+1) = (u{k}(Ey+1:end, Ex+2) - u{k}(Ey+1:end, Ex)) / (2 * dx);

        u_y{k}(Ey, 1:Ex-1) = - (u{k}(Ey-1, 1:Ex-1) - u{k}(Ey, 1:Ex-1)) / dy;
        u_y{k}(Ey-1, 1:Ex-1) = - (u{k}(Ey-2, 1:Ex-1) - u{k}(Ey, 1:Ex-1)) / (2 * dy);

        v_x{k}(Ey+1:end, Ex) = (v{k}(Ey+1:end, Ex+1) - v{k}(Ey+1:end, Ex)) / dx;
        v_x{k}(Ey+1:end, Ex+1) = (v{k}(Ey+1:end, Ex+2) - v{k}(Ey+1:end, Ex)) / (2 * dx);

        v_y{k}(Ey, 1:Ex-1) = - (v{k}(Ey-1, 1:Ex-1) - v{k}(Ey, 1:Ex-1)) / dy;
        v_y{k}(Ey-1, 1:Ex-1) = - (v{k}(Ey-2, 1:Ex-1) - v{k}(Ey, 1:Ex-1)) / (2 * dy);

        if mod(k, 100) == 0
            result = sprintf('ポアソン方程式 %d / %d', k, Nt);
            disp(result)
        end
    end
    
    for k = 1:Nt
        RHS{k} = - 2 .* rho .* (u_x{k} .* v_y{k} - u_y{k} .* v_x{k});
        p{k} = zeros(Ny,Nx);
        pd{k} = zeros(Ny,Nx);
    end

    for k = 1:Nt
        for n = 1:steps
            pd{k} = p{k};
            err2 = err;
            norm2 = 0;
            for j = 2:Ny-1
                for i = 2:Nx-1
                    if (isempty(Ex)~=1) && (isempty(Ey)~=1)
                        if (j == Ey) && (i == Ex)
                            p{k}(j,i) = p{k}(j,i)+omega*((p{k}(j,i+1)+p{k}(j,i-1)+p{k}(j+1,i)+p{k}(j-1,i))/4 - p{k}(j,i) - dx*dx*RHS{k}(j,i)/4);
                        end
                        if (j>Ey-1) && (i<Ex+1)
                            continue
                        end
                    end
                    p{k}(j,i) = p{k}(j,i)+omega*((p{k}(j,i+1)+p{k}(j,i-1)+p{k}(j+1,i)+p{k}(j-1,i))/4 - p{k}(j,i) - dx*dx*RHS{k}(j,i)/4);
                end
            end


    %% 境界条件　円柱
    % ノイマン条件
            % 上
            p{k}(1,:) = p{k}(2,:);
            % 下
            p{k}(Ny,:) = p{k}(Ny-1,:);
            % 左
            p{k}(:,1) = p{k}(:,2);
            % 右
            p{k}(:,Nx) = p{k}(:,Nx-1);
            
            % 追加条件
            % 右端
%             p{k}(1,Nx-1:end) = 0;
%             p{k}(Ny,Nx-1:end) = 0;
            % 左端
            p{k}(1,1:2) = 0;
            p{k}(Ny,1:2) = 0;
            
    %  物体表面（ノイマン条件）
            if (isempty(Ex)~=1) && (isempty(Ey)~=1)
                p{k}(Ey,1:Ex-1) = p{k}(Ey-1,1:Ex-1);
                p{k}(Ey+1:end,Ex) = p{k}(Ey+1:end,Ex+1);
            end

    % ディレクレ条件
            % 上
%             p{k}(1,:) = 0;
            % 下
%             p{k}(Ny,:) = 0;
            % 左
%             p{k}(:,1) = 0;
            % 右
%             p{k}(:,Nx) = 0;



            for j = 1:Ny
                for i = 1:Nx
                    norm2 = norm2 + (p{k}(j,i) - pd{k}(j,i))^2;
                    err = sqrt(norm2);
                end
            end

    %         if mod(n, 1000) == 0
    %             result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
    %             disp(result)
    %         end

            if (isnan(err))
                result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
                disp(result)
                disp('err is NaN')
                return
            end

            if err == err2
                result = sprintf('%d / %d, 反復回数 %d, 誤差 %e', k, Nt, n, err);
                disp(result)
                disp('err is constant')
                return
            end

            if err < terms
    %             if mod(k,100) == 0
                    result = sprintf('(%d) %d / %d, 反復回数 %d, 誤差 %e', num, k, Nt, n, err);
                    disp(result)
    %             end
                break;
            end
        end

        for j = 1:Ny
            for i = 1:Nx
                Cp{k} = p{k} / (1/2 * rho * U^2);
            end
        end

        for j = 1:Ny
            for i = 1:Nx
                if vectorType(j,i) == 0
                    u{k}(j,i) = NaN;
                    v{k}(j,i) = NaN;
                    p{k}(j,i) = NaN;
                    Cp{k}(j,i) = NaN;
                end
            end
        end        
    end

    clearvars err err2 i j k n norm2 pd result RHS steps terms u_x u_y v_x v_y

end





%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************




%% Cylinder new Poisson
% function [Cp, dx, dy, mu, Nt, Nx, Ny, omega, p, Re, rho, u, U, ul, v, vectorType, vl, x, xl, y, yl] = ...
%     sPoisson(x, y, u_component, v_component, typevector_original, L, U, num)
% end