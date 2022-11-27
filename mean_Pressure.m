clear all
close all
clc

disp("meanPressure")
%% 15mm 1.0m/s
disp("(1)Prism 15mm 1.0m/s")
load Prism_15mm_10ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_15mm_10ms_dir_p_mean.mat

clear all
close all

%% 15mm 1.5m/s
disp("(2)Prism 15mm 1.5m/s")
load Prism_15mm_15ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_15mm_15ms_dir_p_mean.mat

clear all
close all

%% 30mm 1.0m/s
disp("(3)Prism 30mm 1.0m/s")
load Prism_30mm_10ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_30mm_10ms_dir_p_mean.mat

clear all
close all

%% 30mm 1.5m/s
disp("(4)Prism 30mm 1.5m/s")
load Prism_30mm_15ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_30mm_15ms_dir_p_mean.mat

clear all
close all

%% 60mm 1.0m/s
disp("(5)Prism 60mm 1.0m/s")
load Prism_60mm_10ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_60mm_10ms_dir_p_mean.mat

clear all
close all

%% 60mm 1.5m/s
disp("(6)Prism 60mm 1.5m/s")
load Prism_60mm_15ms_dir_p.mat
[p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
save Prism_60mm_15ms_dir_p_mean.mat

clear all
close all

%% cylinder
% disp("(7)cylinder 10m/s")
% load cylinder_10ms_leftedge_0_p.mat
% [p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt);
% save cylinder_10ms_leftedge_0_p_mean.mat
% 
% clear all
% close all

disp("Program is end")

function [p_mean, Cp_mean, u_mean, v_mean] = mean_Pressure_Velocity(p, Cp, u, v, Nx, Ny, Nt)
    p_sum = zeros(Ny,Nx);
    Cp_sum = zeros(Ny,Nx);
    u_sum = zeros(Ny,Nx);
    v_sum = zeros(Ny,Nx);
    for k = 1:Nt
        p_sum = p_sum + p{k};
        Cp_sum = Cp_sum + Cp{k};
        u_sum = u_sum + u{k};
        v_sum = v_sum + v{k};
    end

    p_mean = p_sum / Nt;
    Cp_mean = Cp_sum / Nt;
    u_mean = u_sum / Nt;
    v_mean = v_sum / Nt;
end