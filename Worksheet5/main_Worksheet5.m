% Scientific Computing Lab - Worksheet 5 - Group 16
%
% This is the main file used to solve the tasks from the Worksheet 5.
% 
% In this worksheet we study the instationary heat equation
%       T_t = T_xx + T_yy, (1) on the unit square ]0;1[^2,
% with the homogeneous Dirichlet boundary conditions
%       T(x,y,t)=0, (2) for all (x,y) \in \partial ]0;1[^2, t \in ]0;\inf[,
% and the initial condition
%       T(x,y,0)=1, (3) for all (x,y) \in ]0;1[^2.
% 
% The approach to compute the temperatures is to create a grid for the
% spatial-discretization to obtain a system of ODEs. This system of ODEs in
% integrated in time using the Explicit (tasks b)&c)) and Implicit (tasks
% d)&e)) Euler methods.
%
% - For task c), we create a table with the stable cases and 4 windows that
% subplots the computed solutions at 4 different times (t) with 4 different
% grid sizes (Nx=Ny) and 7 different time steps (delta_t). 
%
% The 4 windows are automatically saved in the /pics folder with the names
% EE_t=?%8.png.
% Moreover, the individual plots are saved in the /pics/Explicit_Euler_imgs
% and, if we want to overwrite them, we set to a non-zero value the
% variable save_subplots_EE. The plots have the names t=1%?_Nx=?_dt=1%?.png
% where the ? denotes the particular value for time, grid size or step
% size.
%
% - For task e), we compute our solutions at the same 4 different times (t)
% and 4 different grid sizes (Nx=Ny), but with only 1 time step delta_t. We  
% plot all solutions in the same window and save the image in the /pics
% folder as IE_dt=1%64.png.
%
% - The implementations of the Implicit and Explicit Euler steps for this 
% problem are in the files Implicit_Euler_step.m and Explicit_Euler_step.m.
% - Solvers.m is the class that implements the 2 Euler methods using the
% above functions.
% - Utilities.m is the class using for creating, editing and saving the
% plots and the table.

%%
clear; clc; close all;
save_subplots_EE = 1; % 1(0) - run the code with(out) generating individual plots
% The plots are stored in  the pics folder
%% a) The steady-state solution
% lim_{t->inf} T(x,y,t) = 0

%% b),c) Explicit Euler
t = [1/8, 2/8, 3/8, 4/8]; % the times
dt = [1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048, 1/4096]; % step-sizes
Nx_vec = [3, 7, 15, 31]; % Ox grid size
Ny_vec = [3, 7, 15, 31]; % Oy grid size
no_steps = (t - [0, t(1:end - 1)])'*(1./dt); % no. of steps between two consecutive solutions

% open the 4 figures that will contain the 28 subplots (4*7)
Utilities.Create_Figures_Explicit_Euler(t);

stable_cases = zeros(length(Nx_vec), length(dt));
for i = 1 : length(Nx_vec)
    % initial condition
    T_0 = [zeros(1, Nx_vec(i) + 2); ...
           [zeros(Ny_vec(i), 1), ones(Nx_vec(i), Nx_vec(i)), zeros(Ny_vec(i), 1)]; ...
           zeros(1, Nx_vec(i) + 2)];
    for j = 1 : length(dt)
        T = T_0;
        for k = 1 : length(t)
            % compute the solution
            T = Solvers.Explicit_Euler(Nx_vec(i), Ny_vec(i), dt(j), T, no_steps(k, j));
            % plot the solution
            Utilities.Plot_sol_EE(T, Nx_vec, Ny_vec, dt, t, i, j, k, save_subplots_EE);
        end
        % check if the solution is stable, i.e, the temperatures decrease
        stable_cases(i, j) = double(sum(sum(double(T <= T_0))) == ...
                                    (Nx_vec(i) + 2) * (Ny_vec(i) + 2));
    end
end
Utilities.Edit_Figs_EE(t);
Table_stable_cases = Utilities.Create_Table_err(stable_cases);
fprintf('Stable cases:\n\n');
disp(Table_stable_cases);

%% d),e) Implicit Euler
t = [1/8, 2/8, 3/8, 4/8]; % the times
dt = 1/64; % step-sizes
Nx_vec = [3, 7, 15, 31]; % Ox grid size
Ny_vec = [3, 7, 15, 31]; % Oy grid size
no_steps = (t - [0, t(1:end - 1)])'*(1./dt); % no. of steps between two consecutive solutions

figIE = Utilities.Create_Figure_Implicit_Euler();
for i = 1 : length(Nx_vec)
    % T - initial condition and the solution's vector
    T = [zeros(1, Nx_vec(i) + 2); ...
         [zeros(Ny_vec(i), 1), ones(Nx_vec(i), Nx_vec(i)), zeros(Ny_vec(i), 1)]; ...
         zeros(1, Nx_vec(i) + 2)];
    for k = 1 : length(t)
        % compute the solution
        T = Solvers.Implicit_Euler(Nx_vec(i), Nx_vec(i), dt, T, no_steps(k, 1));
        % plot the solution
        Utilities.Plot_sol_IE(T, figIE, Nx_vec, Ny_vec, t, i, k);
    end
end
Utilities.Edit_Figs_IE(figIE);
