clear; clc; close all;

%% d) Solve the system
Nx_vec = [3, 7, 15, 31, 63];
Ny_vec = Nx_vec;
runtime = zeros(3, length(Nx_vec));
storage = zeros(3, length(Nx_vec));
err = zeros(1, length(Nx_vec) + 1);
full_sol = cell(length(Nx_vec), 1);
sparse_sol = cell(length(Nx_vec), 1);
GS_sol = cell(length(Nx_vec), 1);

for i = 1 : length(Nx_vec)
    % solve the system Ax=b with the 3 methods
    [full_sys, sparse_sys, GS_sys] = Solvers.solve_system (Nx_vec(i), Ny_vec(i));
 
    % data structures
    runtime(:, i) = [full_sys.time, sparse_sys.time, GS_sys.time];
    storage(:, i) = [full_sys.storage, sparse_sys.storage, GS_sys.storage];
    full_sol{i} = full_sys.T;
    sparse_sol{i} = sparse_sys.T;
    GS_sol{i} = GS_sys.T;
 
    % g) Comparison between Gauss-Seidel and analytical solution
    T_analytical = sin(pi * (1:Ny_vec(i)) / (Ny_vec(i) + 1))'*...
                   sin(pi * (1:Nx_vec(i)) / (Nx_vec(i) + 1));
 
    % vector with errors
    err(i) = sqrt(1 / (Nx_vec(i) * Ny_vec(i)) * ...
            sum(sum((GS_sys.T(2:end - 1, 2:end - 1) - T_analytical) .^ 2)));
end

%% e) Visualize the solutions
Utilities.Create_Figures(Nx_vec, Ny_vec, GS_sol);

%% f) Tables - runtime & storage
[TableFullMatrix, TableSparseMatrix, TableGS] = ...
            Utilities.Create_Tables_performance(Nx_vec, runtime, storage);
disp(TableFullMatrix.Properties.Description)
disp(TableFullMatrix)
disp(TableSparseMatrix.Properties.Description)
disp(TableSparseMatrix)
disp(TableGS.Properties.Description)
disp(TableGS)

%% g) GS Vs. Analytical solution
% Nx = Ny = 127
Nx = 127; Ny = Nx;
T_analytical = sin(pi * (1:Ny) / (Ny + 1))'*sin(pi*(1:Nx)/(Nx+1));
% create b vector for this particular case
b = - 2 * pi ^ 2 * reshape(sin(pi * (1:Ny) / (Ny + 1))'*...
                           sin(pi * (1:Nx) / (Nx + 1)), Nx * Ny, 1);
GS_sys = Solvers.GS_solver(b, Nx, Ny);

err(end) = sqrt(1 / (Nx * Ny) * ...
        sum(sum((GS_sys.T(2:end - 1, 2:end - 1) - T_analytical) .^ 2)));
err_red = err(1:end - 1) ./ err(2:end);

Table_err = Utilities.Create_Table_err([7, 15, 31, 63, 127], ...
                                       err(2:end), [- 1, err_red(2:end)]);
disp(Table_err.Properties.Description);
disp(Table_err);
