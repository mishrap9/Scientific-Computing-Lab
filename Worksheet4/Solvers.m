classdef Solvers
    % This class implements the functions for solving the system.
    methods(Static)
        function A_full = create_A(Nx, Ny)
            % This function creates the A matrix in the full form, using
            % the grid coefficients Nx and Ny.
            N = Nx * Ny;
            e = ones(N, 1);
            % It is easy and fast to create A as a sparse and band matrix,
            % with 5 diagonals, and deleting some elements in the end.
            A_sparse = spdiags([(Ny + 1) ^ 2, (Nx + 1) ^ 2, ...
                                -2 * ((Nx + 1) ^ 2 + (Ny + 1) ^ 2), ...
                                (Nx + 1) ^ 2, (Ny + 1) ^ 2] .* e, ...
                                [-Nx, -1:1, Nx], N, N);
         
            % A is a block tridiagonal matrix => we delete from the bands
            % the elements from (i*Nx+1,i*Nx) and (i*Nx,i*Nx+1), i=1:Ny-1,
            % i.e. the bottom and, respectively, right elements of the
            % bottom-right element from each block on the main diagonal.
            pos_make_zeros = (1:Ny-1) * Nx;
            A_sparse(pos_make_zeros, pos_make_zeros + 1) = 0;
            A_sparse(pos_make_zeros + 1, pos_make_zeros) = 0;
            % create the full matrix from the sparse form
            A_full = full(A_sparse);
        end
     
        function [full_sys, sparse_sys, GS_sys] = solve_system (Nx, Ny)
            % This function unifies the three methods used to solve the
            % system.
            %   Inputs:
            %       Nx, Ny - grid constants on Ox, Oy axis
            %   Outputs:
            %       full_sys - the solution in the full matrix case
            %       sparse_sys - the solution in the sparse matrix case
            %       GS_sys - the solution with Gauss-Seidel method
         
            % Setup
            % create A matrix
            A_full = Solvers.create_A(Nx, Ny);
            A_sparse = sparse(A_full);
            % create b vector
            b = -2 * pi ^ 2 * reshape(sin(pi * (1:Ny) / (Ny + 1))'*...
                                sin(pi * (1:Nx) / (Nx + 1)), Nx * Ny, 1);
         
            % Solvers
            % 1) Direct solution using full matrix
            full_sys = Solvers.full_solver(A_full, b, Nx, Ny);
         
            % 2) Direct solution using sparse representation
            sparse_sys = Solvers.sparse_solver(A_sparse, b, Nx, Ny);
         
            % 3) Iterative solution using Gauss-Seidel
            GS_sys = Solvers.GS_solver(b, Nx, Ny);
        end
     
        function full_sys = full_solver(A_full, b, Nx, Ny)
            % This function computes the direct solution using the full
            % matrix form.
            %   Inputs:
            %       A_full - the full Nx*Ny matrix
            %       b      - the vector with data
            %       Nx, Ny - grid constants on Ox, Oy axis
            %   Outputs:
            %       full_sys - a structure that contains the solution, the
            %                  runtime and the storage.
            tic
            x_full = A_full \ b;
            full_sys.time = toc;
            full_sys.T = [zeros(1, Nx + 2); zeros(Ny, 1), ...
                          reshape(x_full, Ny, Nx), ...
                          zeros(Ny, 1); zeros(1, Nx + 2)];
            full_sys.storage = numel(A_full) + numel(b) + numel(x_full);
        end
     
        function sparse_sys = sparse_solver(A_sparse, b, Nx, Ny)
            % This function computes the direct solution using the sparse
            % representation of the matrix.
            %   Inputs:
            %       A_sparse - the sparse representation of Nx*Ny matrix
            %       b      - the vector with data
            %       Nx, Ny - grid constants on Ox, Oy axis
            %   Outputs:
            %       sparse_sys - a structure that contains the solution,
            %       the runtime and the storage.
            tic
            x_sparse = A_sparse \ b;
            sparse_sys.time = toc;
            sparse_sys.T = [zeros(1, Nx + 2); zeros(Ny, 1), ...
                            reshape(x_sparse, Ny, Nx), ...
                            zeros(Ny, 1); zeros(1, Nx + 2)];
            sparse_sys.storage = nzmax(A_sparse) + numel(b) + numel(x_sparse);
        end
     
        function GS_sys = GS_solver(b, Nx, Ny)
            % This function computes the iterative solution using the
            % Gauss-Seidel method.
            %   Inputs:
            %       b      - the vector with data
            %       Nx, Ny - grid constants on Ox, Oy axis
            %   Outputs:
            %       GS_sys - a structure that contains the solution, the
            %                runtime and the storage.
            tic
            GS_sys.T = GS(b, Nx, Ny);
            GS_sys.time = toc;
            % The returned T is boarded with zeros, so these elements
            % should not be considered.
            GS_sys.storage = numel(b) + (numel(GS_sys.T) - (2 * (Nx + Ny) + 4));
        end
     
    end
end