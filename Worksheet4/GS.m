function T = GS (b, Nx, Ny)
% The function GS implements a Gauss-Seidel for the two-dimensional heat
% equation:
% T_xx+T_yy=-2*pi^2*sin(pi*x)*sin(pi*x), on the unit square ]0;1[^2
% with the homogeneous Dirichlet boundary conditions
% T(x,y)=0, for all (x,y) in \partial ]0;1[^2,
% using (Nx+2)*(Ny+2) grid points for computing the temperature.
% The implemented method do not store the A matrix, because of the known
% structure. Each updated element has at most 4 neigbours: left, right, up
% and bottom, that are considered in the computations if they are
% non-zeros. The current solution benefits from the property of null
% boundary conditions.
%
%   Inputs:
%       - b - the vector with data
%       - Nx+2, Ny+2 - the number of discrete points on the Ox, Oy axis
%   Output:
%       - T - the solution consisting in the temperature
 
    tol = 1e-4;
    N = Nx * Ny;
    x = zeros(N, 1); % the unknowns' vector
    R = sqrt(1 / N * sum(b .^ 2)); % the initial residual norm
    ct_ii = - 1 / (2 * ((Nx + 1) ^ 2 + (Ny + 1) ^ 2)); % the constant 1/A(i,i)
 
    while R > tol
        % first row
        % (1,1)
        x(1) = ct_ii * (b(1) - ...
                (Nx + 1) ^ 2 * x(2) - ... % right
                (Ny + 1) ^ 2 * x(Nx + 1)); % bottom
        % (1, 2 : Nx-1)
        for j = 2 : Nx - 1
            x(j) = ct_ii * (b(j) - ...
                   (Nx + 1) ^ 2 * (x(j - 1) + x(j + 1)) - ... % left & right
                   (Ny + 1) ^ 2 * x(Nx + j)); % bottom
        end
        % (1,Nx)
        x(Nx) = ct_ii * (b(Nx) - ...
                    (Nx + 1) ^ 2 * x(Nx - 1) - ... % left
                    (Ny + 1) ^ 2 * x(Nx + Nx)); % bottom
        % second to penultimate row
        for i = 2 : Ny - 1
            % (i,1)
            x(Nx * (i - 1) + 1) = ct_ii * (b(Nx * (i - 1) + 1) - ...
                        (Nx + 1) ^ 2 * x(Nx * (i - 1) + 2) - ... % right
                        (Ny + 1) ^ 2 * (x(Nx * (i - 2) + 1) + x(Nx * i + 1))); % up & bottom
            % (i, 2 : Nx-1)
            for j = 2 : Nx - 1
                x(Nx * (i - 1) + j) = ct_ii * (b(Nx * (i - 1) + j) - ...
                            (Nx + 1) ^ 2 * (x(Nx * (i - 1) + j - 1) + x(Nx * (i - 1) + j + 1)) - ... % left & right
                            (Ny + 1) ^ 2 * (x(Nx * (i - 2) + j) + x(Nx * i + j))); % up & bottom
            end
            % (i,Nx)
            x(Nx * i) = ct_ii * (b(Nx * i) - ...
                (Nx + 1) ^ 2 * x(Nx * i - 1) - ... % left
                (Ny + 1) ^ 2 * (x(Nx * (i - 1)) + x(Nx * (i + 1)))); % up & bottom
        end
        % last row
        % (Ny,1)
        x(Nx * (Ny - 1) + 1) = ct_ii * (b(Nx * (Ny - 1) + 1) - ...
                        (Nx + 1) ^ 2 * x(Nx * (Ny - 1) + 2) - ... % right
                        (Ny + 1) ^ 2 * x(Nx * (Ny - 2) + 1)); % up
        % (Ny, 2 : Nx-1)
        for j = 2 : Nx - 1
            x(Nx * (Ny - 1) + j) = ct_ii * (b(Nx * (Ny - 1) + j) - ...
                        (Nx + 1) ^ 2 * (x(Nx * (Ny - 1) + j - 1) + x(Nx * (Ny - 1) + j + 1)) - ... % left & right
                        (Ny + 1) ^ 2 * x(Nx * (Ny - 2) + j)); % up
        end
        % (Ny,Nx)
        x(Nx * Ny) = ct_ii * (b(Nx * Ny) - ...
                (Nx + 1) ^ 2 * x(Nx * Ny - 1) - ... % left
                (Ny + 1) ^ 2 * x(Nx * (Ny - 1))); % up
     
        % compute the residual norm
        S = 0; % norm(b-A*x)^2
        % first row
        % (1,1)
        S = S + (b(1) - ...
                1 / ct_ii * x(1) - ...
                (Nx + 1) ^ 2 * x(2) - ... % right
                (Ny + 1) ^ 2 * x(Nx + 1)) ^ 2; % bottom
        % (1, 2 : Nx-1)
        for j = 2 : Nx - 1
            S = S + (b(j) - ...
                    1 / ct_ii * x(j) - ...
                    (Nx + 1) ^ 2 * (x(j - 1) + x(j + 1)) - ... % left & right
                    (Ny + 1) ^ 2 * x(Nx + j)) ^ 2; % bottom
        end
        % (1,Nx)
        S = S + (b(Nx) - ...
                1 / ct_ii * x(Nx) - ...
                (Nx + 1) ^ 2 * x(Nx - 1) - ... % left
                (Ny + 1) ^ 2 * x(Nx + Nx)) ^ 2; % bottom
        % second to penultimate row
        for i = 2 : Ny - 1
            % (i,1)
            S = S + (b(Nx * (i - 1) + 1) - ...
                    1 / ct_ii * x(Nx * (i - 1) + 1) - ...
                    (Nx + 1) ^ 2 * x(Nx * (i - 1) + 2) - ... % right
                    (Ny + 1) ^ 2 * (x(Nx * (i - 2) + 1) + x(Nx * i + 1))) ^ 2; % up & bottom
            % (i, 2 : Nx-1)
            for j = 2 : Nx - 1
                S = S + (b(Nx * (i - 1) + j) - ...
                        1 / ct_ii * x(Nx * (i - 1) + j) - ...
                        (Nx + 1) ^ 2 * (x(Nx * (i - 1) + j - 1) + x(Nx * (i - 1) + j + 1)) - ... % left & right
                        (Ny + 1) ^ 2 * (x(Nx * (i - 2) + j) + x(Nx * i + j))) ^ 2; % up & bottom
            end
            % (i,Nx)
            S = S + (b(Nx * i) - ...
                    1 / ct_ii * x(Nx * i) - ...
                    (Nx + 1) ^ 2 * x(Nx * i - 1) - ... % left
                    (Ny + 1) ^ 2 * (x(Nx * (i - 1)) + x(Nx * (i + 1)))) ^ 2; % up & bottom
        end
        % last row
        % (Ny,1)
        S = S + (b(Nx * (Ny - 1) + 1) - ...
                1 / ct_ii * x(Nx * (Ny - 1) + 1) - ...
                (Nx + 1) ^ 2 * x(Nx * (Ny - 1) + 2) - ... % right
                (Ny + 1) ^ 2 * x(Nx * (Ny - 2) + 1)) ^ 2; % up
        % (Ny, 2 : Nx-1)
        for j = 2 : Nx - 1
            S = S + (b(Nx * (Ny - 1) + j) - ...
                    1 / ct_ii * x(Nx * (Ny - 1) + j) - ...
                    (Nx + 1) ^ 2 * (x(Nx * (Ny - 1) + j - 1) + x(Nx * (Ny - 1) + j + 1)) - ... % left & right
                    (Ny + 1) ^ 2 * x(Nx * (Ny - 2) + j)) ^ 2; % up
        end
        % (Ny,Nx)
        S = S + (b(Nx * Ny) - ...
                1 / ct_ii * x(Nx * Ny) - ...
                (Nx + 1) ^ 2 * x(Nx * Ny - 1) - ... % left
                (Ny + 1) ^ 2 * x(Nx * (Ny - 1))) ^ 2; % up
        R = sqrt(1 / N * S);
    end
 
    % bordering the solution with the nule edges
    T = [zeros(1, Nx + 2); ...
         zeros(Ny, 1), reshape(x, Ny, Nx), zeros(Ny, 1); ...
         zeros(1, Nx + 2)];
 
end