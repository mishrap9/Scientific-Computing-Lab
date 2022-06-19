function T_new = Implicit_Euler_step(Nx, Ny, dt, T_old)
%
% This function implements an implicit Euler step for the instationary heat
% equation:
%       T_t = T_xx + T_yy, (1) on the unit square ]0;1[^2,
% with the homogeneous Dirichlet boundary conditions
%       T(x,y,t)=0, (2) for all (x,y) \in \partial ]0;1[^2, t \in ]0;\inf[.
% We used a spatial discretization with (Nx+2) and (Ny+2) grid points on Ox
% and, respectively, Oy axis, where the edges are always 0.
% The Implicit Euler iteration has the form:
%       (I - dt*A) * T_new = T_old (*),
% where A is a sparse matrix obtained by the spatial-discretization and is
% not explicitly computed.
% The linear system (*) is solved iteratively using the Gauss-Seidel (GS)
% method.
% Obs.: In (*), T_new and T_old are viewed, but not explicitely represented
% as long vectors obtained from the respective matrices.
%
% Inputs:
%   - Nx, Ny - grid size constants for Ox and Oy axis
%   - dt - time step for the Implicit Euler method
%   - T_old - the solution at time n*dt; T_old is a (Nx+2)*(Ny+2) matrix,
%             bordered with zeros for an easier implementation. The
%             unknowns are in the positions T_old(2:Nx+1,2:Ny+1).
%
% Outputs:
%   - T_new - the computed solution at time n*dt, with the same structure
%             as T_old
 
    tol = 1e-6; % the tolerance for the residual computed in GS
    R = tol + 1; % the initial residual norm
    itr = 0; % the number of iterations for GS
    T_new = T_old; % the initial guess 
    
    a_ii = - 2 * ((Nx + 1) ^ 2 + (Ny + 1) ^ 2); % the elements over the main diagonal of A
    while R > tol
        itr = itr + 1;
        % compute T_new[itr]
        for i = 2 : Nx + 1
            for j = 2 : Ny + 1
                T_new(i, j) = 1 / (1 - dt * a_ii) * ...
                        (T_old(i, j) + ... % the b-term from GS
                          dt * (Ny + 1) ^ 2 * (T_new(i - 1, j) + ... % up
                                               T_new(i + 1, j)) + ... % down
                          dt * (Nx + 1) ^ 2 * (T_new(i, j - 1) + ... % left
                                               T_new(i, j + 1))); % right
            end
        end
        % compute the residual norm: R = sqrt(1/N*(b-Ax).^2)
        R = 0;
        for i = 2 : Nx + 1
            for j = 2 : Ny + 1
                R = R + (T_old(i, j) - ... % b-term
                  (1 + 2 * dt * ((Nx + 1) ^ 2 + (Ny + 1) ^ 2)) * T_new(i, j) + ...
                  dt * (Ny + 1) ^ 2 * (T_new(i - 1, j) + ... % up
                                     T_new(i + 1, j)) + ... % down
                  dt * (Nx + 1) ^ 2 * (T_new(i, j - 1) + ... % left
                                       T_new(i, j + 1))) .^ 2; % right
            end
        end
        R = sqrt(1 / (Nx * Ny) * R);
    end
end