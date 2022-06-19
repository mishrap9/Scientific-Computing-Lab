function [x, flag] = newton_method(G, dG, x0, eps_newton, Nmax)
    % newton_NL_method implements the iterative Newton method for finding a
    % root of the nonlinear equation G(x)=0.
    % Stop conditions:
    %   - the derivative dG(x) is less than the accuracy limit eps_newton;
    %   - there were more than Nmax steps without finding a concrete result
    %
    %  Inputs:
    %    - G          -> the function that describes the nonlinear equation
    %    - dG         -> the derivative of G
    %    - x0         -> the initial point of searching
    %    - eps_Newton -> the accuracy for the Newton iteration (implicit,
    %                    1e-4)
    %    - Nmax       -> the maximum number of iterations
    %
    %  Outputs:
    %    - x          -> the root found
    %    - flag       -> if no root found
 
    % Default values
    if (nargin < 5)
        Nmax = 100;
    end
    if (nargin < 4)
        Nmax = 100;
        eps_newton = 1e-4;
    end
    flag = 0;
 
    % Implementation of Newton method
    steps = 1;
 
    x = x0 - G(x0) / dG(x0);
    while (norm(x - x0) > eps_newton) && (steps < Nmax)
        x0 = x;
        x = x0 - G(x0) / dG(x0);
        steps = steps + 1;
    end
 
    % Abort the method if there were too many steps without a result
    if (steps >= Nmax)
        flag = 1;
    end
 
end