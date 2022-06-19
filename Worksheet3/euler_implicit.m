function y = euler_implicit(f, df, y0, deltat, tend, eps_newton, Nmax)
    % euler_implicit implements the implicit Euler method for a general
    % initial value problem: dy(t)/d(t) = f(y(t)), y(0)=y0.
    % The function uses the iterative Newton method for finding a root of
    % the nonlinear equation G(x)=0; for the implicit Euler method, the
    % (n+1)-th solution is obtained using Newton method for the equation:
    % 0=G(x_{n+1})=f(x_{n+1})*delta_t - x_{n+1} + x_n.
    %
    %  Inputs:
    %    - f            -> f(t,y) - the function to evaluate
    %    - df           -> df(t,y)/dy - the derivative of the function
    %    - y0           -> the initial value
    %    - deltat       -> the timestep size of the integration
    %    - tend         -> the end time
    %    - [eps_Newton] -> the accuracy for the Newton iteration (implicit,
    %                      1e-4) [optional arg]
    %    - [Nmax]       -> the maximum number of iterations [optional arg]
    %
    %  Outputs:
    %    - y            -> all aproximated values for y
 
    if (nargin == 6)
        Nmax = 20;
    else
        if (nargin == 5)
            Nmax = 100;
            eps_newton = 1e-4;
        else
            if (nargin < 5)
                error('Not enough input arguments');
            end
        end
    end
 
    t = 0 : deltat : tend;
    y = [y0, zeros(1, tend / deltat)];
 
    for i = 1 : length(t) - 1
        % applying Newton method for finding y(i+1)
        G = @(x)(f(t(i), x) * deltat - x + y(i));
        dG = @(x)(df(t(i), x) * deltat - 1);
        [next_val, flag] = newton_method(G, dG, y(i), eps_newton, Nmax);
        if flag % too many steps in the Newton solution without converging
            y = [];
            fprintf('Implicit Euler: The Newton method do not converge! \\delta_t=%.3d\n\n', deltat);
            return;
        else
            y(i + 1) = next_val;
        end
    end
 
end