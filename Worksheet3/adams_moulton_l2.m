function y = adams_moulton_l2(y0, deltat, tend, eps_newton, Nmax)
    % Adams_Moulton implements the second order Adams_Moulton method for
    % the solution of the initial value problem:
    %    dy(t)/d(t) = f(y(t)), y(0)=y0,
    % using the linearisation (l2 in the WS):
    %   y_{n+1}=y_n+delta_t/2*(7*(1-y_n/10)*y_n+7*(1-y_{n}/10)*y_{n+1}).
    %
    %   Inputs:
    %     - y0     -> the initial value
    %     - deltat -> the timestep size of the integration
    %     - tend   -> the end time
    %
    %   Outputs:
    %     - y      -> all aproximated values for y
 
    if (nargin == 4)
        Nmax = 20;
    else
        if (nargin == 3)
            Nmax = 100;
            eps_newton = 1e-4;
        else
            if (nargin < 3)
                error('Not enough input arguments');
            end
        end
    end
 
    t = 0 : deltat : tend;
    y = [y0, zeros(1, tend / deltat)];
 
    for i = 1 : length(t) - 1
        % applying Newton method for finding y(i+1)
        dG = @(x)(1 + 7 / 20 * deltat * y(i) - 7 / 2 * deltat); % the derivative is constant over the Newton iterations
        G = @(x)(x - y(i) - deltat / 2 * 7 * (1 - y(i) / 10) * (y(i) + x));
        [next_val, flag] = newton_method(G, dG, y(i), eps_newton, Nmax);
        if flag % too many steps in the Newton solution without converging
            y = [];
            fprintf('Adams-Moulton - linearisation 1: The Newton method do not converge! \\delta_t=%.3d\n\n', deltat);
            return;
        else
            y(i + 1) = next_val;
        end
    end
 
    % Explicit formula
    % y(i+1) = (1+7/2*deltat*(1-y(i)/10))/(1-7/2*deltat*(1-y(i)/10)) * y(i);
 
end