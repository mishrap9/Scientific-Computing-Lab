classdef problem_functions
    methods(Static)
        function p = sol_p (t)
            % Analytical solution
            p = 200 ./ (20 - 10 * exp(- 7 * t));
        end
     
        function dp = func_p(t, p)
            % defining the ODE to be solved numerically
            dp = 7 * (1 - p / 10) * p;
        end
     
        function ddp = dfunc_p(t, p)
            % defining the function's derivative
            ddp = 7 - 7 / 5 * p;
        end
    end
end