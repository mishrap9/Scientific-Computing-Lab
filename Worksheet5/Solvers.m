classdef Solvers
    % This class implements the functions for calling the Euler Methods.
    methods(Static)
     
        function T = Explicit_Euler(Nx, Ny, dt, T_0, no_steps)
        % no_steps - number of steps to do with explicit Euler with dt
        %            no_steps = (tmax-tmin)/dt (= 1/(8*dt) in our problem)
            T = T_0;
            for i = 1 : no_steps
                T = Explicit_Euler_step(Nx, Ny, dt, T);
            end
        end
     
        function T = Implicit_Euler(Nx, Ny, dt, T_0, no_steps)
        % no_steps - number of steps to do with implicit Euler with the
        %            step size dt
        %            no_steps = (tmax-tmin)/dt (= 1/(8*dt) in our problem)
            T = T_0;
            for i = 1 : no_steps
                T = Implicit_Euler_step(Nx, Ny, dt, T);
            end
        end
     
    end
end