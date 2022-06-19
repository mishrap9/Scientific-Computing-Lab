function y = runge_kutta_4_meth (f, y0, deltat, tend)
    % euler_meth implements the Runge-Kutta method (4th order) for a
    % general initial value problem: dy(t)/d(t) = f(y(t)), y(0)=y0.
    %
    %   Inputs:
    %     - f      -> the function to evaluate
    %     - y0     -> the initial value
    %     - deltat -> the timestep size of the integration
    %     - tend   -> the end time
    %
    %   Outputs:
    %     - y      -> all aproximated values for y
 
    t = 0 : deltat : tend;
    y = [y0, zeros(1, tend / deltat)];
 
    for i = 1 : length(t) - 1
        Y1 = f(t(i), y(i));
        Y2 = f(t(i) + deltat / 2, y(i) + deltat / 2 * Y1);
        Y3 = f(t(i) + deltat / 2, y(i) + deltat / 2 * Y2);
        Y4 = f(t(i + 1), y(i) + deltat * Y3);
        y(i + 1) = y(i) + deltat / 6 * (Y1 + 2 * Y2 + 2 * Y3 + Y4);
    end
 
end