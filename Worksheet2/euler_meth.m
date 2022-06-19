function y = euler_meth (f, y0, deltat, tend)
    % euler_meth implements the explicit Euler method for a general initial
    % value problem: dy(t)/d(t) = f(y(t)), y(0)=y0.
    %
    %   Inputs:
    %     - f      -> f(y) - the function to evaluate
    %     - y0     -> the initial value
    %     - deltat -> the timestep size of the integration
    %     - tend   -> the end time
    %
    %   Outputs:
    %     - y      -> all aproximated values for y
 
    t = 0 : deltat : tend;
    y = [y0, zeros(1, tend / deltat)];
 
    for i = 1 : length(t) - 1
        y(i + 1) = y(i) + deltat * f(t(i), y(i));
    end
 
end