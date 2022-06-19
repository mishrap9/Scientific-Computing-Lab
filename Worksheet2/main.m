% Lab Course Scientific Computing - Worksheet 2, Group 16
clear; clc; close all;
%% Figure properties constants
labelsize = 18;
legFontSize = 18; % legend font size
titleFontSize = 20;
LineWidthAnalytical = 0.5;
MarkerSize = 7;
%% a) Plot the analytical solution
% defining the time interval for computing the analytical solution
t0 = 0;
tend = 10;
tstep = 1e-2;
t = t0 : tstep : tend;
p = sol_p(t);

Fig1 = figure(1); grid on; hold on;
set (Fig1, 'NumberTitle', 'off', ...
  'Name', 'Analytical solution', ...
  'Units', 'normalized', ...
  'OuterPosition', [0 0 1 1]);
p1 = plot(t, p, 'r');
p1.LineWidth = LineWidthAnalytical;
xlab = xlabel('$t$', 'FontSize', labelsize);
set(xlab, 'Interpreter', 'latex');
ylab = ylabel('$p(t)$', 'FontSize', labelsize);
set(ylab, 'Interpreter', 'latex');
leg1 = legend('$p(t)$');
leg1.FontSize = legFontSize;
leg1.Location = 'east';
set(leg1, 'Interpreter', 'latex');
title('a) Graph of $p(t)=\frac{10}{1+9e^{-t}}$', 'interpreter', 'latex', 'FontSize', titleFontSize);
%% c) Approximation error for deltat={1,1/2,1/4,1/8}
tend = 5;
t = 0 : tstep : tend; % time for ploting analytical p
p = sol_p(t);
% setting figures properties
Fig2 = figure(2); % figure for explicit Euler method
set(Fig2, 'NumberTitle', 'off', ...
  'Name', 'Explicit Euler method', ...
  'Units', 'normalized', ...
  'OuterPosition', [0 0 1 1]);
p2 = plot(t, p, 'r');
hold on;
grid on;
p2.LineWidth = LineWidthAnalytical;
xlab = xlabel('t', 'FontSize', labelsize);
set(xlab, 'Interpreter', 'latex');
ylab = ylabel('p(t)', 'FontSize', labelsize);
set(ylab, 'Interpreter', 'latex');
title('c) Explicit Euler with time steps $\delta (t)=\{1,\frac{1}{2},\frac{1}{4},\frac{1}{8}\}$', ...
  'interpreter', 'latex', 'FontSize', titleFontSize);

Fig3 = figure(3); % figure for method of Heun
set(Fig3, 'NumberTitle', 'off', ...
  'Name', 'Heun method', ...
  'Units', 'normalized', ...
  'OuterPosition', [0 0 1 1]);
p3 = plot(t, p, 'r');
hold on;
grid on;
p3.LineWidth = LineWidthAnalytical;
xlab = xlabel('t', 'FontSize', labelsize);
set(xlab, 'Interpreter', 'latex');
ylab = ylabel('p(t)', 'FontSize', labelsize);
set(ylab, 'Interpreter', 'latex');
title('c) Heun with time steps $\delta (t)=\{1,\frac{1}{2},\frac{1}{4},\frac{1}{8}\}$', ...
  'interpreter', 'latex', 'FontSize', titleFontSize);

Fig4 = figure(4); % figure for Runge-Kutta method (4th order)
set(Fig4, 'NumberTitle', 'off', ...
  'Name', 'Runge-Kutta method (4th order)', ...
  'Units', 'normalized', ...
  'OuterPosition', [0 0 1 1]);
p4 = plot(t, p, 'r');
hold on;
grid on;
p4.LineWidth = LineWidthAnalytical;
xlab = xlabel('t', 'FontSize', labelsize);
set(xlab, 'Interpreter', 'latex');
ylab = ylabel('p(t)', 'FontSize', labelsize);
set(ylab, 'Interpreter', 'latex');
title('c) Runge-Kutta method (4th order) with time steps $\delta (t)=\{1,\frac{1}{2},\frac{1}{4},\frac{1}{8}\}$', ...
  'interpreter', 'latex', 'FontSize', titleFontSize);

% Different markers with different colours
markers = {'o', '+', 'd', '*'};
colors = {'b', 'k', 'm', 'g'};

p0 = 1; % initial condition
deltat = [1, 1/2, 1/4, 1/8];
E = zeros(3, length(deltat)); % each line of E contains the approximation error for a specific method
for i = 1 : length(deltat)
    t_deltas = 0 : deltat(i) : tend;
    p_deltas = sol_p(t_deltas); % the analytical solution with time step deltat(i)
 
    y_euler = euler_meth(@func_p, p0, deltat(i), tend); % Euler method
    figure(2);
    p2 = plot(t_deltas, y_euler, markers{i});
    p2.MarkerSize = MarkerSize;
    p2.Color = colors{i};
 
    y_heun = heun_meth(@func_p, p0, deltat(i), tend); % Heun method
    figure(3);
    p3 = plot(t_deltas, y_heun, markers{i});
    p3.MarkerSize = MarkerSize;
    p3.Color = colors{i};
 
    y_RK4 = runge_kutta_4_meth(@func_p, p0, deltat(i), tend); % Runge-Kutta method
    figure(4); p4 = plot(t_deltas, y_RK4, markers{i});
    p4.MarkerSize = MarkerSize;
    p4.Color = colors{i};
    E(:, i) = sqrt(deltat(i) / tend) * sqrt([sum((y_euler - p_deltas) .^ 2); ...
    sum((y_heun - p_deltas) .^ 2); ...
      sum((y_RK4 - p_deltas) .^ 2)]);
end

figure(2);
leg2 = legend('$p(t)$', '$\delta t=1$', '$\delta t=1/2$', '$\delta t=1/4$', '$\delta t=1/8$');
set(leg2, 'Interpreter', 'latex');
leg2.FontSize = legFontSize;
leg2.Location = 'southeast';

figure(3);
leg3 = legend('$p(t)$', '$\delta t=1$', '$\delta t=1/2$', '$\delta t=1/4$', '$\delta t=1/8$');
set(leg3, 'Interpreter', 'latex');
leg3.FontSize = legFontSize;
leg3.Location = 'southeast';

figure(4);
leg4 = legend('$p(t)$', '$\delta t=1$', '$\delta t=1/2$', '$\delta t=1/4$', '$\delta t=1/8$');
set(leg4, 'Interpreter', 'latex');
leg4.FontSize = legFontSize;
leg4.Location = 'southeast';
%% d) Error reduced
deltat_h = deltat/2;
E_red = zeros(3, length(deltat_h)); % each line of E contains the approximation error for a specific method
for i = 1 : length(deltat_h)-1
    E_red(:,i+1) = E(:,i)./E(:,i+1);
end
%% e) Error approximation
deltat_best = deltat(end);
t_deltat_best = 0 : deltat_best : tend;
y_euler_best = euler_meth(@func_p, p0, deltat_best, tend); % Euler method
y_heun_best = heun_meth(@func_p, p0, deltat_best, tend); % Heun method
y_RK4_best = runge_kutta_4_meth(@func_p, p0, deltat_best, tend); % Runge-Kutta method
E_tilde = zeros(3, length(deltat)); % each line of E contains the approximation error for a specific method
for i = 1 : length(deltat)-1
    t_deltat = 0 : deltat(i) : tend;
    y_euler = euler_meth(@func_p, p0, deltat(i), tend); % Euler method
    y_heun = heun_meth(@func_p, p0, deltat(i), tend); % Heun method
    y_RK4 = runge_kutta_4_meth(@func_p, p0, deltat(i), tend); % Runge-Kutta method
    E_tilde(:,i) = sqrt((deltat(i)/tend)* ...
            ([sum((y_euler-y_euler_best(1:deltat(i)/deltat_best:end)).^2); ... % select only the corresponding elements
              sum((y_heun-y_heun_best(1:deltat(i)/deltat_best:end)).^2); ...
              sum((y_RK4-y_RK4_best(1:deltat(i)/deltat_best:end)).^2)]));
end
%% Display tables
format short e; % display numbers in the short engineering form

T1 = table([deltat; E(1,:); E_red(1,:); E_tilde(1,:)],...
    'VariableNames',{'explicit_Euler_method_q_1'},...
    'RowNames',{'deltat' 'error' 'error red.' 'error app'});
T1.Properties.Description='explicit Euler method (q=1)';

T2 = table([deltat; E(2,:); E_red(2,:); E_tilde(2,:)],...
    'VariableNames',{'method_of_Heun_q_2'},...
    'RowNames',{'deltat' 'error' 'error red.' 'error app'});
T2.Properties.Description='method of Heun (q=2)';

T3 = table([deltat; E(3,:); E_red(3,:); E_tilde(3,:)],...
    'VariableNames',{'Runge_Kutta_method_q_4'},...
    'RowNames',{'deltat$' 'error' 'error red.' 'error app'});
T3.Properties.Description='Runge-Kutta method (q=4)';

sprintf(['Table 1: ', T1.Properties.Description]);
disp(T1)
sprintf(['Table 2: ', T2.Properties.Description])
disp(T2)
sprintf(['Table 3: ', T3.Properties.Description])
disp(T3)

format % reset the format
%% Functions
function p = sol_p (t)
    % Analytical solution
    p = 10 ./ (1 + 9 * exp(- t));
end

function dp = func_p(t, p)
    % defining the ODE to be solved numerically
    dp = (1 - p / 10) * p;
end