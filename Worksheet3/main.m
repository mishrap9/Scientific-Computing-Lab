% Lab Course Scientific Computing - Worksheet 3, Group 16
clear; clc; close all;
%% Figure properties constants
fig_prop.LabelSize = 18;
fig_prop.LegFontSize = 18; % legend font size
fig_prop.TitleFontSize = 20;
fig_prop.LineWidthAnalytical = 1;
fig_prop.MarkerSize = 14;
fig_prop.LineStyle = '-.';

%% a) Analytical solution and
%  b) Approximation error for Explicit Methods and deltat={1/2,1/4,1/8,1/16,1/32}
% defining the time interval for computing the analytical solution
tstep = 1e-2;
tend = 5;
t = 0 : tstep : tend;
p = problem_functions.sol_p(t);

% setting figures properties
fig_names = {'Analytical solution', 'Explicit Euler', 'Heun'};
title_names = {'a) Graph of $p(t)=\frac{200}{20-10e^{-7t}}$', 'b) Explicit Euler', 'b) Heun'};
count_prev_figs = 0;
count_curr_figs = Utilities.Open_Figures(fig_names, title_names, t, p, fig_prop, count_prev_figs);

% Different markers with different colours
markers = {'.', '.', '.', '.', '.'};
colors = {[0.5 0.3 0], 'b', [0.5 0.1 1], 'm', [0.1 1 0.1]};

p0 = 20; % initial condition
deltat = [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32];
% matrix with stable cases
stable_cases = zeros(length(deltat), 6);
% each line of EE contains the approximation error for a specific explicit method
EE = zeros(2, length(deltat)); % Results for Explicit Euler and Heun
y_explicit = cell(2, length(deltat));
for i = 1 : length(deltat)
    t_deltas = 0 : deltat(i) : tend;
    p_deltas = problem_functions.sol_p(t_deltas); % the analytical solution with time step deltat(i)
 
    y_explicit{1,i} = euler_explicit(@problem_functions.func_p, p0, deltat(i), tend); % Explicit Euler method
    Utilities.make_plot(2, t_deltas, y_explicit{1,i}, markers{i}, colors{i}, fig_prop);
 
    y_explicit{2,i} = heun(@problem_functions.func_p, p0, deltat(i), tend); % Heun method
    Utilities.make_plot(3, t_deltas, y_explicit{2,i}, markers{i}, colors{i}, fig_prop);
 
    EE(:, i) = sqrt(deltat(i) / tend) * ...
               sqrt([sum((y_explicit{1,i} - p_deltas) .^ 2); ...
                     sum((y_explicit{2,i} - p_deltas) .^ 2)]);
end

leg_text = {'$p(t)$', '$\delta t=1/2$', '$\delta t=1/4$', '$\delta t=1/8$', '$\delta t=1/16$', '$\delta t=1/32$'};
figs_nr = [2, 3];
Utilities.set_legends(figs_nr, leg_text, fig_prop);

%% d),f),g) Approximation error for Implicit Methods and deltat={1/2,1/4,1/8,1/16,1/32}
% setting figures properties
fig_inames = {'Implicit Euler', 'Adams-Moulton', 'Adams-Moulton (l1)', 'Adams-Moulton (l2)'};
title_inames = {'d) Implicit Euler', 'd) Adams-Moulton', 'f) Adams-Moulton (l1)', 'f) Adams-Moulton (l2)'};
count_prev_figs = count_curr_figs;
count_curr_figs = Utilities.Open_Figures(fig_inames, title_inames, t, p, fig_prop, count_prev_figs);

% Different markers with different colours
markers = {'.', '.', '.', '.', '.'};
colors = {[0.5 0.3 0], 'b', [0.5 0.1 1], 'm', [0.1 1 0.1]};

p0 = 20; % initial condition
deltat = [1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 32];
% each line of EI contains the approximation error for a specific implicit method
EI = - ones(4, length(deltat)); % -1 means we don't have a solution
leg_entries = zeros(4,length(deltat)); % 0 - we don't have a solution; 1 - we plot the solution
y_implicit = cell(2, length(deltat));
for i = 1 : length(deltat)
    t_deltas = 0 : deltat(i) : tend;
    p_deltas = problem_functions.sol_p(t_deltas); % the analytical solution with time step deltat(i)
 
    y_implicit{1,i} = euler_implicit(@problem_functions.func_p, @problem_functions.dfunc_p, p0, deltat(i), tend); % Implicit Euler
    if (~isempty(y_implicit{1,i}))
        Utilities.make_plot(count_prev_figs + 1, t_deltas, y_implicit{1,i}, markers{i}, colors{i}, fig_prop);
        leg_entries(1,i) = 1;
        EI(1, i) = sqrt(deltat(i) / tend * sum((y_implicit{1,i} - p_deltas) .^ 2));
    end
 
    y_implicit{2,i} = adams_moulton(@problem_functions.func_p, @problem_functions.dfunc_p, p0, deltat(i), tend); % Adams-Moulton
    if (~isempty(y_implicit{2,i}))
        Utilities.make_plot(count_prev_figs + 2, t_deltas, y_implicit{2,i}, markers{i}, colors{i}, fig_prop);
        leg_entries(2,i) = 1;
        EI(2, i) = sqrt(deltat(i) / tend * sum((y_implicit{2,i} - p_deltas) .^ 2));
    end
 
    y_implicit{3,i} = adams_moulton_l1(p0, deltat(i), tend); % A-M linearisation 1
    if (~isempty(y_implicit{3,i}))
        Utilities.make_plot(count_prev_figs + 3, t_deltas, y_implicit{3,i}, markers{i}, colors{i}, fig_prop);
        leg_entries(3,i) = 1;
        EI(3, i) = sqrt(deltat(i) / tend * sum((y_implicit{3,i} - p_deltas) .^ 2));
    end
 
    y_implicit{4,i} = adams_moulton_l2(p0, deltat(i), tend); % A-M linearisation 2
    if (~isempty(y_implicit{4,i}))
        Utilities.make_plot(count_prev_figs + 4, t_deltas, y_implicit{4,i}, markers{i}, colors{i}, fig_prop);
        leg_entries(4,i) = 1;
        EI(4, i) = sqrt(deltat(i) / tend * sum((y_implicit{4,i} - p_deltas) .^ 2));
    end
end

leg_text = {'$p(t)$', '$\delta t=1/2$', '$\delta t=1/4$', '$\delta t=1/8$', '$\delta t=1/16$', '$\delta t=1/32$'};
figs_nr = count_prev_figs + (1:4);
Utilities.set_legends(figs_nr, leg_text, fig_prop, leg_entries);

%% h) Error reduced
% each line of EE_red contains the approximation errors for a specific explicit method
% each line of EI_red contains the approximation errors for a specific implicit method
EE_red = - ones(2, length(deltat));
EI_red = - ones(4, length(deltat));
for i = 1 : length(deltat) - 1
    EE_red(:, i + 1) = EE(:, i) ./ EE(:, i + 1);
    EI_red(:, i + 1) = EI(:, i) ./ EI(:, i + 1);
end
% the table positions without a meaning
EE_red(EE_red <= 0) = - 1;
EI_red(EI_red <= 0) = - 1;
%% g)&h) Display tables
% format short e; % display numbers in the short engineering format

TableEEuler = table([deltat; EE(1, :); EE_red(1, :)], ...
'VariableNames', {'explicit_Euler'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableEEuler.Properties.Description = 'explicit Euler';

TableHeun = table([deltat; EE(2, :); EE_red(2, :)], ...
  'VariableNames', {'Heun'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableHeun.Properties.Description = 'method of Heun';

TableIEuler = table([deltat; EI(1, :); EI_red(1, :)], ...
  'VariableNames', {'implicit_Euler'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableIEuler.Properties.Description = 'implicit Euler method';

TableAM = table([deltat; EI(2, :); EI_red(2, :)], ...
  'VariableNames', {'Adams_Moulton'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableAM.Properties.Description = 'Adams-Moulton method';

TableAM1 = table([deltat; EI(3, :); EI_red(3, :)], ...
  'VariableNames', {'Adams_Moulton_linearisation_1'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableAM1.Properties.Description = 'Adams-Moulton - linearisation 1';

TableAM2 = table([deltat; EI(4, :); EI_red(4, :)], ...
  'VariableNames', {'Adams_Moulton_linearisation_2'}, ...
  'RowNames', {'deltat' 'error' 'error red.'});
TableAM2.Properties.Description = 'Adams-Moulton - linearisation 2';

disp(TableEEuler)
disp(TableHeun)
disp(TableIEuler)
disp(TableAM)
disp(TableAM1)
disp(TableAM2)

format % reset the format
fprintf('Obs.: The value -1 means that the corresponding cell does not have meaning!\n\n');

%% i) Stability
% For each method, we choose the best obtained result, corresponding to the
% minimum deltat = 1/32; this results is stable
y_results = [y_explicit; y_implicit]'; % cell array containing the y's
table_stable_cases = ones(size(y_results));
eps_sol = 1e-1; % the accuracy between consecutive solution
for j = 1 : size(y_results, 2)
    % each set of results from a method is a column of y_results
    y_best = cell2mat(y_results(end, j));
    y_best_end = y_best(end);
    for i = 1 : size(y_results, 1)
        y_curr = cell2mat(y_results(i, j));
        if isempty(y_curr)
            table_stable_cases(i, j) = 0;
        else
            if (abs(y_curr(end)) > 1 / eps)
                table_stable_cases(i, j) = 0;
            else
                % if the distances between our points and the asympthotic
                % correct solution increase in one of the two last points, and
                % we are not near the solution => the computed solution is
                % unstable
                ybias = abs(y_curr - y_best_end);
                increasing = (ybias(2:end) - ybias(1:end - 1)) > 0;
                if ybias(end) > eps_sol && ...
                    (increasing(end) || increasing(end - 1))
                    table_stable_cases(i, j) = 0;
                end
            end
        end
    end
end

delta_t = {'1/2'; '1/4'; '1/8'; '1/16'; '1/32'};
explicit_Euler = table_stable_cases(:, 1);
Heun = table_stable_cases(:, 2);
implicit_Euler = table_stable_cases(:, 3);
Adams_Moulton = table_stable_cases(:, 4);
Adams_Moulton_l1 = table_stable_cases(:, 5);
Adams_Moulton_l2 = table_stable_cases(:, 6);
TableStability = table(delta_t, explicit_Euler, Heun, implicit_Euler, Adams_Moulton, Adams_Moulton_l1, Adams_Moulton_l2);
TableStability.Properties.Description = 'Stable cases';
fprintf('\nTable with Stable cases: \n\n');
disp(TableStability)
