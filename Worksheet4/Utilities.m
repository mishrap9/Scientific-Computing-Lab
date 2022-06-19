classdef Utilities
    % This class implements function for creating the figures and tables.
    methods(Static)
        function Create_Figures(Nx_vec, Ny_vec, T_cell)
            % Create the figures with coloured surface and contour plot
            % representing the temperatures.
            LabelSize = 18;
            TitleFontSize = 22;
            SubTitleFontSize = 20;
            fig_names = {'Nx = Ny = 3', 'Nx = Ny = 7', 'Nx = Ny = 15', ...
                         'Nx = Ny = 31', 'Nx = Ny = 63'};
            subplot_titles = {'Coloured surface', 'Contour plot'};
            title_name = 'Numerical solution $T(x,y)$ for $N_x = N_y =$ %d';
            for i = 1 : length(Nx_vec)
                Fig = figure(i);
                set(Fig, 'NumberTitle', 'off', ...
                  'Name', fig_names{i}, ...
                  'Units', 'normalized', ...
                  'OuterPosition', [0 0 1 1]);
             
                % 1) coloured surface
                subplot(1, 2, 1);
             
                [X, Y] = meshgrid((0:Nx_vec(i) + 1) / (Nx_vec(i) + 1), ...
                                  (0:Ny_vec(i) + 1) / (Ny_vec(i) + 1));
                h1 = surf(X, Y, T_cell{i});
                h1.EdgeColor = 'none';
                h1.FaceColor = 'texturemap';
                colormap('hsv')
                colorbar;
                axis tight
             
                xlab = xlabel('x', 'FontSize', LabelSize);
                set(xlab, 'Interpreter', 'latex');
                ylab = ylabel('y', 'FontSize', LabelSize);
                set(ylab, 'Interpreter', 'latex');
                zlab = zlabel('$T(x,y)$', 'FontSize', LabelSize);
                set(zlab, 'Interpreter', 'latex');
             
                % 2) contour plot
                subplot(1, 2, 2);
             
                x_contour = (0:Nx_vec(i) + 1) / (Nx_vec(i) + 1);
                y_contour = (0:Ny_vec(i) + 1) / (Ny_vec(i) + 1);
                [~, h2] = contour(x_contour, y_contour, T_cell{i});
                h2.LevelList = (0:Nx_vec(i) + 1) / (Nx_vec(i) + 1);
                h2.LineColor = 'auto';
                h2.LineWidth = 2.5;
                h2.LevelStep = 1 / (Nx_vec(i) + 1);
                colorbar
                axis tight
             
                xlab = xlabel('x', 'FontSize', LabelSize);
                set(xlab, 'Interpreter', 'latex');
                ylab = ylabel('y', 'FontSize', LabelSize);
                set(ylab, 'Interpreter', 'latex');
                check_v = strcmp(version('-release'), '2018b'); % check version
                if check_v || license('test', 'Bioinformatics_Toolbox')
                    xtitle = sprintf(title_name, Nx_vec(i));
                    subplot(1, 2, 1);
                    title(subplot_titles{1}, ...
                      'interpreter', 'latex', 'FontSize', SubTitleFontSize);
                    subplot(1, 2, 2);
                    title(subplot_titles{2}, ...
                      'interpreter', 'latex', 'FontSize', SubTitleFontSize);
                    if check_v
                        sgtitle(xtitle, ...
                          'interpreter', 'latex', 'FontSize', TitleFontSize);
                    else
                        suptitle(xtitle, ...
                          'interpreter', 'latex', 'FontSize', TitleFontSize); % Bioinformatics Toolbox
                    end
                else
                    subplot(1, 2, 1);
                    xtitle = sprintf([subplot_titles{1}, ' of $T(x,y)$ for $N_x = N_y = %d$'], Nx_vec(i));
                    title(xtitle, ...
                      'interpreter', 'latex', 'FontSize', SubTitleFontSize);
                    subplot(1, 2, 2);
                    xtitle = sprintf([subplot_titles{2}, ' of $T(x,y)$ for $N_x = N_y = %d$'], Nx_vec(i));
                    title(xtitle, ...
                      'interpreter', 'latex', 'FontSize', SubTitleFontSize);
                end
            end
        end
     
        function [TableFullMatrix, TableSparseMatrix, TableGS] = ...
                        Create_Tables_performance(Nx_vec, runtime, storage)
            % This function creates the tables with runtime and storage
            TableFullMatrix = table([Nx_vec(2:end); runtime(1, 2:end); storage(1, 2:end)], ...
                'VariableNames', {'Direct_solution_with_full_matrix'}, ...
                'RowNames', {'Nx, Ny', 'runtime [sec]', 'storage [doubles]'});
            TableFullMatrix.Properties.Description = 'Direct solution with full matrix';
         
            TableSparseMatrix = table([Nx_vec(2:end); runtime(2, 2:end); storage(2, 2:end)], ...
                'VariableNames', {'Direct_solution_with_sparse_matrix'}, ...
                'RowNames', {'Nx, Ny', 'runtime [sec]', 'storage [doubles]'});
            TableSparseMatrix.Properties.Description = 'Direct solution with sparse matrix';
         
            TableGS = table([Nx_vec(2:end); runtime(3, 2:end); storage(3, 2:end)], ...
                'VariableNames', {'Iterative_solution_with_Gauss_Seidel'}, ...
                'RowNames', {'Nx, Ny', 'runtime [sec]', 'storage [doubles]'});
            TableGS.Properties.Description = 'Iterative solution with Gauss-Seidel';
        end
     
        function Table_err = Create_Table_err(Nx_vec, err, err_red)
            % This function creates the table with the errors between the
            % analytical solution and the iterative solution obtained using
            % the Gauss-Seidel method.
            Table_err = table([Nx_vec; err; err_red], ...
                           'VariableNames', {'Table_Errors'}, ...
                           'RowNames', {'Nx=Ny', 'error', 'error red.'});
            Table_err.Properties.Description = 'Table Errors';
        end
     
    end
end