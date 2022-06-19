classdef Utilities
    % This class implements functions for creating the figures and the
    % table.
    methods(Static)
        function Create_Figures_Explicit_Euler(t)
        % Creates the windows that will contain the plots for the
        % implementation with the Explicit Euler Method.
            fig_names = {'EE; t = 1/8', 'EE; t = 2/8', ...
                         'EE; t = 3/8', 'EE; t = 4/8'};
            for i = 1 : length(t)
                Fig = figure(i);
                set(Fig, 'NumberTitle', 'off', ...
                  'Name', fig_names{i}, ...
                  'Units', 'normalized', ...
                  'OuterPosition', [0 0 1 1]);
            end
        end
     
        function Plot_sol_EE(T, Nx_vec, Ny_vec, dt, t, i, j, k, save_subplots_EE)
        % Plots the solution T in the corresponding subplot and save it
        % as a .PNG image.
            TitleSize = 10;
            figure(k);
            sb = subplot(length(Nx_vec), length(dt), ...
                         length(dt) * (i - 1) + j); % current subplot
            [X, Y] = meshgrid((0:Nx_vec(i) + 1) / (Nx_vec(i) + 1), ...
                              (0:Ny_vec(i) + 1) / (Ny_vec(i) + 1));
            h = mesh(X, Y, T);
            h.LineWidth = 1;
            axis tight
            Nxtitle = sprintf('$N_x=N_y=%d; ', Nx_vec(i));
            dttitle = sprintf('=1/%d$', 1 / dt(j));
            title([Nxtitle, '\delta t', dttitle], ...
                  'interpreter', 'latex', ...
                  'FontSize', TitleSize);
            % reduce a bit the white space arround the subplots
            outerpos = sb.OuterPosition;
            ti = sb.TightInset;
            left = outerpos(1) + ti(1);
            bottom = outerpos(2) + ti(2);
            sx_width = outerpos(3) - ti(1) - ti(3);
            sb_height = outerpos(4) - ti(2) - ti(4);
            sb.Position = [left bottom sx_width sb_height];
            % save image
            if save_subplots_EE
                LabelSize = 16;
                TitleSize = 18;
                hfig = figure(10);
                h = mesh(X, Y, T);
                h.LineWidth = 1;
                colorbar;
                axis tight
                Nxtitle = sprintf('N_x=N_y=%d; ', Nx_vec(i));
                dttitle = sprintf('=1/%d$', 1 / dt(j));
                ttitle = sprintf('$t=%.3f; ', t(k));
                title([ttitle, Nxtitle, '\delta t', dttitle], ...
                      'interpreter', 'latex', ...
                      'FontSize', TitleSize);
                xlab = xlabel('x', 'FontSize', LabelSize);
                set(xlab, 'Interpreter', 'latex');
                ylab = ylabel('y', 'FontSize', LabelSize);
                set(ylab, 'Interpreter', 'latex');
                zlabel_txt = sprintf('$T(x,y,t=%d/8)$', 8 * t(k));
                zlab = zlabel(zlabel_txt, 'FontSize', LabelSize);
                set(zlab, 'Interpreter', 'latex');
                img_name = sprintf('t=%d%%8_Nx=%d_dt=1%%%d', 8 * t(k), Nx_vec(i), 1 / dt(j));
                cd pics;
                saveas(hfig, fullfile(cd, 'Explicit_Euler_imgs', [img_name, '.png']));
                cd ..
            end
        end
     
        function Edit_Figs_EE(t)
        % Adds Titles to the windows with the solutions from the Explicit 
        % Euler Method implementation and save them as .PNG images.
            % Close the Figure for saving the subplots (figure(10)).
            figure(10);
            close;
            TitleFontSize = 26;
            for k = 1 : length(t)
                hfig = figure(k);
                title_text = sprintf('Temperatures $T(x,y,t=%d/8)$  (Explicit Euler)', 8 * t(k));
                % Add Title if Matlab 2018b or Bioinformatics Toolbox is„
                % installed
                if strcmp(version('-release'), '2018b') % check version
                    gt = sgtitle(title_text, 'Interpreter', 'latex');
                    gt.FontSize = TitleFontSize;
                else
                    if license('test', 'Bioinformatics_Toolbox')
                        suptitle(title_text); % Bioinformatics Toolbox
                    end
                end
                img_window_title = sprintf('EE_t=%d%%8', 8 * t(k));
                cd pics
                saveas(hfig, [img_window_title, '.png']);
                cd ..
            end
        end
     
        function Table_stable_cases = Create_Table_err(stable_cases)
        % Creates the table that shows the stable cases.
        %   0 - unstable; 1 - stable
            Nx_Ny = [3, 7, 15, 31]';
            dt_64 = stable_cases(:, 1);
            dt_128 = stable_cases(:, 2);
            dt_256 = stable_cases(:, 3);
            dt_512 = stable_cases(:, 4);
            dt_1024 = stable_cases(:, 5);
            dt_2048 = stable_cases(:, 6);
            dt_4096 = stable_cases(:, 7);
            Table_stable_cases = table(Nx_Ny, ...
                                       dt_64, dt_128, dt_256, dt_512, ...
                                       dt_1024, dt_2048, dt_4096);
            Table_stable_cases.Properties.Description = 'Table Stability';
        end
     
        function Fig = Create_Figure_Implicit_Euler()
        % Creates the window that will contain the plots from the
        % implementation with the Implicit Euler Method.
            fig_name = 'IE; dt = 1/64';
            Fig = figure(5);
            set(Fig, 'NumberTitle', 'off', ...
                'Name', fig_name, ...
                'Units', 'normalized', ...
                'OuterPosition', [0 0 1 1]);
        end
     
        function Plot_sol_IE(T, Fig, Nx_vec, Ny_vec, t, i, k)
            % Plots the solution T in the corresponding subplot.
            TitleSize = 13;
            LabelSize = 10;
            figure(Fig);
            subplot(length(Nx_vec), length(t), length(t) * (i - 1) + k);
            [X, Y] = meshgrid((0:Nx_vec(i) + 1) / (Nx_vec(i) + 1), ...
                              (0:Ny_vec(i) + 1) / (Ny_vec(i) + 1));
            mesh(X, Y, T);
            Nxtitle = sprintf('$N_x=N_y=%d;', Nx_vec(i));
            ttitle = sprintf(' t=%d/8$', 8 * t(k));
            title([Nxtitle, ttitle], ...
                  'interpreter', 'latex', ...
                  'FontSize', TitleSize);
            zlabel_txt = sprintf('$T(x,y,t=%d/8)$', 8 * t(k));
            zlab = zlabel(zlabel_txt, 'FontSize', LabelSize);
            set(zlab, 'Interpreter', 'latex');
        end
     
        function Edit_Figs_IE(Fig)
        % Add Titles to the window with the solutions from the
        % implementation with the Implicit Euler Method and save it as
        % a .PNG image.
            figure(Fig);
            TitleFontSize = 22;
            title_text = 'Temperatures $T(x,y,t)$ with time step $\delta t=1/64$ (Implicit Euler)';
            % Add Title if Matlab 2018b or Bioinformatics Toolbox
            % installed
            if strcmp(version('-release'), '2018b') % check version
                gt = sgtitle(title_text, 'Interpreter', 'latex');
                gt.FontSize = TitleFontSize;
            else
                if license('test', 'Bioinformatics_Toolbox')
                    suptitle(title_text); % Bioinformatics Toolbox
                end
            end
            img_window_title = sprintf('IE_dt=1%%64');
            cd pics
            saveas(Fig, [img_window_title, '.png']);
            cd ..
        end
     
    end
end