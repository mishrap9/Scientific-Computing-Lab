classdef Utilities
    methods(Static)
        function count_prev_figs = Open_Figures(fig_names, title_names, t, p, fig_prop, count_prev_figs)
            % Open_Figures is used for opening figures and managing their
            % properties
            for i = 1 : length(fig_names)
                count_prev_figs = count_prev_figs + 1;
                Fig = figure(count_prev_figs);
                set(Fig, 'NumberTitle', 'off', ...
                  'Name', fig_names{i}, ...
                  'Units', 'normalized', ...
                  'OuterPosition', [0 0 1 1]);
                ploti = plot(t, p, 'r');
                hold on;
                grid on;
                ploti.LineWidth = fig_prop.LineWidthAnalytical;
                xlab = xlabel('t', 'FontSize', fig_prop.LabelSize);
                set(xlab, 'Interpreter', 'latex');
                ylab = ylabel('p(t)', 'FontSize', fig_prop.LabelSize);
                set(ylab, 'Interpreter', 'latex');
                if i == 1
                    title(title_names{i}, 'interpreter', 'latex', 'FontSize', fig_prop.TitleFontSize);
                    leg = legend('$p(t)$');
                    set(leg, 'Interpreter', 'latex');
                    leg.FontSize = fig_prop.LegFontSize;
                    leg.Location = 'southeast';
                else
                    title([title_names{i}, ' with time steps $\delta (t)=\{\frac{1}{2},\frac{1}{4},\frac{1}{8},\frac{1}{16},\frac{1}{32}\}$'], ...
                      'interpreter', 'latex', 'FontSize', fig_prop.TitleFontSize);
                end
                xlim([0, 5]);
                ylim([0, 20]);
            end
        end
     
        function make_plot(fig_nr, t, y, marker, color, fig_prop)
            % make_plot is used for applying different properties to plots
            figure(fig_nr);
            p = plot(t, y, marker);
            p.MarkerSize = fig_prop.MarkerSize;
            p.Color = color;
            p.LineStyle = fig_prop.LineStyle;
        end
     
        function set_legends(figs_nr, leg_text, fig_prop, leg_entries)
            % set_legends manage the legends of the figures
            %   The input leg_entries contains the logical assignments for
            %   the cases that can be ploted
            if nargin < 4
                leg_entries = ones(length(figs_nr), size(leg_text,2)-1);
            end
            leg_entries = [ones(length(figs_nr),1), leg_entries]; % 1st entry is the analytical solution
            leg_entries = logical(leg_entries);
            for i = 1 : length(figs_nr)
                figure(figs_nr(i));
                leg = legend(leg_text(leg_entries(i,:)));
                set(leg, 'Interpreter', 'latex');
                leg.FontSize = fig_prop.LegFontSize;
                leg.Location = 'southeast';
                xlim([0, 5]);
                ylim([0, 20]);
            end
        end
        
    end
end