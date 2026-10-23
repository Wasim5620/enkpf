function [fig1, fig2] = plot_confidence_intervals_and_data(P, PenkfP, z2, x_hat, x_hatenkfP, Date,Date1,Date2, ~)
    % Calculate confidence intervals
    n_sta_var = size(P, 1);
    n_mes = size(z2, 2);
    n_forecast = size(P, 3) - n_mes;
    
    % Extract the variances from the diagonal of the P matrices
    variances1 = zeros(n_sta_var, n_mes + n_forecast);
     variances2 = zeros(n_sta_var, n_mes + n_forecast);
    for t = 1:(n_mes + n_forecast)
        variances1(:, t) = diag(P(:, :, t));
        variances2(:, t) = diag(PenkfP(:, :, t));
    end
    
    % Calculate the standard deviations
    std_devs1 = sqrt(variances1);
    std_devs2 = sqrt(variances2);
    % Calculate the confidence intervals (95%)
    conf_lower_intr = x_hat - 1.96 * std_devs1;
    conf_upper_intr = x_hat + 1.96 * std_devs1;
    
    % Similarly, you can calculate conf_lower_intr_enkfP and conf_upper_intr_enkfP
    % Assuming you have the std_devs for x_hatenkfP (not provided in the original code)
    % Here, you should calculate or use the appropriate variance matrix for x_hatenkfP
    conf_lower_intr_enkfP = x_hatenkfP - 1.96 * std_devs2; % Placeholder
    conf_upper_intr_enkfP = x_hatenkfP + 1.96 * std_devs2; % Placeholder

    % Plotting (first figure)
    fig1 = figure();
    set(fig1, 'Position', [360 98 1300 550]);
     % pos1 = [0.085 0.27 0.42 0.7];
     % pos2 = [0.575, 0.28, 0.45, 0.7];
     % ax1 = axes('Position', pos1);
    
    % Initialize plot to get colors
    pl = plot(1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10, 1:10);
    c = get(pl, 'Color');
    delete(pl); % Remove temporary plot
    tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
    nexttile;
    % Plot confidence intervals and data
    ciplot(conf_lower_intr(4,:), conf_upper_intr(4,:), Date, c{1}, 0.20);
    hold on;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5, 'YColor', c{1},'YGrid','on');
    plot(Date, z2(1,:), 'Color', c{2}, 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 0.5);
    plot(Date, x_hat(4,:), 'Color', c{1}, 'LineWidth', 1.5, 'LineStyle', '--');
    xlabel('Time (days)', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    ylabel('Active Cases, Q', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    leg = legend('Confidence Interval EnKF', 'True', 'Estimated value', 'Location', 'best', 'fontsize',10);
    set(leg, 'Interpreter', 'tex');
    set(gca, 'TickLabelInterpreter', 'tex');
    
    
    xline(Date1(end), '--k', '2nd wave start','HandleVisibility', 'off');
    xline(Date2(end), '--k', '3rd wave start','HandleVisibility', 'off');
    axis tight;
    % Zoomed-in plot
    % ax3 = axes('Position', [0.355 0.22 0.12 0.24]);
    % box on; % Put box around new pair of axes
    % indexOfInterest = (T < 40) & (T > 30); % Range of t for zoom
    % ciplot(conf_lower_intr(4, indexOfInterest), conf_upper_intr(4, indexOfInterest), Date(indexOfInterest), c{1}, 0.16);
    % hold on;
    % plot(ax3, Date(indexOfInterest), z2(1, indexOfInterest), 'Color', c{2}, 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 0.5);
    % plot(ax3, Date(indexOfInterest), x_hat(4, indexOfInterest), 'Color', c{1}, 'LineWidth', 1.5, 'LineStyle', '--');
    % axis tight;
    % 
    % Second plot
    % ax2 = axes('Position', pos2);
    nexttile;
    ciplot(conf_lower_intr_enkfP(4,:), conf_upper_intr_enkfP(4,:), Date, c{4}, 0.16);
    hold on;
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5, 'YColor', c{4},'YGrid','on');
    plot(Date, z2(1,:), 'Color', c{2}, 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 0.5);
    plot(Date, x_hatenkfP(4,:), 'Color', c{4}, 'LineWidth', 1.5, 'LineStyle', '--');
    xlabel('Time (days)', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    ylabel('Active Cases, Q', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    leg = legend('Confidence Interval EnKPF', 'True', 'Estimated value', 'Location', 'best','fontsize',10);
    set(leg, 'Interpreter', 'tex');
    set(gca, 'TickLabelInterpreter', 'tex');
    
    xline(Date1(end), '--k', '2nd wave start','HandleVisibility', 'off');
    xline(Date2(end), '--k', '3rd wave start','HandleVisibility', 'off');
    axis tight;
    
    % Zoomed-in plot for the second plot
    % ax4 = axes('Position', [0.855 0.22 0.12 0.24]);
    % box on; % Put box around new pair of axes
    %  % indexOfInterest = (T < 10) & (T > 1); % Range of t for zoom
    % ciplot(conf_lower_intr_enkfP(4, indexOfInterest), conf_upper_intr_enkfP(4, indexOfInterest), Date(indexOfInterest), c{4}, 0.16);
    % hold on;
    % plot(ax4, Date(indexOfInterest), z2(1, indexOfInterest), 'Color', c{2}, 'LineWidth', 2, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 0.5);
    % plot(ax4, Date(indexOfInterest), x_hatenkfP(4, indexOfInterest), 'Color', c{4}, 'LineWidth', 1.5, 'LineStyle', '-.');
    % axis tight;
    
   
    % tightfig(gcf)
     hold off;

    %% Second Figure
    fig2 = figure();
    set(fig2, 'Position', [360 200 1300 550]);
    % axes('Position', pos1);
     tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
    nexttile;
    ciplot(conf_lower_intr(5,:), conf_upper_intr(5,:), Date, c{1}, 0.20);
    hold on;
    set(gca, 'FontSize', 11, 'FontWeight', 'bold', 'LineWidth', 1.5, 'YColor', c{1},'YGrid','on');
    plot(Date, z2(2,:), 'Color', c{2}, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 0.5);
    plot(Date, x_hat(5,:), 'Color', c{1}, 'LineWidth', 1.5, 'LineStyle', '--');
    xlabel('Time (days)', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    ylabel('Recovered Cases, R', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    leg = legend('Confidence Interval EnKF', 'True', 'Estimated values', 'Location', 'best','fontsize',10);
    set(leg, 'Interpreter', 'tex');
    set(gca, 'TickLabelInterpreter', 'tex');
     
    xline(Date1(end), '--k', '2nd wave start','HandleVisibility', 'off');
    xline(Date2(end), '--k', '3rd wave start','HandleVisibility', 'off');
    axis tight;
    
    % Zoomed-in plot for Recovered Cases, R
    % axes('Position', [0.355 0.445 0.12 0.2]);
    % box on; % Put box around new pair of axes
    % ciplot(conf_lower_intr(5, indexOfInterest), conf_upper_intr(5, indexOfInterest), Date(indexOfInterest), c{1}, 0.16);
    % hold on;
    % plot(Date(indexOfInterest), z2(2, indexOfInterest), 'Color', c{2}, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 0.5);
    % plot(Date(indexOfInterest), x_hat(5, indexOfInterest), 'Color', c{1}, 'LineWidth', 1.5, 'LineStyle', '--');
    % axis tight;
    
    % Second plot for Recovered Cases, R
    % axes('Position', pos2);
    nexttile;
    ciplot(conf_lower_intr_enkfP(5,:), conf_upper_intr_enkfP(5,:), Date, c{4}, 0.16);
    hold on;
    set(gca, 'FontSize', 11, 'FontWeight', 'bold', 'LineWidth', 1.5, 'YColor', c{4},'YGrid','on');
    plot(Date, z2(2,:), 'Color', c{2}, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 0.5);
    plot(Date, x_hatenkfP(5,:), 'Color', c{4}, 'LineWidth', 1.5, 'LineStyle', '--');
    xlabel('Time (days)', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    ylabel('Recovered Cases, R', 'Interpreter', 'tex', 'FontWeight', 'bold','FontSize',16);
    leg = legend('Confidence Interval EnKPF', 'True', 'Estimated value', 'Location', 'best','fontsize',10);
    set(leg, 'Interpreter', 'tex');
    set(gca, 'TickLabelInterpreter', 'tex');
     
    xline(Date1(end), '--k', '2nd wave start','HandleVisibility', 'off');
    xline(Date2(end), '--k', '3rd wave start','HandleVisibility', 'off');
    axis tight;
    
    % % Zoomed-in plot for Recovered Cases, R
    % axes('Position', [0.855 0.445 0.12 0.2]);
    % box on; % Put box around new pair of axes
    % % indexOfInterest = (T < 20) & (T > 1); % Range of t for zoom
    % ciplot(conf_lower_intr_enkfP(5, indexOfInterest), conf_upper_intr_enkfP(5, indexOfInterest), Date(indexOfInterest), c{4}, 0.16);
    % hold on;
    % plot(Date(indexOfInterest), z2(2, indexOfInterest), 'Color', c{2}, 'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 0.5);
    % plot(Date(indexOfInterest), x_hatenkfP(5, indexOfInterest), 'Color', c{4}, 'LineWidth', 1.5, 'LineStyle', '-.');
    % axis tight;
    
    % tightfig(gcf); % Adjust figure layout
    hold off;
    % Return figure handles fig1 and fig2
end


