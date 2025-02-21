function [Density_1D_1, Density_1D_2] = Plot1D_DensityProfile(Psi, Geometry2D, plot_title)
    % Compute probability densities in 2D
    Psi_1_density = abs(Psi{1}).^2;
    Psi_2_density = abs(Psi{2}).^2;
    
    % Integrate along the Y-direction to get 1D profile along X
    Y_vals = Geometry2D.Y(:,1);

    Density_1D_1 = trapz(Y_vals, Psi_1_density, 1); % Integrate along Y
    Density_1D_2 = trapz(Y_vals, Psi_2_density, 1); % Integrate along Y
    
    % Create figure and set size
    figure;
    set(gcf, 'Position', [100, 100, 800, 500]); % Set figure size

    % Plot 1D density profiles
    plot(Geometry2D.X(1,:), Density_1D_1, 'b-', 'LineWidth', 2); hold on;
    plot(Geometry2D.X(1,:), Density_1D_2, 'r-', 'LineWidth', 2);
    hold off;
    
    % Customize plot
    xlabel('X');
    ylabel('Integrated Density');
    legend('|Psi_1|^2', '|Psi_2|^2');
    title(plot_title);
    grid on;
end