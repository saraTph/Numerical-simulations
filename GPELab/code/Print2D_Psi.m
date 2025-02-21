function [Psi_1_density, Psi_2_density] = Print2D_Psi(Psi, Geometry2D)
    % Compute probability densities
    Psi_1_density = abs(Psi{1}).^2;
    Psi_2_density = abs(Psi{2}).^2;

    % Plot for the first component
    figure;
    set(gcf, 'Position', [100, 100, 800, 300]); % Set figure size
    subplot(1,2,1);
    imagesc(Geometry2D.X(1,:), Geometry2D.Y(:,1), Psi_1_density);
    colorbar;
    xlabel('X'); ylabel('Y'); title('|Psi_1|^2');
    set(gca, 'YDir', 'normal');

    % Plot for the second component
    subplot(1,2,2);
    imagesc(Geometry2D.X(1,:), Geometry2D.Y(:,1), Psi_2_density);
    colorbar;
    xlabel('X'); ylabel('Y'); title('|Psi_2|^2');
    set(gca, 'YDir', 'normal');
end

