function [en] = InteractionEnergy(Psi,g11,g22,g12,Geometry2D)
    Psi1 = Psi{1};
    Psi2 = Psi{2};
    
    % Compute interaction energy density
    interaction_density = (0.5 * g11 * abs(Psi1).^4 + 0.5 * g22 * abs(Psi2).^4 + g12 * abs(Psi1).^2 .* abs(Psi2).^2 );

    % Perform 2D integration using trapz
    en = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), interaction_density, 2), 1);
end