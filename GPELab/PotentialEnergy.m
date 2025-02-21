function PE = PotentialEnergy(Psi,gamma_x,gamma_y,Geometry2D)

    Psi1 = Psi{1};
    Psi2 = Psi{2};

    % Compute potential energy densities for each component
    V = (1/2)*(gamma_x*Geometry2D.X.^2+gamma_y*Geometry2D.Y.^2);  % External potential
    
    % Compute total potential energy density
    PE_density = V .* abs(Psi1).^2 + V .* abs(Psi2).^2;

    % Perform integration along x, then y using trapz
    PE = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), PE_density, 2), 1);

end