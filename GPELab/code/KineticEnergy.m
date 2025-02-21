function KE = KineticEnergy(Psi,Delta,Geometry2D)
    Psi1 = Psi{1};
    Psi2 = Psi{2};  

    dx = Geometry2D.X(1,2) - Geometry2D.X(1,1);  % Step size in the x-direction
    dy = Geometry2D.Y(2,1) - Geometry2D.Y(1,1);  % Step size in the y-direction


    % Compute gradients of Psi1 and Psi2
    [gradPsi1_x, gradPsi1_y] = gradient(Psi1, dx, dy);
    [gradPsi2_x, gradPsi2_y] = gradient(Psi2, dx, dy);
    
    % Kinetic energy density
    KE_density = abs(gradPsi1_x).^2 + abs(gradPsi1_y).^2 + abs(gradPsi2_x).^2 + abs(gradPsi2_y).^2;
    
    % Perform integration along x, then y using trapz
    KE = Delta * trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), KE_density, 2), 1);
end