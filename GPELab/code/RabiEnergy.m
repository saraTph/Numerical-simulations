function [E_rabi] = RabiEnergy(Psi,delta,Omega,Geometry2D)
    Psi1 = Psi{1};
    Psi2 = Psi{2};
    
    % Compute rabi energy density
    Rabi_density = (delta * (abs(Psi1).^2 - abs(Psi2).^2) + Omega * (conj(Psi1).*Psi2 + conj(Psi2).*Psi1));

    % Perform 2D integration using trapz
    E_rabi = 0.5 * trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), Rabi_density, 2), 1);
end