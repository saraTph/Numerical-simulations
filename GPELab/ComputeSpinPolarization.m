%% Computation of spin mixture polarization P
%% P=1: Fully spin-up state
%% P=0: Equal superposition of spin states
%% P=-1: Fully spin-down state

function P_total = ComputeSpinPolarization(Psi, Geometry2D)
    % Compute probability densities
    Psi_1_density = abs(Psi{1}).^2;
    Psi_2_density = abs(Psi{2}).^2;
    
    % Compute numerator and denominator for total spin polarization
    %Y_vals = Geometry2D.Y(:,1);
    num = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), Psi_1_density - Psi_2_density, 2)); % ∫(|Psi_1|^2 - |Psi_2|^2) dA
    den = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), Psi_1_density + Psi_2_density, 2)); % ∫(|Psi_1|^2 + |Psi_2|^2) dA
    
    % Compute total spin polarization
    P_total = num / den;
end