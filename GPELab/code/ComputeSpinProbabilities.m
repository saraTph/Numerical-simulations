function [P_up, P_down] = ComputeSpinProbabilities(Psi, Geometry2D)
    % Compute probability densities
    Psi_1_density = abs(Psi{1}).^2; % Spin-up density
    Psi_2_density = abs(Psi{2}).^2; % Spin-down density
    
    % Compute the integrals
    N_up = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), Psi_1_density, 2)); % ∫ |Psi_1|^2 dA
    N_down = trapz(Geometry2D.Y(:,1), trapz(Geometry2D.X(1,:), Psi_2_density, 2)); % ∫ |Psi_2|^2 dA
    N_total = N_up + N_down; % Total norm ∫ (|Psi_1|^2 + |Psi_2|^2) dA
    
    % Compute probabilities
    P_up = N_up / N_total;   % Probability of being in the up state
    P_down = N_down / N_total; % Probability of being in the down state

    % Handle floating-point errors
    if abs(P_up) < 1e-10  
        P_up = 0;
    end
    if abs(P_down) < 1e-10  
        P_down = 0;
    end
end
