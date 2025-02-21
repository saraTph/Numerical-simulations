% author: Sara Tiengo 
% date: 30/01/2025

%% GROUND STATE COMPUTATION OF A MULTICOMPONENT BEC WITH COUPLED NONLINEARITIES

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP'; % method to solve the Continuous Normalized Gradient Flow (CNGF)
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-6}; % epsilon
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -10; % computational domine dimensions
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^7+1; % # grid points
Ny = 2^7+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);
Deltax = (xmax - xmin)/Nx;

%%

% scans
delta_values = linspace(3,-3, 15);
%delta_values = 0.32;
N_values = linspace(150e3, 150e6, 5);
%N_values = [150e3 150e5];

% initialize output vectors
P = zeros(length(N_values),length(delta_values));
P_up = zeros(length(N_values),length(delta_values));
P_down = zeros(length(N_values),length(delta_values));

energy_1 = zeros(length(N_values),length(delta_values));
energy_2 = zeros(length(N_values),length(delta_values));
energy_tot = zeros(length(N_values),length(delta_values));

chemPot_1 = zeros(length(N_values),length(delta_values));
chemPot_2 = zeros(length(N_values),length(delta_values));
chemPot_tot = zeros(length(N_values),length(delta_values));

IE = zeros(length(N_values),length(delta_values));
RE = zeros(length(N_values),length(delta_values));
PE = zeros(length(N_values),length(delta_values));
KE = zeros(length(N_values),length(delta_values));


size_1 = zeros(length(N_values),length(delta_values));
size_2 = zeros(length(N_values),length(delta_values));

densityProfile1D_comp1 = zeros(length(N_values),Nx);
densityProfile1D_comp2 = zeros(length(N_values),Nx);

j = 1% scan of densities (N)
for N = N_values
    N
    i=1; %scan of detunings
    for delta = delta_values
        %delta 
        %-----------------------------------------------------------
        % Setting the data
        %-----------------------------------------------------------
        
        
        
        %% Setting the physical problem
    
        Om = 5000;  % Rabi frequency
        wr = 200;   % radial trap frequency
        wz = 50;    % axial trap frequncy
    
        %N = 200*10^6;
        a_bohr = 52.917721 * 10^(-11);
        a11 = (96.3657*a_bohr); 
        a22 = (33.0784*a_bohr);
        a12 = (-53.5173*a_bohr);
    
        % adimensional parameters
        T = 1/Om;             % characteristic time
        L = 1/sqrt(wr);       % characteristic legth
        gamma_x = wr/Om;      % adimensional trap frequncies along x and y
        gamma_y = wr/Om;
        Rabi = Om/Om;         % adimensional Rabi frequency = 1 
        Delta = 0.5*T/(L^2);  % kinetic beta = T/L^2
        Beta = 1;             % constant in front of the interaction term
    
        % define INTERACTION Adimensional parameters g22, g22, g12
        g11 = (4*pi)* a11 * (N*T)/(L^3) / (sqrt(2*pi)*(1/sqrt(wz)));   % 4 pi a_s * (N T /L^3)
        g22 = (4*pi)* a22 * (N*T)/(L^3) / (sqrt(2*pi)*(1/sqrt(wz)));
        g12 = (4*pi)* a12 * (N*T)/(L^3) / (sqrt(2*pi)*(1/sqrt(wz)));
    
        %% Define H terms
    
        Potential{1,1} = @(X,Y) (1/2)*(gamma_x*X.^2+gamma_y*Y.^2) + (1/2)*(delta);
        Potential{1,2} = @(X,Y) (1/2)*Rabi;
        Potential{2,1} = @(X,Y) (1/2)*Rabi;
        Potential{2,2} = @(X,Y) (1/2)*(gamma_x*X.^2+gamma_y*Y.^2) - (1/2)*(delta);
    
        CoupledSpinNL{1,1} = @(Phi,X,Y) (g11*(abs(Phi{1})).^2 + g12*(abs(Phi{2})).^2);
        CoupledSpinNL{1,2} = @(Phi,X,Y) 0;
        CoupledSpinNL{2,1} = @(Phi,X,Y) 0;
        CoupledSpinNL{2,2} = @(Phi,X,Y) (g22*(abs(Phi{2})).^2 + g12*(abs(Phi{1})).^2);
    
        Physics2D = Physics2D_Var2d(Method, Delta, Beta); 
        Physics2D = Dispersion_Var2d(Method, Physics2D);
        Physics2D = Potential_Var2d(Method, Physics2D, Potential);
        Physics2D = Nonlinearity_Var2d(Method, Physics2D, CoupledSpinNL);
       
        
        %% Setting the initial data
        %InitialData_Choice = 1; % Gaussian
        X0 = 0;
        Y0 = 0;
        Phi_0 = InitialData_Var2d_custom(Method, Geometry2D, Physics2D, Rabi, delta, X0, Y0, gamma_x, gamma_y);
        %[Psi1_density, Psi2_density] = Print2D_Psi(Phi_0, Geometry2D); % 2D plot
        %[Density_1D_1, Density_1D_2] = Plot1D_DensityProfile(Phi_0, Geometry2D,'inital guess density profile');
        
        %% Setting informations and outputs
        Outputs = OutputsINI_Var2d(Method);
        Printing = 0;
        Evo = 300;
        Draw = 0;
        Print = Print_Var2d(Printing,Evo,Draw);
    
        %-----------------------------------------------------------
        % Launching computation
        %-----------------------------------------------------------
        
        [Phi_1, Outputs] = GPELab2d(Phi_0,Method,Geometry2D,Physics2D,Outputs,[],Print);

        % save outputs
        % spin composition
        P(j,i) = ComputeSpinPolarization(Phi_1, Geometry2D);
        [P_up(j,i), P_down(j,i)] = ComputeSpinProbabilities(Phi_1, Geometry2D);

        % save Energy and Chemical Potential values
        energy_1(j,i) = Outputs.Energy{1}(end);
        energy_2(j,i) = Outputs.Energy{2}(end);
        energy_tot(j,i) = energy_1(j,i) + energy_2(j,i);

        chemPot_1(j,i) = Outputs.Chemical_potential{1}(end);
        chemPot_2(j,i) = Outputs.Chemical_potential{1}(end);
        

        % Calculate Energies
        IE(j,i) = InteractionEnergy(Phi_1,g11,g22,g12,Geometry2D);
        RE(j,i) = RabiEnergy(Phi_1,delta,Rabi,Geometry2D);
        KE(j,i) = KineticEnergy(Phi_1,Delta,Geometry2D);
        PE(j,i) = PotentialEnergy(Phi_1,gamma_x,gamma_y,Geometry2D);


        % save size BEC
        size_1(j,i) = Outputs.x_rms{1}(end);
        size_2(j,i) = Outputs.y_rms{1}(end);

        if delta == 1

            %     1D density profiles
            %     Compute probability densities in 2D ground state solution
            Psi1_1_density = abs(Phi_1{1}).^2;  %component 1
            Psi1_2_density = abs(Phi_1{2}).^2;  %component 2
            
            %     Integrate along the Y-direction to get 1D profile along X
            Y_vals = Geometry2D.Y(:,1);
            Density1_1D_1 = trapz(Y_vals, Psi1_1_density, 1); % Integrate along Y
            Density1_1D_2 = trapz(Y_vals, Psi1_2_density, 1); % Integrate along Y
        end

        
        i=i+1;
    end
    
    if delta == 1
        % save 1D delsity
        densityProfile1D_comp1(j,:) = Density1_1D_1;
        densityProfile1D_comp2(j,:) = Density1_1D_2;
    end

    
    j=j+1;
end

outputFolder = 'C:\Users\sarat\OneDrive\Documenti\InstOptique\Simulations\GPELab\outputs\results\';  % Change this to your desired folder name
fileName = 'output_data.mat';

% Check if the folder exists, if not, create it
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Save the variables into the .mat file inside the folder
save(fullfile(outputFolder, fileName), 'N_values','delta_values', ...
    'densityProfile1D_comp1', 'densityProfile1D_comp2', ...
    'energy_tot','PE','KE',"IE",'RE', ...
    'P','P_down','P_up', ...
    'size_1','size_2');