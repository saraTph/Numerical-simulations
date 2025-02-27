% author: Sara Tiengo 
% date: 30/01/2025

%% GROUND STATE COMPUTATION OF A MULTICOMPONENT BEC WITH COUPLED NONLINEARITIES

%% Setting the method and geometry
Computation = 'Ground';
Ncomponents = 2;
Type = 'BESP'; % method to solve the Continuous Normalized Gradient Flow (CNGF)
Deltat = 1e-1;
Stop_time = [];
Stop_crit = {'MaxNorm',1e-4}; % epsilon
Method = Method_Var2d(Computation, Ncomponents, Type, Deltat, Stop_time, Stop_crit);
xmin = -10; % computational domine dimensions
xmax = 10;
ymin = -10;
ymax = 10;
Nx = 2^7+1; % # grid points
Ny = 2^7+1;
Geometry2D = Geometry2D_Var2d(xmin,xmax,ymin,ymax,Nx,Ny);
Deltax = (xmax - xmin)/Nx;

%% set physical quantities

a_bohr = 0.52917721e-10;
a11 = ((86.4014+1.65)*a_bohr); 
a22 = (33.2755*a_bohr);
a12 = ((-53.1022 + 0* 2 * 0.7843132831992122)*a_bohr);

wr = 169*2*pi;   % radial trap frequency
wz = 26*2*pi;    % axial trap frequncy

% define scan values for delta and Omega or n1D
delta_values = linspace(-2,2,40);
delta_values = [-8 -7 -5 -4 -3 -2.5 delta_values 2.5 3 4 5 7 8];
%delta_values = [-8 0.25 8];
%Omega_values = linspace(4*wr,20*wr,10);
Omega_values = 2810*2*pi;
n_values = [4186709333.3333 3539989333.3333 2411669333.3333 1131989333.3333 196309333.3333];
%n_values = [3894613333.3333, 3293013333.3333,  2243413333.3333,1053013333.3333, 182613333.3333];

%% initialize output vectors
P = zeros(length(Omega_values),length(delta_values));
P_up = zeros(length(Omega_values),length(delta_values));
P_down = zeros(length(Omega_values),length(delta_values));

energy_1 = zeros(length(Omega_values),length(delta_values));
energy_2 = zeros(length(Omega_values),length(delta_values));
energy_tot = zeros(length(Omega_values),length(delta_values));

IE = zeros(length(Omega_values),length(delta_values));
RE = zeros(length(Omega_values),length(delta_values));
PE = zeros(length(Omega_values),length(delta_values));
KE = zeros(length(Omega_values),length(delta_values));


%% evaluate groung state
j = 1;% scan of densities (N)
for n1D = n_values
    n1D
    Om = Omega_values
    i=1; %scan of detunings
    for delta = delta_values
        delta 
        %-----------------------------------------------------------
        % Setting the data
        %-----------------------------------------------------------
      
        %% characteristic lengths
        T = 1/Om;             % characteristic time
        L = sqrt(1/wr);       % characteristic legth
        Lz = sqrt(2/wz);

        % adimensional parameters
        gamma_x = wr/Om;      % adimensional trap frequncies along x and y
        gamma_y = wr/Om;
        Rabi = Om/Om;         % adimensional Rabi frequency = 1 
        Delta = 0.5*T/(L^2);  % kinetic beta = T/L^2
        Beta = 1;             % constant in front of the interaction term
    
        % define INTERACTION Adimensional parameters g22, g22, g12
        g11 = (n1D * (4*pi*a11)) * T/(L^2);
        g22 = (n1D * (4*pi*a22)) * T/(L^2);
        g12 = (n1D * (4*pi*a12)) * T/(L^2);
    
        %% Define H terms
    
        Potential{1,1} = @(X,Y) (1/2)*(gamma_x*X.^2+gamma_y*Y.^2) + (1/2)*(delta);
        Potential{1,2} = @(X,Y) (1/2)*Rabi;
        Potential{2,1} = @(X,Y) (1/2)*Rabi;
        Potential{2,2} = @(X,Y) (1/2)*(gamma_x*X.^2+gamma_y*Y.^2) - (1/2)*(delta);
    
        NL{1,1} = @(Phi,X,Y) (g11*(abs(Phi{1})).^2 + g12*(abs(Phi{2})).^2);
        NL{1,2} = @(Phi,X,Y) 0;
        NL{2,1} = @(Phi,X,Y) 0;
        NL{2,2} = @(Phi,X,Y) (g22*(abs(Phi{2})).^2 + g12*(abs(Phi{1})).^2);

        NLE{1,1} = @(Phi,X,Y) (1/2)*(g11*(abs(Phi{1})).^2 + g12*(abs(Phi{2})).^2);
        NLE{1,2} = @(Phi,X,Y) 0;
        NLE{2,1} = @(Phi,X,Y) 0;
        NLE{2,2} = @(Phi,X,Y) (1/2)*(g22*(abs(Phi{2})).^2 + g12*(abs(Phi{1})).^2);
    
        Physics2D = Physics2D_Var2d(Method, Delta, Beta); 
        Physics2D = Dispersion_Var2d(Method, Physics2D);
        Physics2D = Potential_Var2d(Method, Physics2D, Potential);
        Physics2D = Nonlinearity_Var2d(Method, Physics2D, NL, [], NLE);
       
        
        %% Setting the initial data
        %InitialData_Choice = 1; % Gaussian
        X0 = 0;
        Y0 = 0;
        Phi_0 = InitialData_Var2d_custom(Method, Geometry2D, Physics2D, Rabi, delta, X0, Y0, gamma_x, gamma_y);
        %[Psi1_density, Psi2_density] = Print2D_Psi(Phi_0, Geometry2D); % 2D plot
        %[Density_1D_1, Density_1D_2] = Plot1D_DensityProfile(Phi_0, Geometry2D,'inital guess density profile');
        
        %% Setting informations and outputs
        Outputs = OutputsINI_Var2d(Method);
        Printing = 1;
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
        

        % Calculate Energies
        IE(j,i) = InteractionEnergy(Phi_1,g11,g22,g12,Geometry2D);
        RE(j,i) = RabiEnergy(Phi_1,delta,Rabi,Geometry2D);
        KE(j,i) = KineticEnergy(Phi_1,Delta,Geometry2D);
        PE(j,i) = PotentialEnergy(Phi_1,gamma_x,gamma_y,Geometry2D);

%         if delta == 2
% 
%             %     1D density profiles
%             %     Compute probability densities in 2D ground state solution
%             Psi1_1_density = abs(Phi_1{1}).^2;  %component 1
%             Psi1_2_density = abs(Phi_1{2}).^2;  %component 2
%             
%             %     Integrate along the Y-direction to get 1D profile along X
%             Y_vals = Geometry2D.Y(:,1);
%             Density1_1D_1 = trapz(Y_vals, Psi1_1_density, 1); % Integrate along Y
%             Density1_1D_2 = trapz(Y_vals, Psi1_2_density, 1); % Integrate along Y
%         end

        
        i=i+1;
    end

    
    j=j+1;
end

outputFolder = 'C:\Users\Sarah\Documents\GitHub\Numerical-simulations\GPELab\outputs';  % Change this to your desired folder name
fileName = 'output_data.mat';

% Check if the folder exists, if not, create it
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Save the variables into the .mat file inside the folder
save(fullfile(outputFolder, fileName), 'Omega_values','delta_values', ...
    'energy_tot','PE','KE',"IE",'RE', ...
    'P','P_down','P_up');
beep;